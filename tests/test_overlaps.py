"""Tests for rrecon.overlaps.merge_overlaps.

Pins down the two things Task 1's output is supposed to mean:
  - n_overlaps is the read *depth* at a sub-interval (how many reads
    cover that exact stretch), not the size of a whole transitively-
    chained cluster.
  - start/end bound *only* the overlapping sub-region itself, not the
    full span of any contributing read.
"""

import pandas as pd
import pytest

from rrecon.overlaps import OVERLAP_COLUMNS, merge_overlaps, records_to_frame
from rrecon.bam_io import ReadRecord


def reads(*rows, barcode="bc1", chrom="chr1"):
    """Build a merge_overlaps-input DataFrame from (start, end) pairs, all
    sharing one barcode/chrom unless the row itself overrides them via a
    4-tuple (barcode, chrom, start, end)."""
    out = []
    for row in rows:
        if len(row) == 2:
            start, end = row
            out.append((barcode, chrom, start, end))
        else:
            out.append(row)
    return pd.DataFrame(out, columns=["barcode", "chrom", "start", "end"])


def as_records(df):
    """Round-trip through DataFrame -> tuple set for order-independent
    comparison of merge_overlaps output."""
    return set(
        df[["barcode", "chrom", "start", "end", "n_overlaps"]].itertuples(
            index=False, name=None
        )
    )


def test_empty_input_returns_empty_with_correct_columns():
    out = merge_overlaps(pd.DataFrame(columns=["barcode", "chrom", "start", "end"]))
    assert list(out.columns) == OVERLAP_COLUMNS
    assert len(out) == 0


def test_single_read_has_no_overlap():
    # One read, nothing to overlap with -> no output rows at all.
    out = merge_overlaps(reads((100, 200)))
    assert len(out) == 0


def test_two_overlapping_reads_report_only_the_shared_region():
    # A(100,300), B(150,400) share bases [150,300) only.
    out = merge_overlaps(reads((100, 300), (150, 400)))
    assert as_records(out) == {("bc1", "chr1", 150, 300, 2)}
    # Explicitly not the union of the two reads (100-400).
    row = out.iloc[0]
    assert (row["start"], row["end"]) != (100, 400)


def test_chain_without_triple_overlap_splits_into_two_rows():
    # A(100,200), B(150,300), C(250,400): A overlaps B, B overlaps C, but
    # A and C never share a base. Must NOT collapse into one 100-400
    # cluster with n_overlaps=3.
    out = merge_overlaps(reads((100, 200), (150, 300), (250, 400)))
    assert as_records(out) == {
        ("bc1", "chr1", 150, 200, 2),
        ("bc1", "chr1", 250, 300, 2),
    }


def test_touching_reads_do_not_count_as_overlapping():
    # Half-open [start, end): a read ending at 200 and one starting at
    # 200 share zero bases.
    out = merge_overlaps(reads((100, 200), (200, 300)))
    assert len(out) == 0


def test_nested_read_produces_two_boundaries():
    # B is fully contained inside A: overlap is exactly B's span.
    out = merge_overlaps(reads((100, 500), (200, 300)))
    assert as_records(out) == {("bc1", "chr1", 200, 300, 2)}


def test_triple_overlap_region_reports_depth_three():
    # A(100,300), B(150,350), C(200,250): the middle is covered by all
    # three, flanked by depth-2 regions on either side.
    out = merge_overlaps(reads((100, 300), (150, 350), (200, 250)))
    assert as_records(out) == {
        ("bc1", "chr1", 150, 200, 2),
        ("bc1", "chr1", 200, 250, 3),
        ("bc1", "chr1", 250, 300, 2),
    }


def test_overlaps_do_not_cross_barcode_boundary():
    df = reads(
        ("bc1", "chr1", 100, 300),
        ("bc2", "chr1", 150, 400),  # would overlap bc1's read if barcode ignored
    )
    out = merge_overlaps(df)
    assert len(out) == 0


def test_overlaps_do_not_cross_chrom_boundary():
    df = reads(
        ("bc1", "chr1", 100, 300),
        ("bc1", "chr2", 150, 400),  # would overlap on chr1's read if chrom ignored
    )
    out = merge_overlaps(df)
    assert len(out) == 0


def test_multiple_independent_groups_each_report_correctly():
    df = pd.concat(
        [
            reads((100, 300), (150, 400), barcode="bc1", chrom="chr1"),
            reads((10, 20), barcode="bc2", chrom="chr1"),  # singleton, no overlap
            reads((500, 700), (600, 800), barcode="bc1", chrom="chr2"),
        ],
        ignore_index=True,
    )
    out = merge_overlaps(df)
    assert as_records(out) == {
        ("bc1", "chr1", 150, 300, 2),
        ("bc1", "chr2", 600, 700, 2),
    }


def test_row_order_in_input_does_not_matter():
    df = reads((250, 400), (100, 200), (150, 300))  # deliberately unsorted
    out = merge_overlaps(df)
    assert as_records(out) == {
        ("bc1", "chr1", 150, 200, 2),
        ("bc1", "chr1", 250, 300, 2),
    }


def test_output_dtypes_are_integer():
    out = merge_overlaps(reads((100, 300), (150, 400)))
    assert out["start"].dtype.kind == "i"
    assert out["end"].dtype.kind == "i"
    assert out["n_overlaps"].dtype.kind == "i"


def test_records_to_frame_round_trip_matches_merge_overlaps():
    records = [
        ReadRecord(barcode="bc1", chrom="chr1", start=100, end=300, read_name="r1"),
        ReadRecord(barcode="bc1", chrom="chr1", start=150, end=400, read_name="r2"),
    ]
    df = records_to_frame(records)
    out = merge_overlaps(df)
    assert as_records(out) == {("bc1", "chr1", 150, 300, 2)}


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
