"""Task 1: per-barcode overlapping-read clustering.

Groups aligned reads by (cell barcode, chromosome), merges overlapping
genomic intervals within each group, and counts how many reads contributed
to each merged interval. Equivalent to running `bedtools merge -c 1 -o
count` separately within every cell barcode.
"""

from typing import Iterable, Optional

import pandas as pd

from rrecon.bam_io import ReadRecord, iter_filtered_reads

OVERLAP_COLUMNS = ["barcode", "chrom", "start", "end", "n_overlaps"]


def records_to_frame(records: Iterable[ReadRecord]) -> pd.DataFrame:
    """Convert an iterable of ReadRecord into a flat DataFrame."""
    return pd.DataFrame(
        ((r.barcode, r.chrom, r.start, r.end) for r in records),
        columns=["barcode", "chrom", "start", "end"],
    )


def merge_overlaps(df: pd.DataFrame) -> pd.DataFrame:
    """Merge overlapping intervals within each (barcode, chrom) group.

    Args:
        df: DataFrame with columns barcode, chrom, start, end (one row
            per read).

    Returns:
        DataFrame with columns barcode, chrom, start, end, n_overlaps.
        Each output row is a maximal cluster of mutually/transitively
        overlapping reads that share a barcode and chromosome; n_overlaps
        is how many input reads collapsed into it. Interval size is
        end - start (not stored separately).
    """
    if df.empty:
        return pd.DataFrame(columns=OVERLAP_COLUMNS)

    df = df.sort_values(["barcode", "chrom", "start"]).reset_index(drop=True)
    group_keys = ["barcode", "chrom"]

    # Running max of `end` seen in *earlier* rows of the same group.
    # (cummax includes the current row, so shift(1) drops back to
    # "before this row".)
    prev_end_cummax = df.groupby(group_keys)["end"].cummax().shift(1)

    # True at the first row of a (barcode, chrom) group.
    same_group_as_prev = df[group_keys].eq(df[group_keys].shift(1)).all(axis=1)

    # A new cluster starts when we've entered a new group, or when this
    # row's start doesn't reach back into anything accumulated so far.
    starts_new_cluster = (~same_group_as_prev) | (df["start"] > prev_end_cummax)
    cluster_id = starts_new_cluster.cumsum()

    merged = (
        df.groupby(cluster_id)
        .agg(
            barcode=("barcode", "first"),
            chrom=("chrom", "first"),
            start=("start", "min"),
            end=("end", "max"),
            n_overlaps=("start", "size"),
        )
        .reset_index(drop=True)
    )
    return merged


def process_chromosome(
    bam_path: str,
    chrom: str,
    barcode_tag: str = "CB",
    min_mapq: int = 0,
    hts_threads: int = 1,
    collect_stats: bool = False,
    barcodes=None,
):
    """Extract + merge overlaps for a single chromosome only.

    Safe to run independently/in parallel across chromosomes:
    merge_overlaps groups by (barcode, chrom), and restricting the input
    to one chromosome at a time means it can never merge across a
    chromosome boundary -- results are identical to running the
    whole-file version and then subsetting to this chrom. This is the
    function that gets handed to a process pool, one call per chromosome.

    `barcodes`: optional single barcode or iterable of barcodes -- if
    given, only those cells' reads are kept (see iter_filtered_reads for
    what this does and doesn't save you). None keeps every barcode.

    Returns:
        (merged_df, stats_dict_or_None) -- a tuple so the same function
        works whether called directly (sequential) or via
        ProcessPoolExecutor.submit (parallel); the return value has to be
        picklable either way, and a (DataFrame, dict) tuple is.
    """
    stats = {} if collect_stats else None
    records = iter_filtered_reads(
        bam_path,
        barcode_tag=barcode_tag,
        min_mapq=min_mapq,
        region=chrom,
        threads=hts_threads,
        stats=stats,
        barcodes=barcodes,
    )
    df = records_to_frame(records)
    merged = merge_overlaps(df)
    return merged, stats
