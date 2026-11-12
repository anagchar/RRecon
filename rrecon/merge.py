"""Depth>=2 overlap segments per (barcode, chrom): the original Task 1 sweep-line.

Used by scripts/find_overlaps.py. The read-pair overlaps and Tn5 junctions
live in rrecon.api (rrecon.overlaps() / rrecon.tn5_junctions()).
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
    """Find genuinely-overlapping sub-intervals within each (barcode, chrom) group.

    Args:
        df: DataFrame with columns barcode, chrom, start, end (one row
            per read).

    Returns:
        DataFrame with columns barcode, chrom, start, end, n_overlaps.
    """
    if df.empty:
        return pd.DataFrame(columns=OVERLAP_COLUMNS)

    group_keys = ["barcode", "chrom"]

    # Two sweep-line events per read: +1 where it starts, -1 where it ends.
    starts = df[group_keys + ["start"]].rename(columns={"start": "pos"})
    starts["delta"] = 1
    ends = df[group_keys + ["end"]].rename(columns={"end": "pos"})
    ends["delta"] = -1
    events = pd.concat([starts, ends], ignore_index=True)

    # Net delta at each distinct position in a group -- a read ending
    # exactly where another starts cancels out here (half-open coords:
    # no shared base, so no depth change at that boundary).
    events = (
        events.groupby(group_keys + ["pos"], as_index=False)["delta"]
        .sum()
        .sort_values(group_keys + ["pos"])
        .reset_index(drop=True)
    )

    events["depth"] = events.groupby(group_keys)["delta"].cumsum()

    # A new segment starts whenever depth differs from the previous event
    # in the same group (grouped shift(1) is NaN at each group's first
    # row, which always compares unequal, so group boundaries fall out
    # for free without a separate check).
    prev_depth = events.groupby(group_keys)["depth"].shift(1)
    starts_new_segment = events["depth"] != prev_depth
    segment_id = starts_new_segment.cumsum()

    # A segment's end is wherever the *next* event in the group occurs --
    # depth is constant right up to that position. The final event of
    # each group always has depth 0 (every read has ended by then), so
    # it's dropped by the depth>=2 filter below regardless of its
    # (undefined) end.
    events["seg_end"] = events.groupby(group_keys)["pos"].shift(-1)

    segments = (
        events.groupby(segment_id)
        .agg(
            barcode=("barcode", "first"),
            chrom=("chrom", "first"),
            start=("pos", "first"),
            end=("seg_end", "last"),
            n_overlaps=("depth", "first"),
        )
        .reset_index(drop=True)
    )

    overlaps = segments[segments["n_overlaps"] >= 2].reset_index(drop=True)
    # `end` picks up float64 from the NaN seg_end of each group's final
    # (always-dropped, depth-0) event; every surviving row has a real
    # value at this point, so this is just restoring the int dtype.
    overlaps["end"] = overlaps["end"].astype(int)
    overlaps["n_overlaps"] = overlaps["n_overlaps"].astype(int)
    return overlaps[OVERLAP_COLUMNS]


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

    Returns:
        (merged_df, stats_dict_or_None): a tuple so the same function
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
