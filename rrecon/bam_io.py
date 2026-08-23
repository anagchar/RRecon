"""Generic BAM I/O helpers.

Task-agnostic: this module knows how to pull clean, filtered alignment
records out of a BAM file. It has no opinion about what you do with them
afterwards -- rrecon.overlaps (Task 1) consumes it, and the fragment
reconstruction step (Task 2 proper, not the scatter plot) will reuse it too.
"""

from dataclasses import dataclass
from typing import Optional

import pysam


@dataclass
class ReadRecord:
    """One filtered alignment, reduced to what downstream tasks need."""

    barcode: str
    chrom: str
    start: int  # 0-based, inclusive (pysam reference_start)
    end: int  # 0-based, exclusive (pysam reference_end; CIGAR-aware)
    read_name: str


def iter_filtered_reads(
    bam_path: str,
    barcode_tag: str = "CB",
    min_mapq: int = 0,
    stats: Optional[dict] = None,
    region: Optional[str] = None,
    threads: int = 1,
):
    """Yield a ReadRecord for each alignment that passes the standard filter.

    Kept: primary, mapped, non-duplicate, non-secondary, non-supplementary
    alignments with mapping quality >= min_mapq and a `barcode_tag` tag
    present.

    Args:
        bam_path: path to a coordinate-sorted, indexed BAM file.
        barcode_tag: BAM tag holding the cell barcode. 10x convention is
            "CB" (corrected/filtered barcode), but ONT single-cell demux
            pipelines vary -- confirm this against your actual file with
            `samtools view <bam> | head -1` before trusting the default.
        min_mapq: minimum mapping quality to keep a read.
        stats: optional dict to update in place with filtering counts
            (total, unmapped, secondary_or_supplementary, duplicate,
            low_mapq, missing_barcode, kept). Pass `{}` to collect them.
        region: if given, restrict the scan to this contig only (e.g.
            "chr1"), using pysam's indexed `.fetch()` instead of a full
            linear scan. Requires a `.bai`/`.csi` index next to bam_path.
            This is what lets chromosomes be processed independently/in
            parallel -- leave None to scan the whole file.
        threads: htslib decompression threads for this file handle
            (pysam/htslib-native parallelism, separate from and
            complementary to Python-level multiprocessing across
            chromosomes).

    Yields:
        ReadRecord
    """
    if stats is not None:
        for key in (
            "total",
            "unmapped",
            "secondary_or_supplementary",
            "duplicate",
            "low_mapq",
            "missing_barcode",
            "kept",
        ):
            stats.setdefault(key, 0)

    with pysam.AlignmentFile(bam_path, "rb", threads=threads) as bam:
        read_iter = bam.fetch(region) if region is not None else bam
        for read in read_iter:
            if stats is not None:
                stats["total"] += 1

            if read.is_unmapped:
                if stats is not None:
                    stats["unmapped"] += 1
                continue
            if read.is_secondary or read.is_supplementary:
                if stats is not None:
                    stats["secondary_or_supplementary"] += 1
                continue
            if read.is_duplicate:
                if stats is not None:
                    stats["duplicate"] += 1
                continue
            if read.mapping_quality < min_mapq:
                if stats is not None:
                    stats["low_mapq"] += 1
                continue
            try:
                barcode = read.get_tag(barcode_tag)
            except KeyError:
                if stats is not None:
                    stats["missing_barcode"] += 1
                continue

            if stats is not None:
                stats["kept"] += 1

            yield ReadRecord(
                barcode=barcode,
                chrom=read.reference_name,
                start=read.reference_start,
                end=read.reference_end,
                read_name=read.query_name,
            )
