"""Generic BAM I/O helpers.

Task-agnostic: this module knows how to pull clean, filtered alignment
records out of a BAM file. It has no opinion about what you do with them
afterwards -- rrecon.overlaps (Task 1) consumes it, and the fragment
reconstruction step (Task 2 proper, not the scatter plot) will reuse it too.
"""

from dataclasses import dataclass
from typing import Iterable, Optional, Union

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
    barcodes: Optional[Union[str, Iterable[str]]] = None,
    canonical_chroms: Optional[Iterable[str]] = None,
):
    """Yield a ReadRecord for each alignment that passes the standard filter.

    Kept: primary, mapped, non-duplicate, non-secondary, non-supplementary
    alignments with mapping quality >= min_mapq and a `barcode_tag` tag
    present (and, if `barcodes` is given, matching one of those barcodes).

    Args:
        bam_path: path to a coordinate-sorted, indexed BAM file.
        barcode_tag: BAM tag holding the cell barcode. 10x convention is
            "CB" (corrected/filtered barcode), but ONT single-cell demux
            pipelines vary -- confirm this against your actual file with
            `samtools view <bam> | head -1` before trusting the default.
        min_mapq: minimum mapping quality to keep a read.
        stats: optional dict to update in place with filtering counts
            (total, unmapped, secondary_or_supplementary, duplicate,
            low_mapq, missing_barcode, other_barcode, kept). Pass `{}` to
            collect them.
        region: if given, restrict the scan to this contig only (e.g.
            "chr1"), using pysam's indexed `.fetch()` instead of a full
            linear scan. Requires a `.bai`/`.csi` index next to bam_path.
            This is what lets chromosomes be processed independently/in
            parallel -- leave None to scan the whole file.
        threads: htslib decompression threads for this file handle
            (pysam/htslib-native parallelism, separate from and
            complementary to Python-level multiprocessing across
            chromosomes).
        barcodes: if given (a single barcode string, or any iterable of
            them), only reads carrying one of these barcodes are kept --
            everything else is dropped after its tag is read. NOTE: this
            does NOT skip decompressing/parsing non-matching reads --
            a coordinate-sorted BAM has no index by barcode, so every
            record still has to be inspected. It only reduces what gets
            held in memory and merged downstream. Leave None to keep
            every barcode (current default behaviour).
        canonical_chroms: if given, reads whose chromosome isn't in this
            set are dropped (see rrecon.chroms.DEFAULT_CANONICAL_CHROMS --
            excludes chrM, alt haplotypes, decoys, unplaced scaffolds).
            When scanning with `region` already set to a single canonical
            chromosome (the normal find_overlaps.py path), this is a
            no-op safety net; it matters when scanning without `region`
            (a full linear pass) or calling this generically. Leave None
            to keep every chromosome present in the file.

    Yields:
        ReadRecord
    """
    if isinstance(barcodes, str):
        barcodes = {barcodes}
    elif barcodes is not None:
        barcodes = set(barcodes)

    if stats is not None:
        for key in (
            "total",
            "unmapped",
            "secondary_or_supplementary",
            "duplicate",
            "low_mapq",
            "missing_barcode",
            "other_barcode",
            "non_canonical_chrom",
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
            if canonical_chroms is not None and read.reference_name not in canonical_chroms:
                if stats is not None:
                    stats["non_canonical_chrom"] += 1
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

            if barcodes is not None and barcode not in barcodes:
                if stats is not None:
                    stats["other_barcode"] += 1
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
