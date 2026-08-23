#!/usr/bin/env python
"""Task 1 CLI: build a per-barcode overlapping-read table from a BAM file.

Usage:
    python scripts/find_overlaps.py --bam in.bam --out overlaps.tsv \
        [--barcode-tag CB] [--min-mapq 0] [--processes 4] [--hts-threads 2] \
        [--verbose]

Output TSV columns: barcode, chrom, start, end, n_overlaps
(interval size = end - start, not stored separately)

Parallelism: work is split by chromosome (--processes), since overlaps are
only ever merged within a (barcode, chrom) group -- never across
chromosomes -- so each chromosome can be processed completely
independently. --hts-threads controls htslib's own decompression
parallelism per worker and stacks with --processes. GPUs don't apply here:
this is per-read parsing/tag-lookup work, not bulk numeric compute.
"""

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import pysam

from rrecon.overlaps import OVERLAP_COLUMNS, process_chromosome


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--bam", required=True, help="Input coordinate-sorted BAM file")
    p.add_argument("--out", required=True, help="Output TSV path")
    p.add_argument(
        "--barcode-tag",
        default="CB",
        help="BAM tag holding the cell barcode (default: CB). Check your "
        "file's actual tag names with `samtools view <bam> | head -1` "
        "if this default doesn't match.",
    )
    p.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum mapping quality to keep a read (default: 0)",
    )
    p.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Number of chromosomes to process in parallel (default: 1, "
        "i.e. sequential). A sensible value is min(n_chromosomes, n_cpus).",
    )
    p.add_argument(
        "--hts-threads",
        type=int,
        default=1,
        help="htslib decompression threads per worker (default: 1). Helps "
        "even with --processes 1.",
    )
    p.add_argument(
        "--verbose",
        action="store_true",
        help="Print per-chromosome read-filtering stats to stderr",
    )
    return p.parse_args()


def ensure_indexed(bam_path):
    """Per-chromosome fetch() needs an index; build one if it's missing."""
    if not (os.path.exists(bam_path + ".bai") or os.path.exists(bam_path + ".csi")):
        print(f"No index found for {bam_path} -- building one...", file=sys.stderr)
        pysam.index(bam_path)


def main():
    args = parse_args()
    ensure_indexed(args.bam)

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        chroms = list(bam.references)

    merged_parts = []

    if args.processes <= 1:
        for chrom in chroms:
            merged, stats = process_chromosome(
                args.bam, chrom, args.barcode_tag, args.min_mapq, args.hts_threads, args.verbose
            )
            if args.verbose:
                print(f"{chrom}: {stats}", file=sys.stderr)
            merged_parts.append(merged)
    else:
        with ProcessPoolExecutor(max_workers=args.processes) as pool:
            futures = {
                pool.submit(
                    process_chromosome,
                    args.bam,
                    chrom,
                    args.barcode_tag,
                    args.min_mapq,
                    args.hts_threads,
                    args.verbose,
                ): chrom
                for chrom in chroms
            }
            for fut in as_completed(futures):
                chrom = futures[fut]
                merged, stats = fut.result()
                if args.verbose:
                    print(f"{chrom}: {stats}", file=sys.stderr)
                merged_parts.append(merged)

    merged = (
        pd.concat(merged_parts, ignore_index=True)
        if merged_parts
        else pd.DataFrame(columns=OVERLAP_COLUMNS)
    )

    if merged.empty:
        raise SystemExit(
            f"No reads passed filtering (barcode tag='{args.barcode_tag}', "
            f"min_mapq={args.min_mapq}). Most likely cause: --barcode-tag "
            "doesn't match your BAM's actual tag name. Check with "
            "`samtools view <bam> | head -1` and look at the tag fields. "
            "Re-run with --verbose to see the per-chromosome breakdown."
        )

    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(merged)} merged intervals to {args.out}")


if __name__ == "__main__":
    main()
