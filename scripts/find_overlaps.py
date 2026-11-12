#!/usr/bin/env python
"""Task 1 CLI: build a per-barcode overlapping-read table from a BAM file.

Usage:
    python scripts/find_overlaps.py --bam in.bam --out overlaps.tsv \
        [--barcode-tag CB] [--barcode AAAA [AAAA2 ...]] [--min-mapq 0] \
        [--processes 4] [--hts-threads 2] [--verbose]

Output TSV columns: barcode, chrom, start, end, n_overlaps
start/end bound only the overlapping sub-interval itself (not the full
span of any contributing read); n_overlaps is the read depth there.
Reads are not chained transitively -- A-overlaps-B-overlaps-C with no
base shared by all three yields two rows (A∩B, B∩C), not one A-to-C row.

Parallelism: work is split by chromosome (--processes), since overlaps are
only ever merged within a (barcode, chrom) group -- never across
chromosomes -- so each chromosome can be processed completely
independently. --hts-threads controls htslib's own decompression
parallelism per worker and stacks with --processes. GPUs don't apply here:
this is per-read parsing/tag-lookup work, not bulk numeric compute.

--barcode restricts output to one or more specific cells. NOTE: this saves
memory and merge time, not scan time -- a coordinate-sorted BAM has no
index by barcode, so every alignment record in the (selected) chromosomes
still has to be decompressed and its tag checked either way.

By default, only canonical chromosomes (1-22, X, Y) are scanned at all --
alt haplotypes, decoys, chrM, and unplaced/unlocalized scaffolds are
skipped before any work is dispatched, not filtered out after scanning.
This one DOES save scan time, unlike --barcode. Pass --all-contigs to
process every contig in the BAM instead.
"""

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import pysam

from rrecon.chroms import is_canonical_chrom
from rrecon.merge import OVERLAP_COLUMNS, process_chromosome


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
        "--barcode",
        nargs="+",
        default=None,
        help="Restrict to one or more specific cell barcodes (space-"
        "separated). Default: keep every barcode. Saves memory/merge time, "
        "not scan time -- see module docstring.",
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
        "--all-contigs",
        action="store_true",
        help="Process every contig in the BAM, including alt haplotypes, "
        "decoys, chrM, and unplaced/unlocalized scaffolds. Default: only "
        "canonical chromosomes 1-22, X, Y are scanned.",
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
        all_chroms = list(bam.references)

    if args.all_contigs:
        chroms = all_chroms
    else:
        chroms = [c for c in all_chroms if is_canonical_chrom(c)]
        skipped = [c for c in all_chroms if c not in chroms]
        if args.verbose and skipped:
            preview = ", ".join(skipped[:10]) + (" ..." if len(skipped) > 10 else "")
            print(
                f"Skipping {len(skipped)} non-canonical contig(s): {preview}",
                file=sys.stderr,
            )
        if not chroms:
            raise SystemExit(
                "None of this BAM's contigs matched the canonical chr1-22/X/Y "
                f"naming (found: {all_chroms[:10]}...). Re-run with "
                "--all-contigs if this reference uses different names."
            )

    merged_parts = []

    if args.processes <= 1:
        for chrom in chroms:
            merged, stats = process_chromosome(
                args.bam,
                chrom,
                args.barcode_tag,
                args.min_mapq,
                args.hts_threads,
                args.verbose,
                args.barcode,
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
                    args.barcode,
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
            f"min_mapq={args.min_mapq}, --barcode filter={args.barcode}). "
            "Most likely cause: --barcode-tag doesn't match your BAM's "
            "actual tag name, or --barcode names a barcode that isn't in "
            "this file. Check with `samtools view <bam> | head -1` and "
            "look at the tag fields. Re-run with --verbose to see the "
            "per-chromosome breakdown."
        )

    merged.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(merged)} merged intervals to {args.out}")


if __name__ == "__main__":
    main()
