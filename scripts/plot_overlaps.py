#!/usr/bin/env python
"""Task 2 CLI: genome-wide scatter of one cell's overlap clusters.

Usage:
    python scripts/plot_overlaps.py --overlaps overlaps.tsv --barcode AAAA \
        --out AAAA_overlaps.png [--min-overlap-bp 20] \
        [--bam original.bam | --fai genome.fa.fai]

Reads the TSV produced by scripts/find_overlaps.py. x-axis = genomic
position (chromosomes concatenated); y-axis = n_overlaps; dot size = the
overlap's bp span (end - start). --min-overlap-bp is the easily-tunable
knob for "only show overlaps bigger than Nbp" (e.g. --min-overlap-bp 20).
Underlying filter (primary / non-duplicate / non-secondary /
non-supplementary) was already applied upstream in find_overlaps.py.

Chromosome spacing on the x-axis: pass --bam (reads exact lengths straight
from the BAM header -- no .fai/reference FASTA needed) or --fai (a
samtools .fai index) for accurate spacing. If neither is given, spacing is
inferred from the data itself -- fine visually, not exact. If both are
given, --fai wins.
"""

import argparse

import matplotlib

matplotlib.use("Agg")  # headless: this script writes a file, it doesn't need a display
import matplotlib.pyplot as plt
import pandas as pd

from rrecon.plotting import (
    load_chrom_sizes_from_bam,
    load_chrom_sizes_from_fai,
    plot_barcode_overlaps,
)


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--overlaps", required=True, help="TSV from scripts/find_overlaps.py")
    p.add_argument("--barcode", required=True, help="Which cell barcode to plot")
    p.add_argument("--out", required=True, help="Output image path (.png/.pdf/.svg)")
    p.add_argument(
        "--min-overlap-bp",
        type=int,
        default=0,
        help="Only plot merged intervals with (end-start) >= this many bp "
        "(default: 0, i.e. no filter). The tunable 'overlaps above Nbp' knob.",
    )
    p.add_argument(
        "--bam",
        default=None,
        help="Optional path to the original BAM -- reads exact chromosome "
        "lengths straight from its header for accurate x-axis spacing. No "
        ".bai or reference FASTA needed, just the BAM.",
    )
    p.add_argument(
        "--fai",
        default=None,
        help="Optional samtools .fai index (samtools faidx genome.fa) for "
        "accurate chromosome spacing. Without --bam or --fai, spacing is "
        "inferred from the data itself -- fine visually, not exact. If "
        "both are given, --fai wins.",
    )
    return p.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.overlaps, sep="\t")

    if args.barcode not in df["barcode"].unique():
        top = df["barcode"].value_counts().head(10)
        raise SystemExit(
            f"Barcode '{args.barcode}' not found in {args.overlaps}.\n"
            f"Most frequent barcodes in this file:\n{top.to_string()}"
        )

    if args.fai:
        chrom_sizes = load_chrom_sizes_from_fai(args.fai)
    elif args.bam:
        chrom_sizes = load_chrom_sizes_from_bam(args.bam)
    else:
        chrom_sizes = None

    fig, ax = plt.subplots(figsize=(14, 5))
    plot_barcode_overlaps(
        df,
        args.barcode,
        min_overlap_bp=args.min_overlap_bp,
        chrom_sizes=chrom_sizes,
        ax=ax,
    )
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
