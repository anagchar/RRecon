#!/usr/bin/env python
"""Make the overlapscope demo: a synthetic haplotagged BAM run through `rrecon all`.

Usage:
    python scripts/make_demo.py web/public/demo [--cells 24] [--reads 2500]

Writes <out>/demo.*.parquet, demo.meta.json and index.json (the file list the
GitHub Pages build fetches). Everything is synthetic: human chromosome names
with lengths scaled down 50x, random reads per cell and haplotype, and planted
9 bp Tn5 junction pairs (some with one mismatch).
"""

import argparse
import json
import tempfile
from pathlib import Path

import rrecon
from rrecon.synthetic import random_reads, write_bam

HG38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path)
    ap.add_argument("--cells", type=int, default=24)
    ap.add_argument("--reads", type=int, default=4000, help="reads per cell")
    ap.add_argument("--scale", type=int, default=50, help="shrink chromosome lengths by this")
    args = ap.parse_args()

    chroms = {c: n // args.scale for c, n in HG38.items()}
    args.out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        bam = write_bam(
            Path(tmp) / "demo.bam",
            random_reads(chroms, n_cells=args.cells, reads_per_cell=args.reads,
                         read_len=(60, 4000), junction_frac=0.12, seed=11),
            chroms,
        )
        res = rrecon.run(bam, args.out, prefix="demo", hp_tag="HP", ps_tag="PS",
                         threads=4, bin_size=500_000, formats=("parquet",))
    names = sorted(p.name for p in args.out.iterdir() if p.name.startswith("demo."))
    (args.out / "index.json").write_text(json.dumps(names) + "\n")
    print(f"demo [{res.engine}]: {res.stats['kept']:,} reads, files in {args.out}:")
    for n in names:
        print("  ", n, f"{(args.out / n).stat().st_size / 1e3:.0f} kB")


if __name__ == "__main__":
    main()
