#!/usr/bin/env python
"""Time the Python and Rust engines on a synthetic BAM, with 1 and N threads.

Usage:
    python scripts/benchmark.py [--cells 40] [--reads 5000] [--threads 8] [--repeat 3]

Builds a synthetic haplotagged BAM in a temporary folder (24 chromosomes,
lengths scaled down 20x from GRCh38), runs `rrecon all` with each engine and
thread count, checks that every engine wrote the same tables, and prints the
best wall time of --repeat runs.
"""

import argparse
import gzip
import os
import tempfile
import time
from pathlib import Path

import rrecon
from rrecon.synthetic import random_reads, write_bam

HG38_MB = [249, 242, 198, 190, 182, 171, 159, 145, 138, 134, 135, 133,
           114, 107, 102, 90, 83, 80, 59, 64, 47, 51, 156, 57]


def text(path):
    with gzip.open(path, "rt") as f:
        return f.read()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", type=int, default=40)
    ap.add_argument("--reads", type=int, default=5000, help="reads per cell")
    ap.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 1))
    ap.add_argument("--repeat", type=int, default=3)
    ap.add_argument("--scale", type=int, default=20)
    ap.add_argument("--skip-python", action="store_true", help="time only the Rust engine")
    args = ap.parse_args()

    names = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chroms = {c: mb * 1_000_000 // args.scale for c, mb in zip(names, HG38_MB)}
    runs = [] if args.skip_python else [("python", 1), ("python", args.threads)]
    if rrecon.rust_available():
        runs += [("rust", 1), ("rust", args.threads)]
    else:
        print("Rust extension not built: timing the Python engine only")

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        t0 = time.perf_counter()
        reads = random_reads(chroms, n_cells=args.cells, reads_per_cell=args.reads,
                             read_len=(60, 4000), junction_frac=0.12, seed=3)
        bam = write_bam(tmp / "bench.bam", reads, chroms)
        print(f"synthetic BAM: {len(reads):,} reads, {args.cells} cells, "
              f"{bam.stat().st_size / 1e6:.1f} MB, built in {time.perf_counter() - t0:.1f}s\n")

        results, reference = [], None
        for engine, threads in runs:
            best = float("inf")
            for k in range(args.repeat):
                out = tmp / f"{engine}{threads}_{k}"
                t = time.perf_counter()
                res = rrecon.run(bam, out, engine=engine, threads=threads, hp_tag="HP",
                                 ps_tag="PS", formats=("tsv",))
                best = min(best, time.perf_counter() - t)
            tables = {kind: text(res.path(kind, "tsv")) for kind in res.files}
            if reference is None:
                reference = tables
            same = tables == reference
            results.append((engine, threads, best, same, len(tables["overlaps"].splitlines()) - 1))

    base = results[0][2]
    print(f"{'engine':<8}{'threads':>8}{'seconds':>10}{'speed-up':>10}{'pairs':>10}  identical")
    for engine, threads, sec, same, n in results:
        print(f"{engine:<8}{threads:>8}{sec:>10.2f}{base / sec:>9.1f}x{n:>10,}  {'yes' if same else 'NO'}")
    if not all(r[3] for r in results):
        raise SystemExit("engines disagree")


if __name__ == "__main__":
    main()
