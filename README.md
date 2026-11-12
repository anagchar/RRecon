# RRecon

Read density, pairwise read overlaps and Tn5 9 bp junctions from Tn5-based
single-cell sequencing, including whatshap-haplotagged BAMs.

Tn5 duplicates 9 bp of target sequence at each cut. Two fragments cut from the
same molecule at the same site therefore overlap by exactly those 9 bp. RRecon
finds every overlapping read pair per cell (and per haplotype), and flags the
pairs whose overlap is a 9 bp target-site duplication.

One pass over a coordinate-sorted, indexed BAM gives three outputs:

| output | files | what |
|---|---|---|
| **density** | `density_segments`, `density_bins` | read depth per cell and haplotype, as segments and binned counts |
| **overlaps** | `overlaps` | every overlapping read pair, with the bases of both reads over the overlap and their mismatches |
| **junctions** | `tn5_junctions`, `tn5_summary`, `tn5_kmers` | the 9 bp junction pairs, per-cell counts and the junction 9-mers |

It has a fast Rust core (multithreaded by chromosome), a pure-Python fallback
that gives identical results, a Python API, a CLI, and **overlapscope**, a web
viewer built on the same stack and style as [cnvscope](https://anagchar.github.io/cnvscope/).

## Install

Requirements: Python ≥ 3.10, and a Rust toolchain (<https://rustup.rs>) to build the core.

```bash
git clone https://github.com/anagchar/RRecon && cd RRecon
pip install maturin
maturin develop --release        # builds the Rust core into your environment
```

Or `pip install .`, which builds a wheel the same way. To bundle the web viewer
into the package, build it first (Node.js ≥ 20):

```bash
(cd web && npm ci && npm run build)   # writes rrecon/static
```

Without Rust, the package still works with the Python engine. Check which
engine is in use with `rrecon info`.

## Quick start

On a whatshap-haplotagged BAM (`whatshap haplotag … -o cells.haplotagged.bam`):

```bash
samtools index cells.haplotagged.bam      # rrecon indexes it for you if needed

rrecon all cells.haplotagged.bam -o results/ \
    --hp-tag HP --ps-tag PS \
    --barcode-tag CB --min-mapq 20 \
    --threads 8

rrecon gui --dir results/
```

`results/` then holds `cells.haplotagged.*.parquet`, the same tables as `*.tsv.gz`, and
`cells.haplotagged.meta.json`, which records the parameters, read filtering
counts and chromosome lengths.

### Commands

| command | writes |
|---|---|
| `rrecon all BAM` | all three outputs from one pass over the BAM |
| `rrecon density BAM` | `density_segments`, `density_bins` |
| `rrecon overlaps BAM` | `overlaps` |
| `rrecon junctions BAM` | `tn5_junctions`, `tn5_summary`, `tn5_kmers` |
| `rrecon gui [--dir DIR]` | opens overlapscope in the browser |
| `rrecon info` | version and engine |

### Options (scan commands)

| option | default | |
|---|---|---|
| `-o, --out-dir` | `rrecon_out` | output folder |
| `--prefix` | BAM name | file name prefix |
| `-t, --threads` | 1 | chromosomes processed in parallel |
| `--hts-threads` | 1 | BGZF decompression threads per chromosome |
| `--barcode-tag` | `CB` | cell barcode tag; reads without it are dropped |
| `--barcodes` | – | allow-list: a file (one per line, `.gz` ok) or `A,B,C` |
| `--hp-tag` | – | haplotype tag (`HP`). Reads with value 1 or 2 form haplotypes `1` and `2`; all others are `untagged`. Without it every read is `untagged`. |
| `--ps-tag` | – | phase-set tag (`PS`), reported per pair as `phase_block` |
| `--min-mapq` | 0 | minimum MAPQ |
| `--all-contigs` | off | also scan alt, decoy, unplaced contigs and chrM (default: 1–22, X, Y, with or without `chr`) |
| `--bin-size` | 1000000 | bin size for `density_bins` |
| `--junction-size` | 9 | Tn5 duplication length |
| `--overlap-tolerance` | 0 | accept junctions of junction-size ± this |
| `--max-mismatches` | 1 | most mismatches in a junction (and for `seq_match`) |
| `--max-reads-per-group` | 0 (off) | skip pair-finding for a (cell, haplotype, chromosome) group with more reads; listed in `meta.json` |
| `--max-seq-bp` | 100 | write `seq1`, `seq2` and `mismatch_positions` only for overlaps up to this size (always for junction-sized ones) |
| `--format` | `tsv` and `parquet` | repeat for both: `--format tsv --format parquet` |
| `--engine` | `auto` | `rust`, `python` or `auto` (Rust when built) |

Reads kept: mapped, primary, not duplicate, not secondary or supplementary,
MAPQ ≥ `--min-mapq`, with a barcode (and in the allow-list, if one is given),
and spanning at least one reference base. The BAM must be coordinate-sorted.

### Python

```python
import rrecon

res = rrecon.run("cells.haplotagged.bam", "results/",
                 hp_tag="HP", ps_tag="PS", min_mapq=20, threads=8)
res.engine                    # "rust" or "python"
res.stats                     # read filtering counts
j = res.table("tn5_junctions")  # pandas DataFrame

rrecon.density(bam, out_dir, bin_size=500_000)
rrecon.overlaps(bam, out_dir, barcodes=["AAACGA…-1"])
rrecon.tn5_junctions(bam, out_dir, max_mismatches=0, engine="python")
rrecon.read_table("results/x.overlaps.tsv.gz")  # typed DataFrame from either format
```

## How overlaps and junctions are defined

Reads are grouped by **(barcode, haplotype, chromosome)**: reads of different
cells or haplotypes are never paired. Within a group, two reads overlap when
they share at least one reference base. Coordinates are 0-based and half-open,
as in BED, so reads that only touch do not overlap.

For a pair, `read1` is the read that sorts first by (start, end, name), and

```
overlap = [read2.start, min(read1.end, read2.end))
```

The overlap type is `contained` when one read lies inside the other, and
`end_to_start` otherwise. For a contained read the overlap is the inner read's
own span. (The earlier `find_tn5_junctions.py` took `read1.end - read2.start`
here, which is too long.)

**Bases and mismatches.** Each read's CIGAR is walked to find the base it
places on every reference position of the overlap:

- soft clips are skipped
- a deletion (or `N` skip) gives `-`
- insertions are not on the reference, so they are left out of `seq1`/`seq2`
  but set `has_indel_in_overlap` when they fall strictly inside the overlap

A mismatch is any position where the two reads' bases differ, including a base
against `-`. `mismatch_positions` are offsets from the overlap start (0-based).

**Tn5 junction:** `|overlap_bp − junction_size| ≤ overlap_tolerance` and
`n_mismatches ≤ max_mismatches`. `seq_match` additionally requires no indel.

## Output columns

All tables have a header line. Missing values are `.`.

**overlaps** and **tn5_junctions** (same columns)

| column | |
|---|---|
| `barcode`, `haplotype` | cell, and `1`, `2` or `untagged` |
| `phase_block` | the reads' PS value if they share it, else `ps1/ps2` (`.` without `--ps-tag`) |
| `chrom` | |
| `read1_name`, `read1_start`, `read1_end`, `read1_strand`, `read1_mapq` | first read (start 0-based, end exclusive, `+`/`-`) |
| `read2_name`, `read2_start`, `read2_end`, `read2_strand`, `read2_mapq` | second read |
| `overlap_start`, `overlap_end`, `overlap_bp` | the shared stretch |
| `overlap_type` | `end_to_start` or `contained` |
| `seq1`, `seq2` | each read's bases over the overlap (`-` = deletion); `.` above `--max-seq-bp` |
| `n_mismatches` | positions where `seq1` and `seq2` differ |
| `mismatch_positions` | their offsets, comma-separated; `.` when none |
| `has_indel_in_overlap` | `true` if either read has an insertion or deletion inside the overlap |
| `seq_match` | `n_mismatches ≤ max_mismatches` and no indel |

**density_segments**: `barcode, haplotype, chrom, start, end, depth`. These are
stretches of constant read depth ≥ 1, with neighbouring stretches of equal
depth merged.

**density_bins**: `barcode, haplotype, chrom, bin_start, bin_end, n_reads`.
Reads are counted in the bin holding their start. Empty bins are not written.

**tn5_summary**: `barcode, haplotype, n_reads, n_overlaps, n_end_to_start,
n_contained, n_junctions, pct_junctions`. Each row covers one cell and
haplotype across all chromosomes; `pct_junctions` is the percentage of
overlaps that are junctions.

**tn5_kmers**: `barcode, haplotype, junction_seq, count`, where
`junction_seq` is read 1's bases over the junction.

**meta.json**: version, engine, parameters, chromosome lengths, read
filtering counts, skipped groups, and the file names.

## overlapscope (the viewer)

`rrecon gui --dir results/` serves the viewer and opens a results folder.
You can also drop the files (all at once: `*.parquet` or `*.tsv.gz`, plus
`*.meta.json`) on its start screen. The files are read in the browser and
nothing is uploaded. The same build runs as a static site on GitHub Pages
with a synthetic demo.

The side panel selects a **cell**, a **haplotype** (all, HP1, HP2,
untagged) and a **chromosome or region**. Every view shows the three
haplotypes as separate lanes or series.

- **Overview**
  - Read density along the genome; at ≤ 5 Mb with one cell, read depth
    instead.
  - Overlap locations, with y = overlap size (log scale), coloured by type
    and with junctions highlighted. Click a dot to view its reads.
  - Overlap-size histogram (log y), with the junction size marked.
  - Drag on a track to zoom; double-click or press `Esc` to return to the
    genome.
- **Tn5 junctions**
  - Totals.
  - Junctions per cell, stacked by haplotype. Click a bar to select the cell.
  - Junction 9-mer logo (frequency or bits) and the most common 9-mers.
  - Mismatch counts and mismatch positions of the junction-sized overlaps.
- **Region**: an IGV-style view.
  - A coverage track per haplotype, then the reads packed in rows, with
    strand arrows.
  - Each overlap is shaded on the reads, and mismatches are coloured by base.
  - At base level, the selected pair's bases are drawn on the reads.
  - Scroll to zoom and drag to pan. Below the reads are the pair's details
    and an aligned `seq1`/`seq2`.
- **Table**: every pair in the selection.
  - Sort by any column.
  - Filter by size, mismatches, type, indels, phase block, or junctions only.
  - Export the filtered rows as CSV.
  - Click a row to open it in Region.
- **Dataset**: files, parameters, read filtering and skipped groups.

Keys: `/` searches cells, `↑`/`↓` step through cells, `1`–`5` switch tabs,
`t` toggles the theme, `Esc` returns to the genome view.

The browser holds all rows in memory. For very large runs, make the overlaps
file smaller with `--barcodes` or `--max-seq-bp`, or load only the junction
files.

### Developing the viewer

```bash
cd web
npm ci
npm run demo                 # synthetic demo into web/public/demo (needs rrecon)
npm run dev -- --mode pages  # http://localhost:5173, the static site with the demo
npm run dev                  # against `rrecon gui --dir … --port 8766`
npm run lint && npm test
npm run build                # into rrecon/static, for the Python package
npm run build:pages          # the GitHub Pages site into web/dist-pages
```

The *Pages* workflow (`.github/workflows/pages.yml`, run by hand) publishes the
static site.

## Tests and benchmark

```bash
cd rust && cargo test && cd ..   # Rust unit tests
pytest                           # both engines on a synthetic BAM, and Rust == Python
python scripts/benchmark.py      # Python vs Rust, 1 and N threads
```

The tests build a BAM with pysam in a temporary folder (`rrecon/synthetic.py`)
covering:

- contained reads
- an exact 9 bp junction
- 1 and 2 mismatches
- a deletion and an insertion inside the overlap
- reads that only touch
- HP1, HP2 and untagged reads kept apart
- duplicate, supplementary and secondary reads
- missing barcodes
- soft clips
- the group limit

The parity tests compare every output file of the two engines byte for byte.
No real data is used.

`scripts/benchmark.py` on a 16-core Mac (40 cells, 224k reads of 60–4000 bp,
84 MB BAM):

| engine | threads | seconds | speed-up |
|---|---|---|---|
| python | 1 | 5.28 | 1.0× |
| python | 8 | 1.33 | 4.0× |
| rust | 1 | 1.36 | 3.9× |
| rust | 8 | 0.52 | 10.2× |

CI (`.github/workflows/ci.yml`) runs `cargo test` and clippy, pytest on Linux
and macOS, the web lint, tests and build, and builds maturin wheels for Linux
and macOS (x86_64 and arm64) with the viewer bundled.

## Earlier tools

`scripts/find_overlaps.py` and `scripts/plot_overlaps.py` still work. They
use depth ≥ 2 sweep-line segments per barcode (`rrecon.merge`, formerly
`rrecon.overlaps`), not read pairs.
