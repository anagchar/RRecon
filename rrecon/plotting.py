"""Task 2: genome-wide scatter of one cell's overlap clusters.

Consumes the TSV produced by scripts/find_overlaps.py (Task 1) and, for a
single cell barcode, plots each merged overlap interval as a dot:

    x = genomic position (chromosomes concatenated, Manhattan-style)
    y = n_overlaps (how many reads contributed to that merged interval)
    dot size = the overlap's bp span (end - start)
    color = chromosome (visual separation only)

`min_overlap_bp` is the tunable knob for "only show overlaps above Nbp".
"""

from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd


def _chrom_sort_key(chrom: str):
    """Natural chromosome order: chr1, chr2, ..., chr10, ..., chrX, chrY, chrM."""
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name))
    tail_order = {"X": 1, "Y": 2, "M": 3, "MT": 3}
    return (1, tail_order.get(name.upper(), 99), name)


def load_chrom_sizes_from_fai(fai_path: str) -> dict:
    """Parse a samtools .fai index (chrom, length, ...) into {chrom: length}."""
    sizes = {}
    with open(fai_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            sizes[fields[0]] = int(fields[1])
    return sizes


def _chrom_extents(df: pd.DataFrame, chrom_sizes: Optional[dict] = None) -> dict:
    """{chrom: extent}. Uses true lengths from chrom_sizes when given;
    otherwise infers each chromosome's extent as the max `end` seen for it
    in df -- a rough stand-in that's fine for visual spacing but not for
    reporting exact coordinates."""
    chroms = df["chrom"].unique()
    if chrom_sizes:
        return {
            c: chrom_sizes.get(c, df.loc[df["chrom"] == c, "end"].max()) for c in chroms
        }
    return df.groupby("chrom")["end"].max().to_dict()


def compute_genome_offsets(df: pd.DataFrame, chrom_sizes: Optional[dict] = None) -> dict:
    """{chrom: cumulative_offset}, laying chromosomes out left-to-right in
    natural order so per-chromosome positions become one genome-wide axis."""
    extents = _chrom_extents(df, chrom_sizes)
    chroms = sorted(extents, key=_chrom_sort_key)
    offsets = {}
    cursor = 0
    for chrom in chroms:
        offsets[chrom] = cursor
        cursor += extents[chrom]
    return offsets


def plot_barcode_overlaps(
    df: pd.DataFrame,
    barcode: str,
    min_overlap_bp: int = 0,
    chrom_sizes: Optional[dict] = None,
    ax=None,
):
    """Scatter plot of one barcode's merged overlap intervals, genome-wide.

    Args:
        df: Task 1 overlaps table (columns: barcode, chrom, start, end, n_overlaps).
        barcode: which cell to plot.
        min_overlap_bp: drop merged intervals with (end - start) below this.
            The "show overlaps above 20bp" knob -- change the default here,
            or pass --min-overlap-bp on the CLI without touching this file.
        chrom_sizes: optional {chrom: length} (e.g. from a .fai) for
            accurate chromosome spacing; inferred from the data if omitted.
        ax: existing matplotlib Axes to draw into; creates one if None.

    Returns:
        the matplotlib Axes drawn into.
    """
    sub = df[df["barcode"] == barcode].copy()
    if sub.empty:
        raise ValueError(f"No rows for barcode '{barcode}' in the overlaps table.")

    sub["overlap_bp"] = sub["end"] - sub["start"]
    sub = sub[sub["overlap_bp"] >= min_overlap_bp]
    if sub.empty:
        raise ValueError(
            f"Barcode '{barcode}' has no merged intervals with overlap_bp >= "
            f"{min_overlap_bp}. Try lowering min_overlap_bp."
        )

    extents = _chrom_extents(df, chrom_sizes)
    offsets = compute_genome_offsets(df, chrom_sizes)
    sub["genome_x"] = sub["chrom"].map(offsets) + (sub["start"] + sub["end"]) / 2

    chroms_in_order = sorted(offsets, key=lambda c: offsets[c])
    palette = plt.get_cmap("tab20").colors
    color_of = {c: palette[i % len(palette)] for i, c in enumerate(chroms_in_order)}

    if ax is None:
        _, ax = plt.subplots(figsize=(14, 5))

    # Scale dot *area* (not radius) so size differences read honestly.
    size_scale = 3000 / max(sub["overlap_bp"].max(), 1)
    ax.scatter(
        sub["genome_x"],
        sub["n_overlaps"],
        s=sub["overlap_bp"] * size_scale + 10,
        c=[color_of[c] for c in sub["chrom"]],
        alpha=0.6,
        edgecolors="none",
    )

    for boundary in sorted(offsets.values())[1:]:
        ax.axvline(boundary, color="grey", linewidth=0.5, alpha=0.5)

    tick_positions = [offsets[c] + extents[c] / 2 for c in chroms_in_order]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(chroms_in_order, rotation=90, fontsize=7)

    ax.set_xlabel("Genomic position (chromosomes concatenated)")
    ax.set_ylabel("Number of overlapping reads (n_overlaps)")
    ax.set_title(
        f"Overlap clusters for barcode {barcode} "
        f"(dot size = overlap bp-span, min {min_overlap_bp}bp)"
    )
    return ax
