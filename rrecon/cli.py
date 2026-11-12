"""The `rrecon` command line.

    rrecon all       cells.bam -o results/ --hp-tag HP --ps-tag PS --threads 8
    rrecon density   cells.bam -o results/
    rrecon overlaps  cells.bam -o results/
    rrecon junctions cells.bam -o results/ --max-mismatches 1
    rrecon gui --dir results/
"""

from __future__ import annotations

import sys
from pathlib import Path

import click

from . import __version__
from .api import ALL_OUTPUTS, Params, rust_available, run


def _barcodes(value: str | None):
    """A file with one barcode per line, or a comma-separated list."""
    if not value:
        return None
    p = Path(value)
    if p.is_file():
        import gzip

        opener = gzip.open if p.suffix == ".gz" else open
        with opener(p, "rt") as f:
            return [ln.strip().split("\t")[0] for ln in f if ln.strip()]
    return [b for b in value.split(",") if b]


def scan_options(f):
    opts = [
        click.argument("bam", type=click.Path(exists=True, dir_okay=False, path_type=Path)),
        click.option("-o", "--out-dir", type=click.Path(file_okay=False, path_type=Path),
                     default=Path("rrecon_out"), show_default=True, help="Output folder."),
        click.option("--prefix", default=None, help="File name prefix [default: BAM name]."),
        click.option("--threads", "-t", type=int, default=1, show_default=True,
                     help="Chromosomes processed in parallel."),
        click.option("--hts-threads", type=int, default=1, show_default=True,
                     help="BGZF decompression threads per chromosome."),
        click.option("--barcode-tag", default="CB", show_default=True,
                     help="BAM tag with the cell barcode."),
        click.option("--hp-tag", default=None,
                     help="Haplotype tag from whatshap haplotag (usually HP). Reads with HP 1 or 2 "
                          "form haplotypes 1 and 2; all others are 'untagged'."),
        click.option("--ps-tag", default=None, help="Phase-set tag (usually PS)."),
        click.option("--min-mapq", type=int, default=0, show_default=True),
        click.option("--barcodes", default=None,
                     help="Allow-list: a file (one barcode per line, .gz ok) or a comma list."),
        click.option("--chroms", default=None,
                     help="Only these contigs, comma-separated (e.g. chr22 for a quick test)."),
        click.option("--all-contigs", is_flag=True,
                     help="Scan every contig, not only chromosomes 1-22, X, Y."),
        click.option("--bin-size", type=int, default=1_000_000, show_default=True,
                     help="Bin size for density_bins."),
        click.option("--junction-size", type=int, default=9, show_default=True,
                     help="Tn5 target-site duplication length."),
        click.option("--overlap-tolerance", type=int, default=0, show_default=True,
                     help="Accept junctions of junction-size ± this many bp."),
        click.option("--max-mismatches", type=int, default=1, show_default=True,
                     help="Most mismatches allowed in a junction (and for seq_match)."),
        click.option("--max-reads-per-group", type=int, default=0, show_default=True,
                     help="Skip pair-finding for (cell, haplotype, chrom) groups with more reads "
                          "than this (0 = no limit). Skipped groups are listed in meta.json."),
        click.option("--max-seq-bp", type=int, default=100, show_default=True,
                     help="Write seq1/seq2/mismatch_positions only for overlaps up to this size."),
        click.option("--format", "formats", type=click.Choice(["tsv", "parquet"]), multiple=True,
                     default=("tsv", "parquet"), show_default=True,
                     help="Output formats (repeat the option for both)."),
        click.option("--engine", type=click.Choice(["auto", "rust", "python"]), default="auto",
                     show_default=True, help="Rust core or the pure-Python fallback."),
    ]
    for o in reversed(opts):
        f = o(f)
    return f


def _run(outputs, bam, out_dir, prefix, engine, barcodes, chroms, **kw):
    params = Params(barcodes=_barcodes(barcodes), **kw)
    contigs = [c for c in chroms.split(",") if c] if chroms else None
    try:
        res = run(bam, out_dir, outputs, prefix=prefix, engine=engine, params=params,
                  chroms=contigs)
    except (ValueError, RuntimeError, OSError) as e:
        raise click.ClickException(str(e)) from e
    s = res.stats
    click.echo(
        f"rrecon {__version__} [{res.engine}] {res.seconds:.1f}s - "
        f"{s['kept']:,} of {s['total']:,} reads kept",
        err=True,
    )
    if res.skipped_groups:
        click.echo(f"  {len(res.skipped_groups)} group(s) over --max-reads-per-group skipped "
                   "(see meta.json)", err=True)
    for kind, e in res.files.items():
        for fmt, path in e.items():
            click.echo(f"  {kind:<17} {fmt:<8} {path}")
    click.echo(f"  {'meta':<17} {'json':<8} {res.meta_path}")
    return res


def _command(name, outputs, help_):
    @click.command(name=name, help=help_)
    @scan_options
    def cmd(**kw):
        _run(outputs, **kw)

    return cmd


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="rrecon")
def main():
    """Read density, read overlaps and Tn5 9 bp junctions from a (haplotagged) BAM."""


main.add_command(_command("density", ("density",),
                          "Depth segments and binned read counts per cell and haplotype."))
main.add_command(_command("overlaps", ("overlaps",),
                          "Every overlapping read pair per cell, haplotype and chromosome."))
main.add_command(_command("junctions", ("junctions",),
                          "Tn5 junction overlaps (9 bp by default) and per-cell summaries."))
main.add_command(_command("all", ALL_OUTPUTS, "All three outputs from one pass over the BAM."))


@main.command()
@click.option("--dir", "folder", type=click.Path(exists=True, file_okay=False, path_type=Path),
              default=None, help="Folder with rrecon outputs to open.")
@click.option("--host", default="127.0.0.1", show_default=True)
@click.option("--port", type=int, default=8766, show_default=True)
@click.option("--no-open", is_flag=True, help="Do not open a browser.")
def gui(folder, host, port, no_open):
    """Open overlapscope, the web viewer, on a results folder."""
    from .server import serve

    serve(folder, host=host, port=port, open_browser=not no_open)


@main.command()
def info():
    """Show the version and which engine is available."""
    click.echo(f"rrecon {__version__}")
    click.echo(f"rust core: {'yes' if rust_available() else 'no (pure-Python fallback)'}")
    click.echo(f"python: {sys.version.split()[0]}")


if __name__ == "__main__":
    main()
