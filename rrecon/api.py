"""Public API: one BAM pass producing read density, read overlaps and Tn5 junctions.

    import rrecon
    res = rrecon.run("cells.haplotagged.bam", "results/", hp_tag="HP", ps_tag="PS", threads=8)
    res.table("tn5_junctions")          # pandas DataFrame

    rrecon.density(bam, out_dir, ...)       # density_segments + density_bins
    rrecon.overlaps(bam, out_dir, ...)      # overlaps
    rrecon.tn5_junctions(bam, out_dir, ...) # tn5_junctions + tn5_summary + tn5_kmers

The work is done by the Rust extension (rrecon._core) when it is built, else by
the pure-Python engine (rrecon._pyengine); both write identical tables.
"""

from __future__ import annotations

import json
import os
import time
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path

import pysam

from . import _pyengine
from .chroms import is_canonical_chrom
from .schema import BOOL_COLUMNS, COLUMNS, FLOAT_COLUMNS, INT_COLUMNS, OUTPUT_KINDS

try:  # the compiled core is optional
    from . import _core  # type: ignore[attr-defined]
except ImportError:  # pragma: no cover - depends on the build
    _core = None

__all__ = [
    "Params",
    "RunResult",
    "rust_available",
    "run",
    "density",
    "overlaps",
    "tn5_junctions",
    "read_table",
]

ALL_OUTPUTS = ("density", "overlaps", "junctions")


def rust_available() -> bool:
    return _core is not None


@dataclass
class Params:
    """Options shared by every command (see the README for what each does)."""

    barcode_tag: str = "CB"
    hp_tag: str | None = None
    ps_tag: str | None = None
    min_mapq: int = 0
    barcodes: Iterable[str] | None = None
    all_contigs: bool = False
    bin_size: int = 1_000_000
    junction_size: int = 9
    overlap_tolerance: int = 0
    max_mismatches: int = 1
    max_reads_per_group: int = 0
    max_seq_bp: int = 100
    threads: int = 1
    hts_threads: int = 1
    formats: Sequence[str] = ("tsv", "parquet")

    def validate(self) -> None:
        for name in ("barcode_tag", "hp_tag", "ps_tag"):
            v = getattr(self, name)
            if v is not None and len(v) != 2:
                raise ValueError(f"{name} must be a 2-character BAM tag, got {v!r}")
        if self.bin_size <= 0:
            raise ValueError("bin_size must be positive")
        if not 0 <= self.min_mapq <= 255:
            raise ValueError("min_mapq must be 0-255")
        for name in ("overlap_tolerance", "max_mismatches", "max_reads_per_group", "max_seq_bp"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be >= 0")
        bad = set(self.formats) - {"tsv", "parquet"}
        if bad or not self.formats:
            raise ValueError(f"formats must be tsv and/or parquet, got {list(self.formats)}")


@dataclass
class RunResult:
    out_dir: Path
    prefix: str
    engine: str
    stats: dict
    skipped_groups: list
    files: dict = field(default_factory=dict)  # kind -> {"tsv": Path, "parquet": Path}
    meta_path: Path | None = None
    seconds: float = 0.0

    def path(self, kind: str, fmt: str | None = None) -> Path:
        f = self.files[kind]
        if fmt:
            return f[fmt]
        return f.get("parquet") or f["tsv"]

    def table(self, kind: str):
        """Load one output table as a pandas DataFrame."""
        return read_table(self.path(kind))


def read_table(path):
    """Read an rrecon output (.parquet or .tsv.gz) into a DataFrame with proper types."""
    import pandas as pd

    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    kind = _kind_of(path)
    dtypes = {c: _pandas_type(c) for c in COLUMNS.get(kind, [])}
    df = pd.read_csv(
        path,
        sep="\t",
        dtype={c: t for c, t in dtypes.items() if t != "bool"},
        keep_default_na=False,
        na_values=[],
        quoting=3,
    )
    for c in BOOL_COLUMNS & set(df.columns):
        df[c] = df[c].map({"true": True, "false": False, True: True, False: False})
    return df


def _pandas_type(c: str) -> str:
    if c in INT_COLUMNS:
        return "int64"
    if c in FLOAT_COLUMNS:
        return "float64"
    if c in BOOL_COLUMNS:
        return "bool"
    return "str"


def _kind_of(path: Path) -> str:
    name = path.name
    for kind in COLUMNS:
        if f".{kind}." in name:
            return kind
    return ""


def ensure_index(bam: str | os.PathLike) -> None:
    b = str(bam)
    if not any(os.path.exists(b + ext) for ext in (".bai", ".csi")) and not os.path.exists(
        str(Path(b).with_suffix(".bai"))
    ):
        pysam.index(b)


def resolve_chroms(bam, all_contigs: bool = False, chroms: Iterable[str] | None = None):
    """(name, length) of the contigs to scan, in BAM header order."""
    with pysam.AlignmentFile(str(bam), "rb") as f:
        header = list(zip(f.references, f.lengths))
    if chroms is not None:
        wanted = set(chroms)
        out = [(c, n) for c, n in header if c in wanted]
        missing = wanted - {c for c, _ in out}
        if missing:
            raise ValueError(f"contigs not in the BAM header: {sorted(missing)}")
        return out
    if all_contigs:
        return header
    out = [(c, n) for c, n in header if is_canonical_chrom(c)]
    if not out:
        raise ValueError(
            "none of the BAM's contigs are canonical chromosomes (1-22, X, Y, with or "
            f"without 'chr'); found {[c for c, _ in header[:8]]}... Use all_contigs=True."
        )
    return out


def _pick_engine(engine: str) -> str:
    if engine == "auto":
        return "rust" if _core is not None else "python"
    if engine == "rust" and _core is None:
        raise RuntimeError(
            "the Rust extension is not built; run `maturin develop --release` or use "
            "engine='python'"
        )
    if engine not in ("rust", "python"):
        raise ValueError("engine must be auto, rust or python")
    return engine


def to_parquet(tsv_gz: Path, dest: Path, kind: str) -> Path:
    """Stream a TSV.gz output into Parquet with explicit column types."""
    import pyarrow as pa
    import pyarrow.csv as pv
    import pyarrow.parquet as pq

    cols = COLUMNS[kind]

    def arrow_type(c):
        if c in INT_COLUMNS:
            return pa.int64()
        if c in FLOAT_COLUMNS:
            return pa.float64()
        if c in BOOL_COLUMNS:
            return pa.bool_()
        return pa.string()

    schema = pa.schema([(c, arrow_type(c)) for c in cols])
    reader = pv.open_csv(
        tsv_gz,
        read_options=pv.ReadOptions(block_size=8 << 20),
        parse_options=pv.ParseOptions(delimiter="\t", quote_char=False),
        convert_options=pv.ConvertOptions(
            column_types={c: schema.field(c).type for c in cols},
            null_values=[],
            strings_can_be_null=False,
            true_values=["true"],
            false_values=["false"],
        ),
    )
    with pq.ParquetWriter(dest, schema, compression="zstd") as w:
        for batch in reader:
            w.write_batch(batch)
    return dest


def run(
    bam,
    out_dir,
    outputs: Iterable[str] = ALL_OUTPUTS,
    *,
    prefix: str | None = None,
    engine: str = "auto",
    chroms: Iterable[str] | None = None,
    params: Params | None = None,
    **kwargs,
) -> RunResult:
    """Scan `bam` once and write the requested outputs to `out_dir`.

    outputs: any of "density", "overlaps", "junctions".
    Keyword arguments set fields of Params (e.g. hp_tag="HP", threads=8).
    """
    p = params or Params()
    if kwargs:
        p = Params(**{**asdict(p), **kwargs})
    p.validate()
    outputs = list(dict.fromkeys(outputs))
    unknown = set(outputs) - set(ALL_OUTPUTS)
    if unknown or not outputs:
        raise ValueError(f"outputs must be among {ALL_OUTPUTS}, got {outputs}")
    which = _pick_engine(engine)

    bam = Path(bam)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = prefix or bam.name.removesuffix(".bam")
    ensure_index(bam)
    contigs = resolve_chroms(bam, p.all_contigs, chroms)
    barcodes = sorted(set(p.barcodes)) if p.barcodes is not None else None

    t0 = time.perf_counter()
    common = dict(
        barcode_tag=p.barcode_tag,
        hp_tag=p.hp_tag,
        ps_tag=p.ps_tag,
        min_mapq=p.min_mapq,
        bin_size=p.bin_size,
        junction_size=p.junction_size,
        overlap_tolerance=p.overlap_tolerance,
        max_mismatches=p.max_mismatches,
        max_reads_per_group=p.max_reads_per_group,
        max_seq_bp=p.max_seq_bp,
        write_density="density" in outputs,
        write_overlaps="overlaps" in outputs,
        write_junctions="junctions" in outputs,
        hts_threads=max(1, p.hts_threads),
    )
    if which == "rust":
        stats, skipped, files = _core.run(
            str(bam),
            str(out_dir),
            prefix,
            contigs,
            barcodes=barcodes,
            threads=max(1, p.threads),
            **common,
        )
    else:
        ep = _pyengine.EngineParams(
            barcodes=frozenset(barcodes) if barcodes is not None else None, **common
        )
        stats, skipped, files = _pyengine.run(bam, out_dir, prefix, contigs, ep, max(1, p.threads))
    seconds = time.perf_counter() - t0

    result = RunResult(
        out_dir=out_dir,
        prefix=prefix,
        engine=which,
        stats=dict(stats),
        skipped_groups=[
            {"barcode": b, "haplotype": h, "chrom": c, "n_reads": n} for b, h, c, n in skipped
        ],
        seconds=seconds,
    )
    for kind, tsv in files:
        tsv = Path(tsv)
        entry = {}
        if "parquet" in p.formats:
            entry["parquet"] = to_parquet(tsv, out_dir / f"{prefix}.{kind}.parquet", kind)
        if "tsv" in p.formats:
            entry["tsv"] = tsv
        else:
            tsv.unlink()
        result.files[kind] = entry
    result.meta_path = _write_meta(result, bam, outputs, contigs, p)
    return result


def _write_meta(res: RunResult, bam: Path, outputs, contigs, p: Params) -> Path:
    from . import __version__

    params = asdict(p)
    params["barcodes"] = None if p.barcodes is None else len(set(p.barcodes))
    params["formats"] = list(p.formats)
    meta = {
        "tool": "rrecon",
        "version": __version__,
        "engine": res.engine,
        "bam": bam.name,
        "prefix": res.prefix,
        "outputs": list(outputs),
        "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "seconds": round(res.seconds, 3),
        "params": params,
        "chroms": [{"name": c, "length": n} for c, n in contigs],
        "stats": res.stats,
        "skipped_groups": res.skipped_groups,
        "files": {k: {f: Path(v).name for f, v in e.items()} for k, e in res.files.items()},
    }
    path = res.out_dir / f"{res.prefix}.meta.json"
    path.write_text(json.dumps(meta, indent=2) + "\n")
    return path


def density(bam, out_dir, **kwargs) -> RunResult:
    """Depth segments (depth >= 1) and binned read counts per cell and haplotype."""
    return run(bam, out_dir, ("density",), **kwargs)


def overlaps(bam, out_dir, **kwargs) -> RunResult:
    """Every overlapping read pair within a (barcode, haplotype, chrom) group."""
    return run(bam, out_dir, ("overlaps",), **kwargs)


def tn5_junctions(bam, out_dir, **kwargs) -> RunResult:
    """Junction-sized overlaps (9 bp by default) with few mismatches, plus summaries."""
    return run(bam, out_dir, ("junctions",), **kwargs)
