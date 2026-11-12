"""RRecon: read density, read overlaps and Tn5 junctions for Tn5-based sequencing.

Tn5 duplicates 9 bp of target sequence at each cut, so two fragments cut from
the same molecule at the same site share a 9 bp target-site duplication (TSD)
at their junction. RRecon finds every overlapping read pair per cell (and per
haplotype, for a whatshap-haplotagged BAM) and flags those junctions.

    import rrecon
    res = rrecon.run("cells.bam", "results/", hp_tag="HP", ps_tag="PS", threads=8)

Modules:
    api         -- run(), density(), overlaps(), tn5_junctions()
    schema      -- output tables and columns
    _core       -- the Rust engine (when built), _pyengine the Python fallback
    bam_io, overlaps, plotting -- the earlier per-barcode overlap tools
"""

__version__ = "0.2.0"

from .api import (  # noqa: E402
    Params,
    RunResult,
    density,
    overlaps,
    read_table,
    run,
    rust_available,
    tn5_junctions,
)

__all__ = [
    "Params",
    "RunResult",
    "density",
    "overlaps",
    "read_table",
    "run",
    "rust_available",
    "tn5_junctions",
]
