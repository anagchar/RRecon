"""Pure-Python engine: the fallback used when the Rust extension is not built.

It mirrors rust/src/{align,core,scan}.rs rule for rule, so both engines write
byte-identical (uncompressed) tables. Keep the two in step.
"""

from __future__ import annotations

import gzip
import os
import shutil
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import pysam

from .schema import (
    BIN_COLUMNS,
    KMER_COLUMNS,
    MISSING,
    OVERLAP_COLUMNS,
    SEGMENT_COLUMNS,
    SUMMARY_COLUMNS,
    UNTAGGED,
)

STAT_KEYS = [
    "total",
    "unmapped",
    "secondary_or_supplementary",
    "duplicate",
    "low_mapq",
    "missing_barcode",
    "other_barcode",
    "no_alignment_span",
    "kept",
]

# CIGAR codes: M I D N S H P = X
_REF_OPS = {0, 2, 3, 7, 8}
_ALIGNED = {0, 7, 8}


@dataclass
class EngineParams:
    barcode_tag: str = "CB"
    hp_tag: str | None = None
    ps_tag: str | None = None
    min_mapq: int = 0
    barcodes: frozenset | None = None
    bin_size: int = 1_000_000
    junction_size: int = 9
    overlap_tolerance: int = 0
    max_mismatches: int = 1
    max_reads_per_group: int = 0
    max_seq_bp: int = 100
    write_density: bool = True
    write_overlaps: bool = True
    write_junctions: bool = True
    hts_threads: int = 1

    def is_junction_size(self, bp: int) -> bool:
        return abs(bp - self.junction_size) <= self.overlap_tolerance

    def seq_limit(self) -> int:
        return max(self.max_seq_bp, self.junction_size + self.overlap_tolerance)


@dataclass
class Read:
    name: str
    start: int
    end: int
    reverse: bool
    mapq: int
    ps: str
    cigar: list = field(default_factory=list)
    seq: str = ""

    def key(self):
        return (self.start, self.end, self.name)


# --------------------------------------------------------------------------- alignment


def reference_span(cigar) -> int:
    return sum(n for op, n in cigar if op in _REF_OPS)


def aligned_window(ref_start: int, cigar, seq: str, ws: int, we: int) -> tuple[str, bool]:
    """The read's bases on reference window [ws, we): see rust/src/align.rs."""
    n = max(we - ws, 0)
    out = ["-"] * n
    insertion = False
    rpos, qpos = ref_start, 0
    for op, length in cigar:
        if rpos >= we:
            break
        if op in _ALIGNED:
            lo, hi = max(rpos, ws), min(rpos + length, we)
            for p in range(lo, hi):
                q = qpos + (p - rpos)
                out[p - ws] = seq[q].upper() if q < len(seq) else "N"
            rpos += length
            qpos += length
        elif op == 1:
            if ws < rpos < we:
                insertion = True
            qpos += length
        elif op in (2, 3):
            rpos += length
        elif op == 4:
            qpos += length
    return "".join(out), insertion


def overlap_window(r1: Read, r2: Read) -> tuple[int, int, bool]:
    """[r2.start, min(r1.end, r2.end)); r1 sorts first by (start, end, name)."""
    os_ = r2.start
    oe = min(r1.end, r2.end)
    contained = r1.start == r2.start or r2.end <= r1.end
    return os_, oe, contained


def pair_detail(r1: Read, r2: Read, os_: int, oe: int):
    s1, ins1 = aligned_window(r1.start, r1.cigar, r1.seq, os_, oe)
    s2, ins2 = aligned_window(r2.start, r2.cigar, r2.seq, os_, oe)
    positions = [k for k, (a, b) in enumerate(zip(s1, s2)) if a != b]
    has_indel = ins1 or ins2 or "-" in s1 or "-" in s2
    return s1, s2, positions, has_indel


def phase_block(ps1: str, ps2: str) -> str:
    return ps1 if ps1 == ps2 else f"{ps1}/{ps2}"


def _b(v: bool) -> str:
    return "true" if v else "false"


def format_pair(p, group, chrom, r1, r2, os_, oe, contained, s1, s2, positions, has_indel):
    bp = oe - os_
    show = bp <= p.seq_limit()
    nmm = len(positions)
    seq_match = nmm <= p.max_mismatches and not has_indel
    fields = [
        group[0],
        group[1],
        phase_block(r1.ps, r2.ps),
        chrom,
        r1.name,
        r1.start,
        r1.end,
        "-" if r1.reverse else "+",
        r1.mapq,
        r2.name,
        r2.start,
        r2.end,
        "-" if r2.reverse else "+",
        r2.mapq,
        os_,
        oe,
        bp,
        "contained" if contained else "end_to_start",
        s1 if show else MISSING,
        s2 if show else MISSING,
        nmm,
        ",".join(map(str, positions)) if show and positions else MISSING,
        _b(has_indel),
        _b(seq_match),
    ]
    return "\t".join(map(str, fields)) + "\n"


# --------------------------------------------------------------------------- density


def depth_segments(intervals):
    net: dict[int, int] = defaultdict(int)
    for s, e in intervals:
        net[s] += 1
        net[e] -= 1
    events = sorted((pos, d) for pos, d in net.items() if d != 0)
    out, depth = [], 0
    for i, (pos, d) in enumerate(events):
        depth += d
        if depth > 0 and i + 1 < len(events):
            out.append((pos, events[i + 1][0], depth))
    return out


def format_pct(n: int, d: int) -> str:
    return "0.0000" if d == 0 else f"{100.0 * n / d:.4f}"


# --------------------------------------------------------------------------- reading


def _tag(read, tag):
    try:
        v = read.get_tag(tag)
    except KeyError:
        return None
    if isinstance(v, (str, int, float)) and not isinstance(v, bool):
        return str(v)
    return None


def extract(read, p: EngineParams, stats: dict, with_seq: bool):
    stats["total"] += 1
    if read.is_unmapped:
        stats["unmapped"] += 1
        return None
    if read.is_secondary or read.is_supplementary:
        stats["secondary_or_supplementary"] += 1
        return None
    if read.is_duplicate:
        stats["duplicate"] += 1
        return None
    mapq = read.mapping_quality
    if mapq < p.min_mapq:
        stats["low_mapq"] += 1
        return None
    barcode = _tag(read, p.barcode_tag)
    if barcode is None:
        stats["missing_barcode"] += 1
        return None
    if p.barcodes is not None and barcode not in p.barcodes:
        stats["other_barcode"] += 1
        return None
    if read.reference_start < 0:
        stats["unmapped"] += 1
        return None
    cigar = read.cigartuples or []
    span = reference_span(cigar)
    if span <= 0:
        stats["no_alignment_span"] += 1
        return None
    stats["kept"] += 1
    hap = _tag(read, p.hp_tag) if p.hp_tag else None
    if hap not in ("1", "2"):
        hap = UNTAGGED
    ps = (_tag(read, p.ps_tag) if p.ps_tag else None) or MISSING
    start = read.reference_start
    return (barcode, hap), Read(
        name=read.query_name or "*",
        start=start,
        end=start + span,
        reverse=read.is_reverse,
        mapq=mapq,
        ps=ps,
        cigar=cigar,
        seq=(read.query_sequence or "") if with_seq else "",
    )


def _part(parts_dir: Path, idx: int, kind: str) -> Path:
    return parts_dir / f"{idx:05d}.{kind}.tsv"


def process_chrom(bam_path, chrom, idx, chrom_len, p: EngineParams, parts_dir):
    """One chromosome; returns (stats, groups, kmers, skipped)."""
    parts_dir = Path(parts_dir)
    stats = dict.fromkeys(STAT_KEYS, 0)
    groups: dict = {}
    kmers: dict = defaultdict(int)
    skipped = []
    pairs_on = p.write_overlaps or p.write_junctions

    def summary(g):
        s = groups.get(g)
        if s is None:
            s = groups[g] = dict.fromkeys(
                ["n_reads", "n_overlaps", "n_end_to_start", "n_contained", "n_junctions"], 0
            )
        return s

    with pysam.AlignmentFile(bam_path, "rb", threads=p.hts_threads) as bam:
        over = set()
        if pairs_on and p.max_reads_per_group > 0:
            counts: dict = defaultdict(int)
            scratch = dict.fromkeys(STAT_KEYS, 0)
            for rec in bam.fetch(contig=chrom):
                got = extract(rec, p, scratch, False)
                if got:
                    counts[got[0]] += 1
            for g, n in counts.items():
                if n > p.max_reads_per_group:
                    skipped.append((g[0], g[1], n))
                    over.add(g)
            skipped.sort()

        w_ov = open(_part(parts_dir, idx, "overlaps"), "w") if p.write_overlaps else None
        w_jn = open(_part(parts_dir, idx, "tn5_junctions"), "w") if p.write_junctions else None
        active: dict = defaultdict(list)
        intervals: dict = defaultdict(list)
        bins: dict = defaultdict(lambda: defaultdict(int))
        last_start = None
        try:
            for rec in bam.fetch(contig=chrom):
                got = extract(rec, p, stats, pairs_on)
                if got is None:
                    continue
                g, read = got
                if last_start is not None and read.start < last_start:
                    raise ValueError(f"{chrom}: BAM is not coordinate-sorted (samtools sort first)")
                last_start = read.start
                summary(g)["n_reads"] += 1
                if p.write_density:
                    intervals[g].append((read.start, read.end))
                    bins[g][read.start // p.bin_size] += 1
                if not pairs_on or g in over:
                    continue
                lst = [a for a in active[g] if a.end > read.start]
                for a in lst:
                    _pair(p, chrom, g, a, read, summary(g), kmers, w_ov, w_jn)
                lst.append(read)
                active[g] = lst
        finally:
            for w in (w_ov, w_jn):
                if w:
                    w.close()

    if p.write_density:
        with open(_part(parts_dir, idx, "density_segments"), "w") as w:
            for g in sorted(intervals):
                for s, e, d in depth_segments(intervals[g]):
                    w.write(f"{g[0]}\t{g[1]}\t{chrom}\t{s}\t{e}\t{d}\n")
        with open(_part(parts_dir, idx, "density_bins"), "w") as w:
            for g in sorted(bins):
                for b in sorted(bins[g]):
                    s = b * p.bin_size
                    e = min((b + 1) * p.bin_size, max(chrom_len, s + 1))
                    w.write(f"{g[0]}\t{g[1]}\t{chrom}\t{s}\t{e}\t{bins[g][b]}\n")
    return stats, groups, dict(kmers), skipped


def _pair(p, chrom, g, a, b, gs, kmers, w_ov, w_jn):
    r1, r2 = (a, b) if a.key() <= b.key() else (b, a)
    os_, oe, contained = overlap_window(r1, r2)
    bp = oe - os_
    jsize = p.is_junction_size(bp)
    need = w_ov is not None or (jsize and w_jn is not None)
    detail = pair_detail(r1, r2, os_, oe) if need else None
    is_junction = jsize and detail is not None and len(detail[2]) <= p.max_mismatches
    gs["n_overlaps"] += 1
    gs["n_contained" if contained else "n_end_to_start"] += 1
    if is_junction:
        gs["n_junctions"] += 1
    if detail is None:
        return
    line = format_pair(p, g, chrom, r1, r2, os_, oe, contained, *detail)
    if w_ov is not None:
        w_ov.write(line)
    if is_junction:
        if w_jn is not None:
            w_jn.write(line)
        kmers[(g[0], g[1], detail[0])] += 1


def _join(parts_dir: Path, n: int, kind: str, columns, dest: Path):
    with gzip.open(dest, "wt", compresslevel=6) as out:
        out.write("\t".join(columns) + "\n")
        for i in range(n):
            part = _part(parts_dir, i, kind)
            if part.exists():
                with open(part) as f:
                    shutil.copyfileobj(f, out, 1 << 20)
                part.unlink()


def run(bam_path, out_dir, prefix, chroms, p: EngineParams, threads: int = 1):
    """Same contract as rrecon._core.run: returns (stats, skipped, files)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    parts_dir = out_dir / f".{prefix}.parts.{os.getpid()}"
    parts_dir.mkdir(parents=True, exist_ok=True)
    try:
        jobs = [(str(bam_path), c, i, n, p, str(parts_dir)) for i, (c, n) in enumerate(chroms)]
        if threads > 1 and len(jobs) > 1:
            with ProcessPoolExecutor(max_workers=min(threads, len(jobs))) as pool:
                results = list(pool.map(process_chrom, *zip(*jobs)))
        else:
            results = [process_chrom(*j) for j in jobs]

        stats = dict.fromkeys(STAT_KEYS, 0)
        groups: dict = {}
        kmers: dict = defaultdict(int)
        skipped = []
        for (chrom, _), (st, gr, km, sk) in zip(chroms, results):
            for k in STAT_KEYS:
                stats[k] += st[k]
            for g, s in gr.items():
                acc = groups.setdefault(g, dict.fromkeys(s, 0))
                for k, v in s.items():
                    acc[k] += v
            for k, v in km.items():
                kmers[k] += v
            skipped += [(b, h, chrom, n) for b, h, n in sk]

        files = []
        n = len(chroms)

        def emit(kind, cols):
            dest = out_dir / f"{prefix}.{kind}.tsv.gz"
            _join(parts_dir, n, kind, cols, dest)
            files.append((kind, str(dest)))

        if p.write_density:
            emit("density_segments", SEGMENT_COLUMNS)
            emit("density_bins", BIN_COLUMNS)
        if p.write_overlaps:
            emit("overlaps", OVERLAP_COLUMNS)
        if p.write_junctions:
            emit("tn5_junctions", OVERLAP_COLUMNS)
            dest = out_dir / f"{prefix}.tn5_summary.tsv.gz"
            with gzip.open(dest, "wt", compresslevel=6) as w:
                w.write("\t".join(SUMMARY_COLUMNS) + "\n")
                for g in sorted(groups):
                    s = groups[g]
                    w.write(
                        f"{g[0]}\t{g[1]}\t{s['n_reads']}\t{s['n_overlaps']}\t"
                        f"{s['n_end_to_start']}\t{s['n_contained']}\t{s['n_junctions']}\t"
                        f"{format_pct(s['n_junctions'], s['n_overlaps'])}\n"
                    )
            files.append(("tn5_summary", str(dest)))
            dest = out_dir / f"{prefix}.tn5_kmers.tsv.gz"
            with gzip.open(dest, "wt", compresslevel=6) as w:
                w.write("\t".join(KMER_COLUMNS) + "\n")
                for k in sorted(kmers):
                    w.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{kmers[k]}\n")
            files.append(("tn5_kmers", str(dest)))
        return [(k, stats[k]) for k in STAT_KEYS], skipped, files
    finally:
        shutil.rmtree(parts_dir, ignore_errors=True)
