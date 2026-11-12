"""Synthetic BAMs for tests, the benchmark and the GUI demo. No real data needed.

`SynthRead` describes one alignment against a random reference; `write_bam`
turns a list of them into a sorted, indexed BAM. `case_reads()` is the fixed
set of edge cases the tests check; `random_reads()` makes a large, realistic
mix (with Tn5 junctions planted) for benchmarking and the demo.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from pathlib import Path

import pysam

CIGAR = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}


def random_genome(chroms: dict[str, int], seed: int = 1) -> dict[str, str]:
    rng = random.Random(seed)
    return {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in chroms.items()}


@dataclass
class SynthRead:
    name: str
    chrom: str
    start: int
    cigar: list[tuple[str, int]]
    barcode: str | None = "AAAC"
    hp: int | None = None
    ps: int | None = None
    reverse: bool = False
    mapq: int = 60
    duplicate: bool = False
    supplementary: bool = False
    secondary: bool = False
    # 0-based offsets into the read sequence to replace with another base
    edits: list[int] = field(default_factory=list)

    def sequence(self, genome: dict[str, str]) -> str:
        ref = genome[self.chrom]
        out, rpos = [], self.start
        rng = random.Random(hash((self.name, self.start)) & 0xFFFF)
        for op, n in self.cigar:
            if op in "M=X":
                out.append(ref[rpos : rpos + n])
                rpos += n
            elif op in "IS":
                out.append("".join(rng.choice("ACGT") for _ in range(n)))
            elif op in "DN":
                rpos += n
        seq = list("".join(out))
        for k in self.edits:
            seq[k] = {"A": "C", "C": "G", "G": "T", "T": "A"}[seq[k]]
        return "".join(seq)


def write_bam(path, reads: list[SynthRead], chroms: dict[str, int], genome=None) -> Path:
    """Write a coordinate-sorted, indexed BAM; returns its path."""
    path = Path(path)
    genome = genome or random_genome(chroms)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": n} for c, n in chroms.items()],
    }
    tid = {c: i for i, c in enumerate(chroms)}
    unsorted = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for r in reads:
            a = pysam.AlignedSegment(out.header)
            a.query_name = r.name
            a.reference_id = tid[r.chrom]
            a.reference_start = r.start
            a.cigartuples = [(CIGAR[op], n) for op, n in r.cigar]
            a.query_sequence = r.sequence(genome)
            a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
            a.mapping_quality = r.mapq
            flag = 0
            if r.reverse:
                flag |= 0x10
            if r.secondary:
                flag |= 0x100
            if r.duplicate:
                flag |= 0x400
            if r.supplementary:
                flag |= 0x800
            a.flag = flag
            tags = []
            if r.barcode is not None:
                tags.append(("CB", r.barcode))
            if r.hp is not None:
                tags.append(("HP", r.hp))
            if r.ps is not None:
                tags.append(("PS", r.ps))
            a.set_tags(tags)
            out.write(a)
    pysam.sort("-o", str(path), str(unsorted))
    unsorted.unlink()
    pysam.index(str(path))
    return path


def M(n):
    return [("M", n)]


CASE_CHROMS = {"chr1": 20_000, "chr2": 5_000, "chrUn_test": 2_000}


def case_reads() -> list[SynthRead]:
    """The edge cases the tests check, all on cell AAAC unless stated."""
    R = SynthRead
    return [
        # 1 contained: B lies inside A -> overlap = B (50 bp), not A's tail
        R("contained_outer", "chr1", 1000, M(200), hp=1, ps=900),
        R("contained_inner", "chr1", 1050, M(50), hp=1, ps=900, reverse=True),
        # 2 exact 9 bp end-to-start junction
        R("j9_left", "chr1", 2000, M(100), hp=1, ps=900),
        R("j9_right", "chr1", 2091, M(100), hp=1, ps=900),
        # 3 one mismatch in the junction (offset 4 of the overlap)
        R("mm1_left", "chr1", 3000, M(100), hp=1, ps=900),
        R("mm1_right", "chr1", 3091, M(100), hp=1, ps=900, edits=[4]),
        # 4 two mismatches: still an overlap, not a junction (max 1)
        R("mm2_left", "chr1", 4000, M(100), hp=1, ps=900),
        R("mm2_right", "chr1", 4091, M(100), hp=1, ps=900, edits=[1, 6]),
        # 5 deletion inside the overlap [5080, 5100): ref 5090-5091 deleted
        R("del_left", "chr1", 5000, M(100), hp=1),
        R("del_right", "chr1", 5080, [("M", 10), ("D", 2), ("M", 50)], hp=1),
        # 6 insertion inside the overlap [6080, 6100)
        R("ins_left", "chr1", 6000, M(100), hp=1),
        R("ins_right", "chr1", 6080, [("M", 10), ("I", 3), ("M", 50)], hp=1),
        # 7 reads that only touch: no overlap
        R("touch_left", "chr1", 7000, M(100), hp=1),
        R("touch_right", "chr1", 7100, M(100), hp=1),
        # 8 haplotypes kept apart: HP1 and HP2 reads overlap but are never paired;
        #   the two untagged reads are paired with each other only
        R("hp1_read", "chr1", 8000, M(100), hp=1, ps=7000),
        R("hp2_read", "chr1", 8091, M(100), hp=2, ps=7000),
        R("untagged_a", "chr1", 8095, M(55)),
        R("untagged_b", "chr1", 8141, M(60)),
        # 9 duplicate and supplementary reads are dropped
        R("dup_left", "chr1", 9000, M(100), hp=1, duplicate=True),
        R("dup_right", "chr1", 9091, M(100), hp=1, duplicate=True),
        R("supp_read", "chr1", 9050, M(100), hp=1, supplementary=True),
        R("secondary_read", "chr1", 9060, M(100), hp=1, secondary=True),
        # 10 a missing barcode tag drops the read
        R("nobc_left", "chr1", 9500, M(100), barcode=None, hp=1),
        R("nobc_right", "chr1", 9591, M(100), barcode=None, hp=1),
        # 11 soft clips are skipped when reading the junction sequence
        R("clip_left", "chr1", 10000, M(100) + [("S", 7)], hp=1),
        R("clip_right", "chr1", 10091, [("S", 5)] + M(100), hp=1),
        # low MAPQ (kept by default, dropped with --min-mapq 30)
        R("lowq_left", "chr1", 12000, M(100), hp=1, mapq=5),
        R("lowq_right", "chr1", 12091, M(100), hp=1, mapq=5),
        # another cell, haplotype 2: a 9 bp junction
        R("cell2_left", "chr1", 11000, M(100), barcode="GGGT", hp=2, ps=10500),
        R("cell2_right", "chr1", 11091, M(100), barcode="GGGT", hp=2, ps=10500),
        # chr2: one junction, and a non-canonical contig that is skipped by default
        R("chr2_left", "chr2", 500, M(100), hp=1, ps=400),
        R("chr2_right", "chr2", 591, M(100), hp=1, ps=400),
        R("unplaced_a", "chrUn_test", 100, M(100)),
        R("unplaced_b", "chrUn_test", 191, M(100)),
    ]


def case_bam(path) -> Path:
    return write_bam(path, case_reads(), CASE_CHROMS)


def random_reads(
    chroms: dict[str, int],
    n_cells: int = 20,
    reads_per_cell: int = 2000,
    read_len: tuple[int, int] = (80, 3000),
    junction_frac: float = 0.15,
    seed: int = 7,
) -> list[SynthRead]:
    """A realistic mix: random reads per cell and haplotype, with a fraction
    followed by a planted neighbour that overlaps it by 9 bp (a Tn5 junction,
    sometimes with a mismatch) or by a random 1-60 bp."""
    rng = random.Random(seed)
    names = list(chroms)
    weights = [chroms[c] for c in names]
    cells = ["".join(rng.choice("ACGT") for _ in range(16)) + "-1" for _ in range(n_cells)]
    reads: list[SynthRead] = []
    k = 0
    for cell in cells:
        for _ in range(reads_per_cell):
            chrom = rng.choices(names, weights)[0]
            n = rng.randint(*read_len)
            start = rng.randrange(0, max(1, chroms[chrom] - 2 * n))
            hp = rng.choice([1, 2, None])
            ps = (start // 1_000_000) * 1_000_000 + 1 if hp else None
            mapq = rng.choice([60, 60, 60, 30, 0])
            k += 1
            reads.append(SynthRead(f"r{k}", chrom, start, M(n), cell, hp, ps,
                                   reverse=rng.random() < 0.5, mapq=mapq))
            if rng.random() < junction_frac:
                # a planted neighbour: mostly a 9 bp Tn5 junction (some with a
                # mismatch), otherwise a short chance overlap of 1-60 bp
                size = 9 if rng.random() < 0.7 else rng.randint(1, 60)
                s2 = start + n - size
                n2 = min(rng.randint(max(read_len[0], size + 1), read_len[1]), chroms[chrom] - s2)
                if n2 <= size:
                    continue
                k += 1
                edits = [rng.randrange(size)] if size == 9 and rng.random() < 0.2 else []
                reads.append(SynthRead(f"r{k}", chrom, s2, M(n2), cell, hp, ps,
                                       reverse=rng.random() < 0.5, mapq=mapq, edits=edits))
    return reads
