"""Behaviour of rrecon.run on the synthetic edge-case BAM, for both engines,
and Rust == Python parity on every output file."""

import gzip
import json

import pandas as pd
import pytest

import rrecon
from rrecon.schema import COLUMNS

ENGINES = ["python"] + (["rust"] if rrecon.rust_available() else [])
HAP = dict(hp_tag="HP", ps_tag="PS")


def run(bam, tmp_path, engine, outputs=("density", "overlaps", "junctions"), **kw):
    return rrecon.run(bam, tmp_path / engine, outputs, engine=engine, **{**HAP, **kw})


def pair(df, a, b):
    m = df[(df.read1_name == a) & (df.read2_name == b)]
    assert len(m) == 1, f"expected one pair {a}-{b}, got {len(m)}"
    return m.iloc[0]


def names(df):
    return {(r.read1_name, r.read2_name) for r in df.itertuples()}


@pytest.fixture(scope="module", params=ENGINES)
def result(request, bam, tmp_path_factory):
    return run(bam, tmp_path_factory.mktemp("out"), request.param)


def test_files_and_columns(result):
    for kind, cols in COLUMNS.items():
        for fmt in ("tsv", "parquet"):
            df = rrecon.read_table(result.path(kind, fmt))
            assert list(df.columns) == cols, (kind, fmt)
    meta = json.loads(result.meta_path.read_text())
    assert [c["name"] for c in meta["chroms"]] == ["chr1", "chr2"]  # chrUn skipped
    assert meta["engine"] == result.engine


def test_contained_read_overlap_is_inner_read(result):
    r = pair(result.table("overlaps"), "contained_outer", "contained_inner")
    assert (r.overlap_start, r.overlap_end, r.overlap_bp) == (1050, 1100, 50)
    assert r.overlap_type == "contained"
    assert r.read2_strand == "-"
    assert r.seq1 == r.seq2 and len(r.seq1) == 50


def test_exact_9bp_junction(result):
    r = pair(result.table("tn5_junctions"), "j9_left", "j9_right")
    assert (r.overlap_start, r.overlap_end, r.overlap_bp) == (2091, 2100, 9)
    assert r.overlap_type == "end_to_start"
    assert r.n_mismatches == 0 and r.mismatch_positions == "."
    assert r.seq_match and not r.has_indel_in_overlap
    assert r.phase_block == "900"


def test_one_mismatch_is_a_junction(result):
    r = pair(result.table("tn5_junctions"), "mm1_left", "mm1_right")
    assert r.n_mismatches == 1 and r.mismatch_positions == "4"
    assert r.seq1[4] != r.seq2[4]


def test_two_mismatches_overlap_but_not_junction(result):
    r = pair(result.table("overlaps"), "mm2_left", "mm2_right")
    assert r.n_mismatches == 2 and r.mismatch_positions == "1,6"
    assert not r.seq_match
    assert ("mm2_left", "mm2_right") not in names(result.table("tn5_junctions"))


def test_deletion_in_overlap(result):
    r = pair(result.table("overlaps"), "del_left", "del_right")
    assert (r.overlap_start, r.overlap_end) == (5080, 5100)
    assert r.has_indel_in_overlap
    assert r.seq2[10:12] == "--" and "-" not in r.seq1
    assert r.mismatch_positions == "10,11"


def test_insertion_in_overlap(result):
    r = pair(result.table("overlaps"), "ins_left", "ins_right")
    assert r.has_indel_in_overlap
    assert r.n_mismatches == 0  # inserted bases are not on the reference window
    assert not r.seq_match


def test_touching_reads_do_not_overlap(result):
    df = result.table("overlaps")
    assert not ({"touch_left", "touch_right"} & (set(df.read1_name) | set(df.read2_name)))


def test_haplotypes_kept_apart(result):
    df = result.table("overlaps")
    hp_reads = {"hp1_read", "hp2_read"}
    assert not (hp_reads & (set(df.read1_name) | set(df.read2_name)))
    r = pair(df, "untagged_a", "untagged_b")
    assert r.haplotype == "untagged" and r.phase_block == "."
    j = result.table("tn5_junctions")
    assert set(j.haplotype) == {"1", "2", "untagged"}


def test_filtered_reads_absent(result):
    df = result.table("overlaps")
    seen = set(df.read1_name) | set(df.read2_name)
    for n in ("dup_left", "dup_right", "supp_read", "secondary_read", "nobc_left", "nobc_right"):
        assert n not in seen
    s = result.stats
    assert s["duplicate"] == 2
    assert s["secondary_or_supplementary"] == 2
    assert s["missing_barcode"] == 2


def test_soft_clips_skipped(result):
    r = pair(result.table("tn5_junctions"), "clip_left", "clip_right")
    assert r.overlap_bp == 9 and r.n_mismatches == 0


def test_expected_junction_set(result):
    assert names(result.table("tn5_junctions")) == {
        ("j9_left", "j9_right"),
        ("mm1_left", "mm1_right"),
        ("untagged_a", "untagged_b"),
        ("clip_left", "clip_right"),
        ("lowq_left", "lowq_right"),
        ("cell2_left", "cell2_right"),
        ("chr2_left", "chr2_right"),
    }


def test_summary_and_kmers(result):
    s = result.table("tn5_summary").set_index(["barcode", "haplotype"])
    hp1 = s.loc[("AAAC", "1")]
    assert hp1.n_junctions == 5
    assert hp1.n_overlaps == hp1.n_contained + hp1.n_end_to_start
    assert hp1.pct_junctions == pytest.approx(100 * hp1.n_junctions / hp1.n_overlaps, abs=1e-4)
    assert s.loc[("GGGT", "2")].n_junctions == 1
    k = result.table("tn5_kmers")
    assert k["count"].sum() == 7
    assert (k.junction_seq.str.len() == 9).all()


def test_density(result):
    seg = result.table("density_segments")
    a = seg[(seg.barcode == "AAAC") & (seg.haplotype == "1") & (seg.chrom == "chr1")]
    # contained pair: depth 1 over [1000,1050), 2 over [1050,1100), 1 over [1100,1200)
    assert a[(a.start >= 1000) & (a.end <= 1200)][["start", "end", "depth"]].values.tolist() == [
        [1000, 1050, 1],
        [1050, 1100, 2],
        [1100, 1200, 1],
    ]
    bins = result.table("density_bins")
    tot = bins.groupby(["barcode", "haplotype"])["n_reads"].sum()
    summ = result.table("tn5_summary").set_index(["barcode", "haplotype"])["n_reads"]
    assert tot.sort_index().tolist() == summ.sort_index().tolist()


def test_min_mapq_and_allow_list(bam, tmp_path):
    for engine in ENGINES:
        r = run(bam, tmp_path / "q", engine, ("junctions",), min_mapq=30, barcodes=["GGGT"])
        j = r.table("tn5_junctions")
        assert names(j) == {("cell2_left", "cell2_right")}


def test_no_hp_tag_means_all_untagged(bam, tmp_path):
    for engine in ENGINES:
        r = rrecon.run(bam, tmp_path / engine, ("overlaps",), engine=engine)
        df = r.table("overlaps")
        assert set(df.haplotype) == {"untagged"}
        # without haplotypes, the HP1/HP2 reads now pair up
        assert ("hp1_read", "hp2_read") in names(df)


def test_max_reads_per_group(bam, tmp_path):
    for engine in ENGINES:
        r = run(bam, tmp_path, engine, ("junctions",), max_reads_per_group=3)
        skipped = {(g["barcode"], g["haplotype"], g["chrom"]) for g in r.skipped_groups}
        assert ("AAAC", "1", "chr1") in skipped
        assert "j9_left" not in set(r.table("tn5_junctions").read1_name)


def test_all_contigs(bam, tmp_path):
    r = run(bam, tmp_path / "all", ENGINES[-1], ("overlaps",), all_contigs=True)
    assert "chrUn_test" in set(r.table("overlaps").chrom)


# --------------------------------------------------------------------------- parity


def _text(path):
    with gzip.open(path, "rt") as f:
        return f.read()


@pytest.mark.skipif(not rrecon.rust_available(), reason="Rust extension not built")
@pytest.mark.parametrize(
    "kw",
    [
        {},
        {"threads": 3},
        {"min_mapq": 30, "max_mismatches": 0, "overlap_tolerance": 2, "bin_size": 1000},
        {"max_reads_per_group": 3, "max_seq_bp": 10},
        {"hp_tag": None, "ps_tag": None, "all_contigs": True},
    ],
)
def test_rust_python_identical(bam, tmp_path, kw):
    py = run(bam, tmp_path, "python", **kw)
    rs = run(bam, tmp_path, "rust", **kw)
    assert py.stats == rs.stats
    assert py.skipped_groups == rs.skipped_groups
    assert py.files.keys() == rs.files.keys()
    for kind in py.files:
        assert _text(py.path(kind, "tsv")) == _text(rs.path(kind, "tsv")), kind
        a = pd.read_parquet(py.path(kind, "parquet"))
        b = pd.read_parquet(rs.path(kind, "parquet"))
        pd.testing.assert_frame_equal(a, b)


@pytest.mark.skipif(not rrecon.rust_available(), reason="Rust extension not built")
def test_parity_on_random_bam(tmp_path):
    from rrecon.synthetic import random_reads, write_bam

    chroms = {"chr1": 300_000, "chr2": 200_000, "chrX": 100_000}
    bam = write_bam(tmp_path / "rand.bam", random_reads(chroms, n_cells=4, reads_per_cell=400,
                                                        read_len=(50, 2000)), chroms)
    py = run(bam, tmp_path, "python", formats=("tsv",))
    rs = run(bam, tmp_path, "rust", formats=("tsv",), threads=4, hts_threads=2)
    for kind in py.files:
        assert _text(py.path(kind, "tsv")) == _text(rs.path(kind, "tsv")), kind
    assert len(py.table("overlaps")) > 100
    assert len(py.table("tn5_junctions")) > 10
