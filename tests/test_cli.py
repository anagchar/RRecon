import json

from click.testing import CliRunner
from fastapi.testclient import TestClient

import rrecon
from rrecon.cli import main
from rrecon.server import create_app


def test_all_command(bam, tmp_path):
    out = tmp_path / "res"
    r = CliRunner().invoke(
        main,
        ["all", str(bam), "-o", str(out), "--hp-tag", "HP", "--ps-tag", "PS",
         "--threads", "2", "--format", "parquet", "--prefix", "s1"],
    )
    assert r.exit_code == 0, r.output
    names = sorted(p.name for p in out.iterdir())
    assert names == sorted(
        [f"s1.{k}.parquet" for k in ("density_segments", "density_bins", "overlaps",
                                      "tn5_junctions", "tn5_summary", "tn5_kmers")]
        + ["s1.meta.json"]
    )
    meta = json.loads((out / "s1.meta.json").read_text())
    assert meta["params"]["hp_tag"] == "HP"
    assert meta["files"]["overlaps"] == {"parquet": "s1.overlaps.parquet"}


def test_subcommands_write_only_their_files(bam, tmp_path):
    for cmd, kinds in [("density", {"density_segments", "density_bins"}),
                       ("overlaps", {"overlaps"}),
                       ("junctions", {"tn5_junctions", "tn5_summary", "tn5_kmers"})]:
        out = tmp_path / cmd
        r = CliRunner().invoke(main, [cmd, str(bam), "-o", str(out), "--format", "tsv",
                                      "--engine", "python"])
        assert r.exit_code == 0, r.output
        got = {p.name.split(".")[1] for p in out.glob("*.tsv.gz")}
        assert got == kinds


def test_barcode_file(bam, tmp_path):
    f = tmp_path / "keep.txt"
    f.write_text("GGGT\n")
    out = tmp_path / "bc"
    r = CliRunner().invoke(main, ["junctions", str(bam), "-o", str(out), "--hp-tag", "HP",
                                  "--barcodes", str(f), "--format", "tsv"])
    assert r.exit_code == 0, r.output
    df = rrecon.read_table(out / "cases.tn5_junctions.tsv.gz")
    assert set(df.barcode) == {"GGGT"}


def test_bad_tag_is_a_clean_error(bam, tmp_path):
    r = CliRunner().invoke(main, ["all", str(bam), "-o", str(tmp_path), "--hp-tag", "HAP"])
    assert r.exit_code != 0
    assert "2-character" in r.output


def test_server_lists_and_serves_outputs(bam, tmp_path):
    rrecon.run(bam, tmp_path, ("junctions",), hp_tag="HP")
    (tmp_path / "unrelated.txt").write_text("x")
    c = TestClient(create_app(tmp_path))
    files = {f["name"] for f in c.get("/api/files").json()["files"]}
    assert "cases.tn5_junctions.parquet" in files and "cases.meta.json" in files
    assert "unrelated.txt" not in files
    assert c.get("/api/file/cases.meta.json").json()["tool"] == "rrecon"
    assert c.get("/api/file/unrelated.txt").status_code == 404
    assert c.get("/api/file/..%2Fetc%2Fpasswd").status_code == 404


def test_server_without_folder():
    c = TestClient(create_app(None))
    assert c.get("/api/files").json()["files"] == []
