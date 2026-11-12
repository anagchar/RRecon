import { describe, expect, it } from "vitest";
import { assemble, chooseRun } from "./dataset";
import { formatRegion, makeGenome, parseRegion } from "./genome";
import { baseCounts, defaultRowFilter, filterRows, logoHeights, mismatchStats, selectRows, sizeHistogram, sortRows } from "./stats";
import { kindOf, parseTsv, prefixOf, toDelimited, where } from "./table";

const OVERLAPS = [
  "barcode\thaplotype\tphase_block\tchrom\tread1_name\tread1_start\tread1_end\tread1_strand\tread1_mapq\tread2_name\tread2_start\tread2_end\tread2_strand\tread2_mapq\toverlap_start\toverlap_end\toverlap_bp\toverlap_type\tseq1\tseq2\tn_mismatches\tmismatch_positions\thas_indel_in_overlap\tseq_match",
  "AAAC\t1\t900\tchr1\ta\t100\t200\t+\t60\tb\t191\t300\t-\t60\t191\t200\t9\tend_to_start\tACGTACGTA\tACGTACGTA\t0\t.\tfalse\ttrue",
  "AAAC\t1\t900\tchr1\tc\t1000\t1200\t+\t60\td\t1050\t1100\t+\t60\t1050\t1100\t50\tcontained\t.\t.\t3\t.\tfalse\tfalse",
  "AAAC\t2\t.\tchr2\te\t10\t110\t+\t60\tf\t101\t200\t+\t60\t101\t110\t9\tend_to_start\tAAAAAAAAA\tAAAACAAAA\t1\t4\tfalse\ttrue",
  "GGGT\tuntagged\t.\tchr1\tg\t500\t600\t+\t60\th\t591\t700\t+\t60\t591\t600\t9\tend_to_start\tCCCCCCCCC\tCCGCCCCCA\t2\t2,8\ttrue\tfalse",
].join("\n");

const genome = makeGenome([["chr2", 5000], ["chr1", 20000], ["chrX", 1000]]);
const t = parseTsv(OVERLAPS + "\n");
const all = { barcode: null, hap: "all" as const, region: { chrom: null, start: 0, end: genome.length } };
const pos: [string, string] = ["overlap_start", "overlap_end"];

describe("table", () => {
  it("parses typed columns", () => {
    expect(t.n).toBe(4);
    expect(t.cols.overlap_bp).toBeInstanceOf(Float64Array);
    expect(Array.from(t.cols.has_indel_in_overlap as Uint8Array)).toEqual([0, 0, 0, 1]);
    expect((t.cols.phase_block as string[])[0]).toBe("900");
  });
  it("names files", () => {
    expect(kindOf("s1.tn5_junctions.parquet")).toBe("tn5_junctions");
    expect(kindOf("s1.meta.json")).toBe("meta");
    expect(kindOf("notes.txt")).toBeNull();
    expect(prefixOf("my.sample.overlaps.tsv.gz")).toBe("my.sample");
  });
  it("writes CSV with quoting", () => {
    const csv = toDelimited(t, [3], ",");
    expect(csv.split("\n")[1]).toContain('"2,8"');
    expect(csv.split("\n")[1]).toContain(",true,false");
  });
  it("chooses one run and prefers parquet", () => {
    const { prefix, picked } = chooseRun(["a.overlaps.tsv.gz", "a.overlaps.parquet", "a.meta.json", "b.overlaps.parquet"]);
    expect(prefix).toBe("a");
    expect(picked.get("overlaps")).toBe("a.overlaps.parquet");
  });
});

describe("genome", () => {
  it("orders chromosomes naturally", () => {
    expect(genome.chroms.map((c) => c.name)).toEqual(["chr1", "chr2", "chrX"]);
    expect(genome.byName.get("chr2")!.offset).toBe(20000);
  });
  it("parses regions", () => {
    expect(parseRegion(genome, "chr1:1,001-2,000")).toEqual({ chrom: "chr1", start: 1000, end: 2000 });
    expect(parseRegion(genome, "2")).toEqual({ chrom: "chr2", start: 0, end: 5000 });
    expect(parseRegion(genome, "chr1:1k-2k")).toEqual({ chrom: "chr1", start: 999, end: 2000 });
    expect(parseRegion(genome, "chr9")).toBeNull();
    expect(parseRegion(genome, "genome")!.chrom).toBeNull();
    expect(formatRegion({ chrom: "chr1", start: 999, end: 2000 })).toBe("chr1:1,000-2,000");
  });
});

describe("stats", () => {
  it("selects by cell, haplotype and region", () => {
    expect(selectRows(t, genome, all, pos).length).toBe(4);
    expect(Array.from(selectRows(t, genome, { ...all, barcode: "AAAC", hap: "1" }, pos))).toEqual([0, 1]);
    expect(Array.from(selectRows(t, genome, { ...all, region: { chrom: "chr1", start: 0, end: 600 } }, pos))).toEqual([0, 3]);
  });
  it("histograms sizes per haplotype", () => {
    const h = sizeHistogram(t, where(t, () => true), 30);
    expect(h[0][9]).toBe(1);
    expect(h[0][30]).toBe(1); // 50 bp lands in the last bin
    expect(h[1][9]).toBe(1);
    expect(h[2][9]).toBe(1);
  });
  it("counts mismatches and their positions", () => {
    const m = mismatchStats(t, Int32Array.from([0, 2, 3]), 3, 9);
    expect(m.counts[0][0]).toBe(1);
    expect(m.counts[1][1]).toBe(1);
    expect(m.positions[1][4]).toBe(1);
    expect(m.positions[2][2]).toBe(1);
    expect(m.positions[2][8]).toBe(1);
  });
  it("filters and sorts rows", () => {
    const rows = where(t, () => true);
    expect(Array.from(filterRows(t, rows, { ...defaultRowFilter(), junctionsOnly: true }, 9, 0, 1))).toEqual([0, 2]);
    expect(Array.from(filterRows(t, rows, { ...defaultRowFilter(), type: "contained" }, 9, 0, 1))).toEqual([1]);
    expect(Array.from(filterRows(t, rows, { ...defaultRowFilter(), indel: "with" }, 9, 0, 1))).toEqual([3]);
    expect(Array.from(filterRows(t, rows, { ...defaultRowFilter(), phase: "900" }, 9, 0, 1))).toEqual([0, 1]);
    expect(Array.from(sortRows(t, rows, "n_mismatches", true))).toEqual([1, 3, 2, 0]);
    expect(Array.from(sortRows(t, rows, "barcode", false))).toEqual([0, 1, 2, 3]);
  });
  it("builds logo heights", () => {
    const c = baseCounts(["AC", "AG", "AT", "xx"], [2, 1, 1, 5], 2);
    expect(c[0]).toEqual([4, 0, 0, 0]);
    const h = logoHeights(c);
    expect(h[0][0]).toBeCloseTo(2);
    expect(h[1].reduce((a, b) => a + b, 0)).toBeLessThan(1);
  });
});

describe("dataset", () => {
  it("assembles cells and the junction rule", () => {
    const ds = assemble("x", { prefix: "x", params: { junction_size: 9, max_mismatches: 1 }, chroms: [{ name: "chr1", length: 20000 }, { name: "chr2", length: 5000 }] }, { overlaps: t });
    expect(ds.cells.map((c) => c.barcode)).toEqual(["AAAC", "GGGT"]);
    expect(ds.cells[0].overlaps).toEqual([2, 1, 0]);
    expect(ds.genome.length).toBe(25000);
  });
});
