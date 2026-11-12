// One loaded rrecon run: its tables, genome layout and per-cell totals.

import { makeGenome, type Genome } from "./genome";
import { HAPS, hapIndex } from "./stats";
import { kindOf, num, prefixOf, readTable, str, type Kind, type Table } from "./table";

export interface RunMeta {
  tool?: string;
  version?: string;
  engine?: string;
  bam?: string;
  prefix?: string;
  created?: string;
  seconds?: number;
  outputs?: string[];
  params?: Record<string, unknown>;
  chroms?: { name: string; length: number }[];
  stats?: Record<string, number>;
  skipped_groups?: { barcode: string; haplotype: string; chrom: string; n_reads: number }[];
}

export interface CellInfo {
  barcode: string;
  /** per haplotype (1, 2, untagged) */
  reads: number[];
  overlaps: number[];
  junctions: number[];
}

export interface Dataset {
  name: string;
  meta: RunMeta | null;
  genome: Genome;
  tables: Partial<Record<Kind, Table>>;
  files: string[];
  cells: CellInfo[];
  junctionSize: number;
  tolerance: number;
  maxMismatches: number;
}

export interface InFile {
  name: string;
  data: ArrayBuffer;
}

/** Pick one run's files out of what was dropped (the prefix with most outputs). */
export function chooseRun(names: string[]): { prefix: string; picked: Map<Kind | "meta", string> } {
  const byPrefix = new Map<string, Map<Kind | "meta", string>>();
  for (const n of names) {
    const k = kindOf(n);
    if (!k) continue;
    const p = prefixOf(n);
    if (!byPrefix.has(p)) byPrefix.set(p, new Map());
    const m = byPrefix.get(p)!;
    // Parquet over TSV when both are there
    if (!m.has(k) || n.endsWith(".parquet")) m.set(k, n);
  }
  let best = "", bestMap = new Map<Kind | "meta", string>();
  for (const [p, m] of byPrefix) if (m.size > bestMap.size) { best = p; bestMap = m; }
  return { prefix: best, picked: bestMap };
}

export async function loadDataset(files: InFile[]): Promise<Dataset> {
  const { prefix, picked } = chooseRun(files.map((f) => f.name));
  if (!picked.size) {
    throw new Error(
      "No rrecon outputs found. Expected files like <prefix>.overlaps.parquet, " +
        "<prefix>.tn5_junctions.tsv.gz or <prefix>.meta.json.",
    );
  }
  const byName = new Map(files.map((f) => [f.name, f]));
  let meta: RunMeta | null = null;
  const tables: Partial<Record<Kind, Table>> = {};
  for (const [kind, name] of picked) {
    const f = byName.get(name)!;
    try {
      if (kind === "meta") meta = JSON.parse(new TextDecoder().decode(f.data));
      else tables[kind] = await readTable(name, f.data);
    } catch (e) {
      throw new Error(`Could not read ${name}: ${e instanceof Error ? e.message : e}`, { cause: e });
    }
  }
  return assemble(prefix, meta, tables, [...picked.values()]);
}

function inferLengths(tables: Partial<Record<Kind, Table>>): [string, number][] {
  const max = new Map<string, number>();
  const add = (t: Table | undefined, endCol: string) => {
    if (!t || !t.cols[endCol]) return;
    const c = str(t, "chrom"), e = num(t, endCol);
    for (let i = 0; i < t.n; i++) max.set(c[i], Math.max(max.get(c[i]) ?? 0, e[i]));
  };
  add(tables.density_bins, "bin_end");
  add(tables.density_segments, "end");
  add(tables.overlaps, "read2_end");
  add(tables.tn5_junctions, "read2_end");
  return [...max];
}

export function assemble(
  name: string,
  meta: RunMeta | null,
  tables: Partial<Record<Kind, Table>>,
  files: string[] = [],
): Dataset {
  const lengths: [string, number][] = meta?.chroms?.length
    ? meta.chroms.map((c) => [c.name, c.length])
    : inferLengths(tables);
  const p = meta?.params ?? {};
  const n = (k: string, d: number) => (typeof p[k] === "number" ? (p[k] as number) : d);
  return {
    name: meta?.prefix ?? name,
    meta,
    genome: makeGenome(lengths),
    tables,
    files,
    cells: cellTotals(tables),
    junctionSize: n("junction_size", 9),
    tolerance: n("overlap_tolerance", 0),
    maxMismatches: n("max_mismatches", 1),
  };
}

/** Reads, overlaps and junctions per cell and haplotype, from whichever tables exist. */
export function cellTotals(tables: Partial<Record<Kind, Table>>): CellInfo[] {
  const cells = new Map<string, CellInfo>();
  const get = (b: string) => {
    let c = cells.get(b);
    if (!c) cells.set(b, (c = { barcode: b, reads: [0, 0, 0], overlaps: [0, 0, 0], junctions: [0, 0, 0] }));
    return c;
  };
  const s = tables.tn5_summary;
  if (s) {
    const b = str(s, "barcode"), h = str(s, "haplotype");
    const r = num(s, "n_reads"), o = num(s, "n_overlaps"), j = num(s, "n_junctions");
    for (let i = 0; i < s.n; i++) {
      const c = get(b[i]), k = hapIndex(h[i]);
      c.reads[k] += r[i]; c.overlaps[k] += o[i]; c.junctions[k] += j[i];
    }
  } else {
    const bins = tables.density_bins;
    if (bins) {
      const b = str(bins, "barcode"), h = str(bins, "haplotype"), r = num(bins, "n_reads");
      for (let i = 0; i < bins.n; i++) get(b[i]).reads[hapIndex(h[i])] += r[i];
    }
    const count = (t: Table | undefined, key: "overlaps" | "junctions") => {
      if (!t) return;
      const b = str(t, "barcode"), h = str(t, "haplotype");
      for (let i = 0; i < t.n; i++) get(b[i])[key][hapIndex(h[i])] += 1;
    };
    count(tables.overlaps, "overlaps");
    count(tables.tn5_junctions, "junctions");
  }
  return [...cells.values()].sort((a, b) => sum(b.junctions) - sum(a.junctions) || sum(b.reads) - sum(a.reads) || a.barcode.localeCompare(b.barcode));
}

export const sum = (v: number[]) => v.reduce((a, b) => a + b, 0);
export { HAPS };
