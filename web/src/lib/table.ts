// Columnar tables read in the browser from rrecon outputs (.tsv.gz or .parquet).
// Numbers become Float64Array, booleans Uint8Array, everything else string[].

import { parquetMetadata, parquetRead } from "hyparquet";

export type Column = string[] | Float64Array | Uint8Array;
export interface Table {
  n: number;
  cols: Record<string, Column>;
}

export const KINDS = [
  "density_segments",
  "density_bins",
  "overlaps",
  "tn5_junctions",
  "tn5_summary",
  "tn5_kmers",
] as const;
export type Kind = (typeof KINDS)[number];

// keep in step with rrecon/schema.py
const NUM = new Set([
  "read1_start", "read1_end", "read1_mapq", "read2_start", "read2_end", "read2_mapq",
  "overlap_start", "overlap_end", "overlap_bp", "n_mismatches", "start", "end", "depth",
  "bin_start", "bin_end", "n_reads", "n_overlaps", "n_end_to_start", "n_contained",
  "n_junctions", "count", "pct_junctions",
]);
const BOOL = new Set(["has_indel_in_overlap", "seq_match"]);

export const str = (t: Table, c: string) => t.cols[c] as string[];
export const num = (t: Table, c: string) => t.cols[c] as Float64Array;
export const bool = (t: Table, c: string) => t.cols[c] as Uint8Array;

/** Which output a file holds, from its name (`<prefix>.<kind>.<ext>`). */
export function kindOf(name: string): Kind | "meta" | null {
  if (name.endsWith(".meta.json")) return "meta";
  for (const k of KINDS) if (name.includes(`.${k}.`)) return k;
  return null;
}
export function prefixOf(name: string): string {
  const k = kindOf(name);
  if (!k) return name;
  return name.slice(0, name.indexOf(k === "meta" ? ".meta.json" : `.${k}.`));
}

function makeColumn(name: string, n: number): Column {
  if (NUM.has(name)) return new Float64Array(n);
  if (BOOL.has(name)) return new Uint8Array(n);
  return new Array<string>(n);
}

/** Parse tab-separated text (header line first) into a table. */
export function parseTsv(text: string): Table {
  const lines = text.split("\n");
  while (lines.length && lines[lines.length - 1] === "") lines.pop();
  const header = (lines[0] ?? "").split("\t");
  const n = Math.max(0, lines.length - 1);
  const cols: Record<string, Column> = {};
  for (const h of header) cols[h] = makeColumn(h, n);
  const setters = header.map((h) => {
    const c = cols[h];
    if (c instanceof Float64Array) return (i: number, v: string) => { c[i] = +v; };
    if (c instanceof Uint8Array) return (i: number, v: string) => { c[i] = v === "true" ? 1 : 0; };
    const s = c as string[];
    return (i: number, v: string) => { s[i] = v; };
  });
  for (let i = 0; i < n; i++) {
    const f = lines[i + 1].split("\t");
    for (let k = 0; k < setters.length; k++) setters[k](i, f[k] ?? "");
  }
  return { n, cols };
}

async function gunzipText(buf: ArrayBuffer): Promise<string> {
  const stream = new Blob([buf]).stream().pipeThrough(new DecompressionStream("gzip"));
  return await new Response(stream).text();
}

export async function readParquet(buf: ArrayBuffer): Promise<Table> {
  const meta = parquetMetadata(buf);
  const n = Number(meta.num_rows);
  const names = meta.schema.slice(1).map((s) => s.name);
  const cols: Record<string, Column> = {};
  for (const c of names) cols[c] = makeColumn(c, n);
  await parquetRead({
    file: buf,
    metadata: meta,
    onChunk: ({ columnName, columnData, rowStart, rowEnd }) => {
      const col = cols[columnName];
      if (!col) return;
      const lo = Math.max(0, rowStart), hi = Math.min(n, rowEnd);
      for (let i = lo; i < hi; i++) {
        const v = columnData[i - rowStart];
        if (col instanceof Float64Array) col[i] = typeof v === "bigint" ? Number(v) : (v as number);
        else if (col instanceof Uint8Array) col[i] = v ? 1 : 0;
        else (col as string[])[i] = v == null ? "" : String(v);
      }
    },
  });
  return { n, cols };
}

export async function readTable(name: string, buf: ArrayBuffer): Promise<Table> {
  if (name.endsWith(".parquet")) return readParquet(buf);
  const text = name.endsWith(".gz") ? await gunzipText(buf) : new TextDecoder().decode(buf);
  return parseTsv(text);
}

/** Row indices of `t` for which `keep` is true. */
export function where(t: Table, keep: (i: number) => boolean): Int32Array {
  const out = new Int32Array(t.n);
  let k = 0;
  for (let i = 0; i < t.n; i++) if (keep(i)) out[k++] = i;
  return out.subarray(0, k);
}

/** Rows of `t` as tab- or comma-separated text (with header), for downloads. */
export function toDelimited(t: Table, rows: ArrayLike<number>, sep = ","): string {
  const names = Object.keys(t.cols);
  const q = (v: string) => (sep === "," && /[",\n]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v);
  const out = [names.join(sep)];
  for (let r = 0; r < rows.length; r++) {
    const i = rows[r];
    out.push(names.map((c) => {
      const col = t.cols[c];
      if (col instanceof Uint8Array) return col[i] ? "true" : "false";
      return q(String(col[i]));
    }).join(sep));
  }
  return out.join("\n") + "\n";
}
