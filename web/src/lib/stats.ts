// Selections and summaries computed from the tables (pure, unit tested).

import { axisPos, type Genome, type Region } from "./genome";
import { bool, num, str, where, type Table } from "./table";

export const HAPS = ["1", "2", "untagged"] as const;
export type Hap = (typeof HAPS)[number];
export type HapSel = Hap | "all";
export const HAP_LABEL: Record<Hap, string> = { "1": "HP1", "2": "HP2", untagged: "untagged" };
export const HAP_COLOR: Record<Hap, string> = { "1": "--hp1", "2": "--hp2", untagged: "--hpu" };

export function hapIndex(h: string): number {
  return h === "1" ? 0 : h === "2" ? 1 : 2;
}
export function hapsShown(sel: HapSel): Hap[] {
  return sel === "all" ? [...HAPS] : [sel];
}

export interface Selection {
  barcode: string | null; // null = all cells
  hap: HapSel;
  region: Region;
}

/** Rows of a table in the selection; `pos` names the coordinate column(s) to test. */
export function selectRows(t: Table | undefined, g: Genome, sel: Selection, pos: [string, string]): Int32Array {
  if (!t) return new Int32Array(0);
  const b = str(t, "barcode"), h = str(t, "haplotype"), c = str(t, "chrom");
  const s = num(t, pos[0]), e = num(t, pos[1]);
  const r = sel.region;
  return where(t, (i) =>
    (sel.barcode == null || b[i] === sel.barcode) &&
    (sel.hap === "all" || h[i] === sel.hap) &&
    (r.chrom == null ? g.byName.has(c[i]) : c[i] === r.chrom && e[i] > r.start && s[i] < r.end),
  );
}

/** Histogram of overlap sizes, 1-bp bins from 1 to maxBp (the last bin collects larger), per haplotype. */
export function sizeHistogram(t: Table, rows: ArrayLike<number>, maxBp: number): number[][] {
  const out = HAPS.map(() => new Array<number>(maxBp + 1).fill(0));
  const bp = num(t, "overlap_bp"), h = str(t, "haplotype");
  for (let k = 0; k < rows.length; k++) {
    const i = rows[k];
    out[hapIndex(h[i])][Math.min(maxBp, Math.max(0, bp[i]))] += 1;
  }
  return out;
}

/** Counts of n_mismatches (0..max, last collects larger) and of mismatch offsets (0..size-1), per haplotype. */
export function mismatchStats(t: Table, rows: ArrayLike<number>, maxCount: number, size: number) {
  const counts = HAPS.map(() => new Array<number>(maxCount + 1).fill(0));
  const positions = HAPS.map(() => new Array<number>(size).fill(0));
  const nm = num(t, "n_mismatches"), mp = str(t, "mismatch_positions"), h = str(t, "haplotype");
  for (let k = 0; k < rows.length; k++) {
    const i = rows[k], hi = hapIndex(h[i]);
    counts[hi][Math.min(maxCount, nm[i])] += 1;
    if (mp[i] && mp[i] !== ".") for (const p of mp[i].split(",")) {
      const x = +p;
      if (x >= 0 && x < size) positions[hi][x] += 1;
    }
  }
  return { counts, positions };
}

/** Rows that are junctions under the run's rule. */
export function isJunction(t: Table, i: number, size: number, tol: number, maxMm: number): boolean {
  return Math.abs(num(t, "overlap_bp")[i] - size) <= tol && num(t, "n_mismatches")[i] <= maxMm;
}

/** Per-position base counts (A, C, G, T) of junction sequences, weighted by count. */
export function baseCounts(seqs: string[], weights: ArrayLike<number>, length: number): number[][] {
  const idx: Record<string, number> = { A: 0, C: 1, G: 2, T: 3 };
  const out = Array.from({ length }, () => [0, 0, 0, 0]);
  seqs.forEach((s, k) => {
    if (s.length !== length) return;
    for (let p = 0; p < length; p++) {
      const b = idx[s[p]];
      if (b != null) out[p][b] += weights[k];
    }
  });
  return out;
}

/** Information content (bits) of each base at each position, for a sequence logo. */
export function logoHeights(counts: number[][]): number[][] {
  return counts.map((c) => {
    const n = c.reduce((a, b) => a + b, 0);
    if (!n) return [0, 0, 0, 0];
    const f = c.map((x) => x / n);
    const h = -f.reduce((a, p) => a + (p > 0 ? p * Math.log2(p) : 0), 0);
    const ic = Math.max(0, 2 - h);
    return f.map((p) => p * ic);
  });
}

/** Points of the overlap track: axis x, size, type and junction flag. */
export function overlapPoints(t: Table, rows: ArrayLike<number>, g: Genome, r: Region, size: number, tol: number, maxMm: number) {
  const c = str(t, "chrom"), os = num(t, "overlap_start"), oe = num(t, "overlap_end"), bp = num(t, "overlap_bp");
  const ty = str(t, "overlap_type"), h = str(t, "haplotype");
  const pts: { i: number; x: number; bp: number; contained: boolean; junction: boolean; hap: number }[] = [];
  for (let k = 0; k < rows.length; k++) {
    const i = rows[k];
    const x = axisPos(g, r, c[i], (os[i] + oe[i]) / 2);
    if (x == null) continue;
    pts.push({ i, x, bp: bp[i], contained: ty[i] === "contained", junction: isJunction(t, i, size, tol, maxMm), hap: hapIndex(h[i]) });
  }
  return pts;
}

export interface RowFilter {
  minBp: number | null;
  maxBp: number | null;
  maxMm: number | null;
  type: "all" | "end_to_start" | "contained";
  junctionsOnly: boolean;
  phase: string;
  indel: "all" | "with" | "without";
}
export const defaultRowFilter = (): RowFilter => ({ minBp: null, maxBp: null, maxMm: null, type: "all", junctionsOnly: false, phase: "", indel: "all" });

export function filterRows(t: Table, rows: Int32Array, f: RowFilter, size: number, tol: number, maxMm: number): Int32Array {
  const bp = num(t, "overlap_bp"), nm = num(t, "n_mismatches"), ty = str(t, "overlap_type");
  const pb = str(t, "phase_block"), ind = bool(t, "has_indel_in_overlap");
  const q = f.phase.trim();
  const out = new Int32Array(rows.length);
  let k = 0;
  for (let r = 0; r < rows.length; r++) {
    const i = rows[r];
    if (f.minBp != null && bp[i] < f.minBp) continue;
    if (f.maxBp != null && bp[i] > f.maxBp) continue;
    if (f.maxMm != null && nm[i] > f.maxMm) continue;
    if (f.type !== "all" && ty[i] !== f.type) continue;
    if (f.indel !== "all" && (ind[i] === 1) !== (f.indel === "with")) continue;
    if (q && !pb[i].includes(q)) continue;
    if (f.junctionsOnly && !isJunction(t, i, size, tol, maxMm)) continue;
    out[k++] = i;
  }
  return out.subarray(0, k);
}

/** Sort row indices by a column (stable; strings or numbers). */
export function sortRows(t: Table, rows: Int32Array, col: string, desc: boolean): Int32Array {
  const c = t.cols[col];
  const out = Array.from(rows);
  const d = desc ? -1 : 1;
  if (Array.isArray(c)) out.sort((a, b) => d * (c[a] < c[b] ? -1 : c[a] > c[b] ? 1 : 0) || a - b);
  else out.sort((a, b) => d * (c[a] - c[b]) || a - b);
  return Int32Array.from(out);
}
