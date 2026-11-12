// Chromosome layout on one genome-wide axis, and region parsing.

export interface Chrom {
  name: string;
  length: number;
  offset: number;
}
export interface Genome {
  chroms: Chrom[];
  byName: Map<string, Chrom>;
  length: number;
}
/** A window on the axis: one chromosome (or the whole genome) and a range in it. */
export interface Region {
  chrom: string | null; // null = whole genome
  start: number;
  end: number;
}

export function chromKey(c: string): [number, number, string] {
  const n = c.replace(/^chr/i, "");
  if (/^\d+$/.test(n)) return [0, +n, ""];
  const tail: Record<string, number> = { X: 1, Y: 2, M: 3, MT: 3 };
  return [1, tail[n.toUpperCase()] ?? 99, n];
}
export function chromCompare(a: string, b: string): number {
  const x = chromKey(a), y = chromKey(b);
  return x[0] - y[0] || x[1] - y[1] || x[2].localeCompare(y[2]);
}

/** Build the layout; `lengths` in the BAM header order (from meta.json) or inferred. */
export function makeGenome(lengths: [string, number][]): Genome {
  const sorted = [...lengths].sort((a, b) => chromCompare(a[0], b[0]));
  let off = 0;
  const chroms = sorted.map(([name, length]) => {
    const c = { name, length, offset: off };
    off += length;
    return c;
  });
  return { chroms, byName: new Map(chroms.map((c) => [c.name, c])), length: off };
}

export function regionOf(g: Genome, chrom: string | null): Region {
  if (chrom == null || !g.byName.has(chrom)) return { chrom: null, start: 0, end: g.length };
  return { chrom, start: 0, end: g.byName.get(chrom)!.length };
}

/** Axis coordinate of a position: genome-wide when the region is the genome. */
export function axisPos(g: Genome, r: Region, chrom: string, pos: number): number | null {
  if (r.chrom == null) {
    const c = g.byName.get(chrom);
    return c ? c.offset + pos : null;
  }
  return chrom === r.chrom ? pos : null;
}

/** "chr1", "chr1:1,000-2,000", "1:5M-6M" → Region (1-based inclusive in, 0-based half-open out). */
export function parseRegion(g: Genome, text: string): Region | null {
  const t = text.trim().replace(/,/g, "");
  if (!t || /^(all|genome)$/i.test(t)) return regionOf(g, null);
  const m = /^([^:]+?)(?::([\d.]+[kKmM]?)\s*-\s*([\d.]+[kKmM]?))?$/.exec(t);
  if (!m) return null;
  const name = g.byName.has(m[1]) ? m[1] : g.byName.has("chr" + m[1]) ? "chr" + m[1] : g.byName.has(m[1].replace(/^chr/, "")) ? m[1].replace(/^chr/, "") : null;
  if (!name) return null;
  const len = g.byName.get(name)!.length;
  if (!m[2]) return { chrom: name, start: 0, end: len };
  const v = (s: string) => {
    const u = s.slice(-1).toLowerCase();
    const x = parseFloat(s);
    return u === "k" ? x * 1e3 : u === "m" ? x * 1e6 : x;
  };
  const start = Math.max(0, Math.round(v(m[2])) - 1), end = Math.min(len, Math.round(v(m[3])));
  return end > start ? { chrom: name, start, end } : null;
}

export function formatRegion(r: Region): string {
  if (r.chrom == null) return "genome";
  return `${r.chrom}:${(r.start + 1).toLocaleString("en-US")}-${r.end.toLocaleString("en-US")}`;
}
