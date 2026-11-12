// Shared pieces of the genome-position charts: x scale, axis, haplotype lanes, drag-to-zoom.
import { chrLabel, css, FONT, niceTicks } from "../lib/canvas";
import type { Genome, Region } from "../lib/genome";
import type { MouseEvent } from "react";
import { useStore } from "../state/store";

type CanvasMouse = MouseEvent<HTMLCanvasElement>;

export interface XScale {
  x0: number;
  x1: number;
  l: number;
  w: number;
  px: (v: number) => number;
  inv: (p: number) => number;
}

export function xScale(r: Region, l: number, w: number): XScale {
  const x0 = r.start, x1 = Math.max(r.end, r.start + 1);
  return {
    x0, x1, l, w,
    px: (v) => l + ((v - x0) / (x1 - x0)) * w,
    inv: (p) => x0 + ((p - l) / w) * (x1 - x0),
  };
}

/** Axis positions are genome-wide when the region is the whole genome. */
export function regionAxis(g: Genome, r: Region): Region {
  return r.chrom == null ? { chrom: null, start: 0, end: g.length } : r;
}

export function fmtPos(v: number, span: number): string {
  if (span > 5e6) return (v / 1e6).toFixed(span > 5e7 ? 0 : 1) + " Mb";
  if (span > 5e3) return (v / 1e3).toFixed(span > 5e4 ? 0 : 1) + " kb";
  return Math.round(v).toLocaleString("en-US");
}

/** Chromosome boundaries and labels (genome) or position ticks (one chromosome). */
export function drawXAxis(g: CanvasRenderingContext2D, gen: Genome, r: Region, sc: XScale, top: number, bottom: number) {
  g.font = `11px ${FONT}`;
  g.textBaseline = "top";
  g.textAlign = "center";
  if (r.chrom == null) {
    let flip = false;
    for (const c of gen.chroms) {
      const a = sc.px(c.offset), b = sc.px(c.offset + c.length);
      if (flip) {
        g.fillStyle = css("--grid");
        g.globalAlpha = 0.45;
        g.fillRect(a, top, b - a, bottom - top);
        g.globalAlpha = 1;
      }
      flip = !flip;
      if (b - a > 12) {
        g.fillStyle = css("--ink3");
        g.fillText(chrLabel(c.name), (a + b) / 2, bottom + 5);
      }
    }
  } else {
    const span = sc.x1 - sc.x0;
    g.strokeStyle = css("--grid");
    g.lineWidth = 1;
    for (const v of niceTicks(sc.x0, sc.x1, Math.max(3, Math.floor(sc.w / 110)))) {
      const x = sc.px(v);
      g.beginPath(); g.moveTo(x, top); g.lineTo(x, bottom); g.stroke();
      g.fillStyle = css("--ink3");
      g.fillText(fmtPos(v, span), x, bottom + 5);
    }
  }
}

/** A region from two axis coordinates (drag-to-zoom); null when the drag spans chromosomes. */
export function zoomRegion(gen: Genome, r: Region, a: number, b: number): Region | null {
  const lo = Math.min(a, b), hi = Math.max(a, b);
  if (r.chrom != null) {
    const len = gen.byName.get(r.chrom)?.length ?? hi;
    const s = Math.max(0, Math.floor(lo)), e = Math.min(len, Math.ceil(hi));
    return e - s >= 20 ? { chrom: r.chrom, start: s, end: e } : null;
  }
  const c = gen.chroms.find((c) => lo >= c.offset && lo < c.offset + c.length);
  if (!c) return null;
  const e = Math.min(hi, c.offset + c.length);
  return { chrom: c.name, start: Math.max(0, Math.floor(lo - c.offset)), end: Math.ceil(e - c.offset) };
}

/**
 * Mouse handlers for drag-to-zoom on a canvas whose scale is held in `scale`.
 * Draws the selection band on `overlay` (a positioned div) while dragging.
 */
export function dragZoom(gen: Genome, getRegion: () => Region, scale: { current: XScale | null }, overlay: { current: HTMLDivElement | null }) {
  let start: number | null = null;
  const x = (e: CanvasMouse) => e.clientX - e.currentTarget.getBoundingClientRect().left;
  const band = (a: number | null, b = 0) => {
    const o = overlay.current;
    if (!o) return;
    if (a == null) { o.style.display = "none"; return; }
    o.style.display = "block";
    o.style.left = Math.min(a, b) + "px";
    o.style.width = Math.abs(b - a) + "px";
  };
  return {
    onMouseDown: (e: CanvasMouse) => {
      const sc = scale.current;
      if (!sc || e.button !== 0) return;
      const p = x(e);
      if (p < sc.l || p > sc.l + sc.w) return;
      start = p;
    },
    onMouseMoveDrag: (e: CanvasMouse) => {
      if (start != null) band(start, x(e));
      return start != null;
    },
    onMouseUp: (e: CanvasMouse) => {
      const sc = scale.current;
      if (start == null || !sc) return false;
      const p = x(e), s = start;
      start = null;
      band(null);
      if (Math.abs(p - s) < 5) return false;
      const z = zoomRegion(gen, getRegion(), sc.inv(s), sc.inv(p));
      if (z) useStore.getState().set({ region: z });
      return true;
    },
    onLeave: () => { start = null; band(null); },
  };
}
