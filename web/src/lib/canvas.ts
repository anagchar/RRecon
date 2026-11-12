// Canvas and formatting helpers shared by every chart (from cnvscope, itself ported from the ASCAT.sc
// Explorer's chart.js), so charts read the theme from the CSS tokens.

export const FONT = "'Helvetica Neue', Helvetica, Arial, sans-serif";

export const css = (token: string) =>
  token.startsWith("--") ? getComputedStyle(document.documentElement).getPropertyValue(token).trim() : token;

/** Resolve a token or colour to RGB, whatever CSS colour syntax it uses. */
let probe: CanvasRenderingContext2D | null = null;
export function rgbOf(token: string): [number, number, number] {
  probe ??= document.createElement("canvas").getContext("2d", { willReadFrequently: true })!;
  probe.clearRect(0, 0, 1, 1);
  probe.fillStyle = css(token) || "#888";
  probe.fillRect(0, 0, 1, 1);
  const d = probe.getImageData(0, 0, 1, 1).data;
  return [d[0], d[1], d[2]];
}
export const hexRGB = (h: string): [number, number, number] => [
  parseInt(h.slice(1, 3), 16),
  parseInt(h.slice(3, 5), 16),
  parseInt(h.slice(5, 7), 16),
];
export const mixRGB = (a: number[], b: number[], t: number) =>
  `rgb(${a.map((v, k) => Math.round(v + (b[k] - v) * t)).join(",")})`;

/** Size a canvas for the device pixel ratio; returns its 2D context in CSS pixels. */
export function fit(cv: HTMLCanvasElement, w: number, h: number): CanvasRenderingContext2D {
  const dpr = Math.max(1, window.devicePixelRatio || 1);
  cv.width = Math.round(w * dpr);
  cv.height = Math.round(h * dpr);
  cv.style.width = w + "px";
  cv.style.height = h + "px";
  const g = cv.getContext("2d")!;
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  return g;
}

type Num = number | null | undefined;
const bad = (v: Num): v is null | undefined => v == null || !Number.isFinite(v);
export const fmt = (v: Num, d = 2) => (bad(v) ? "–" : v.toFixed(d));
export const fmtInt = (v: Num) => (bad(v) ? "–" : Math.round(v).toLocaleString());
/** compact counts: 950, 12.0k, 3.62M */
export const fmtCount = (v: Num) =>
  bad(v) ? "–" : v < 1e3 ? String(Math.round(v)) : v < 1e6 ? (v / 1e3).toFixed(v < 1e4 ? 2 : 1) + "k" : (v / 1e6).toFixed(v < 1e7 ? 2 : 1) + "M";
/** axis labels for log scales: 1, 10, 100, 1k, 10k, 1M */
export const fmtTick = (v: number) =>
  v >= 1e6 ? +(v / 1e6).toFixed(1) + "M" : v >= 1e3 ? +(v / 1e3).toFixed(1) + "k" : String(+v.toFixed(2));
export const fmtBp = (b: number) => (b >= 1e6 ? (b / 1e6).toFixed(b % 1e6 ? 2 : 0) + " Mb" : (b / 1e3).toFixed(0) + " kb");
export const chrLabel = (c: string) => String(c).replace(/^chr/, "");
export const esc = (s: string) =>
  String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]!);

export function quantile(sorted: ArrayLike<number>, p: number): number {
  if (!sorted.length) return NaN;
  const i = (sorted.length - 1) * p, a = Math.floor(i), b = Math.ceil(i);
  return sorted[a] + (sorted[b] - sorted[a]) * (i - a);
}
export function niceTicks(a: number, b: number, n = 5): number[] {
  if (!(b > a)) return [a];
  const s = Math.pow(10, Math.floor(Math.log10((b - a) / n)));
  const st = [1, 2, 5, 10].map((m) => m * s).find((v) => (b - a) / v <= n) || s * 10;
  const out: number[] = [];
  for (let v = Math.ceil(a / st) * st; v <= b + 1e-12; v += st) out.push(+v.toFixed(10));
  return out;
}
export function logTicks(a: number, b: number): number[] {
  const out: number[] = [];
  for (let e = Math.floor(Math.log10(a)); e <= Math.ceil(Math.log10(b)); e++)
    for (const m of [1, 2, 5]) {
      const v = m * Math.pow(10, e);
      if (v >= a && v <= b) out.push(v);
    }
  return out;
}
export function axisTicks(a: number, b: number) {
  const t = logTicks(a, b);
  return t.length >= 3 ? t : niceTicks(a, b, 5).filter((v) => v > 0);
}

export interface Box { l: number; t: number; w: number; h: number }
export interface AxesOpt {
  x0: number; x1: number; y0: number; y1: number;
  xlog?: boolean; ylog?: boolean; xlab?: string; ylab?: string;
  xfmt?: (v: number) => string; yfmt?: (v: number) => string;
}
/** Gridlines, tick labels and axis titles; returns the value → pixel maps. */
export function axes(g: CanvasRenderingContext2D, box: Box, xs: number[], ys: number[], opt: AxesOpt) {
  const { x0, x1, y0, y1, xlog, ylog } = opt;
  const tx = (v: number) => box.l + ((xlog ? Math.log10(v) - Math.log10(x0) : v - x0) / (xlog ? Math.log10(x1) - Math.log10(x0) : x1 - x0)) * box.w;
  const ty = (v: number) => box.t + box.h - ((ylog ? Math.log10(v) - Math.log10(y0) : v - y0) / (ylog ? Math.log10(y1) - Math.log10(y0) : y1 - y0)) * box.h;
  g.font = `11px ${FONT}`;
  g.lineWidth = 1;
  g.strokeStyle = css("--grid");
  g.fillStyle = css("--ink3");
  g.textAlign = "right";
  g.textBaseline = "middle";
  for (const v of ys) {
    const y = ty(v);
    g.beginPath(); g.moveTo(box.l, y); g.lineTo(box.l + box.w, y); g.stroke();
    g.fillText(opt.yfmt ? opt.yfmt(v) : String(v), box.l - 6, y);
  }
  g.textAlign = "center";
  g.textBaseline = "top";
  for (const v of xs) g.fillText(opt.xfmt ? opt.xfmt(v) : String(v), tx(v), box.t + box.h + 6);
  g.fillStyle = css("--ink2");
  if (opt.xlab) g.fillText(opt.xlab, box.l + box.w / 2, box.t + box.h + 22);
  if (opt.ylab) {
    g.save(); g.translate(12, box.t + box.h / 2); g.rotate(-Math.PI / 2); g.textBaseline = "middle"; g.fillText(opt.ylab, 0, 0); g.restore();
  }
  return { tx, ty };
}

/** A row of [colour token, label, shape] legend items drawn on the canvas. */
export function legendRow(g: CanvasRenderingContext2D, x: number, y: number, items: [string, string, "dot" | "tri" | "line"][]) {
  g.font = `11.5px ${FONT}`;
  g.textBaseline = "middle";
  g.textAlign = "left";
  for (const [c, t, shape] of items) {
    g.fillStyle = css(c);
    g.beginPath();
    if (shape === "tri") { g.moveTo(x + 4, y - 4); g.lineTo(x + 8, y + 3); g.lineTo(x, y + 3); }
    else if (shape === "line") g.fillRect(x, y - 1.5, 14, 3);
    else g.arc(x + 4, y, 3.5, 0, 7);
    g.fill();
    g.fillStyle = css("--ink2");
    g.fillText(t, x + (shape === "line" ? 20 : 12), y);
    x += g.measureText(t).width + (shape === "line" ? 34 : 26);
  }
}

export function roundTop(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, r: number) {
  if (h <= 0) return;
  r = Math.min(r, w / 2, h);
  g.beginPath(); g.moveTo(x, y + h); g.lineTo(x, y + r); g.quadraticCurveTo(x, y, x + r, y); g.lineTo(x + w - r, y);
  g.quadraticCurveTo(x + w, y, x + w, y + r); g.lineTo(x + w, y + h); g.fill();
}

export function emptyMsg(g: CanvasRenderingContext2D, W: number, H: number, text: string) {
  g.fillStyle = css("--ink3"); g.font = `13px ${FONT}`; g.textAlign = "center"; g.textBaseline = "middle"; g.fillText(text, W / 2, H / 2);
}


/* ---------- tooltip: one fixed element shared by every chart ---------- */
let tipEl: HTMLDivElement | null = null;
export function tip(e: { clientX: number; clientY: number } | null, html?: string | null) {
  if (!tipEl) {
    tipEl = document.createElement("div");
    tipEl.className = "tip";
    document.body.appendChild(tipEl);
  }
  if (!e || !html) { tipEl.style.display = "none"; return; }
  tipEl.innerHTML = html;
  tipEl.style.display = "block";
  const r = tipEl.getBoundingClientRect();
  let x = e.clientX + 14, y = e.clientY + 14;
  if (x + r.width > innerWidth - 8) x = e.clientX - r.width - 14;
  if (y + r.height > innerHeight - 8) y = e.clientY - r.height - 14;
  tipEl.style.left = x + "px";
  tipEl.style.top = y + "px";
}
export const swatch = (c: string) => `<span class="swatch" style="background:${c}"></span>`;

/** Save canvases side by side as PNG or PDF, on the card's background colour. */
export async function saveCanvas(canvases: (HTMLCanvasElement | null | undefined)[], basename: string, format: "png" | "pdf" = "png") {
  const list = canvases.filter((c): c is HTMLCanvasElement => !!c);
  if (!list.length) return;
  const W = list.reduce((s, c) => s + c.width, 0), H = Math.max(...list.map((c) => c.height));
  const out = document.createElement("canvas");
  out.width = W; out.height = H;
  const g = out.getContext("2d")!;
  g.fillStyle = css("--panel");
  g.fillRect(0, 0, W, H);
  let x = 0;
  for (const c of list) { g.drawImage(c, x, 0); x += c.width; }
  if (format === "pdf") {
    const { jsPDF } = await import("jspdf");
    const cssW = list.reduce((s, c) => s + c.offsetWidth, 0), cssH = Math.max(...list.map((c) => c.offsetHeight));
    const pdf = new jsPDF({ orientation: cssW > cssH ? "landscape" : "portrait", unit: "px", format: [cssW, cssH], hotfixes: ["px_scaling"] });
    pdf.addImage(out.toDataURL("image/png"), "PNG", 0, 0, cssW, cssH);
    pdf.save(basename + ".pdf");
  } else {
    const a = document.createElement("a");
    a.download = basename + ".png";
    a.href = out.toDataURL("image/png");
    a.click();
  }
}
