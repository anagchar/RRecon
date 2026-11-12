// Tn5 junctions: counts per cell and haplotype, the junction sequence logo and 9-mer counts,
// and the mismatch counts and positions of junction-sized overlaps.
import { useEffect, useMemo, useRef, useState, type MouseEvent, type RefObject } from "react";
import { axes, css, emptyMsg, fit, fmtCount, FONT, niceTicks, saveCanvas, tip } from "../lib/canvas";
import { sum } from "../lib/dataset";
import { baseCounts, HAP_COLOR, HAP_LABEL, hapIndex, HAPS, hapsShown, logoHeights, mismatchStats, selectRows } from "../lib/stats";
import { num, str, where } from "../lib/table";
import { pairs, useStore } from "../state/store";
import { useSelection } from "./OverviewTab";
import { Card, DownloadMenu, HapKey, Seg, useThemeKey, useWidth } from "./widgets";

const BASE_COLOR = ["--base-a", "--base-c", "--base-g", "--base-t"];
const BASES = "ACGT";

function Tiles() {
  const { ds, sel } = useSelection();
  const lanes = hapsShown(sel.hap).map(hapIndex);
  const cells = sel.barcode == null ? ds.cells : ds.cells.filter((c) => c.barcode === sel.barcode);
  const j = cells.reduce((a, c) => a + lanes.reduce((s, k) => s + c.junctions[k], 0), 0);
  const o = cells.reduce((a, c) => a + lanes.reduce((s, k) => s + c.overlaps[k], 0), 0);
  const withJ = ds.cells.filter((c) => lanes.some((k) => c.junctions[k] > 0)).length;
  const top = useTopKmers(1)[0];
  return (
    <div className="tiles">
      <div className="tile"><div className="k">Junctions</div><div className="v">{fmtCount(j)}</div><div className="s">{sel.barcode ?? "all cells"} · genome-wide</div></div>
      <div className="tile"><div className="k">Share of overlaps</div><div className="v">{o ? ((100 * j) / o).toFixed(2) + "%" : "–"}</div><div className="s">{fmtCount(o)} overlapping pairs</div></div>
      <div className="tile"><div className="k">Cells with junctions</div><div className="v">{withJ.toLocaleString()}</div><div className="s">of {ds.cells.length.toLocaleString()} cells</div></div>
      <div className="tile"><div className="k">Most common {ds.junctionSize}-mer</div><div className="v long mono">{top ? top[0] : "–"}</div><div className="s">{top ? `${top[1].toLocaleString()} junction${top[1] === 1 ? "" : "s"}` : "no junctions"}</div></div>
    </div>
  );
}

function useTopKmers(limit: number): [string, number][] {
  const { ds, sel } = useSelection();
  return useMemo(() => {
    const t = ds.tables.tn5_kmers;
    if (!t) return [];
    const b = str(t, "barcode"), h = str(t, "haplotype"), s = str(t, "junction_seq"), n = num(t, "count");
    const acc = new Map<string, number>();
    for (let i = 0; i < t.n; i++) {
      if (sel.barcode != null && b[i] !== sel.barcode) continue;
      if (sel.hap !== "all" && h[i] !== sel.hap) continue;
      acc.set(s[i], (acc.get(s[i]) ?? 0) + n[i]);
    }
    return [...acc].sort((x, y) => y[1] - x[1] || x[0].localeCompare(y[0])).slice(0, limit);
  }, [ds, sel.barcode, sel.hap, limit]);
}

function PerCellCard() {
  const { ds, sel } = useSelection();
  const theme = useThemeKey();
  const wrap = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  const W = useWidth(wrap);
  const lanes = hapsShown(sel.hap);
  const N = Math.max(10, Math.min(60, Math.floor((W - 70) / 14)));
  const cells = useMemo(() => {
    const ks = lanes.map(hapIndex);
    return [...ds.cells].sort((a, b) => ks.reduce((s, k) => s + b.junctions[k], 0) - ks.reduce((s, k) => s + a.junctions[k], 0)).slice(0, N);
  }, [ds, lanes, N]);
  const geo = useRef<{ l: number; bw: number } | null>(null);
  const H = 260;
  useEffect(() => {
    const g = fit(cv.current!, W, H);
    if (!cells.length) return emptyMsg(g, W, H, "No cells");
    const max = Math.max(1, ...cells.map((c) => lanes.reduce((s, h) => s + c.junctions[hapIndex(h)], 0)));
    const box = { l: 56, t: 10, w: W - 56 - 12, h: H - 10 - 34 };
    const { ty } = axes(g, box, [], niceTicks(0, max, 4), { x0: 0, x1: 1, y0: 0, y1: max * 1.08, ylab: "junctions", yfmt: fmtCount });
    const bw = box.w / cells.length;
    geo.current = { l: box.l, bw };
    cells.forEach((c, k) => {
      let base = 0;
      const x = box.l + k * bw + 1;
      for (const h of lanes) {
        const v = c.junctions[hapIndex(h)];
        if (!v) continue;
        g.fillStyle = css(HAP_COLOR[h]);
        const y = ty(base + v);
        g.fillRect(x, y, Math.max(1, bw - 2), ty(base) - y);
        base += v;
      }
      if (c.barcode === sel.barcode) {
        g.strokeStyle = css("--ink"); g.lineWidth = 1.5;
        g.strokeRect(x - 1, ty(base) - 1, bw, ty(0) - ty(base) + 2);
      }
    });
    g.fillStyle = css("--ink3"); g.font = `11px ${FONT}`; g.textAlign = "center"; g.textBaseline = "top";
    g.fillText(`top ${cells.length} cells by junction count · click a bar to select the cell`, box.l + box.w / 2, box.t + box.h + 10);
  }, [W, cells, lanes, theme, sel.barcode]);
  const at = (e: MouseEvent<HTMLCanvasElement>) => {
    const gm = geo.current;
    if (!gm) return null;
    const k = Math.floor((e.clientX - e.currentTarget.getBoundingClientRect().left - gm.l) / gm.bw);
    return cells[k] ?? null;
  };
  return (
    <Card title="Junctions per cell" sub="stacked by haplotype" right={<><HapKey hap={sel.hap} /><DownloadMenu onSave={(f) => saveCanvas([cv.current], "junctions_per_cell", f)} /></>}>
      <div ref={wrap}>
        <canvas
          ref={cv}
          style={{ cursor: "pointer" }}
          onClick={(e) => { const c = at(e); if (c) useStore.getState().set({ barcode: c.barcode }); }}
          onMouseLeave={() => tip(null)}
          onMouseMove={(e) => {
            const c = at(e);
            tip(e, c ? `<b class="mono">${c.barcode}</b><br>` + HAPS.map((h) => `${HAP_LABEL[h]}: <b>${c.junctions[hapIndex(h)]}</b> junctions of ${c.overlaps[hapIndex(h)].toLocaleString()} overlaps`).join("<br>") + `<br><span class="k">reads</span> ${sum(c.reads).toLocaleString()}` : null);
          }}
        />
      </div>
    </Card>
  );
}

function LogoCard() {
  const { ds, sel } = useSelection();
  const theme = useThemeKey();
  const wrap = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  const W = useWidth(wrap, 400);
  const top = useTopKmers(12);
  const [mode, setMode] = useState<"freq" | "bits">("freq");
  const counts = useMemo(() => {
    const t = ds.tables.tn5_kmers;
    if (!t) return null;
    const rows = where(t, (i) => (sel.barcode == null || str(t, "barcode")[i] === sel.barcode) && (sel.hap === "all" || str(t, "haplotype")[i] === sel.hap));
    const seqs = Array.from(rows, (i) => str(t, "junction_seq")[i]);
    const w = Array.from(rows, (i) => num(t, "count")[i]);
    return baseCounts(seqs, w, ds.junctionSize);
  }, [ds, sel.barcode, sel.hap]);
  const H = 200;
  useEffect(() => {
    const g = fit(cv.current!, W, H);
    if (!counts) return emptyMsg(g, W, H, "No tn5_kmers file loaded");
    if (!counts.some((c) => c.some((v) => v))) return emptyMsg(g, W, H, "No junctions in this selection");
    const bits = mode === "bits";
    const heights = bits ? logoHeights(counts) : counts.map((c) => { const n = c.reduce((a, b) => a + b, 0) || 1; return c.map((v) => v / n); });
    const box = { l: 40, t: 10, w: W - 40 - 10, h: H - 10 - 30 };
    const { ty } = axes(g, box, [], bits ? [0, 1, 2] : [0, 0.5, 1], { x0: 0, x1: 1, y0: 0, y1: bits ? 2 : 1, ylab: bits ? "bits" : "frequency" });
    const cw = box.w / counts.length;
    heights.forEach((hs, p) => {
      const order = [0, 1, 2, 3].sort((a, b) => hs[a] - hs[b]);
      let y = ty(0);
      for (const b of order) {
        const hpx = ty(0) - ty(hs[b]);
        if (hpx < 0.5) continue;
        drawLetter(g, BASES[b], box.l + p * cw + 2, y - hpx, cw - 4, hpx, css(BASE_COLOR[b]));
        y -= hpx;
      }
      g.fillStyle = css("--ink3"); g.font = `11px ${FONT}`; g.textAlign = "center"; g.textBaseline = "top";
      g.fillText(String(p + 1), box.l + (p + 0.5) * cw, box.t + box.h + 6);
    });
  }, [W, counts, theme, mode]);
  return (
    <Card title={`Junction ${ds.junctionSize}-mers`} sub={`sequence logo (read 1's bases) · ${sel.barcode ?? "all cells"}`}
      right={<><Seg label="Logo scale" value={mode} onChange={setMode} options={[["freq", "frequency"], ["bits", "bits"]]} /><DownloadMenu onSave={(f) => saveCanvas([cv.current], "junction_logo", f)} /></>}>
      <div ref={wrap}><canvas ref={cv} /></div>
      {top.length > 0 && (
        <div className="tablewrap kmers">
          <table>
            <thead><tr><th>{ds.junctionSize}-mer</th><th>junctions</th><th /></tr></thead>
            <tbody>
              {top.map(([k, n]) => (
                <tr key={k}>
                  <td className="mono seq">{k.split("").map((b, i) => <span key={i} style={{ color: `var(${BASE_COLOR[BASES.indexOf(b)] ?? "--ink2"})` }}>{b}</span>)}</td>
                  <td>{n.toLocaleString()}</td>
                  <td style={{ width: "40%" }}><span className="bar" style={{ width: `${(100 * n) / top[0][1]}%` }} /></td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </Card>
  );
}

/** A letter stretched to fill a box (the usual sequence-logo glyph). */
function drawLetter(g: CanvasRenderingContext2D, ch: string, x: number, y: number, w: number, h: number, color: string) {
  g.save();
  g.font = `bold 100px ${FONT}`;
  const m = g.measureText(ch);
  const asc = m.actualBoundingBoxAscent, desc = m.actualBoundingBoxDescent;
  g.translate(x, y);
  g.scale(w / m.width, h / (asc + desc));
  g.fillStyle = color;
  g.textBaseline = "alphabetic";
  g.textAlign = "left";
  g.fillText(ch, 0, asc);
  g.restore();
}

function MismatchCards() {
  const { ds, sel } = useSelection();
  const theme = useThemeKey();
  const t = pairs(ds);
  const lanes = hapsShown(sel.hap);
  const MAXC = Math.max(3, ds.maxMismatches + 2);
  const stats = useMemo(() => {
    if (!t) return null;
    const bp = num(t, "overlap_bp");
    const rows = selectRows(t, ds.genome, sel, ["overlap_start", "overlap_end"]).filter((i) => Math.abs(bp[i] - ds.junctionSize) <= ds.tolerance);
    return { n: rows.length, ...mismatchStats(t, rows, MAXC, ds.junctionSize + ds.tolerance) };
  }, [t, ds, sel.barcode, sel.hap, sel.region, MAXC]); // eslint-disable-line react-hooks/exhaustive-deps
  const onlyJunctions = !ds.tables.overlaps && !!ds.tables.tn5_junctions;
  const cvCounts = useRef<HTMLCanvasElement>(null), cvPos = useRef<HTMLCanvasElement>(null);
  return (
    <div className="grid2">
      <Card title="Mismatches per junction" sub={`${ds.junctionSize}${ds.tolerance ? ` ± ${ds.tolerance}` : ""} bp overlaps${onlyJunctions ? " (junctions only: load overlaps for all)" : ""}`}
        right={<><HapKey hap={sel.hap} /><DownloadMenu onSave={(f) => saveCanvas([cvCounts.current], "mismatch_counts", f)} /></>}>
        <Bars cv={cvCounts} theme={theme} lanes={lanes} data={stats?.counts ?? null} labels={Array.from({ length: MAXC + 1 }, (_, k) => (k === MAXC ? `≥${k}` : String(k)))} xlab="mismatches" cutoff={ds.maxMismatches} />
      </Card>
      <Card title="Mismatch positions" sub="offset within the junction (1 = its first base)"
        right={<><HapKey hap={sel.hap} /><DownloadMenu onSave={(f) => saveCanvas([cvPos.current], "mismatch_positions", f)} /></>}>
        <Bars cv={cvPos} theme={theme} lanes={lanes} data={stats?.positions ?? null} labels={Array.from({ length: ds.junctionSize + ds.tolerance }, (_, k) => String(k + 1))} xlab="position" />
      </Card>
    </div>
  );
}

/** Grouped bars, one group per label, one bar per haplotype. */
function Bars({ lanes, data, labels, xlab, cutoff, cv, theme }: { lanes: (typeof HAPS)[number][]; data: number[][] | null; labels: string[]; xlab: string; cutoff?: number; cv: RefObject<HTMLCanvasElement | null>; theme: string }) {
  const wrap = useRef<HTMLDivElement>(null);
  const W = useWidth(wrap, 400);
  const H = 220;
  useEffect(() => {
    const g = fit(cv.current!, W, H);
    if (!data) return emptyMsg(g, W, H, "No overlaps file loaded");
    const max = Math.max(...lanes.map((h) => Math.max(...data[hapIndex(h)])));
    if (!max) return emptyMsg(g, W, H, "Nothing to show");
    const box = { l: 52, t: 10, w: W - 52 - 10, h: H - 10 - 40 };
    const { ty } = axes(g, box, [], niceTicks(0, max, 4), { x0: 0, x1: 1, y0: 0, y1: max * 1.1, xlab, yfmt: fmtCount });
    const gw = box.w / labels.length, bw = (gw - 6) / lanes.length;
    labels.forEach((lab, k) => {
      lanes.forEach((h, li) => {
        const v = data[hapIndex(h)][k];
        g.fillStyle = css(HAP_COLOR[h]);
        const y = ty(v);
        g.fillRect(box.l + k * gw + 3 + li * bw, y, Math.max(1, bw - 1), ty(0) - y);
      });
      g.fillStyle = css("--ink3"); g.font = `11px ${FONT}`; g.textAlign = "center"; g.textBaseline = "top";
      g.fillText(lab, box.l + (k + 0.5) * gw, box.t + box.h + 5);
    });
    if (cutoff != null) {
      // between the last allowed mismatch count and the next
      const x = box.l + (cutoff + 0.5) * gw;
      g.strokeStyle = css("--junction"); g.setLineDash([4, 3]); g.lineWidth = 1.5;
      g.beginPath(); g.moveTo(x, box.t); g.lineTo(x, box.t + box.h); g.stroke(); g.setLineDash([]);
      g.fillStyle = css("--junction"); g.font = `600 11px ${FONT}`; g.textAlign = "left"; g.textBaseline = "top";
      g.fillText("junction cutoff", x + 4, box.t);
    }
  }, [W, data, lanes, labels, xlab, cutoff, theme, cv]);
  return (
    <div ref={wrap}>
      <canvas ref={cv} />
    </div>
  );
}

export function Tn5Tab() {
  return (
    <>
      <Tiles />
      <PerCellCard />
      <LogoCard />
      <MismatchCards />
    </>
  );
}
