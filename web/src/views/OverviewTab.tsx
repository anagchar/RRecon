// Overview: density track, overlap locations and overlap sizes, each split by haplotype.
import { useEffect, useMemo, useRef, type MouseEvent } from "react";
import { axes, css, emptyMsg, fit, fmtCount, fmtTick, FONT, logTicks, saveCanvas, tip } from "../lib/canvas";
import { axisPos, formatRegion } from "../lib/genome";
import { HAP_COLOR, HAP_LABEL, hapIndex, HAPS, hapsShown, overlapPoints, selectRows, sizeHistogram } from "../lib/stats";
import { num, str } from "../lib/table";
import { pairs, useStore } from "../state/store";
import { dragZoom, drawXAxis, fmtPos, xScale, type XScale } from "./track";
import { Card, DownloadMenu, HapKey, Seg, useThemeKey, useWidth } from "./widgets";

const L = 84, R = 14, LANE = 64, GAP = 10, TOP = 6, AXIS = 22;

export function useSelection() {
  const ds = useStore((s) => s.ds)!;
  const barcode = useStore((s) => s.barcode);
  const hap = useStore((s) => s.hap);
  const region = useStore((s) => s.region);
  return { ds, sel: { barcode, hap, region } };
}

function laneLabel(g: CanvasRenderingContext2D, h: (typeof HAPS)[number], y: number, sub: string) {
  g.fillStyle = css(HAP_COLOR[h]);
  g.beginPath(); g.arc(10, y + 9, 4, 0, 7); g.fill();
  g.fillStyle = css("--ink"); g.font = `600 12px ${FONT}`; g.textAlign = "left"; g.textBaseline = "middle";
  g.fillText(HAP_LABEL[h], 20, y + 9);
  g.fillStyle = css("--ink3"); g.font = `11px ${FONT}`;
  g.fillText(sub, 20, y + 25);
}

function ResetZoom() {
  const ds = useStore((s) => s.ds)!;
  const region = useStore((s) => s.region);
  if (region.chrom == null) return null;
  const whole = region.start === 0 && region.end === ds.genome.byName.get(region.chrom)?.length;
  return (
    <button className="zoomchip" onClick={() => useStore.getState().set({ region: whole ? { chrom: null, start: 0, end: ds.genome.length } : { chrom: region.chrom, start: 0, end: ds.genome.byName.get(region.chrom!)!.length } })}>
      {formatRegion(region)} · {whole ? "back to genome" : "whole chromosome"} ✕
    </button>
  );
}

export function DensityCard() {
  const { ds, sel } = useSelection();
  const theme = useThemeKey();
  const wrap = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  const band = useRef<HTMLDivElement>(null);
  const scale = useRef<XScale | null>(null);
  const hits = useRef<{ lane: number; top: number; items: [number, number, string][] }[]>([]);
  const W = useWidth(wrap);
  const lanes = hapsShown(sel.hap);
  const span = sel.region.end - sel.region.start;
  const useSegs = sel.barcode != null && sel.region.chrom != null && span <= 5e6 && !!ds.tables.density_segments;

  const data = useMemo(() => {
    const out = lanes.map(() => [] as [number, number, number, string][]);
    if (useSegs) {
      const t = ds.tables.density_segments!;
      const rows = selectRows(t, ds.genome, sel, ["start", "end"]);
      const h = str(t, "haplotype"), s = num(t, "start"), e = num(t, "end"), d = num(t, "depth");
      for (const i of rows) {
        const li = lanes.indexOf(h[i] as never);
        if (li >= 0) out[li].push([s[i], e[i], d[i], `${s[i].toLocaleString()}–${e[i].toLocaleString()}`]);
      }
      return { out, unit: "depth" };
    }
    const t = ds.tables.density_bins;
    if (!t) return { out, unit: "reads" };
    const rows = selectRows(t, ds.genome, sel, ["bin_start", "bin_end"]);
    const h = str(t, "haplotype"), c = str(t, "chrom"), s = num(t, "bin_start"), e = num(t, "bin_end"), n = num(t, "n_reads");
    const acc = lanes.map(() => new Map<number, [number, number, number, string]>());
    for (const i of rows) {
      const li = lanes.indexOf(h[i] as never);
      if (li < 0) continue;
      const x = axisPos(ds.genome, sel.region, c[i], s[i])!;
      const cur = acc[li].get(x);
      if (cur) cur[2] += n[i];
      else acc[li].set(x, [x, x + (e[i] - s[i]), n[i], `${c[i]}:${fmtPos(s[i], 1e9)}–${fmtPos(e[i], 1e9)}`]);
    }
    return { out: acc.map((m) => [...m.values()]), unit: "reads / bin" };
  }, [ds, sel.barcode, sel.hap, sel.region, useSegs]); // eslint-disable-line react-hooks/exhaustive-deps

  const H = TOP + lanes.length * (LANE + GAP) + AXIS;
  useEffect(() => {
    const g = fit(cv.current!, W, H);
    const sc = xScale(sel.region, L, W - L - R);
    scale.current = sc;
    const max = Math.max(1, ...data.out.flat().map((d) => d[2]));
    drawXAxis(g, ds.genome, sel.region, sc, TOP, H - AXIS);
    hits.current = [];
    lanes.forEach((h, li) => {
      const y0 = TOP + li * (LANE + GAP);
      laneLabel(g, h, y0, `max ${fmtCount(max)}`);
      g.fillStyle = css(HAP_COLOR[h]);
      const items: [number, number, string][] = [];
      for (const [a, b, v, label] of data.out[li]) {
        if (b < sc.x0 || a > sc.x1) continue;
        const x = sc.px(Math.max(a, sc.x0)), x2 = sc.px(Math.min(b, sc.x1));
        const hh = (v / max) * LANE;
        g.fillRect(x, y0 + LANE - hh, Math.max(1, x2 - x), hh);
        items.push([x, Math.max(x + 1, x2), `${label}<br><span class="k">${data.unit}</span> <b>${v.toLocaleString()}</b>`]);
      }
      g.fillStyle = css("--line2");
      g.fillRect(L, y0 + LANE, W - L - R, 1);
      hits.current.push({ lane: li, top: y0, items });
    });
    if (!data.out.some((d) => d.length)) emptyMsg(g, W, H - AXIS, ds.tables.density_bins ? "No reads in this selection" : "No density file loaded");
  }, [W, H, data, theme, sel.region, ds, lanes]);

  const z = dragZoom(ds.genome, () => useStore.getState().region, scale, band);
  return (
    <Card
      title="Read density"
      sub={<>{useSegs ? "Read depth" : "Reads per bin"} for {sel.barcode ?? "all cells"} · drag to zoom</>}
      right={<><HapKey hap={sel.hap} /><ResetZoom /><DownloadMenu onSave={(f) => saveCanvas([cv.current], "density", f)} /></>}
    >
      <div ref={wrap} className="hmwrap">
        <canvas
          ref={cv}
          onMouseDown={z.onMouseDown}
          onMouseUp={z.onMouseUp}
          onDoubleClick={() => useStore.getState().set({ region: { chrom: null, start: 0, end: ds.genome.length } })}
          onMouseLeave={() => { z.onLeave(); tip(null); }}
          onMouseMove={(e) => {
            if (z.onMouseMoveDrag(e)) return tip(null);
            const r = e.currentTarget.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
            const lane = hits.current.find((l) => y >= l.top && y <= l.top + LANE);
            const it = lane?.items.find(([a, b]) => x >= a - 0.5 && x <= b + 0.5);
            tip(e, it ? `<b>${HAP_LABEL[lanes[lane!.lane]]}</b><br>${it[2]}` : null);
          }}
        />
        <div ref={band} className="zoomband" />
      </div>
    </Card>
  );
}

export function OverlapTrackCard() {
  const { ds, sel } = useSelection();
  const theme = useThemeKey();
  const wrap = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  const band = useRef<HTMLDivElement>(null);
  const scale = useRef<XScale | null>(null);
  const hitPts = useRef<{ x: number; y: number; i: number }[]>([]);
  const W = useWidth(wrap);
  const t = pairs(ds);
  const lanes = hapsShown(sel.hap);
  const pts = useMemo(() => {
    if (!t) return [];
    const rows = selectRows(t, ds.genome, sel, ["overlap_start", "overlap_end"]);
    return overlapPoints(t, rows, ds.genome, sel.region, ds.junctionSize, ds.tolerance, ds.maxMismatches);
  }, [t, ds, sel.barcode, sel.hap, sel.region]); // eslint-disable-line react-hooks/exhaustive-deps
  const LH = 92;
  const H = TOP + lanes.length * (LH + GAP) + AXIS;

  useEffect(() => {
    const g = fit(cv.current!, W, H);
    const sc = xScale(sel.region, L, W - L - R);
    scale.current = sc;
    drawXAxis(g, ds.genome, sel.region, sc, TOP, H - AXIS);
    const maxBp = Math.max(10, ...pts.map((p) => p.bp));
    const ly = (bp: number, y0: number) => y0 + LH - 4 - (Math.log10(Math.max(1, bp)) / Math.log10(maxBp)) * (LH - 10);
    const hit: { x: number; y: number; i: number }[] = [];
    lanes.forEach((h, li) => {
      const y0 = TOP + li * (LH + GAP);
      const own = pts.filter((p) => p.hap === hapIndex(h));
      laneLabel(g, h, y0, `${fmtCount(own.length)} pairs`);
      g.font = `10px ${FONT}`; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
      for (const v of logTicks(1, maxBp).filter((v) => [1, 10, 100, 1000, 10000, 100000].includes(v))) {
        const y = ly(v, y0);
        g.fillText(fmtTick(v), L - 4, y);
        g.fillStyle = css("--grid"); g.fillRect(L, y, W - L - R, 1); g.fillStyle = css("--ink3");
      }
      const jy = ly(ds.junctionSize, y0);
      g.fillStyle = css("--junction"); g.globalAlpha = 0.12; g.fillRect(L, jy - 3, W - L - R, 6); g.globalAlpha = 1;
      // plain overlaps first, junctions on top
      for (const pass of [0, 1])
        for (const p of own) {
          if (+p.junction !== pass) continue;
          const x = sc.px(p.x), y = ly(p.bp, y0);
          if (x < L || x > W - R) continue;
          if (p.junction) {
            g.fillStyle = css("--junction");
            g.beginPath(); g.arc(x, y, 3.6, 0, 7); g.fill();
          } else {
            g.fillStyle = css(p.contained ? "--contained" : "--endstart");
            g.globalAlpha = 0.65;
            g.beginPath(); g.arc(x, y, 2.4, 0, 7); g.fill();
            g.globalAlpha = 1;
          }
          hit.push({ x, y, i: p.i });
        }
      g.fillStyle = css("--line2"); g.fillRect(L, y0 + LH, W - L - R, 1);
    });
    hitPts.current = hit;
    if (!t) emptyMsg(g, W, H - AXIS, "No overlaps file loaded");
    else if (!pts.length) emptyMsg(g, W, H - AXIS, "No overlaps in this selection");
  }, [W, H, pts, theme, sel.region, ds, lanes, t]);

  const nearest = (e: MouseEvent<HTMLCanvasElement>) => {
    const r = e.currentTarget.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
    let best = -1, bd = 49;
    for (const p of hitPts.current) {
      const d = (p.x - x) ** 2 + (p.y - y) ** 2;
      if (d < bd) { bd = d; best = p.i; }
    }
    return best;
  };
  const z = dragZoom(ds.genome, () => useStore.getState().region, scale, band);
  return (
    <Card
      title="Overlap locations"
      sub={<>y = overlap size (log) · <span className="dotkey"><i style={{ background: "var(--endstart)" }} />end-to-start</span> <span className="dotkey"><i style={{ background: "var(--contained)" }} />contained</span> <span className="dotkey"><i style={{ background: "var(--junction)" }} />{ds.junctionSize} bp junction</span> · click a dot to view the reads</>}
      right={<><ResetZoom /><DownloadMenu onSave={(f) => saveCanvas([cv.current], "overlap_locations", f)} /></>}
    >
      <div ref={wrap} className="hmwrap">
        <canvas
          ref={cv}
          onMouseDown={z.onMouseDown}
          onMouseUp={(e) => {
            if (z.onMouseUp(e)) return;
            const i = nearest(e);
            if (i >= 0) useStore.getState().focusRow(i);
          }}
          onMouseLeave={() => { z.onLeave(); tip(null); }}
          onMouseMove={(e) => {
            if (z.onMouseMoveDrag(e)) return tip(null);
            const i = nearest(e);
            e.currentTarget.style.cursor = i >= 0 ? "pointer" : "crosshair";
            tip(e, i >= 0 && t ? pairTip(t, i) : null);
          }}
        />
        <div ref={band} className="zoomband" />
      </div>
    </Card>
  );
}

export function pairTip(t: NonNullable<ReturnType<typeof pairs>>, i: number) {
  const s = (c: string) => str(t, c)[i], n = (c: string) => num(t, c)[i];
  return `<b>${s("chrom")}:${(n("overlap_start") + 1).toLocaleString()}-${n("overlap_end").toLocaleString()}</b><br>` +
    `<span class="k">size</span> ${n("overlap_bp")} bp · ${s("overlap_type").replace(/_/g, "-")}<br>` +
    `<span class="k">mismatches</span> ${n("n_mismatches")}${s("mismatch_positions") !== "." ? ` (at ${s("mismatch_positions")})` : ""}<br>` +
    `<span class="k">cell</span> <span class="mono">${s("barcode")}</span> · ${HAP_LABEL[s("haplotype") as "1"] ?? s("haplotype")}` +
    (s("phase_block") !== "." ? ` · PS ${s("phase_block")}` : "") +
    (s("seq1") !== "." && n("overlap_bp") <= 30 ? `<br><span class="mono">${s("seq1")}<br>${s("seq2")}</span>` : "");
}

export function SizeHistCard() {
  const { ds, sel } = useSelection();
  const theme = useThemeKey();
  const histMax = useStore((s) => s.histMax);
  const wrap = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  const W = useWidth(wrap);
  const t = pairs(ds);
  const lanes = hapsShown(sel.hap);
  const hist = useMemo(() => {
    if (!t) return null;
    return sizeHistogram(t, selectRows(t, ds.genome, sel, ["overlap_start", "overlap_end"]), histMax);
  }, [t, ds, sel.barcode, sel.hap, sel.region, histMax]); // eslint-disable-line react-hooks/exhaustive-deps
  const H = 280;
  const geo = useRef<{ tx: (v: number) => number; box: { l: number; w: number } } | null>(null);

  useEffect(() => {
    const g = fit(cv.current!, W, H);
    if (!hist) return emptyMsg(g, W, H, "No overlaps file loaded");
    const ymax = Math.max(1, ...lanes.map((h) => Math.max(...hist[hapIndex(h)])));
    const box = { l: 56, t: 12, w: W - 56 - 16, h: H - 12 - 44 };
    const { tx, ty } = axes(g, box, [1, ...niceX(histMax)], logTicks(1, ymax * 1.2), {
      x0: 0.5, x1: histMax + 0.5, y0: 0.8, y1: ymax * 1.4, ylog: true,
      xlab: `overlap size (bp)${histMax ? ` — last bar: ≥ ${histMax}` : ""}`, ylab: "read pairs",
      yfmt: fmtTick, xfmt: (v) => String(v),
    });
    geo.current = { tx, box };
    // the junction size band
    const j0 = tx(ds.junctionSize - ds.tolerance - 0.5), j1 = tx(ds.junctionSize + ds.tolerance + 0.5);
    g.fillStyle = css("--junction"); g.globalAlpha = 0.14; g.fillRect(j0, box.t, j1 - j0, box.h); g.globalAlpha = 1;
    g.fillStyle = css("--junction"); g.font = `600 11px ${FONT}`; g.textAlign = "center"; g.textBaseline = "bottom";
    g.fillText(`${ds.junctionSize} bp`, (j0 + j1) / 2, box.t + 12);
    const bw = box.w / histMax;
    const shown = lanes.length;
    lanes.forEach((h, li) => {
      const counts = hist[hapIndex(h)];
      g.fillStyle = css(HAP_COLOR[h]);
      for (let bp = 1; bp <= histMax; bp++) {
        const v = counts[bp];
        if (!v) continue;
        const x = tx(bp - 0.5) + (li * bw) / shown, y = ty(v);
        g.fillRect(x, y, Math.max(1, bw / shown - (bw > 6 ? 1 : 0)), box.t + box.h - y);
      }
    });
  }, [W, hist, theme, lanes, histMax, ds]);

  return (
    <Card
      title="Overlap sizes"
      sub={<>{sel.barcode ?? "all cells"} · {formatRegion(sel.region)} · log y</>}
      right={<>
        <HapKey hap={sel.hap} />
        <Seg label="Largest size shown" value={histMax} options={[[30, "30"], [60, "60"], [200, "200"], [1000, "1k"]]} onChange={(v) => useStore.getState().set({ histMax: v })} />
        <DownloadMenu onSave={(f) => saveCanvas([cv.current], "overlap_sizes", f)} />
      </>}
    >
      <div ref={wrap}>
        <canvas
          ref={cv}
          onMouseLeave={() => tip(null)}
          onMouseMove={(e) => {
            const gm = geo.current;
            if (!gm || !hist) return;
            const x = e.clientX - e.currentTarget.getBoundingClientRect().left;
            const bp = Math.round(0.5 + ((x - gm.box.l) / gm.box.w) * histMax);
            if (bp < 1 || bp > histMax) return tip(null);
            tip(e, `<b>${bp === histMax ? "≥ " : ""}${bp} bp</b><br>` + lanes.map((h) => `${HAP_LABEL[h]}: <b>${hist[hapIndex(h)][bp].toLocaleString()}</b>`).join("<br>"));
          }}
        />
      </div>
    </Card>
  );
}

function niceX(max: number): number[] {
  const step = max <= 30 ? 5 : max <= 60 ? 10 : max <= 200 ? 25 : 100;
  const out = [];
  for (let v = step; v <= max; v += step) out.push(v);
  return out;
}

export function OverviewTab() {
  return (
    <>
      <DensityCard />
      <OverlapTrackCard />
      <SizeHistCard />
    </>
  );
}
