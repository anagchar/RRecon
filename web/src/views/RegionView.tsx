// Zoomable region view: the reads of the selected cell, packed in rows per haplotype, with
// each overlap shaded and mismatches ticked; at base level the focused pair shows its bases.
import { useEffect, useMemo, useRef, useState } from "react";
import { css, emptyMsg, fit, FONT, niceTicks, saveCanvas, tip } from "../lib/canvas";
import { HAP_COLOR, HAP_LABEL, hapIndex, hapsShown, isJunction } from "../lib/stats";
import { bool, num, str, type Table } from "../lib/table";
import { pairs, useStore } from "../state/store";
import { useSelection } from "./OverviewTab";
import { fmtPos } from "./track";
import { Card, DownloadMenu, HapKey, useThemeKey, useWidth } from "./widgets";

const L = 84, R = 14, ROW = 14, ROWGAP = 4, AXIS = 26, MAX_ROWS = 24, MAX_SPAN = 200_000, COV = 30;
// IGV's base colours, as theme tokens
const BASE_TOKEN: Record<string, string> = { A: "--base-a", C: "--base-c", G: "--base-g", T: "--base-t" };

interface ReadBox { key: string; name: string; start: number; end: number; reverse: boolean; hap: number; row: number; pairs: number[] }
interface Win { chrom: string; start: number; end: number }

function windowFor(t: Table, i: number, tight: boolean): Win {
  const os = num(t, "overlap_start")[i], oe = num(t, "overlap_end")[i];
  const a = num(t, "read1_start")[i], b = Math.max(num(t, "read1_end")[i], num(t, "read2_end")[i]);
  if (tight || b - a > 4000) {
    const pad = Math.max(25, Math.round((oe - os) * 0.6));
    return { chrom: str(t, "chrom")[i], start: os - pad, end: oe + pad };
  }
  const pad = Math.round((b - a) * 0.05) + 5;
  return { chrom: str(t, "chrom")[i], start: a - pad, end: b + pad };
}

export function RegionTab() {
  const { ds, sel } = useSelection();
  const focus = useStore((s) => s.focus);
  const t = pairs(ds);
  const [win, setWin] = useState<Win | null>(null);

  // a newly focused pair, or a narrow region from the side panel, sets the window
  useEffect(() => {
    if (t && focus != null && focus < t.n) setWin(windowFor(t, focus, num(t, "overlap_bp")[focus] < 200));
  }, [t, focus]);
  useEffect(() => {
    const r = sel.region;
    if (r.chrom != null && r.end - r.start <= MAX_SPAN) setWin({ chrom: r.chrom, start: r.start, end: r.end });
  }, [sel.region]);

  if (!t) return <div className="card"><div className="empty">Load an overlaps or tn5_junctions file to see reads.</div></div>;
  if (!win) {
    return (
      <div className="card"><div className="empty">
        Click an overlap in <b>Overview</b> or a row in <b>Table</b>, or set a region of at most {fmtPos(MAX_SPAN, MAX_SPAN)} in the side panel.
      </div></div>
    );
  }
  return (
    <>
      <ReadsCard win={win} setWin={setWin} />
      {focus != null && focus < t.n && <PairDetail t={t} i={focus} />}
    </>
  );
}

function ReadsCard({ win, setWin }: { win: Win; setWin: (w: Win) => void }) {
  const { ds, sel } = useSelection();
  const theme = useThemeKey();
  const focus = useStore((s) => s.focus);
  const t = pairs(ds)!;
  const wrap = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  const W = useWidth(wrap);
  const lanes = hapsShown(sel.hap);
  const chromLen = ds.genome.byName.get(win.chrom)?.length ?? win.end;

  // reads and pairs overlapping a generous window around the view, packed per lane
  const { reads, laneRows } = useMemo(() => {
    const span = win.end - win.start;
    const lo = win.start - span, hi = win.end + span;
    const b = str(t, "barcode"), h = str(t, "haplotype"), c = str(t, "chrom");
    const s1 = num(t, "read1_start"), e1 = num(t, "read1_end"), s2 = num(t, "read2_start"), e2 = num(t, "read2_end");
    const n1 = str(t, "read1_name"), n2 = str(t, "read2_name"), st1 = str(t, "read1_strand"), st2 = str(t, "read2_strand");
    const map = new Map<string, ReadBox>();
    const add = (name: string, s: number, e: number, strand: string, hap: number, i: number) => {
      const key = `${hap}|${name}|${s}`;
      let r = map.get(key);
      if (!r) map.set(key, (r = { key, name, start: s, end: e, reverse: strand === "-", hap, row: -1, pairs: [] }));
      r.pairs.push(i);
    };
    for (let i = 0; i < t.n; i++) {
      if (c[i] !== win.chrom) continue;
      if (sel.barcode != null && b[i] !== sel.barcode && i !== focus) continue;
      if (sel.hap !== "all" && h[i] !== sel.hap && i !== focus) continue;
      if (Math.max(e1[i], e2[i]) < lo || s1[i] > hi) continue;
      const hp = hapIndex(h[i]);
      add(n1[i], s1[i], e1[i], st1[i], hp, i);
      add(n2[i], s2[i], e2[i], st2[i], hp, i);
    }
    const reads = [...map.values()].sort((a, b) => a.start - b.start || a.end - b.end);
    const laneRows = [0, 0, 0];
    const ends: number[][] = [[], [], []];
    for (const r of reads) {
      const rows = ends[r.hap];
      let k = rows.findIndex((e) => e + 2 <= r.start);
      if (k < 0) { k = rows.length; rows.push(0); }
      rows[k] = r.end;
      r.row = k;
      laneRows[r.hap] = Math.max(laneRows[r.hap], k + 1);
    }
    return { reads, laneRows };
  }, [t, win, sel.barcode, sel.hap, focus]);

  // IGV-style coverage from the density segments (per cell, so only with one cell selected)
  const cov = useMemo(() => {
    const out: [number, number, number][][] = [[], [], []];
    const d = ds.tables.density_segments;
    if (!d || sel.barcode == null) return null;
    const b = str(d, "barcode"), h = str(d, "haplotype"), c = str(d, "chrom"), s = num(d, "start"), e = num(d, "end"), dp = num(d, "depth");
    for (let i = 0; i < d.n; i++) {
      if (b[i] !== sel.barcode || c[i] !== win.chrom || e[i] <= win.start || s[i] >= win.end) continue;
      out[hapIndex(h[i])].push([s[i], e[i], dp[i]]);
    }
    return out;
  }, [ds, sel.barcode, win]);
  const covH = cov ? COV : 0;
  const laneH = (k: number) => covH + Math.max(1, Math.min(MAX_ROWS, laneRows[k])) * (ROW + ROWGAP) + 26;
  const H = AXIS + lanes.reduce((s, h) => s + laneH(hapIndex(h)) + 8, 0) + 4;
  const hits = useRef<{ x0: number; x1: number; y: number; r: ReadBox }[]>([]);

  useEffect(() => {
    const g = fit(cv.current!, W, H);
    const x0 = win.start, x1 = win.end, w = W - L - R;
    const px = (v: number) => L + ((v - x0) / (x1 - x0)) * w;
    const ppb = w / (x1 - x0);
    // axis
    g.font = `11px ${FONT}`; g.textAlign = "center"; g.textBaseline = "bottom";
    for (const v of niceTicks(x0, x1, Math.max(3, Math.floor(w / 110)))) {
      g.fillStyle = css("--grid"); g.fillRect(px(v), AXIS - 4, 1, H - AXIS);
      g.fillStyle = css("--ink3"); g.fillText(fmtPos(v, x1 - x0), px(v), AXIS - 8);
    }
    g.textAlign = "left"; g.fillStyle = css("--ink2");
    g.fillText(win.chrom, 8, AXIS - 8);

    const showBases = ppb >= 7;
    const os = num(t, "overlap_start"), oe = num(t, "overlap_end"), ty = str(t, "overlap_type");
    const mp = str(t, "mismatch_positions"), seq1 = str(t, "seq1"), seq2 = str(t, "seq2");
    const n1 = str(t, "read1_name"), s1 = num(t, "read1_start");
    const hit: typeof hits.current = [];
    let y0 = AXIS;
    for (const hap of lanes) {
      const k = hapIndex(hap);
      g.fillStyle = css(HAP_COLOR[hap]); g.beginPath(); g.arc(10, y0 + 10, 4, 0, 7); g.fill();
      g.fillStyle = css("--ink"); g.font = `600 12px ${FONT}`; g.textAlign = "left"; g.textBaseline = "middle";
      g.fillText(HAP_LABEL[hap], 20, y0 + 10);
      if (cov) {
        const segs = cov[k], max = Math.max(1, ...segs.map((q) => q[2]));
        const cy = y0 + 20;
        g.fillStyle = css("--field"); g.fillRect(L, cy, W - L - R, COV - 6);
        g.fillStyle = css(HAP_COLOR[hap]); g.globalAlpha = 0.6;
        for (const [a, b, dp] of segs) {
          const xa = Math.max(L, px(a)), xb = Math.min(W - R, px(b)), hh = ((COV - 6) * dp) / max;
          g.fillRect(xa, cy + COV - 6 - hh, Math.max(1, xb - xa), hh);
        }
        g.globalAlpha = 1;
        g.fillStyle = css("--ink3"); g.font = `10px ${FONT}`; g.textAlign = "right"; g.textBaseline = "middle";
        g.fillText(segs.length ? `depth ${max}` : "depth 0", L - 4, cy + (COV - 6) / 2);
      }
      const lane = reads.filter((r) => r.hap === k);
      if (!lane.length) { g.fillStyle = css("--ink3"); g.font = `11px ${FONT}`; g.textAlign = "left"; g.textBaseline = "middle"; g.fillText("no overlapping reads here", L + 8, y0 + 20 + covH + ROW / 2); }
      for (const r of lane) {
        if (r.row >= MAX_ROWS || r.end < x0 || r.start > x1) continue;
        const y = y0 + 20 + covH + r.row * (ROW + ROWGAP);
        const a = Math.max(L, px(r.start)), b = Math.min(W - R, px(r.end));
        const focused = focus != null && r.pairs.includes(focus) && (r.name === n1[focus] && r.start === s1[focus] || r.name === str(t, "read2_name")[focus] && r.start === num(t, "read2_start")[focus]);
        g.fillStyle = css(HAP_COLOR[hap]);
        g.globalAlpha = focused ? 0.55 : 0.28;
        g.fillRect(a, y, Math.max(1, b - a), ROW);
        g.globalAlpha = 1;
        // strand arrow at the read's 3' end
        if (b - a > 10) {
          g.fillStyle = css(HAP_COLOR[hap]);
          const tipX = r.reverse ? px(r.start) : px(r.end);
          if (tipX >= L && tipX <= W - R) {
            g.beginPath();
            if (r.reverse) { g.moveTo(tipX, y + ROW / 2); g.lineTo(tipX + 6, y); g.lineTo(tipX + 6, y + ROW); }
            else { g.moveTo(tipX, y + ROW / 2); g.lineTo(tipX - 6, y); g.lineTo(tipX - 6, y + ROW); }
            g.fill();
          }
        }
        if (focused) { g.strokeStyle = css("--ink"); g.lineWidth = 1.5; g.strokeRect(a, y + 0.5, Math.max(1, b - a), ROW - 1); }
        // overlaps on this read, and mismatches
        for (const i of r.pairs) {
          const oa = Math.max(L, px(os[i])), ob = Math.min(W - R, px(oe[i]));
          if (ob < L || oa > W - R) continue;
          g.fillStyle = css(isJunction(t, i, ds.junctionSize, ds.tolerance, ds.maxMismatches) ? "--junction" : ty[i] === "contained" ? "--contained" : "--endstart");
          g.globalAlpha = 0.85;
          g.fillRect(oa, y + ROW - 4, Math.max(1.5, ob - oa), 4);
          g.globalAlpha = 1;
          if (mp[i] !== ".") {
            const mine = r.name === n1[i] && r.start === s1[i] ? seq1[i] : seq2[i];
            const bw = Math.max(2, ppb);
            for (const p of mp[i].split(",")) {
              const x = px(os[i] + +p);
              if (x < L - bw || x > W - R) continue;
              g.fillStyle = css(BASE_TOKEN[mine[+p]] ?? "--critical");
              g.fillRect(x, y - 2, bw, ROW + 4);
            }
          }
        }
        // bases of the focused pair
        if (showBases && focused && focus != null && seq1[focus] !== ".") {
          const isR1 = r.name === n1[focus] && r.start === s1[focus];
          const sq = isR1 ? seq1[focus] : seq2[focus], other = isR1 ? seq2[focus] : seq1[focus];
          g.font = `bold ${Math.min(12, ppb * 0.9)}px ${FONT}`; g.textAlign = "center"; g.textBaseline = "middle";
          for (let p = 0; p < sq.length; p++) {
            const x = px(os[focus] + p + 0.5);
            if (x < L || x > W - R) continue;
            g.fillStyle = sq[p] !== other[p] ? "#fff" : css("--ink");
            g.fillText(sq[p], x, y + ROW / 2 - 1);
          }
        }
        hit.push({ x0: a, x1: b, y, r });
      }
      if (laneRows[k] > MAX_ROWS) {
        g.fillStyle = css("--ink3"); g.font = `11px ${FONT}`; g.textAlign = "right"; g.textBaseline = "bottom";
        g.fillText(`+${laneRows[k] - MAX_ROWS} more rows`, W - R, y0 + laneH(k));
      }
      y0 += laneH(k) + 8;
      g.fillStyle = css("--line"); g.fillRect(L, y0 - 5, W - L - R, 1);
    }
    hits.current = hit;
    if (!reads.length) emptyMsg(g, W, H, "No overlapping reads here for this selection");
  }, [W, H, reads, win, theme, focus, lanes, ds, t, laneRows, cov, covH]); // eslint-disable-line react-hooks/exhaustive-deps

  // wheel zoom around the cursor, drag to pan
  const drag = useRef<{ x: number; win: Win } | null>(null);
  const clamp = (s: number, e: number): Win => {
    const span = Math.max(20, Math.min(MAX_SPAN * 5, e - s));
    let a = s, b = s + span;
    if (a < 0) { b -= a; a = 0; }
    if (b > chromLen) { a -= b - chromLen; b = chromLen; }
    return { chrom: win.chrom, start: Math.max(0, Math.round(a)), end: Math.round(b) };
  };
  useEffect(() => {
    const el = cv.current!;
    const onWheel = (e: WheelEvent) => {
      e.preventDefault();
      const x = e.clientX - el.getBoundingClientRect().left;
      const f = Math.exp(e.deltaY * 0.0015);
      const at = win.start + ((x - L) / (W - L - R)) * (win.end - win.start);
      setWin(clamp(at - (at - win.start) * f, at + (win.end - at) * f));
    };
    el.addEventListener("wheel", onWheel, { passive: false });
    return () => el.removeEventListener("wheel", onWheel);
  });

  const zoom = (f: number) => {
    const mid = (win.start + win.end) / 2, half = ((win.end - win.start) * f) / 2;
    setWin(clamp(mid - half, mid + half));
  };
  const span = win.end - win.start;
  return (
    <Card
      title="Reads"
      sub={<>{win.chrom}:{(win.start + 1).toLocaleString()}-{win.end.toLocaleString()} ({span.toLocaleString()} bp) · {sel.barcode ?? "all cells (pick one for coverage)"} · scroll to zoom, drag to move{span > 300 ? ", zoom in to see bases" : ""}</>}
      right={<>
        <HapKey hap={sel.hap} />
        <button className="iconbtn" onClick={() => zoom(0.5)} title="Zoom in">＋</button>
        <button className="iconbtn" onClick={() => zoom(2)} title="Zoom out">－</button>
        {focus != null && <button className="iconbtn" onClick={() => setWin(windowFor(t, focus, true))}>Fit overlap</button>}
        {focus != null && <button className="iconbtn" onClick={() => setWin(windowFor(t, focus, false))}>Fit reads</button>}
        <DownloadMenu onSave={(f) => saveCanvas([cv.current], "region", f)} />
      </>}
    >
      <div ref={wrap}>
        <canvas
          ref={cv}
          style={{ cursor: drag.current ? "grabbing" : "grab" }}
          onMouseDown={(e) => { drag.current = { x: e.clientX, win }; }}
          onMouseUp={(e) => {
            const d = drag.current;
            drag.current = null;
            if (d && Math.abs(e.clientX - d.x) < 3) {
              const r = e.currentTarget.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
              const h = hits.current.find((h) => x >= h.x0 && x <= h.x1 && y >= h.y && y <= h.y + ROW);
              if (h) useStore.getState().set({ focus: h.r.pairs[0] });
            }
          }}
          onMouseLeave={() => { drag.current = null; tip(null); }}
          onMouseMove={(e) => {
            const d = drag.current;
            if (d) {
              const dx = ((e.clientX - d.x) / (W - L - R)) * (d.win.end - d.win.start);
              setWin(clamp(d.win.start - dx, d.win.end - dx));
              return tip(null);
            }
            const r = e.currentTarget.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
            const h = hits.current.find((h) => x >= h.x0 && x <= h.x1 && y >= h.y && y <= h.y + ROW);
            tip(e, h ? `<b class="mono">${h.r.name}</b><br>${win.chrom}:${(h.r.start + 1).toLocaleString()}-${h.r.end.toLocaleString()} (${h.r.reverse ? "−" : "+"})<br><span class="k">overlaps</span> ${h.r.pairs.length} · click to inspect` : null);
          }}
        />
      </div>
    </Card>
  );
}

function PairDetail({ t, i }: { t: Table; i: number }) {
  const s = (c: string) => str(t, c)[i], n = (c: string) => num(t, c)[i], b = (c: string) => bool(t, c)[i] === 1;
  const ds = useStore((st) => st.ds)!;
  const junction = isJunction(t, i, ds.junctionSize, ds.tolerance, ds.maxMismatches);
  const mm = new Set(s("mismatch_positions") === "." ? [] : s("mismatch_positions").split(",").map(Number));
  const seq1 = s("seq1"), seq2 = s("seq2");
  return (
    <Card title="Read pair" sub={<>{junction ? <span className="badge junction">Tn5 junction</span> : <span className="badge">{s("overlap_type").replace(/_/g, "-")}</span>} <span className="mono">{s("barcode")}</span> · {HAP_LABEL[s("haplotype") as "1"] ?? s("haplotype")}{s("phase_block") !== "." ? ` · phase block ${s("phase_block")}` : ""}</>}>
      <div className="pills">
        <div className="pill"><span className="k">overlap</span><span className="v">{n("overlap_bp")} bp</span></div>
        <div className="pill"><span className="k">position</span><span className="v mono">{s("chrom")}:{(n("overlap_start") + 1).toLocaleString()}-{n("overlap_end").toLocaleString()}</span></div>
        <div className={"pill" + (n("n_mismatches") > ds.maxMismatches ? " bad" : "")}><span className="k">mismatches</span><span className="v">{n("n_mismatches")}</span></div>
        <div className={"pill" + (b("has_indel_in_overlap") ? " bad" : "")}><span className="k">indel in overlap</span><span className="v">{b("has_indel_in_overlap") ? "yes" : "no"}</span></div>
        <div className="pill"><span className="k">sequences match</span><span className="v">{b("seq_match") ? "yes" : "no"}</span></div>
      </div>
      <div className="tablewrap">
        <table>
          <thead><tr><th>read</th><th>name</th><th>start</th><th>end</th><th>strand</th><th>MAPQ</th></tr></thead>
          <tbody>
            {[1, 2].map((k) => (
              <tr key={k}><td>{k}</td><td className="mono">{s(`read${k}_name`)}</td><td>{(n(`read${k}_start`) + 1).toLocaleString()}</td><td>{n(`read${k}_end`).toLocaleString()}</td><td>{s(`read${k}_strand`)}</td><td>{n(`read${k}_mapq`)}</td></tr>
            ))}
          </tbody>
        </table>
      </div>
      {seq1 !== "." ? (
        <div className="seqblock mono">
          {[seq1, seq2].map((q, k) => (
            <div key={k}><span className="k">read {k + 1}</span>{q.split("").map((c, p) => <span key={p} className={mm.has(p) ? "mm" : ""}>{c}</span>)}</div>
          ))}
          <div><span className="k" />{seq1.split("").map((_, p) => <span key={p} className="mk">{mm.has(p) ? "^" : " "}</span>)}</div>
        </div>
      ) : (
        <div className="hint-keys">Sequences are written only for overlaps up to rrecon's --max-seq-bp.</div>
      )}
    </Card>
  );
}
