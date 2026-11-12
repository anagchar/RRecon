import { forwardRef, useEffect, useMemo, useState } from "react";
import { sum } from "../lib/dataset";
import { formatRegion, parseRegion, regionOf } from "../lib/genome";
import { HAP_COLOR, HAP_LABEL, hapIndex, HAPS, hapsShown, type HapSel } from "../lib/stats";
import { useStore } from "../state/store";
import { fmtCount } from "../lib/canvas";

const LIST_CAP = 400;

/** Step through the cell list (keyboard ↑/↓). */
export function navigate(dir: number) {
  const { barcode, set } = useStore.getState();
  if (!currentList.length) return;
  const i = barcode == null ? -1 : currentList.indexOf(barcode);
  set({ barcode: currentList[Math.max(0, Math.min(currentList.length - 1, i < 0 ? 0 : i + dir))] });
}
let currentList: string[] = [];

export const Sidebar = forwardRef<HTMLInputElement>(function Sidebar(_, searchRef) {
  const ds = useStore((s) => s.ds)!;
  const s = useStore();
  const ks = hapsShown(s.hap).map(hapIndex);
  const list = useMemo(() => {
    const q = s.search.trim().toLowerCase();
    return q ? ds.cells.filter((c) => c.barcode.toLowerCase().includes(q)) : ds.cells;
  }, [ds, s.search]);
  currentList = list.map((c) => c.barcode);
  useEffect(() => { document.querySelector(".cell.on")?.scrollIntoView({ block: "nearest" }); }, [s.barcode]);

  const [regionText, setRegionText] = useState(formatRegion(s.region));
  const [bad, setBad] = useState(false);
  useEffect(() => { setRegionText(formatRegion(s.region)); setBad(false); }, [s.region]);
  const applyRegion = () => {
    const r = parseRegion(ds.genome, regionText);
    if (r) s.set({ region: r });
    setBad(!r);
  };
  const hapCount = (h: HapSel) => {
    const cells = s.barcode == null ? ds.cells : ds.cells.filter((c) => c.barcode === s.barcode);
    const k = h === "all" ? [0, 1, 2] : [hapIndex(h)];
    return cells.reduce((a, c) => a + k.reduce((x, i) => x + c.reads[i], 0), 0);
  };

  return (
    <aside><div className="aside-scroll">
      <div className="group">
        <h3>Haplotype</h3>
        <div className="seg wrap" role="radiogroup" aria-label="Haplotype">
          {(["all", ...HAPS] as HapSel[]).map((h) => (
            <button key={h} role="radio" aria-checked={s.hap === h} className={s.hap === h ? "on" : ""} onClick={() => s.set({ hap: h })} title={`${hapCount(h).toLocaleString()} reads`}>
              {h !== "all" && <i className="hapdot" style={{ background: `var(${HAP_COLOR[h]})` }} />}
              {h === "all" ? "All" : HAP_LABEL[h]}
            </button>
          ))}
        </div>
        {ds.cells.length > 0 && ds.cells.every((c) => !c.reads[0] && !c.reads[1]) && (
          <small className="note">No HP-tagged reads: run rrecon with <code>--hp-tag HP</code> on a whatshap-haplotagged BAM.</small>
        )}
      </div>

      <div className="group">
        <h3>Region</h3>
        <select className="chromsel" value={s.region.chrom ?? ""} onChange={(e) => s.set({ region: regionOf(ds.genome, e.target.value || null) })} aria-label="Chromosome">
          <option value="">Whole genome</option>
          {ds.genome.chroms.map((c) => <option key={c.name} value={c.name}>{c.name}</option>)}
        </select>
        <input
          className={"search region" + (bad ? " bad" : "")}
          value={regionText}
          spellCheck={false}
          onChange={(e) => setRegionText(e.target.value)}
          onKeyDown={(e) => { if (e.key === "Enter") applyRegion(); }}
          onBlur={applyRegion}
          aria-label="Region, e.g. chr1:1,000,000-2,000,000"
        />
        <small className="note">e.g. <code>chr1:1.2M-1.5M</code>; drag on a track to zoom, double-click for the genome</small>
      </div>

      <div className="group">
        <h3>Cell <span className="count">{ds.cells.length.toLocaleString()}</span></h3>
        <input ref={searchRef} className="search" placeholder="Search barcodes  /" value={s.search} onChange={(e) => s.set({ search: e.target.value })} aria-label="Search barcodes" />
        <div className="cells" role="listbox" aria-label="Cells">
          <button className={"cell" + (s.barcode == null ? " on" : "")} role="option" aria-selected={s.barcode == null} onClick={() => s.set({ barcode: null })}>
            <span className="nm">All cells</span>
          </button>
          {list.slice(0, LIST_CAP).map((c) => {
            const j = ks.reduce((a, k) => a + c.junctions[k], 0), o = ks.reduce((a, k) => a + c.overlaps[k], 0);
            return (
              <button key={c.barcode} className={"cell" + (s.barcode === c.barcode ? " on" : "")} role="option" aria-selected={s.barcode === c.barcode} onClick={() => s.set({ barcode: c.barcode })}
                title={`${sum(c.reads).toLocaleString()} reads · ${o.toLocaleString()} overlaps · ${j.toLocaleString()} junctions`}>
                <span className="nm">{c.barcode}</span>
                <span className="pl" title="junctions / overlaps">{fmtCount(j)}<span className="ct"> / {fmtCount(o)}</span></span>
              </button>
            );
          })}
          {list.length > LIST_CAP && <div className="more">{(list.length - LIST_CAP).toLocaleString()} more: search to find them</div>}
          {!list.length && <div className="more">No cell matches</div>}
        </div>
        <small className="note">sorted by junctions · numbers: junctions / overlaps · <kbd>↑</kbd> <kbd>↓</kbd> to step</small>
      </div>

      <div className="group">
        <h3>Key</h3>
        <div className="keylist">
          {HAPS.map((h) => <span key={h}><i style={{ background: `var(${HAP_COLOR[h]})` }} />{HAP_LABEL[h]}</span>)}
          <span><i className="round" style={{ background: "var(--endstart)" }} />end-to-start overlap</span>
          <span><i className="round" style={{ background: "var(--contained)" }} />contained read</span>
          <span><i className="round" style={{ background: "var(--junction)" }} />{ds.junctionSize} bp Tn5 junction</span>
          <span><i style={{ background: "var(--critical)", width: 3 }} />mismatch</span>
        </div>
      </div>
      <div className="hint-keys"><kbd>/</kbd> search · <kbd>t</kbd> theme · <kbd>1</kbd>–<kbd>5</kbd> tabs · <kbd>Esc</kbd> genome</div>
    </div></aside>
  );
});
