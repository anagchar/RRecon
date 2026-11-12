import { useEffect, useRef } from "react";
import { sum } from "./lib/dataset";
import { useStore, type Tab } from "./state/store";
import { DatasetTab } from "./views/DatasetTab";
import { LogoIcon } from "./views/Logo";
import { OverviewTab } from "./views/OverviewTab";
import { RegionTab } from "./views/RegionView";
import { navigate, Sidebar } from "./views/Sidebar";
import { TableTab } from "./views/TableTab";
import { Tn5Tab } from "./views/Tn5Tab";
import { UploadScreen } from "./views/UploadScreen";
import { TAB_ICONS, ThemeIcon } from "./views/widgets";

const TABS: [Tab, string, string][] = [
  ["overview", "Overview", "Density, overlap locations and sizes"],
  ["tn5", "Tn5 junctions", "Junctions per cell, 9-mers and mismatches"],
  ["region", "Region", "The reads behind an overlap"],
  ["table", "Table", "Every read pair, sortable and filterable"],
  ["data", "Dataset", "What was loaded"],
];

export function App() {
  const ds = useStore((s) => s.ds);
  const booting = useStore((s) => s.booting);
  const picking = useStore((s) => s.picking);
  const loading = useStore((s) => s.loading);
  const boot = useStore((s) => s.boot);
  useEffect(() => { boot(); }, [boot]);
  if (booting && loading) return <div className="upload"><div className="inner"><div className="spinner" /><div className="empty">{loading}</div></div></div>;
  if (booting) return null;
  if (!ds || picking) return <UploadScreen />;
  return <Explorer key={ds.name + ds.files.join()} />;
}

function Explorer() {
  const s = useStore();
  const ds = s.ds!;
  const searchRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    const h = (e: KeyboardEvent) => {
      const t = e.target as HTMLElement;
      if (t.tagName === "INPUT" || t.tagName === "SELECT") { if (e.key === "Escape") t.blur(); return; }
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      const st = useStore.getState();
      if (e.key === "ArrowUp" || e.key === "k") { e.preventDefault(); navigate(-1); }
      else if (e.key === "ArrowDown" || e.key === "j") { e.preventDefault(); navigate(1); }
      else if (e.key === "Escape") st.set({ region: { chrom: null, start: 0, end: st.ds!.genome.length } });
      else if (e.key === "/") { e.preventDefault(); st.set({ panel: true }); setTimeout(() => searchRef.current?.focus(), 0); }
      else if (e.key === "t") st.toggleTheme();
      else if (/^[1-5]$/.test(e.key)) st.set({ tab: TABS[+e.key - 1][0] });
    };
    addEventListener("keydown", h);
    return () => removeEventListener("keydown", h);
  }, []);

  const nPairs = ds.tables.overlaps?.n;
  const nJ = ds.cells.reduce((a, c) => a + sum(c.junctions), 0);
  return (
    <div className={"app" + (s.panel ? "" : " nopanel")}>
      <header>
        <div className="brand"><LogoIcon size={30} />overlapscope</div>
        <div className="project" title={ds.name}>{ds.name}</div>
        <div className="chips">
          <span className="chip">{ds.cells.length.toLocaleString()} cells</span>
          {nPairs != null && <span className="chip">{nPairs.toLocaleString()} overlaps</span>}
          <span className="chip jn">{nJ.toLocaleString()} junctions</span>
          {ds.meta?.engine && <span className="chip">rrecon {ds.meta.version} · {ds.meta.engine}</span>}
        </div>
        <div className="spacer" />
        <button className="iconbtn" onClick={() => s.set({ panel: !s.panel })} title={s.panel ? "Hide the side panel" : "Show the side panel"} aria-pressed={s.panel}>
          <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5"><rect x="1.5" y="2.5" width="13" height="11" rx="2" /><path d="M6 2.5v11" /></svg><span className="btnlabel">Panel</span>
        </button>
        <button className="iconbtn" onClick={() => s.set({ picking: true, error: null })} title="Open other files">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5-5 5 5M12 5v12" /></svg><span className="btnlabel">Open</span>
        </button>
        <button className="iconbtn" onClick={s.toggleTheme} title="Light or dark (t)" aria-label="Switch between light and dark"><ThemeIcon /></button>
      </header>

      <Sidebar ref={searchRef} />

      <main>
        <div className="tabs" role="tablist">
          {TABS.map(([id, label, title]) => (
            <button key={id} className={"tab" + (s.tab === id ? " on" : "")} role="tab" aria-selected={s.tab === id} title={title} onClick={() => s.set({ tab: id })}>{TAB_ICONS[id]}{label}</button>
          ))}
        </div>
        {s.error && <div className="error">{s.error}</div>}
        <section>
          {s.tab === "overview" && <OverviewTab />}
          {s.tab === "tn5" && <Tn5Tab />}
          {s.tab === "region" && <RegionTab />}
          {s.tab === "table" && <TableTab />}
          {s.tab === "data" && <DatasetTab />}
        </section>
      </main>
    </div>
  );
}
