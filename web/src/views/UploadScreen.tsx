import { useRef, useState } from "react";
import { isBrowserOnly, useStore } from "../state/store";
import { LogoIcon } from "./Logo";
import { ThemeIcon } from "./widgets";

export function UploadScreen() {
  const { openFiles, demo, loading, error, ds, toggleTheme, set } = useStore();
  const [over, setOver] = useState(false);
  const input = useRef<HTMLInputElement>(null);
  const take = (list?: FileList | null) => { if (list?.length) openFiles([...list]); };
  return (
    <div className="upload">
      <button className="iconbtn" style={{ position: "fixed", top: 16, right: 16 }} onClick={toggleTheme} title="Light or dark (t)" aria-label="Switch between light and dark"><ThemeIcon /></button>
      <div className="inner">
        <div className="hero">
          <div className="mark"><LogoIcon size={110} /></div>
          <h1>overlapscope</h1>
          <p>Read overlaps and Tn5 junctions, per cell and haplotype</p>
        </div>
        <button
          className={"card drop" + (over ? " over" : "")}
          onDragOver={(e) => { e.preventDefault(); setOver(true); }}
          onDragLeave={() => setOver(false)}
          onDrop={(e) => { e.preventDefault(); setOver(false); take(e.dataTransfer.files); }}
          onClick={() => input.current?.click()}
        >
          <input ref={input} type="file" multiple accept=".parquet,.gz,.json" hidden onChange={(e) => take(e.target.files)} />
          {loading ? <><div className="spinner" /><div className="small">{loading}</div></> : <>
            <svg width="34" height="34" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.6"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5-5 5 5M12 5v12" /></svg>
            <div className="big">Drop the files from an <code>rrecon</code> run</div>
            <div className="small">all of them at once: <code>*.parquet</code> or <code>*.tsv.gz</code>, and <code>*.meta.json</code></div>
            <div className="small">Files are read in your browser; nothing is uploaded.</div>
          </>}
        </button>
        {error && <div className="err">{error}</div>}
        <div className="actions">
          {isBrowserOnly() && <button className="iconbtn" onClick={demo} disabled={!!loading}>Open the demo</button>}
          {ds && <button className="iconbtn" onClick={() => set({ picking: false, error: null })}>Back to {ds.name}</button>}
        </div>
        <div className="card how">
          <pre>{`pip install git+https://github.com/anagchar/RRecon
rrecon all cells.haplotagged.bam -o results/ \\
    --hp-tag HP --ps-tag PS --threads 8
rrecon gui --dir results/`}</pre>
        </div>
      </div>
    </div>
  );
}
