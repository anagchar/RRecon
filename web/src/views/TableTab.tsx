// Sortable, filterable table of read pairs, with CSV export; a row opens the region view.
import { useMemo, useState } from "react";
import { HAP_LABEL, filterRows, isJunction, selectRows, sortRows, type RowFilter } from "../lib/stats";
import { bool, num, str, toDelimited } from "../lib/table";
import { pairs, useStore } from "../state/store";
import { useSelection } from "./OverviewTab";
import { Card, DlIcon, NumField, Seg, Switch, downloadText } from "./widgets";

const PAGE = 100;
const COLS: [string, string][] = [
  ["barcode", "cell"], ["haplotype", "hap"], ["phase_block", "phase block"], ["chrom", "chrom"],
  ["overlap_start", "start"], ["overlap_end", "end"], ["overlap_bp", "size"], ["overlap_type", "type"],
  ["n_mismatches", "mm"], ["mismatch_positions", "mm at"], ["has_indel_in_overlap", "indel"],
  ["read1_name", "read 1"], ["read2_name", "read 2"], ["seq1", "seq 1"], ["seq2", "seq 2"],
];

export function TableTab() {
  const { ds, sel } = useSelection();
  const t = pairs(ds);
  const f = useStore((s) => s.rowFilter);
  const sort = useStore((s) => s.sort);
  const focus = useStore((s) => s.focus);
  const [page, setPage] = useState(0);
  const setF = (p: Partial<RowFilter>) => { useStore.getState().set({ rowFilter: { ...f, ...p } }); setPage(0); };

  const rows = useMemo(() => {
    if (!t) return new Int32Array(0);
    const base = selectRows(t, ds.genome, sel, ["overlap_start", "overlap_end"]);
    return sortRows(t, filterRows(t, base, f, ds.junctionSize, ds.tolerance, ds.maxMismatches), sort.col, sort.desc);
  }, [t, ds, sel.barcode, sel.hap, sel.region, f, sort]); // eslint-disable-line react-hooks/exhaustive-deps

  if (!t) return <div className="card"><div className="empty">Load an overlaps or tn5_junctions file to list read pairs.</div></div>;
  const pages = Math.max(1, Math.ceil(rows.length / PAGE));
  const shown = rows.subarray(page * PAGE, page * PAGE + PAGE);
  const cell = (c: string, i: number) => {
    if (c === "has_indel_in_overlap") return bool(t, c)[i] ? "yes" : "";
    if (c === "haplotype") return HAP_LABEL[str(t, c)[i] as "1"] ?? str(t, c)[i];
    if (c === "overlap_start") return (num(t, c)[i] + 1).toLocaleString();
    const col = t.cols[c];
    return col instanceof Float64Array ? col[i].toLocaleString() : String(col[i]);
  };
  const exportCsv = () => {
    const safe = (ds.name || "rrecon").replace(/[^\w.-]+/g, "_");
    downloadText(toDelimited(t, rows, ","), `${safe}_${sel.barcode ?? "all"}_${sel.hap}_pairs.csv`);
  };
  return (
    <Card
      title="Read pairs"
      sub={<>{rows.length.toLocaleString()} rows · {sel.barcode ?? "all cells"} · {sel.hap === "all" ? "all haplotypes" : HAP_LABEL[sel.hap]}{!ds.tables.overlaps ? " · junctions only" : ""}</>}
      right={<button className="iconbtn" onClick={exportCsv} disabled={!rows.length}><DlIcon />CSV</button>}
    >
      <div className="filterbar">
        <label>size <NumField label="Smallest overlap" value={f.minBp} onChange={(v) => setF({ minBp: v })} placeholder="min" /> – <NumField label="Largest overlap" value={f.maxBp} onChange={(v) => setF({ maxBp: v })} placeholder="max" suffix="bp" /></label>
        <label>mismatches ≤ <NumField label="Most mismatches" value={f.maxMm} onChange={(v) => setF({ maxMm: v })} /></label>
        <Seg label="Overlap type" value={f.type} onChange={(v) => setF({ type: v })} options={[["all", "all types"], ["end_to_start", "end-to-start"], ["contained", "contained"]]} />
        <Seg label="Indels" value={f.indel} onChange={(v) => setF({ indel: v })} options={[["all", "any"], ["without", "no indel"], ["with", "indel"]]} />
        <label>phase block <input className="tinput" value={f.phase} placeholder="any" onChange={(e) => setF({ phase: e.target.value })} /></label>
        <label className="sw">junctions only <Switch label="Junctions only" on={f.junctionsOnly} onChange={(v) => setF({ junctionsOnly: v })} /></label>
      </div>
      <div className="tablewrap pairs">
        <table>
          <thead>
            <tr>
              {COLS.map(([c, label]) => (
                <th key={c} className={"sortable" + (sort.col === c ? " on" : "")} onClick={() => useStore.getState().set({ sort: { col: c, desc: sort.col === c ? !sort.desc : false } })}>
                  {label}{sort.col === c ? (sort.desc ? " ↓" : " ↑") : ""}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {Array.from(shown, (i) => (
              <tr key={i} className={"clickable" + (i === focus ? " on" : "") + (isJunction(t, i, ds.junctionSize, ds.tolerance, ds.maxMismatches) ? " junction" : "")} onClick={() => useStore.getState().focusRow(i)} title="Show these reads">
                {COLS.map(([c]) => <td key={c} className={c.startsWith("seq") || c.endsWith("name") || c === "barcode" ? "mono" : ""}>{cell(c, i)}</td>)}
              </tr>
            ))}
          </tbody>
        </table>
        {!rows.length && <div className="empty">No read pairs match.</div>}
      </div>
      {pages > 1 && (
        <div className="pager">
          <button className="iconbtn" disabled={page === 0} onClick={() => setPage(page - 1)}>‹ Previous</button>
          <span>page {page + 1} of {pages.toLocaleString()}</span>
          <button className="iconbtn" disabled={page >= pages - 1} onClick={() => setPage(page + 1)}>Next ›</button>
        </div>
      )}
    </Card>
  );
}
