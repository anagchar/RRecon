// What was loaded: the run's files, parameters, read filtering and skipped groups.
import { useStore } from "../state/store";
import { Card } from "./widgets";

const STAT_LABEL: Record<string, string> = {
  total: "reads scanned", unmapped: "unmapped", secondary_or_supplementary: "secondary or supplementary",
  duplicate: "duplicates", low_mapq: "below --min-mapq", missing_barcode: "no barcode tag",
  other_barcode: "not in the barcode list", no_alignment_span: "no aligned bases", kept: "kept",
};

export function DatasetTab() {
  const ds = useStore((s) => s.ds)!;
  const m = ds.meta;
  return (
    <>
      <div className="tiles">
        <div className="tile"><div className="k">Cells</div><div className="v">{ds.cells.length.toLocaleString()}</div></div>
        <div className="tile"><div className="k">Chromosomes</div><div className="v">{ds.genome.chroms.length}</div></div>
        <div className="tile"><div className="k">Read pairs loaded</div><div className="v">{(ds.tables.overlaps?.n ?? ds.tables.tn5_junctions?.n ?? 0).toLocaleString()}</div></div>
        <div className="tile"><div className="k">Junction rule</div><div className="v long">{ds.junctionSize}{ds.tolerance ? ` ± ${ds.tolerance}` : ""} bp, ≤ {ds.maxMismatches} mm</div></div>
      </div>
      <div className="grid2">
        <Card title="Files" sub={ds.name}>
          <table>
            <tbody>{ds.files.map((f) => <tr key={f}><td className="mono">{f}</td></tr>)}</tbody>
          </table>
          {!m && <p className="hint-keys">No meta.json was loaded, so chromosome lengths are estimated from the data and the junction rule is the default.</p>}
        </Card>
        {m?.stats && (
          <Card title="Read filtering" sub={m.bam ? `from ${m.bam}` : undefined}>
            <table>
              <tbody>
                {Object.entries(m.stats).map(([k, v]) => (
                  <tr key={k}><td>{STAT_LABEL[k] ?? k}</td><td>{v.toLocaleString()}</td><td>{m.stats!.total ? ((100 * v) / m.stats!.total).toFixed(1) + "%" : ""}</td></tr>
                ))}
              </tbody>
            </table>
          </Card>
        )}
      </div>
      {m && (
        <Card title="Run" sub={`rrecon ${m.version ?? ""} · ${m.engine ?? ""} engine${m.seconds != null ? ` · ${m.seconds}s` : ""}${m.created ? ` · ${m.created}` : ""}`}>
          <table>
            <tbody>
              {Object.entries(m.params ?? {}).map(([k, v]) => <tr key={k}><td>{k}</td><td className="mono">{v == null ? "–" : String(v)}</td></tr>)}
            </tbody>
          </table>
        </Card>
      )}
      {m?.skipped_groups?.length ? (
        <Card title="Skipped groups" sub="over --max-reads-per-group: no pairs were computed for these">
          <table>
            <thead><tr><th>cell</th><th>haplotype</th><th>chrom</th><th>reads</th></tr></thead>
            <tbody>{m.skipped_groups.map((g, k) => <tr key={k}><td className="mono">{g.barcode}</td><td>{g.haplotype}</td><td>{g.chrom}</td><td>{g.n_reads.toLocaleString()}</td></tr>)}</tbody>
          </table>
        </Card>
      ) : null}
    </>
  );
}
