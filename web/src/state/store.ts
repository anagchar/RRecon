import { create } from "zustand";
import { loadDataset, type Dataset, type InFile } from "../lib/dataset";
import { regionOf, type Region } from "../lib/genome";
import { defaultRowFilter, type HapSel, type RowFilter } from "../lib/stats";
import type { Table } from "../lib/table";

export type Theme = "light" | "dark";
export type Tab = "overview" | "tn5" | "region" | "table" | "data";

/** Above this, warn: the browser holds every row in memory. */
export const BIG_FILE = 400 * 1024 * 1024;

export const isBrowserOnly = () => import.meta.env.VITE_STATIC === "1";

interface State {
  ds: Dataset | null;
  booting: boolean;
  loading: string | null;
  picking: boolean;
  error: string | null;
  serverFolder: string | null;

  barcode: string | null;
  hap: HapSel;
  region: Region;
  tab: Tab;
  panel: boolean;
  theme: Theme;
  search: string;
  /** row of pairs() shown in the region view */
  focus: number | null;
  rowFilter: RowFilter;
  sort: { col: string; desc: boolean };
  histMax: number;

  boot: () => Promise<void>;
  openFiles: (files: File[]) => Promise<void>;
  demo: () => Promise<void>;
  setDataset: (ds: Dataset) => void;
  set: (p: Partial<State>) => void;
  toggleTheme: () => void;
  focusRow: (i: number) => void;
}

const initialTheme = (): Theme => (document.documentElement.dataset.theme === "dark" ? "dark" : "light");

/** The read-pair table the views draw: all overlaps, or only junctions if that is what was loaded. */
export function pairs(ds: Dataset): Table | undefined {
  return ds.tables.overlaps ?? ds.tables.tn5_junctions;
}

async function fetchFiles(base: string, names: string[], progress: (s: string) => void): Promise<InFile[]> {
  const out: InFile[] = [];
  for (const name of names) {
    progress(`Fetching ${name}…`);
    const r = await fetch(base + encodeURIComponent(name));
    if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`);
    out.push({ name, data: await r.arrayBuffer() });
  }
  return out;
}

export const useStore = create<State>((set, get) => ({
  ds: null,
  booting: true,
  loading: null,
  picking: false,
  error: null,
  serverFolder: null,
  barcode: null,
  hap: "all",
  region: { chrom: null, start: 0, end: 1 },
  tab: "overview",
  panel: true,
  theme: initialTheme(),
  search: "",
  focus: null,
  rowFilter: defaultRowFilter(),
  sort: { col: "overlap_bp", desc: false },
  histMax: 60,

  boot: async () => {
    if (isBrowserOnly()) {
      set({ booting: false });
      return;
    }
    try {
      const r = await fetch("api/files");
      if (r.ok) {
        const { folder, files } = (await r.json()) as { folder: string | null; files: { name: string; size: number }[] };
        set({ serverFolder: folder });
        if (files.length) {
          set({ loading: "Loading…" });
          const big = files.filter((f) => f.size > BIG_FILE).map((f) => f.name);
          if (big.length) console.warn("overlapscope: large files, loading may be slow:", big);
          const all = await fetchFiles("api/file/", files.map((f) => f.name), (s) => set({ loading: s }));
          set({ loading: "Reading tables…" });
          get().setDataset(await loadDataset(all));
        }
      }
    } catch (e) {
      set({ error: e instanceof Error ? e.message : String(e) });
    } finally {
      set({ booting: false, loading: null });
    }
  },

  openFiles: async (files) => {
    const big = files.filter((f) => f.size > BIG_FILE);
    set({
      loading: "Reading files…",
      error: big.length
        ? `Large file${big.length > 1 ? "s" : ""} (${big.map((f) => f.name).join(", ")}): the browser keeps every row in memory, so this can be slow. Run rrecon with --barcodes or on fewer contigs to make smaller outputs.`
        : null,
    });
    try {
      const inputs = await Promise.all(files.map(async (f) => ({ name: f.name, data: await f.arrayBuffer() })));
      set({ loading: "Reading tables…" });
      get().setDataset(await loadDataset(inputs));
    } catch (e) {
      set({ error: e instanceof Error ? e.message : String(e) });
    } finally {
      set({ loading: null });
    }
  },

  demo: async () => {
    set({ loading: "Fetching the demo…", error: null });
    try {
      const r = await fetch("demo/index.json");
      if (!r.ok) throw new Error("the demo is not part of this build (npm run demo)");
      const names = (await r.json()) as string[];
      const files = await fetchFiles("demo/", names, (s) => set({ loading: s }));
      get().setDataset(await loadDataset(files));
    } catch (e) {
      set({ error: e instanceof Error ? e.message : String(e) });
    } finally {
      set({ loading: null });
    }
  },

  setDataset: (ds) =>
    set({
      ds,
      picking: false,
      error: null,
      barcode: ds.cells[0]?.barcode ?? null,
      hap: "all",
      region: regionOf(ds.genome, null),
      tab: "overview",
      search: "",
      focus: null,
      rowFilter: defaultRowFilter(),
      sort: { col: "overlap_bp", desc: false },
      histMax: Math.max(30, Math.min(200, (ds.junctionSize + ds.tolerance) * 6)),
    }),

  set: (p) => set(p),

  toggleTheme: () => {
    const theme: Theme = get().theme === "dark" ? "light" : "dark";
    document.documentElement.dataset.theme = theme;
    try {
      localStorage.setItem("overlapscope-theme", theme);
    } catch {
      /* storage blocked */
    }
    set({ theme });
  },

  focusRow: (i) => set({ focus: i, tab: "region" }),
}));
