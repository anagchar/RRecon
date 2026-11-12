// Small shared widgets (the size hooks, icons and menus are cnvscope's).
import { useEffect, useRef, useState, type CSSProperties, type ReactNode, type RefObject } from "react";
import { HAP_COLOR, HAP_LABEL, HAPS, type HapSel } from "../lib/stats";
import { useStore } from "../state/store";

/** Width of an element, followed as it resizes. */
export function useWidth(ref: RefObject<HTMLElement | null>, fallback = 600): number {
  const [w, setW] = useState(fallback);
  useEffect(() => {
    if (!ref.current) return;
    const ro = new ResizeObserver((e) => {
      const v = e[0]?.contentRect.width;
      if (v > 0) setW(Math.round(v));
    });
    ro.observe(ref.current);
    return () => ro.disconnect();
  }, [ref]);
  return w;
}

export function ThemeIcon() {
  return (
    <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5">
      <circle cx="8" cy="8" r="6" />
      <path d="M8 2a6 6 0 0 0 0 12z" fill="currentColor" />
    </svg>
  );
}

export function DlIcon() {
  return (
    <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
      <path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5 5 5-5M12 15V3" />
    </svg>
  );
}

export function Switch({ on, onChange, label }: { on: boolean; onChange: (v: boolean) => void; label: string }) {
  return <button className="switch" role="switch" aria-checked={on} aria-label={label} onClick={() => onChange(!on)} />;
}

export function DownloadMenu({ onSave }: { onSave: (f: "png" | "pdf") => void }) {
  const [open, setOpen] = useState(false);
  const ref = useRef<HTMLDivElement>(null);
  useEffect(() => {
    if (!open) return;
    const h = (e: MouseEvent) => {
      if (!ref.current?.contains(e.target as Node)) setOpen(false);
    };
    document.addEventListener("mousedown", h);
    return () => document.removeEventListener("mousedown", h);
  }, [open]);
  return (
    <div className="menu" ref={ref}>
      <button className="iconbtn" onClick={() => setOpen((o) => !o)} title="Save the figure" aria-haspopup="menu" aria-expanded={open}>
        <DlIcon />
        Save
      </button>
      {open && (
        <div className="pop" role="menu">
          {(["png", "pdf"] as const).map((f) => (
            <button key={f} role="menuitem" onClick={() => { onSave(f); setOpen(false); }}>
              {f.toUpperCase()}
            </button>
          ))}
        </div>
      )}
    </div>
  );
}

export function Seg<T extends string | number>({ value, options, onChange, label }: {
  value: T; options: [T, ReactNode][]; onChange: (v: T) => void; label: string;
}) {
  return (
    <div className="seg" role="radiogroup" aria-label={label}>
      {options.map(([v, l]) => (
        <button key={String(v)} className={v === value ? "on" : ""} role="radio" aria-checked={v === value} onClick={() => onChange(v)}>{l}</button>
      ))}
    </div>
  );
}

export function NumField({ value, onChange, placeholder, suffix, label }: {
  value: number | null; onChange: (v: number | null) => void; placeholder?: string; suffix?: string; label: string;
}) {
  return (
    <label className="numfield" title={label}>
      <input
        type="number" min={0} aria-label={label} value={value ?? ""} placeholder={placeholder ?? "–"}
        onChange={(e) => onChange(e.target.value === "" ? null : Math.max(0, +e.target.value))}
      />
      {suffix && <span>{suffix}</span>}
    </label>
  );
}

/** Haplotype key: coloured swatches for the haplotypes shown. */
export function HapKey({ hap }: { hap: HapSel }) {
  return (
    <span className="hapkey">
      {HAPS.filter((h) => hap === "all" || h === hap).map((h) => (
        <span key={h}><i style={{ background: `var(${HAP_COLOR[h]})` } as CSSProperties} />{HAP_LABEL[h]}</span>
      ))}
    </span>
  );
}

/** Redraw when the theme changes (colours are read from CSS tokens). */
export const useThemeKey = () => useStore((s) => s.theme);

export function Card({ title, sub, right, children, id }: { title: string; sub?: ReactNode; right?: ReactNode; children: ReactNode; id?: string }) {
  return (
    <div className="card" id={id}>
      <div className="hd">
        <h2>{title}</h2>
        {sub != null && <span className="sub">{sub}</span>}
        {right}
      </div>
      <div className="bd">{children}</div>
    </div>
  );
}

export const TAB_ICONS = {
  overview: (
    <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
      <path d="M1.5 12.5v-3M4.5 12.5v-6M7.5 12.5v-4M10.5 12.5v-8M13.5 12.5v-5" />
      <path d="M1 14.5h14" opacity=".45" />
    </svg>
  ),
  tn5: (
    <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" aria-hidden="true">
      <path d="M1.5 6h8M6.5 10h8" />
      <rect x="6.5" y="4.5" width="3" height="7" rx="1" opacity=".5" />
    </svg>
  ),
  region: (
    <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" aria-hidden="true">
      <circle cx="7" cy="7" r="4.6" /><path d="M10.5 10.5 14 14M5 7h4M7 5v4" />
    </svg>
  ),
  table: (
    <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5" aria-hidden="true">
      <rect x="1.5" y="2.5" width="13" height="11" rx="2" /><path d="M1.5 6.5h13M1.5 10h13M6 6.5v7" />
    </svg>
  ),
  data: (
    <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
      <rect x="2.5" y="1.5" width="11" height="13" rx="2" />
      <path d="M5.5 5.5h5M5.5 8.5h5M5.5 11.5h3" />
    </svg>
  ),
};

export function downloadText(text: string, name: string, type = "text/csv") {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], { type }));
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
