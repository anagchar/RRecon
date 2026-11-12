// The overlapscope mark: two reads drawn as cnvscope draws chromosomes (background fill,
// ink outline), offset so they overlap, with the shared 9 bp stretch in the junction colour.
// Colours: --logo-* and --junction in theme.css.

const H = 10, SW = 1.8;

export function LogoIcon({ size = 32, title = "overlapscope" }: { size?: number; title?: string }) {
  return (
    <svg className="logo-icon" height={size} width={size * 1.25} viewBox="0 0 80 64" role="img" aria-label={title}>
      <rect x="4" y="16" width="46" height={H} rx={H / 2} fill="var(--logo-fill)" stroke="var(--logo-line)" strokeWidth={SW} />
      <rect x="30" y="38" width="46" height={H} rx={H / 2} fill="var(--logo-fill)" stroke="var(--logo-line)" strokeWidth={SW} />
      <rect x="31.9" y="17.9" width="16.2" height={H - 3.8} rx="1.6" fill="var(--junction)" />
      <rect x="31.9" y="39.9" width="16.2" height={H - 3.8} rx="1.6" fill="var(--junction)" />
      <path d="M30 30v4M50 30v4" stroke="var(--logo-line)" strokeWidth="1.4" strokeLinecap="round" opacity=".6" />
    </svg>
  );
}
