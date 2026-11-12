import { rmSync } from "node:fs";
import { resolve } from "node:path";
import { defineConfig, type Plugin } from "vitest/config";
import react from "@vitejs/plugin-react";

// The built UI ships inside the Python wheel; `rrecon gui` serves it.
// `vite build --mode pages` (npm run build:pages) makes the GitHub Pages site
// instead: no server, files are read in the browser (.env.pages: VITE_STATIC=1),
// with the synthetic demo from `npm run demo` (public/demo).
const outDir = (mode: string) => (mode === "pages" ? "dist-pages" : "../rrecon/static");

/** Keep the demo data out of the wheel: only the Pages site serves it. */
const dropDemo = (mode: string): Plugin => ({
  name: "drop-demo",
  apply: "build",
  closeBundle() {
    if (mode !== "pages") rmSync(resolve(__dirname, outDir(mode), "demo"), { recursive: true, force: true });
  },
});

export default defineConfig(({ mode }) => ({
  plugins: [react(), dropDemo(mode)],
  base: "./",
  build: { outDir: outDir(mode), emptyOutDir: true },
  server: { proxy: { "/api": "http://127.0.0.1:8766" } },
  test: { environment: "node" },
}));
