"""Serve overlapscope (the built web UI in rrecon/static) and one results folder.

- ``GET /api/files``        the folder's rrecon outputs: {"files": [{name, size}]}
                            ({"files": []} with no folder; the UI then shows its drop screen)
- ``GET /api/file/{name}``  one of those files, as is
"""

from __future__ import annotations

import socket
import threading
import webbrowser
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles

from . import __version__

STATIC = Path(__file__).parent / "static"
SUFFIXES = (".tsv.gz", ".parquet", ".meta.json")


def output_files(folder: Path | None) -> list[Path]:
    if folder is None:
        return []
    return sorted(p for p in folder.iterdir() if p.is_file() and p.name.endswith(SUFFIXES))


def create_app(folder: Path | None = None) -> FastAPI:
    app = FastAPI(title="overlapscope", version=__version__)
    folder = folder.resolve() if folder else None

    @app.get("/api/files")
    def files() -> dict:
        return {
            "folder": folder.name if folder else None,
            "files": [{"name": p.name, "size": p.stat().st_size} for p in output_files(folder)],
        }

    @app.get("/api/file/{name}")
    def file(name: str):
        allowed = {p.name: p for p in output_files(folder)}
        if name not in allowed:
            raise HTTPException(404, "no such output file")
        return FileResponse(allowed[name])

    @app.middleware("http")
    async def revalidate_page(request, call_next):
        response = await call_next(request)
        if not request.url.path.startswith(("/assets/", "/api/")):
            response.headers["Cache-Control"] = "no-cache"
        return response

    if (STATIC / "index.html").exists():
        app.mount("/", StaticFiles(directory=STATIC, html=True), name="ui")
    else:

        @app.get("/", response_class=HTMLResponse)
        def no_ui() -> str:
            return (
                "<h1>rrecon API is running</h1><p>The overlapscope UI is not built. Run "
                "<code>cd web &amp;&amp; npm install &amp;&amp; npm run build</code>, or "
                "<code>npm run dev</code> and open http://localhost:5173.</p>"
            )

    return app


def free_port(host: str, port: int, tries: int = 50) -> int:
    for p in range(port, port + tries):
        with socket.socket() as s:
            try:
                s.bind((host, p))
                return p
            except OSError:
                continue
    raise OSError(f"no free port between {port} and {port + tries - 1}")


def serve(folder: Path | None, host: str = "127.0.0.1", port: int = 8766, open_browser=True):
    import uvicorn

    port = free_port(host, port)
    url = f"http://{host}:{port}/"
    print(f"overlapscope: {url}" + (f"  ({folder})" if folder else ""))
    if open_browser:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    uvicorn.run(create_app(folder), host=host, port=port, log_level="warning")
