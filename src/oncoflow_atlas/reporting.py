"""One offline report renderer for CLI outputs and GitHub Pages."""

from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path

from .analysis import clean_json

WEB = Path(__file__).parent / "web"


def write_json(path: str | Path, value: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(clean_json(value), indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def render_report(reports: list[dict], output: str | Path) -> None:
    template = (WEB / "index.html").read_text(encoding="utf-8")
    css = (WEB / "styles.css").read_text(encoding="utf-8")
    scripts = [(WEB / name).read_text(encoding="utf-8") for name in ("analysis.js", "app.js")]
    payload = json.dumps(reports, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
    payload = (
        payload.replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )
    hashes = ["'sha256-" + base64.b64encode(hashlib.sha256(s.encode()).digest()).decode() + "'" for s in scripts]
    policy = (
        "default-src 'none'; script-src "
        + " ".join(hashes)
        + "; style-src 'unsafe-inline'; img-src data: blob:; connect-src 'none'; base-uri 'none'; form-action 'none'; object-src 'none'"
    )
    html = (
        template.replace("__CSP__", policy)
        .replace("__STYLES__", css)
        .replace("__PAYLOAD__", payload)
        .replace("__ANALYSIS__", scripts[0])
        .replace("__APP__", scripts[1])
    )
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(html, encoding="utf-8")
    temporary.replace(path)
