#!/usr/bin/env python3
"""Rebuild the self-contained WordPress HTML from its template.

Inlines every {{BASE64:<name>}} placeholder in
v2_followup_wordpress.template.html with the base64 data URI of
v2_figures/<name>.png and writes v2_followup_wordpress.html.

The template already carries the design-rule SVG inline and the article
body; only the molecule figures need to be embedded at build time so the
template stays readable and diffable.

Usage:
    python blog/build_wordpress_html.py
"""
import base64
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TEMPLATE = ROOT / "v2_followup_wordpress.template.html"
OUTPUT = ROOT / "v2_followup_wordpress.html"
FIGDIR = ROOT / "v2_figures"


def build() -> str:
    html = TEMPLATE.read_text(encoding="utf-8")

    def inline(match: re.Match) -> str:
        name = match.group(1)
        path = FIGDIR / f"{name}.png"
        if not path.exists():
            raise FileNotFoundError(f"figure for {name}: {path}")
        data = base64.b64encode(path.read_bytes()).decode("ascii")
        return f"data:image/png;base64,{data}"

    html, n = re.subn(r"\{\{BASE64:([a-z0-9_]+)\}\}", inline, html)
    if n == 0:
        raise RuntimeError("no {{BASE64:...}} placeholders found in template")
    return html


def main() -> None:
    html = build()
    OUTPUT.write_text(html, encoding="utf-8")
    print(f"wrote {OUTPUT} ({len(html):,} bytes)")


if __name__ == "__main__":
    main()
