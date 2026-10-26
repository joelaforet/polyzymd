#!/usr/bin/env python
"""Check that the docs search puts the right page first for what users type.

Each query in ``docs_search_queries.yaml`` (next to this script) runs through
the real Sphinx search of a built HTML site: a headless Chromium opens
``search.html?q=<query>``, the page's own JavaScript searches
``searchindex.js``, and the script reads the result list. It reports the rank
of the best expected page (1 = top) and the top-1 and top-3 hit rates.

Run it before a release, after a docs build::

    sphinx-build -E -b html docs/source /tmp/docs_html
    python scripts/check_docs_search.py /tmp/docs_html

or let it build the site into a temporary folder::

    python scripts/check_docs_search.py --build

It needs no network. It needs a Chromium or Chrome binary: set ``CHROME`` to
its path, or have ``chromium``, ``google-chrome`` or ``chrome-headless-shell``
on ``PATH`` (``playwright install chromium-headless-shell`` provides one under
``~/.cache/ms-playwright``). ``--fail-under-top3 0.9`` makes it exit non-zero
when fewer than 90% of the queries have an expected page in the top three.
"""

from __future__ import annotations

import argparse
import html
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.parse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
QUERIES = Path(__file__).with_name("docs_search_queries.yaml")

_RESULTS_RE = re.compile(r'<ul[^>]*class="search"[^>]*>(.*?)</ul>', re.S)
_LINK_RE = re.compile(r'<li[^>]*>\s*<a href="([^"#?]+)[^"]*"[^>]*>(.*?)</a>', re.S)


def find_chrome() -> str:
    """Return a headless-capable Chromium binary, or exit with a message."""
    if os.environ.get("CHROME"):
        return os.environ["CHROME"]
    for name in ("chrome-headless-shell", "chromium", "chromium-browser", "google-chrome"):
        path = shutil.which(name)
        if path:
            return path
    cache = Path.home() / ".cache" / "ms-playwright"
    found = sorted(cache.glob("chromium*/chrome*/chrome-headless-shell")) + sorted(
        cache.glob("chromium-*/chrome-linux*/chrome")
    )
    if found:
        return str(found[-1])
    sys.exit("No Chromium found: set CHROME to a chrome or chrome-headless-shell binary.")


def search(chrome: str, html_dir: Path, query: str) -> list[str]:
    """Return the docnames of the search results for ``query``, best first."""
    url = (html_dir / "search.html").resolve().as_uri() + "?" + urllib.parse.urlencode({"q": query})
    dom = subprocess.run(
        [chrome, "--headless", "--no-sandbox", "--disable-gpu",
         "--virtual-time-budget=8000", "--dump-dom", url],
        capture_output=True, text=True, timeout=120, check=False,
    ).stdout
    match = _RESULTS_RE.search(dom)
    if not match:
        return []
    docnames: list[str] = []
    for href, _title in _LINK_RE.findall(match.group(1)):
        docname = html.unescape(href).removesuffix(".html")
        if docname not in docnames:
            docnames.append(docname)
    return docnames


def build(out: Path) -> None:
    subprocess.run(
        [sys.executable, "-m", "sphinx", "-E", "-b", "html", "-q",
         str(ROOT / "docs" / "source"), str(out)],
        check=True,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("html_dir", nargs="?", type=Path, help="a built HTML site")
    parser.add_argument("--build", action="store_true", help="build the site first")
    parser.add_argument("--queries", type=Path, default=QUERIES)
    parser.add_argument("--show", type=int, default=3, help="results to print per query")
    parser.add_argument("--fail-under-top3", type=float, default=None)
    args = parser.parse_args(argv)

    with tempfile.TemporaryDirectory() as tmp:
        html_dir = args.html_dir
        if args.build:
            html_dir = Path(tmp) / "html"
            build(html_dir)
        if html_dir is None or not (html_dir / "searchindex.js").exists():
            parser.error("give a built HTML folder (with searchindex.js) or --build")

        cases = yaml.safe_load(args.queries.read_text(encoding="utf-8"))
        chrome = find_chrome()
        with ThreadPoolExecutor(max_workers=min(8, os.cpu_count() or 1)) as pool:
            results = list(pool.map(lambda c: search(chrome, html_dir, c["query"]), cases))

    top1 = top3 = 0
    width = max(len(c["query"]) for c in cases)
    for case, found in zip(cases, results):
        ranks = [found.index(page) + 1 for page in case["expect"] if page in found]
        rank = min(ranks) if ranks else None
        top1 += rank == 1
        top3 += rank is not None and rank <= 3
        shown = ", ".join(found[: args.show]) or "(no results)"
        print(f"{case['query']:<{width}}  {rank if rank else '-':>3}  {shown}")
    n = len(cases)
    print(f"\ntop-1: {top1}/{n} ({top1 / n:.0%})   top-3: {top3}/{n} ({top3 / n:.0%})")
    if args.fail_under_top3 is not None and top3 / n < args.fail_under_top3:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
