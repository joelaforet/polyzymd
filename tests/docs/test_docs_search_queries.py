"""The docs search query list names real pages, and index entries are searchable.

The ranking itself needs a built site: run ``scripts/check_docs_search.py``.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs" / "source"
CASES = yaml.safe_load((ROOT / "scripts" / "docs_search_queries.yaml").read_text(encoding="utf-8"))


@pytest.mark.parametrize("case", CASES, ids=[case["query"] for case in CASES])
def test_expected_pages_exist(case: dict) -> None:
    assert case["expect"], case["query"]
    for docname in case["expect"]:
        assert (DOCS / f"{docname}.md").exists() or (DOCS / f"{docname}.rst").exists(), docname


def test_queries_are_unique() -> None:
    queries = [case["query"] for case in CASES]
    assert len(queries) == len(set(queries))


def test_index_entries_are_lowercase() -> None:
    # Sphinx search compares the lower-cased query with the entry as written,
    # so an entry with a capital letter never matches.
    for path in DOCS.rglob("*.md"):
        for line in re.findall(r"^```\{index\} (.*)$", path.read_text(encoding="utf-8"), re.M):
            assert line == line.lower(), f"{path.relative_to(DOCS)}: {line}"
