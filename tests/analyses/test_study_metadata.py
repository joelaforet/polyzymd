"""Checks of the publishing metadata block."""

from __future__ import annotations

from polyzymd.analyses.study_metadata import check_metadata


def test_todo_placeholders_are_a_gap() -> None:
    """Authors and a paper title left as the template's TODO each give a warning."""
    raw = {
        "title": "Lipase in water",
        "authors": [{"family-names": "TODO", "given-names": "TODO", "orcid": "TODO"}],
        "related": {"paper": {"title": "TODO", "doi": "TODO"}},
    }
    _, warnings = check_metadata(raw)
    (todo,) = [w for w in warnings if w.startswith("TODO placeholders")]
    assert "metadata.authors[0].family-names" in todo
    assert "metadata.related.paper.title" in todo
    assert "doi" not in todo
    assert not any("placeholders" in w for w in check_metadata({"title": "Lipase"})[1])
