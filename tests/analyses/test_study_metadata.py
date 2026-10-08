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


def test_placeholder_dois_of_the_documentation_are_gaps() -> None:
    """DOIs such as 10.5281/zenodo.NNNNNNN are warned about and never written as DOIs."""
    from polyzymd.analyses.study_metadata import is_placeholder

    for doi in ("10.5281/zenodo.NNNNNNN", "10.5281/zenodo.0000000", "10.1021/...", "10.xxxx/y"):
        assert is_placeholder(doi), doi
    for doi in ("10.5281/zenodo.1234567", "10.1021/acs.jctc.0c00100", "10.1038/nnano.2010.1"):
        assert not is_placeholder(doi), doi
    meta, warnings = check_metadata(
        {
            "doi": "10.5281/zenodo.NNNNNNN",
            "funding": [{"funder": "NSF", "award": "1", "funder_doi": "10.13039/..."}],
        }
    )
    assert meta["doi"] is None
    assert any(w.startswith("metadata.doi is not set") for w in warnings)
    assert "funder_doi" not in meta["funding"][0]
    assert any("funder_doi" in w for w in warnings)
