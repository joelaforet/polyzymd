"""How to cite PolyzyMD.

The repository's ``CITATION.cff`` is the source of these values, and a test
checks that they agree with it. ``polyzymd study check`` prints
:func:`citation_line`, so people and agents who use a study folder know
which framework produced it and how to credit it, and ``polyzymd study
freeze`` writes :func:`software_reference` and :func:`paper_reference` into
each study's ``CITATION.cff`` and ``.zenodo.json``.

When the PolyzyMD paper is published, set :data:`PAPER` (and ``DOI`` once
PolyzyMD has a Zenodo DOI) here and in ``CITATION.cff``; every study frozen
afterwards cites it.
"""

from __future__ import annotations

from typing import Any

TITLE = "PolyzyMD: build, run and analyse molecular dynamics of complex biomolecular systems"
AUTHORS = "Laforet, Joseph R., Jr."
REPOSITORY = "https://github.com/joelaforet/polyzymd"
#: The authors, as CITATION.cff persons.
PEOPLE: list[dict[str, str]] = [
    {
        "family-names": "Laforet",
        "given-names": "Joseph R.",
        "name-suffix": "Jr.",
        "affiliation": "Shirts group, Department of Chemical and Biological Engineering, "
        "University of Colorado Boulder",
    }
]
#: PolyzyMD's own DOI, once it has one.
DOI: str | None = None
#: The PolyzyMD paper, once it exists: a CITATION.cff reference with at least
#: ``title`` and ``status`` (or ``doi``). ``None`` leaves it out.
PAPER: dict[str, Any] | None = None


def citation_line() -> str:
    """Return a one-line citation of the installed PolyzyMD version."""
    from polyzymd import __version__

    return f"{AUTHORS} {TITLE} (version {__version__}). {REPOSITORY}"


def software_reference() -> dict[str, Any]:
    """Return the installed PolyzyMD version as a CITATION.cff reference."""
    from polyzymd import __version__

    reference: dict[str, Any] = {
        "type": "software",
        "title": TITLE,
        "authors": [dict(person) for person in PEOPLE],
        "version": __version__,
        "repository-code": REPOSITORY,
    }
    if DOI:
        reference["doi"] = DOI
    return reference


def paper_reference() -> dict[str, Any] | None:
    """Return the PolyzyMD paper as a CITATION.cff reference, or ``None`` before it exists."""
    if PAPER is None:
        return None
    return {"type": "article", "authors": [dict(p) for p in PEOPLE], **PAPER}
