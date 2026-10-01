"""How to cite PolyzyMD.

The repository's ``CITATION.cff`` is the source of these values, and a test
checks that they agree with it. ``polyzymd study check`` prints
:func:`citation_line`, so people and agents who use a study folder know
which framework produced it and how to credit it.
"""

from __future__ import annotations

TITLE = "PolyzyMD: build, run and analyse molecular dynamics of complex biomolecular systems"
AUTHORS = "Laforet, Joseph R., Jr."
REPOSITORY = "https://github.com/joelaforet/polyzymd"


def citation_line() -> str:
    """Return a one-line citation of the installed PolyzyMD version."""
    from polyzymd import __version__

    return f"{AUTHORS} {TITLE} (version {__version__}). {REPOSITORY}"
