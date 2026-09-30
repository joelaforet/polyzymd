"""Figure checks shared by the analyses tests: does a figure draw an uncertainty, and say what it is."""

from __future__ import annotations

from typing import Any


def figure_draws_uncertainty(fig: Any) -> bool:
    """Return whether any axes on *fig* draws an error bar or a shaded band."""

    from matplotlib.collections import PolyCollection

    return any(
        any(
            getattr(container, "has_yerr", False) or getattr(container, "has_xerr", False)
            for container in ax.containers
        )
        or any(isinstance(collection, PolyCollection) for collection in ax.collections)
        for ax in fig.axes
    )


def figure_has_uncertainty_footnote(fig: Any) -> bool:
    """Return whether *fig* carries a text naming what its uncertainty is."""

    return any(
        ("95%" in text.get_text() or "SEM" in text.get_text()) and "replicates" in text.get_text()
        for text in fig.texts
    )
