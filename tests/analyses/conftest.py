"""Shared fixtures for the analyses tests, and the footnote audit of every saved figure."""

from __future__ import annotations

from typing import Any

import pytest


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


@pytest.fixture(autouse=True)
def audit_plot_uncertainty_footnotes(monkeypatch: pytest.MonkeyPatch) -> None:
    """Fail any test that saves a figure drawing an uncertainty with no footnote.

    Grossfield et al. (2018, LiveCoMS 1:5067) require every figure to describe
    the meaning and basis of its uncertainties. Every study figure in
    :mod:`polyzymd.analyses.figures` is saved through
    :func:`polyzymd.analyses.shared.plotting.save_figure`, which this fixture
    wraps, so each test that draws a figure also checks its footnote. A test
    that replaces ``save_figure`` with a stub that does not call the one it
    found bypasses the audit.
    """

    pytest.importorskip("matplotlib")
    from polyzymd.analyses.shared import plotting

    original = plotting.save_figure

    def checked(fig: Any, *args: Any, **kwargs: Any) -> Any:
        if figure_draws_uncertainty(fig) and not figure_has_uncertainty_footnote(fig):
            raise AssertionError(
                "A figure that draws an uncertainty was saved without a footnote saying "
                "what it is and what it is computed across"
            )
        return original(fig, *args, **kwargs)

    monkeypatch.setattr(plotting, "save_figure", checked)
