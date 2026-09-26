"""Shared fixtures for the analyses test suite."""

from __future__ import annotations

import importlib
from typing import Any

import pytest

_PLOTTER_MODULES = tuple(
    f"polyzymd.analyses.{name}._plotters"
    for name in (
        "rmsd",
        "rmsf",
        "sasa",
        "contacts",
        "distances",
        "hydrogen_bonds",
        "secondary_structure",
        "catalytic_triad",
    )
)


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
    """Fail any test whose plotter saves an uncertainty figure with no footnote.

    Grossfield et al. (2018) require every figure to describe the meaning and
    basis of its uncertainties. Wrapping ``save_figure`` in every plotter module
    turns the existing body of plot tests into that check on real rendered
    figures. A test that installs its own ``save_figure`` stub bypasses the
    audit; ``tests/analyses/scientific/test_uncertainty_plots.py`` renders those
    paths directly instead.
    """

    pytest.importorskip("matplotlib")

    for module_name in _PLOTTER_MODULES:
        module = importlib.import_module(module_name)
        original = getattr(module, "save_figure", None)
        if original is None:
            continue

        def _checked(
            fig: Any,
            *args: Any,
            _original: Any = original,
            _name: str = module_name,
            **kwargs: Any,
        ) -> Any:
            if figure_draws_uncertainty(fig) and not figure_has_uncertainty_footnote(fig):
                raise AssertionError(
                    f"{_name} saved a figure that draws an uncertainty "
                    "without a footnote saying what it is"
                )
            return _original(fig, *args, **kwargs)

        monkeypatch.setattr(module, "save_figure", _checked)
