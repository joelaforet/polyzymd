"""Known-answer tests for how figures draw and describe uncertainty.

Grossfield et al. (2018, LiveCoMS 1:5067) ask that authors graph 95 percent
confidence intervals and describe the meaning and basis of the uncertainties in
every figure. These tests pin the interval arithmetic behind an error bar or a
shaded band, and check that a figure which draws either one also carries the
footnote saying what it is.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from polyzymd.analyses._framework.comparison_models import BasePlotSettings  # noqa: E402
from polyzymd.analyses.exceptions import StatisticsError  # noqa: E402
from polyzymd.analyses.shared.plotting import (  # noqa: E402
    add_uncertainty_footnote,
    band_half_widths,
    error_bar_half_widths,
    plugin_plot_settings,
    resolve_error_bar,
)
from tests.analyses.conftest import (  # noqa: E402
    figure_draws_uncertainty,
    figure_has_uncertainty_footnote,
)

T_FACTOR_N3 = 4.302652729749462

PLUGINS_WITH_PLOT_SETTINGS = (
    "rmsd",
    "rmsf",
    "rg",
    "sasa",
    "contacts",
    "distances",
    "secondary_structure",
    "catalytic_triad",
)


def _footnote_texts(fig: "plt.Figure") -> list[str]:
    """Return the figure-level texts that read as an uncertainty footnote.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to inspect.

    Returns
    -------
    list of str
        Texts naming either a 95 percent interval or a standard error across
        replicates.
    """
    texts = [text.get_text() for text in fig.texts]
    return [text for text in texts if ("95%" in text or "SEM" in text) and "replicates" in text]


def _draws_uncertainty(fig: "plt.Figure") -> bool:
    """Return whether any axes on *fig* draws an error bar or a shaded band.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to inspect.

    Returns
    -------
    bool
        True when an errorbar container or a fill_between collection is present.
    """
    from matplotlib.collections import PolyCollection

    for ax in fig.axes:
        if len(ax.containers) and any(
            getattr(container, "has_yerr", False) or getattr(container, "has_xerr", False)
            for container in ax.containers
        ):
            return True
        if any(isinstance(collection, PolyCollection) for collection in ax.collections):
            return True
    return False


class TestErrorBarWidths:
    """A drawn error bar must span the interval the figure claims."""

    def test_bar_half_widths_scale_the_sem_by_the_coverage_factor(self) -> None:
        """A three-replicate bar spans 4.303 standard errors either side."""

        half_widths = error_bar_half_widths(
            [0.1, 0.2],
            [[1.0, 1.1, 1.2], [2.0, 2.2, 2.4]],
            error_bar="ci95",
        )

        assert half_widths == pytest.approx([0.1 * T_FACTOR_N3, 0.2 * T_FACTOR_N3])

    def test_sem_mode_draws_one_standard_error(self) -> None:
        """Asking for SEM bars draws exactly one standard error."""

        assert error_bar_half_widths([0.1], [[1.0, 1.1, 1.2]], error_bar="sem") == pytest.approx(
            [0.1]
        )

    def test_each_bar_uses_its_own_replicate_count(self) -> None:
        """Coverage factors differ per bar when replicate counts differ."""

        half_widths = error_bar_half_widths(
            [0.1, 0.1],
            [[1.0, 1.1, 1.2], [2.0, 2.1, 2.2, 2.3, 2.4]],
            error_bar="ci95",
        )

        assert half_widths[0] > half_widths[1]

    def test_singleton_bars_get_no_error_bar(self) -> None:
        """A bar backed by one replicate shows no interval at all."""

        assert error_bar_half_widths([0.1], [[1.0]], error_bar="ci95") is None

    def test_band_half_widths_use_the_coverage_factor(self) -> None:
        """A shaded band over replicate traces is also a 95 percent interval."""

        matrix = np.array([[2.0, 2.0], [2.2, 2.2], [2.4, 2.4]])
        band = band_half_widths(matrix, error_bar="ci95")

        expected = T_FACTOR_N3 * 0.2 / math.sqrt(3.0)
        assert band == pytest.approx([expected, expected])

    def test_band_is_none_for_one_replicate(self) -> None:
        """One trace gives no band."""

        assert band_half_widths(np.array([[2.0, 2.0]]), error_bar="ci95") is None

    @pytest.mark.parametrize("bad", ["ci", "stdev", ""])
    def test_unknown_error_bar_is_a_typed_error(self, bad: str) -> None:
        """An unsupported interval name fails loudly with a typed error."""

        with pytest.raises(StatisticsError):
            error_bar_half_widths([0.1], [[1.0, 1.1]], error_bar=bad)
        with pytest.raises(StatisticsError):
            band_half_widths(np.array([[1.0], [1.1]]), error_bar=bad)


class TestErrorBarSetting:
    """The plugin's own setting must reach the plotter."""

    def test_default_is_the_confidence_interval(self) -> None:
        """A plugin that declares no preference gets the 95 percent interval."""

        assert BasePlotSettings().error_bar == "ci95"

    def test_plugin_settings_win_over_the_global_object(self) -> None:
        """A user's per-analysis choice is honoured, not silently dropped."""

        plugin = BasePlotSettings(error_bar="sem")

        assert resolve_error_bar(plugin, object()) == "sem"

    def test_missing_settings_fall_back_to_the_default(self) -> None:
        """A plugin with no plot settings model still gets a valid choice."""

        assert resolve_error_bar(None, object()) == "ci95"

    @pytest.mark.parametrize("analysis_name", PLUGINS_WITH_PLOT_SETTINGS)
    def test_every_plugin_setting_is_reachable_from_the_global_object(
        self, analysis_name: str
    ) -> None:
        """The global plot settings expose each plugin's error_bar choice."""

        from polyzymd.config.comparison import PlotSettings

        settings = plugin_plot_settings(PlotSettings(), analysis_name)

        assert settings is not None, f"{analysis_name} settings not reachable"
        assert getattr(settings, "error_bar", None) == "ci95"


class TestFootnote:
    """Every figure that draws an uncertainty says what it is."""

    def test_footnote_states_the_interval_and_the_sampling_unit(self) -> None:
        """The footnote names the coverage, the replicates and the window."""

        fig, _ax = plt.subplots()
        try:
            add_uncertainty_footnote(fig, error_bar="ci95", n_replicates=3, equilibration="10ns")
            footnotes = _footnote_texts(fig)
        finally:
            plt.close(fig)

        assert footnotes
        assert "95%" in footnotes[0]
        assert "n = 3" in footnotes[0]
        assert "10ns" in footnotes[0]

    def test_sem_footnote_warns_that_it_is_not_a_95_percent_interval(self) -> None:
        """Choosing SEM bars must not let a reader assume 95 percent coverage."""

        fig, _ax = plt.subplots()
        try:
            add_uncertainty_footnote(fig, error_bar="sem", n_replicates=3)
            footnote = _footnote_texts(fig)[0]
        finally:
            plt.close(fig)

        assert "SEM" in footnote
        assert "not a 95% interval" in footnote

    def test_detector_sees_an_error_bar(self) -> None:
        """The uncertainty detector recognises a bar chart with yerr."""

        fig, ax = plt.subplots()
        try:
            ax.bar([0, 1], [1.0, 2.0], yerr=[0.1, 0.2])
            assert figure_draws_uncertainty(fig)
        finally:
            plt.close(fig)

    def test_detector_sees_a_shaded_band(self) -> None:
        """The uncertainty detector recognises a fill_between band."""

        fig, ax = plt.subplots()
        try:
            ax.fill_between([0, 1], [0.9, 1.9], [1.1, 2.1])
            assert figure_draws_uncertainty(fig)
        finally:
            plt.close(fig)

    def test_detector_ignores_a_plain_line(self) -> None:
        """A figure with no uncertainty is not asked for a footnote."""

        fig, ax = plt.subplots()
        try:
            ax.plot([0, 1], [1.0, 2.0])
            assert not figure_draws_uncertainty(fig)
        finally:
            plt.close(fig)


class TestEveryPluginFigureIsAudited:
    """Each plugin must render at least one audited uncertainty figure.

    The conftest audit only sees figures that reach the real ``save_figure``.
    A test that installs its own stub bypasses it, so this records which plugin
    modules the audit actually inspected during the session and fails if any
    plugin that draws uncertainty is never covered.
    """

    def test_conftest_audit_covers_the_plotter_modules(self) -> None:
        """The audit must be installed on every plotter module."""

        import importlib

        from tests.analyses.conftest import _PLOTTER_MODULES

        for module_name in _PLOTTER_MODULES:
            module = importlib.import_module(module_name)
            assert hasattr(module, "save_figure"), module_name
            assert module.save_figure.__name__ == "_checked", (
                f"{module_name} save_figure is not wrapped by the footnote audit; "
                "a test in this session replaced it"
            )

    @pytest.mark.parametrize("analysis_name", PLUGINS_WITH_PLOT_SETTINGS + ("hydrogen_bonds",))
    def test_plotter_module_annotates_every_uncertainty_figure(
        self, analysis_name: str, tmp_path: "Path"
    ) -> None:
        """Rendering a bar chart through a plugin's save_figure must be audited.

        Builds a minimal figure that draws an error bar, routes it through that
        plugin's own ``save_figure`` reference, and checks that the audit
        rejects it without a footnote and accepts it with one.
        """

        import importlib

        module = importlib.import_module(f"polyzymd.analyses.{analysis_name}._plotters")

        from polyzymd.config.comparison import PlotSettings

        settings = PlotSettings(output_dir=tmp_path)

        fig, ax = plt.subplots()
        try:
            ax.bar([0], [1.0], yerr=[0.1])
            with pytest.raises(AssertionError, match="without a footnote"):
                module.save_figure(
                    fig, tmp_path / f"{analysis_name}_bad.png", settings, close=False
                )

            add_uncertainty_footnote(fig, error_bar="ci95", n_replicates=3)
            module.save_figure(fig, tmp_path / f"{analysis_name}_good.png", settings, close=False)
        finally:
            plt.close(fig)

        assert (tmp_path / f"{analysis_name}_good.png").exists()
        assert not (tmp_path / f"{analysis_name}_bad.png").exists()


class TestRealPlottersCarryTheFootnote:
    """Render the bar plotters whose own tests stub save_figure.

    Those tests replace ``save_figure`` before the conftest audit can see the
    figure, so deleting a footnote call in RMSD or Rg would otherwise go
    unnoticed. These render the same plotters with the real ``save_figure`` and
    read the footnote back off the saved figure.
    """

    @staticmethod
    def _captured_figures(monkeypatch: "pytest.MonkeyPatch") -> list:
        """Record every figure handed to save_figure without closing it."""

        from polyzymd.analyses.shared import plotting

        captured: list = []
        original = plotting.save_figure

        def _capture(fig, output_path, plot_settings, **kwargs):
            captured.append(fig)
            kwargs["close"] = False
            return original(fig, output_path, plot_settings, **kwargs)

        for module_name in ("polyzymd.analyses.rmsd._plotters", "polyzymd.analyses.rg._plotters"):
            import importlib

            monkeypatch.setattr(importlib.import_module(module_name), "save_figure", _capture)
        return captured

    def test_rg_comparison_bars_are_footnoted(
        self, tmp_path: Path, monkeypatch: "pytest.MonkeyPatch"
    ) -> None:
        """The Rg comparison bar chart carries the footnote when rendered."""

        from datetime import datetime

        from polyzymd.analyses.base import PlotContext
        from polyzymd.analyses.rg._comparison_results import (
            RgComparisonResult,
            RgConditionSummary,
            RgRunSummary,
        )
        from polyzymd.analyses.rg._plotters import plot_rg_comparison_bars
        from polyzymd.config.comparison import PlotSettings

        captured = self._captured_figures(monkeypatch)
        comparison = RgComparisonResult(
            metric="mean_rg",
            name="rg_compare",
            n_runs=1,
            run_labels=["protein_rg"],
            control_label="Control",
            conditions=[
                RgConditionSummary(
                    label="Control",
                    config_path="/fake/control.yaml",
                    n_replicates=3,
                    run_summaries=[
                        RgRunSummary(
                            label="protein_rg",
                            selection="protein",
                            mean_rg=10.5,
                            sem_rg=0.5,
                            per_replicate_means=[10.0, 10.5, 11.0],
                            replicates=[1, 2, 3],
                            n_replicates=3,
                        )
                    ],
                )
            ],
            pairwise_comparisons=[],
            anova_by_run=None,
            ranking_by_run={"protein_rg": ["Control"]},
            equilibration_time="10ns",
            created_at=datetime.now(),
            polyzymd_version="test",
        )
        ctx = PlotContext(
            conditions=[],
            analysis_dirs={},
            results_dir=tmp_path,
            output_dir=tmp_path / "figures",
            settings=None,
            plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
            equilibration="10ns",
        )

        plot_rg_comparison_bars(ctx, comparison)

        assert captured, "no figure was rendered"
        try:
            for fig in captured:
                assert figure_draws_uncertainty(fig)
                assert figure_has_uncertainty_footnote(fig), "rendered Rg bars lack a footnote"
                assert "n = 3" in _footnote_texts(fig)[0]
                assert "10ns" in _footnote_texts(fig)[0]
        finally:
            for fig in captured:
                plt.close(fig)
