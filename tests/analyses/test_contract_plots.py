"""Figures the contract runner generates from an observable's kind."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, ClassVar, Sequence

import pytest
from pydantic import BaseModel

from polyzymd.analyses import contract_plots
from polyzymd.analyses.contract import Observable, ObservableAggregate, iter_frames
from polyzymd.analyses.contract_plots import ContractPlotSettings
from polyzymd.analyses.contract_runner import contract_analysis
from tests.analyses.conftest import make_simulation_config, make_synthetic_universe

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

STUDENT_T_AT_THREE = 4.302652729911275


class SyntheticSettings(BaseModel):
    """Which kind of observable the synthetic plugin should report."""

    kind: str = "mean_of_timeseries"
    n_index: int = 0
    integer_index: bool = True


class Synthetic:
    """Reports one observable of a chosen kind, scaled by the replicate."""

    name: ClassVar[str] = "synthetic"
    Settings: ClassVar[type[BaseModel]] = SyntheticSettings
    references: ClassVar[tuple[str, ...]] = ()

    def compute(
        self, universe: Any, frames: Any, settings: SyntheticSettings
    ) -> Sequence[Observable]:
        """Return one observable whose values depend on the replicate's scale."""
        import numpy as np

        scale = float(universe.select_atoms("all").radius_of_gyration())
        n_frames = sum(1 for _ in iter_frames(universe, frames))
        if settings.kind == "profile":
            index = (
                np.arange(settings.n_index, dtype=float)
                if settings.integer_index
                else np.linspace(0.25, 0.25 + settings.n_index, settings.n_index)
            )
            return [
                Observable(
                    name="occupancy",
                    kind="profile",
                    unit="A",
                    values=(scale + 0.1 * index).tolist(),
                    index=index.tolist(),
                )
            ]
        if settings.kind == "fraction":
            filled = min(n_frames, int(round(scale * 4)))
            return [
                Observable(
                    name="bound",
                    kind="fraction",
                    unit="fraction",
                    values=[1.0] * filled + [0.0] * (n_frames - filled),
                )
            ]
        if settings.kind == "fluctuation":
            values = (scale * np.sin(np.arange(n_frames))).tolist()
        else:
            values = (scale + 0.05 * np.arange(n_frames)).tolist()
        return [Observable(name="size", kind=settings.kind, unit="A", values=values)]


SyntheticAnalysis = contract_analysis(Synthetic)


def _universes(base: float) -> Callable[[int], Any]:
    """Replicate factory whose radius of gyration grows by 0.5 per replicate."""
    return lambda replicate: make_synthetic_universe(scale=base + 0.5 * replicate, n_frames=24)


@pytest.fixture
def rendered(monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    """Keep every saved figure by file name, leaving the footnote audit in place."""
    figures: dict[str, Any] = {}
    audited = contract_plots.save_figure

    def _record(fig: Any, output_path: Path, plot_settings: Any, **kwargs: Any) -> Path:
        figures[Path(output_path).name] = fig
        return audited(fig, output_path, plot_settings, **kwargs)

    monkeypatch.setattr(contract_plots, "save_figure", _record)
    return figures


@pytest.fixture
def render(tmp_path: Path, run_contract_analysis: Any) -> Callable[..., list[Path]]:
    """Run the synthetic plugin over two conditions and plot the result."""

    def _render(settings: SyntheticSettings, labels: Sequence[str] = ("A", "B")) -> list[Path]:
        from polyzymd.analyses._framework.contexts import Condition, PlotContext
        from polyzymd.config.comparison import PlotSettings

        for index, label in enumerate(labels):
            run_contract_analysis(
                SyntheticAnalysis,
                settings,
                _universes(1.0 + 2.0 * index),
                label=label,
                root=tmp_path,
            )
        conditions = [
            Condition(
                label=label,
                config_path=tmp_path / f"{label}.yaml",
                replicates=(1, 2, 3),
                sim_config=make_simulation_config(label),
            )
            for label in labels
        ]
        return SyntheticAnalysis().plot(
            PlotContext(
                conditions=conditions,
                analysis_dirs={
                    label: tmp_path / "analysis" / label / "synthetic" for label in labels
                },
                results_dir=tmp_path / "results",
                output_dir=tmp_path / "figures",
                settings=settings,
                plot_settings=PlotSettings(),
                control_label=labels[0],
                equilibration="5ns",
            )
        )

    return _render


def _aggregate(tmp_path: Path, label: str, name: str) -> ObservableAggregate:
    """Read one observable's aggregate back off disk."""
    from polyzymd.analyses.mda.store import ArtifactStore

    directory = tmp_path / "analysis" / label / "synthetic" / "aggregated"
    artifact = ArtifactStore(directory).read_condition_result("result.json")
    payloads = {
        payload["name"]: payload for payload in artifact.payload["observables"]  # type: ignore[index]
    }
    return ObservableAggregate.model_validate(payloads[name])


def _error_half_widths(ax: Any) -> list[float]:
    """Half height of every error bar drawn on an axes, read from the artists."""
    half_widths = []
    for container in ax.containers:
        errorbar = getattr(container, "errorbar", None)
        if errorbar is None:
            continue
        for collection in errorbar[2]:
            for segment in collection.get_segments():
                half_widths.append((float(segment[1][1]) - float(segment[0][1])) / 2.0)
    return half_widths


def _footnotes(fig: Any) -> str:
    """Every figure-level text joined, so a footnote can be matched in one go."""
    return " | ".join(text.get_text() for text in fig.texts)


@pytest.mark.parametrize("kind", ["mean_of_timeseries", "fluctuation"])
def test_time_series_kinds_get_bars_and_a_series_panel(
    kind: str, tmp_path: Path, render: Callable[..., list[Path]]
) -> None:
    """Both time-series kinds produce a comparison chart and a per-frame panel."""
    paths = render(SyntheticSettings(kind=kind))

    assert [path.name for path in paths] == [
        "synthetic_size_comparison.png",
        "synthetic_size_timeseries.png",
    ]
    assert all(path.exists() for path in paths)


def test_bars_span_the_student_t_interval_across_replicates(
    tmp_path: Path, render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """Each error bar is the 95 percent Student t interval, 4.303 SEM at n = 3."""
    render(SyntheticSettings(kind="mean_of_timeseries"))

    figure = rendered["synthetic_size_comparison.png"]
    half_widths = _error_half_widths(figure.axes[0])
    expected = [
        _aggregate(tmp_path, label, "size").sem * STUDENT_T_AT_THREE for label in ("A", "B")
    ]

    assert half_widths == pytest.approx(expected, rel=1e-6)
    assert all(half_width > 0.0 for half_width in half_widths)


def test_the_comparison_figure_says_what_its_error_bars_are(
    render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """The footnote names the interval, the replicate count and the window."""
    render(SyntheticSettings(kind="mean_of_timeseries"))

    footnote = _footnotes(rendered["synthetic_size_comparison.png"])

    assert "95% CI (Student t)" in footnote
    assert "n = 3 replicates" in footnote
    assert "t >= 5ns" in footnote


def test_a_fraction_is_drawn_on_a_clamped_axis(
    render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """A fraction gets bars only, on an axis pinned to [0, 1]."""
    paths = render(SyntheticSettings(kind="fraction"))

    assert [path.name for path in paths] == ["synthetic_bound_comparison.png"]
    figure = rendered["synthetic_bound_comparison.png"]
    assert figure.axes[0].get_ylim() == (0.0, 1.0)
    assert figure.axes[0].get_ylabel() == "bound (fraction)"


def test_a_categorical_profile_is_drawn_as_bars_up_to_thirty_categories(
    render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """Thirty residue labels still fit as grouped bars, one group per label."""
    from matplotlib.container import BarContainer

    render(SyntheticSettings(kind="profile", n_index=30))

    axes = rendered["synthetic_occupancy_comparison.png"].axes[0]
    bars = [container for container in axes.containers if isinstance(container, BarContainer)]

    assert [len(container) for container in bars] == [30, 30]
    assert len(axes.get_xticks()) == 30
    assert _error_half_widths(axes)


def test_a_longer_profile_falls_back_to_lines_with_a_band(
    render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """Thirty one categories is past the limit, so the profile becomes lines."""
    from matplotlib.collections import PolyCollection

    render(SyntheticSettings(kind="profile", n_index=31))

    figure = rendered["synthetic_occupancy_comparison.png"]
    axes = figure.axes[0]

    assert not axes.containers
    assert len(axes.lines) == 8
    assert any(isinstance(collection, PolyCollection) for collection in axes.collections)
    assert "95% CI (Student t)" in _footnotes(figure)


def test_a_continuous_profile_index_is_never_drawn_as_bars(
    render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """Histogram bin centres are a coordinate, not a category, so they get a line."""
    render(SyntheticSettings(kind="profile", n_index=12, integer_index=False))

    assert not rendered["synthetic_occupancy_comparison.png"].axes[0].containers


def test_the_time_series_panel_draws_every_replicate_and_the_mean(
    render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """Three faint traces plus one mean line per condition, and no footnote to claim."""
    render(SyntheticSettings(kind="mean_of_timeseries"))

    figure = rendered["synthetic_size_timeseries.png"]
    lines = figure.axes[0].lines

    assert len(lines) == 8
    assert [line.get_label() for line in lines if line.get_label().startswith("A")] == ["A (n = 3)"]
    assert "Error bars" not in _footnotes(figure)


def test_plot_settings_turn_the_replicate_overlay_off() -> None:
    """The default settings model carries the three fields plugins may override."""
    settings = ContractPlotSettings(show_replicates=False, max_categories_for_bars=8)

    assert settings.error_bar == "ci95"
    assert settings.figsize == (10.0, 6.0)
    assert (settings.show_replicates, settings.max_categories_for_bars) == (False, 8)


def test_a_contract_plugin_gets_the_default_plot_settings_model() -> None:
    """A plugin that declares nothing still exposes ContractPlotSettings."""
    from polyzymd.analyses.rg_contract import Rg2Analysis

    assert Rg2Analysis.PlotSettingsModel is ContractPlotSettings
    assert SyntheticAnalysis.PlotSettingsModel is ContractPlotSettings


def test_rg2_writes_its_figures_through_run_comparison(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A real plugin gets figures from the full comparison lifecycle, unedited."""
    from types import SimpleNamespace

    from polyzymd.analyses._framework.contexts import Condition
    from polyzymd.analyses._framework.lifecycle import AnalysisLifecycle
    from polyzymd.analyses.rg_contract import Rg2Analysis, RgSettings
    from polyzymd.config.comparison import PlotSettings
    from tests.analyses.conftest import _stubbed

    settings = RgSettings(runs=[{"label": "protein", "selection": "all"}])
    monkeypatch.setattr(
        Condition,
        "from_condition_config",
        staticmethod(
            lambda cfg: Condition(
                cfg.label, cfg.config, tuple(cfg.replicates), make_simulation_config(cfg.label)
            )
        ),
    )
    config = SimpleNamespace(
        name="project",
        source_path=tmp_path / "comparison.yaml",
        defaults=SimpleNamespace(equilibration_time="0ns"),
        control="A",
        conditions=[
            SimpleNamespace(label=label, config=tmp_path / f"{label}.yaml", replicates=[1, 2, 3])
            for label in ("A", "B")
        ],
        plugins=SimpleNamespace(get=lambda name: None),
        plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
    )
    config.model_copy = lambda deep=True: config
    stubbed = _stubbed(
        Rg2Analysis,
        lambda replicate: make_synthetic_universe(scale=1.0 + 0.5 * replicate, n_frames=8),
        (),
    )

    result = AnalysisLifecycle(stubbed(), settings_resolver=lambda _a, _c: settings).run_comparison(
        config
    )

    assert [path.name for path in result["plots"]] == [
        "rg2_protein_comparison.png",
        "rg2_protein_timeseries.png",
    ]
    assert all(path.exists() for path in result["plots"])
