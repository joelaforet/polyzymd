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
    index_label: str | None = None
    high_occupancy: bool = False


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
                    index_label=settings.index_label,
                )
            ]
        if settings.kind == "fraction":
            filled = (
                n_frames - int(round((2.5 - scale) * 4.8))
                if settings.high_occupancy
                else min(n_frames, int(round(scale * 4)))
            )
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

    def _render(
        settings: SyntheticSettings,
        labels: Sequence[str] = ("A", "B"),
        mutate: Callable[[str, dict[str, Any]], None] | None = None,
    ) -> list[Path]:
        import json

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
            if mutate is not None:
                path = tmp_path / "analysis" / label / "synthetic" / "aggregated" / "result.json"
                payload = json.loads(path.read_text())
                mutate(label, payload["payload"]["observables"][0])
                path.write_text(json.dumps(payload))
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


def _bars(ax: Any) -> list[Any]:
    """Every bar container on an axes, in the order the series were drawn."""
    from matplotlib.container import BarContainer

    return [container for container in ax.containers if isinstance(container, BarContainer)]


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


def test_a_fraction_is_drawn_with_its_physical_bound_marked(
    render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """A fraction gets bars only, with 1.0 drawn as a line rather than as a clip."""
    paths = render(SyntheticSettings(kind="fraction"))

    assert [path.name for path in paths] == ["synthetic_bound_comparison.png"]
    axes = rendered["synthetic_bound_comparison.png"].axes[0]
    assert axes.get_ylim()[0] == 0.0
    assert axes.get_ylim()[1] >= 1.0
    assert axes.get_ylabel() == "bound (fraction)"
    assert "Physical bound (fraction = 1)" in [line.get_label() for line in axes.lines]


def test_a_fraction_interval_past_one_is_shown_whole(
    render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """An interval whose upper arm passes 1.0 is never cut off at the bound."""
    render(SyntheticSettings(kind="fraction", high_occupancy=True))

    axes = rendered["synthetic_bound_comparison.png"].axes[0]
    arms = _error_half_widths(axes)
    tops = [
        container.patches[0].get_height() + arm
        for container, arm in zip(_bars(axes), arms, strict=True)
    ]

    assert max(tops) > 1.0
    assert axes.get_ylim()[1] >= max(tops)


def test_a_categorical_profile_is_drawn_as_bars_up_to_thirty_categories(
    render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """Thirty residue labels still fit as grouped bars, one group per label."""
    render(SyntheticSettings(kind="profile", n_index=30))

    axes = rendered["synthetic_occupancy_comparison.png"].axes[0]

    assert [len(container) for container in _bars(axes)] == [30, 30]
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
    from polyzymd.analyses.rg import RgAnalysis

    assert RgAnalysis.PlotSettingsModel is ContractPlotSettings
    assert SyntheticAnalysis.PlotSettingsModel is ContractPlotSettings


def test_profiles_are_aligned_on_the_index_every_condition_reports(
    render: Callable[..., list[Path]], rendered: dict[str, Any], caplog: pytest.LogCaptureFixture
) -> None:
    """A residue missing from one condition is dropped, and the drop is logged."""

    def _drop_last(label: str, observable: dict[str, Any]) -> None:
        if label != "B":
            return
        for key in ("index", "profile_mean", "profile_sem"):
            observable[key] = observable[key][:-1]

    with caplog.at_level("WARNING", logger="polyzymd.analyses.contract_plots"):
        render(SyntheticSettings(kind="profile", n_index=10), mutate=_drop_last)

    axes = rendered["synthetic_occupancy_comparison.png"].axes[0]

    assert [len(container) for container in _bars(axes)] == [9, 9]
    assert "keeps the 9 index entries every condition reports and drops [9.0]" in caplog.text


def test_incompatible_profile_indices_raise_a_named_contract_error(
    render: Callable[..., list[Path]],
) -> None:
    """Conditions that share no index entry are an error, not a numpy traceback."""
    from polyzymd.analyses.exceptions import PluginContractError

    def _shift(label: str, observable: dict[str, Any]) -> None:
        if label == "B":
            observable["index"] = [value + 100.0 for value in observable["index"]]

    with pytest.raises(PluginContractError, match="occupancy.*incompatible"):
        render(SyntheticSettings(kind="profile", n_index=10), mutate=_shift)


def test_a_condition_with_no_estimate_leaves_a_gap(
    render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """A null mean is hatched out and named in the tick label, never drawn as zero."""

    def _blank(label: str, observable: dict[str, Any]) -> None:
        if label == "B":
            observable.update({"mean": None, "sem": None, "ci95_low": None, "ci95_high": None})

    render(SyntheticSettings(kind="mean_of_timeseries"), mutate=_blank)

    axes = rendered["synthetic_size_comparison.png"].axes[0]
    hatched = [
        patch for container in _bars(axes) for patch in container.patches if patch.get_hatch()
    ]

    assert len(hatched) == 1
    assert hatched[0].get_height() < 0.1 * float(axes.get_ylim()[1])
    assert axes.get_xticklabels()[0].get_text() == "size (n/a: B)"


def test_unequal_replicate_counts_get_their_own_t_factor(
    tmp_path: Path, render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """With no replicate points to count from, each bar still uses its own n."""

    def _thin(label: str, observable: dict[str, Any]) -> None:
        if label != "B":
            return
        observable["replicate_values"] = observable["replicate_values"][:2]
        observable["n_replicates"] = 2

    render(SyntheticSettings(kind="mean_of_timeseries"), mutate=_thin)

    half_widths = _error_half_widths(rendered["synthetic_size_comparison.png"].axes[0])
    factors = [
        half / _aggregate(tmp_path, label, "size").sem
        for half, label in zip(half_widths, ("A", "B"), strict=True)
    ]

    assert factors == pytest.approx([STUDENT_T_AT_THREE, 12.706204736432095], rel=1e-6)


def test_the_band_comes_from_the_aggregate_not_the_sidecars(
    render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """Zeroing the stored profile SEM removes that condition's band on its own.

    The replicate profiles in the sidecars still differ, so a band recomputed
    from them would survive. Only a band read from the aggregate disappears,
    which is what keeps the figure and the text report saying the same thing.
    """
    from matplotlib.collections import PolyCollection

    def _flatten(label: str, observable: dict[str, Any]) -> None:
        if label == "B":
            observable["profile_sem"] = [0.0] * len(observable["profile_sem"])

    render(SyntheticSettings(kind="profile", n_index=40))
    both = _band_count(rendered["synthetic_occupancy_comparison.png"], PolyCollection)
    rendered.clear()
    render(SyntheticSettings(kind="profile", n_index=40), mutate=_flatten)
    one = _band_count(rendered["synthetic_occupancy_comparison.png"], PolyCollection)

    assert (both, one) == (2, 1)


def _band_count(fig: Any, poly: type) -> int:
    """How many shaded bands the first axes of a figure carries."""
    return sum(isinstance(collection, poly) for collection in fig.axes[0].collections)


def test_a_profile_labels_its_axis_from_the_observable(
    render: Callable[..., list[Path]], rendered: dict[str, Any]
) -> None:
    """index_label reaches the figure, and "Index" stands in when it is absent."""
    render(SyntheticSettings(kind="profile", n_index=40, index_label="Residue"))
    assert rendered["synthetic_occupancy_comparison.png"].axes[0].get_xlabel() == "Residue"

    rendered.clear()
    render(SyntheticSettings(kind="profile", n_index=40))
    assert rendered["synthetic_occupancy_comparison.png"].axes[0].get_xlabel() == "Index"


def test_the_footnote_promises_points_only_when_points_are_drawn(
    tmp_path: Path, run_contract_analysis: Any
) -> None:
    """Turning the replicate overlay off drops the sentence about the points."""
    from polyzymd.analyses.shared.plotting import add_uncertainty_footnote

    figure = matplotlib.figure.Figure()

    assert "Points are per-replicate values." in add_uncertainty_footnote(
        figure, n_replicates=3, points=True
    )
    assert "Points" not in add_uncertainty_footnote(figure, n_replicates=3, points=False)


def test_rg_writes_its_figures_through_run_comparison(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A real plugin gets figures from the full comparison lifecycle, unedited."""
    from types import SimpleNamespace

    from polyzymd.analyses._framework.contexts import Condition
    from polyzymd.analyses._framework.lifecycle import AnalysisLifecycle
    from polyzymd.analyses.rg import RgAnalysis, RgSettings
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
        RgAnalysis,
        lambda replicate: make_synthetic_universe(scale=1.0 + 0.5 * replicate, n_frames=8),
        (),
    )

    result = AnalysisLifecycle(stubbed(), settings_resolver=lambda _a, _c: settings).run_comparison(
        config
    )

    assert [path.name for path in result["plots"]] == [
        "rg_rg_protein_comparison.png",
        "rg_rg_protein_timeseries.png",
    ]
    assert all(path.exists() for path in result["plots"])
