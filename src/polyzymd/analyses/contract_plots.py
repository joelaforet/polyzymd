"""Generic comparison figures for contract plugins, keyed on observable kind.

A plugin under the observable contract declares a kind for every quantity it
reports, and the kind already fixes how the framework reduces and compares the
numbers. It fixes the figure too, so a plugin gets its plots from here instead
of carrying a ``_plotters`` module.

Each kind maps to a fixed set of figures.

======================= =================================================
kind                    figures
======================= =================================================
``mean_of_timeseries``  comparison bars across conditions, time series
``fluctuation``         comparison bars across conditions, time series
``fraction``            comparison bars clamped to [0, 1]
``profile``             one line per condition over the index with a band,
                        or grouped bars when the index is categorical
======================= =================================================

Every bar and every band is the interval named by the plugin's ``error_bar``
setting across replicates, never across frames, and every figure that draws one
carries the footnote saying what it is. Bars and bands come from the
condition-level aggregate; replicate traces and per-frame series come from the
NPZ sidecar the runner writes beside each replicate artifact.

References
----------
Grossfield, A., Patrone, P. N., Roe, D. R., Schultz, A. J., Siderius, D. W. &
Zuckerman, D. M. (2018). Best practices for quantifying the uncertainty in
molecular simulations. *Living Journal of Computational Molecular Science*,
1(1), 5067. doi:10.33011/livecoms.1.1.5067
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterator, Sequence

from polyzymd.analyses._framework.comparison_models import BasePlotSettings
from polyzymd.analyses.contract import ObservableAggregate
from polyzymd.analyses.shared.plotting import (
    annotate_uncertainty,
    apply_axis_style,
    apply_legend,
    band_half_widths,
    get_condition_colors,
    get_output_path,
    grouped_bars,
    load_canonical_plot_artifacts,
    order_condition_labels,
    plugin_plot_settings,
    resolve_error_bar,
    save_figure,
)

if TYPE_CHECKING:
    import numpy as np

    from polyzymd.analyses.base import PlotContext

logger = logging.getLogger(__name__)

#: Kinds whose per-frame series is worth a time-series panel of its own.
_TIMESERIES_KINDS = ("mean_of_timeseries", "fluctuation")


class ContractPlotSettings(BasePlotSettings):
    """Plot settings every contract plugin gets, and any plugin may extend.

    Attributes
    ----------
    figsize : tuple of float
        Width and height of every generated figure, in inches.
    show_replicates : bool
        Draw the per-replicate points on bars and the per-replicate traces on
        lines. Turn it off for a condition with many replicates.
    max_categories_for_bars : int
        Longest categorical profile still drawn as grouped bars. A longer one,
        or one indexed by a non-integer coordinate such as a histogram bin
        centre, is drawn as a line.
    """

    figsize: tuple[float, float] = (10.0, 6.0)
    show_replicates: bool = True
    max_categories_for_bars: int = 30


@dataclass(frozen=True)
class _ConditionData:
    """One condition's aggregates plus the per-replicate series behind them."""

    label: str
    aggregates: dict[str, ObservableAggregate]
    series: dict[str, list[Any]]


def plot_observables(analysis_name: str, ctx: PlotContext) -> list[Path]:
    """Render every figure the observables of one analysis call for.

    Parameters
    ----------
    analysis_name : str
        Name of the analysis, used for the plugin plot settings and as the
        first part of every file name.
    ctx : PlotContext
        Framework plot context. ``ctx.output_dir`` is the figures directory.

    Returns
    -------
    list[Path]
        Paths of the figures written, in the order the plugin reported its
        observables.
    """
    conditions = [
        data
        for data in (_load_condition(ctx, label) for label in _ordered_labels(ctx))
        if data is not None
    ]
    if not conditions:
        logger.warning("%s: no condition aggregate on disk, skipping plots", analysis_name)
        return []

    settings = plugin_plot_settings(ctx.plot_settings, analysis_name) or ContractPlotSettings()
    error_bar = resolve_error_bar(settings, ctx.plot_settings)
    labels = [data.label for data in conditions]
    palette = dict(
        zip(
            labels,
            get_condition_colors(labels, ctx.plot_settings, control_label=ctx.control_label),
            strict=True,
        )
    )
    paths: list[Path] = []
    for name, kind in _observable_order(conditions):
        present = [data for data in conditions if name in data.aggregates]
        shown = [palette[data.label] for data in present]
        if kind == "profile":
            paths.append(
                _profile_figure(analysis_name, ctx, settings, error_bar, name, present, shown)
            )
            continue
        paths.append(_bar_figure(analysis_name, ctx, settings, error_bar, name, present, shown))
        if kind in _TIMESERIES_KINDS:
            series = _timeseries_figure(analysis_name, ctx, settings, name, present, shown)
            if series is not None:
                paths.append(series)
    return paths


def _bar_figure(
    analysis_name: str,
    ctx: PlotContext,
    settings: Any,
    error_bar: str,
    name: str,
    conditions: Sequence[_ConditionData],
    colors: Sequence[Any],
) -> Path:
    """One bar per condition for a scalar observable, with replicate points."""
    import matplotlib.pyplot as plt
    import numpy as np

    aggregates = [data.aggregates[name] for data in conditions]
    replicate_values = [[list(aggregate.replicate_values)] for aggregate in aggregates]
    fig, ax = plt.subplots(figsize=tuple(settings.figsize))
    grouped_bars(
        ax,
        np.arange(1),
        [
            (data.label, [aggregate.mean or 0.0], [aggregate.sem])
            for data, aggregate in zip(conditions, aggregates, strict=True)
        ],
        colors,
        ctx.plot_settings,
        error_bar=error_bar,
        reference_line=None,
        replicate_values=replicate_values if settings.show_replicates else None,
        n_replicates=min(aggregate.n_replicates for aggregate in aggregates),
    )
    ax.set_xticks([0])
    ax.set_xticklabels([name])
    apply_axis_style(ax, ctx.plot_settings, title=name, ylabel=_axis_label(aggregates[0]))
    if aggregates[0].kind == "fraction":
        ax.set_ylim(0.0, 1.0)
    apply_legend(ax, ctx.plot_settings)
    _footnote(fig, ctx, analysis_name, aggregates)
    return save_figure(
        fig,
        get_output_path(ctx.output_dir, f"{analysis_name}_{name}_comparison", ctx.plot_settings),
        ctx.plot_settings,
    )


def _timeseries_figure(
    analysis_name: str,
    ctx: PlotContext,
    settings: Any,
    name: str,
    conditions: Sequence[_ConditionData],
    colors: Sequence[Any],
) -> Path | None:
    """One faint line per replicate plus the condition mean, over frame index."""
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=tuple(settings.figsize))
    drawn = False
    for data, color in zip(conditions, colors, strict=True):
        matrix = _stacked(data.series.get(name))
        if matrix is None:
            continue
        drawn = True
        frames = np.arange(matrix.shape[1])
        if settings.show_replicates:
            for trace in matrix:
                ax.plot(frames, trace, color=color, alpha=0.25, linewidth=0.8)
        ax.plot(
            frames,
            matrix.mean(axis=0),
            color=color,
            linewidth=1.8,
            label=f"{data.label} (n = {matrix.shape[0]})",
        )
    if not drawn:
        plt.close(fig)
        logger.warning("%s: no per-frame series for %r, skipping its panel", analysis_name, name)
        return None
    apply_axis_style(
        ax,
        ctx.plot_settings,
        title=f"{name} per frame, production window t >= {ctx.equilibration}",
        xlabel="Frame in production window",
        ylabel=_axis_label(conditions[0].aggregates[name]),
    )
    apply_legend(ax, ctx.plot_settings)
    return save_figure(
        fig,
        get_output_path(ctx.output_dir, f"{analysis_name}_{name}_timeseries", ctx.plot_settings),
        ctx.plot_settings,
    )


def _profile_figure(
    analysis_name: str,
    ctx: PlotContext,
    settings: Any,
    error_bar: str,
    name: str,
    conditions: Sequence[_ConditionData],
    colors: Sequence[Any],
) -> Path:
    """A profile as lines with a band, or as grouped bars over few categories."""
    import matplotlib.pyplot as plt
    import numpy as np

    aggregates = [data.aggregates[name] for data in conditions]
    index = np.asarray(aggregates[0].index or [], dtype=float)
    fig, ax = plt.subplots(figsize=tuple(settings.figsize))
    categorical = _is_categorical(index, int(settings.max_categories_for_bars))
    if categorical:
        _profile_bars(ax, ctx, settings, error_bar, index, conditions, aggregates, colors)
    else:
        _profile_lines(ax, ctx, settings, error_bar, index, name, conditions, aggregates, colors)
    apply_axis_style(
        ax,
        ctx.plot_settings,
        title=name,
        xlabel="Index" if not categorical else "Category",
        ylabel=_axis_label(aggregates[0]),
    )
    if aggregates[0].kind == "fraction":
        ax.set_ylim(0.0, 1.0)
    apply_legend(ax, ctx.plot_settings)
    _footnote(fig, ctx, analysis_name, aggregates)
    return save_figure(
        fig,
        get_output_path(ctx.output_dir, f"{analysis_name}_{name}_comparison", ctx.plot_settings),
        ctx.plot_settings,
    )


def _profile_bars(
    ax: Any,
    ctx: PlotContext,
    settings: Any,
    error_bar: str,
    index: np.ndarray,
    conditions: Sequence[_ConditionData],
    aggregates: Sequence[ObservableAggregate],
    colors: Sequence[Any],
) -> None:
    """Draw a short categorical profile as one bar group per category."""
    import numpy as np

    replicate_values = []
    for data, aggregate in zip(conditions, aggregates, strict=True):
        matrix = _stacked(data.series.get(aggregate.name))
        if matrix is None or matrix.shape[1] != index.size:
            replicate_values = []
            break
        replicate_values.append([matrix[:, column].tolist() for column in range(index.size)])
    grouped_bars(
        ax,
        np.arange(index.size),
        [
            (data.label, list(aggregate.profile_mean or []), list(aggregate.profile_sem or []))
            for data, aggregate in zip(conditions, aggregates, strict=True)
        ],
        colors,
        ctx.plot_settings,
        error_bar=error_bar,
        reference_line=None,
        replicate_values=(
            replicate_values if replicate_values and settings.show_replicates else None
        ),
        n_replicates=min(aggregate.n_replicates for aggregate in aggregates),
    )
    ax.set_xticks(np.arange(index.size))
    ax.set_xticklabels([f"{value:g}" for value in index], rotation=90, fontsize=7)


def _profile_lines(
    ax: Any,
    ctx: PlotContext,
    settings: Any,
    error_bar: str,
    index: np.ndarray,
    name: str,
    conditions: Sequence[_ConditionData],
    aggregates: Sequence[ObservableAggregate],
    colors: Sequence[Any],
) -> None:
    """Draw a long or continuous profile as one line per condition with a band."""
    import numpy as np

    for data, aggregate, color in zip(conditions, aggregates, colors, strict=True):
        mean = np.asarray(aggregate.profile_mean or [], dtype=float)
        matrix = _stacked(data.series.get(name))
        if settings.show_replicates and matrix is not None and matrix.shape[1] == index.size:
            for trace in matrix:
                ax.plot(index, trace, color=color, alpha=0.2, linewidth=0.7)
        half = None if matrix is None else band_half_widths(matrix, error_bar=error_bar)
        if half is not None and half.size == mean.size:
            ax.fill_between(index, mean - half, mean + half, color=color, alpha=0.2, linewidth=0)
        ax.plot(
            index,
            mean,
            color=color,
            linewidth=1.8,
            label=f"{data.label} (n = {aggregate.n_replicates})",
        )


def _footnote(
    fig: Any, ctx: PlotContext, analysis_name: str, aggregates: Sequence[ObservableAggregate]
) -> None:
    """Say what the figure's interval is, unless one replicate left none to draw."""
    n_replicates = min(aggregate.n_replicates for aggregate in aggregates)
    if n_replicates < 2:
        return
    annotate_uncertainty(
        fig,
        ctx.plot_settings,
        analysis_name,
        n_replicates=n_replicates,
        equilibration=ctx.equilibration,
    )


def _ordered_labels(ctx: PlotContext) -> list[str]:
    """Condition labels in plot order, keeping only those with a directory."""
    labels = [
        condition.label for condition in ctx.conditions if condition.label in ctx.analysis_dirs
    ]
    return order_condition_labels(labels, ctx.plot_settings)


def _load_condition(ctx: PlotContext, label: str) -> _ConditionData | None:
    """Read one condition's aggregate and NPZ sidecars, or None when absent."""
    replicates = next(
        (list(condition.replicates) for condition in ctx.conditions if condition.label == label), []
    )
    artifacts = load_canonical_plot_artifacts(
        ctx.analysis_dirs[label], replicates, require_replicates=False
    )
    if artifacts.condition_artifact is None:
        return None
    aggregates = {
        aggregate.name: aggregate
        for aggregate in (
            ObservableAggregate.model_validate(payload)
            for payload in artifacts.condition_artifact.payload.get("observables", [])
        )
    }
    series: dict[str, list[Any]] = {}
    for replicate, artifact in sorted(artifacts.replicate_artifacts.items()):
        for name, values in _sidecar_arrays(artifacts.run_dirs[replicate], artifact):
            series.setdefault(name, []).append(values)
    return _ConditionData(label=label, aggregates=aggregates, series=series)


def _sidecar_arrays(run_dir: Path, artifact: Any) -> Iterator[tuple[str, Any]]:
    """Yield every named array in the NPZ sidecars of one replicate artifact."""
    import numpy as np

    for reference in getattr(artifact, "sidecars", []):
        path = run_dir / getattr(reference, "path", "")
        if path.suffix != ".npz" or not path.exists():
            continue
        with np.load(path) as payload:
            for name in payload.files:
                yield name, np.asarray(payload[name], dtype=float)


def _observable_order(conditions: Sequence[_ConditionData]) -> list[tuple[str, str]]:
    """Observable names and kinds, in the order the first condition reports them."""
    ordered: dict[str, str] = {}
    for data in conditions:
        for name, aggregate in data.aggregates.items():
            ordered.setdefault(name, aggregate.kind)
    return list(ordered.items())


def _stacked(traces: Sequence[Any] | None) -> np.ndarray | None:
    """Replicate traces as one (n_replicates, n_points) array, or None."""
    import numpy as np

    if not traces:
        return None
    length = min(int(np.asarray(trace).size) for trace in traces)
    if length == 0:
        return None
    return np.asarray([np.asarray(trace, dtype=float)[:length] for trace in traces])


def _is_categorical(index: np.ndarray, limit: int) -> bool:
    """Whether a profile index names discrete items few enough to draw as bars."""
    import numpy as np

    return bool(index.size and index.size <= limit and np.all(index == np.round(index)))


def _axis_label(aggregate: ObservableAggregate) -> str:
    """Observable name with its unit, or bare when the quantity has none."""
    return aggregate.name if not aggregate.unit else f"{aggregate.name} ({aggregate.unit})"
