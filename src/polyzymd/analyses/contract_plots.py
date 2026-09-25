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
``fraction``            comparison bars with 1.0 marked as the bound
``profile``             one line per condition over the index with a band,
                        or grouped bars when the index is categorical
======================= =================================================

Every bar and every band is the interval named by the plugin's ``error_bar``
setting across replicates, never across frames, and every figure that draws one
carries the footnote saying what it is. Every mean, SEM and interval is read
from the condition-level aggregate, so a figure and the text report can never
disagree; the NPZ sidecar the runner writes beside each replicate artifact is
read only for the faint per-replicate traces and the per-frame panels.

A figure never hides a value. A fraction interval that reaches past 1.0 is
drawn in full with a dashed line marking the physical bound, and a condition
with no estimate leaves a hatched gap rather than a zero bar.

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
from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.shared.plotting import (
    annotate_uncertainty,
    apply_axis_style,
    apply_legend,
    get_condition_colors,
    get_output_path,
    grouped_bars,
    load_canonical_plot_artifacts,
    order_condition_labels,
    plugin_plot_settings,
    resolve_error_bar,
    save_figure,
    shared_count_half_widths,
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
    timeseries_figsize : tuple of float
        Width and height of the per-frame panels, which are usually wider.
    show_replicates : bool
        Draw the per-replicate points on bars and the per-replicate traces on
        lines. Turn it off for a condition with many replicates.
    max_categories_for_bars : int
        Longest categorical profile still drawn as grouped bars. A longer one,
        or one indexed by a non-integer coordinate such as a histogram bin
        centre, is drawn as a line.
    """

    figsize: tuple[float, float] = (10.0, 6.0)
    timeseries_figsize: tuple[float, float] = (12.0, 5.0)
    show_replicates: bool = True
    max_categories_for_bars: int = 30


@dataclass(frozen=True)
class _ConditionData:
    """One condition's aggregates plus the per-replicate series behind them."""

    label: str
    aggregates: dict[str, ObservableAggregate]
    series: dict[str, list[Any]]
    completeness: dict[str, Any] | None = None


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
    from matplotlib.container import BarContainer

    aggregates = [data.aggregates[name] for data in conditions]
    counts = [aggregate.n_replicates for aggregate in aggregates]
    points = bool(settings.show_replicates)
    fig, ax = plt.subplots(figsize=tuple(settings.figsize))
    grouped_bars(
        ax,
        np.arange(1),
        [
            (data.label, [aggregate.mean if aggregate.mean is not None else 0.0], [aggregate.sem])
            for data, aggregate in zip(conditions, aggregates, strict=True)
        ],
        colors,
        ctx.plot_settings,
        bar_width=min(0.8 / max(len(aggregates), 1), 0.3),
        error_bar=error_bar,
        reference_line=None,
        replicate_values=(
            [[list(aggregate.replicate_values)] for aggregate in aggregates] if points else None
        ),
        n_replicates=counts,
    )
    ax.set_xlim(-0.5, 0.5)
    ax.set_xticks([0])
    missing = _hatch_missing(
        ax,
        [c for c in ax.containers if isinstance(c, BarContainer)],
        conditions,
        aggregates,
        analysis_name,
    )
    ax.set_xticklabels([name if not missing else f"{name} (n/a: {', '.join(missing)})"])
    apply_axis_style(ax, ctx.plot_settings, title=name, ylabel=_axis_label(aggregates[0]))
    if aggregates[0].kind == "fraction":
        _mark_fraction_bound(ax)
    apply_legend(ax, ctx.plot_settings)
    _footnote(fig, ctx, analysis_name, counts, points=points, conditions=conditions)
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

    fig, ax = plt.subplots(figsize=tuple(settings.timeseries_figsize))
    drawn = False
    for data, color in zip(conditions, colors, strict=True):
        matrix = _stacked(data.series.get(name))
        if matrix is None:
            logger.warning(
                "%s: condition %r has no %r sidecar series, leaving it off the panel",
                analysis_name,
                data.label,
                name,
            )
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
    _partial_footnote(fig, conditions)
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

    aggregates = [data.aggregates[name] for data in conditions]
    index, columns = _aligned_index(analysis_name, name, conditions, aggregates)
    counts = [aggregate.n_replicates for aggregate in aggregates]
    fig, ax = plt.subplots(figsize=tuple(settings.figsize))
    categorical = _is_categorical(index, int(settings.max_categories_for_bars))
    if categorical:
        points = _profile_bars(
            ax, ctx, settings, error_bar, index, columns, conditions, aggregates, colors
        )
    else:
        points = _profile_lines(
            ax, settings, error_bar, index, columns, name, conditions, aggregates, colors
        )
    apply_axis_style(
        ax,
        ctx.plot_settings,
        title=name,
        xlabel=aggregates[0].index_label or ("Category" if categorical else "Index"),
        ylabel=_axis_label(aggregates[0]),
    )
    if aggregates[0].kind == "fraction":
        _mark_fraction_bound(ax)
    apply_legend(ax, ctx.plot_settings)
    _footnote(fig, ctx, analysis_name, counts, points=points, conditions=conditions)
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
    columns: Sequence[np.ndarray],
    conditions: Sequence[_ConditionData],
    aggregates: Sequence[ObservableAggregate],
    colors: Sequence[Any],
) -> bool:
    """Draw a short categorical profile as bars, saying whether it drew points."""
    import numpy as np

    replicate_values: list[Any] = []
    for data, aggregate, column in zip(conditions, aggregates, columns, strict=True):
        matrix = _replicate_profiles(data, aggregate, column)
        if matrix is None:
            replicate_values = []
            break
        replicate_values.append([matrix[:, position].tolist() for position in range(index.size)])
    points = bool(replicate_values and settings.show_replicates)
    grouped_bars(
        ax,
        np.arange(index.size),
        [
            (data.label, _at(aggregate.profile_mean, column), _at(aggregate.profile_sem, column))
            for data, aggregate, column in zip(conditions, aggregates, columns, strict=True)
        ],
        colors,
        ctx.plot_settings,
        error_bar=error_bar,
        reference_line=None,
        replicate_values=replicate_values if points else None,
        n_replicates=[aggregate.n_replicates for aggregate in aggregates],
    )
    ax.set_xticks(np.arange(index.size))
    ax.set_xticklabels([f"{value:g}" for value in index], rotation=90, fontsize=7)
    return points


def _profile_lines(
    ax: Any,
    settings: Any,
    error_bar: str,
    index: np.ndarray,
    columns: Sequence[np.ndarray],
    name: str,
    conditions: Sequence[_ConditionData],
    aggregates: Sequence[ObservableAggregate],
    colors: Sequence[Any],
) -> bool:
    """Draw a long or continuous profile as a line and a band per condition.

    The band is the interval on the aggregate's own ``profile_sem``, so it says
    the same thing as the text report whether or not the sidecars are present.
    The sidecars supply only the faint replicate traces.
    """
    import numpy as np

    drew_traces = False
    for data, aggregate, column, color in zip(conditions, aggregates, columns, colors, strict=True):
        mean = np.asarray(_at(aggregate.profile_mean, column), dtype=float)
        matrix = _replicate_profiles(data, aggregate, column)
        if settings.show_replicates and matrix is None:
            logger.warning(
                "condition %r has no %r sidecar profiles, drawing no replicate traces",
                data.label,
                name,
            )
        elif settings.show_replicates:
            drew_traces = True
            for trace in matrix:
                ax.plot(index, trace, color=color, alpha=0.2, linewidth=0.7)
        half = shared_count_half_widths(
            np.asarray(_at(aggregate.profile_sem, column), dtype=float),
            aggregate.n_replicates,
            error_bar=error_bar,
        )
        if half.size == mean.size and bool(np.any(half > 0.0)):
            ax.fill_between(index, mean - half, mean + half, color=color, alpha=0.2, linewidth=0)
        ax.plot(
            index,
            mean,
            color=color,
            linewidth=1.8,
            label=f"{data.label} (n = {aggregate.n_replicates})",
        )
    return drew_traces


def _footnote(
    fig: Any,
    ctx: PlotContext,
    analysis_name: str,
    counts: Sequence[int],
    *,
    points: bool,
    conditions: Sequence[_ConditionData] = (),
) -> None:
    """Say what the figure's interval is, and whether its data are partial.

    The interval sentence is left out when one replicate leaves none to draw;
    the partial note never is.
    """
    partial = _partial_text(conditions)
    if min(counts, default=0) < 2:
        _partial_footnote(fig, conditions)
        return
    annotate_uncertainty(
        fig,
        ctx.plot_settings,
        analysis_name,
        n_replicates=min(counts),
        equilibration=ctx.equilibration,
        points=points,
        partial=partial,
    )


def _partial_text(conditions: Sequence[_ConditionData]) -> str | None:
    """``Partial: ...`` naming what each partial condition in a figure is missing."""
    from polyzymd.analyses.completeness import summaries

    lines = summaries({data.label: data.completeness for data in conditions})
    return f"Partial: {'; '.join(lines)}." if lines else None


def _partial_footnote(fig: Any, conditions: Sequence[_ConditionData]) -> None:
    """Write the partial note alone, on a figure with no interval sentence."""
    text = _partial_text(conditions)
    if text:
        fig.text(0.01, 0.01, text, fontsize=7, color="dimgray", ha="left", va="bottom")


def _mark_fraction_bound(ax: Any) -> None:
    """Mark 1.0 on a fraction axis without cutting an interval that passes it."""
    bottom, top = ax.get_ylim()
    ax.set_ylim(min(0.0, float(bottom)), max(1.05, float(top)))
    ax.axhline(
        1.0, color="dimgray", linestyle="--", linewidth=1.0, label="Physical bound (fraction = 1)"
    )


def _hatch_missing(
    ax: Any,
    containers: Sequence[Any],
    conditions: Sequence[_ConditionData],
    aggregates: Sequence[ObservableAggregate],
    analysis_name: str,
) -> list[str]:
    """Turn the bar of a condition with no estimate into a hatched gap."""
    missing = []
    stub = 0.02 * float(ax.get_ylim()[1])
    for container, data, aggregate in zip(containers, conditions, aggregates, strict=False):
        if aggregate.mean is not None:
            continue
        missing.append(data.label)
        for patch in container.patches:
            patch.set_height(stub)
            patch.set_facecolor("none")
            patch.set_edgecolor("dimgray")
            patch.set_hatch("///")
    if missing:
        logger.warning(
            "%s: conditions %s have no estimate for %r, drawing them as gaps",
            analysis_name,
            missing,
            aggregates[0].name,
        )
    return missing


def _aligned_index(
    analysis_name: str,
    name: str,
    conditions: Sequence[_ConditionData],
    aggregates: Sequence[ObservableAggregate],
) -> tuple[np.ndarray, list[np.ndarray]]:
    """Index entries every condition reports, and where each one holds them.

    Conditions may disagree on the index when a residue or a pair is absent
    from one system. The figure keeps the intersection, which is the part every
    condition actually measured, and names the dropped entries in a warning.

    Raises
    ------
    PluginContractError
        If a condition reports the profile without an index, or if the
        conditions share no index entry at all.
    """
    import numpy as np

    indices = []
    for data, aggregate in zip(conditions, aggregates, strict=True):
        if not aggregate.index:
            raise PluginContractError(
                f"{analysis_name}: profile observable {name!r} has no index for condition "
                f"{data.label!r}; every condition must report the same kind of index"
            )
        indices.append(np.asarray(aggregate.index, dtype=float))

    shared = indices[0]
    for other in indices[1:]:
        shared = np.intersect1d(shared, other)
    if shared.size == 0:
        raise PluginContractError(
            f"{analysis_name}: profile observable {name!r} has no index entry shared by the "
            f"conditions {[data.label for data in conditions]}; their indices are incompatible"
        )
    dropped = sorted(set(np.concatenate(indices).tolist()) - set(shared.tolist()))
    if dropped:
        logger.warning(
            "%s: profile %r keeps the %d index entries every condition reports and drops %s",
            analysis_name,
            name,
            int(shared.size),
            dropped,
        )
    return shared, [
        np.asarray([int(np.flatnonzero(index == value)[0]) for value in shared])
        for index in indices
    ]


def _at(values: Sequence[float] | None, column: np.ndarray) -> list[float]:
    """The entries of a per-index vector that survived the index alignment."""
    import numpy as np

    if not values:
        return []
    return np.asarray(values, dtype=float)[column].tolist()


def _replicate_profiles(
    data: _ConditionData, aggregate: ObservableAggregate, column: np.ndarray
) -> np.ndarray | None:
    """Per-replicate profiles from the sidecars, aligned to the shared index."""
    matrix = _stacked(data.series.get(aggregate.name))
    if matrix is None or matrix.shape[1] <= int(column.max()):
        return None
    return matrix[:, column]


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
    return _ConditionData(
        label=label,
        aggregates=aggregates,
        series=series,
        completeness=artifacts.condition_artifact.metadata.get("completeness"),
    )


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
