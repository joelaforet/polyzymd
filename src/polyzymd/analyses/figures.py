"""Draw figures of stored study results with matplotlib, without reading a trajectory.

:func:`plot_timeseries`, :func:`plot_distribution` and :func:`plot_values`
are reached as ``Timeseries.plot``, ``Timeseries.plot_distribution`` and
``ReplicateValues.plot``. They read the values already held by those objects,
style each figure through :mod:`polyzymd.analyses.shared.plotting` and a
:class:`~polyzymd.config.comparison.PlotSettings`, and save it with
:func:`~polyzymd.analyses.shared.plotting.save_figure`. Every condition and
every replicate is drawn. matplotlib is imported only when a figure is drawn.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from polyzymd.analyses.timeseries import ReplicateValues, Timeseries
    from polyzymd.config.comparison import PlotSettings


def _setup(source: Timeseries, plot_settings: PlotSettings | None) -> tuple[Any, list, dict]:
    """Return the plot settings, the condition labels in plot order and their colours."""
    from polyzymd.analyses.shared.plotting import get_condition_color_map, order_condition_labels
    from polyzymd.config.comparison import PlotSettings

    settings = plot_settings or PlotSettings()
    labels = list(source.series)
    colors = get_condition_color_map(labels, settings, control_label=source.study.control)
    return settings, order_condition_labels(labels, settings), colors


def _label(quantity: str, unit: str | None) -> str:
    """Write an axis label from a quantity name and its unit, with Å for ``A``."""
    text = quantity.replace("_", " ")
    return f"{text} ({'Å' if unit == 'A' else unit})" if unit else text


def _footnote(fig: Any, counts: set[int], equilibration: str, **wording: str) -> None:
    """Footnote the interval with n, or say that no interval is drawn with one replicate each."""
    from polyzymd.analyses.shared.plotting import add_uncertainty_footnote

    if max(counts) < 2:
        _note(
            fig,
            f"No interval: every condition has one replicate; production window t >= {equilibration}.",
        )
        return
    n = counts.pop() if len(counts) == 1 else None
    add_uncertainty_footnote(fig, n_replicates=n, equilibration=equilibration, **wording)


def _note(fig: Any, text: str) -> None:
    """Write ``text`` at the bottom left of ``fig`` in the footnote style."""
    fig.text(0.01, 0.01, text, fontsize=7, color="dimgray", ha="left", va="bottom")


def _save(fig: Any, output_dir: Path, name: str, settings: Any) -> Path:
    """Save ``fig`` as ``<output_dir>/<name>.<format>`` and close it."""
    from polyzymd.analyses.shared.plotting import get_output_path, save_figure
    from polyzymd.analyses.timeseries import _safe

    return save_figure(fig, get_output_path(Path(output_dir), _safe(name), settings), settings)


def plot_timeseries(
    source: Timeseries,
    output_dir: str | Path,
    name: str,
    plot_settings: PlotSettings | None = None,
) -> Path:
    """Draw every replicate's series against simulation time, coloured by condition.

    Each replicate is a thin line. Where every replicate of a condition has
    the same times over its first frames, a thick line gives their mean over
    those frames, with a band of the 95 percent Student t interval across
    replicates at each time, as the legacy rg and rmsd plugins drew it. The
    equilibration window, which is not stored, is shaded from the start of
    each trajectory, found from the frame indices and times, to the end of
    the window.

    Returns
    -------
    Path
        The saved figure file.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    from polyzymd.analyses.shared.loader import convert_time, parse_time_string
    from polyzymd.analyses.shared.plotting import apply_axis_style, apply_legend, band_half_widths

    settings, labels, colors = _setup(source, plot_settings)
    equilibration = source.study[labels[0]].equilibration
    window_ns = convert_time(*parse_time_string(equilibration), "ns")
    fig, ax = plt.subplots(figsize=(12, 5))
    starts, counts = [], set()
    for label in labels:
        items, color = source.series[label], colors[label]
        counts.add(len(items))
        for item in items:
            ax.plot(item.times, item.values, color=color, linewidth=0.8, alpha=0.35, zorder=1)
            if len(item.frames) > 1:
                step = (item.times[1] - item.times[0]) / (item.frames[1] - item.frames[0])
                starts.append(item.times[0] - item.frames[0] * step)
        common = min(len(item.values) for item in items)
        times = items[0].times[:common]
        if all(np.allclose(item.times[:common], times) for item in items):
            matrix = np.vstack([item.values[:common] for item in items])
            mean, band = matrix.mean(axis=0), band_half_widths(matrix)
            ax.plot(
                times,
                mean,
                color=color,
                linewidth=2.0,
                zorder=3,
                label=f"{label} (n = {len(items)})",
            )
            if band is not None:
                ax.fill_between(times, mean - band, mean + band, color=color, alpha=0.2, zorder=2)
        else:
            ax.plot(
                [],
                [],
                color=color,
                linewidth=0.8,
                label=f"{label} (n = {len(items)}, no common times)",
            )
    if window_ns > 0 and starts:
        start = min(starts)
        ax.axvspan(start, start + window_ns, color="0.85", zorder=0, label="equilibration window")
    apply_axis_style(
        ax, settings, title=source.name, xlabel="Time (ns)", ylabel=_label(source.name, source.unit)
    )
    apply_legend(ax, settings, loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
    fig.tight_layout(rect=[0, 0.04, 0.78, 1])
    _footnote(
        fig,
        counts,
        equilibration,
        drawn="Band",
        points="Thin lines are per-replicate series; the thick line is their mean",
    )
    return _save(fig, output_dir, name, settings)


def _density(ax: Any, values: Any, **style: Any) -> None:
    """Draw a Gaussian KDE of ``values``, or a vertical line when they are all equal.

    The KDE is ``scipy.stats.gaussian_kde`` with Scott's bandwidth, drawn over
    the data range extended by three bandwidths, as seaborn's ``kdeplot``
    draws it.
    """
    import numpy as np
    from scipy.stats import gaussian_kde

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2 or np.allclose(values, values[0]):
        ax.axvline(values[0], **style)
        return
    kde = gaussian_kde(values)
    width = 3.0 * kde.factor * float(np.std(values, ddof=1))
    grid = np.linspace(values.min() - width, values.max() + width, 200)
    ax.plot(grid, kde(grid), **style)


def plot_distribution(
    source: Timeseries,
    output_dir: str | Path,
    name: str,
    threshold: float | None = None,
    title: str | None = None,
    plot_settings: PlotSettings | None = None,
) -> Path:
    """Draw the distribution of the per-frame values of each condition.

    For each condition, a thick line is the Gaussian KDE of every production
    frame of every replicate pooled, as the legacy distances and catalytic
    triad plugins drew it, and a thin line of the same colour is the KDE of
    each replicate alone, so the spread between replicates shows. A
    ``threshold`` is drawn as a red dashed vertical line. ``title``, by
    default the series name, also names the quantity on the x axis.

    Returns
    -------
    Path
        The saved figure file.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    from polyzymd.analyses.shared.plotting import apply_axis_style, apply_legend

    settings, labels, colors = _setup(source, plot_settings)
    fig, ax = plt.subplots(figsize=(10, 6))
    for label in labels:
        items, color = source.series[label], colors[label]
        for item in items:
            _density(ax, item.values, color=color, linewidth=0.8, alpha=0.4)
        pooled = np.concatenate([item.values for item in items])
        _density(ax, pooled, color=color, linewidth=2.0, label=f"{label} (n = {len(items)})")
    if threshold is not None:
        ax.axvline(
            threshold,
            color="red",
            linestyle=settings.theme.reference_line_style,
            linewidth=settings.theme.reference_line_width,
            label=f"threshold {threshold:g} {'Å' if source.unit == 'A' else source.unit or ''}",
        )
    title = title or source.name
    apply_axis_style(ax, settings, title=title, xlabel=_label(title, source.unit), ylabel="Density")
    apply_legend(ax, settings, loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
    fig.tight_layout(rect=[0, 0.04, 0.78, 1])
    equilibration = source.study[labels[0]].equilibration
    _note(
        fig,
        "Thick lines: Gaussian KDE of all replicates' frames pooled; thin lines: each "
        f"replicate; n = replicates; production window t >= {equilibration}.",
    )
    return _save(fig, output_dir, name, settings)


def plot_values(
    values: ReplicateValues,
    output_dir: str | Path,
    name: str,
    title: str | None = None,
    plot_settings: PlotSettings | None = None,
) -> Path:
    """Draw one bar per condition at the mean, with its 95 percent interval and every replicate.

    The bar height is the mean and the error bar the Student t interval of
    :func:`~polyzymd.analyses.shared.statistics.mean_sem_ci`, the values that
    ``summary()`` reports. Every replicate value is a point on its bar, as
    Grossfield et al. (2018) recommend for fewer than 10 samples. Each tick
    label gives the condition's number of replicates, and a condition with
    one replicate gets no error bar.

    Returns
    -------
    Path
        The saved figure file.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    from polyzymd.analyses.shared.plotting import (
        apply_axis_style,
        error_bar_half_widths,
        scatter_replicate_values,
    )
    from polyzymd.analyses.shared.statistics import mean_sem_ci

    settings, labels, colors = _setup(values.source, plot_settings)
    data = [values.values[label] for label in labels]
    stats = [mean_sem_ci(entry) for entry in data]
    positions, theme = np.arange(len(labels)), settings.theme
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(
        positions,
        [item.mean for item in stats],
        yerr=error_bar_half_widths([item.sem for item in stats], data),
        color=[colors[label] for label in labels],
        edgecolor=theme.bar_edgecolor,
        linewidth=theme.bar_linewidth,
        capsize=theme.bar_capsize,
        alpha=theme.bar_alpha,
    )
    scatter_replicate_values(ax, positions, data, settings)
    ax.set_xticks(positions)
    ax.set_xticklabels(
        [f"{label}\nn = {len(entry)}" for label, entry in zip(labels, data, strict=True)],
        rotation=30,
        ha="right",
    )
    if values.is_fraction:
        ax.set_ylim(0, 1.05)
    apply_axis_style(
        ax, settings, title=title or values.source.name, ylabel=_label(values.metric, values.unit)
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    _footnote(fig, {len(entry) for entry in data}, values.source.study[labels[0]].equilibration)
    return _save(fig, output_dir, name, settings)
