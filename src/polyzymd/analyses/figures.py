"""Draw figures of stored study results with matplotlib, without reading a trajectory.

:func:`plot_timeseries`, :func:`plot_distribution` and
:func:`plot_condition_values` are reached as ``Timeseries.plot``,
``Timeseries.plot_distribution`` and ``ReplicateValues.plot``, and
:func:`plot_profile` as ``ReplicateValues.plot`` of labelled values.
:func:`plot_values` and :func:`plot_distributions`, reached as
``pz.plot_values`` and ``pz.plot_distributions``, draw several results in
one figure. They read the values already held by those objects,
style each figure through :mod:`polyzymd.analyses.shared.plotting` and a
:class:`~polyzymd.config.analysis_settings.PlotSettings`, and save it with
:func:`~polyzymd.analyses.shared.plotting.save_figure`. Every condition and
every replicate is drawn. matplotlib is imported only when a figure is drawn.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

if TYPE_CHECKING:
    from polyzymd.analyses.timeseries import ReplicateValues, Timeseries
    from polyzymd.config.analysis_settings import PlotSettings


def _setup(source: Timeseries, plot_settings: PlotSettings | None) -> tuple[Any, list, dict]:
    """Return the plot settings, the condition labels in plot order and their colours."""
    from polyzymd.analyses.shared.plotting import get_condition_color_map, order_condition_labels
    from polyzymd.config.analysis_settings import PlotSettings

    settings = plot_settings or PlotSettings()
    labels = list(source.series)
    colors = get_condition_color_map(labels, settings, control_label=source.study.control)
    return settings, order_condition_labels(labels, settings), colors


def _label(quantity: str, unit: str | None) -> str:
    """Write an axis label from a quantity name and its unit, with Å for ``A``."""
    text = quantity.replace("_", " ")
    return f"{text} ({'Å' if unit == 'A' else unit})" if unit else text


def _footnote(
    fig: Any, replicates: set[int], equilibration: str | None, window: str, **wording: str
) -> None:
    """Footnote what the interval is and its n, or say that no interval is drawn.

    ``replicates`` holds the replicate counts of the conditions drawn. With one
    count the footnote gives that n, and with several it says that n is per
    condition. A condition with one replicate has no interval, and the
    footnote says so.
    """
    from polyzymd.analyses.shared.plotting import add_uncertainty_footnote

    if max(replicates) < 2:
        _note(fig, f"No interval: every condition has one replicate; {window}.")
        return
    if min(replicates) < 2:
        wording["points"] += "; a condition with one replicate has no interval"
    n = next(iter(replicates)) if len(replicates) == 1 else None
    add_uncertainty_footnote(fig, n_replicates=n, equilibration=equilibration, **wording)


def _note(fig: Any, text: str) -> None:
    """Write ``text`` under the axes of ``fig`` in the footnote style."""
    from polyzymd.analyses.shared.plotting import add_figure_note

    add_figure_note(fig, text)


#: How the footnote of :func:`plot_differences` names the test of each interval.
_TESTS = {"welch_t": "Welch t", "student_t": "Student t (pooled variance)"}


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
    """Draw every replicate's series against simulation time, colored by condition.

    Each replicate is a thin line. Where every replicate of a condition has
    the same times over its first frames, a thick line gives their mean over
    those frames, with a band of the 95 percent Student t interval across
    replicates at each time. The
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
    fig.tight_layout(rect=[0, 0, 0.78, 1])
    _footnote(
        fig,
        counts,
        equilibration,
        f"production window t >= {equilibration}",
        drawn="Band",
        of="the condition mean at each time",
        counts="n per condition in the legend",
        points="Thin lines: per-replicate series; thick lines: condition means; grey: the "
        "equilibration window, whose frames are not in the series",
    )
    return _save(fig, output_dir, name, settings)


def reflected_kde(values: Any, bounds: tuple = (None, None), points: int = 200) -> tuple:
    """Evaluate a Gaussian KDE of ``values`` that does not put density outside ``bounds``.

    ``f`` is ``scipy.stats.gaussian_kde`` with Scott's bandwidth ``h``. The
    grid runs from ``min(values) - 3h`` to ``max(values) + 3h``, as seaborn's
    ``kdeplot`` draws it, cut at each finite bound. At a lower bound ``a``
    the density is ``f(x) + f(2a - x)``, and at an upper bound ``b`` it gains
    ``f(2b - x)``, the reflection method of Schuster (1985) and Silverman
    (1986, section 2.10), so the curve integrates to 1 over the support and
    equals ``f`` far from a bound.

    Returns
    -------
    tuple of numpy.ndarray
        The grid and the density on it.
    """
    import numpy as np
    from scipy.stats import gaussian_kde

    values = np.asarray(values, dtype=float)
    kde = gaussian_kde(values)
    width = 3.0 * kde.factor * float(np.std(values, ddof=1))
    low, high = bounds
    start, stop = values.min() - width, values.max() + width
    grid = np.linspace(
        start if low is None else max(low, start), stop if high is None else min(high, stop), points
    )
    density = kde(grid)
    for bound in (low, high):
        if bound is not None:
            density += kde(2.0 * bound - grid)
    return grid, density


def _density(ax: Any, values: Any, bounds: tuple, **style: Any) -> None:
    """Draw :func:`reflected_kde` of ``values``, or a vertical line when they are all equal."""
    import numpy as np

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2 or np.allclose(values, values[0]):
        ax.axvline(values[0], **style)
        return
    ax.plot(*reflected_kde(values, bounds), **style)


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
    frame of every replicate pooled, and a thin line of the same colour is the KDE of
    each replicate alone, so the spread between replicates shows. A
    ``threshold`` is drawn as a red dashed vertical line. ``title``, by
    default the series name, also names the quantity on the x axis.

    Returns
    -------
    Path
        The saved figure file.
    """
    title = title or source.name
    return plot_distributions(
        [source], [threshold], [title], output_dir, name, title, plot_settings
    )


def _distribution_axes(ax: Any, source: Timeseries, threshold: float | None, style: tuple) -> None:
    """Draw the pooled and per-replicate KDEs of every condition and the threshold on ``ax``."""
    import numpy as np

    settings, labels, colors = style
    for label in labels:
        items, color = source.series[label], colors[label]
        for item in items:
            _density(ax, item.values, source.bounds, color=color, linewidth=0.8, alpha=0.4)
        pooled = np.concatenate([item.values for item in items])
        _density(
            ax,
            pooled,
            source.bounds,
            color=color,
            linewidth=2.0,
            label=f"{label} (n = {len(items)})",
        )
    if threshold is not None:
        ax.axvline(
            threshold,
            color="red",
            linestyle=settings.theme.reference_line_style,
            linewidth=settings.theme.reference_line_width,
            label=f"threshold {threshold:g} {'Å' if source.unit == 'A' else source.unit or ''}",
        )


def _one_unit(results: list, kind: str) -> str | None:
    """Return the unit shared by ``results``, refusing an empty list or mixed units."""
    from polyzymd.analyses.exceptions import ProtocolError

    units = {item.unit for item in results}
    if len(units) != 1:
        raise ProtocolError(
            f"{kind} in one figure must share one unit; got {sorted(map(str, units)) or 'none'}.",
            hint="Plot results with different units in separate figures.",
        )
    return units.pop()


def plot_distributions(
    series: Sequence[Timeseries],
    thresholds: Sequence[float | None] | None = None,
    titles: Sequence[str] | None = None,
    output_dir: str | Path | None = None,
    name: str = "distributions",
    quantity: str = "value",
    plot_settings: PlotSettings | None = None,
) -> Path:
    """Draw the distribution of each series in its own panel of one figure.

    Each panel shows, for every condition, the Gaussian KDE of every frame
    of every replicate pooled as a thick line and the KDE of each replicate
    as a thin line, with the series' threshold from ``thresholds`` as a red
    dashed line.
    The panels share the x axis, labelled ``quantity`` and the unit of the
    series. ``titles`` default to the series names, and ``output_dir`` to
    the ``figures`` folder next to ``polyzymd_results``.

    Returns
    -------
    Path
        The saved figure file.

    Raises
    ------
    ProtocolError
        If the series do not share one unit.
    """
    import matplotlib.pyplot as plt

    from polyzymd.analyses.shared.plotting import apply_axis_style, apply_legend

    series = list(series)
    unit = _one_unit(series, "Series")
    thresholds = list(thresholds) if thresholds is not None else [None] * len(series)
    titles = list(titles) if titles is not None else [item.name for item in series]
    style = _setup(series[0], plot_settings)
    settings, labels = style[0], style[1]
    fig, axes = plt.subplots(
        len(series),
        1,
        figsize=(10, 6 if len(series) == 1 else 3 * len(series)),
        sharex=True,
        squeeze=False,
    )
    for ax, source, threshold, title in zip(axes[:, 0], series, thresholds, titles, strict=True):
        _distribution_axes(ax, source, threshold, style)
        apply_axis_style(ax, settings, title=title, ylabel="Density")
        apply_legend(ax, settings, loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
    apply_axis_style(axes[-1, 0], settings, xlabel=_label(quantity, unit))
    fig.tight_layout(rect=[0, 0, 0.78, 1])
    _note(
        fig,
        "Thick lines: Gaussian KDE of all replicates' frames pooled; thin lines: each "
        f"replicate; n = replicates; production window t >= "
        f"{series[0].study[labels[0]].equilibration}."
        + (
            "\nGaussian KDE (Scott bandwidth), evaluated only within the physical support; "
            "near a bound the estimate is corrected by reflection (Schuster 1985; Silverman 1986)."
            if any(bound is not None for item in series for bound in item.bounds)
            else ""
        ),
    )
    folder = output_dir or series[0].path.parent.parent / "figures"
    return _save(fig, folder, name, settings)


def plot_condition_values(
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
    labels = [label for label in labels if values.values.get(label)]
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
    fig.tight_layout()
    equilibration = values.source.study[labels[0]].equilibration
    _footnote(
        fig,
        {len(entry) for entry in data},
        equilibration,
        f"production window t >= {equilibration}",
        of="the condition mean",
        counts="n per condition under each bar",
        points="Bars: condition means; points: per-replicate values",
    )
    return _save(fig, output_dir, name, settings)


def plot_profile(
    values: ReplicateValues,
    output_dir: str | Path,
    name: str,
    title: str | None = None,
    plot_settings: PlotSettings | None = None,
    highlight: Sequence = (),
    xlabel: str = "label",
) -> Path:
    """Draw a labelled result along its labels, one line per replicate and per condition.

    Each replicate's value at every label is a thin line. A thick line gives
    each condition's mean at every label, with a band of the 95 percent
    Student t interval across replicates, the values that ``summary()``
    reports.
    Numeric labels, such as residue IDs, are placed at their value on the x
    axis and other labels in order. Each ``highlight`` label is marked with a
    red dashed vertical line.

    Returns
    -------
    Path
        The saved figure file.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    from polyzymd.analyses.shared.plotting import apply_axis_style, apply_legend, band_half_widths

    settings, labels, colors = _setup(values.source, plot_settings)
    labels = [label for label in labels if values.values.get(label)]
    try:
        x = np.asarray(values.labels, dtype=float)
        ticks = None
    except (TypeError, ValueError):
        x, ticks = np.arange(len(values.labels), dtype=float), [str(k) for k in values.labels]
    fig, ax = plt.subplots(figsize=(14, 5))
    counts = set()
    for label in labels:
        matrix, color = np.vstack(values.values[label]), colors[label]
        counts.add(len(matrix))
        for row in matrix:
            ax.plot(x, row, color=color, linewidth=0.6, alpha=0.35, zorder=1)
        mean, band = matrix.mean(axis=0), band_half_widths(matrix)
        ax.plot(x, mean, color=color, linewidth=1.8, zorder=3, label=f"{label} (n = {len(matrix)})")
        if band is not None:
            ax.fill_between(x, mean - band, mean + band, color=color, alpha=0.2, zorder=2)
    where = {str(k): position for k, position in zip(values.labels, x, strict=True)}
    for key in highlight:
        if str(key) in where:
            ax.axvline(where[str(key)], color="red", linestyle="--", linewidth=1, alpha=0.6)
    if ticks is not None:
        ax.set_xticks(x)
        ax.set_xticklabels(ticks, rotation=90)
    apply_axis_style(
        ax,
        settings,
        title=title or values.metric,
        xlabel=xlabel,
        ylabel=_label(values.metric, values.unit),
    )
    apply_legend(ax, settings, loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
    fig.tight_layout(rect=[0, 0, 0.84, 1])
    equilibration = values.source.study[labels[0]].equilibration
    _footnote(
        fig,
        counts,
        equilibration,
        f"production window t >= {equilibration}",
        drawn="Band",
        of=f"the condition mean at each {xlabel.lower()}",
        counts="n per condition in the legend",
        points="Thin lines: per-replicate values; thick lines: condition means",
    )
    return _save(fig, output_dir, name, settings)


def plot_decomposition(
    parts: dict[str, ReplicateValues],
    output_dir: str | Path,
    name: str,
    title: str | None = None,
    plot_settings: PlotSettings | None = None,
    xlabel: str = "label",
) -> Path:
    """Draw several labelled results of one study together, one panel per condition.

    Each panel holds one line per result, its mean over the condition's
    replicates at every label, with the band of its 95 percent Student t
    interval, the values that ``summary()`` reports. For
    :func:`~polyzymd.analyses.functions.rms_decomposition` the lines are the
    RMS deviation from the reference, the RMSF and the offset of the mean
    position from the reference. The results must share one unit.

    Returns
    -------
    Path
        The saved figure file.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    from polyzymd.analyses.shared.plotting import (
        apply_axis_style,
        apply_legend,
        band_half_widths,
        get_palette_colors,
    )

    results = list(parts.values())
    unit = _one_unit(results, "Results")
    settings, labels, _ = _setup(results[0].source, plot_settings)
    colors = get_palette_colors(len(parts), settings)
    x = np.asarray(results[0].labels, dtype=float)
    fig, axes = plt.subplots(
        len(labels), 1, figsize=(14, 2.6 * len(labels) + 1), sharex=True, squeeze=False
    )
    counts = set()
    for ax, label in zip(axes[:, 0], labels, strict=True):
        for (part, values), color in zip(parts.items(), colors, strict=True):
            matrix = np.vstack(values.values[label])
            counts.add(len(matrix))
            mean, band = matrix.mean(axis=0), band_half_widths(matrix)
            ax.plot(x, mean, color=color, linewidth=1.5, label=part.replace("_", " "))
            if band is not None:
                ax.fill_between(x, mean - band, mean + band, color=color, alpha=0.2)
        apply_axis_style(
            ax, settings, title=f"{label} (n = {len(matrix)})", ylabel=_label("value", unit)
        )
        apply_legend(ax, settings, loc="center left", bbox_to_anchor=(1.01, 0.5), borderaxespad=0)
    apply_axis_style(axes[-1, 0], settings, xlabel=xlabel)
    if title:
        fig.suptitle(title)
    fig.tight_layout(rect=[0, 0, 0.86, 1])
    equilibration = results[0].source.study[labels[0]].equilibration
    _footnote(
        fig,
        counts,
        equilibration,
        f"production window t >= {equilibration}",
        drawn="Bands",
        of=f"each line's condition mean at each {xlabel.lower()}",
        counts="n in each panel title",
        points="Lines: condition means over replicates",
    )
    return _save(fig, output_dir, name, settings)


def plot_differences(
    values: ReplicateValues,
    report: Any,
    output_dir: str | Path,
    name: str,
    title: str | None = None,
    plot_settings: PlotSettings | None = None,
    xlabel: str = "label",
) -> Path:
    """Draw each condition's per-label difference from the control, one panel per condition.

    ``report`` is the ``compare()`` report of the labelled ``values``. Each
    panel draws, along the labels, ``delta`` of every comparison row, the
    mean of the condition minus the mean of the control, with a band of its
    95 percent interval from the report's test, Welch's by default, and a
    point on each label that is significant after the Benjamini-Hochberg
    correction. The panels share the y axis, and a label with no interval
    leaves a gap in the band.

    Returns
    -------
    Path
        The saved figure file.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    from polyzymd.analyses.shared.plotting import apply_axis_style, uncertainty_footnote_text

    settings, _, colors = _setup(values.source, plot_settings)
    try:
        where = {str(k): float(k) for k in values.labels}
    except (TypeError, ValueError):
        where = {str(k): float(i) for i, k in enumerate(values.labels)}
    conditions = list(dict.fromkeys(row.b for row in report.pairwise))
    fig, axes = plt.subplots(
        len(conditions),
        1,
        figsize=(14, 2.6 * len(conditions) + 1),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    for ax, label in zip(axes[:, 0], conditions, strict=True):
        rows = [row for row in report.pairwise if row.b == label]
        x = np.array([where[row.entry] for row in rows])
        delta = np.array([row.delta for row in rows])
        low, high = (
            np.array([row.delta_ci95[i] if row.delta_ci95 else np.nan for row in rows])
            for i in (0, 1)
        )
        color = colors[label]
        ax.axhline(0.0, color="0.5", linewidth=0.8)
        ax.plot(x, delta, color=color, linewidth=1.5, label="difference")
        ax.fill_between(x, low, high, color=color, alpha=0.25, label="95% interval")
        marked = [row.significant for row in rows]
        ax.scatter(x[marked], delta[marked], color="red", s=14, zorder=3, label="significant")
        n = f"n = {len(values.values[label])} vs {len(values.values[rows[0].a])}"
        apply_axis_style(
            ax,
            settings,
            title=f"{label} minus {rows[0].a} ({n})",
            ylabel=_label(f"delta {values.metric}", values.unit),
        )
    axes[0, 0].legend(loc="upper right", fontsize=8)
    apply_axis_style(axes[-1, 0], settings, xlabel=xlabel)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    family = next((row.family_size for row in report.pairwise if row.family_size), None)
    _note(
        fig,
        uncertainty_footnote_text(
            equilibration=report.equilibration,
            drawn="Bands",
            of=f"the difference of condition means at each {xlabel.lower()}",
            method=_TESTS.get(report.pairwise[0].test, report.pairwise[0].test),
            counts="n of both conditions in each panel title",
            points="Lines: condition mean minus control mean; red points: significant after "
            f"the Benjamini-Hochberg correction over a family of {family} tests; a label with "
            "fewer than two replicates in either condition, or no variance, has no interval",
        ),
    )
    return _save(fig, output_dir, name, settings)


def plot_values(
    results: Sequence[ReplicateValues],
    labels: Sequence[str] | None = None,
    output_dir: str | Path | None = None,
    name: str = "values",
    title: str | None = None,
    plot_settings: PlotSettings | None = None,
) -> Path:
    """Draw several results in one figure, one group per result and one bar per condition.

    Each bar is the condition's mean with its 95 percent Student t interval
    from :func:`~polyzymd.analyses.shared.statistics.mean_sem_ci`, drawn by
    :func:`~polyzymd.analyses.shared.plotting.grouped_bars` with every
    replicate value as a point. ``labels`` name the
    groups and default to the source names. ``output_dir`` defaults to the
    ``figures`` folder next to ``polyzymd_results``.

    Returns
    -------
    Path
        The saved figure file.

    Raises
    ------
    ProtocolError
        If the results do not share one unit and one set of conditions.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.shared.plotting import apply_axis_style, apply_legend, grouped_bars
    from polyzymd.analyses.shared.statistics import mean_sem_ci

    results = list(results)
    unit = _one_unit(results, "Results")
    labels = list(labels) if labels is not None else [item.source.name for item in results]
    settings, conditions, colors = _setup(results[0].source, plot_settings)
    if len(labels) != len(results) or any(set(item.values) != set(conditions) for item in results):
        raise ProtocolError(
            "Results in one figure need one label each and the same conditions.",
            hint="Plot results of the same study, with one label per result.",
        )
    data = [[item.values[label] for item in results] for label in conditions]
    stats = [[mean_sem_ci(entry) for entry in row] for row in data]
    series = [
        (
            f"{label} (n = {len(row[0])})",
            [item.mean for item in cells],
            [item.sem for item in cells],
        )
        for label, row, cells in zip(conditions, data, stats, strict=True)
    ]
    fig, ax = plt.subplots(figsize=(10, 6))
    positions = np.arange(len(results))
    grouped_bars(
        ax,
        positions,
        series,
        [colors[label] for label in conditions],
        settings,
        reference_line=None,
        replicate_values=data,
    )
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    metrics = {item.metric for item in results}
    fraction = all(item.is_fraction for item in results)
    quantity = metrics.pop() if len(metrics) == 1 else "fraction of frames" if fraction else "value"
    if fraction:
        ax.set_ylim(0, 1.05)
    apply_axis_style(ax, settings, title=title, ylabel=_label(quantity, unit))
    apply_legend(ax, settings, loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0)
    fig.tight_layout(rect=[0, 0, 0.78, 1])
    counts = {len(entry) for row in data for entry in row}
    equilibration = results[0].source.study[conditions[0]].equilibration
    _footnote(
        fig,
        counts,
        equilibration,
        f"production window t >= {equilibration}",
        of="the condition mean",
        counts="n per condition in the legend",
        points="Bars: condition means; points: per-replicate values",
    )
    folder = output_dir or results[0].source.path.parent.parent / "figures"
    return _save(fig, folder, name, settings)
