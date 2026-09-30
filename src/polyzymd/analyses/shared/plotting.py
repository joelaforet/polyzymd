"""Plotting helpers that the figures of :mod:`polyzymd.analyses.figures` are drawn with.

They style axes and legends, pick and order condition colours, draw grouped
bars with every replicate value, work out error bar and band half-widths,
add the uncertainty footnote, and save figures with the PolyzyMD watermark.

All functions accept a ``PlotSettings`` object (from
``polyzymd.config.analysis_settings``) so that user-configured themes, palettes,
and DPI settings are respected automatically.

Examples
--------
>>> from polyzymd.analyses.shared.plotting import (
...     apply_axis_style, apply_legend, get_palette_colors, save_figure,
... )
>>>
>>> fig, ax = plt.subplots()
>>> colors = get_palette_colors(3, plot_settings)
>>> ax.bar(x, y, color=colors[0])
>>> apply_axis_style(ax, plot_settings, title="My Plot", ylabel="Value (Å)")
>>> apply_legend(ax, plot_settings)
>>> save_figure(fig, output_dir / "my_plot.png", plot_settings)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

from polyzymd.core.branding import PLOT_WATERMARK

if TYPE_CHECKING:
    import numpy as np
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    from polyzymd.config.analysis_settings import PlotSettings

from polyzymd.analyses.exceptions import StatisticsError

logger = logging.getLogger(__name__)

_UNSET = object()  # sentinel for apply_legend defaults


# ---------------------------------------------------------------------------
# Axis styling
# ---------------------------------------------------------------------------


def apply_axis_style(
    ax: "Axes",
    plot_settings: "PlotSettings",
    *,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
) -> None:
    """Apply standard axis chrome from the theme.

    Hides spines according to theme settings, sizes tick labels, and
    optionally sets title / xlabel / ylabel with themed font sizes.

    Parameters
    ----------
    ax : matplotlib Axes
        Target axes to style.
    plot_settings : PlotSettings
        Global plot settings.
    title : str, optional
        Axes title.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    """
    t = plot_settings.theme
    if t.hide_top_spine:
        ax.spines["top"].set_visible(False)
    if t.hide_right_spine:
        ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=t.tick_fontsize)
    if title is not None:
        ax.set_title(title, fontsize=t.title_fontsize, fontweight=t.title_fontweight)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=t.label_fontsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=t.label_fontsize)


def apply_legend(
    ax: "Axes",
    plot_settings: "PlotSettings",
    *,
    loc: str | None = None,
    bbox_to_anchor: tuple[float, float] | None | object = _UNSET,
    fontsize: int | None = None,
    **kwargs: Any,
) -> None:
    """Apply legend with themed defaults.

    Uses ``theme.legend_loc`` and ``theme.legend_bbox`` unless overridden
    by the caller.  Extra *kwargs* are forwarded to ``ax.legend()``.

    Parameters
    ----------
    ax : matplotlib Axes
        Target axes.
    plot_settings : PlotSettings
        Global plot settings.
    loc : str, optional
        Override ``theme.legend_loc``.
    bbox_to_anchor : tuple of float or None, optional
        Override ``theme.legend_bbox``.  Pass ``None`` explicitly to
        suppress the bbox (e.g. for inside-axes placement).
    fontsize : int, optional
        Override ``theme.legend_fontsize``.
    **kwargs
        Forwarded to ``ax.legend()``.
    """
    t = plot_settings.theme
    resolved_loc = loc or t.legend_loc
    resolved_fs = fontsize or t.legend_fontsize

    legend_kwargs: dict[str, Any] = {"loc": resolved_loc, "fontsize": resolved_fs}

    # Resolve bbox_to_anchor: _UNSET → theme default, None → omit
    if bbox_to_anchor is _UNSET:
        legend_kwargs["bbox_to_anchor"] = t.legend_bbox
    elif bbox_to_anchor is not None:
        legend_kwargs["bbox_to_anchor"] = bbox_to_anchor

    legend_kwargs.update(kwargs)
    ax.legend(**legend_kwargs)


# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------


def get_palette_colors(n: int, plot_settings: "PlotSettings") -> list:
    """Get *n* distinct colors from the configured palette.

    Tries seaborn first (richer palette support), falls back to a
    matplotlib colormap sampled at evenly-spaced intervals.

    Parameters
    ----------
    n : int
        Number of colors needed.
    plot_settings : PlotSettings
        Global plot settings (carries ``color_palette``).

    Returns
    -------
    list
        List of color values (RGB tuples or matplotlib color specs).
    """
    palette = plot_settings.color_palette
    try:
        import seaborn as sns

        try:
            return list(sns.color_palette(palette, n))
        except ValueError:
            logger.warning(
                "Color palette %r is not a valid seaborn palette. Falling back to 'tab10'.",
                palette,
            )
            return _matplotlib_colormap_colors("tab10", n)
    except ImportError:
        pass

    try:
        return _matplotlib_colormap_colors(palette, n)
    except (KeyError, ValueError):
        logger.warning(
            "Color palette %r is not a valid matplotlib colormap "
            "(seaborn is not installed). Falling back to 'tab10'.",
            palette,
        )
        return _matplotlib_colormap_colors("tab10", n)


def _matplotlib_colormap_colors(colormap_name: str, n: int) -> list:
    """Sample colors from a matplotlib colormap.

    Parameters
    ----------
    colormap_name : str
        Name of the matplotlib colormap to sample.
    n : int
        Number of colors to return.

    Returns
    -------
    list
        RGBA colors sampled across the colormap range.
    """
    import matplotlib.pyplot as plt

    cmap = plt.colormaps[colormap_name]
    return [cmap(i / max(1, n - 1)) for i in range(n)]


def order_condition_labels(labels: Sequence[str], plot_settings: "PlotSettings") -> list[str]:
    """Return condition labels in semantic plot order when enabled.

    Ordering only affects plot display order. It does not alter comparison
    statistics, rankings, or condition result files.

    Parameters
    ----------
    labels : sequence of str
        Condition labels in their original order.
    plot_settings : PlotSettings
        Global plot settings carrying optional semantic color settings.

    Returns
    -------
    list of str
        Ordered labels for plotting.
    """
    label_list = list(labels)
    semantic = getattr(plot_settings, "semantic_colors", None)
    if semantic is None or not semantic.enabled:
        return label_list

    remaining = list(label_list)
    ordered: list[str] = []
    for label in semantic.order:
        if label in remaining:
            ordered.append(label)
            remaining.remove(label)

    indexed_remaining = list(enumerate(remaining))
    with_order: list[tuple[int, str, int]] = []
    without_order: list[tuple[int, str]] = []
    for relative_index, label in indexed_remaining:
        condition = semantic.conditions.get(label)
        if condition is not None and condition.order is not None:
            with_order.append((condition.order, label, relative_index))
        else:
            without_order.append((relative_index, label))

    ordered.extend(label for _, label, _ in sorted(with_order, key=lambda item: (item[0], item[2])))
    ordered.extend(label for _, label in without_order)
    return ordered


def get_condition_color_map(
    labels: Sequence[str],
    plot_settings: "PlotSettings",
    *,
    control_label: str | None = None,
) -> dict[str, Any]:
    """Return a label-to-color map using semantic condition color rules.

    Resolution precedence is manual color, condition color, control color,
    family/value color, missing metadata fallback, then existing palette fallback.
    Invalid color or colormap values warn and continue to a safe fallback.

    Parameters
    ----------
    labels : sequence of str
        Condition labels in their original palette-alignment order.
    plot_settings : PlotSettings
        Global plot settings carrying optional semantic color settings.
    control_label : str, optional
        Label that should use the configured semantic control color.

    Returns
    -------
    dict of str to Any
        Mapping from each label to its resolved matplotlib-compatible color.
    """
    label_list = list(labels)
    palette_colors = get_palette_colors(len(label_list), plot_settings)
    palette_by_label = dict(zip(label_list, palette_colors))
    semantic = getattr(plot_settings, "semantic_colors", None)
    if semantic is None or not semantic.enabled:
        return dict(palette_by_label)

    observed_values = _collect_family_values(label_list, semantic.conditions)
    color_map: dict[str, Any] = {}
    for label in label_list:
        color_map[label] = _resolve_condition_color(
            label,
            semantic,
            palette_by_label,
            observed_values,
            control_label=control_label,
        )
    return color_map


def _resolve_condition_color(
    label: str,
    semantic: Any,
    palette_by_label: dict[str, Any],
    observed_values: dict[str, list[Any]],
    *,
    control_label: str | None,
) -> Any:
    """Resolve one condition color using semantic precedence."""
    manual_color = _validated_color(
        semantic.manual_colors.get(label), f"manual color for {label!r}"
    )
    if manual_color is not None:
        return manual_color

    condition = semantic.conditions.get(label)
    if condition is None:
        default_color = _validated_color(semantic.default_color, "semantic default_color")
        return default_color if default_color is not None else palette_by_label[label]

    condition_color = _validated_color(condition.color, f"condition color for {label!r}")
    if condition_color is not None:
        return condition_color

    is_control = control_label is not None and label == control_label
    is_control = is_control or condition.role == "control"
    if is_control:
        control_color = _validated_color(semantic.control_color, "semantic control_color")
        if control_color is not None:
            return control_color

    family_color = _resolve_family_color(condition, semantic.families, observed_values, label)
    if family_color is not None:
        return family_color

    missing_color = _validated_color(semantic.missing_color, "semantic missing_color")
    return missing_color if missing_color is not None else palette_by_label[label]


def _collect_family_values(
    labels: Sequence[str], conditions: dict[str, Any]
) -> dict[str, list[Any]]:
    """Collect observed semantic values by family in label order."""
    observed_values: dict[str, list[Any]] = {}
    for label in labels:
        condition = conditions.get(label)
        if condition is None or condition.family is None or condition.value is None:
            continue
        values = observed_values.setdefault(condition.family, [])
        if condition.value not in values:
            values.append(condition.value)
    return observed_values


def _validated_color(color: Any, context: str) -> Any | None:
    """Return a color when matplotlib accepts it, otherwise warn."""
    if color is None:
        return None

    from matplotlib.colors import is_color_like

    if is_color_like(color):
        return color
    logger.warning("Invalid %s %r. Falling back to the next available color rule.", context, color)
    return None


def _resolve_family_color(
    condition: Any,
    families: dict[str, Any],
    observed_values: dict[str, list[Any]],
    label: str,
) -> Any | None:
    """Resolve a condition color from its semantic family and value."""
    if condition.family is None or condition.value is None:
        return None

    family = families.get(condition.family)
    if family is None:
        logger.warning(
            "Condition %r references unknown semantic color family %r.",
            label,
            condition.family,
        )
        return None

    value_color = _resolve_explicit_value_color(family.value_colors, condition.value, label)
    if value_color is not None:
        return value_color

    cmap = _get_colormap(family.colormap, label)
    if cmap is None:
        return None

    if family.scale == "ordinal":
        fraction = _ordinal_fraction(
            condition.value, family, observed_values.get(condition.family, [])
        )
    else:
        fraction = _linear_fraction(
            condition.value, family, observed_values.get(condition.family, [])
        )

    if fraction is None:
        return None
    if family.reverse:
        fraction = 1.0 - fraction
    low, high = family.colormap_range
    return cmap(low + ((high - low) * fraction))


def _resolve_explicit_value_color(
    value_colors: dict[str, Any], value: Any, label: str
) -> Any | None:
    """Resolve an exact semantic value color if configured."""
    if value in value_colors:
        color = value_colors[value]
    else:
        color = value_colors.get(str(value))
    return _validated_color(color, f"value color for {label!r}")


def _get_colormap(colormap_name: str, label: str) -> Any | None:
    """Return a matplotlib colormap or warn and return ``None``."""
    import matplotlib.pyplot as plt

    try:
        return plt.colormaps[colormap_name]
    except (KeyError, ValueError):
        logger.warning(
            "Semantic color colormap %r for condition %r is invalid. Falling back.",
            colormap_name,
            label,
        )
        return None


def _ordinal_fraction(value: Any, family: Any, observed_values: list[Any]) -> float | None:
    """Return an ordinal colormap fraction for a semantic value."""
    value_order = family.value_order or observed_values
    if value not in value_order:
        logger.warning("Semantic ordinal value %r is not present in value_order.", value)
        return None
    if len(value_order) <= 1:
        return 0.5
    return value_order.index(value) / (len(value_order) - 1)


def _linear_fraction(value: Any, family: Any, observed_values: list[Any]) -> float | None:
    """Return a linear colormap fraction for a semantic numeric value."""
    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        logger.warning("Semantic linear value %r is not numeric.", value)
        return None

    numeric_values: list[float] = []
    for observed_value in observed_values:
        try:
            numeric_values.append(float(observed_value))
        except (TypeError, ValueError):
            continue

    vmin = family.vmin if family.vmin is not None else min(numeric_values, default=numeric_value)
    vmax = family.vmax if family.vmax is not None else max(numeric_values, default=numeric_value)
    if vmax == vmin:
        return 0.5
    return min(1.0, max(0.0, (numeric_value - vmin) / (vmax - vmin)))


# ---------------------------------------------------------------------------
# Figure saving
# ---------------------------------------------------------------------------


def get_output_path(output_dir: Path, name: str, plot_settings: "PlotSettings") -> Path:
    """Generate output file path with correct format extension.

    Parameters
    ----------
    output_dir : Path
        Output directory.
    name : str
        Base filename (without extension).
    plot_settings : PlotSettings
        Global plot settings (carries ``format``).

    Returns
    -------
    Path
        Full output path with extension.
    """
    return output_dir / f"{name}.{plot_settings.format}"


def save_figure(
    fig: "Figure",
    output_path: Path,
    plot_settings: "PlotSettings",
    *,
    experimental_features: Sequence[str] | None = None,
    close: bool = True,
) -> Path:
    """Save figure with DPI, watermark, and optional experimental stamp.

    Parameters
    ----------
    fig : matplotlib Figure
        Figure to save.
    output_path : Path
        Output file path.
    plot_settings : PlotSettings
        Global plot settings (carries ``dpi`` and ``theme``).
    experimental_features : sequence of str or None, optional
        Experimental feature ids to stamp onto the figure.
    close : bool, optional
        If True, close the figure after saving. Set False when the caller
        needs to keep using the figure object.

    Returns
    -------
    Path
        Path to saved figure.
    """
    import matplotlib.pyplot as plt

    from polyzymd.core.experimental import annotate_experimental_figure

    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        if experimental_features:
            annotate_experimental_figure(fig, experimental_features)

        # Add watermark if enabled
        if plot_settings.theme.show_watermark:
            fig.text(
                0.99,
                0.01,
                PLOT_WATERMARK,
                fontsize=7,
                color="dimgray",
                ha="right",
                va="bottom",
                alpha=0.85,
                style="italic",
            )

        fig.savefig(
            output_path,
            dpi=plot_settings.dpi,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )
    finally:
        if close:
            plt.close(fig)
    logger.info(f"Saved plot: {output_path}")
    return output_path


# ---------------------------------------------------------------------------
# Grouped bar chart
# ---------------------------------------------------------------------------


def finite_numeric_values(values: Any) -> "np.ndarray":
    """Return finite numeric values as a one-dimensional float array.

    Parameters
    ----------
    values : Any
        Candidate scalar or sequence of replicate values. ``None`` and
        non-numeric inputs are treated as missing data.

    Returns
    -------
    numpy.ndarray
        One-dimensional array containing only finite floats. The array is
        empty when no finite numeric values are available.
    """
    import numpy as np

    if values is None:
        return np.array([], dtype=float)

    try:
        value_array = np.asarray(values, dtype=float)
    except (TypeError, ValueError):
        if isinstance(values, str | bytes):
            return np.array([], dtype=float)
        try:
            iterator = iter(values)
        except TypeError:
            try:
                scalar_value = float(values)
            except (TypeError, ValueError):
                return np.array([], dtype=float)
            if np.isfinite(scalar_value):
                return np.array([scalar_value], dtype=float)
            return np.array([], dtype=float)

        finite_values: list[float] = []
        for value in iterator:
            try:
                numeric_value = float(value)
            except (TypeError, ValueError):
                continue
            if np.isfinite(numeric_value):
                finite_values.append(numeric_value)
        return np.array(finite_values, dtype=float)

    if value_array.ndim == 0:
        value_array = value_array.reshape(1)

    value_array = value_array.ravel()
    return value_array[np.isfinite(value_array)]


def replicate_jitter_offsets(n_values: int, bar_width: float) -> "np.ndarray":
    """Return deterministic offsets for replicate dot overlays.

    Offsets are centred on the corresponding bar position so overlays are
    reproducible across runs and independent of random-number state.

    Parameters
    ----------
    n_values : int
        Number of replicate values to display for one bar.
    bar_width : float
        Width or height of the corresponding bar, depending on orientation.

    Returns
    -------
    numpy.ndarray
        Jitter offsets centred around zero.
    """
    import numpy as np

    if n_values <= 0:
        return np.array([], dtype=float)
    if n_values == 1:
        return np.array([0.0], dtype=float)

    max_jitter = bar_width * 0.25
    return np.linspace(-max_jitter, max_jitter, n_values)


def scatter_replicate_values(
    ax: "Axes",
    bar_positions: "Sequence[float] | np.ndarray",
    replicate_values: "Sequence[Any]",
    plot_settings: "PlotSettings",
    *,
    orientation: str = "vertical",
    bar_width: float = 0.8,
    dot_color: Any | None = None,
    dot_size: float | None = None,
    dot_alpha: float | None = None,
    zorder: float = 5,
) -> int:
    """Overlay jittered per-replicate values on bars.

    For vertical bars, bar positions are x-coordinates, jitter is applied in
    x, and replicate values are plotted on y. For horizontal bars, bar
    positions are y-coordinates, jitter is applied in y, and replicate values
    are plotted on x.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes containing the bar chart.
    bar_positions : sequence of float or numpy.ndarray
        Bar centre positions aligned to ``replicate_values``.
    replicate_values : sequence of Any
        Per-bar replicate values. Each item may be a scalar or sequence; only
        finite numeric values are plotted.
    plot_settings : PlotSettings
        Global plot settings whose theme provides default dot styling.
    orientation : {"vertical", "horizontal"}, optional
        Bar orientation, by default ``"vertical"``.
    bar_width : float, optional
        Width or height of the bars, by default ``0.8``.
    dot_color : Any, optional
        Override for theme dot colour.
    dot_size : float, optional
        Override for theme dot size. Dots are skipped when non-positive.
    dot_alpha : float, optional
        Override for theme dot alpha. Dots are skipped when non-positive.
    zorder : float, optional
        Matplotlib z-order for dot overlays, by default ``5``.

    Returns
    -------
    int
        Number of scatter calls emitted.

    Raises
    ------
    ValueError
        If ``orientation`` is not ``"vertical"`` or ``"horizontal"``, or if
        ``bar_positions`` and ``replicate_values`` are not the same length.
    """
    import numpy as np

    if orientation not in {"vertical", "horizontal"}:
        raise ValueError("orientation must be 'vertical' or 'horizontal'")

    theme = plot_settings.theme
    resolved_size = theme.dot_size if dot_size is None else dot_size
    resolved_alpha = theme.dot_alpha if dot_alpha is None else dot_alpha
    if resolved_size <= 0 or resolved_alpha <= 0:
        return 0

    resolved_color = theme.dot_color if dot_color is None else dot_color
    positions = np.asarray(bar_positions, dtype=float)
    if len(replicate_values) != len(positions):
        raise ValueError(
            "replicate_values length must match bar_positions length "
            f"({len(replicate_values)} != {len(positions)})"
        )

    n_scattered = 0
    for idx, values in enumerate(replicate_values):
        rep_arr = finite_numeric_values(values)
        if rep_arr.size == 0:
            continue

        jitter = replicate_jitter_offsets(rep_arr.size, bar_width)
        position_arr = np.full(rep_arr.shape, float(positions[idx]), dtype=float) + jitter
        if orientation == "vertical":
            x_values = position_arr
            y_values = rep_arr
        else:
            x_values = rep_arr
            y_values = position_arr

        ax.scatter(
            x_values,
            y_values,
            color=resolved_color,
            s=resolved_size,
            zorder=zorder,
            alpha=resolved_alpha,
            edgecolors="none",
        )
        n_scattered += 1

    return n_scattered


def grouped_bars(
    ax: "Axes",
    x: "np.ndarray",
    series: "Sequence[tuple[str, Sequence[float], Sequence[float]]]",
    colors: "Sequence",
    plot_settings: "PlotSettings",
    *,
    bar_width: float | None = None,
    show_error: bool = True,
    error_bar: str = "ci95",
    n_replicates: int | None = None,
    reference_line: float | None = 0.0,
    reference_label: str = "Neutral (0)",
    replicate_values: "Sequence[Sequence[Sequence[float]]] | None" = None,
    **style_overrides,
) -> None:
    """Render grouped bars with optional error bars and reference line.

    Style values (alpha, capsize, edgecolor, linewidth, dot_size, etc.)
    are read from ``plot_settings.theme``.  Callers can override any of
    them via ``**style_overrides`` using the theme field names as keys.

    ``errors`` holds the standard error of each bar. It is converted to the
    interval named by ``error_bar`` using the per-bar replicate counts from
    ``replicate_values``, or the shared count in ``n_replicates`` when the
    per-bar values are not available. Without either, no error bar is drawn, so
    a figure never shows a raw SEM under a footnote claiming an interval.

    Parameters
    ----------
    ax : matplotlib Axes
        Target axes.
    x : np.ndarray
        1-D array of group centre positions (e.g. ``np.arange(n_groups)``).
    series : sequence of (label, means, errors)
        One tuple per condition.  *means* and *errors* must have the same
        length as *x*.
    colors : sequence
        One colour per condition (same length as *series*).
    plot_settings : PlotSettings
        Global plot settings.
    bar_width : float | None, optional
        Width of each individual bar.  When ``None`` (default) the width
        is computed as ``0.8 / len(series)``.
    show_error : bool, optional
        If ``False``, error bars are suppressed, by default ``True``.
    reference_line : float | None, optional
        Y-value for a horizontal reference line.  Set to ``None`` to
        skip, by default ``0.0``.
    reference_label : str, optional
        Legend label for the reference line, by default ``"Neutral (0)"``.
    replicate_values : sequence or None, optional
        Per-replicate values for jittered dot overlay.  Indexed as
        ``replicate_values[condition_idx][group_idx]`` -> sequence of
        floats (one per replicate).  When ``None`` (default), no dots
        are drawn.
    **style_overrides
        Override any theme field for this call only.  Accepted keys:
        ``bar_alpha``, ``bar_capsize``, ``bar_edgecolor``,
        ``bar_linewidth``, ``dot_size``, ``dot_alpha``, ``dot_color``,
        ``reference_line_color``, ``reference_line_style``,
        ``reference_line_width``.
    """
    import numpy as np

    t = plot_settings.theme
    alpha = style_overrides.get("bar_alpha", t.bar_alpha)
    capsize = style_overrides.get("bar_capsize", t.bar_capsize)
    edgecolor = style_overrides.get("bar_edgecolor", t.bar_edgecolor)
    linewidth = style_overrides.get("bar_linewidth", t.bar_linewidth)
    dot_s = style_overrides.get("dot_size", t.dot_size)
    dot_a = style_overrides.get("dot_alpha", t.dot_alpha)
    dot_c = style_overrides.get("dot_color", t.dot_color)
    ref_color = style_overrides.get("reference_line_color", t.reference_line_color)
    ref_style = style_overrides.get("reference_line_style", t.reference_line_style)
    ref_width = style_overrides.get("reference_line_width", t.reference_line_width)

    n = len(series)
    w = bar_width if bar_width is not None else 0.8 / max(n, 1)

    if replicate_values is not None:
        if len(replicate_values) != n:
            raise ValueError(
                f"replicate_values length must match series length ({len(replicate_values)} != {n})"
            )
        for idx, series_replicates in enumerate(replicate_values):
            if len(series_replicates) != len(x):
                raise ValueError(
                    "replicate_values entries must match x length "
                    f"for series {idx} ({len(series_replicates)} != {len(x)})"
                )

    for i, (label, means, errors) in enumerate(series):
        offset = (i - n / 2 + 0.5) * w
        bar_kwargs: dict = {
            "width": w,
            "label": label,
            "color": colors[i],
            "alpha": alpha,
            "capsize": capsize,
            "edgecolor": edgecolor,
            "linewidth": linewidth,
        }
        if show_error:
            if replicate_values is not None:
                errors_for_plot = error_bar_half_widths(
                    list(errors), replicate_values[i], error_bar=error_bar
                )
            else:
                errors_for_plot = error_bar_half_widths(
                    list(errors), error_bar=error_bar, n_replicates=n_replicates
                )
            if errors_for_plot is not None:
                bar_kwargs["yerr"] = errors_for_plot
        bar_positions = np.asarray(x) + offset
        ax.bar(bar_positions, means, **bar_kwargs)

        # Overlay jittered replicate dots
        if replicate_values is not None:
            scatter_replicate_values(
                ax,
                bar_positions,
                replicate_values[i],
                plot_settings,
                orientation="vertical",
                bar_width=w,
                dot_color=dot_c,
                dot_size=dot_s,
                dot_alpha=dot_a,
            )

    if reference_line is not None:
        ax.axhline(
            y=reference_line,
            color=ref_color,
            linestyle=ref_style,
            linewidth=ref_width,
            label=reference_label,
        )


# ---------------------------------------------------------------------------
# Uncertainty intervals on figures
# ---------------------------------------------------------------------------


def error_bar_half_widths(
    sems: Sequence[float | None],
    replicate_values: Sequence[Any] | None = None,
    *,
    error_bar: str = "ci95",
    n_replicates: int | None = None,
) -> list[float] | None:
    """Return the half width of the error bar to draw on each bar.

    With ``"ci95"`` each standard error is multiplied by the Student t factor
    for that bar's own replicate count, so the bar spans the 95 percent
    interval; with ``"sem"`` it is drawn unchanged. ``replicate_values`` gives
    the per-bar counts, or ``n_replicates`` a count shared by every bar. Bars
    backed by one replicate get no error bar, and the result is ``None`` when
    no bar has displayable uncertainty.
    """
    from polyzymd.analyses.shared.statistics import student_t_coverage_factor

    if error_bar not in ("ci95", "sem"):
        raise StatisticsError(f"error_bar must be 'ci95' or 'sem', got {error_bar!r}")

    counts: list[int] = []
    for index in range(len(sems)):
        if replicate_values is not None and index < len(replicate_values):
            counts.append(int(finite_numeric_values(replicate_values[index]).size))
        else:
            counts.append(int(n_replicates or 0))

    half_widths: list[float] = []
    any_uncertainty = False
    for sem, count in zip(sems, counts, strict=False):
        if sem is None or count < 2:
            half_widths.append(0.0)
            continue
        factor = 1.0
        if error_bar == "ci95":
            factor = student_t_coverage_factor(count) or 1.0
        half_widths.append(float(sem) * factor)
        any_uncertainty = True

    return half_widths if any_uncertainty else None


def band_half_widths(
    matrix: Any,
    *,
    error_bar: str = "ci95",
) -> "np.ndarray | None":
    """Return the half width of a shaded band over per-replicate traces.

    ``matrix`` is shaped ``(n_replicates, n_points)``. The band is centred on
    the mean across replicates and is ``None`` when fewer than two replicates
    make it inestimable.
    """
    import numpy as np

    from polyzymd.analyses.shared.statistics import student_t_coverage_factor

    if error_bar not in ("ci95", "sem"):
        raise StatisticsError(f"error_bar must be 'ci95' or 'sem', got {error_bar!r}")

    array = np.asarray(matrix, dtype=float)
    if array.ndim != 2 or array.shape[0] < 2:
        return None
    n = int(array.shape[0])
    sem = np.std(array, axis=0, ddof=1) / np.sqrt(float(n))
    if error_bar == "sem":
        return sem
    return sem * float(student_t_coverage_factor(n) or 1.0)


def add_uncertainty_footnote(
    fig: "Figure",
    *,
    error_bar: str = "ci95",
    n_replicates: int | None = None,
    equilibration: str | None = None,
    drawn: str = "Error bars",
    points: str = "Points are per-replicate values",
) -> str:
    """Write the sentence saying what a figure's error bars mean, and return it.

    Grossfield et al. (2018) ask that every figure describe the meaning and
    basis of its uncertainties. This is that sentence. ``drawn`` names the
    mark that shows the interval, such as ``"Band"``, and ``points`` says
    what the per-replicate marks are.
    """
    what = (
        f"{drawn}: 1 SEM (not a 95% interval)"
        if error_bar == "sem"
        else f"{drawn}: 95% CI (Student t)"
    )
    across = (
        f" across n = {n_replicates} replicates"
        if n_replicates and n_replicates >= 2
        else " across replicates"
    )
    window = f"; production window t >= {equilibration}" if equilibration else ""
    text = f"{what}{across}{window}. {points}."
    fig.text(
        0.01,
        0.01,
        text,
        fontsize=7,
        color="dimgray",
        ha="left",
        va="bottom",
    )
    return text
