"""Shared plotting utilities for analysis plugins.

This module provides reusable plotting helper functions extracted from the
plotter infrastructure.  Analysis plugins import these functions to apply
consistent styling, save figures with watermarks, and render common chart
elements (grouped bars, heatmap annotations, etc.) without inheriting from
a base class.

All functions accept a ``PlotSettings`` object (from
``polyzymd.compare.config``) so that user-configured themes, palettes,
and DPI settings are respected automatically.

Examples
--------
>>> from polyzymd.analyses.shared.plotting import (
...     apply_axis_style, apply_legend, get_colors, save_figure,
... )
>>>
>>> fig, ax = plt.subplots()
>>> colors = get_colors(3, plot_settings)
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

    from polyzymd.compare.config import PlotSettings, PlotTheme

logger = logging.getLogger(__name__)

_UNSET = object()  # sentinel for apply_legend defaults


# ---------------------------------------------------------------------------
# Theme access
# ---------------------------------------------------------------------------


def get_theme(plot_settings: "PlotSettings") -> "PlotTheme":
    """Return the resolved ``PlotTheme`` from *plot_settings*.

    Parameters
    ----------
    plot_settings : PlotSettings
        Global plot settings (carries a ``.theme`` property).

    Returns
    -------
    PlotTheme
    """
    return plot_settings.theme


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


def get_colors(n: int, plot_settings: "PlotSettings") -> list:
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
    try:
        import seaborn as sns

        return list(sns.color_palette(plot_settings.color_palette, n))
    except ImportError:
        pass

    import matplotlib.pyplot as plt

    palette = plot_settings.color_palette
    try:
        cmap = plt.colormaps[palette]
    except (KeyError, ValueError):
        logger.warning(
            "Color palette %r is not a valid matplotlib colormap "
            "(seaborn is not installed). Falling back to 'tab10'.",
            palette,
        )
        cmap = plt.colormaps["tab10"]
    return [cmap(i / max(1, n - 1)) for i in range(n)]


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

    Returns
    -------
    Path
        Path to saved figure.
    """
    import matplotlib.pyplot as plt

    from polyzymd.core.experimental import annotate_experimental_figure

    output_path.parent.mkdir(parents=True, exist_ok=True)

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

    try:
        fig.savefig(
            output_path,
            dpi=plot_settings.dpi,
            bbox_inches="tight",
            facecolor="white",
            edgecolor="none",
        )
    finally:
        plt.close(fig)
    logger.info(f"Saved plot: {output_path}")
    return output_path


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------


def find_json(
    directory: Path,
    preferred: str,
    glob_pattern: str = "*.json",
) -> Path | None:
    """Locate a JSON result file inside *directory*.

    Tries the *preferred* filename first; if it does not exist, falls
    back to the first file matching *glob_pattern* (sorted
    lexicographically so results are deterministic).

    Parameters
    ----------
    directory : Path
        Directory to search.
    preferred : str
        Exact filename to try first (e.g. ``"rmsf_aggregated.json"``).
    glob_pattern : str, optional
        Glob to use as fallback, by default ``"*.json"``.

    Returns
    -------
    Path | None
        Path to the located file, or ``None`` if nothing was found.
    """
    exact = directory / preferred
    if exact.is_file():
        return exact
    candidates = sorted(directory.glob(glob_pattern))
    return candidates[0] if candidates else None


# ---------------------------------------------------------------------------
# Grouped bar chart
# ---------------------------------------------------------------------------


def grouped_bars(
    ax: "Axes",
    x: "np.ndarray",
    series: "Sequence[tuple[str, Sequence[float], Sequence[float]]]",
    colors: "Sequence",
    plot_settings: "PlotSettings",
    *,
    bar_width: float | None = None,
    show_error: bool = True,
    reference_line: float | None = 0.0,
    reference_label: str = "Neutral (0)",
    replicate_values: "Sequence[Sequence[Sequence[float]]] | None" = None,
    **style_overrides,
) -> None:
    """Render grouped bars with optional error bars and reference line.

    Style values (alpha, capsize, edgecolor, linewidth, dot_size, etc.)
    are read from ``plot_settings.theme``.  Callers can override any of
    them via ``**style_overrides`` using the theme field names as keys.

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
            bar_kwargs["yerr"] = errors
        ax.bar(np.asarray(x) + offset, means, **bar_kwargs)

        # Overlay jittered replicate dots
        if replicate_values is not None and i < len(replicate_values):
            rng = np.random.default_rng(seed=42 + i)
            cond_reps = replicate_values[i]
            for j in range(len(x)):
                if j < len(cond_reps) and cond_reps[j] is not None and len(cond_reps[j]) > 0:
                    rep_vals = np.asarray(cond_reps[j], dtype=float)
                    jitter = rng.uniform(-w * 0.3, w * 0.3, size=len(rep_vals))
                    ax.scatter(
                        np.full_like(rep_vals, float(x[j]) + offset) + jitter,
                        rep_vals,
                        color=dot_c,
                        s=dot_s,
                        zorder=5,
                        alpha=dot_a,
                        edgecolors="none",
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
# Heatmap cell annotation
# ---------------------------------------------------------------------------


def annotate_cells(
    ax: "Axes",
    matrix: "np.ndarray",
    plot_settings: "PlotSettings",
    *,
    fmt: str = ".2f",
    fontsize: int | None = None,
    threshold: float = 0.3,
    sem_matrix: "np.ndarray | None" = None,
    show_sign: bool = True,
    linespacing: float | None = None,
) -> None:
    """Annotate heatmap cells with formatted values.

    Iterates over every element of *matrix* and places a text label at
    the corresponding (col, row) position on *ax*.  NaN cells are
    skipped.  Text colour flips between black and white depending on
    the background intensity (controlled by *threshold*).

    Parameters
    ----------
    ax : matplotlib Axes
        The axes containing the heatmap image.
    matrix : np.ndarray
        2-D array of values (rows x cols) matching the heatmap.
    plot_settings : PlotSettings
        Global plot settings.
    fmt : str, optional
        Format spec for the value, by default ``".2f"``.
    fontsize : int | None, optional
        Annotation font size.  When ``None`` (default), uses
        ``plot_settings.theme.annotation_fontsize``.
    threshold : float, optional
        Absolute-value threshold above which text turns white.
    sem_matrix : np.ndarray | None, optional
        If provided, a second line ``±{sem}`` is appended when the
        SEM value is finite.
    show_sign : bool, optional
        Prefix positive values with ``"+"`` , by default ``True``.
    linespacing : float | None, optional
        Passed to ``ax.text(linespacing=...)`` when SEM is shown.
    """
    import numpy as np

    fs = fontsize if fontsize is not None else plot_settings.theme.annotation_fontsize

    n_rows, n_cols = matrix.shape
    for i in range(n_rows):
        for j in range(n_cols):
            val = matrix[i, j]
            if not np.isfinite(val):
                continue
            text_color = "white" if abs(val) > threshold else "black"
            sign = "+" if show_sign and val > 0 else ""
            label_str = f"{sign}{val:{fmt}}"
            if sem_matrix is not None:
                sem = sem_matrix[i, j]
                if not np.isnan(sem):
                    label_str = f"{label_str}\n\u00b1{sem:{fmt}}"
            kwargs: dict = {
                "ha": "center",
                "va": "center",
                "color": text_color,
                "fontsize": fs,
            }
            if linespacing is not None:
                kwargs["linespacing"] = linespacing
            ax.text(j, i, label_str, **kwargs)


# ---------------------------------------------------------------------------
# Symmetric colour-limit helper
# ---------------------------------------------------------------------------


def symmetric_clim(
    values: "Sequence[float] | np.ndarray",
    pad: float = 0.1,
) -> tuple[float, float]:
    """Compute symmetric colour limits centred on zero.

    Parameters
    ----------
    values : sequence of float or ndarray
        Finite data values to derive limits from.
    pad : float, optional
        Extra padding added to both sides, by default ``0.1``.

    Returns
    -------
    tuple[float, float]
        ``(vmin, vmax)`` with ``vmin == -vmax`` (before padding).
    """
    import numpy as np

    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return (-pad, pad)
    max_abs = float(max(abs(arr.min()), abs(arr.max())))
    return (-(max_abs + pad), max_abs + pad)
