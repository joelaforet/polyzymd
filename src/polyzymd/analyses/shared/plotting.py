"""Shared plotting utilities for analysis plugins.

This module provides reusable plotting helper functions extracted from the
plotter infrastructure.  Analysis plugins import these functions to apply
consistent styling, save figures with watermarks, and render common chart
elements (grouped bars, heatmap annotations, etc.) without inheriting from
a base class.

All functions accept a ``PlotSettings`` object (from
``polyzymd.config.comparison``) so that user-configured themes, palettes,
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
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

from polyzymd.core.branding import PLOT_WATERMARK

if TYPE_CHECKING:
    import numpy as np
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    from polyzymd.config.comparison import PlotSettings, PlotTheme

logger = logging.getLogger(__name__)

_UNSET = object()  # sentinel for apply_legend defaults


@dataclass(frozen=True)
class ArtifactPlotData:
    """Canonical artifacts loaded for plot-time data access."""

    analysis_dir: Path
    condition_artifact: Any | None
    replicate_artifacts: dict[int, Any]
    aggregated_dir: Path
    run_dirs: dict[int, Path]


def load_canonical_plot_artifacts(
    analysis_dir: Path,
    replicates: Sequence[int],
    *,
    require_condition: bool = False,
    require_replicates: bool = True,
) -> ArtifactPlotData:
    """Load plot inputs from canonical MDAnalysis artifacts only.

    The loader reads ``aggregated/result.json`` and the configured
    ``run_N/result.json`` files through :class:`ArtifactStore`. It never scans
    directories, opens non-canonical JSON files, or imports trajectory packages.

    Parameters
    ----------
    analysis_dir : Path
        Condition-level analysis directory containing ``aggregated`` and
        ``run_N`` subdirectories.
    replicates : sequence of int
        Configured replicate IDs to load. Extra run directories are ignored.
    require_condition : bool, optional
        Raise when ``aggregated/result.json`` is absent, by default False.
    require_replicates : bool, optional
        Raise when any configured ``run_N/result.json`` is absent, by default
        True.

    Returns
    -------
    ArtifactPlotData
        Loaded canonical condition and replicate artifacts.
    """

    from polyzymd.analyses.mda import ArtifactStore, ArtifactStoreError

    root = Path(analysis_dir)
    aggregated_dir = root / "aggregated"
    condition_artifact = None
    condition_path = aggregated_dir / "result.json"
    if condition_path.exists():
        condition_artifact = ArtifactStore(aggregated_dir).read_condition_result("result.json")
    elif require_condition:
        raise ArtifactStoreError(f"Missing canonical condition artifact: {condition_path}")

    replicate_artifacts: dict[int, Any] = {}
    run_dirs: dict[int, Path] = {}
    for replicate in replicates:
        replicate_id = int(replicate)
        run_dir = root / f"run_{replicate_id}"
        run_dirs[replicate_id] = run_dir
        replicate_path = run_dir / "result.json"
        if replicate_path.exists():
            replicate_artifacts[replicate_id] = ArtifactStore(run_dir).read_replicate_result(
                "result.json"
            )
        elif require_replicates:
            raise ArtifactStoreError(f"Missing canonical replicate artifact: {replicate_path}")

    return ArtifactPlotData(
        analysis_dir=root,
        condition_artifact=condition_artifact,
        replicate_artifacts=replicate_artifacts,
        aggregated_dir=aggregated_dir,
        run_dirs=run_dirs,
    )


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


def get_colors(n: int, plot_settings: "PlotSettings") -> list:
    """Get *n* existing palette colors.

    This backward-compatible alias preserves the historical behavior of
    ``get_colors()`` while semantic condition color helpers opt in separately.

    Parameters
    ----------
    n : int
        Number of colors needed.
    plot_settings : PlotSettings
        Global plot settings.

    Returns
    -------
    list
        List of color values from the configured palette.
    """
    return get_palette_colors(n, plot_settings)


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


def get_condition_colors(
    labels: Sequence[str],
    plot_settings: "PlotSettings",
    *,
    control_label: str | None = None,
) -> list:
    """Return colors for condition labels using semantic settings if enabled.

    Parameters
    ----------
    labels : sequence of str
        Condition labels in plot order.
    plot_settings : PlotSettings
        Global plot settings carrying optional semantic color settings.
    control_label : str, optional
        Label that should use the configured semantic control color.

    Returns
    -------
    list
        Color values aligned to ``labels``.
    """
    color_map = get_condition_color_map(labels, plot_settings, control_label=control_label)
    return [color_map[label] for label in labels]


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


def has_replicate_uncertainty(
    replicate_values: Any = None,
    *,
    n_replicates: int | None = None,
) -> bool:
    """Return whether replicate-level uncertainty can be displayed.

    Parameters
    ----------
    replicate_values : Any, optional
        Per-condition or per-bar replicate values. Finite numeric entries are
        counted after coercion.
    n_replicates : int or None, optional
        Explicit replicate count when the raw replicate values are not
        available.

    Returns
    -------
    bool
        True when at least two finite independent replicate values are present.
    """

    if n_replicates is not None:
        return n_replicates >= 2
    return finite_numeric_values(replicate_values).size >= 2


def suppress_singleton_errors(
    errors: Sequence[float],
    replicate_values: Sequence[Any] | None,
) -> list[float] | None:
    """Return errors with singleton replicate uncertainties suppressed.

    Parameters
    ----------
    errors : sequence of float
        SEM or uncertainty values aligned to ``replicate_values``.
    replicate_values : sequence or None
        Per-bar replicate values used to decide whether an error bar is
        statistically displayable.

    Returns
    -------
    list of float or None
        Sanitized error values. Returns ``None`` when no bar has replicate
        uncertainty, allowing callers to omit error bars entirely.
    """

    if replicate_values is None:
        return list(errors)

    sanitized: list[float] = []
    has_any_uncertainty = False
    for error, values in zip(errors, replicate_values, strict=False):
        if has_replicate_uncertainty(values):
            sanitized.append(float(error))
            has_any_uncertainty = True
        else:
            sanitized.append(0.0)

    return sanitized if has_any_uncertainty else None


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


def scatter_stacked_segment_replicates(
    ax: "Axes",
    x_position: float,
    bottom_value: float,
    replicate_values: Sequence[Any],
    plot_settings: "PlotSettings",
    *,
    replicate_base_values: Sequence[Any] | None = None,
    positive_base_values: Sequence[Any] | None = None,
    negative_base_values: Sequence[Any] | None = None,
    bar_width: float = 0.8,
    dot_color: Any | None = None,
    dot_size: float | None = None,
    dot_alpha: float | None = None,
    placement: str = "center",
    zorder: float = 5,
) -> int:
    """Overlay replicate dots on stacked segments.

    The per-component replicate value is a segment height, not an absolute
    stacked coordinate. Plotting at ``base + replicate / 2`` places each dot at
    the center of the component-specific replicate segment. Callers should pass
    replicate-specific bases when earlier stacked components vary by replicate.
    Signed stacks may pass separate positive and negative bases so each dot is
    placed on the same sign stack as its own replicate value.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes containing the stacked bar chart.
    x_position : float
        Center x-coordinate of the condition bar.
    bottom_value : float
        Aggregate stack baseline for the current segment.
    replicate_values : sequence of Any
        Component-specific per-replicate segment heights.
    plot_settings : PlotSettings
        Plot configuration used for dot styling.
    replicate_base_values : sequence of Any, optional
        Per-replicate cumulative stack bases for unsigned stacks. When omitted,
        ``bottom_value`` is used for every replicate for backward compatibility.
    positive_base_values : sequence of Any, optional
        Per-replicate cumulative positive stack bases for signed stacks.
    negative_base_values : sequence of Any, optional
        Per-replicate cumulative negative stack bases for signed stacks.
    bar_width : float, optional
        Width used for deterministic jitter, by default ``0.8``.
    dot_color : Any, optional
        Override for theme dot colour.
    dot_size : float, optional
        Override for theme dot size.
    dot_alpha : float, optional
        Override for theme dot alpha.
    placement : {"center", "end"}, optional
        Dot placement within each replicate segment. ``"center"`` uses
        ``base + replicate / 2`` and ``"end"`` uses ``base + replicate``.
    zorder : float, optional
        Matplotlib z-order for dot overlays, by default ``5``.

    Returns
    -------
    int
        Number of scatter calls emitted.

    Raises
    ------
    ValueError
        If replicate base arrays do not align with ``replicate_values``.
    """
    import math

    import numpy as np

    if placement not in {"center", "end"}:
        raise ValueError("placement must be 'center' or 'end'")

    raw_values = list(replicate_values)
    if positive_base_values is not None or negative_base_values is not None:
        if positive_base_values is None or negative_base_values is None:
            raise ValueError(
                "positive_base_values and negative_base_values must be provided together"
            )
        positive_bases = list(positive_base_values)
        negative_bases = list(negative_base_values)
        if len(positive_bases) != len(raw_values) or len(negative_bases) != len(raw_values):
            raise ValueError("signed replicate base lengths must match replicate_values length")
        base_values: list[float] = []
        segment_values_list: list[float] = []
        for value, positive_base, negative_base in zip(raw_values, positive_bases, negative_bases):
            try:
                segment_value = float(value)
                base_value = float(positive_base if segment_value >= 0.0 else negative_base)
            except (TypeError, ValueError):
                continue
            if math.isfinite(segment_value) and math.isfinite(base_value):
                segment_values_list.append(segment_value)
                base_values.append(base_value)
        segment_values = np.asarray(segment_values_list, dtype=float)
        bases = np.asarray(base_values, dtype=float)
    elif replicate_base_values is not None:
        raw_bases = list(replicate_base_values)
        if len(raw_bases) != len(raw_values):
            raise ValueError("replicate_base_values length must match replicate_values length")
        base_values = []
        segment_values_list = []
        for value, base in zip(raw_values, raw_bases):
            try:
                segment_value = float(value)
                base_value = float(base)
            except (TypeError, ValueError):
                continue
            if math.isfinite(segment_value) and math.isfinite(base_value):
                segment_values_list.append(segment_value)
                base_values.append(base_value)
        segment_values = np.asarray(segment_values_list, dtype=float)
        bases = np.asarray(base_values, dtype=float)
    else:
        segment_values = finite_numeric_values(raw_values)
        bases = np.full(segment_values.shape, float(bottom_value), dtype=float)

    if segment_values.size == 0:
        return 0

    divisor = 2.0 if placement == "center" else 1.0
    segment_positions = [
        float(base) + float(value) / divisor for base, value in zip(bases, segment_values)
    ]
    return scatter_replicate_values(
        ax,
        [x_position],
        [segment_positions],
        plot_settings,
        orientation="vertical",
        bar_width=bar_width,
        dot_color=dot_color,
        dot_size=dot_size,
        dot_alpha=dot_alpha,
        zorder=zorder,
    )


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
            errors_for_plot = errors
            if replicate_values is not None:
                errors_for_plot = suppress_singleton_errors(errors, replicate_values[i])
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
