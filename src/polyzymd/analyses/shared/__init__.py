"""Reusable shared utilities for analysis plugins.

This package provides contributor-facing utilities that are broadly reusable
across analysis plugins. Framework internals, CLI helpers, and plugin-private
artifact helpers live with their owning packages.

In particular, selectors and custom selections are not re-exported from this
package root and should be imported from:

- ``polyzymd.analyses.shared.selectors``
- ``polyzymd.analyses.shared.selections``

Sub-modules
-----------
loader
    Trajectory loading, time parsing, frame conversion.
alignment
    Trajectory alignment (centroid/average/frame/external).
centroid
    K-Means centroid frame finding, reference mode dispatch.
statistics
    SEM, per-residue/region aggregation, weighted mean.
autocorrelation
    ACF, correlation time, statistical inefficiency.
selections
    Extended selection syntax (midpoint, COM), position retrieval.
diagnostics
    Selection diagnostics, equilibration validation.
window
    Centralized trajectory window resolution for MDAnalysis job lifecycles.
plotting
    Shared plotting helpers (axis styling, figure saving, grouped bars, etc.).
"""

from __future__ import annotations

# Re-export the most commonly used symbols for convenience.
# Plugins can do:  from polyzymd.analyses.shared import TrajectoryLoader, AlignmentConfig
# or import specific sub-modules directly.
from polyzymd.analyses.shared.alignment import (
    AlignmentConfig,
    ReferenceMode,
    align_trajectory,
    get_alignment_description,
)
from polyzymd.analyses.shared.autocorrelation import (
    MIN_RECOMMENDED_N_INDEPENDENT,
    ACFResult,
    CorrelationTimeResult,
    check_statistical_reliability,
    compute_acf,
    estimate_correlation_time,
    n_effective,
    statistical_inefficiency,
    statistical_inefficiency_multiple,
)
from polyzymd.analyses.shared.loader import (
    TrajectoryInfo,
    TrajectoryLoader,
    convert_time,
    parse_time_string,
    time_to_frame,
)
from polyzymd.analyses.shared.plotting import (
    add_uncertainty_footnote,
    annotate_cells,
    apply_axis_style,
    apply_legend,
    band_half_widths,
    error_bar_half_widths,
    get_condition_color_map,
    get_condition_colors,
    get_output_path,
    get_palette_colors,
    get_theme,
    grouped_bars,
    order_condition_labels,
    plugin_plot_settings,
    resolve_error_bar,
    save_figure,
    symmetric_clim,
)
from polyzymd.analyses.shared.statistics import (
    CI_METHOD_STUDENT_T,
    UNCERTAINTY_KIND_REPLICATE_SEM,
    MeanSemCI,
    PerResidueStats,
    StatResult,
    aggregate_per_residue_stats,
    aggregate_region_stats,
    compute_sem,
    mean_sem_ci,
    student_t_coverage_factor,
    uncertainty_block,
    weighted_mean_with_sem,
)
from polyzymd.analyses.shared.window import (
    TrajectoryWindow,
    resolve_replicate_trajectory_window,
    resolve_trajectory_window,
)

__all__ = [
    # Loader
    "TrajectoryInfo",
    "TrajectoryLoader",
    "parse_time_string",
    "convert_time",
    "time_to_frame",
    "TrajectoryWindow",
    "resolve_trajectory_window",
    "resolve_replicate_trajectory_window",
    # Alignment
    "AlignmentConfig",
    "ReferenceMode",
    "align_trajectory",
    "get_alignment_description",
    # Statistics
    "CI_METHOD_STUDENT_T",
    "UNCERTAINTY_KIND_REPLICATE_SEM",
    "MeanSemCI",
    "StatResult",
    "PerResidueStats",
    "compute_sem",
    "mean_sem_ci",
    "student_t_coverage_factor",
    "uncertainty_block",
    "aggregate_per_residue_stats",
    "aggregate_region_stats",
    "weighted_mean_with_sem",
    # Autocorrelation
    "ACFResult",
    "CorrelationTimeResult",
    "MIN_RECOMMENDED_N_INDEPENDENT",
    "compute_acf",
    "estimate_correlation_time",
    "statistical_inefficiency",
    "statistical_inefficiency_multiple",
    "n_effective",
    "check_statistical_reliability",
    # Plotting
    "get_theme",
    "apply_axis_style",
    "apply_legend",
    "get_palette_colors",
    "order_condition_labels",
    "get_condition_colors",
    "get_condition_color_map",
    "add_uncertainty_footnote",
    "band_half_widths",
    "error_bar_half_widths",
    "plugin_plot_settings",
    "resolve_error_bar",
    "get_output_path",
    "save_figure",
    "grouped_bars",
    "annotate_cells",
    "symmetric_clim",
    # Plot settings (lazily re-exported from config.comparison)
    "PlotSettings",
]


def __getattr__(name: str):
    """Lazily expose ``PlotSettings`` without creating an import cycle."""
    if name == "PlotSettings":
        from polyzymd.config.comparison import PlotSettings

        return PlotSettings
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
