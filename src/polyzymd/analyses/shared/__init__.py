"""Building blocks the analyses are made of: loading, windows, statistics and plotting.

Sub-modules
-----------
loader
    Trajectory loading, time parsing, frame conversion.
window
    The production window of a replicate after the equilibration time.
selections
    Extended selection syntax (midpoint, COM), position retrieval.
diagnostics
    Selection diagnostics, equilibration validation.
centroid
    The frame closest to the iterative average structure, for centroid references.
topology
    Checks that a topology carries the bonds an analysis needs.
statistics
    Mean, standard error and Student t interval of replicate values.
inferential_statistics
    t tests, effect sizes and the Benjamini-Hochberg correction.
autocorrelation
    Statistical inefficiency, effective sample size and detected equilibration.
aa_classification
    Maximum accessible surface area of each amino acid.
groupings
    Physicochemical classes of amino acids.
plotting
    Axis styling, condition colours, grouped bars, uncertainty bands, figure saving.
"""

from __future__ import annotations

# Re-export the most commonly used symbols, so that
# ``from polyzymd.analyses.shared import TrajectoryLoader`` works.
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
    apply_axis_style,
    apply_legend,
    band_half_widths,
    error_bar_half_widths,
    get_condition_color_map,
    get_output_path,
    get_palette_colors,
    grouped_bars,
    order_condition_labels,
    save_figure,
)
from polyzymd.analyses.shared.statistics import (
    CI_METHOD_STUDENT_T,
    MeanSemCI,
    mean_sem_ci,
    student_t_coverage_factor,
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
    # Statistics
    "CI_METHOD_STUDENT_T",
    "MeanSemCI",
    "mean_sem_ci",
    "student_t_coverage_factor",
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
    "apply_axis_style",
    "apply_legend",
    "get_palette_colors",
    "order_condition_labels",
    "get_condition_color_map",
    "add_uncertainty_footnote",
    "band_half_widths",
    "error_bar_half_widths",
    "get_output_path",
    "save_figure",
    "grouped_bars",
    # Plot settings (lazily re-exported from config.analysis_settings)
    "PlotSettings",
]


def __getattr__(name: str):
    """Lazily expose ``PlotSettings`` without creating an import cycle."""
    if name == "PlotSettings":
        from polyzymd.config.analysis_settings import PlotSettings

        return PlotSettings
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
