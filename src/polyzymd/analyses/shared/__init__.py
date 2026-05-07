"""Shared infrastructure for analysis plugins.

This package provides a curated set of commonly used shared utilities.
For advanced or less common helpers, import directly from the relevant
shared submodules.

In particular, selectors, custom selections, and result I/O helpers are not
re-exported from this package root and should be imported from:

- ``polyzymd.analyses.shared.selectors``
- ``polyzymd.analyses.shared.selections``
- ``polyzymd.analyses.shared.result_io``

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
aggregation
    Replicate collection, distance pair aggregation.
autocorrelation
    ACF, correlation time, statistical inefficiency.
pbc
    Minimum image distance, PBC-aware distance matrix.
selections
    Extended selection syntax (midpoint, COM), position retrieval.
diagnostics
    Selection diagnostics, equilibration validation.
window
    Centralized trajectory window resolution for runner-based analyses.
config_hash
    Config hashing for cache validation.
constants
    Shared default values (cutoffs, thresholds).
defaults
    AnalysisDefaults model (equilibration_time).
logging_utils
    Colored terminal logging.
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
    get_independent_indices,
    n_effective,
    statistical_inefficiency,
    statistical_inefficiency_multiple,
)
from polyzymd.analyses.shared.config_hash import (
    compute_config_hash,
    validate_config_hash,
)
from polyzymd.analyses.shared.constants import (
    DEFAULT_CONTACT_CUTOFF,
    DEFAULT_DISTANCE_THRESHOLD,
)
from polyzymd.analyses.shared.convergence import (
    ConvergenceResult,
    find_convergence_time,
)
from polyzymd.analyses.shared.defaults import AnalysisDefaults
from polyzymd.analyses.shared.loader import (
    TrajectoryInfo,
    TrajectoryLoader,
    convert_time,
    parse_time_string,
    time_to_frame,
)
from polyzymd.analyses.shared.pbc import (
    minimum_image_distance,
    pairwise_distances_pbc,
)
from polyzymd.analyses.shared.plotting import (
    annotate_cells,
    apply_axis_style,
    apply_legend,
    find_json,
    get_colors,
    get_output_path,
    get_theme,
    grouped_bars,
    save_figure,
    symmetric_clim,
)
from polyzymd.analyses.shared.statistics import (
    PerResidueStats,
    StatResult,
    aggregate_per_residue_stats,
    aggregate_region_stats,
    compute_sem,
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
    "StatResult",
    "PerResidueStats",
    "compute_sem",
    "aggregate_per_residue_stats",
    "aggregate_region_stats",
    "weighted_mean_with_sem",
    # Autocorrelation
    "ACFResult",
    "CorrelationTimeResult",
    "MIN_RECOMMENDED_N_INDEPENDENT",
    "compute_acf",
    "estimate_correlation_time",
    "get_independent_indices",
    "statistical_inefficiency",
    "statistical_inefficiency_multiple",
    "n_effective",
    "check_statistical_reliability",
    # PBC
    "minimum_image_distance",
    "pairwise_distances_pbc",
    # Constants
    "DEFAULT_CONTACT_CUTOFF",
    "DEFAULT_DISTANCE_THRESHOLD",
    # Defaults
    "AnalysisDefaults",
    # Config hash
    "compute_config_hash",
    "validate_config_hash",
    # Convergence
    "ConvergenceResult",
    "find_convergence_time",
    # Plotting
    "get_theme",
    "apply_axis_style",
    "apply_legend",
    "get_colors",
    "get_output_path",
    "save_figure",
    "grouped_bars",
    "find_json",
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
