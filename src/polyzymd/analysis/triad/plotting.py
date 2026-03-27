"""Backward-compatibility shim — canonical location is ``polyzymd.analyses._triad_plotting``.

All public names are re-exported so that existing ``from polyzymd.analysis.triad.plotting
import …`` statements continue to work during the migration period.
"""

from polyzymd.analyses._triad_plotting import (  # noqa: F401
    plot_triad_2d_kde,
    plot_triad_2d_kde_comparison,
    plot_triad_kde_panel,
    plot_triad_kde_panel_pooled,
    plot_triad_threshold_bars,
)

__all__ = [
    "plot_triad_kde_panel",
    "plot_triad_kde_panel_pooled",
    "plot_triad_threshold_bars",
    "plot_triad_2d_kde",
    "plot_triad_2d_kde_comparison",
]
