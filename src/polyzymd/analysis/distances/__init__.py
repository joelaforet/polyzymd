"""Distance analysis module.

This module provides tools for computing inter-atomic distances
from MD trajectories with proper statistical handling.

Classes
-------
DistanceCalculator
    Main class for distance computation and aggregation

Functions
---------
plot_distance_histogram
    Plot distance distribution
plot_distance_comparison
    Compare multiple distance distributions
plot_distance_timeseries
    Plot distance over time
plot_contact_fraction_bar
    Bar chart of contact fractions
"""

from polyzymd.analysis.distances.calculator import DistanceCalculator
from polyzymd.analysis.distances.plotting import (
    plot_contact_fraction_bar,
    plot_distance_comparison,
    plot_distance_histogram,
    plot_distance_timeseries,
    save_distance_plot,
)

__all__ = [
    "DistanceCalculator",
    "plot_distance_histogram",
    "plot_distance_comparison",
    "plot_distance_timeseries",
    "plot_contact_fraction_bar",
    "save_distance_plot",
]
