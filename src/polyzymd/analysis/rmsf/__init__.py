"""RMSF analysis module.

This module provides tools for computing Root Mean Square Fluctuation (RMSF)
from MD trajectories with proper statistical handling.

Classes
-------
RMSFCalculator
    Main class for RMSF computation and aggregation

Functions
---------
plot_rmsf
    Plot per-residue RMSF
plot_rmsf_comparison
    Compare multiple RMSF results
plot_rmsf_heatmap
    Heatmap of RMSF across replicates
"""

from polyzymd.analysis.rmsf.calculator import RMSFCalculator
from polyzymd.analysis.rmsf.plotting import (
    plot_rmsf,
    plot_rmsf_comparison,
    plot_rmsf_heatmap,
    save_rmsf_plot,
)

__all__ = [
    "RMSFCalculator",
    "plot_rmsf",
    "plot_rmsf_comparison",
    "plot_rmsf_heatmap",
    "save_rmsf_plot",
]
