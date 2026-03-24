"""Catalytic triad/active site analysis module.

This module provides tools for analyzing catalytic triad and active site
geometry from MD trajectories.

Classes
-------
CatalyticTriadAnalyzer
    Main analyzer class for computing triad distances and contact fractions

Functions
---------
plot_triad_kde_panel
    Create multi-row KDE panel comparing triad distances across conditions
plot_triad_kde_panel_pooled
    Create KDE panel from pre-pooled distance arrays
plot_triad_threshold_bars
    Create grouped bar chart of threshold fractions
plot_triad_2d_kde
    Create 2D joint KDE plot of two triad pairs
plot_triad_2d_kde_comparison
    Create 2D joint KDE comparing multiple conditions
"""

from polyzymd.analysis.triad.analyzer import CatalyticTriadAnalyzer
from polyzymd.analysis.triad.plotting import (
    plot_triad_2d_kde,
    plot_triad_2d_kde_comparison,
    plot_triad_kde_panel,
    plot_triad_kde_panel_pooled,
    plot_triad_threshold_bars,
)

__all__ = [
    "CatalyticTriadAnalyzer",
    "plot_triad_kde_panel",
    "plot_triad_kde_panel_pooled",
    "plot_triad_threshold_bars",
    "plot_triad_2d_kde",
    "plot_triad_2d_kde_comparison",
]
