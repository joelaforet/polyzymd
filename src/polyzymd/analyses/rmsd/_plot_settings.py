"""RMSD plot settings — registered with PlotSettingsRegistry at import time."""

from __future__ import annotations

from polyzymd.compare.registries import BasePlotSettings, PlotSettingsRegistry


@PlotSettingsRegistry.register("rmsd")
class RMSDPlotSettings(BasePlotSettings):
    """RMSD analysis plot customization.

    Attributes
    ----------
    show_per_replicate : bool
        Overlay individual replicate traces on timeseries plots.
    figsize : tuple[float, float]
        Default figure size for RMSD plots.
    timeseries_figsize : tuple[float, float]
        Figure size for RMSD vs time plots (often wider).
    show_convergence_plots : bool
        Generate per-replicate convergence diagnostic plots.
    convergence_figsize : tuple[float, float]
        Figure size for convergence diagnostic panels.
    """

    show_per_replicate: bool = False
    figsize: tuple[float, float] = (10, 6)
    timeseries_figsize: tuple[float, float] = (12, 5)
    show_convergence_plots: bool = False
    convergence_figsize: tuple[float, float] = (12.0, 5.0)
