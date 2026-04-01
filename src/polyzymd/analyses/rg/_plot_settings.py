"""Rg plot settings — registered with PlotSettingsRegistry at import time."""

from __future__ import annotations

from polyzymd.compare.registries import BasePlotSettings, PlotSettingsRegistry


@PlotSettingsRegistry.register("rg")
class RgPlotSettings(BasePlotSettings):
    """Rg analysis plot customization.

    Attributes
    ----------
    show_per_replicate : bool
        Overlay individual replicate traces on timeseries plots.
    figsize : tuple[float, float]
        Default figure size for Rg plots.
    timeseries_figsize : tuple[float, float]
        Figure size for Rg vs time plots (often wider).
    """

    show_per_replicate: bool = False
    figsize: tuple[float, float] = (10, 6)
    timeseries_figsize: tuple[float, float] = (12, 5)
