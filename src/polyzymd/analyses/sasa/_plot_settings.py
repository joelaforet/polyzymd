"""SASA plot settings — registered with PlotSettingsRegistry at import time."""

from __future__ import annotations

from polyzymd.compare.registries import BasePlotSettings, PlotSettingsRegistry


@PlotSettingsRegistry.register("sasa")
class SASAPlotSettings(BasePlotSettings):
    """SASA analysis plot customization."""

    show_per_replicate: bool = False
    figsize: tuple[float, float] = (10, 6)
    timeseries_figsize: tuple[float, float] = (12, 5)
    profile_figsize: tuple[float, float] = (12, 5)
