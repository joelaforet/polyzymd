"""Polymer affinity plot settings."""

from __future__ import annotations

from polyzymd.analyses.base import BasePlotSettings


class AffinityPlotSettings(BasePlotSettings):
    """Polymer affinity score plot customization.

    Attributes
    ----------
    generate_stacked_bars : bool
        Generate stacked bar chart of total score by condition, broken
        down by polymer type. Default True.
    generate_group_bars : bool
        Generate grouped bar chart showing per-group contributions
        across conditions. Default True.
    figsize_stacked : tuple[float, float]
        Figure size for stacked bar chart.
    figsize_group_bars : tuple[float, float]
        Figure size for grouped bar charts.
    show_error_bars : bool
        Show SEM error bars on plots. Default True.
    """

    generate_stacked_bars: bool = True
    generate_group_bars: bool = True
    figsize_stacked: tuple[float, float] = (10, 6)
    figsize_group_bars: tuple[float, float] = (10, 6)
    show_error_bars: bool = True
