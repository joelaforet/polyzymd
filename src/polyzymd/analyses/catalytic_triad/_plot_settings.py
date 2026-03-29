"""Catalytic triad plot settings — registered with PlotSettingsRegistry at import time."""

from __future__ import annotations

from polyzymd.compare.registries import BasePlotSettings, PlotSettingsRegistry


@PlotSettingsRegistry.register("triad")
class TriadPlotSettings(BasePlotSettings):
    """Triad-specific plot customization.

    Attributes
    ----------
    generate_kde_panel : bool
        Generate multi-row KDE panel plot (default True)
    generate_bars : bool
        Generate grouped threshold bar chart (default True)
    generate_2d_kde : bool
        Generate 2D joint KDE plot (default False, more specialized)
    threshold_line_color : str
        Color for threshold vertical line
    kde_fill_alpha : float
        Transparency for KDE fill (0-1)
    kde_xlim : tuple[float, float]
        X-axis limits for KDE panel in Angstroms (default ``(0, 7)``).
    figsize_kde_panel : tuple[float, float] | None
        Figure size for KDE panel (auto-calculated if None)
    figsize_bars : tuple[float, float]
        Figure size for bar chart
    """

    generate_kde_panel: bool = True
    generate_bars: bool = True
    generate_2d_kde: bool = False
    threshold_line_color: str = "red"
    kde_fill_alpha: float = 0.7
    kde_xlim: tuple[float, float] = (0.0, 7.0)
    figsize_kde_panel: tuple[float, float] | None = None
    figsize_bars: tuple[float, float] = (10, 6)
