"""Binding free energy plot settings."""

from __future__ import annotations

from polyzymd.analyses.base import BasePlotSettings


class BFEPlotSettings(BasePlotSettings):
    """Binding free energy plot customization.

    Attributes
    ----------
    generate_heatmap : bool
        Generate ΔG_sel heatmap (rows = AA groups, columns = conditions). Default True.
    generate_bars : bool
        Generate ΔG_sel grouped bar chart (one bar per condition per AA group). Default True.
    figsize_heatmap : tuple[float, float] | None
        Figure size for ΔG_sel heatmap (auto-calculated if None).
    figsize_bars : tuple[float, float]
        Figure size for ΔG_sel bar charts.
    colormap : str
        Diverging colormap for heatmap (default "RdBu_r": red = avoidance, blue = preference).
    show_error_bars : bool
        Show SEM error bars on bar charts. Default True.
    annotate_heatmap : bool
        Annotate each heatmap cell with its ΔG_sel value. Default True.
    """

    generate_heatmap: bool = True
    generate_bars: bool = True
    figsize_heatmap: tuple[float, float] | None = None
    figsize_bars: tuple[float, float] = (10, 6)
    colormap: str = "RdBu_r"
    show_error_bars: bool = True
    annotate_heatmap: bool = True
