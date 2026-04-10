"""Secondary structure plot settings."""

from __future__ import annotations

from polyzymd.analyses.base import BasePlotSettings


class SSPlotSettings(BasePlotSettings):
    """Secondary structure plot customization.

    Attributes
    ----------
    generate_timeline : bool
        Generate per-condition residue x time SS heatmap. Default True.
    generate_content_bars : bool
        Generate grouped bar chart of helix/strand/coil fractions. Default True.
    generate_individual_bars : bool
        Generate one bar chart per SS type (helix, beta-sheet, no-SS). Default True.
    generate_diff_heatmap : bool
        Generate condition x residue persistence difference heatmap. Default True.
    figsize_timeline : tuple[float, float]
        Figure size for timeline heatmap.
    figsize_content_bars : tuple[float, float]
        Figure size for content bar chart.
    figsize_diff_heatmap : tuple[float, float] | None
        Figure size for difference heatmap (auto-calculated if None).
    diff_colormap : str
        Diverging colormap for difference heatmap.
    """

    generate_timeline: bool = True
    generate_content_bars: bool = True
    generate_individual_bars: bool = True
    generate_diff_heatmap: bool = True
    figsize_timeline: tuple[float, float] = (14, 6)
    figsize_content_bars: tuple[float, float] = (10, 6)
    figsize_diff_heatmap: tuple[float, float] | None = None
    diff_colormap: str = "RdBu_r"
