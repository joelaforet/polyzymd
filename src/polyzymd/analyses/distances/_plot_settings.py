"""Distances plot settings."""

from __future__ import annotations

from polyzymd.analyses.base import BasePlotSettings


class DistancesPlotSettings(BasePlotSettings):
    """Distance analysis plot customization.

    Attributes
    ----------
    show_threshold : bool
        Show threshold line on distribution plots
    use_kde : bool
        Use KDE instead of histogram for distributions
    generate_state_bars : bool
        Generate per-pair state bar charts (above/below threshold).
        Each pair gets its own figure showing the fraction of frames
        in each state per condition. Default True.
    figsize : tuple[float, float]
        Default figure size for distance plots
    """

    show_threshold: bool = True
    use_kde: bool = True
    generate_state_bars: bool = True
    figsize: tuple[float, float] = (10, 6)
