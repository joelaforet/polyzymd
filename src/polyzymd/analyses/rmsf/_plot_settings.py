"""RMSF plot settings."""

from __future__ import annotations

from pydantic import Field

from polyzymd.analyses.base import BasePlotSettings


class RMSFPlotSettings(BasePlotSettings):
    """RMSF-specific plot customization.

    Attributes
    ----------
    show_error : bool
        Show error bands/bars on plots (default True)
    highlight_residues : list[int]
        Residue numbers to highlight with vertical lines (e.g., active site)
    figsize_profile : tuple[float, float]
        Figure size for per-residue profile plots
    figsize_comparison : tuple[float, float]
        Figure size for bar comparison plots
    """

    show_error: bool = True
    highlight_residues: list[int] = Field(default_factory=list)
    figsize_profile: tuple[float, float] = (14, 4)
    figsize_comparison: tuple[float, float] = (8, 6)
