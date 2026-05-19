"""RMSF plot settings."""

from __future__ import annotations

from pydantic import Field

from polyzymd.analyses.base import BasePlotSettings


class RMSFPlotSettings(BasePlotSettings):
    """RMSF-specific plot customization."""

    show_error: bool = True
    show_reference_secondary_structure: bool = True
    highlight_residues: list[int] = Field(default_factory=list)
    figsize_profile: tuple[float, float] = (14, 4)
    figsize_comparison: tuple[float, float] = (8, 6)
