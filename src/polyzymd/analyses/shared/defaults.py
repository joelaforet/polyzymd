"""Default parameters shared across all analyses.

This module provides ``AnalysisDefaults``, a small Pydantic model that
holds parameters applied to every analysis type (e.g., equilibration time).

Provides shared defaults so that ``compare/`` can import common
parameters without depending on plugin internals.
"""

from __future__ import annotations

from pydantic import BaseModel, Field


class AnalysisDefaults(BaseModel):
    """Default parameters applied to all analyses.

    Attributes
    ----------
    equilibration_time : str
        Time to skip for equilibration (e.g., "10ns", "5000ps")
    fdr_alpha : float
        False discovery rate threshold for Benjamini-Hochberg correction
        in default scalar pairwise comparisons.
    """

    equilibration_time: str = "10ns"
    fdr_alpha: float = Field(default=0.05, gt=0.0, le=1.0)
