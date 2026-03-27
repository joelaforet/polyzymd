"""Default parameters shared across all analyses.

This module provides ``AnalysisDefaults``, a small Pydantic model that
holds parameters applied to every analysis type (e.g., equilibration time).

Previously lived in ``polyzymd.analysis.config``.  Relocated here during
the Phase 3 OCP-compliance refactor so that ``compare/`` can import it
without depending on the legacy ``analysis/`` package.
"""

from __future__ import annotations

from pydantic import BaseModel


class AnalysisDefaults(BaseModel):
    """Default parameters applied to all analyses.

    Attributes
    ----------
    equilibration_time : str
        Time to skip for equilibration (e.g., "10ns", "5000ps")
    """

    equilibration_time: str = "10ns"
