"""Registry pattern for analysis types — backward-compatibility shim.

All classes have been relocated to :mod:`polyzymd.compare.registries`.
This module re-exports them so that existing ``analysis/`` calculators
and tests continue to work until they are migrated or deleted.

.. deprecated::
    Import from ``polyzymd.compare.registries`` instead.
"""

from __future__ import annotations

# Re-export everything from the new authoritative location
from polyzymd.compare.registries import (  # noqa: F401
    AnalysisSettingsRegistry,
    AnalyzerRegistry,
    BaseAnalysisSettings,
    BaseAnalyzer,
    BaseComparisonSettings,
    BasePlotSettings,
    ComparisonSettingsRegistry,
    PlotSettingsRegistry,
)

__all__ = [
    "AnalysisSettingsRegistry",
    "AnalyzerRegistry",
    "BaseAnalysisSettings",
    "BaseAnalyzer",
    "BaseComparisonSettings",
    "BasePlotSettings",
    "ComparisonSettingsRegistry",
    "PlotSettingsRegistry",
]
