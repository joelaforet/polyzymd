"""Structured exceptions for the analysis orchestration lifecycle.

These exceptions provide explicit failure categories while preserving
the existing behavior where expected per-condition failures can be
handled gracefully by higher-level orchestration.
"""


class AnalysisError(Exception):
    """Base class for analysis lifecycle errors."""


class PluginContractError(AnalysisError):
    """Raised when a plugin violates the Analysis contract."""


class ReplicateError(AnalysisError):
    """Raised when per-replicate computation fails unexpectedly."""


class AggregationError(AnalysisError):
    """Raised when condition-level aggregation fails unexpectedly."""


class ComparisonError(AnalysisError):
    """Raised when cross-condition comparison fails."""


class PlotError(AnalysisError):
    """Raised when plot generation fails."""
