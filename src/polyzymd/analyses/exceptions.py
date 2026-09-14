"""Structured exceptions for the analysis orchestration lifecycle.

These exceptions provide explicit failure categories while preserving
the existing behavior where expected per-condition failures can be
handled gracefully by higher-level orchestration.
"""


class AnalysisError(Exception):
    """Base class for analysis lifecycle errors."""


class PluginContractError(AnalysisError):
    """Raised when a plugin violates the Analysis contract."""


class ReplicateSkippedError(AnalysisError):
    """Raised when a replicate is skipped for a known recoverable reason."""


class ReplicateError(AnalysisError):
    """Raised when per-replicate computation fails unexpectedly."""


class AggregationError(AnalysisError):
    """Raised when condition-level aggregation fails unexpectedly."""


class ComparisonError(AnalysisError):
    """Raised when cross-condition comparison fails."""


class PlotError(AnalysisError):
    """Raised when plot generation fails."""


class DependencyError(AnalysisError):
    """Raised when declared analysis dependencies are invalid or missing."""


class StatisticsError(AnalysisError, ValueError):
    """Raised when a statistical estimator is given an input it cannot use.

    An invalid sample or an unsupported option is a typed failure, not a
    silently degraded zero. It also subclasses ``ValueError`` so that callers
    written before the typed error existed keep working.
    """

class SelectionError(AnalysisError):
    """Raised when an analysis selection is empty or cannot be resolved."""
