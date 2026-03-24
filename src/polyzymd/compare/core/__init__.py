"""Core infrastructure for comparison analysis.

This module provides base classes and registries for building comparators
following the Open-Closed Principle and Template Method pattern.

Classes
-------
BaseComparator
    Abstract base class for all comparators.
BaseComparisonResult
    Abstract base class for comparison results.
BaseConditionSummary
    Abstract base class for condition summaries.
PairwiseComparison
    Shared model for pairwise statistical comparisons.
ANOVASummary
    Shared model for ANOVA results.
ComparatorRegistry
    Registry for dynamically creating comparators.

Example
-------
Creating a custom comparator:

>>> from polyzymd.compare.core import (
...     BaseComparator,
...     BaseComparisonResult,
...     BaseConditionSummary,
...     ComparatorRegistry,
... )
>>>
>>> @ComparatorRegistry.register("my_metric")
... class MyComparator(BaseComparator):
...     ...
"""

from polyzymd.compare.core.base import (
    ANOVASummary,
    BaseComparator,
    BaseComparisonResult,
    BaseConditionSummary,
    MetricType,
    PairwiseComparison,
)
from polyzymd.compare.core.registry import ComparatorRegistry

__all__ = [
    "ANOVASummary",
    "BaseComparator",
    "BaseComparisonResult",
    "BaseConditionSummary",
    "ComparatorRegistry",
    "MetricType",
    "PairwiseComparison",
]
