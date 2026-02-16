"""Comparator implementations.

This module provides comparator classes that inherit from BaseComparator.
Each comparator is automatically registered with ComparatorRegistry.
"""

from polyzymd.compare.comparators.rmsf import RMSFComparator

__all__ = [
    "RMSFComparator",
]
