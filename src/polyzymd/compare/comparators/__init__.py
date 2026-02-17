"""Comparator implementations.

This module provides comparator classes that inherit from BaseComparator.
Each comparator is automatically registered with ComparatorRegistry.
"""

from polyzymd.compare.comparators.contacts import ContactsComparator
from polyzymd.compare.comparators.distances import DistancesComparator
from polyzymd.compare.comparators.rmsf import RMSFComparator
from polyzymd.compare.comparators.triad import TriadComparator

__all__ = [
    "ContactsComparator",
    "DistancesComparator",
    "RMSFComparator",
    "TriadComparator",
]
