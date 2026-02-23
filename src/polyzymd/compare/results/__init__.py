"""Result models for comparison analysis.

This module provides structured result models for each comparison type.
All new result classes inherit from BaseComparisonResult.

Submodules
----------
rmsf : RMSF comparison results (new OOP-compliant)
rmsf_legacy : Legacy RMSF results (for backward compatibility)
triad : Catalytic triad comparison results
contacts : Polymer-protein contacts comparison results
binding_free_energy : Binding free energy (ΔΔG) comparison results
"""

# Binding free energy result classes
from polyzymd.compare.results.binding_free_energy import (
    BindingFreeEnergyResult,
    FreeEnergyConditionSummary,
    FreeEnergyEntry,
    FreeEnergyPairwiseEntry,
)

# New OOP-compliant result classes
# Contacts result classes
from polyzymd.compare.results.contacts import (
    AggregateComparisonResult,
    ContactsANOVASummary,
    ContactsComparisonResult,
    ContactsConditionSummary,
    ContactsPairwiseComparison,
)

# Distance result classes
from polyzymd.compare.results.distances import (
    DistanceComparisonResult,
    DistanceConditionSummary,
    DistancePairANOVA,
    DistancePairSummary,
    DistancePairwiseComparison,
)
from polyzymd.compare.results.rmsf import RMSFComparisonResult, RMSFConditionSummary

# Legacy RMSF result classes (for backward compatibility with old comparator)
from polyzymd.compare.results.rmsf_legacy import (
    ANOVASummary,
    ComparisonResult,
    ConditionSummary,
    PairwiseComparison,
)

# Triad result classes
from polyzymd.compare.results.triad import (
    TriadANOVASummary,
    TriadComparisonResult,
    TriadConditionSummary,
    TriadPairSummary,
    TriadPairwiseComparison,
)

__all__ = [
    # New OOP RMSF
    "RMSFComparisonResult",
    "RMSFConditionSummary",
    # Legacy RMSF
    "ANOVASummary",
    "ComparisonResult",
    "ConditionSummary",
    "PairwiseComparison",
    # Triad
    "TriadANOVASummary",
    "TriadComparisonResult",
    "TriadConditionSummary",
    "TriadPairSummary",
    "TriadPairwiseComparison",
    # Contacts
    "AggregateComparisonResult",
    "ContactsANOVASummary",
    "ContactsComparisonResult",
    "ContactsConditionSummary",
    "ContactsPairwiseComparison",
    # Distances
    "DistanceComparisonResult",
    "DistanceConditionSummary",
    "DistancePairANOVA",
    "DistancePairSummary",
    "DistancePairwiseComparison",
    # Binding free energy
    "BindingFreeEnergyResult",
    "FreeEnergyConditionSummary",
    "FreeEnergyEntry",
    "FreeEnergyPairwiseEntry",
]
