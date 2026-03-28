"""Result models for comparison analysis.

.. versionchanged:: 1.3.0
    All plugin-specific result models now live in their respective
    analysis plugin packages under ``polyzymd.analyses.<plugin>._comparison_results``.
    This module re-exports them for convenience.
"""

# Binding free energy result classes
from polyzymd.analyses.binding_free_energy._comparison_results import (
    BindingFreeEnergyResult,
    FreeEnergyConditionSummary,
    FreeEnergyEntry,
    FreeEnergyPairwiseEntry,
)

# Contacts result classes
from polyzymd.analyses.contacts._comparison_results import (
    AggregateComparisonResult,
    ContactsANOVASummary,
    ContactsComparisonResult,
    ContactsConditionSummary,
    ContactsPairwiseComparison,
)

# Distance result classes
from polyzymd.analyses.distances._comparison_results import (
    DistanceComparisonResult,
    DistanceConditionSummary,
    DistancePairANOVA,
    DistancePairSummary,
    DistancePairwiseComparison,
)

# Exposure result classes
from polyzymd.analyses.exposure._comparison_results import (
    ExposureComparisonResult,
    ExposureConditionSummary,
)

# Polymer affinity score result classes
from polyzymd.analyses.polymer_affinity._comparison_results import (
    AffinityScoreConditionSummary,
    AffinityScoreEntry,
    AffinityScorePairwiseEntry,
    PolymerAffinityScoreResult,
    PolymerTypeScore,
)

# RMSF result classes (new OOP-compliant)
from polyzymd.analyses.rmsf._comparison_results import (
    RMSFComparisonResult,
    RMSFConditionSummary,
)

# Secondary structure result classes
from polyzymd.analyses.secondary_structure._comparison_results import (
    SSComparisonResult,
    SSConditionSummary,
)

# Triad result classes
from polyzymd.analyses.catalytic_triad._comparison_results import (
    TriadANOVASummary,
    TriadComparisonResult,
    TriadConditionSummary,
    TriadPairSummary,
    TriadPairwiseComparison,
)

__all__ = [
    # RMSF
    "RMSFComparisonResult",
    "RMSFConditionSummary",
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
    # Exposure
    "ExposureComparisonResult",
    "ExposureConditionSummary",
    # Binding free energy
    "BindingFreeEnergyResult",
    "FreeEnergyConditionSummary",
    "FreeEnergyEntry",
    "FreeEnergyPairwiseEntry",
    # Polymer affinity score
    "AffinityScoreConditionSummary",
    "AffinityScoreEntry",
    "AffinityScorePairwiseEntry",
    "PolymerAffinityScoreResult",
    "PolymerTypeScore",
    # Secondary structure
    "SSComparisonResult",
    "SSConditionSummary",
]
