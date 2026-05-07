"""Binding preference metrics for polymer-protein contacts.

Shared compute layer consumed by contacts binding-preference analysis.
"""

from __future__ import annotations

from ._aggregate import (
    aggregate_binding_preference,
    aggregate_polymer_binding_preference,
    aggregate_system_coverage,
)
from ._compute import compute_binding_preference
from ._models import (
    AggregatedBindingPreferenceEntry,
    AggregatedBindingPreferenceResult,
    AggregatedPartitionBindingEntry,
    AggregatedPartitionBindingResult,
    AggregatedPartitionCoverageEntry,
    AggregatedPartitionCoverageResult,
    AggregatedPolymerBindingPreferenceResult,
    AggregatedSystemCoverageResult,
    BindingPreferenceEntry,
    BindingPreferenceResult,
    PartitionBindingEntry,
    PartitionBindingResult,
    PartitionCoverageEntry,
    PartitionCoverageResult,
    PolymerBindingPreferenceResult,
    PolymerComposition,
    SystemCoverageEntry,
    SystemCoverageResult,
)
from ._orchestration import (
    compute_binding_preference_from_config,
    compute_condition_binding_preference,
)
from ._polymer import extract_polymer_composition
from ._resolution import (
    resolve_polymer_type_selections,
    resolve_protein_group_selections,
    resolve_protein_groups_from_surface_exposure,
)

aggregate_binding_preference_results = aggregate_binding_preference

__all__ = [
    "PolymerComposition",
    "BindingPreferenceEntry",
    "BindingPreferenceResult",
    "compute_binding_preference",
    "aggregate_binding_preference",
    "aggregate_binding_preference_results",
    "AggregatedBindingPreferenceEntry",
    "AggregatedBindingPreferenceResult",
    "SystemCoverageEntry",
    "PartitionCoverageEntry",
    "PartitionCoverageResult",
    "SystemCoverageResult",
    "PartitionBindingEntry",
    "PartitionBindingResult",
    "PolymerBindingPreferenceResult",
    "AggregatedPartitionBindingEntry",
    "AggregatedPartitionBindingResult",
    "AggregatedPolymerBindingPreferenceResult",
    "aggregate_polymer_binding_preference",
    "AggregatedPartitionCoverageEntry",
    "AggregatedPartitionCoverageResult",
    "AggregatedSystemCoverageResult",
    "aggregate_system_coverage",
    "resolve_protein_group_selections",
    "resolve_protein_groups_from_surface_exposure",
    "resolve_polymer_type_selections",
    "extract_polymer_composition",
    "compute_binding_preference_from_config",
    "compute_condition_binding_preference",
]
