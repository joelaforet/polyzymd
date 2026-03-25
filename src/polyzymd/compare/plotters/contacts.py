"""Contacts analysis plotters — backward-compatible re-export shim.

All plotter implementations have been split into focused modules:
- contacts_binding_preference: BindingPreferenceHeatmapPlotter, BindingPreferenceBarPlotter
- contacts_coverage: SystemCoverageHeatmapPlotter, SystemCoverageBarPlotter,
  UserPartitionBarPlotter
- contacts_profiles: ContactFractionProfilePlotter, ResidenceTimeProfilePlotter
- contacts_grouped_bars: CF/RT by AA class and partition bar plotters

Importing this module registers all contacts plotters with PlotterRegistry.
"""

from __future__ import annotations

from polyzymd.compare.plotters._contacts_shared import (
    _get_enrichment_value,
    _get_enrichment_with_sem,
    _get_polymer_types_and_aa_classes,
    _load_aggregated_contact_results,
    _load_binding_preference_results,
    _load_partition_definitions,
    _load_system_coverage_results,
)
from polyzymd.compare.plotters.contacts_binding_preference import (
    BindingPreferenceBarPlotter,
    BindingPreferenceHeatmapPlotter,
)
from polyzymd.compare.plotters.contacts_coverage import (
    SystemCoverageBarPlotter,
    SystemCoverageHeatmapPlotter,
    UserPartitionBarPlotter,
)
from polyzymd.compare.plotters.contacts_grouped_bars import (
    ContactFractionByAAClassBarPlotter,
    ContactFractionByPartitionBarPlotter,
    ResidenceTimeByAAClassBarPlotter,
    ResidenceTimeByPartitionBarPlotter,
)
from polyzymd.compare.plotters.contacts_profiles import (
    ContactFractionProfilePlotter,
    ResidenceTimeProfilePlotter,
)

__all__ = [
    "BindingPreferenceHeatmapPlotter",
    "BindingPreferenceBarPlotter",
    "SystemCoverageHeatmapPlotter",
    "SystemCoverageBarPlotter",
    "UserPartitionBarPlotter",
    "ContactFractionProfilePlotter",
    "ResidenceTimeProfilePlotter",
    "ContactFractionByAAClassBarPlotter",
    "ContactFractionByPartitionBarPlotter",
    "ResidenceTimeByAAClassBarPlotter",
    "ResidenceTimeByPartitionBarPlotter",
    "_get_polymer_types_and_aa_classes",
    "_get_enrichment_value",
    "_get_enrichment_with_sem",
    "_load_binding_preference_results",
    "_load_system_coverage_results",
    "_load_aggregated_contact_results",
    "_load_partition_definitions",
]
