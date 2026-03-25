"""Tests for split contacts plotter modules and compatibility shim."""

from __future__ import annotations

from polyzymd.compare.config import PlotSettings
from polyzymd.compare.plotter import BasePlotter, PlotterRegistry
from polyzymd.compare.plotters.contacts import (
    BindingPreferenceBarPlotter,
    BindingPreferenceHeatmapPlotter,
    ContactFractionByAAClassBarPlotter,
    ContactFractionByPartitionBarPlotter,
    ContactFractionProfilePlotter,
    ResidenceTimeByAAClassBarPlotter,
    ResidenceTimeByPartitionBarPlotter,
    ResidenceTimeProfilePlotter,
    SystemCoverageBarPlotter,
    SystemCoverageHeatmapPlotter,
    UserPartitionBarPlotter,
    _get_enrichment_value,
    _get_enrichment_with_sem,
    _get_polymer_types_and_aa_classes,
    _load_aggregated_contact_results,
    _load_binding_preference_results,
    _load_partition_definitions,
    _load_system_coverage_results,
)
from polyzymd.compare.plotters.contacts_binding_preference import (
    BindingPreferenceBarPlotter as BindingPreferenceBarPlotterDirect,
)
from polyzymd.compare.plotters.contacts_binding_preference import (
    BindingPreferenceHeatmapPlotter as BindingPreferenceHeatmapPlotterDirect,
)
from polyzymd.compare.plotters.contacts_coverage import (
    SystemCoverageBarPlotter as SystemCoverageBarPlotterDirect,
)
from polyzymd.compare.plotters.contacts_coverage import (
    SystemCoverageHeatmapPlotter as SystemCoverageHeatmapPlotterDirect,
)
from polyzymd.compare.plotters.contacts_coverage import (
    UserPartitionBarPlotter as UserPartitionBarPlotterDirect,
)
from polyzymd.compare.plotters.contacts_grouped_bars import (
    ContactFractionByAAClassBarPlotter as ContactFractionByAAClassBarPlotterDirect,
)
from polyzymd.compare.plotters.contacts_grouped_bars import (
    ContactFractionByPartitionBarPlotter as ContactFractionByPartitionBarPlotterDirect,
)
from polyzymd.compare.plotters.contacts_grouped_bars import (
    ResidenceTimeByAAClassBarPlotter as ResidenceTimeByAAClassBarPlotterDirect,
)
from polyzymd.compare.plotters.contacts_grouped_bars import (
    ResidenceTimeByPartitionBarPlotter as ResidenceTimeByPartitionBarPlotterDirect,
)
from polyzymd.compare.plotters.contacts_profiles import (
    ContactFractionProfilePlotter as ContactFractionProfilePlotterDirect,
)
from polyzymd.compare.plotters.contacts_profiles import (
    ResidenceTimeProfilePlotter as ResidenceTimeProfilePlotterDirect,
)

EXPECTED_CONTACTS_PLOT_TYPES = [
    "binding_preference_heatmap",
    "binding_preference_bars",
    "system_coverage_heatmap",
    "system_coverage_bars",
    "user_partition_bars",
    "contact_fraction_profile",
    "residence_time_profile",
    "cf_by_aa_class_bars",
    "cf_by_partition_bars",
    "rt_by_aa_class_bars",
    "rt_by_partition_bars",
]

PLOTTER_CLASSES = [
    BindingPreferenceHeatmapPlotter,
    BindingPreferenceBarPlotter,
    SystemCoverageHeatmapPlotter,
    SystemCoverageBarPlotter,
    UserPartitionBarPlotter,
    ContactFractionProfilePlotter,
    ResidenceTimeProfilePlotter,
    ContactFractionByAAClassBarPlotter,
    ContactFractionByPartitionBarPlotter,
    ResidenceTimeByAAClassBarPlotter,
    ResidenceTimeByPartitionBarPlotter,
]

DIRECT_IMPORT_CLASSES = [
    BindingPreferenceHeatmapPlotterDirect,
    BindingPreferenceBarPlotterDirect,
    SystemCoverageHeatmapPlotterDirect,
    SystemCoverageBarPlotterDirect,
    UserPartitionBarPlotterDirect,
    ContactFractionProfilePlotterDirect,
    ResidenceTimeProfilePlotterDirect,
    ContactFractionByAAClassBarPlotterDirect,
    ContactFractionByPartitionBarPlotterDirect,
    ResidenceTimeByAAClassBarPlotterDirect,
    ResidenceTimeByPartitionBarPlotterDirect,
]


def test_all_contacts_plot_types_are_registered() -> None:
    """All contacts plotter keys should remain registered in PlotterRegistry."""
    for plot_type in EXPECTED_CONTACTS_PLOT_TYPES:
        assert PlotterRegistry.is_registered(plot_type), (
            f"Expected registered plot type: {plot_type}"
        )


def test_contacts_shim_re_exports_plotter_classes_and_helpers() -> None:
    """Contacts shim should re-export all classes and helper functions."""
    exports = [
        BindingPreferenceHeatmapPlotter,
        BindingPreferenceBarPlotter,
        SystemCoverageHeatmapPlotter,
        SystemCoverageBarPlotter,
        UserPartitionBarPlotter,
        ContactFractionProfilePlotter,
        ResidenceTimeProfilePlotter,
        ContactFractionByAAClassBarPlotter,
        ContactFractionByPartitionBarPlotter,
        ResidenceTimeByAAClassBarPlotter,
        ResidenceTimeByPartitionBarPlotter,
        _get_polymer_types_and_aa_classes,
        _get_enrichment_value,
        _get_enrichment_with_sem,
        _load_binding_preference_results,
        _load_system_coverage_results,
        _load_aggregated_contact_results,
        _load_partition_definitions,
    ]
    assert all(export is not None for export in exports)


def test_contacts_plotter_classes_are_subclasses_of_base_plotter() -> None:
    """All contacts plotter classes should inherit from BasePlotter."""
    for cls in PLOTTER_CLASSES:
        assert issubclass(cls, BasePlotter)


def test_contacts_plot_type_strings_match_expected() -> None:
    """Each contacts plotter class should report the expected plot type string."""
    observed = [cls.plot_type() for cls in PLOTTER_CLASSES]
    assert observed == EXPECTED_CONTACTS_PLOT_TYPES


def test_contacts_plotter_can_plot_guards_non_contacts_analysis() -> None:
    """Contacts plotters should reject non-contacts analysis types."""
    settings = PlotSettings()
    for cls in PLOTTER_CLASSES:
        plotter = cls(settings)
        assert not plotter.can_plot(comparison_config=None, analysis_type="rmsf")


def test_direct_module_imports_expose_same_plotter_classes() -> None:
    """Direct imports from split modules should expose all plotter classes."""
    assert len(DIRECT_IMPORT_CLASSES) == len(PLOTTER_CLASSES)
    for direct_cls, shim_cls in zip(DIRECT_IMPORT_CLASSES, PLOTTER_CLASSES):
        assert direct_cls is shim_cls
