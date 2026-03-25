"""Tests for BaseComparator factory hooks and ComparatorRegistry config-driven creation.

Phase 3 of OCP refactoring: verifies that comparison_type_name() is concrete,
config_settings_key() resolves correctly, and create_from_config() works for
all registered comparators.
"""

from __future__ import annotations

import importlib
import sys
from typing import Any, ClassVar
from unittest.mock import MagicMock

import pytest

from polyzymd.analysis.core.metric_type import MetricType
from polyzymd.compare.core.base import BaseComparator, BaseComparisonResult, BaseConditionSummary
from polyzymd.compare.core.registry import ComparatorRegistry


class _FakeConditionSummary(BaseConditionSummary):
    """Minimal concrete condition summary for testing."""

    @property
    def primary_metric_value(self) -> float:
        return 0.0

    @property
    def primary_metric_sem(self) -> float:
        return 0.0


class _FakeResult(BaseComparisonResult):
    comparison_type: ClassVar[str] = "fake"


class _FakeComparatorSimple(BaseComparator):
    """Comparator with defaults and no comparison settings."""

    comparison_type: ClassVar[str] = "fake_simple"

    @property
    def metric_type(self) -> MetricType:
        return MetricType.MEAN_BASED

    def _load_or_compute(self, cond, recompute):
        return {}

    def _build_condition_summary(self, cond, data):
        return _FakeConditionSummary(
            label=cond.label,
            config_path="",
            n_replicates=1,
            replicate_values=[1.0],
        )


class _FakeComparatorWithAlias(BaseComparator):
    """Comparator with analysis_settings_key override like triad."""

    comparison_type: ClassVar[str] = "fake_alias"
    analysis_settings_key: ClassVar[str] = "fake_catalytic"

    @property
    def metric_type(self) -> MetricType:
        return MetricType.MEAN_BASED

    def _load_or_compute(self, cond, recompute):
        return {}

    def _build_condition_summary(self, cond, data):
        return _FakeConditionSummary(
            label=cond.label,
            config_path="",
            n_replicates=1,
            replicate_values=[1.0],
        )


class _FakeComparatorWithCompSettings(BaseComparator):
    """Comparator that uses comparison_settings like contacts."""

    comparison_type: ClassVar[str] = "fake_comp_settings"
    uses_comparison_settings: ClassVar[bool] = True

    def __init__(self, config, analysis_settings, equilibration=None, comparison_settings=None):
        super().__init__(config, analysis_settings, equilibration)
        self.comparison_settings = comparison_settings

    @property
    def metric_type(self) -> MetricType:
        return MetricType.MEAN_BASED

    def _load_or_compute(self, cond, recompute):
        return {}

    def _build_condition_summary(self, cond, data):
        return _FakeConditionSummary(
            label=cond.label,
            config_path="",
            n_replicates=1,
            replicate_values=[1.0],
        )


class TestComparisonTypeNameConcrete:
    """Verify comparison_type_name() returns comparison_type without override."""

    def test_returns_comparison_type(self):
        assert _FakeComparatorSimple.comparison_type_name() == "fake_simple"

    def test_alias_comparator_returns_comparison_type(self):
        """comparison_type_name() returns comparison_type, not analysis_settings_key."""
        assert _FakeComparatorWithAlias.comparison_type_name() == "fake_alias"


class TestConfigSettingsKey:
    """Verify config_settings_key() resolves correctly."""

    def test_default_falls_back_to_comparison_type(self):
        assert _FakeComparatorSimple.config_settings_key() == "fake_simple"

    def test_alias_uses_analysis_settings_key(self):
        assert _FakeComparatorWithAlias.config_settings_key() == "fake_catalytic"

    def test_comp_settings_falls_back(self):
        assert _FakeComparatorWithCompSettings.config_settings_key() == "fake_comp_settings"


class TestRealComparatorMetadata:
    """Verify all registered comparators have correct metadata after Phase 3."""

    def test_triad_analysis_settings_key(self):
        from polyzymd.compare.comparators.triad import TriadComparator

        assert TriadComparator.config_settings_key() == "catalytic_triad"
        assert TriadComparator.comparison_type_name() == "triad"

    def test_rmsf_config_settings_key(self):
        from polyzymd.compare.comparators.rmsf import RMSFComparator

        assert RMSFComparator.config_settings_key() == "rmsf"
        assert RMSFComparator.comparison_type_name() == "rmsf"

    def test_contacts_uses_comparison_settings(self):
        from polyzymd.compare.comparators.contacts import ContactsComparator

        assert ContactsComparator.uses_comparison_settings is True
        assert ContactsComparator.config_settings_key() == "contacts"

    def test_exposure_uses_comparison_settings(self):
        from polyzymd.compare.comparators.exposure import ExposureDynamicsComparator

        assert ExposureDynamicsComparator.uses_comparison_settings is True

    def test_bfe_uses_comparison_settings(self):
        from polyzymd.compare.comparators.binding_free_energy import BindingFreeEnergyComparator

        assert BindingFreeEnergyComparator.uses_comparison_settings is True

    def test_polymer_affinity_uses_comparison_settings(self):
        from polyzymd.compare.comparators.polymer_affinity import PolymerAffinityScoreComparator

        assert PolymerAffinityScoreComparator.uses_comparison_settings is True

    def test_distances_defaults(self):
        from polyzymd.compare.comparators.distances import DistancesComparator

        assert DistancesComparator.uses_comparison_settings is False
        assert DistancesComparator.config_settings_key() == "distances"

    def test_secondary_structure_defaults(self):
        from polyzymd.compare.comparators.secondary_structure import SecondaryStructureComparator

        assert SecondaryStructureComparator.uses_comparison_settings is False
        assert SecondaryStructureComparator.config_settings_key() == "secondary_structure"


class TestRegistryGetForAnalysisType:
    """Verify analysis-type lookup works for direct and aliased keys."""

    def test_direct_key(self):
        from polyzymd.compare.comparators.rmsf import RMSFComparator

        cls = ComparatorRegistry.get_for_analysis_type("rmsf")
        assert cls is RMSFComparator

    def test_alias_key_catalytic_triad(self):
        """catalytic_triad should resolve to TriadComparator."""
        from polyzymd.compare.comparators.triad import TriadComparator

        cls = ComparatorRegistry.get_for_analysis_type("catalytic_triad")
        assert cls is TriadComparator

    def test_unknown_key_raises(self):
        with pytest.raises(ValueError, match="No comparator for analysis settings key"):
            ComparatorRegistry.get_for_analysis_type("nonexistent_type")


class TestRegistryCreateFromConfig:
    """Verify config-driven factory creation."""

    def _make_mock_config(
        self,
        settings_key: str,
        analysis_data: Any = None,
        comparison_data: Any = None,
    ):
        """Create a mock ComparisonConfig with the given settings."""
        mock_config = MagicMock()
        mock_config.defaults.equilibration_time = "10ns"

        def analysis_get(key: str):
            if key == settings_key:
                return analysis_data or MagicMock()
            return None

        mock_config.analysis_settings.get = analysis_get

        def comparison_get(key: str):
            if key == settings_key:
                return comparison_data
            return None

        mock_config.comparison_settings.get = comparison_get

        return mock_config

    def test_create_rmsf_from_config(self):
        """RMSF should be created without comparison_settings."""
        from polyzymd.compare.comparators.rmsf import RMSFComparator
        from polyzymd.compare.settings import RMSFAnalysisSettings

        mock_settings = MagicMock(spec=RMSFAnalysisSettings)
        mock_settings.selection = "protein and name CA"
        mock_settings.reference_mode = "centroid"
        mock_settings.reference_frame = None
        mock_settings.reference_file = None

        config = self._make_mock_config("rmsf", analysis_data=mock_settings)
        comparator = ComparatorRegistry.create_from_config("rmsf", config, equilibration="5ns")
        assert isinstance(comparator, RMSFComparator)

    def test_create_triad_via_alias(self):
        """catalytic_triad settings key should resolve to TriadComparator."""
        from polyzymd.compare.comparators.triad import TriadComparator

        config = self._make_mock_config("catalytic_triad")
        comparator = ComparatorRegistry.create_from_config("catalytic_triad", config)
        assert isinstance(comparator, TriadComparator)

    def test_missing_settings_raises(self):
        """Should raise ValueError when analysis_settings section is missing."""
        config = self._make_mock_config("nonexistent")
        with pytest.raises(ValueError, match="No comparator for analysis settings key"):
            ComparatorRegistry.create_from_config("nonexistent", config)

    def test_comparison_settings_injected(self):
        """Comparators with uses_comparison_settings=True get comparison_settings."""
        from polyzymd.compare.comparators.contacts import ContactsComparator

        mock_comp_settings = MagicMock()
        config = self._make_mock_config("contacts", comparison_data=mock_comp_settings)
        comparator = ComparatorRegistry.create_from_config("contacts", config)
        assert isinstance(comparator, ContactsComparator)
        assert comparator.comparison_settings is mock_comp_settings

    def test_no_comparison_settings_for_simple(self):
        """Comparators without uses_comparison_settings don't get it."""
        from polyzymd.compare.comparators.distances import DistancesComparator
        from polyzymd.compare.settings import DistancesAnalysisSettings

        mock_settings = MagicMock(spec=DistancesAnalysisSettings)
        mock_settings.model_dump.return_value = {}
        config = self._make_mock_config("distances", analysis_data=mock_settings)
        comparator = ComparatorRegistry.create_from_config("distances", config)
        assert isinstance(comparator, DistancesComparator)
        assert not hasattr(comparator, "comparison_settings")


class TestEnsureAllComparatorsRegistered:
    """Verify that the CLI bootstrapper populates the comparator registry."""

    def test_all_expected_comparators_registered(self):
        """Bootstrapper should register all expected comparison types."""
        from polyzymd.compare.cli import _ensure_all_comparators_registered

        ComparatorRegistry.clear()

        comparator_modules = [
            mod_name
            for mod_name in list(sys.modules)
            if mod_name == "polyzymd.compare.comparators"
            or mod_name.startswith("polyzymd.compare.comparators.")
        ]
        for mod_name in comparator_modules:
            sys.modules.pop(mod_name, None)

        importlib.invalidate_caches()
        _ensure_all_comparators_registered()

        available = ComparatorRegistry.list_available()
        expected = {
            "rmsf",
            "triad",
            "contacts",
            "distances",
            "secondary_structure",
            "exposure",
            "binding_free_energy",
            "polymer_affinity",
        }
        assert expected.issubset(set(available)), (
            f"Missing comparators: {expected - set(available)}"
        )

    def test_catalytic_triad_resolvable_after_registration(self):
        """catalytic_triad should resolve to TriadComparator after bootstrap."""
        from polyzymd.compare.cli import _ensure_all_comparators_registered

        _ensure_all_comparators_registered()

        cls = ComparatorRegistry.get_for_analysis_type("catalytic_triad")
        assert cls.comparison_type_name() == "triad"
