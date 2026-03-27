"""Tests for analyzer registry activation and compatibility.

NOTE: Many tests in this file validated the old ``analysis/`` registry
bootstrap pattern which has been replaced by the ``analyses/`` plugin
system. Those tests are now updated to import from canonical locations.
Tests that validate the old AnalyzerRegistry bootstrap are retained with
imports from the migrated ``analyses/_calculator_*`` modules and
``compare.registries``.
"""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

from polyzymd.analyses._calculator_distances import DistanceCalculator
from polyzymd.analyses._calculator_rmsf import RMSFCalculator
from polyzymd.analyses._calculator_triad import CatalyticTriadAnalyzer
from polyzymd.compare.registries import AnalyzerRegistry, BaseAnalyzer


class TestAnalyzerRegistry:
    """Test analyzer registration, discovery, and factory behavior."""

    def test_all_analyzers_registered_after_bootstrap(self) -> None:
        """Bootstrap should register all expected analyzer keys."""
        # Import triggers @register decorators on calculator classes
        from polyzymd.analyses._calculator_distances import DistanceCalculator  # noqa: F811
        from polyzymd.analyses._calculator_rmsf import RMSFCalculator  # noqa: F811
        from polyzymd.analyses._calculator_triad import CatalyticTriadAnalyzer  # noqa: F811

        expected = {
            "rmsf",
            "distances",
            "catalytic_triad",
        }
        registered = set(AnalyzerRegistry.list_available())
        assert expected.issubset(registered)

    def test_get_returns_expected_classes(self) -> None:
        """Registry lookup should return the concrete analyzer classes."""
        assert AnalyzerRegistry.get("rmsf") is RMSFCalculator
        assert AnalyzerRegistry.get("distances") is DistanceCalculator
        assert AnalyzerRegistry.get("catalytic_triad") is CatalyticTriadAnalyzer

    def test_is_registered_for_known_and_unknown(self) -> None:
        """Registration checks should be correct for known and unknown keys."""
        assert AnalyzerRegistry.is_registered("rmsf")
        assert not AnalyzerRegistry.is_registered("unknown")

    def test_registered_analyzers_are_baseanalyzer_subclasses(self) -> None:
        """All migrated analyzers should satisfy the BaseAnalyzer contract."""
        keys = ["rmsf", "distances", "catalytic_triad"]
        for key in keys:
            analyzer_class = AnalyzerRegistry.get(key)
            assert issubclass(analyzer_class, BaseAnalyzer)

    def test_analysis_type_values(self) -> None:
        """Each analyzer should report its expected analysis type."""
        assert RMSFCalculator.analysis_type() == "rmsf"
        assert DistanceCalculator.analysis_type() == "distances"
        assert CatalyticTriadAnalyzer.analysis_type() == "catalytic_triad"

    def test_baseanalyzer_default_label(self) -> None:
        """Default label should mirror analysis_type when not overridden."""

        class _DummyAnalyzer(BaseAnalyzer):
            @classmethod
            def analysis_type(cls) -> str:
                return "dummy"

            @classmethod
            def from_config(cls, analysis_settings, sim_config, equilibration="0ns"):
                return cls()

            def compute(self, replicate: int, **kwargs):
                return None

            def compute_aggregated(self, replicates, **kwargs):
                return None

        analyzer = _DummyAnalyzer()
        assert analyzer.label == "dummy"

    def test_from_config_builds_instances(self) -> None:
        """Factory methods should map settings to constructor kwargs."""
        sim_config = object()

        rmsf_settings = SimpleNamespace(
            selection="protein and name CA",
            reference_mode="centroid",
            reference_frame=None,
            reference_file=None,
        )
        with patch.object(RMSFCalculator, "__init__", return_value=None) as init_mock:
            result = RMSFCalculator.from_config(rmsf_settings, sim_config, equilibration="10ns")
            assert isinstance(result, RMSFCalculator)
            init_mock.assert_called_once_with(
                config=sim_config,
                selection="protein and name CA",
                equilibration="10ns",
                reference_mode="centroid",
                reference_frame=None,
                reference_file=None,
            )

        distances_settings = SimpleNamespace(
            pairs=[
                SimpleNamespace(selection_a="protein and resid 1", selection_b="chainID C"),
            ]
        )
        with patch.object(DistanceCalculator, "__init__", return_value=None) as init_mock:
            result = DistanceCalculator.from_config(
                distances_settings, sim_config, equilibration="20ns"
            )
            assert isinstance(result, DistanceCalculator)
            init_mock.assert_called_once_with(
                config=sim_config,
                pairs=[("protein and resid 1", "chainID C")],
                equilibration="20ns",
            )

        triad_settings = SimpleNamespace(
            name="triad",
            threshold=3.5,
            pairs=[
                SimpleNamespace(
                    label="a-b",
                    selection_a="protein and resid 77 and name OG",
                    selection_b="protein and resid 156 and name NE2",
                )
            ],
        )
        with patch.object(CatalyticTriadAnalyzer, "__init__", return_value=None) as init_mock:
            result = CatalyticTriadAnalyzer.from_config(
                triad_settings,
                sim_config,
                equilibration="30ns",
            )
            assert isinstance(result, CatalyticTriadAnalyzer)
            init_mock.assert_called_once()
            called_kwargs = init_mock.call_args.kwargs
            assert called_kwargs["config"] is sim_config
            assert called_kwargs["equilibration"] == "30ns"
            assert called_kwargs["triad_config"].name == "triad"

    def test_backward_compatible_analysis_imports(self) -> None:
        """Calculator class names should remain correct."""
        assert RMSFCalculator.__name__ == "RMSFCalculator"
        assert DistanceCalculator.__name__ == "DistanceCalculator"
        assert CatalyticTriadAnalyzer.__name__ == "CatalyticTriadAnalyzer"


class TestAnalysisConfigSecondaryStructure:
    """Test bridge support for legacy analysis config."""

    def test_secondary_structure_in_analysis_config(self) -> None:
        """AnalysisConfig should include SecondaryStructureConfig field."""
        from polyzymd.analyses._analysis_config import AnalysisConfig, SecondaryStructureConfig

        config = AnalysisConfig()

        assert isinstance(config.secondary_structure, SecondaryStructureConfig)
        assert config.secondary_structure.enabled is False
        assert config.secondary_structure.chain_id == "A"

    def test_get_enabled_analyses_includes_secondary_structure(self) -> None:
        """Enabled analyses should include secondary_structure when enabled."""
        from polyzymd.analyses._analysis_config import AnalysisConfig

        config = AnalysisConfig(secondary_structure={"enabled": True, "chain_id": "A"})
        enabled = config.get_enabled_analyses()

        assert "secondary_structure" in enabled
