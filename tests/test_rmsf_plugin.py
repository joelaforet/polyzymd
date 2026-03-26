"""Tests for the RMSF analysis plugin (Phase B).

Tests the RMSFAnalysis class: discovery, settings, compute_replicate,
aggregate, extract_metrics, _deserialize_result, and the full lifecycle
via the orchestrator.

Heavy dependencies (MDAnalysis, trajectories) are mocked.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, ClassVar
from unittest.mock import MagicMock, patch

import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import (
    AggregateContext,
    ComparisonContext,
    Condition,
    MetricValue,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.rmsf import RMSFAnalysis, RMSFSettings


# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def rmsf_analysis():
    """Return a fresh RMSFAnalysis instance."""
    return RMSFAnalysis()


@pytest.fixture
def default_settings():
    """Return default RMSFSettings."""
    return RMSFSettings()


@pytest.fixture
def condition():
    """Return a mock Condition."""
    return Condition(
        label="No Polymer",
        config_path=Path("/fake/config.yaml"),
        replicates=(1, 2, 3),
        sim_config=MagicMock(),
    )


def _make_mock_rmsf_result(replicate: int, mean_rmsf: float = 1.5) -> MagicMock:
    """Create a mock RMSFResult with realistic fields."""
    result = MagicMock()
    result.replicate = replicate
    result.rmsf_values = [1.0, 1.5, 2.0, 1.2, 1.8]
    result.residue_ids = [1, 2, 3, 4, 5]
    result.residue_names = ["ALA", "GLY", "VAL", "LEU", "ILE"]
    result.mean_rmsf = mean_rmsf
    result.std_rmsf = 0.35
    result.min_rmsf = 1.0
    result.max_rmsf = 2.0
    result.config_hash = "abc123"
    result.equilibration_time = 10.0
    result.equilibration_unit = "ns"
    result.selection_string = "protein and name CA"
    result.n_independent_frames = 50
    return result


# ============================================================================
# Test: Discovery
# ============================================================================


class TestRMSFDiscovery:
    """Test that RMSFAnalysis is auto-discovered by the plugin system."""

    def test_discovery_finds_rmsf(self):
        from polyzymd.analyses.discovery import clear_cache, get_analysis, list_analyses

        clear_cache()
        analyses = list_analyses()
        assert "rmsf" in analyses
        assert analyses["rmsf"] is RMSFAnalysis

    def test_get_analysis_returns_rmsf(self):
        from polyzymd.analyses.discovery import clear_cache, get_analysis

        clear_cache()
        cls = get_analysis("rmsf")
        assert cls is RMSFAnalysis


# ============================================================================
# Test: Class variables and Settings
# ============================================================================


class TestRMSFClassVars:
    """Test class variables and settings."""

    def test_name(self, rmsf_analysis):
        assert rmsf_analysis.name == "rmsf"

    def test_aliases_empty(self, rmsf_analysis):
        assert rmsf_analysis.aliases == ()

    def test_dependencies_empty(self, rmsf_analysis):
        assert rmsf_analysis.dependencies == ()

    def test_min_replicates(self, rmsf_analysis):
        assert rmsf_analysis.min_replicates == 2

    def test_repr(self, rmsf_analysis):
        assert repr(rmsf_analysis) == "<RMSFAnalysis(name='rmsf')>"


class TestRMSFSettings:
    """Test RMSFSettings validation and defaults."""

    def test_defaults(self):
        s = RMSFSettings()
        assert s.selection == "protein and name CA"
        assert s.reference_mode == "centroid"
        assert s.reference_frame is None
        assert s.reference_file is None

    def test_custom_selection(self):
        s = RMSFSettings(selection="protein and name N")
        assert s.selection == "protein and name N"

    def test_frame_mode(self):
        s = RMSFSettings(reference_mode="frame", reference_frame=500)
        assert s.reference_mode == "frame"
        assert s.reference_frame == 500

    def test_external_mode(self):
        s = RMSFSettings(reference_mode="external", reference_file="/path/to/ref.pdb")
        assert s.reference_mode == "external"

    def test_invalid_reference_mode(self):
        with pytest.raises(ValueError, match="reference_mode must be one of"):
            RMSFSettings(reference_mode="invalid")

    def test_average_mode(self):
        s = RMSFSettings(reference_mode="average")
        assert s.reference_mode == "average"


# ============================================================================
# Test: compute_replicate
# ============================================================================


class TestComputeReplicate:
    """Test RMSFAnalysis.compute_replicate delegates to RMSFCalculator."""

    @patch("polyzymd.analysis.rmsf.calculator.RMSFCalculator")
    def test_delegates_to_calculator(self, MockCalc, rmsf_analysis, condition, tmp_path):
        mock_result = _make_mock_rmsf_result(replicate=1)
        mock_calc_instance = MagicMock()
        mock_calc_instance.compute.return_value = mock_result
        MockCalc.return_value = mock_calc_instance

        settings = RMSFSettings()
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=settings,
        )

        result = rmsf_analysis.compute_replicate(ctx, 1)

        assert result is mock_result
        MockCalc.assert_called_once_with(
            config=condition.sim_config,
            selection="protein and name CA",
            equilibration="10ns",
            reference_mode="centroid",
            reference_frame=None,
            reference_file=None,
        )
        mock_calc_instance.compute.assert_called_once_with(
            replicate=1,
            save=True,
            output_dir=tmp_path / "run_1",
            recompute=False,
        )

    @patch("polyzymd.analysis.rmsf.calculator.RMSFCalculator")
    def test_passes_custom_settings(self, MockCalc, rmsf_analysis, condition, tmp_path):
        MockCalc.return_value.compute.return_value = _make_mock_rmsf_result(1)

        settings = RMSFSettings(
            selection="backbone",
            reference_mode="frame",
            reference_frame=100,
        )
        ctx = ReplicateContext(
            condition=condition,
            replicate=2,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_2",
            equilibration="5ns",
            recompute=True,
            settings=settings,
        )

        rmsf_analysis.compute_replicate(ctx, 2)

        MockCalc.assert_called_once_with(
            config=condition.sim_config,
            selection="backbone",
            equilibration="5ns",
            reference_mode="frame",
            reference_frame=100,
            reference_file=None,
        )
        MockCalc.return_value.compute.assert_called_once_with(
            replicate=2,
            save=True,
            output_dir=tmp_path / "run_2",
            recompute=True,
        )

    @patch("polyzymd.analysis.rmsf.calculator.RMSFCalculator")
    def test_handles_legacy_settings(self, MockCalc, rmsf_analysis, condition, tmp_path):
        """Legacy RMSFAnalysisSettings from old config should work via getattr."""
        MockCalc.return_value.compute.return_value = _make_mock_rmsf_result(1)

        # Simulate legacy settings with same attributes
        legacy_settings = MagicMock()
        legacy_settings.selection = "protein and name CA"
        legacy_settings.reference_mode = "average"
        legacy_settings.reference_frame = None
        legacy_settings.reference_file = None

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=legacy_settings,
        )

        rmsf_analysis.compute_replicate(ctx, 1)

        MockCalc.assert_called_once_with(
            config=condition.sim_config,
            selection="protein and name CA",
            equilibration="10ns",
            reference_mode="average",
            reference_frame=None,
            reference_file=None,
        )


# ============================================================================
# Test: aggregate
# ============================================================================


class TestAggregate:
    """Test RMSFAnalysis.aggregate."""

    def test_aggregates_results(self, rmsf_analysis, condition, tmp_path):
        """Test aggregate produces an RMSFAggregatedResult with correct stats."""
        results = [_make_mock_rmsf_result(i, mean_rmsf=1.0 + 0.1 * i) for i in range(1, 4)]
        agg_dir = tmp_path / "aggregated"

        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=agg_dir,
            equilibration="10ns",
            settings=RMSFSettings(),
        )

        # Patch only get_polyzymd_version (non-essential) and let real functions run
        with patch("polyzymd.analysis.results.base.get_polyzymd_version", return_value="1.2.1"):
            result = rmsf_analysis.aggregate(ctx, results)

        # Verify the result has correct attributes
        assert result.n_replicates == 3
        assert result.replicates == [1, 2, 3]
        assert len(result.mean_rmsf_per_residue) == 5
        assert len(result.sem_rmsf_per_residue) == 5
        assert len(result.per_replicate_mean_rmsf) == 3
        assert result.overall_mean_rmsf == pytest.approx(1.2, abs=0.01)  # mean of 1.1, 1.2, 1.3

    def test_aggregate_saves_result_file(self, rmsf_analysis, condition, tmp_path):
        """Verify aggregated result is saved to the output directory."""
        results = [_make_mock_rmsf_result(i, mean_rmsf=1.5) for i in range(1, 3)]
        agg_dir = tmp_path / "aggregated"

        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=agg_dir,
            equilibration="10ns",
            settings=RMSFSettings(),
        )

        with patch("polyzymd.analysis.results.base.get_polyzymd_version", return_value="1.2.1"):
            rmsf_analysis.aggregate(ctx, results)

        # Check that a file was saved
        json_files = list(agg_dir.glob("*.json"))
        assert len(json_files) == 1
        assert "rmsf_reps1-2_eq10ns.json" in json_files[0].name

    def test_aggregate_uses_correct_replicate_means(self, rmsf_analysis, condition, tmp_path):
        """Verify that per_replicate_mean_rmsf contains the correct values."""
        results = [_make_mock_rmsf_result(i, mean_rmsf=1.0 + 0.1 * i) for i in range(1, 4)]

        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=RMSFSettings(),
        )

        with patch("polyzymd.analysis.results.base.get_polyzymd_version", return_value="1.2.1"):
            result = rmsf_analysis.aggregate(ctx, results)

        assert result.per_replicate_mean_rmsf == [1.1, 1.2, 1.3]


# ============================================================================
# Test: extract_metrics
# ============================================================================


class TestExtractMetrics:
    """Test RMSFAnalysis.extract_metrics for default comparison pipeline."""

    def test_returns_mean_rmsf_metric(self, rmsf_analysis):
        summary = MagicMock()
        summary.overall_mean_rmsf = 1.45
        summary.overall_sem_rmsf = 0.08
        summary.per_replicate_mean_rmsf = [1.4, 1.5, 1.45]

        metrics = rmsf_analysis.extract_metrics(summary)

        assert "mean_rmsf" in metrics
        m = metrics["mean_rmsf"]
        assert isinstance(m, MetricValue)
        assert m.name == "mean_rmsf"
        assert m.mean == 1.45
        assert m.sem == 0.08
        assert m.replicate_values == [1.4, 1.5, 1.45]
        assert m.higher_is_better is False
        assert m.direction_labels == ("stabilizing", "unchanged", "destabilizing")

    def test_returns_single_metric(self, rmsf_analysis):
        """RMSF should return exactly one metric."""
        summary = MagicMock()
        summary.overall_mean_rmsf = 1.0
        summary.overall_sem_rmsf = 0.05
        summary.per_replicate_mean_rmsf = [1.0, 1.0]

        metrics = rmsf_analysis.extract_metrics(summary)
        assert len(metrics) == 1


# ============================================================================
# Test: _deserialize_result
# ============================================================================


class TestDeserializeResult:
    """Test RMSFAnalysis._deserialize_result."""

    def test_loads_aggregated_result(self, rmsf_analysis, tmp_path):
        """Test that _deserialize_result loads via RMSFAggregatedResult.load()."""
        with patch("polyzymd.analysis.results.rmsf.RMSFAggregatedResult") as MockAgg:
            mock_loaded = MagicMock()
            MockAgg.load.return_value = mock_loaded

            result = rmsf_analysis._deserialize_result(tmp_path / "test.json")

            MockAgg.load.assert_called_once_with(tmp_path / "test.json")
            assert result is mock_loaded


# ============================================================================
# Test: filter_conditions (default behavior)
# ============================================================================


class TestFilterConditions:
    """RMSF applies to all conditions — filter should keep all."""

    def test_keeps_all_conditions(self, rmsf_analysis):
        conditions = [
            Condition(
                label=f"Cond {i}",
                config_path=Path(f"/fake/config_{i}.yaml"),
                replicates=(1, 2),
                sim_config=MagicMock(),
            )
            for i in range(4)
        ]

        filtered = rmsf_analysis.filter_conditions(conditions)
        assert len(filtered) == 4


# ============================================================================
# Test: _make_aggregated_filename
# ============================================================================


class TestMakeAggregatedFilename:
    """Test backward-compatible filename generation."""

    def test_contiguous_replicates(self):
        result = MagicMock()
        result.equilibration_time = 10.0
        result.equilibration_unit = "ns"

        filename = RMSFAnalysis._make_aggregated_filename((1, 2, 3, 4, 5), result)
        assert filename == "rmsf_reps1-5_eq10ns.json"

    def test_non_contiguous_replicates(self):
        result = MagicMock()
        result.equilibration_time = 100.0
        result.equilibration_unit = "ns"

        filename = RMSFAnalysis._make_aggregated_filename((1, 3, 5), result)
        assert filename == "rmsf_reps1_3_5_eq100ns.json"

    def test_single_pair(self):
        result = MagicMock()
        result.equilibration_time = 0.0
        result.equilibration_unit = "ns"

        filename = RMSFAnalysis._make_aggregated_filename((1, 2), result)
        assert filename == "rmsf_reps1-2_eq0ns.json"


# ============================================================================
# Test: plot
# ============================================================================


class TestPlot:
    """Test RMSFAnalysis.plot delegates to existing plotters."""

    def test_plot_returns_empty_on_no_conditions(self, rmsf_analysis, tmp_path):
        ctx = PlotContext(
            conditions=[],
            analysis_dirs={},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=RMSFSettings(),
            plot_settings=None,
        )

        plots = rmsf_analysis.plot(ctx)
        assert plots == []

    @patch("polyzymd.compare.plotters.rmsf.RMSFProfilePlotter")
    @patch("polyzymd.compare.plotters.rmsf.RMSFComparisonPlotter")
    def test_plot_delegates_to_plotters(
        self,
        MockCompPlotter,
        MockProfPlotter,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        # Setup mock plotters
        mock_comp = MagicMock()
        mock_comp.plot.return_value = [tmp_path / "figures" / "rmsf_comparison.png"]
        MockCompPlotter.return_value = mock_comp

        mock_prof = MagicMock()
        mock_prof.plot.return_value = [tmp_path / "figures" / "rmsf_profile.png"]
        MockProfPlotter.return_value = mock_prof

        analysis_dir = tmp_path / "analysis" / "no_polymer" / "rmsf"
        analysis_dir.mkdir(parents=True)

        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"No Polymer": analysis_dir},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=RMSFSettings(),
            plot_settings=None,
        )

        with patch("polyzymd.compare.config.PlotSettings") as MockPlotSettings:
            MockPlotSettings.return_value = MagicMock()
            plots = rmsf_analysis.plot(ctx)

        assert len(plots) == 2
        mock_comp.plot.assert_called_once()
        mock_prof.plot.assert_called_once()

    @patch("polyzymd.compare.plotters.rmsf.RMSFProfilePlotter", side_effect=Exception("plot error"))
    @patch(
        "polyzymd.compare.plotters.rmsf.RMSFComparisonPlotter", side_effect=Exception("plot error")
    )
    def test_plot_catches_exceptions(
        self,
        MockCompPlotter,
        MockProfPlotter,
        rmsf_analysis,
        condition,
        tmp_path,
    ):
        """Plotter failures should be caught, not crash the pipeline."""
        analysis_dir = tmp_path / "analysis" / "no_polymer" / "rmsf"
        analysis_dir.mkdir(parents=True)

        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"No Polymer": analysis_dir},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=RMSFSettings(),
            plot_settings=None,
        )

        with patch("polyzymd.compare.config.PlotSettings") as MockPlotSettings:
            MockPlotSettings.return_value = MagicMock()
            plots = rmsf_analysis.plot(ctx)

        # Should not crash, just return empty
        assert plots == []


# ============================================================================
# Test: Full lifecycle via orchestrator (mocked compute)
# ============================================================================


class TestRMSFLifecycle:
    """Test the full compute -> aggregate -> compare lifecycle.

    Uses the orchestrator's run_analysis() and verifies the data flow
    through the RMSF plugin. Heavy deps are mocked.
    """

    @patch("polyzymd.analysis.rmsf.calculator.RMSFCalculator")
    def test_run_analysis_lifecycle(self, MockCalc, rmsf_analysis, condition, tmp_path):
        """Test compute_replicate -> aggregate via run_analysis()."""
        from polyzymd.analyses.orchestrator import run_analysis

        # Create mock results for each replicate
        mock_results = {
            rep: _make_mock_rmsf_result(rep, mean_rmsf=1.0 + 0.1 * rep) for rep in (1, 2, 3)
        }

        def mock_compute(replicate, save=True, output_dir=None, recompute=False):
            return mock_results[replicate]

        MockCalc.return_value.compute = mock_compute

        # Let real aggregation functions run, only mock version
        with patch("polyzymd.analysis.results.base.get_polyzymd_version", return_value="1.2.1"):
            output_dir = tmp_path / "analysis" / "no_polymer" / "rmsf"
            result = run_analysis(
                rmsf_analysis,
                condition,
                settings=RMSFSettings(),
                equilibration="10ns",
                output_dir=output_dir,
            )

        # Verify we got an aggregated result
        assert result.n_replicates == 3
        assert result.per_replicate_mean_rmsf == [1.1, 1.2, 1.3]
        # Calculator should have been called 3 times (one per replicate)
        assert MockCalc.call_count == 3

    def test_extract_metrics_feeds_default_compare(self, rmsf_analysis):
        """Verify extract_metrics output is compatible with default_scalar_comparison."""
        from polyzymd.analyses.stats import default_scalar_comparison

        # Simulate two conditions with extracted metrics
        metrics_by_condition = {}
        for label, mean, sem, vals in [
            ("No Polymer", 1.8, 0.12, [1.7, 1.9, 1.8]),
            ("100% SBMA", 1.3, 0.08, [1.25, 1.35, 1.3]),
        ]:
            summary = MagicMock()
            summary.overall_mean_rmsf = mean
            summary.overall_sem_rmsf = sem
            summary.per_replicate_mean_rmsf = vals
            metrics_by_condition[label] = rmsf_analysis.extract_metrics(summary)

        # Run default comparison — should not crash
        result = default_scalar_comparison(
            analysis_name="rmsf",
            project_name="test_project",
            metrics_by_condition=metrics_by_condition,
            control_label="No Polymer",
            equilibration="10ns",
        )

        assert result["analysis_type"] == "rmsf"
        assert len(result["conditions"]) == 2
        assert len(result["pairwise_comparisons"]) >= 1
        # Ranking: lower RMSF first (100% SBMA should rank before No Polymer)
        assert result["ranking"][0] == "100% SBMA"
