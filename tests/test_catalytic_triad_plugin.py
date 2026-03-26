"""Tests for the catalytic triad analysis plugin.

Tests the CatalyticTriadAnalysis class: discovery, aliases, settings,
compute_replicate, aggregate, extract_metrics, _deserialize_result,
_settings_to_config, _make_aggregated_filename, plot delegation, and the
full lifecycle via the orchestrator.

Heavy dependencies (MDAnalysis, trajectories) are mocked.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any
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
from polyzymd.analyses.catalytic_triad import (
    CatalyticTriadAnalysis,
    CatalyticTriadSettings,
    TriadPairSettings,
)


# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def triad_analysis():
    """Return a fresh CatalyticTriadAnalysis instance."""
    return CatalyticTriadAnalysis()


@pytest.fixture
def default_settings():
    """Return default CatalyticTriadSettings with two pairs."""
    return CatalyticTriadSettings(
        name="LipA_triad",
        pairs=[
            TriadPairSettings(
                label="Asp133-His156",
                selection_a="resid 133 and name OD1",
                selection_b="resid 156 and name ND1",
            ),
            TriadPairSettings(
                label="His156-Ser77",
                selection_a="resid 156 and name NE2",
                selection_b="resid 77 and name OG",
            ),
        ],
        threshold=3.5,
        description="Ser-His-Asp catalytic triad",
    )


@pytest.fixture
def condition():
    """Return a mock Condition."""
    return Condition(
        label="No Polymer",
        config_path=Path("/fake/config.yaml"),
        replicates=(1, 2, 3),
        sim_config=MagicMock(),
    )


def _make_mock_pair_result(pair_idx: int, replicate: int) -> MagicMock:
    """Create a mock DistancePairResult for one pair in one replicate."""
    pr = MagicMock()
    pr.pair_label = f"Pair_{pair_idx}"
    pr.selection1 = f"resid {100 + pair_idx} and name OD1"
    pr.selection2 = f"resid {200 + pair_idx} and name ND1"
    pr.mean_distance = 3.0 + 0.1 * pair_idx
    pr.std_distance = 0.5
    pr.sem_distance = 0.1
    pr.median_distance = 3.0 + 0.1 * pair_idx
    pr.min_distance = 2.0
    pr.max_distance = 5.0
    pr.replicate = replicate
    pr.threshold = 3.5
    pr.fraction_below_threshold = 0.7 - 0.05 * pair_idx
    pr.kde_peak = 2.8 + 0.1 * pair_idx
    pr.distances = [3.0, 3.2, 2.8, 3.5, 2.9]  # Needed for aggregation stats
    return pr


def _make_mock_triad_result(replicate: int, sim_contact: float = 0.65) -> MagicMock:
    """Create a mock TriadResult with realistic fields."""
    result = MagicMock()
    result.replicate = replicate
    result.triad_name = "LipA_triad"
    result.triad_description = "Ser-His-Asp catalytic triad"
    result.pair_results = [
        _make_mock_pair_result(0, replicate),
        _make_mock_pair_result(1, replicate),
    ]
    result.threshold = 3.5
    result.simultaneous_contact_fraction = sim_contact
    result.n_frames_simultaneous = int(sim_contact * 1000)
    result.n_frames_total = 1200
    result.n_frames_used = 1000
    result.config_hash = "hash123"
    result.equilibration_time = 10.0
    result.equilibration_unit = "ns"
    result.selection_string = "(resid 133 : resid 156); (resid 156 : resid 77)"
    return result


# ============================================================================
# Test: Discovery
# ============================================================================


class TestTriadDiscovery:
    """Test that CatalyticTriadAnalysis is auto-discovered by the plugin system."""

    def test_discovery_finds_catalytic_triad(self):
        from polyzymd.analyses.discovery import clear_cache, list_analyses

        clear_cache()
        analyses = list_analyses()
        assert "catalytic_triad" in analyses
        assert analyses["catalytic_triad"] is CatalyticTriadAnalysis

    def test_get_analysis_by_name(self):
        from polyzymd.analyses.discovery import clear_cache, get_analysis

        clear_cache()
        cls = get_analysis("catalytic_triad")
        assert cls is CatalyticTriadAnalysis

    def test_get_analysis_by_alias(self):
        from polyzymd.analyses.discovery import clear_cache, get_analysis

        clear_cache()
        cls = get_analysis("triad")
        assert cls is CatalyticTriadAnalysis

    def test_list_all_names_includes_alias(self):
        from polyzymd.analyses.discovery import clear_cache, list_all_names

        clear_cache()
        all_names = list_all_names()
        assert "catalytic_triad" in all_names
        assert "triad" in all_names


# ============================================================================
# Test: Class variables
# ============================================================================


class TestTriadClassVars:
    """Test class variables and settings."""

    def test_name(self, triad_analysis):
        assert triad_analysis.name == "catalytic_triad"

    def test_aliases(self, triad_analysis):
        assert triad_analysis.aliases == ("triad",)

    def test_dependencies_empty(self, triad_analysis):
        assert triad_analysis.dependencies == ()

    def test_min_replicates(self, triad_analysis):
        assert triad_analysis.min_replicates == 2

    def test_repr(self, triad_analysis):
        assert repr(triad_analysis) == "<CatalyticTriadAnalysis(name='catalytic_triad')>"

    def test_settings_type(self, triad_analysis):
        assert triad_analysis.Settings is CatalyticTriadSettings


# ============================================================================
# Test: Settings validation
# ============================================================================


class TestCatalyticTriadSettings:
    """Test CatalyticTriadSettings validation and defaults."""

    def test_valid_settings(self, default_settings):
        assert default_settings.name == "LipA_triad"
        assert len(default_settings.pairs) == 2
        assert default_settings.threshold == 3.5
        assert default_settings.description == "Ser-His-Asp catalytic triad"

    def test_default_threshold(self):
        s = CatalyticTriadSettings(
            pairs=[
                TriadPairSettings(label="pair1", selection_a="resid 1", selection_b="resid 2"),
            ],
        )
        assert s.threshold == 3.5

    def test_default_name(self):
        s = CatalyticTriadSettings(
            pairs=[
                TriadPairSettings(label="pair1", selection_a="resid 1", selection_b="resid 2"),
            ],
        )
        assert s.name == "catalytic_triad"

    def test_empty_pairs_raises(self):
        with pytest.raises(ValueError, match="At least one distance pair"):
            CatalyticTriadSettings(pairs=[])

    def test_n_pairs(self, default_settings):
        assert default_settings.n_pairs == 2

    def test_get_pair_selections(self, default_settings):
        sels = default_settings.get_pair_selections()
        assert len(sels) == 2
        assert sels[0] == ("resid 133 and name OD1", "resid 156 and name ND1")

    def test_get_pair_labels(self, default_settings):
        labels = default_settings.get_pair_labels()
        assert labels == ["Asp133-His156", "His156-Ser77"]

    def test_custom_threshold(self):
        s = CatalyticTriadSettings(
            name="test",
            threshold=4.0,
            pairs=[
                TriadPairSettings(label="pair1", selection_a="resid 1", selection_b="resid 2"),
            ],
        )
        assert s.threshold == 4.0

    def test_no_description(self):
        s = CatalyticTriadSettings(
            pairs=[
                TriadPairSettings(label="pair1", selection_a="resid 1", selection_b="resid 2"),
            ],
        )
        assert s.description is None


# ============================================================================
# Test: _settings_to_config
# ============================================================================


class TestSettingsToConfig:
    """Test converting plugin settings to CatalyticTriadConfig."""

    def test_converts_plugin_settings(self, default_settings):
        from polyzymd.compare.config import CatalyticTriadConfig

        config = CatalyticTriadAnalysis._settings_to_config(default_settings)

        assert isinstance(config, CatalyticTriadConfig)
        assert config.name == "LipA_triad"
        assert config.threshold == 3.5
        assert len(config.pairs) == 2
        assert config.pairs[0].label == "Asp133-His156"
        assert config.description == "Ser-His-Asp catalytic triad"

    def test_passes_through_config(self):
        """If settings is already a CatalyticTriadConfig, pass through."""
        from polyzymd.compare.config import CatalyticTriadConfig, TriadPairConfig

        config = CatalyticTriadConfig(
            name="test",
            pairs=[
                TriadPairConfig(label="p1", selection_a="resid 1", selection_b="resid 2"),
            ],
        )

        result = CatalyticTriadAnalysis._settings_to_config(config)
        assert result is config


# ============================================================================
# Test: compute_replicate
# ============================================================================


class TestComputeReplicate:
    """Test CatalyticTriadAnalysis.compute_replicate delegates to analyzer."""

    @patch("polyzymd.analysis.triad.analyzer.CatalyticTriadAnalyzer")
    def test_delegates_to_analyzer(
        self, MockAnalyzer, triad_analysis, condition, tmp_path, default_settings
    ):
        mock_result = _make_mock_triad_result(replicate=1)
        mock_analyzer_inst = MagicMock()
        mock_analyzer_inst.compute.return_value = mock_result
        MockAnalyzer.return_value = mock_analyzer_inst

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=default_settings,
        )

        result = triad_analysis.compute_replicate(ctx, 1)

        assert result is mock_result
        MockAnalyzer.assert_called_once()
        mock_analyzer_inst.compute.assert_called_once_with(
            replicate=1,
            save=True,
            output_dir=tmp_path / "run_1",
            recompute=False,
        )

    @patch("polyzymd.analysis.triad.analyzer.CatalyticTriadAnalyzer")
    def test_passes_recompute_flag(
        self, MockAnalyzer, triad_analysis, condition, tmp_path, default_settings
    ):
        MockAnalyzer.return_value.compute.return_value = _make_mock_triad_result(2)

        ctx = ReplicateContext(
            condition=condition,
            replicate=2,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_2",
            equilibration="5ns",
            recompute=True,
            settings=default_settings,
        )

        triad_analysis.compute_replicate(ctx, 2)

        MockAnalyzer.return_value.compute.assert_called_once_with(
            replicate=2,
            save=True,
            output_dir=tmp_path / "run_2",
            recompute=True,
        )

    @patch("polyzymd.analysis.triad.analyzer.CatalyticTriadAnalyzer")
    def test_handles_legacy_settings(self, MockAnalyzer, triad_analysis, condition, tmp_path):
        """Legacy CatalyticTriadAnalysisSettings should work via _settings_to_config."""
        mock_result = _make_mock_triad_result(1)
        MockAnalyzer.return_value.compute.return_value = mock_result

        # Simulate legacy settings object
        legacy = MagicMock()
        legacy.name = "LipA_triad"
        legacy.threshold = 3.5
        legacy.description = None
        legacy_pair = MagicMock()
        legacy_pair.label = "pair1"
        legacy_pair.selection_a = "resid 1"
        legacy_pair.selection_b = "resid 2"
        legacy.pairs = [legacy_pair]

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=legacy,
        )

        result = triad_analysis.compute_replicate(ctx, 1)
        assert result is mock_result


# ============================================================================
# Test: aggregate
# ============================================================================


class TestAggregate:
    """Test CatalyticTriadAnalysis.aggregate."""

    def test_aggregates_results(self, triad_analysis, condition, tmp_path, default_settings):
        """Test aggregate produces a TriadAggregatedResult with correct stats."""
        results = [_make_mock_triad_result(i, sim_contact=0.6 + 0.05 * i) for i in range(1, 4)]
        agg_dir = tmp_path / "aggregated"

        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=agg_dir,
            equilibration="10ns",
            settings=default_settings,
        )

        with patch("polyzymd.analysis.results.base.get_polyzymd_version", return_value="1.2.1"):
            result = triad_analysis.aggregate(ctx, results)

        assert result.n_replicates == 3
        assert result.replicates == [1, 2, 3]
        assert result.triad_name == "LipA_triad"
        assert len(result.pair_results) == 2
        assert len(result.per_replicate_simultaneous) == 3
        # Mean of 0.65, 0.70, 0.75
        assert result.overall_simultaneous_contact == pytest.approx(0.7, abs=0.01)
        assert result.sem_simultaneous_contact > 0

    def test_aggregate_saves_result_file(
        self, triad_analysis, condition, tmp_path, default_settings
    ):
        """Verify aggregated result is saved to the output directory."""
        results = [_make_mock_triad_result(i, sim_contact=0.65) for i in range(1, 3)]
        agg_dir = tmp_path / "aggregated"

        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=agg_dir,
            equilibration="10ns",
            settings=default_settings,
        )

        with patch("polyzymd.analysis.results.base.get_polyzymd_version", return_value="1.2.1"):
            triad_analysis.aggregate(ctx, results)

        json_files = list(agg_dir.glob("*.json"))
        assert len(json_files) == 1
        assert "triad_LipA_triad_reps1-2_eq10ns.json" in json_files[0].name

    def test_aggregate_preserves_per_replicate_contact(
        self, triad_analysis, condition, tmp_path, default_settings
    ):
        """Verify per_replicate_simultaneous has the correct values."""
        results = [_make_mock_triad_result(i, sim_contact=0.5 + 0.1 * i) for i in range(1, 4)]

        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2, 3),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=default_settings,
        )

        with patch("polyzymd.analysis.results.base.get_polyzymd_version", return_value="1.2.1"):
            result = triad_analysis.aggregate(ctx, results)

        assert result.per_replicate_simultaneous == pytest.approx([0.6, 0.7, 0.8])


# ============================================================================
# Test: extract_metrics
# ============================================================================


class TestExtractMetrics:
    """Test CatalyticTriadAnalysis.extract_metrics for default comparison pipeline."""

    def test_returns_simultaneous_contact_metric(self, triad_analysis):
        summary = MagicMock()
        summary.overall_simultaneous_contact = 0.72
        summary.sem_simultaneous_contact = 0.04
        summary.per_replicate_simultaneous = [0.70, 0.72, 0.74]

        metrics = triad_analysis.extract_metrics(summary)

        assert "simultaneous_contact_fraction" in metrics
        m = metrics["simultaneous_contact_fraction"]
        assert isinstance(m, MetricValue)
        assert m.name == "simultaneous_contact_fraction"
        assert m.mean == 0.72
        assert m.sem == 0.04
        assert m.replicate_values == [0.70, 0.72, 0.74]
        assert m.higher_is_better is True
        assert m.direction_labels == ("worsening", "unchanged", "improving")

    def test_returns_single_metric(self, triad_analysis):
        summary = MagicMock()
        summary.overall_simultaneous_contact = 0.5
        summary.sem_simultaneous_contact = 0.1
        summary.per_replicate_simultaneous = [0.4, 0.6]

        metrics = triad_analysis.extract_metrics(summary)
        assert len(metrics) == 1


# ============================================================================
# Test: _deserialize_result
# ============================================================================


class TestDeserializeResult:
    """Test CatalyticTriadAnalysis._deserialize_result."""

    def test_loads_aggregated_result(self, triad_analysis, tmp_path):
        with patch("polyzymd.analysis.results.triad.TriadAggregatedResult") as MockAgg:
            mock_loaded = MagicMock()
            MockAgg.load.return_value = mock_loaded

            result = triad_analysis._deserialize_result(tmp_path / "test.json")

            MockAgg.load.assert_called_once_with(tmp_path / "test.json")
            assert result is mock_loaded


# ============================================================================
# Test: filter_conditions (default behavior — keeps all)
# ============================================================================


class TestFilterConditions:
    """Catalytic triad applies to all conditions — filter should keep all."""

    def test_keeps_all_conditions(self, triad_analysis):
        conditions = [
            Condition(
                label=f"Cond {i}",
                config_path=Path(f"/fake/config_{i}.yaml"),
                replicates=(1, 2),
                sim_config=MagicMock(),
            )
            for i in range(4)
        ]

        filtered = triad_analysis.filter_conditions(conditions)
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
        result.triad_name = "LipA_triad"

        filename = CatalyticTriadAnalysis._make_aggregated_filename((1, 2, 3, 4, 5), result)
        assert filename == "triad_LipA_triad_reps1-5_eq10ns.json"

    def test_non_contiguous_replicates(self):
        result = MagicMock()
        result.equilibration_time = 100.0
        result.equilibration_unit = "ns"
        result.triad_name = "LipA_triad"

        filename = CatalyticTriadAnalysis._make_aggregated_filename((1, 3, 5), result)
        assert filename == "triad_LipA_triad_reps1_3_5_eq100ns.json"

    def test_single_pair(self):
        result = MagicMock()
        result.equilibration_time = 0.0
        result.equilibration_unit = "ns"
        result.triad_name = "LipA_triad"

        filename = CatalyticTriadAnalysis._make_aggregated_filename((1, 2), result)
        assert filename == "triad_LipA_triad_reps1-2_eq0ns.json"

    def test_name_sanitization(self):
        result = MagicMock()
        result.equilibration_time = 10.0
        result.equilibration_unit = "ns"
        result.triad_name = "Ser His/Asp triad"

        filename = CatalyticTriadAnalysis._make_aggregated_filename((1, 2, 3), result)
        assert filename == "triad_Ser_His-Asp_triad_reps1-3_eq10ns.json"


# ============================================================================
# Test: plot
# ============================================================================


class TestPlot:
    """Test CatalyticTriadAnalysis.plot delegates to existing plotters."""

    def test_plot_returns_empty_on_no_conditions(self, triad_analysis, tmp_path, default_settings):
        ctx = PlotContext(
            conditions=[],
            analysis_dirs={},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=default_settings,
            plot_settings=None,
        )

        plots = triad_analysis.plot(ctx)
        assert plots == []

    @patch("polyzymd.compare.plotters.triad.TriadThresholdBarsPlotter")
    @patch("polyzymd.compare.plotters.triad.TriadKDEPanelPlotter")
    def test_plot_delegates_to_plotters(
        self,
        MockKDEPlotter,
        MockBarsPlotter,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        # Setup mock plotters
        mock_kde = MagicMock()
        mock_kde.plot.return_value = [tmp_path / "figures" / "triad_kde_panel.png"]
        MockKDEPlotter.return_value = mock_kde

        mock_bars = MagicMock()
        mock_bars.plot.return_value = [tmp_path / "figures" / "triad_threshold_bars.png"]
        MockBarsPlotter.return_value = mock_bars

        analysis_dir = tmp_path / "analysis" / "no_polymer" / "catalytic_triad"
        analysis_dir.mkdir(parents=True)

        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"No Polymer": analysis_dir},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=default_settings,
            plot_settings=None,
        )

        with patch("polyzymd.compare.config.PlotSettings") as MockPlotSettings:
            MockPlotSettings.return_value = MagicMock()
            plots = triad_analysis.plot(ctx)

        assert len(plots) == 2
        mock_kde.plot.assert_called_once()
        mock_bars.plot.assert_called_once()

    @patch(
        "polyzymd.compare.plotters.triad.TriadThresholdBarsPlotter",
        side_effect=Exception("plot error"),
    )
    @patch(
        "polyzymd.compare.plotters.triad.TriadKDEPanelPlotter", side_effect=Exception("plot error")
    )
    def test_plot_catches_exceptions(
        self,
        MockKDEPlotter,
        MockBarsPlotter,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        """Plotter failures should be caught, not crash the pipeline."""
        analysis_dir = tmp_path / "analysis" / "no_polymer" / "catalytic_triad"
        analysis_dir.mkdir(parents=True)

        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"No Polymer": analysis_dir},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=default_settings,
            plot_settings=None,
        )

        with patch("polyzymd.compare.config.PlotSettings") as MockPlotSettings:
            MockPlotSettings.return_value = MagicMock()
            plots = triad_analysis.plot(ctx)

        assert plots == []

    def test_plot_passes_data_and_labels(self, triad_analysis, tmp_path, default_settings):
        """Verify the data dict passed to plotters has the expected shape."""
        cond1 = Condition(
            label="No Polymer",
            config_path=Path("/fake/config1.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        cond2 = Condition(
            label="100% SBMA",
            config_path=Path("/fake/config2.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )

        for label in ("no_polymer", "sbma"):
            (tmp_path / "analysis" / label / "catalytic_triad").mkdir(parents=True)

        captured_data = {}
        captured_labels = []

        with patch("polyzymd.compare.plotters.triad.TriadKDEPanelPlotter") as MockKDE:
            mock_kde_inst = MagicMock()
            mock_kde_inst.plot.return_value = []
            MockKDE.return_value = mock_kde_inst

            # Capture what's passed
            def capture_plot(data, labels, output_dir, **kwargs):
                captured_data.update(data)
                captured_labels.extend(labels)
                return []

            mock_kde_inst.plot.side_effect = capture_plot

            with patch("polyzymd.compare.plotters.triad.TriadThresholdBarsPlotter") as MockBars:
                MockBars.return_value.plot.return_value = []

                ctx = PlotContext(
                    conditions=[cond1, cond2],
                    analysis_dirs={
                        "No Polymer": tmp_path / "analysis" / "no_polymer" / "catalytic_triad",
                        "100% SBMA": tmp_path / "analysis" / "sbma" / "catalytic_triad",
                    },
                    results_dir=tmp_path / "comparison",
                    output_dir=tmp_path / "figures",
                    settings=default_settings,
                    plot_settings=MagicMock(),  # Already provided
                )

                triad_analysis.plot(ctx)

        assert "No Polymer" in captured_data
        assert "100% SBMA" in captured_data
        assert captured_labels == ["No Polymer", "100% SBMA"]
        # Each condition entry should have analysis_dir, aggregated_dir, replicates
        for label in ("No Polymer", "100% SBMA"):
            entry = captured_data[label]
            assert "analysis_dir" in entry
            assert "aggregated_dir" in entry
            assert "replicates" in entry


# ============================================================================
# Test: Full lifecycle via orchestrator (mocked compute)
# ============================================================================


class TestTriadLifecycle:
    """Test the full compute -> aggregate -> compare lifecycle."""

    @patch("polyzymd.analysis.triad.analyzer.CatalyticTriadAnalyzer")
    def test_run_analysis_lifecycle(
        self, MockAnalyzer, triad_analysis, condition, tmp_path, default_settings
    ):
        """Test compute_replicate -> aggregate via run_analysis()."""
        from polyzymd.analyses.orchestrator import run_analysis

        mock_results = {
            rep: _make_mock_triad_result(rep, sim_contact=0.6 + 0.05 * rep) for rep in (1, 2, 3)
        }

        def mock_compute(replicate, save=True, output_dir=None, recompute=False):
            return mock_results[replicate]

        MockAnalyzer.return_value.compute = mock_compute

        with patch("polyzymd.analysis.results.base.get_polyzymd_version", return_value="1.2.1"):
            output_dir = tmp_path / "analysis" / "no_polymer" / "catalytic_triad"
            result = run_analysis(
                triad_analysis,
                condition,
                settings=default_settings,
                equilibration="10ns",
                output_dir=output_dir,
            )

        assert result.n_replicates == 3
        assert result.per_replicate_simultaneous == pytest.approx([0.65, 0.70, 0.75])
        assert MockAnalyzer.call_count == 3

    def test_extract_metrics_feeds_default_compare(self, triad_analysis):
        """Verify extract_metrics output is compatible with default_scalar_comparison."""
        from polyzymd.analyses.stats import default_scalar_comparison

        metrics_by_condition = {}
        for label, mean, sem, vals in [
            ("No Polymer", 0.55, 0.04, [0.50, 0.55, 0.60]),
            ("100% SBMA", 0.75, 0.03, [0.72, 0.75, 0.78]),
        ]:
            summary = MagicMock()
            summary.overall_simultaneous_contact = mean
            summary.sem_simultaneous_contact = sem
            summary.per_replicate_simultaneous = vals
            metrics_by_condition[label] = triad_analysis.extract_metrics(summary)

        result = default_scalar_comparison(
            analysis_name="catalytic_triad",
            project_name="test_project",
            metrics_by_condition=metrics_by_condition,
            control_label="No Polymer",
            equilibration="10ns",
        )

        assert result.analysis_type == "catalytic_triad"
        assert len(result.conditions) == 2
        assert len(result.pairwise_comparisons) >= 1
        # Higher triad score = better, so 100% SBMA should rank first
        assert result.ranking[0] == "100% SBMA"
