"""Tests for the catalytic triad analysis plugin.

Tests the CatalyticTriadAnalysis class: discovery, aliases, settings,
compute_replicate, aggregate, extract_metrics, AggregatedResultClass,
_settings_to_config, _make_aggregated_filename, plot delegation, and the
full lifecycle via the orchestrator.

Heavy dependencies (MDAnalysis, trajectories) are mocked.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from polyzymd.analyses.base import (
    AggregateContext,
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
# Test: compute_replicate
# ============================================================================


def _make_mock_distance_result(n_pairs: int = 2, n_frames: int = 1000):
    """Create a mock DistanceResult with *real* DistancePairResult objects.

    Uses genuine ``DistancePairResult`` instances so that ``model_copy()``
    works correctly and ``TriadResult`` Pydantic validation passes.
    """
    import numpy as np

    from polyzymd.analyses.distances._results import DistancePairResult

    mock_result = MagicMock()
    mock_result.n_frames_total = n_frames + 200
    mock_result.n_frames_used = n_frames

    pair_results = []
    for i in range(n_pairs):
        rng = np.random.default_rng(42 + i)
        dists = rng.normal(3.2, 0.8, n_frames).tolist()

        pr = DistancePairResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string=f"resid {100 + i} and name OD1 : resid {200 + i} and name ND1",
            pair_label=f"orig_pair_{i}",
            selection1=f"resid {100 + i} and name OD1",
            selection2=f"resid {200 + i} and name ND1",
            distances=dists,
            mean_distance=float(np.mean(dists)),
            std_distance=float(np.std(dists)),
            median_distance=float(np.median(dists)),
            min_distance=float(np.min(dists)),
            max_distance=float(np.max(dists)),
            sem_distance=0.1,
            threshold=3.5,
            fraction_below_threshold=float(np.mean(np.array(dists) < 3.5)),
            kde_peak=2.8 + 0.1 * i,
            n_frames_total=n_frames + 200,
            n_frames_used=n_frames,
        )
        pair_results.append(pr)

    mock_result.pair_results = pair_results
    return mock_result


class TestComputeReplicate:
    """Test CatalyticTriadAnalysis.compute_replicate with inlined computation."""

    @patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1")
    @patch("polyzymd.analyses.shared.config_hash.validate_config_hash")
    @patch("polyzymd.analyses.shared.config_hash.compute_config_hash", return_value="hash123")
    @patch("polyzymd.analyses.shared.TrajectoryLoader")
    @patch("polyzymd.analyses.distances.DistanceCalculator")
    def test_computes_triad_inline(
        self,
        MockDistCalc,
        MockLoader,
        mock_hash,
        mock_validate_hash,
        mock_version,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        """compute_replicate should use DistanceCalculator and compute simultaneous contact."""

        # Mock TrajectoryLoader
        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst
        mock_u = MagicMock()
        # Make select_atoms return non-empty for validation
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=1)
        mock_u.select_atoms.return_value = mock_atoms
        mock_loader_inst.load_universe.return_value = mock_u
        mock_loader_inst.get_timestep.return_value = 10.0

        # Mock DistanceCalculator
        mock_dist_result = _make_mock_distance_result(n_pairs=2, n_frames=1000)
        mock_dist_inst = MagicMock()
        mock_dist_inst.compute.return_value = mock_dist_result
        MockDistCalc.return_value = mock_dist_inst

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,  # skip cache check
            settings=default_settings,
        )

        result = triad_analysis.compute_replicate(ctx, 1)

        # Verify DistanceCalculator was created with correct args
        MockDistCalc.assert_called_once()
        call_kwargs = MockDistCalc.call_args
        assert call_kwargs.kwargs["config"] is condition.sim_config
        assert len(call_kwargs.kwargs["pairs"]) == 2

        # Verify result has expected triad fields
        assert result.replicate == 1
        assert result.triad_name == "LipA_triad"
        assert result.threshold == 3.5
        assert len(result.pair_results) == 2
        # Pair labels should be updated from settings
        assert result.pair_results[0].pair_label == "Asp133-His156"
        assert result.pair_results[1].pair_label == "His156-Ser77"
        # Simultaneous contact fraction should be computed
        assert 0.0 <= result.simultaneous_contact_fraction <= 1.0
        assert result.n_frames_used == 1000

    @patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1")
    @patch("polyzymd.analyses.shared.config_hash.validate_config_hash")
    @patch("polyzymd.analyses.shared.config_hash.compute_config_hash", return_value="hash123")
    @patch("polyzymd.analyses.shared.TrajectoryLoader")
    @patch("polyzymd.analyses.distances.DistanceCalculator")
    def test_passes_recompute_flag(
        self,
        MockDistCalc,
        MockLoader,
        mock_hash,
        mock_validate_hash,
        mock_version,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        """recompute flag should be passed to DistanceCalculator.compute()."""
        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst
        mock_u = MagicMock()
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=1)
        mock_u.select_atoms.return_value = mock_atoms
        mock_loader_inst.load_universe.return_value = mock_u
        mock_loader_inst.get_timestep.return_value = 10.0

        mock_dist_result = _make_mock_distance_result(n_pairs=2, n_frames=500)
        mock_dist_inst = MagicMock()
        mock_dist_inst.compute.return_value = mock_dist_result
        MockDistCalc.return_value = mock_dist_inst

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

        mock_dist_inst.compute.assert_called_once_with(
            replicate=2,
            save=False,
            recompute=True,
            store_distributions=True,
        )

    @patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1")
    @patch("polyzymd.analyses.shared.config_hash.validate_config_hash")
    @patch("polyzymd.analyses.shared.config_hash.compute_config_hash", return_value="hash123")
    @patch("polyzymd.analyses.shared.TrajectoryLoader")
    @patch("polyzymd.analyses.distances.DistanceCalculator")
    def test_caches_result_file(
        self,
        MockDistCalc,
        MockLoader,
        mock_hash,
        mock_validate_hash,
        mock_version,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        """Result should be saved to the output directory."""
        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst
        mock_u = MagicMock()
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=1)
        mock_u.select_atoms.return_value = mock_atoms
        mock_loader_inst.load_universe.return_value = mock_u
        mock_loader_inst.get_timestep.return_value = 10.0

        mock_dist_result = _make_mock_distance_result(n_pairs=2, n_frames=100)
        mock_dist_inst = MagicMock()
        mock_dist_inst.compute.return_value = mock_dist_result
        MockDistCalc.return_value = mock_dist_inst

        output_dir = tmp_path / "run_1"
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=output_dir,
            equilibration="10ns",
            recompute=True,
            settings=default_settings,
        )

        triad_analysis.compute_replicate(ctx, 1)

        # Check that result JSON was written
        json_files = list(output_dir.glob("*.json"))
        assert len(json_files) == 1
        assert "triad_LipA_triad_eq10.00ns.json" in json_files[0].name


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

        with patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1"):
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

        with patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1"):
            triad_analysis.aggregate(ctx, results)

        json_files = list(agg_dir.glob("*.json"))
        assert len(json_files) == 1
        assert "triad_LipA_triad_reps1-2_eq10.00ns.json" in json_files[0].name

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

        with patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1"):
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
        assert m.mean == 72.0
        assert m.sem == 4.0
        assert m.replicate_values == [70.0, 72.0, 74.0]
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
# Test: AggregatedResultClass and _deserialize_result
# ============================================================================


class TestDeserializeResult:
    """Test CatalyticTriadAnalysis.AggregatedResultClass and _deserialize_result."""

    def test_aggregated_result_class_set(self, triad_analysis):
        """AggregatedResultClass should be TriadAggregatedResult."""
        from polyzymd.analyses.catalytic_triad._results import TriadAggregatedResult

        assert triad_analysis.AggregatedResultClass is TriadAggregatedResult

    def test_loads_aggregated_result(self, triad_analysis, tmp_path):
        mock_loaded = MagicMock()
        with patch.object(
            triad_analysis.AggregatedResultClass, "load", return_value=mock_loaded
        ) as mock_load:
            result = triad_analysis._deserialize_result(tmp_path / "test.json")

            mock_load.assert_called_once_with(tmp_path / "test.json")
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
        assert filename == "triad_LipA_triad_reps1-5_eq10.00ns.json"

    def test_non_contiguous_replicates(self):
        result = MagicMock()
        result.equilibration_time = 100.0
        result.equilibration_unit = "ns"
        result.triad_name = "LipA_triad"

        filename = CatalyticTriadAnalysis._make_aggregated_filename((1, 3, 5), result)
        assert filename == "triad_LipA_triad_reps1_3_5_eq100.00ns.json"

    def test_single_pair(self):
        result = MagicMock()
        result.equilibration_time = 0.0
        result.equilibration_unit = "ns"
        result.triad_name = "LipA_triad"

        filename = CatalyticTriadAnalysis._make_aggregated_filename((1, 2), result)
        assert filename == "triad_LipA_triad_reps1-2_eq0.00ns.json"

    def test_name_sanitization(self):
        result = MagicMock()
        result.equilibration_time = 10.0
        result.equilibration_unit = "ns"
        result.triad_name = "Ser His/Asp triad"

        filename = CatalyticTriadAnalysis._make_aggregated_filename((1, 2, 3), result)
        assert filename == "triad_Ser_His-Asp_triad_reps1-3_eq10.00ns.json"


# ============================================================================
# Test: plot
# ============================================================================


class TestPlot:
    """Test CatalyticTriadAnalysis.plot delegates to inlined plotting helpers."""

    def test_plot_returns_empty_on_no_conditions(self, triad_analysis, tmp_path, default_settings):
        from polyzymd.config.comparison import PlotSettings

        ctx = PlotContext(
            conditions=[],
            analysis_dirs={},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=default_settings,
            plot_settings=PlotSettings(),
        )

        plots = triad_analysis.plot(ctx)
        assert plots == []

    @patch("polyzymd.analyses.catalytic_triad.plot_triad_threshold_bars_from_data")
    @patch("polyzymd.analyses.catalytic_triad.plot_triad_kde_panel_from_data")
    def test_plot_delegates_to_helpers(
        self,
        mock_kde_fn,
        mock_bars_fn,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        from polyzymd.config.comparison import PlotSettings

        # Setup mock return values
        mock_kde_fn.return_value = [tmp_path / "figures" / "triad_kde_panel.png"]
        mock_bars_fn.return_value = [tmp_path / "figures" / "triad_threshold_bars.png"]

        analysis_dir = tmp_path / "analysis" / "no_polymer" / "catalytic_triad"
        analysis_dir.mkdir(parents=True)

        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"No Polymer": analysis_dir},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=default_settings,
            plot_settings=PlotSettings(),
        )

        with patch("polyzymd.config.comparison.PlotSettings") as MockPlotSettings:
            MockPlotSettings.return_value = MagicMock()
            plots = triad_analysis.plot(ctx)

        assert len(plots) == 2
        mock_kde_fn.assert_called_once()
        mock_bars_fn.assert_called_once()

        # Verify PlotSettings is passed to both functions
        kde_args = mock_kde_fn.call_args[0]
        bars_args = mock_bars_fn.call_args[0]
        assert len(kde_args) == 4  # data, labels, output_dir, plot_settings
        assert len(bars_args) == 4

    @patch(
        "polyzymd.analyses.catalytic_triad.plot_triad_threshold_bars_from_data",
        side_effect=Exception("plot error"),
    )
    @patch(
        "polyzymd.analyses.catalytic_triad.plot_triad_kde_panel_from_data",
        side_effect=Exception("plot error"),
    )
    def test_plot_propagates_exceptions(
        self,
        mock_kde_fn,
        mock_bars_fn,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        """Plotter failures should propagate to orchestrator."""
        from polyzymd.config.comparison import PlotSettings

        analysis_dir = tmp_path / "analysis" / "no_polymer" / "catalytic_triad"
        analysis_dir.mkdir(parents=True)

        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={"No Polymer": analysis_dir},
            results_dir=tmp_path / "comparison",
            output_dir=tmp_path / "figures",
            settings=default_settings,
            plot_settings=PlotSettings(),
        )

        with patch("polyzymd.config.comparison.PlotSettings") as MockPlotSettings:
            MockPlotSettings.return_value = MagicMock()
            with pytest.raises(Exception, match="plot error"):
                triad_analysis.plot(ctx)

    def test_plot_passes_data_and_labels(self, triad_analysis, tmp_path, default_settings):
        """Verify the data dict passed to helpers has the expected shape."""
        from polyzymd.config.comparison import PlotSettings

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

        def capture_kde(data, labels, output_dir, plot_settings):
            captured_data.update(data)
            captured_labels.extend(labels)
            return []

        with patch(
            "polyzymd.analyses.catalytic_triad.plot_triad_kde_panel_from_data",
            side_effect=capture_kde,
        ) as mock_kde_fn:
            with patch(
                "polyzymd.analyses.catalytic_triad.plot_triad_threshold_bars_from_data"
            ) as mock_bars_fn:
                mock_bars_fn.return_value = []

                ctx = PlotContext(
                    conditions=[cond1, cond2],
                    analysis_dirs={
                        "No Polymer": tmp_path / "analysis" / "no_polymer" / "catalytic_triad",
                        "100% SBMA": tmp_path / "analysis" / "sbma" / "catalytic_triad",
                    },
                    results_dir=tmp_path / "comparison",
                    output_dir=tmp_path / "figures",
                    settings=default_settings,
                    plot_settings=PlotSettings(),
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

    @patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1")
    @patch("polyzymd.analyses.shared.config_hash.validate_config_hash")
    @patch("polyzymd.analyses.shared.config_hash.compute_config_hash", return_value="hash123")
    @patch("polyzymd.analyses.shared.TrajectoryLoader")
    @patch("polyzymd.analyses.distances.DistanceCalculator")
    def test_run_analysis_lifecycle(
        self,
        MockDistCalc,
        MockLoader,
        mock_hash,
        mock_validate_hash,
        mock_version,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        """Test compute_replicate -> aggregate via run_analysis()."""
        from polyzymd.analyses.orchestrator import run_analysis

        # Mock TrajectoryLoader (shared by all replicates)
        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst
        mock_u = MagicMock()
        mock_atoms = MagicMock()
        mock_atoms.__len__ = MagicMock(return_value=1)
        mock_u.select_atoms.return_value = mock_atoms
        mock_loader_inst.load_universe.return_value = mock_u
        mock_loader_inst.get_timestep.return_value = 10.0

        def make_dist_result_for_replicate(*args, **kwargs):
            return _make_mock_distance_result(n_pairs=2, n_frames=100)

        mock_dist_inst = MagicMock()
        mock_dist_inst.compute.side_effect = make_dist_result_for_replicate
        MockDistCalc.return_value = mock_dist_inst

        output_dir = tmp_path / "analysis" / "no_polymer" / "catalytic_triad"
        result = run_analysis(
            triad_analysis,
            condition,
            settings=default_settings,
            equilibration="10ns",
            output_dir=output_dir,
        )

        # Verify orchestrator called compute for each replicate
        assert mock_dist_inst.compute.call_count == 3

        # Verify aggregate produced a valid result
        assert result.n_replicates == 3
        assert result.replicates == [1, 2, 3]
        assert result.triad_name == "LipA_triad"
        assert len(result.pair_results) == 2
        assert len(result.per_replicate_simultaneous) == 3
        # All replicates used the same random seed → same fraction, so all should be equal
        for val in result.per_replicate_simultaneous:
            assert 0.0 <= val <= 1.0
        assert result.overall_simultaneous_contact == pytest.approx(
            sum(result.per_replicate_simultaneous) / 3, abs=0.01
        )

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
