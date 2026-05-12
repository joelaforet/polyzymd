"""Tests for the catalytic triad analysis plugin.

Tests the CatalyticTriadAnalysis class: discovery, aliases, settings,
run_replicate, aggregate, extract_metrics, AggregatedResultClass,
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
    AggregateValidationError,
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
from polyzymd.analyses.catalytic_triad._results import TriadResult
from polyzymd.analyses.distances._results import DistancePairResult
from polyzymd.analyses.shared.config_hash import settings_fingerprint

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


def _make_mock_pair_result(pair_idx: int, replicate: int) -> DistancePairResult:
    """Create a valid DistancePairResult for one pair in one replicate."""
    pair_schemas = [
        ("Asp133-His156", "resid 133 and name OD1", "resid 156 and name ND1"),
        ("His156-Ser77", "resid 156 and name NE2", "resid 77 and name OG"),
    ]
    pair_label, selection1, selection2 = pair_schemas[pair_idx]
    return DistancePairResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=replicate,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string=f"{selection1} : {selection2}",
        pair_label=pair_label,
        selection1=selection1,
        selection2=selection2,
        mean_distance=3.0 + 0.1 * pair_idx,
        std_distance=0.5,
        sem_distance=0.1,
        median_distance=3.0 + 0.1 * pair_idx,
        min_distance=2.0,
        max_distance=5.0,
        threshold=3.5,
        fraction_below_threshold=0.7 - 0.05 * pair_idx,
        kde_peak=2.8 + 0.1 * pair_idx,
        distances=[3.0, 3.2, 2.8, 3.5, 2.9],
        n_frames_total=1200,
        n_frames_used=1000,
    )


def _make_mock_triad_result(
    replicate: int,
    sim_contact: float = 0.65,
    *,
    settings_fp: str | None = None,
) -> TriadResult:
    """Create a valid TriadResult with realistic fields."""
    return TriadResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=replicate,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="(resid 133 : resid 156); (resid 156 : resid 77)",
        triad_name="LipA_triad",
        triad_description="Ser-His-Asp catalytic triad",
        pair_results=[
            _make_mock_pair_result(0, replicate),
            _make_mock_pair_result(1, replicate),
        ],
        threshold=3.5,
        simultaneous_contact_fraction=sim_contact,
        n_frames_simultaneous=int(sim_contact * 1000),
        n_frames_total=1200,
        n_frames_used=1000,
        settings_fingerprint=settings_fp,
    )


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
        assert triad_analysis.min_replicates == 1

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
# Test: run_replicate
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


class TestRunReplicate:
    """Test catalytic-triad runner seam behavior."""

    def test_run_replicate_delegates_to_runner_seam(
        self,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        from polyzymd.analyses.base import Analysis as AnalysisBase

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,
            settings=default_settings,
        )
        mock_result = MagicMock()

        with patch.object(
            AnalysisBase, "_run_replicate_default", return_value=mock_result
        ) as mock_super:
            result = triad_analysis.run_replicate(ctx, 1)

        assert result is mock_result
        mock_super.assert_called_once_with(ctx, 1)
        mock_result.save.assert_called_once()
        assert mock_result.save.call_args[0][0].name == (
            f"triad_LipA_triad_eq10ns_s{settings_fingerprint(default_settings)}.json"
        )

    def test_run_replicate_cache_filename_changes_when_settings_change(
        self,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        """Replicate cache filenames should change when triad settings change."""

        seen_paths: list[Path] = []

        def _fake_check_cache(result_cls, cache_path, **kwargs):  # noqa: ARG001
            seen_paths.append(cache_path)
            return {"cached": True}

        threshold_settings = default_settings.model_copy(update={"threshold": 4.0})
        triad_analysis._check_cache = _fake_check_cache

        ctx_a = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1a",
            equilibration="10ns",
            recompute=False,
            settings=default_settings,
        )
        ctx_b = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1b",
            equilibration="10ns",
            recompute=False,
            settings=threshold_settings,
        )

        assert triad_analysis.run_replicate(ctx_a, 1) == {"cached": True}
        assert triad_analysis.run_replicate(ctx_b, 1) == {"cached": True}
        assert len(seen_paths) == 2
        assert seen_paths[0].name != seen_paths[1].name

    def test_build_runner_returns_triad_runner(self, triad_analysis, condition, default_settings):
        from polyzymd.analyses.catalytic_triad._runner import CatalyticTriadReplicateRunner

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=Path("/tmp/run_1"),
            equilibration="10ns",
            recompute=True,
            settings=default_settings,
        )

        runner = triad_analysis.build_runner(ctx, 1, MagicMock(), MagicMock(timestep_ps=10.0))

        assert isinstance(runner, CatalyticTriadReplicateRunner)

    @patch("polyzymd.analyses.shared.config_hash.compute_config_hash", return_value="hash123")
    @patch("polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.2.1")
    def test_summarize_replicate_computes_simultaneous_contact(
        self,
        mock_version,
        mock_hash,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        import numpy as np

        from polyzymd.analyses.distances._runner import DistancePairPayload, DistancesRunnerPayload

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,
            settings=default_settings,
        )
        runner = MagicMock(
            results=DistancesRunnerPayload(
                n_frames_total=12,
                n_frames_used=4,
                pair_payloads=[
                    DistancePairPayload(
                        pair_label="auto0",
                        selection1="resid 133 and name OD1",
                        selection2="resid 156 and name ND1",
                        distances=np.asarray([3.0, 3.8, 3.1, 3.2], dtype=np.float64),
                        mean_distance=3.275,
                        std_distance=0.0,
                        median_distance=3.15,
                        min_distance=3.0,
                        max_distance=3.8,
                        sem_distance=None,
                        correlation_time=None,
                        correlation_time_unit=None,
                        n_independent_frames=None,
                        statistical_inefficiency=None,
                        autocorrelation_warning=None,
                        threshold=3.5,
                        fraction_below_threshold=0.75,
                        histogram_edges=np.asarray([3.0, 3.4, 3.8], dtype=np.float64),
                        histogram_counts=np.asarray([3, 1], dtype=np.int64),
                        kde_x=None,
                        kde_y=None,
                        kde_peak=None,
                        kde_bandwidth=None,
                        n_frames_total=12,
                        n_frames_used=4,
                    ),
                    DistancePairPayload(
                        pair_label="auto1",
                        selection1="resid 156 and name NE2",
                        selection2="resid 77 and name OG",
                        distances=np.asarray([3.1, 3.0, 3.7, 3.2], dtype=np.float64),
                        mean_distance=3.25,
                        std_distance=0.0,
                        median_distance=3.15,
                        min_distance=3.0,
                        max_distance=3.7,
                        sem_distance=None,
                        correlation_time=None,
                        correlation_time_unit=None,
                        n_independent_frames=None,
                        statistical_inefficiency=None,
                        autocorrelation_warning=None,
                        threshold=3.5,
                        fraction_below_threshold=0.75,
                        histogram_edges=np.asarray([3.0, 3.35, 3.7], dtype=np.float64),
                        histogram_counts=np.asarray([3, 1], dtype=np.int64),
                        kde_x=None,
                        kde_y=None,
                        kde_peak=None,
                        kde_bandwidth=None,
                        n_frames_total=12,
                        n_frames_used=4,
                    ),
                ],
            )
        )
        window = MagicMock(timestep_ps=10.0, step=1)

        result = triad_analysis.summarize_replicate(ctx, 1, runner, window)

        assert result.replicate == 1
        assert result.triad_name == "LipA_triad"
        assert result.pair_results[0].pair_label == "Asp133-His156"
        assert result.pair_results[1].pair_label == "His156-Ser77"
        assert result.pair_results[0].replicate == 1
        assert result.pair_results[1].replicate == 1
        assert result.simultaneous_contact_fraction == pytest.approx(0.5)
        assert result.n_frames_simultaneous == 2
        assert result.settings_fingerprint == settings_fingerprint(default_settings)

    def test_summarize_replicate_falls_back_on_expected_autocorrelation_failure(
        self,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        """Expected autocorrelation failures should use the naive SEM fallback."""
        import numpy as np

        from polyzymd.analyses.distances._runner import DistancePairPayload, DistancesRunnerPayload

        distances = np.full(20, 3.0, dtype=np.float64)
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,
            settings=default_settings,
        )
        runner = MagicMock(
            results=DistancesRunnerPayload(
                n_frames_total=20,
                n_frames_used=20,
                pair_payloads=[
                    DistancePairPayload(
                        pair_label="auto0",
                        selection1="resid 133 and name OD1",
                        selection2="resid 156 and name ND1",
                        distances=distances,
                        mean_distance=3.0,
                        std_distance=0.0,
                        median_distance=3.0,
                        min_distance=3.0,
                        max_distance=3.0,
                        sem_distance=None,
                        correlation_time=None,
                        correlation_time_unit=None,
                        n_independent_frames=None,
                        statistical_inefficiency=None,
                        autocorrelation_warning=None,
                        threshold=3.5,
                        fraction_below_threshold=1.0,
                        histogram_edges=np.asarray([2.5, 3.0, 3.5], dtype=np.float64),
                        histogram_counts=np.asarray([0, 20], dtype=np.int64),
                        kde_x=None,
                        kde_y=None,
                        kde_peak=None,
                        kde_bandwidth=None,
                        n_frames_total=20,
                        n_frames_used=20,
                    ),
                    DistancePairPayload(
                        pair_label="auto1",
                        selection1="resid 156 and name NE2",
                        selection2="resid 77 and name OG",
                        distances=distances,
                        mean_distance=3.0,
                        std_distance=0.0,
                        median_distance=3.0,
                        min_distance=3.0,
                        max_distance=3.0,
                        sem_distance=None,
                        correlation_time=None,
                        correlation_time_unit=None,
                        n_independent_frames=None,
                        statistical_inefficiency=None,
                        autocorrelation_warning=None,
                        threshold=3.5,
                        fraction_below_threshold=1.0,
                        histogram_edges=np.asarray([2.5, 3.0, 3.5], dtype=np.float64),
                        histogram_counts=np.asarray([0, 20], dtype=np.int64),
                        kde_x=None,
                        kde_y=None,
                        kde_peak=None,
                        kde_bandwidth=None,
                        n_frames_total=20,
                        n_frames_used=20,
                    ),
                ],
            )
        )
        window = MagicMock(timestep_ps=10.0, step=1)

        with patch(
            "polyzymd.analyses.shared.autocorrelation.estimate_correlation_time",
            side_effect=ValueError("autocorrelation failed"),
        ):
            result = triad_analysis.summarize_replicate(ctx, 1, runner, window)

        assert result.sim_contact_sem == pytest.approx(0.0)

    def test_summarize_replicate_surfaces_unexpected_autocorrelation_errors(
        self,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        """Unexpected autocorrelation errors should not be masked."""
        import numpy as np

        from polyzymd.analyses.distances._runner import DistancePairPayload, DistancesRunnerPayload

        distances = np.full(20, 3.0, dtype=np.float64)
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,
            settings=default_settings,
        )
        runner = MagicMock(
            results=DistancesRunnerPayload(
                n_frames_total=20,
                n_frames_used=20,
                pair_payloads=[
                    DistancePairPayload(
                        pair_label="auto0",
                        selection1="resid 133 and name OD1",
                        selection2="resid 156 and name ND1",
                        distances=distances,
                        mean_distance=3.0,
                        std_distance=0.0,
                        median_distance=3.0,
                        min_distance=3.0,
                        max_distance=3.0,
                        sem_distance=None,
                        correlation_time=None,
                        correlation_time_unit=None,
                        n_independent_frames=None,
                        statistical_inefficiency=None,
                        autocorrelation_warning=None,
                        threshold=3.5,
                        fraction_below_threshold=1.0,
                        histogram_edges=np.asarray([2.5, 3.0, 3.5], dtype=np.float64),
                        histogram_counts=np.asarray([0, 20], dtype=np.int64),
                        kde_x=None,
                        kde_y=None,
                        kde_peak=None,
                        kde_bandwidth=None,
                        n_frames_total=20,
                        n_frames_used=20,
                    ),
                    DistancePairPayload(
                        pair_label="auto1",
                        selection1="resid 156 and name NE2",
                        selection2="resid 77 and name OG",
                        distances=distances,
                        mean_distance=3.0,
                        std_distance=0.0,
                        median_distance=3.0,
                        min_distance=3.0,
                        max_distance=3.0,
                        sem_distance=None,
                        correlation_time=None,
                        correlation_time_unit=None,
                        n_independent_frames=None,
                        statistical_inefficiency=None,
                        autocorrelation_warning=None,
                        threshold=3.5,
                        fraction_below_threshold=1.0,
                        histogram_edges=np.asarray([2.5, 3.0, 3.5], dtype=np.float64),
                        histogram_counts=np.asarray([0, 20], dtype=np.int64),
                        kde_x=None,
                        kde_y=None,
                        kde_peak=None,
                        kde_bandwidth=None,
                        n_frames_total=20,
                        n_frames_used=20,
                    ),
                ],
            )
        )
        window = MagicMock(timestep_ps=10.0, step=1)

        with patch(
            "polyzymd.analyses.shared.autocorrelation.estimate_correlation_time",
            side_effect=RuntimeError("unexpected failure"),
        ):
            with pytest.raises(RuntimeError, match="unexpected failure"):
                triad_analysis.summarize_replicate(ctx, 1, runner, window)


# ============================================================================
# Test: aggregate
# ============================================================================


class TestAggregate:
    """Test CatalyticTriadAnalysis.aggregate."""

    def test_aggregates_results(self, triad_analysis, condition, tmp_path, default_settings):
        """Test aggregate produces a TriadAggregatedResult with correct stats."""
        settings_fp = settings_fingerprint(default_settings)
        results = [
            _make_mock_triad_result(i, sim_contact=0.6 + 0.05 * i, settings_fp=settings_fp)
            for i in range(1, 4)
        ]
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
        assert [pair.pair_label for pair in result.pair_results] == [
            "Asp133-His156",
            "His156-Ser77",
        ]
        assert len(result.per_replicate_simultaneous) == 3
        # Mean of 0.65, 0.70, 0.75
        assert result.overall_simultaneous_contact == pytest.approx(0.7, abs=0.01)
        assert result.sem_simultaneous_contact > 0

        stale = result.model_copy(
            update={"config_hash": "unknown", "settings_fingerprint": "deadbeef"}
        )
        with pytest.raises(AggregateValidationError, match="settings fingerprint mismatch"):
            triad_analysis.validate_aggregated_result(
                stale,
                condition=condition,
                settings=default_settings,
                equilibration="10ns",
                expected_replicates=(1, 2, 3),
            )

    def test_aggregate_saves_result_file(
        self, triad_analysis, condition, tmp_path, default_settings
    ):
        """Verify aggregated result is saved to the output directory."""
        settings_fp = settings_fingerprint(default_settings)
        results = [
            _make_mock_triad_result(i, sim_contact=0.65, settings_fp=settings_fp)
            for i in range(1, 3)
        ]
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
        assert "triad_LipA_triad_reps1-2_eq10ns.json" in json_files[0].name

    def test_aggregate_preserves_per_replicate_contact(
        self, triad_analysis, condition, tmp_path, default_settings
    ):
        """Verify per_replicate_simultaneous has the correct values."""
        settings_fp = settings_fingerprint(default_settings)
        results = [
            _make_mock_triad_result(i, sim_contact=0.5 + 0.1 * i, settings_fp=settings_fp)
            for i in range(1, 4)
        ]

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

    def test_aggregate_rejects_legacy_results_missing_settings_fingerprint(
        self, triad_analysis, condition, tmp_path, default_settings
    ):
        """Aggregation should reject legacy triad results without settings identity."""

        results = [_make_mock_triad_result(i, sim_contact=0.65) for i in range(1, 3)]
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=default_settings,
        )

        with pytest.raises(ValueError, match="missing settings fingerprints"):
            triad_analysis.aggregate(ctx, results)

    def test_aggregate_rejects_pair_count_mismatch(
        self, triad_analysis, condition, tmp_path, default_settings
    ):
        """Aggregation should reject replicate results with missing triad pairs."""

        settings_fp = settings_fingerprint(default_settings)
        results = [
            _make_mock_triad_result(1, sim_contact=0.65, settings_fp=settings_fp),
            _make_mock_triad_result(2, sim_contact=0.70, settings_fp=settings_fp),
        ]
        results[1] = results[1].model_copy(update={"pair_results": results[1].pair_results[:1]})
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=default_settings,
        )

        with pytest.raises(ValueError, match="configured pair schema"):
            triad_analysis.aggregate(ctx, results)

    def test_aggregate_rejects_pair_order_and_threshold_mismatch(
        self, triad_analysis, condition, tmp_path, default_settings
    ):
        """Aggregation should reject pair ordering or threshold drift before aggregation."""

        settings_fp = settings_fingerprint(default_settings)
        results = [
            _make_mock_triad_result(1, sim_contact=0.65, settings_fp=settings_fp),
            _make_mock_triad_result(2, sim_contact=0.70, settings_fp=settings_fp),
        ]
        mismatched_pairs = [
            results[1].pair_results[1].model_copy(update={"threshold": 4.0}),
            results[1].pair_results[0],
        ]
        results[1] = results[1].model_copy(update={"pair_results": mismatched_pairs})
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=default_settings,
        )

        with pytest.raises(ValueError, match="configured pair schema"):
            triad_analysis.aggregate(ctx, results)


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
    """Test CatalyticTriadAnalysis.plot delegates to inlined plotting helpers."""

    def test_threshold_bars_overlay_pair_and_simultaneous_replicates(self):
        """Threshold bars should pass per-replicate percentages to shared scatter."""
        import matplotlib.pyplot as plt

        from polyzymd.analyses.catalytic_triad._plotters import plot_triad_threshold_bars
        from polyzymd.config.comparison import PlotSettings

        pair_a = MagicMock(
            pair_label="Asp-His",
            selection1="resid 1",
            selection2="resid 2",
            overall_fraction_below=0.5,
            sem_fraction_below=0.05,
            per_replicate_fractions_below=[0.4, 0.6],
        )
        pair_b = MagicMock(
            pair_label="His-Ser",
            selection1="resid 2",
            selection2="resid 3",
            overall_fraction_below=0.25,
            sem_fraction_below=0.02,
            per_replicate_fractions_below=[0.2, 0.3],
        )
        result = MagicMock(
            pair_results=[pair_a, pair_b],
            overall_simultaneous_contact=0.1,
            sem_simultaneous_contact=0.01,
            per_replicate_simultaneous=[0.08, 0.12],
            threshold=3.5,
        )
        result.get_pair_labels.return_value = ["Asp-His", "His-Ser"]

        with patch(
            "polyzymd.analyses.catalytic_triad._plotters.scatter_replicate_values"
        ) as mock_scatter:
            fig = plot_triad_threshold_bars(
                [result],
                ["Control"],
                colors=["blue"],
                plot_settings=PlotSettings(),
            )

        mock_scatter.assert_called_once()
        assert mock_scatter.call_args.args[2] == [[40.0, 60.0], [20.0, 30.0], [8.0, 12.0]]
        plt.close(fig)

    def test_threshold_bars_aligns_out_of_order_sparse_pair_replicates(self):
        """Threshold dots should align by pair identity and reserve the final slot."""
        import matplotlib.pyplot as plt

        from polyzymd.analyses.catalytic_triad._plotters import plot_triad_threshold_bars
        from polyzymd.config.comparison import PlotSettings

        ref_pair_a = MagicMock(
            pair_label="Asp-His",
            selection1="resid 1",
            selection2="resid 2",
            overall_fraction_below=0.5,
            sem_fraction_below=0.05,
            per_replicate_fractions_below=[0.4, 0.6],
        )
        ref_pair_b = MagicMock(
            pair_label="His-Ser",
            selection1="resid 2",
            selection2="resid 3",
            overall_fraction_below=0.25,
            sem_fraction_below=0.02,
            per_replicate_fractions_below=[0.2, 0.3],
        )
        treatment_pair_b = MagicMock(
            pair_label="His-Ser",
            selection1="resid 2",
            selection2="resid 3",
            overall_fraction_below=0.35,
            sem_fraction_below=0.03,
            per_replicate_fractions_below=[0.3, 0.4],
        )
        treatment_extra_pair = MagicMock(
            pair_label="Extra-Pair",
            selection1="resid 4",
            selection2="resid 5",
            overall_fraction_below=0.95,
            sem_fraction_below=0.01,
            per_replicate_fractions_below=[0.9, 1.0],
        )
        control = MagicMock(
            pair_results=[ref_pair_a, ref_pair_b],
            overall_simultaneous_contact=0.1,
            sem_simultaneous_contact=0.01,
            per_replicate_simultaneous=[0.08, 0.12],
            threshold=3.5,
        )
        control.get_pair_labels.return_value = ["Asp-His", "His-Ser"]
        treatment = MagicMock(
            pair_results=[treatment_pair_b, treatment_extra_pair],
            overall_simultaneous_contact=0.2,
            sem_simultaneous_contact=0.02,
            per_replicate_simultaneous=[0.18, 0.22],
            threshold=3.5,
        )
        treatment.get_pair_labels.return_value = ["His-Ser", "Extra-Pair"]

        with patch(
            "polyzymd.analyses.catalytic_triad._plotters.scatter_replicate_values"
        ) as mock_scatter:
            fig = plot_triad_threshold_bars(
                [control, treatment],
                ["Control", "Treatment"],
                colors=["blue", "orange"],
                plot_settings=PlotSettings(),
            )

        treatment_replicates = mock_scatter.call_args_list[1].args[2]
        assert treatment_replicates == [[], [30.0, 40.0], [18.0, 22.0]]
        assert [90.0, 100.0] not in treatment_replicates
        treatment_heights = [patch_obj.get_height() for patch_obj in fig.axes[0].patches[3:]]
        assert treatment_heights == [0.0, 35.0, 20.0]
        plt.close(fig)

    def test_threshold_bars_from_data_passes_plot_settings(self, tmp_path):
        """Data-backed threshold bars should preserve caller plot settings."""
        from polyzymd.analyses.catalytic_triad._plotters import (
            plot_triad_threshold_bars_from_data,
        )
        from polyzymd.config.comparison import PlotSettings

        plot_settings = PlotSettings()
        mock_result = MagicMock()
        mock_figure = MagicMock()

        with (
            patch(
                "polyzymd.analyses.catalytic_triad._plotters._load_aggregated_results",
                return_value={"Control": mock_result},
            ),
            patch(
                "polyzymd.analyses.catalytic_triad._plotters.plot_triad_threshold_bars",
                return_value=mock_figure,
            ) as mock_plot,
            patch(
                "polyzymd.analyses.shared.plotting.save_figure", return_value=tmp_path / "bars.png"
            ),
        ):
            plot_triad_threshold_bars_from_data({}, ["Control"], tmp_path, plot_settings)

        assert mock_plot.call_args.kwargs["plot_settings"] is plot_settings

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
    @patch("polyzymd.analyses.catalytic_triad.TrajectoryLoader")
    def test_run_analysis_lifecycle(
        self,
        MockLoader,
        mock_hash,
        mock_validate_hash,
        mock_version,
        triad_analysis,
        condition,
        tmp_path,
        default_settings,
    ):
        """Test run_replicate -> aggregate via run_analysis()."""
        from polyzymd.analyses.orchestrator import run_analysis

        # Mock TrajectoryLoader for the runner seam
        mock_loader_inst = MagicMock()
        MockLoader.return_value = mock_loader_inst
        mock_u = MagicMock()
        mock_traj = MagicMock()
        mock_traj.__len__ = MagicMock(return_value=2000)
        mock_u.trajectory = mock_traj
        mock_loader_inst.load_universe.return_value = mock_u
        mock_loader_inst.get_timestep.return_value = 10.0

        mock_runner = MagicMock()
        mock_runner.results = object()
        mock_runner.run.return_value = mock_runner

        def make_triad_result(_ctx, replicate, _runner, _window):
            return _make_mock_triad_result(
                replicate,
                sim_contact=0.60 + 0.05 * (replicate - 1),
                settings_fp=settings_fingerprint(default_settings),
            )

        output_dir = tmp_path / "analysis" / "no_polymer" / "catalytic_triad"
        with (
            patch.object(
                triad_analysis, "build_runner", return_value=mock_runner
            ) as mock_build_runner,
            patch.object(
                triad_analysis,
                "summarize_replicate",
                side_effect=make_triad_result,
            ) as mock_summarize,
        ):
            result = run_analysis(
                triad_analysis,
                condition,
                settings=default_settings,
                equilibration="10ns",
                output_dir=output_dir,
            )

        assert MockLoader.call_count == 3
        assert mock_build_runner.call_count == 3
        assert mock_summarize.call_count == 3

        # Verify aggregate produced a valid result
        assert result.n_replicates == 3
        assert result.replicates == [1, 2, 3]
        assert result.triad_name == "LipA_triad"
        assert len(result.pair_results) == 2
        assert len(result.per_replicate_simultaneous) == 3
        # Mocked replicate summaries should remain valid contact fractions
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
