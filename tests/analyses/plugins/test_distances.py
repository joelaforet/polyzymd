"""Tests for the distances analysis plugin.

Covers discovery, settings, compute-stage dispatch, aggregate, compare (full override),
plot delegation, AggregatedResultClass, artifact deserialization, and lifecycle.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


class TestDiscovery:
    """The plugin is found by the automatic discovery system."""

    def test_discovery_by_name(self):
        from polyzymd.analyses.discovery import get_analysis
        from polyzymd.analyses.distances import DistancesAnalysis

        cls = get_analysis("distances")
        assert cls is DistancesAnalysis
        assert cls.name == "distances"

    def test_listed(self):
        from polyzymd.analyses.discovery import list_analyses

        analyses = list_analyses()
        names = list(analyses.keys())
        assert "distances" in names

    def test_all_names(self):
        from polyzymd.analyses.discovery import list_all_names

        names = list_all_names()
        assert "distances" in names


# ---------------------------------------------------------------------------
# Class Attributes
# ---------------------------------------------------------------------------


class TestClassAttributes:
    """Verify ClassVar declarations on DistancesAnalysis."""

    def test_name(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        assert DistancesAnalysis.name == "distances"

    def test_settings_type(self):
        from polyzymd.analyses.distances import DistancesAnalysis, DistancesSettings

        assert DistancesAnalysis.Settings is DistancesSettings

    def test_aliases(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        assert DistancesAnalysis.aliases == ()

    def test_dependencies(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        assert DistancesAnalysis.dependencies == ()

    def test_min_replicates(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        assert DistancesAnalysis.min_replicates == 1


def test_singleton_distance_pairwise_not_testable() -> None:
    """Distances singleton pairwise metrics should be marked not testable."""
    from polyzymd.analyses.distances import DistancesAnalysis

    pair_a = SimpleNamespace(
        mean_distance=3.0,
        fraction_below_threshold=0.4,
        per_replicate_means=[3.0],
        per_replicate_fractions=[0.4],
    )
    pair_b = SimpleNamespace(
        mean_distance=2.5,
        fraction_below_threshold=0.7,
        per_replicate_means=[2.5],
        per_replicate_fractions=[0.7],
    )

    result = DistancesAnalysis._compare_pair("Pair", "A", "B", pair_a, pair_b)

    assert result.distance_testable is False
    assert result.distance_significant is False
    assert result.distance_note is not None
    assert result.fraction_testable is False
    assert result.fraction_significant is False
    assert result.fraction_note is not None


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class TestSettings:
    """Validate DistancesSettings and DistancePairSettings."""

    def test_minimal_valid(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        s = DistancesSettings(
            pairs=[
                DistancePairSettings(
                    label="Pair1",
                    selection_a="resid 10 and name CA",
                    selection_b="resid 20 and name CA",
                )
            ]
        )
        assert s.threshold == 3.5
        assert s.use_pbc is True
        assert s.align_trajectory is True
        assert len(s.pairs) == 1

    def test_empty_pairs_rejected(self):
        from polyzymd.analyses.distances import DistancesSettings

        with pytest.raises(ValueError):
            DistancesSettings(pairs=[])

    def test_get_pair_selections(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        s = DistancesSettings(
            pairs=[
                DistancePairSettings(label="A", selection_a="sel1", selection_b="sel2"),
                DistancePairSettings(label="B", selection_a="sel3", selection_b="sel4"),
            ]
        )
        assert s.get_pair_selections() == [("sel1", "sel2"), ("sel3", "sel4")]

    def test_get_pair_labels(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        s = DistancesSettings(
            pairs=[
                DistancePairSettings(label="Alpha", selection_a="a", selection_b="b"),
                DistancePairSettings(label="Beta", selection_a="c", selection_b="d"),
            ]
        )
        assert s.get_pair_labels() == ["Alpha", "Beta"]

    def test_get_pair_thresholds_fallback(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        s = DistancesSettings(
            threshold=5.0,
            pairs=[
                DistancePairSettings(label="A", selection_a="a", selection_b="b", threshold=3.0),
                DistancePairSettings(
                    label="B", selection_a="c", selection_b="d"
                ),  # No per-pair threshold
            ],
        )
        thresholds = s.get_pair_thresholds()
        assert thresholds == [3.0, 5.0]  # per-pair, then fallback

    def test_alignment_mode_validation(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        with pytest.raises(ValueError):
            DistancesSettings(
                pairs=[DistancePairSettings(label="A", selection_a="a", selection_b="b")],
                alignment_mode="invalid",
            )

    def test_alignment_frame_required(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        with pytest.raises(ValueError):
            DistancesSettings(
                pairs=[DistancePairSettings(label="A", selection_a="a", selection_b="b")],
                alignment_mode="frame",
                # alignment_frame not provided
            )

    def test_alignment_config(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        s = DistancesSettings(
            pairs=[DistancePairSettings(label="A", selection_a="a", selection_b="b")],
            align_trajectory=True,
            alignment_selection="backbone",
            alignment_mode="centroid",
        )
        cfg = s.get_alignment_config()
        assert cfg.enabled is True
        assert cfg.selection == "backbone"

    def test_serialization_roundtrip(self):
        from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings

        s = DistancesSettings(
            threshold=4.0,
            pairs=[
                DistancePairSettings(
                    label="X",
                    selection_a="resid 1",
                    selection_b="resid 2",
                    threshold=3.0,
                    below_label="Bound",
                    above_label="Unbound",
                )
            ],
            use_pbc=False,
        )
        d = s.model_dump()
        s2 = DistancesSettings.model_validate(d)
        assert s2.threshold == 4.0
        assert s2.pairs[0].below_label == "Bound"
        assert s2.use_pbc is False


# ---------------------------------------------------------------------------
# run_replicate
# ---------------------------------------------------------------------------


def _make_mock_distance_result(
    replicate: int = 1,
    n_pairs: int = 2,
    pair_schemas: list[tuple[str, str, str]] | None = None,
):
    """Create a mock DistanceResult."""
    mock = MagicMock()
    mock.replicate = replicate
    mock.config_hash = "abc123"
    mock.equilibration_time = 100.0
    mock.equilibration_unit = "ns"

    pair_results = []
    for i in range(n_pairs):
        pr = MagicMock()
        if pair_schemas is None:
            pair_label = f"pair{i}"
            selection1 = f"sel_a_{i}"
            selection2 = f"sel_b_{i}"
        else:
            pair_label, selection1, selection2 = pair_schemas[i]

        pr.pair_label = pair_label
        pr.selection1 = selection1
        pr.selection2 = selection2
        pr.mean_distance = 3.5 + i * 0.5
        pr.std_distance = 0.5
        pr.median_distance = 3.4 + i * 0.5
        pr.min_distance = 2.0 + i * 0.3
        pr.max_distance = 6.0 + i * 0.3
        pr.sem_distance = 0.1
        pr.threshold = 3.5
        pr.fraction_below_threshold = 0.6 - i * 0.1
        pr.kde_peak = 3.3 + i * 0.5
        pr.per_replicate_means = [3.5 + i * 0.5]
        pair_results.append(pr)

    mock.pair_results = pair_results
    return mock


def _make_distance_cache_result(
    *,
    config_hash: str,
    pairs: list[tuple[str, str]],
    thresholds: list[float | None],
):
    """Create a concrete ``DistanceResult`` for cache-identity tests."""
    from polyzymd.analyses.distances._results import DistancePairResult, DistanceResult

    pair_results = []
    for idx, ((selection1, selection2), threshold) in enumerate(
        zip(pairs, thresholds, strict=True)
    ):
        pair_results.append(
            DistancePairResult(
                config_hash=config_hash,
                polyzymd_version="1.0.0",
                replicate=1,
                equilibration_time=0.0,
                equilibration_unit="ns",
                selection_string=f"{selection1} : {selection2}",
                pair_label=f"pair{idx}",
                selection1=selection1,
                selection2=selection2,
                distances=[3.0, 3.1, 3.2],
                mean_distance=3.1,
                std_distance=0.1,
                median_distance=3.1,
                min_distance=3.0,
                max_distance=3.2,
                sem_distance=0.05,
                correlation_time=None,
                correlation_time_unit=None,
                n_independent_frames=None,
                statistical_inefficiency=None,
                autocorrelation_warning=None,
                threshold=threshold,
                fraction_below_threshold=1.0,
                histogram_edges=[3.0, 3.1, 3.2],
                histogram_counts=[1, 2],
                kde_x=None,
                kde_y=None,
                kde_peak=None,
                kde_bandwidth=None,
                n_frames_total=3,
                n_frames_used=3,
            )
        )

    return DistanceResult(
        config_hash=config_hash,
        polyzymd_version="1.0.0",
        replicate=1,
        equilibration_time=0.0,
        equilibration_unit="ns",
        selection_string="; ".join(f"({a} : {b})" for a, b in pairs),
        pair_results=pair_results,
        n_frames_total=3,
        n_frames_used=3,
        trajectory_files=["/fake/traj.dcd"],
    )


def _make_distance_artifacts(tmp_path, condition_label, settings, n_reps: int = 3):
    """Create canonical distance replicate artifacts with NPZ sidecars."""
    import numpy as np

    from polyzymd.analyses._framework.cache_identity import settings_fingerprint
    from polyzymd.analyses.mda import ArtifactStore, ReplicateArtifact

    analysis_dir = tmp_path
    artifacts = []
    thresholds = settings.get_pair_thresholds()
    for rep in range(1, n_reps + 1):
        run_dir = analysis_dir / f"run_{rep}"
        store = ArtifactStore(run_dir)
        pair_payloads = []
        matrix_rows = []
        for pair_idx, (pair_setting, threshold) in enumerate(
            zip(settings.pairs, thresholds, strict=True)
        ):
            distances = np.asarray(
                [3.0 + pair_idx + rep * 0.01, 3.2 + pair_idx + rep * 0.01],
                dtype=np.float64,
            )
            matrix_rows.append(distances)
            pair_payloads.append(
                {
                    "pair_label": pair_setting.label,
                    "selection1": pair_setting.selection_a,
                    "selection2": pair_setting.selection_b,
                    "mean_distance": float(np.mean(distances)),
                    "std_distance": float(np.std(distances)),
                    "median_distance": float(np.median(distances)),
                    "min_distance": float(np.min(distances)),
                    "max_distance": float(np.max(distances)),
                    "sem_distance": 0.1,
                    "correlation_time": None,
                    "correlation_time_unit": None,
                    "n_independent_frames": None,
                    "statistical_inefficiency": None,
                    "autocorrelation_warning": None,
                    "threshold": threshold,
                    "fraction_below_threshold": 0.6 - pair_idx * 0.1,
                    "histogram_edges": [3.0, 3.5, 4.0],
                    "histogram_counts": [1, 1],
                    "kde_x": None,
                    "kde_y": None,
                    "kde_peak": 3.3 + pair_idx + rep * 0.01,
                    "kde_bandwidth": None,
                    "n_frames_total": 2,
                    "n_frames_used": 2,
                }
            )
        sidecar = store.write_npz_sidecar(
            "sidecars/00_distances.npz",
            distance_matrix=np.vstack(matrix_rows),
            frames=np.asarray([0, 1], dtype=np.int64),
            time_ns=np.asarray([0.0, 0.01], dtype=np.float64),
            pair_labels=np.asarray([pair.label for pair in settings.pairs]),
            metadata={"kind": "distance_matrix"},
        )
        artifact = ReplicateArtifact(
            analysis_name="distances",
            condition_label=condition_label,
            replicate=rep,
            payload={
                "pairs": pair_payloads,
                "pair_results": pair_payloads,
                "metrics": {
                    f"pair{idx}.mean_distance": pair["mean_distance"]
                    for idx, pair in enumerate(pair_payloads)
                },
                "replicate_metrics": {},
                "n_frames_total": 2,
                "n_frames_used": 2,
            },
            sidecars=[sidecar],
            provenance={"frame_selection": {"start": 0, "stop": 2, "step": 1}},
            metadata={
                "settings_fingerprint": settings_fingerprint(settings),
                "config_hash": "hash123",
                "polyzymd_version": "1.0.0-test",
                "equilibration_time": 100.0,
                "equilibration_unit": "ns",
                "selection_string": "; ".join(
                    f"({pair.selection_a} : {pair.selection_b})" for pair in settings.pairs
                ),
            },
        )
        artifacts.append(artifact)
    return artifacts


class TestDistanceResultFactories:
    """Factory helpers preserve the flat distance result schema."""

    @staticmethod
    def _make_payload():
        """Create a representative runner payload for factory tests."""
        import numpy as np

        from polyzymd.analyses.distances._mda import DistancePairPayload

        return DistancePairPayload(
            pair_label="pair0",
            selection1="sel_a",
            selection2="sel_b",
            distances=np.asarray([3.0, 3.2, 3.4], dtype=np.float64),
            mean_distance=3.2,
            std_distance=0.1,
            median_distance=3.2,
            min_distance=3.0,
            max_distance=3.4,
            sem_distance=0.05,
            correlation_time=20.0,
            correlation_time_unit="ps",
            n_independent_frames=12,
            statistical_inefficiency=1.5,
            autocorrelation_warning=None,
            threshold=3.5,
            fraction_below_threshold=1.0,
            histogram_edges=np.asarray([3.0, 3.2, 3.4], dtype=np.float64),
            histogram_counts=np.asarray([1, 2], dtype=np.int64),
            kde_x=np.asarray([3.0, 3.1, 3.2], dtype=np.float64),
            kde_y=np.asarray([0.1, 0.2, 0.1], dtype=np.float64),
            kde_peak=3.1,
            kde_bandwidth=0.2,
            n_frames_total=100,
            n_frames_used=90,
        )

    @staticmethod
    def _make_metadata(replicate: int | None = 1):
        """Create common result metadata for factory tests."""
        from polyzymd.analyses.distances._results import DistanceResultMetadata

        return DistanceResultMetadata(
            config_hash="hash123",
            polyzymd_version="1.3.0-test",
            replicate=replicate,
            equilibration_time=10.0,
            equilibration_unit="ns",
        )

    def test_pair_factory_matches_direct_constructor_schema(self):
        from polyzymd.analyses.distances._results import DistancePairResult

        payload = self._make_payload()
        metadata = self._make_metadata(replicate=None)

        factory_result = DistancePairResult.from_runner_payload(metadata, payload)
        direct_result = DistancePairResult(
            config_hash="hash123",
            polyzymd_version="1.3.0-test",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="sel_a : sel_b",
            pair_label="pair0",
            selection1="sel_a",
            selection2="sel_b",
            distances=[3.0, 3.2, 3.4],
            mean_distance=3.2,
            std_distance=0.1,
            median_distance=3.2,
            min_distance=3.0,
            max_distance=3.4,
            sem_distance=0.05,
            correlation_time=20.0,
            correlation_time_unit="ps",
            n_independent_frames=12,
            statistical_inefficiency=1.5,
            autocorrelation_warning=None,
            threshold=3.5,
            fraction_below_threshold=1.0,
            histogram_edges=[3.0, 3.2, 3.4],
            histogram_counts=[1, 2],
            kde_x=[3.0, 3.1, 3.2],
            kde_y=[0.1, 0.2, 0.1],
            kde_peak=3.1,
            kde_bandwidth=0.2,
            n_frames_total=100,
            n_frames_used=90,
        )

        assert factory_result.model_dump(exclude={"created_at"}) == direct_result.model_dump(
            exclude={"created_at"}
        )

    def test_pair_factory_can_omit_distributions(self):
        from polyzymd.analyses.distances._results import DistancePairResult

        result = DistancePairResult.from_runner_payload(
            self._make_metadata(replicate=None),
            self._make_payload(),
            store_distributions=False,
        )

        assert result.distances is None
        assert result.histogram_counts == [1, 2]

    def test_factory_results_have_flat_keys(self):
        from polyzymd.analyses.distances._results import DistancePairResult

        result = DistancePairResult.from_runner_payload(
            self._make_metadata(replicate=None),
            self._make_payload(),
        )
        dumped = result.model_dump()

        assert "metadata" not in dumped
        assert "payload" not in dumped
        assert "stats" not in dumped
        assert dumped["config_hash"] == "hash123"
        assert dumped["pair_label"] == "pair0"

    def test_result_factory_stringifies_paths(self):
        from polyzymd.analyses.distances._results import DistancePairResult, DistanceResult

        pair_result = DistancePairResult.from_runner_payload(
            self._make_metadata(replicate=None),
            self._make_payload(),
        )

        result = DistanceResult.from_pair_results(
            self._make_metadata(replicate=1),
            [pair_result],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=[Path("/fake/traj.dcd"), "relative.xtc"],
        )

        assert result.trajectory_files == ["/fake/traj.dcd", "relative.xtc"]

    def test_aggregate_factory_preserves_flat_schema(self):
        from types import SimpleNamespace

        from polyzymd.analyses.distances._results import (
            DistanceAggregatedResult,
            DistancePairAggregatedResult,
        )

        source_pair = SimpleNamespace(pair_label="pair0", selection1="sel_a", selection2="sel_b")
        stats = SimpleNamespace(
            mean_stats=SimpleNamespace(mean=3.2, sem=0.05),
            median_stats=SimpleNamespace(mean=3.1, sem=0.04),
            fraction_stats=SimpleNamespace(mean=0.8, sem=0.02),
            kde_peak_stats=SimpleNamespace(mean=3.0, sem=0.03),
            per_rep_means=[3.1, 3.3],
            per_rep_stds=[0.1, 0.2],
            per_rep_medians=[3.0, 3.2],
            per_rep_fractions=[0.75, 0.85],
            per_rep_kde_peaks=[2.9, 3.1],
        )

        pair_result = DistancePairAggregatedResult.from_aggregated_stats(
            self._make_metadata(replicate=None),
            source_pair,
            stats,
            replicates=(1, 2),
            threshold=3.5,
        )
        result = DistanceAggregatedResult.from_pair_results(
            self._make_metadata(replicate=None),
            [pair_result],
            replicates=(1, 2),
            source_result_files=[Path("/fake/run_1.json")],
        )

        dumped = result.model_dump()
        assert dumped["replicates"] == [1, 2]
        assert dumped["source_result_files"] == ["/fake/run_1.json"]
        assert dumped["pair_results"][0]["overall_mean"] == 3.2
        assert "metadata" not in dumped
        assert "stats" not in dumped

    def test_direct_constructors_still_supported(self):
        from polyzymd.analyses.distances._results import DistancePairResult, DistanceResult

        pair_result = DistancePairResult(
            config_hash="hash123",
            polyzymd_version="1.3.0-test",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="sel_a : sel_b",
            pair_label="pair0",
            selection1="sel_a",
            selection2="sel_b",
            distances=None,
            mean_distance=3.2,
            std_distance=0.1,
            median_distance=3.2,
            min_distance=3.0,
            max_distance=3.4,
            n_frames_total=100,
            n_frames_used=90,
        )
        result = DistanceResult(
            config_hash="hash123",
            polyzymd_version="1.3.0-test",
            replicate=1,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="(sel_a : sel_b)",
            pair_results=[pair_result],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=[],
        )

        assert result.pair_results[0].pair_label == "pair0"

    def test_distance_pair_factory_replicate_is_none_for_distances(self):
        from polyzymd.analyses.distances._results import DistancePairResult

        result = DistancePairResult.from_runner_payload(
            self._make_metadata(replicate=1).with_replicate(None),
            self._make_payload(),
        )

        assert result.replicate is None


class TestRunReplicate:
    """The compute-stage dispatcher delegates to the MDAnalysis job seam."""

    def test_delegates_to_mda_job_seam(self, tmp_path):
        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )
        from polyzymd.analyses.mda import ReplicateArtifact

        analysis = DistancesAnalysis()
        settings = DistancesSettings(
            pairs=[
                DistancePairSettings(
                    label="Configured Label",
                    selection_a="sel_a",
                    selection_b="sel_b",
                    threshold=3.2,
                )
            ]
        )
        mock_sim_config = MagicMock()
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        ctx = ReplicateContext(
            condition=cond,
            replicate=1,
            sim_config=cond.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=False,
            settings=settings,
        )
        artifact = ReplicateArtifact(
            analysis_name="distances",
            condition_label="test",
            replicate=1,
            payload={"metrics": {"distance.mean": 1.0}},
        )
        with patch(
            "polyzymd.analyses.mda.lifecycle.run_mda_replicate_jobs",
            return_value=artifact,
        ) as mock_run:
            result = analysis._run_compute_stage(ctx, 1)

        assert result is artifact
        mock_run.assert_called_once_with(analysis, ctx, 1)

    def test_build_mda_jobs_delegates_to_builder(self):
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()
        settings = DistancesSettings(
            pairs=[DistancePairSettings(label="P1", selection_a="resid 10", selection_b="resid 20")]
        )
        ctx = MagicMock(settings=settings)
        expected_jobs = [MagicMock()]

        with patch("polyzymd.analyses.distances.build_distance_jobs", return_value=expected_jobs):
            jobs = analysis.build_mda_jobs(ctx)

        assert jobs is expected_jobs

    def test_compute_distance_payloads_uses_effective_timestep_for_autocorrelation(
        self, monkeypatch
    ):
        """Distance autocorrelation should use raw timestep multiplied by frame stride."""
        from types import SimpleNamespace

        import numpy as np

        from polyzymd.analyses.distances import _mda as runner_module
        from polyzymd.analyses.distances._mda import compute_distance_payloads

        captured: dict[str, float] = {}

        class _FakeDistanceAnalysis:
            def __init__(self, **kwargs):
                self.results = SimpleNamespace(
                    distance_matrix=[], frames=[], times_ps=[], warnings=[]
                )

            def run(self, *, start: int, stop: int, step: int):
                del start, stop, step
                self.results.distance_matrix = np.asarray(
                    [np.linspace(1.0, 2.0, 20, dtype=np.float64)]
                )
                self.results.frames = np.arange(0, 60, 3, dtype=np.int64)
                return self

        monkeypatch.setattr(
            runner_module, "build_pair_distance_analysis", lambda **kwargs: _FakeDistanceAnalysis()
        )
        monkeypatch.setattr(
            runner_module,
            "resolve_distance_pairs",
            lambda **kwargs: [
                SimpleNamespace(
                    label="P1",
                    selection_a="sel_a",
                    selection_b="sel_b",
                    threshold=None,
                )
            ],
        )

        def _fake_estimate_correlation_time(_series, **kwargs):
            captured["timestep"] = kwargs["timestep"]
            return SimpleNamespace(
                tau=1.0,
                tau_unit="ps",
                n_independent=20,
                statistical_inefficiency=1.0,
                warning=None,
            )

        monkeypatch.setattr(
            "polyzymd.analyses.shared.autocorrelation.estimate_correlation_time",
            _fake_estimate_correlation_time,
        )
        universe = MagicMock(trajectory=list(range(30)))

        compute_distance_payloads(
            universe=universe,
            pairs=[("sel_a", "sel_b")],
            thresholds=[None],
            start=0,
            stop=30,
            step=3,
            timestep_ps=10.0,
            use_pbc=False,
            alignment=SimpleNamespace(enabled=False),
            pair_label_func=lambda _selection_a, _selection_b: "P1",
        )

        assert captured["timestep"] == pytest.approx(30.0)

    def test_collector_writes_summary_json_and_npz_sidecar(self, tmp_path):
        import numpy as np

        from polyzymd.analyses.base import Condition, ReplicateContext
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesSettings,
        )
        from polyzymd.analyses.distances._mda import DistanceArtifactCollector
        from polyzymd.analyses.mda import ArtifactStore, FrameSelection, MDACollectorContext
        from polyzymd.analyses.mda.job import MDABackendPolicy, MDAJobResult, MDAUniversePolicy

        settings = DistancesSettings(
            pairs=[
                DistancePairSettings(
                    label="Configured Label",
                    selection_a="sel_a",
                    selection_b="sel_b",
                    threshold=3.2,
                )
            ]
        )
        condition = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=MagicMock(),
        )
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,
            settings=settings,
        )
        analysis_obj = MagicMock()
        analysis_obj.results.distance_matrix = np.asarray([[3.0, 3.2, 3.4]], dtype=np.float64)
        analysis_obj.results.frames = np.asarray([0, 1, 2], dtype=np.int64)
        analysis_obj.results.times_ps = np.asarray([0.0, 10.0, 20.0], dtype=np.float64)
        analysis_obj.results.warnings = []
        job = MDAJobResult(
            name="pair_distances",
            analysis=analysis_obj,
            results=analysis_obj.results,
            run_kwargs={"start": 0, "stop": 3, "step": 1},
            frame_selection=FrameSelection(start=0, stop=3, step=1, n_frames_total=3),
            backend_policy=MDABackendPolicy(),
            universe_policy=MDAUniversePolicy(),
        )
        collector_ctx = MDACollectorContext(
            analysis_name="distances",
            replicate_context=ctx,
            frame_selection=FrameSelection(start=0, stop=3, step=1, n_frames_total=3),
            universe_policy=MDAUniversePolicy(condition_label="test", replicate=1),
            artifact_store=ArtifactStore(tmp_path / "run_1"),
            settings_fingerprint="abc123",
        )

        result = DistanceArtifactCollector()(collector_ctx, [job])

        assert result.replicate == 1
        assert result.payload["n_frames_used"] == 3
        assert result.payload["pairs"][0]["pair_label"] == "Configured Label"
        assert result.payload["pairs"][0]["mean_distance"] == pytest.approx(3.2)
        assert result.payload["pairs"][0]["fraction_below_threshold"] == pytest.approx(1.0 / 3.0)
        assert "distances" not in result.payload["pairs"][0]
        assert result.sidecars[0].metadata["kind"] == "distance_matrix"
        sidecar_path = ArtifactStore(tmp_path / "run_1").validate_sidecar(result.sidecars[0])
        with np.load(sidecar_path) as npz_data:
            assert npz_data["distance_matrix"].shape == (1, 3)
            assert npz_data["pair_labels"].tolist() == ["Configured Label"]
            assert npz_data["thresholds"].tolist() == pytest.approx([3.2])


# ---------------------------------------------------------------------------
# aggregate
# ---------------------------------------------------------------------------


class TestAggregate:
    """aggregate produces a DistanceAggregatedResult."""

    def _make_mock_results(self, settings, n_reps: int = 3):
        """Create mock per-replicate results for aggregation."""
        from polyzymd.analyses.distances import _make_pair_label

        thresholds = settings.get_pair_thresholds()
        results = []
        for rep in range(1, n_reps + 1):
            mock = MagicMock()
            mock.replicate = rep
            mock.config_hash = "hash123"
            mock.equilibration_time = 100.0
            mock.equilibration_unit = "ns"

            pair_results = []
            for pair_idx, (pair_setting, threshold) in enumerate(
                zip(settings.pairs, thresholds, strict=True)
            ):
                pr = MagicMock()
                pr.pair_label = _make_pair_label(pair_setting.selection_a, pair_setting.selection_b)
                pr.selection1 = pair_setting.selection_a
                pr.selection2 = pair_setting.selection_b
                pr.mean_distance = 3.5 + pair_idx * 0.5 + rep * 0.01
                pr.std_distance = 0.5
                pr.median_distance = 3.4 + pair_idx * 0.5 + rep * 0.01
                pr.threshold = threshold
                pr.fraction_below_threshold = 0.6 - pair_idx * 0.1
                pr.kde_peak = 3.3 + pair_idx * 0.5 + rep * 0.01
                pair_results.append(pr)

            mock.pair_results = pair_results
            results.append(mock)
        return results

    def test_aggregate_produces_result(self, tmp_path):
        from polyzymd.analyses._framework.cache_identity import settings_fingerprint
        from polyzymd.analyses.base import AggregateContext, AggregateValidationError, Condition
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()
        settings = DistancesSettings(
            pairs=[
                DistancePairSettings(label="P0", selection_a="sel_a_0", selection_b="sel_b_0"),
                DistancePairSettings(label="P1", selection_a="sel_a_1", selection_b="sel_b_1"),
            ]
        )

        output_dir = tmp_path / "aggregated"
        mock_sim_config = MagicMock()
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=mock_sim_config,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2, 3),
            output_dir=output_dir,
            equilibration="100ns",
            settings=settings,
        )

        mock_results = _make_distance_artifacts(tmp_path, "test", settings, n_reps=3)

        with patch(
            "polyzymd.analyses._framework.results_base.get_polyzymd_version",
            return_value="1.0.0-test",
        ):
            result = analysis.aggregate(ctx, mock_results)

        assert result is not None
        assert result.metadata["n_replicates"] == 3
        assert len(result.payload["pair_results"]) == 2
        assert [pair["pair_label"] for pair in result.payload["pair_results"]] == ["P0", "P1"]
        assert result.metadata["settings_fingerprint"] == settings_fingerprint(settings)

        stale = result.model_copy(
            update={"metadata": {**result.metadata, "settings_fingerprint": "deadbeef"}}
        )
        with pytest.raises(AggregateValidationError, match="settings fingerprint mismatch"):
            analysis.validate_aggregated_result(
                stale,
                condition=cond,
                settings=settings,
                equilibration="100ns",
                expected_replicates=(1, 2, 3),
            )

    def test_aggregate_returns_artifact_without_canonical_result(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()
        settings = DistancesSettings(
            pairs=[
                DistancePairSettings(label="P0", selection_a="a", selection_b="b"),
            ]
        )

        output_dir = tmp_path / "aggregated"
        mock_sim_config = MagicMock()
        cond = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=mock_sim_config,
        )
        ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=output_dir,
            equilibration="100ns",
            settings=settings,
        )

        mock_results = _make_distance_artifacts(tmp_path, "test", settings, n_reps=2)

        with patch(
            "polyzymd.analyses._framework.results_base.get_polyzymd_version",
            return_value="1.0.0-test",
        ):
            result = analysis.aggregate(ctx, mock_results)

        assert result is not None
        json_files = list(output_dir.glob("*.json"))
        assert json_files == []

    def test_aggregate_rejects_pair_count_mismatch(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()
        settings = DistancesSettings(
            pairs=[
                DistancePairSettings(label="P0", selection_a="sel_a_0", selection_b="sel_b_0"),
                DistancePairSettings(label="P1", selection_a="sel_a_1", selection_b="sel_b_1"),
            ]
        )
        condition = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="100ns",
            settings=settings,
        )
        results = _make_distance_artifacts(tmp_path, "test", settings, n_reps=2)
        results[1].payload["pairs"] = results[1].payload["pairs"][:1]
        results[1].payload["pair_results"] = results[1].payload["pairs"]

        with pytest.raises(ValueError, match="expected 2"):
            analysis.aggregate(ctx, results)

    def test_aggregate_rejects_pair_order_and_threshold_mismatch(self, tmp_path):
        from polyzymd.analyses.base import AggregateContext, Condition
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()
        settings = DistancesSettings(
            pairs=[
                DistancePairSettings(
                    label="P0",
                    selection_a="sel_a_0",
                    selection_b="sel_b_0",
                    threshold=3.0,
                ),
                DistancePairSettings(
                    label="P1",
                    selection_a="sel_a_1",
                    selection_b="sel_b_1",
                    threshold=4.0,
                ),
            ]
        )
        condition = Condition(
            label="test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="100ns",
            settings=settings,
        )
        results = _make_distance_artifacts(tmp_path, "test", settings, n_reps=2)
        results[1].payload["pairs"] = list(reversed(results[1].payload["pairs"]))
        results[1].payload["pair_results"] = results[1].payload["pairs"]

        with pytest.raises(ValueError, match="selection mismatch"):
            analysis.aggregate(ctx, results)


# ---------------------------------------------------------------------------
# compare (full override)
# ---------------------------------------------------------------------------


def _make_mock_agg_result(n_pairs: int = 2, n_reps: int = 3, offset: float = 0.0):
    """Create a canonical condition artifact for comparison tests."""
    from polyzymd.analyses._framework.cache_identity import settings_fingerprint
    from polyzymd.analyses.distances import DistancePairSettings, DistancesSettings
    from polyzymd.analyses.mda import ConditionArtifact

    labels = ["Alpha", "Beta", "Gamma", "Delta"]
    settings = DistancesSettings(
        pairs=[
            DistancePairSettings(
                label=labels[i] if i < len(labels) else f"pair{i}",
                selection_a=f"sel_a_{i}",
                selection_b=f"sel_b_{i}",
            )
            for i in range(n_pairs)
        ]
    )
    pair_results = []
    for i in range(n_pairs):
        pair_results.append(
            {
                "pair_label": f"pair{i}",
                "selection1": f"sel_a_{i}",
                "selection2": f"sel_b_{i}",
                "overall_mean": 3.5 + i * 0.5 + offset,
                "overall_sem": 0.05,
                "threshold": 3.5,
                "overall_fraction_below": 0.6 - i * 0.1 - offset * 0.1,
                "sem_fraction_below": 0.02,
                "per_replicate_means": [3.5 + i * 0.5 + offset + j * 0.01 for j in range(n_reps)],
                "per_replicate_fractions_below": [
                    0.6 - i * 0.1 - offset * 0.1 + j * 0.005 for j in range(n_reps)
                ],
                "replicates": list(range(1, n_reps + 1)),
                "n_replicates": n_reps,
            }
        )

    return ConditionArtifact(
        analysis_name="distances",
        condition_label="mock",
        replicates=list(range(1, n_reps + 1)),
        payload={"pair_results": pair_results, "pairs": pair_results, "n_replicates": n_reps},
        metadata={
            "settings_fingerprint": settings_fingerprint(settings),
            "config_hash": "unknown",
            "equilibration_time": 100.0,
            "equilibration_unit": "ns",
        },
    )


def _distance_condition_artifact_from_pair_objects(
    pairs: list[object], replicates: list[int], settings_fingerprint: str
):
    """Create a canonical distances condition artifact from pair-like objects."""

    from polyzymd.analyses.mda import ConditionArtifact

    pair_payloads = []
    for pair in pairs:
        pair_payloads.append(
            {
                "pair_label": pair.pair_label,
                "selection1": pair.selection1,
                "selection2": pair.selection2,
                "threshold": pair.threshold,
                "overall_mean": pair.overall_mean,
                "overall_sem": pair.overall_sem,
                "overall_fraction_below": pair.overall_fraction_below,
                "sem_fraction_below": pair.sem_fraction_below,
                "per_replicate_means": list(pair.per_replicate_means),
                "per_replicate_fractions_below": list(pair.per_replicate_fractions_below),
                "replicates": replicates,
                "n_replicates": len(replicates),
            }
        )
    return ConditionArtifact(
        analysis_name="distances",
        condition_label="Control",
        replicates=replicates,
        payload={
            "pair_results": pair_payloads,
            "pairs": pair_payloads,
            "n_replicates": len(replicates),
        },
        metadata={
            "settings_fingerprint": settings_fingerprint,
            "config_hash": "unknown",
            "equilibration_time": 100.0,
            "equilibration_unit": "ns",
        },
    )


class TestCompare:
    """compare() produces a DistanceComparisonResult with per-pair stats."""

    def _make_analysis_and_ctx(self, tmp_path, n_conditions: int = 3, control: str = "Control"):
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()

        settings = DistancesSettings(
            pairs=[
                DistancePairSettings(
                    label="Alpha",
                    selection_a="sel_a_0",
                    selection_b="sel_b_0",
                ),
                DistancePairSettings(
                    label="Beta",
                    selection_a="sel_a_1",
                    selection_b="sel_b_1",
                ),
            ]
        )

        conditions = []
        analysis_dirs = {}
        labels = ["Control", "Treatment_A", "Treatment_B"][:n_conditions]

        for label in labels:
            cond = Condition(
                label=label,
                config_path=Path(f"/tmp/{label}/config.yaml"),
                replicates=(1, 2, 3),
                sim_config=MagicMock(),
            )
            conditions.append(cond)
            agg_dir = tmp_path / label / "distances"
            agg_dir.mkdir(parents=True, exist_ok=True)
            analysis_dirs[label] = agg_dir

        ctx = ComparisonContext(
            name="test_comparison",
            conditions=conditions,
            excluded_conditions=[],
            control_label=control,
            analysis_dirs=analysis_dirs,
            results_dir=tmp_path / "results",
            equilibration="100ns",
            settings=settings,
            recompute=False,
        )

        return analysis, ctx, labels

    def test_compare_returns_result(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path)

        # Mock _load_aggregated_result to return different values per condition
        # Path structure: tmp_path / label / "distances" / "aggregated"
        # So agg_dir.parent.parent.name gives the condition label
        offsets = {"Control": 0.0, "Treatment_A": 0.5, "Treatment_B": -0.3}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        assert result is not None
        assert result.metric == "mean_distance"
        assert result.name == "test_comparison"
        assert result.n_pairs == 2
        assert len(result.pair_labels) == 2
        assert result.pair_labels == ["Alpha", "Beta"]
        assert len(result.conditions) == 3
        assert result.control_label == "Control"

    def test_compare_has_per_pair_rankings(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path)

        offsets = {"Control": 0.0, "Treatment_A": 0.5, "Treatment_B": -0.3}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        assert "Alpha" in result.ranking_by_pair
        assert "Beta" in result.ranking_by_pair
        # Treatment_B has lowest offset (-0.3) so should rank first for each pair
        assert result.ranking_by_pair["Alpha"][0] == "Treatment_B"

    def test_compare_pairwise_comparisons(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path)

        offsets = {"Control": 0.0, "Treatment_A": 0.5, "Treatment_B": -0.3}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        # With control, should have (n_conditions - 1) * n_pairs comparisons
        assert len(result.pairwise_comparisons) == 4  # 2 treatments * 2 pairs

        # Each comparison should have distance and fraction fields
        for comp in result.pairwise_comparisons:
            assert comp.condition_a == "Control"
            assert hasattr(comp, "distance_t_statistic")
            assert hasattr(comp, "distance_direction")

    def test_compare_anova_with_3_conditions(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path, n_conditions=3)

        offsets = {"Control": 0.0, "Treatment_A": 0.5, "Treatment_B": -0.3}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        # With 3+ conditions, ANOVA should be computed
        assert result.anova_by_pair is not None
        assert len(result.anova_by_pair) == 2  # one per pair

    def test_compare_no_anova_with_2_conditions(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path, n_conditions=2)

        offsets = {"Control": 0.0, "Treatment_A": 0.5}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        assert result.anova_by_pair is None

    def test_compare_returns_result_with_single_condition(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path, n_conditions=1)

        def mock_load(agg_dir):
            return _make_mock_agg_result(n_pairs=2)

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        assert result is not None
        assert len(result.conditions) == 1
        assert result.pairwise_comparisons == []

    def test_compare_no_control(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path, n_conditions=3, control=None)
        # Override control_label to None
        from polyzymd.analyses.base import ComparisonContext

        ctx = ComparisonContext(
            name=ctx.name,
            conditions=ctx.conditions,
            excluded_conditions=ctx.excluded_conditions,
            control_label=None,
            analysis_dirs=ctx.analysis_dirs,
            results_dir=ctx.results_dir,
            equilibration=ctx.equilibration,
            settings=ctx.settings,
            recompute=ctx.recompute,
        )

        offsets = {"Control": 0.0, "Treatment_A": 0.5, "Treatment_B": -0.3}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        assert result is not None
        # Without control, should have C(3,2) * n_pairs = 3 * 2 = 6 comparisons
        assert len(result.pairwise_comparisons) == 6

    def test_compare_fraction_ranking(self, tmp_path):
        analysis, ctx, labels = self._make_analysis_and_ctx(tmp_path, n_conditions=2)

        offsets = {"Control": 0.0, "Treatment_A": 0.5}

        def mock_load(agg_dir):
            label = agg_dir.parent.parent.name
            return _make_mock_agg_result(n_pairs=2, offset=offsets.get(label, 0.0))

        with patch.object(analysis, "_load_aggregated_result", side_effect=mock_load):
            result = analysis.compare(ctx)

        # fraction_ranking_by_pair should be present since mock has fractions
        assert result.fraction_ranking_by_pair is not None

    def test_compare_preserves_duplicate_selection_pair_labels_by_index(self, tmp_path):
        from polyzymd.analyses._framework.cache_identity import settings_fingerprint
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()
        settings = DistancesSettings(
            pairs=[
                DistancePairSettings(
                    label="Inner duplicate",
                    selection_a="resid 10 and name CA",
                    selection_b="resid 20 and name CA",
                    threshold=3.0,
                ),
                DistancePairSettings(
                    label="Outer duplicate",
                    selection_a="resid 10 and name CA",
                    selection_b="resid 20 and name CA",
                    threshold=5.0,
                ),
            ]
        )
        condition = Condition(
            label="Control",
            config_path=Path("/tmp/Control/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        ctx = ComparisonContext(
            name="duplicate_selection_labels",
            conditions=[condition],
            excluded_conditions=[],
            control_label="Control",
            analysis_dirs={"Control": tmp_path / "Control" / "distances"},
            results_dir=tmp_path / "results",
            equilibration="100ns",
            settings=settings,
            recompute=False,
        )
        aggregate = MagicMock()
        aggregate.config_hash = "unknown"
        aggregate.equilibration_time = 100.0
        aggregate.equilibration_unit = "ns"
        aggregate.settings_fingerprint = settings_fingerprint(settings)
        aggregate.replicates = [1, 2, 3]
        aggregate.n_replicates = 3
        aggregate.pair_results = []
        for label, threshold, offset in (
            ("Inner duplicate", 3.0, 0.0),
            ("Outer duplicate", 5.0, 1.0),
        ):
            pair_result = MagicMock()
            pair_result.pair_label = label
            pair_result.selection1 = "resid 10 and name CA"
            pair_result.selection2 = "resid 20 and name CA"
            pair_result.threshold = threshold
            pair_result.overall_mean = 4.0 + offset
            pair_result.overall_sem = 0.1
            pair_result.overall_fraction_below = 0.5 - offset * 0.1
            pair_result.sem_fraction_below = 0.02
            pair_result.per_replicate_means = [4.0 + offset, 4.1 + offset, 3.9 + offset]
            pair_result.per_replicate_fractions_below = [0.5, 0.55, 0.45]
            aggregate.pair_results.append(pair_result)

        artifact = _distance_condition_artifact_from_pair_objects(
            aggregate.pair_results, aggregate.replicates, aggregate.settings_fingerprint
        )
        with patch.object(analysis, "_load_aggregated_result", return_value=artifact):
            result = analysis.compare(ctx)

        assert result is not None
        summary = result.conditions[0]
        assert [pair.label for pair in summary.pair_summaries] == [
            "Inner duplicate",
            "Outer duplicate",
        ]
        assert summary.get_pair("Inner duplicate").threshold == pytest.approx(3.0)
        assert summary.get_pair("Outer duplicate").threshold == pytest.approx(5.0)
        assert result.ranking_by_pair["Inner duplicate"] == ["Control"]
        assert result.ranking_by_pair["Outer duplicate"] == ["Control"]

    def test_compare_remaps_duplicate_legacy_auto_labels_by_index(self, tmp_path):
        from polyzymd.analyses._framework.cache_identity import settings_fingerprint
        from polyzymd.analyses.base import ComparisonContext, Condition
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
            _make_pair_label,
        )

        selection_a = "resid 10 and name CA"
        selection_b = "resid 20 and name CA"
        auto_label = _make_pair_label(selection_a, selection_b)
        analysis = DistancesAnalysis()
        settings = DistancesSettings(
            pairs=[
                DistancePairSettings(
                    label=auto_label,
                    selection_a=selection_a,
                    selection_b=selection_b,
                    threshold=3.0,
                ),
                DistancePairSettings(
                    label="Distinct duplicate",
                    selection_a=selection_a,
                    selection_b=selection_b,
                    threshold=5.0,
                ),
            ]
        )
        condition = Condition(
            label="Control",
            config_path=Path("/tmp/Control/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )
        ctx = ComparisonContext(
            name="duplicate_legacy_auto_labels",
            conditions=[condition],
            excluded_conditions=[],
            control_label="Control",
            analysis_dirs={"Control": tmp_path / "Control" / "distances"},
            results_dir=tmp_path / "results",
            equilibration="100ns",
            settings=settings,
            recompute=False,
        )
        aggregate = MagicMock()
        aggregate.config_hash = "unknown"
        aggregate.equilibration_time = 100.0
        aggregate.equilibration_unit = "ns"
        aggregate.settings_fingerprint = settings_fingerprint(settings)
        aggregate.replicates = [1, 2, 3]
        aggregate.n_replicates = 3
        aggregate.pair_results = []
        for threshold, offset in ((3.0, 0.0), (5.0, 1.0)):
            pair_result = MagicMock()
            pair_result.pair_label = auto_label
            pair_result.selection1 = selection_a
            pair_result.selection2 = selection_b
            pair_result.threshold = threshold
            pair_result.overall_mean = 4.0 + offset
            pair_result.overall_sem = 0.1
            pair_result.overall_fraction_below = 0.5 - offset * 0.1
            pair_result.sem_fraction_below = 0.02
            pair_result.per_replicate_means = [4.0 + offset, 4.1 + offset, 3.9 + offset]
            pair_result.per_replicate_fractions_below = [0.5, 0.55, 0.45]
            aggregate.pair_results.append(pair_result)

        artifact = _distance_condition_artifact_from_pair_objects(
            aggregate.pair_results, aggregate.replicates, aggregate.settings_fingerprint
        )
        with patch.object(analysis, "_load_aggregated_result", return_value=artifact):
            result = analysis.compare(ctx)

        assert result is not None
        summary = result.conditions[0]
        assert [pair.label for pair in summary.pair_summaries] == [
            auto_label,
            "Distinct duplicate",
        ]
        assert summary.get_pair(auto_label).threshold == pytest.approx(3.0)
        assert summary.get_pair("Distinct duplicate").threshold == pytest.approx(5.0)
        assert result.ranking_by_pair[auto_label] == ["Control"]
        assert result.ranking_by_pair["Distinct duplicate"] == ["Control"]


# ---------------------------------------------------------------------------
# _compare_pair
# ---------------------------------------------------------------------------


class TestComparePair:
    """_compare_pair performs correct statistical comparison for one pair."""

    def test_distance_direction_closer(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        pair_a = MagicMock()
        pair_a.per_replicate_means = [4.0, 4.1, 3.9]
        pair_a.mean_distance = 4.0
        pair_a.per_replicate_fractions = None

        pair_b = MagicMock()
        pair_b.per_replicate_means = [3.0, 3.1, 2.9]
        pair_b.mean_distance = 3.0
        pair_b.per_replicate_fractions = None

        comp = DistancesAnalysis._compare_pair("test_pair", "A", "B", pair_a, pair_b)
        assert comp.distance_direction == "closer"
        assert comp.distance_percent_change < 0

    def test_distance_direction_farther(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        pair_a = MagicMock()
        pair_a.per_replicate_means = [3.0, 3.1, 2.9]
        pair_a.mean_distance = 3.0
        pair_a.per_replicate_fractions = None

        pair_b = MagicMock()
        pair_b.per_replicate_means = [5.0, 5.1, 4.9]
        pair_b.mean_distance = 5.0
        pair_b.per_replicate_fractions = None

        comp = DistancesAnalysis._compare_pair("test_pair", "A", "B", pair_a, pair_b)
        assert comp.distance_direction == "farther"

    def test_fraction_comparison_when_present(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        pair_a = MagicMock()
        pair_a.per_replicate_means = [4.0, 4.1, 3.9]
        pair_a.mean_distance = 4.0
        pair_a.per_replicate_fractions = [0.5, 0.55, 0.45]
        pair_a.fraction_below_threshold = 0.5

        pair_b = MagicMock()
        pair_b.per_replicate_means = [3.0, 3.1, 2.9]
        pair_b.mean_distance = 3.0
        pair_b.per_replicate_fractions = [0.7, 0.75, 0.65]
        pair_b.fraction_below_threshold = 0.7

        comp = DistancesAnalysis._compare_pair("test_pair", "A", "B", pair_a, pair_b)
        assert comp.fraction_t_statistic is not None
        assert comp.fraction_direction == "more_contact"

    def test_no_fraction_when_absent(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        pair_a = MagicMock()
        pair_a.per_replicate_means = [4.0, 4.1]
        pair_a.mean_distance = 4.0
        pair_a.per_replicate_fractions = None

        pair_b = MagicMock()
        pair_b.per_replicate_means = [3.0, 3.1]
        pair_b.mean_distance = 3.0
        pair_b.per_replicate_fractions = None

        comp = DistancesAnalysis._compare_pair("test_pair", "A", "B", pair_a, pair_b)
        assert comp.fraction_t_statistic is None


# ---------------------------------------------------------------------------
# AggregatedResultClass and _deserialize_result
# ---------------------------------------------------------------------------


class TestDeserializeResult:
    """Aggregated results load from canonical condition artifacts only."""

    def test_aggregated_result_class_set(self):
        from polyzymd.analyses.distances import DistancesAnalysis

        assert DistancesAnalysis.AggregatedResultClass is None

    def test_rejects_legacy_aggregate_json(self, tmp_path):
        from polyzymd.analyses.distances import DistancesAnalysis
        from polyzymd.analyses.exceptions import PluginContractError

        analysis = DistancesAnalysis()
        legacy_path = tmp_path / "distances_legacy.json"
        legacy_path.write_text('{"analysis_type": "distances_aggregated"}', encoding="utf-8")

        with pytest.raises(PluginContractError, match="canonical MDAnalysis condition artifact"):
            analysis._deserialize_result(legacy_path)


class TestSettingsCacheTag:
    """Settings fingerprinting should differentiate cache identities."""

    def test_cache_tag_changes_with_use_pbc(self):
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        settings_true = DistancesSettings(
            use_pbc=True,
            pairs=[DistancePairSettings(label="A", selection_a="resid 1", selection_b="resid 2")],
        )
        settings_false = DistancesSettings(
            use_pbc=False,
            pairs=[DistancePairSettings(label="A", selection_a="resid 1", selection_b="resid 2")],
        )

        tag_true = DistancesAnalysis._make_settings_cache_tag(settings_true)
        tag_false = DistancesAnalysis._make_settings_cache_tag(settings_false)

        assert tag_true != tag_false

    def test_cache_tag_changes_with_pair_definitions(self):
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        settings_a = DistancesSettings(
            pairs=[DistancePairSettings(label="A", selection_a="resid 1", selection_b="resid 2")],
        )
        settings_b = DistancesSettings(
            pairs=[DistancePairSettings(label="B", selection_a="resid 10", selection_b="resid 20")],
        )

        tag_a = DistancesAnalysis._make_settings_cache_tag(settings_a)
        tag_b = DistancesAnalysis._make_settings_cache_tag(settings_b)

        assert tag_a != tag_b


class TestDistanceCalculatorCacheIdentity:
    """DistanceCalculator should enforce strict cache identity."""

    def test_default_cache_filename_changes_when_only_later_pair_changes(self):
        from polyzymd.analyses.distances import DistanceCalculator

        config = MagicMock()
        with (
            patch("polyzymd.analyses.shared.loader._require_mdanalysis"),
            patch("polyzymd.analyses.distances.TrajectoryLoader"),
        ):
            calc_a = DistanceCalculator(
                config,
                pairs=[("resid 1 and name CA", "resid 2 and name CA"), ("resid 3", "resid 4")],
            )
            calc_b = DistanceCalculator(
                config,
                pairs=[("resid 1 and name CA", "resid 2 and name CA"), ("resid 30", "resid 40")],
            )

        assert calc_a._settings_tag != "legacy"
        assert calc_a._make_result_filename() != calc_b._make_result_filename()

    def test_cached_result_rejected_when_config_hash_validation_fails(self, tmp_path):
        from polyzymd.analyses.distances import DistanceCalculator

        config = MagicMock()
        pairs = [("resid 1 and name CA", "resid 2 and name CA")]
        thresholds = [3.5]

        with (
            patch("polyzymd.analyses.shared.loader._require_mdanalysis"),
            patch("polyzymd.analyses.distances.TrajectoryLoader"),
        ):
            calc = DistanceCalculator(
                config,
                pairs=pairs,
                thresholds=thresholds,
            )

        result_file = tmp_path / calc._make_result_filename()
        cached_result = _make_distance_cache_result(
            config_hash="stale-config-hash",
            pairs=pairs,
            thresholds=thresholds,
        )
        cached_result.save(result_file)
        calc._write_cache_metadata(result_file)

        with patch(
            "polyzymd.analyses._framework.cache_identity.validate_config_hash",
            return_value=False,
        ) as mock_validate:
            reused = calc._load_cached_result(result_file)

        assert reused is None
        mock_validate.assert_called_once_with("stale-config-hash", config)

    def test_stale_legacy_named_cache_is_rejected_when_settings_change(self, tmp_path):
        from polyzymd.analyses.distances import DistanceCalculator

        config = MagicMock()
        original_pairs = [
            ("resid 1 and name CA", "resid 2 and name CA"),
            ("resid 3 and name CA", "resid 4 and name CA"),
        ]
        changed_pairs = [
            ("resid 1 and name CA", "resid 2 and name CA"),
            ("resid 30 and name CA", "resid 40 and name CA"),
        ]
        thresholds = [3.5, 4.0]

        with (
            patch("polyzymd.analyses.shared.loader._require_mdanalysis"),
            patch("polyzymd.analyses.distances.TrajectoryLoader"),
        ):
            calc_old = DistanceCalculator(
                config,
                pairs=original_pairs,
                thresholds=thresholds,
                settings_tag="legacy",
            )
            calc_new = DistanceCalculator(
                config,
                pairs=changed_pairs,
                thresholds=thresholds,
                settings_tag="legacy",
            )

        result_file = tmp_path / calc_old._make_result_filename()
        assert result_file.name == calc_new._make_result_filename()

        cached_result = _make_distance_cache_result(
            config_hash=calc_old._config_hash,
            pairs=original_pairs,
            thresholds=thresholds,
        )
        cached_result.save(result_file)
        calc_old._write_cache_metadata(result_file)

        with patch(
            "polyzymd.analyses._framework.cache_identity.validate_config_hash",
            return_value=True,
        ):
            reused = calc_new._load_cached_result(result_file)

        assert reused is None


# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------


class TestPlot:
    """plot() delegates to all three inlined plotting functions."""

    def test_delegates_to_three_plotters(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.distances import DistancesAnalysis
        from polyzymd.config.comparison import PlotSettings

        analysis = DistancesAnalysis()

        cond = Condition(
            label="Cond1",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2, 3),
            sim_config=MagicMock(),
        )

        ctx = PlotContext(
            conditions=[cond],
            analysis_dirs={"Cond1": tmp_path / "analysis"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=MagicMock(),
            plot_settings=PlotSettings(),
        )

        with (
            patch("polyzymd.analyses.distances._plot_distance_kde") as MockKDE,
            patch("polyzymd.analyses.distances._plot_distance_threshold_bars") as MockThreshold,
            patch("polyzymd.analyses.distances._plot_distance_state_bars") as MockState,
        ):
            MockKDE.return_value = [Path("kde.png")]
            MockThreshold.return_value = [Path("bars.png")]
            MockState.return_value = [Path("state.png")]

            plots = analysis.plot(ctx)

        assert len(plots) == 3
        MockKDE.assert_called_once()
        MockThreshold.assert_called_once()
        MockState.assert_called_once()

    def test_plot_handles_failure(self, tmp_path):
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.distances import DistancesAnalysis
        from polyzymd.config.comparison import PlotSettings

        analysis = DistancesAnalysis()

        cond = Condition(
            label="Cond1",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=MagicMock(),
        )

        ctx = PlotContext(
            conditions=[cond],
            analysis_dirs={"Cond1": tmp_path / "analysis"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=MagicMock(),
            plot_settings=PlotSettings(),
        )

        with (
            patch(
                "polyzymd.analyses.distances._plot_distance_kde",
            ) as MockKDE,
            patch("polyzymd.analyses.distances._plot_distance_threshold_bars") as MockThreshold,
            patch("polyzymd.analyses.distances._plot_distance_state_bars") as MockState,
        ):
            MockKDE.side_effect = RuntimeError("KDE failed")
            MockThreshold.return_value = [Path("bars.png")]
            MockState.return_value = [Path("state.png")]

            plots = analysis.plot(ctx)

        # KDE failed, but threshold and state should still produce
        assert len(plots) == 2

    def test_plot_empty_conditions(self, tmp_path):
        from polyzymd.analyses.base import PlotContext
        from polyzymd.analyses.distances import DistancesAnalysis
        from polyzymd.config.comparison import PlotSettings

        analysis = DistancesAnalysis()

        ctx = PlotContext(
            conditions=[],
            analysis_dirs={},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=MagicMock(),
            plot_settings=PlotSettings(),
        )

        plots = analysis.plot(ctx)
        assert plots == []

    def test_plot_creates_default_plot_settings(self, tmp_path):
        """Orchestrator guarantees non-None plot_settings; plugin passes it through."""
        from polyzymd.analyses.base import Condition, PlotContext
        from polyzymd.analyses.distances import DistancesAnalysis
        from polyzymd.config.comparison import PlotSettings

        analysis = DistancesAnalysis()

        cond = Condition(
            label="Cond1",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1,),
            sim_config=MagicMock(),
        )

        plot_settings = PlotSettings()
        ctx = PlotContext(
            conditions=[cond],
            analysis_dirs={"Cond1": tmp_path / "analysis"},
            results_dir=tmp_path / "results",
            output_dir=tmp_path / "plots",
            settings=MagicMock(),
            plot_settings=plot_settings,
        )

        with (
            patch("polyzymd.analyses.distances._plot_distance_kde") as MockKDE,
            patch("polyzymd.analyses.distances._plot_distance_threshold_bars") as MockThreshold,
            patch("polyzymd.analyses.distances._plot_distance_state_bars") as MockState,
        ):
            MockKDE.return_value = []
            MockThreshold.return_value = []
            MockState.return_value = []

            analysis.plot(ctx)

        # Verify PlotSettings was passed through
        MockKDE.assert_called_once()
        call_args = MockKDE.call_args
        assert isinstance(call_args[0][3], PlotSettings)

    def test_semantic_colors_apply_to_condition_plots_but_not_state_categories(self, tmp_path):
        """Distance condition plots should use semantic colors while state bars stay categorical."""

        import numpy as np

        from polyzymd.analyses.distances._plotters import (
            DistancePlotData,
            _plot_distance_kde,
            _plot_distance_state_bars,
            _plot_distance_threshold_bars,
        )
        from polyzymd.config.comparison import PlotSettings

        plot_settings = PlotSettings(
            semantic_colors={
                "enabled": True,
                "order": ["Control", "Treatment"],
                "conditions": {
                    "Control": {"role": "control"},
                    "Treatment": {"color": "#ff7f0e"},
                },
                "control_color": "#111111",
            }
        )
        pair_result = {
            "pair_label": "Pair 1",
            "overall_fraction_below": 0.25,
            "sem_fraction_below": 0.05,
            "per_replicate_fractions_below": [0.2, 0.3],
            "threshold": 3.5,
        }
        plot_data = DistancePlotData(
            pooled_distances={
                "Pair 1": {
                    "Treatment": {"distances": np.asarray([3.0, 3.2, 3.4]), "threshold": 3.5},
                    "Control": {"distances": np.asarray([2.7, 2.9, 3.1]), "threshold": 3.5},
                }
            },
            aggregated_results={
                "Treatment": {"pair_results": [pair_result]},
                "Control": {"pair_results": [pair_result]},
            },
            control_label="Control",
        )

        with (
            patch(
                "polyzymd.analyses.distances._plotters.get_condition_colors",
                return_value=["#111111", "#ff7f0e"],
            ) as colors,
            patch(
                "polyzymd.analyses.distances._plotters.save_figure",
                side_effect=lambda fig, path, _: path,
            ),
            patch("polyzymd.analyses.distances._plotters.grouped_bars") as grouped,
        ):
            kde_paths = _plot_distance_kde(
                plot_data,
                ["Treatment", "Control"],
                tmp_path,
                plot_settings,
            )
            threshold_paths = _plot_distance_threshold_bars(
                plot_data,
                ["Treatment", "Control"],
                tmp_path,
                plot_settings,
            )
            state_paths = _plot_distance_state_bars(
                plot_data,
                ["Treatment", "Control"],
                tmp_path,
                plot_settings,
            )

        assert kde_paths == [tmp_path / "distance_kde_pair_1.png"]
        assert threshold_paths == [tmp_path / "distance_threshold_bars.png"]
        assert state_paths == [tmp_path / "distance_state_pair_1.png"]
        assert colors.call_args_list[0].args[0] == ["Control", "Treatment"]
        assert grouped.call_args_list[0].args[2][0][0] == "Control"
        assert grouped.call_args_list[0].args[3] == ["#111111", "#ff7f0e"]
        assert grouped.call_args_list[1].args[2][0][0].startswith("Below")
        assert grouped.call_args_list[1].args[3] != ["#111111", "#ff7f0e"]


# ---------------------------------------------------------------------------
# Plot loader cache filtering
# ---------------------------------------------------------------------------


class TestPlotLoadersCacheFiltering:
    """Plot data loaders should only read distances cache files."""

    def test_pooled_loader_ignores_unrelated_json_files(self, tmp_path):
        import numpy as np

        from polyzymd.analyses.distances._plotters import _load_pooled_distances
        from polyzymd.analyses.mda import ArtifactStore, ReplicateArtifact

        analysis_dir = tmp_path / "condition" / "distances"
        run_dir = analysis_dir / "run_1"
        run_dir.mkdir(parents=True)
        sidecar = ArtifactStore(run_dir).write_npz_sidecar(
            "sidecars/distances.npz",
            distance_matrix=np.asarray([[3.0, 3.2, 3.4]], dtype=np.float64),
            metadata={"kind": "distance_matrix"},
        )

        ArtifactStore(run_dir).write_replicate_result(
            ReplicateArtifact(
                analysis_name="distances",
                condition_label="Cond1",
                replicate=1,
                payload={"pairs": [{"pair_label": "Distance Pair", "threshold": 3.5}]},
                sidecars=[sidecar],
            )
        )
        (run_dir / "contacts_result.json").write_text(
            '{"pair_results": [{"pair_label": "Contacts Pair"}]}', encoding="utf-8"
        )
        (run_dir / "notes.json").write_text('{"note": "ignore me"}', encoding="utf-8")

        pooled = _load_pooled_distances(analysis_dir, [1])

        assert set(pooled.keys()) == {"Distance Pair"}
        assert pooled["Distance Pair"]["threshold"] == pytest.approx(3.5)

    def test_aggregated_loader_ignores_unrelated_json_files(self, tmp_path):
        from polyzymd.analyses.distances._plotters import _load_distance_aggregated_results
        from polyzymd.analyses.mda import ArtifactStore, ConditionArtifact

        aggregated_dir = tmp_path / "condition" / "distances" / "aggregated"
        aggregated_dir.mkdir(parents=True)

        ArtifactStore(aggregated_dir).write_condition_result(
            ConditionArtifact(
                analysis_name="distances",
                condition_label="Cond1",
                replicates=[1],
                payload={"pair_results": [{"pair_label": "Distance Pair", "overall_mean": 3.3}]},
            )
        )
        (aggregated_dir / "contacts_result.json").write_text(
            '{"pair_results": [{"pair_label": "Contacts Pair", "overall_mean": 9.9}]}',
            encoding="utf-8",
        )
        (aggregated_dir / "notes.json").write_text('{"note": "ignore me"}', encoding="utf-8")

        data = {"Cond1": {"aggregated_dir": aggregated_dir}}
        loaded = _load_distance_aggregated_results(data, ["Cond1"])

        assert "Cond1" in loaded
        pair_results = loaded["Cond1"]["pair_results"]
        assert len(pair_results) == 1
        assert pair_results[0]["pair_label"] == "Distance Pair"


# ---------------------------------------------------------------------------
# filter_conditions (default keeps all)
# ---------------------------------------------------------------------------


class TestFilterConditions:
    """Default filter_conditions keeps all conditions."""

    def test_keeps_all_conditions(self):
        from polyzymd.analyses.base import Condition
        from polyzymd.analyses.distances import DistancesAnalysis

        analysis = DistancesAnalysis()

        conditions = [
            Condition(
                label="A",
                config_path=Path("/tmp/a.yaml"),
                replicates=(1,),
                sim_config=MagicMock(),
            ),
            Condition(
                label="B",
                config_path=Path("/tmp/b.yaml"),
                replicates=(1,),
                sim_config=MagicMock(),
            ),
        ]

        result = analysis.filter_conditions(conditions)
        assert len(result) == 2


# ---------------------------------------------------------------------------
# Lifecycle test
# ---------------------------------------------------------------------------


class TestLifecycle:
    """End-to-end lifecycle: settings -> compute -> aggregate -> compare."""

    def test_settings_to_compute_to_aggregate_to_compare(self, tmp_path):
        """Verify the full lifecycle flow works end-to-end with mocks."""
        from polyzymd.analyses.base import (
            AggregateContext,
            Condition,
            ReplicateContext,
        )
        from polyzymd.analyses.distances import (
            DistancePairSettings,
            DistancesAnalysis,
            DistancesSettings,
        )

        analysis = DistancesAnalysis()
        settings = DistancesSettings(
            pairs=[DistancePairSettings(label="Pair1", selection_a="sel_a", selection_b="sel_b")]
        )

        mock_sim_config = MagicMock()
        cond = Condition(
            label="Test",
            config_path=Path("/tmp/config.yaml"),
            replicates=(1, 2),
            sim_config=mock_sim_config,
        )

        # 1. Compute-stage output is now canonical MDA replicate artifacts.
        rep_results = _make_distance_artifacts(tmp_path, "Test", settings, n_reps=2)
        assert len(rep_results) == 2

        # 2. Aggregate
        agg_ctx = AggregateContext(
            condition=cond,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=settings,
        )

        with patch(
            "polyzymd.analyses._framework.results_base.get_polyzymd_version",
            return_value="1.0.0-test",
        ):
            agg_result = analysis.aggregate(agg_ctx, rep_results)

        assert agg_result is not None
        assert agg_result.metadata["n_replicates"] == 2

    def test_extract_metrics_not_needed(self):
        """Since compare() is overridden, extract_metrics returns empty."""
        from polyzymd.analyses.distances import DistancesAnalysis

        analysis = DistancesAnalysis()
        # extract_metrics uses the default (returns {})
        assert analysis.extract_metrics(MagicMock()) == {}
