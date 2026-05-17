"""Tests for the MDAnalysis-native catalytic-triad plugin."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from polyzymd.analyses._framework.cache_identity import settings_fingerprint
from polyzymd.analyses.base import (
    AggregateContext,
    ComparisonContext,
    Condition,
    MetricValue,
    PlotContext,
    PluginContractError,
    ReplicateContext,
)
from polyzymd.analyses.catalytic_triad import (
    CatalyticTriadAnalysis,
    CatalyticTriadSettings,
    TriadPairSettings,
)
from polyzymd.analyses.catalytic_triad._mda import (
    SIMULTANEOUS_CONTACT_METRIC,
    TriadArtifactCollector,
    artifact_to_triad_result,
    condition_artifact_to_legacy_result,
)
from polyzymd.analyses.catalytic_triad._results import TriadAggregatedResult
from polyzymd.analyses.distances._mda import DistancePairPayload, DistanceReplicatePayload
from polyzymd.analyses.mda import (
    ArtifactStore,
    ComparisonArtifact,
    FrameSelection,
    MDACollectorContext,
    MDAJobResult,
    MDAUniversePolicy,
    ReplicateArtifact,
)
from polyzymd.analyses.mda.job import MDABackendPolicy


@pytest.fixture
def triad_analysis() -> CatalyticTriadAnalysis:
    """Return a fresh catalytic-triad analysis."""

    return CatalyticTriadAnalysis()


@pytest.fixture
def default_settings() -> CatalyticTriadSettings:
    """Return default catalytic-triad settings with two pairs."""

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
def condition() -> Condition:
    """Return a representative comparison condition."""

    return Condition(
        label="No Polymer",
        config_path=Path("/fake/config.yaml"),
        replicates=(1, 2, 3),
        sim_config=MagicMock(),
    )


def _collector_context(
    tmp_path: Path,
    condition: Condition,
    settings: CatalyticTriadSettings,
    *,
    replicate: int = 1,
) -> MDACollectorContext:
    """Build an MDA collector context for synthetic triad jobs."""

    from polyzymd.analyses.base import ReplicateContext

    output_dir = tmp_path / f"run_{replicate}"
    return MDACollectorContext(
        analysis_name="catalytic_triad",
        replicate_context=ReplicateContext(
            condition=condition,
            replicate=replicate,
            sim_config=condition.sim_config,
            output_dir=output_dir,
            equilibration="10ns",
            recompute=True,
            settings=settings,
        ),
        frame_selection=FrameSelection(
            start=0,
            stop=4,
            step=1,
            n_frames_total=4,
            timestep_ps=10.0,
        ),
        universe_policy=MDAUniversePolicy(condition_label=condition.label, replicate=replicate),
        artifact_store=ArtifactStore(output_dir),
        settings_fingerprint=settings_fingerprint(settings),
    )


def _completed_pair_distance_job(distance_matrix: np.ndarray) -> MDAJobResult:
    """Return a completed job with pair-distance analysis results."""

    results = SimpleNamespace(
        distance_matrix=np.asarray(distance_matrix, dtype=np.float64),
        frames=np.arange(distance_matrix.shape[1], dtype=np.int64),
        times_ps=np.arange(distance_matrix.shape[1], dtype=np.float64) * 10.0,
        warnings=[],
    )
    analysis = SimpleNamespace(results=results, frames=results.frames, times=results.times_ps)
    return MDAJobResult(
        name="triad_pair_distances",
        analysis=analysis,
        results=results,
        run_kwargs={},
        frame_selection=FrameSelection(start=0, stop=distance_matrix.shape[1], step=1),
        backend_policy=MDABackendPolicy(),
        universe_policy=MDAUniversePolicy(condition_label="No Polymer", replicate=1),
    )


def _replicate_artifact(
    tmp_path: Path,
    condition: Condition,
    settings: CatalyticTriadSettings,
    distance_matrix: np.ndarray,
    *,
    replicate: int,
) -> ReplicateArtifact:
    """Collect and persist one synthetic triad replicate artifact."""

    ctx = _collector_context(tmp_path, condition, settings, replicate=replicate)
    artifact = TriadArtifactCollector()(ctx, [_completed_pair_distance_job(distance_matrix)])
    ArtifactStore(ctx.output_dir).write_replicate_result(artifact)
    return artifact


class TestTriadDiscoveryAndSettings:
    """Tests for discovery, class variables, and settings."""

    def test_discovery_finds_catalytic_triad(self) -> None:
        """Catalytic triad should be auto-discovered by plugin discovery."""

        from polyzymd.analyses.discovery import clear_cache, get_analysis, list_analyses

        clear_cache()
        assert list_analyses()["catalytic_triad"] is CatalyticTriadAnalysis
        assert get_analysis("triad") is CatalyticTriadAnalysis

    def test_class_vars_use_mda_replicate_artifacts(
        self, triad_analysis: CatalyticTriadAnalysis
    ) -> None:
        """The migrated plugin should expose MDA jobs and no legacy replicate model."""

        assert triad_analysis.name == "catalytic_triad"
        assert triad_analysis.aliases == ("triad",)
        assert triad_analysis.ReplicateResultClass is None
        assert triad_analysis.AggregatedResultClass is TriadAggregatedResult

    def test_extract_metrics_metadata_preserves_primary_metric(
        self, triad_analysis: CatalyticTriadAnalysis
    ) -> None:
        """Metric extraction should preserve the primary comparison metadata."""

        summary = MagicMock(
            overall_simultaneous_contact=0.72,
            sem_simultaneous_contact=0.04,
            per_replicate_simultaneous=[0.70, 0.72, 0.74],
        )
        metric = triad_analysis.extract_metrics(summary)[SIMULTANEOUS_CONTACT_METRIC]

        assert metric.name == SIMULTANEOUS_CONTACT_METRIC
        assert metric.higher_is_better is True
        assert metric.direction_labels == ("worsening", "unchanged", "improving")

    def test_settings_validation(self, default_settings: CatalyticTriadSettings) -> None:
        """Settings should validate required pairs and expose labels."""

        assert default_settings.n_pairs == 2
        assert default_settings.get_pair_labels() == ["Asp133-His156", "His156-Ser77"]
        with pytest.raises(ValueError, match="At least one distance pair"):
            CatalyticTriadSettings(pairs=[])


class TestTriadMDACollector:
    """Tests for triad pair-distance artifact collection."""

    def test_collector_uses_strict_threshold_and_writes_sidecar(
        self, tmp_path: Path, condition: Condition, default_settings: CatalyticTriadSettings
    ) -> None:
        """Simultaneous contact should require every pair to be strictly below threshold."""

        matrix = np.asarray(
            [
                [3.0, 3.5, 3.49, 3.2],
                [3.1, 3.0, 3.6, 3.4],
            ],
            dtype=np.float64,
        )
        ctx = _collector_context(tmp_path, condition, default_settings)

        with patch(
            "polyzymd.analyses._framework.results_base.get_polyzymd_version",
            return_value="1.3.0",
        ):
            artifact = TriadArtifactCollector()(ctx, [_completed_pair_distance_job(matrix)])

        assert artifact.analysis_name == "catalytic_triad"
        assert artifact.payload["simultaneous_contact_fraction"] == pytest.approx(0.5)
        assert artifact.payload["n_frames_simultaneous"] == 2
        assert artifact.payload["metrics"][SIMULTANEOUS_CONTACT_METRIC] == pytest.approx(50.0)
        assert len(artifact.sidecars) == 1

        sidecar_path = ArtifactStore(ctx.output_dir).validate_sidecar(artifact.sidecars[0])
        with np.load(sidecar_path) as npz_data:
            np.testing.assert_allclose(npz_data["distance_matrix"], matrix)
            np.testing.assert_array_equal(
                npz_data["simultaneous_contact"], np.asarray([True, False, False, True])
            )
            assert npz_data["thresholds"].tolist() == [3.5, 3.5]
            assert npz_data["pair_labels"].tolist() == ["Asp133-His156", "His156-Ser77"]

    def test_artifact_to_legacy_result_preserves_pair_labels(
        self, tmp_path: Path, condition: Condition, default_settings: CatalyticTriadSettings
    ) -> None:
        """Artifact adapter should preserve established in-memory result fields."""

        artifact = _replicate_artifact(
            tmp_path,
            condition,
            default_settings,
            np.asarray([[3.0, 3.2], [3.1, 3.3]], dtype=np.float64),
            replicate=1,
        )

        result = artifact_to_triad_result(artifact)

        assert result.replicate == 1
        assert result.triad_name == "LipA_triad"
        assert result.get_pair_labels() == ["Asp133-His156", "His156-Ser77"]
        assert result.simultaneous_contact_fraction == 1.0
        assert result.settings_fingerprint == settings_fingerprint(default_settings)


class TestTriadAggregationAndComparison:
    """Tests for artifact aggregation and scalar comparison compatibility."""

    def test_mda_not_configured_fails_without_scalar_fallback(
        self,
        triad_analysis: CatalyticTriadAnalysis,
        condition: Condition,
        tmp_path: Path,
        default_settings: CatalyticTriadSettings,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """A missing MDA job path should fail closed before scalar fallback."""

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="10ns",
            recompute=True,
            settings=default_settings,
        )
        monkeypatch.setattr(
            "polyzymd.analyses.mda.lifecycle.run_mda_replicate_jobs",
            lambda *_args: None,
        )

        with pytest.raises(PluginContractError, match=r"build_mda_jobs\(\) returned None"):
            triad_analysis.run_replicate(ctx, replicate=1)

    def test_direct_scalar_runner_path_is_removed(
        self,
        triad_analysis: CatalyticTriadAnalysis,
    ) -> None:
        """The old runner hook should not exist on catalytic triad."""

        assert "build_runner" not in type(triad_analysis).__dict__
        assert not hasattr(triad_analysis, "_run_replicate_via_runner")

    def test_aggregate_requires_replicate_artifacts(
        self,
        triad_analysis: CatalyticTriadAnalysis,
        condition: Condition,
        tmp_path: Path,
        default_settings: CatalyticTriadSettings,
    ) -> None:
        """Legacy replicate results should be rejected with recompute guidance."""

        ctx = AggregateContext(
            condition=condition,
            replicates=(1,),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=default_settings,
        )

        with pytest.raises(TypeError, match="recompute the condition"):
            triad_analysis.aggregate(ctx, [object()])

    def test_aggregate_writes_condition_artifact_and_adapter(
        self,
        triad_analysis: CatalyticTriadAnalysis,
        condition: Condition,
        tmp_path: Path,
        default_settings: CatalyticTriadSettings,
    ) -> None:
        """Aggregation should produce canonical artifacts and legacy-compatible summaries."""

        artifacts = [
            _replicate_artifact(
                tmp_path,
                condition,
                default_settings,
                np.asarray([[3.0, 3.8, 3.1, 3.2], [3.1, 3.0, 3.7, 3.2]]),
                replicate=1,
            ),
            _replicate_artifact(
                tmp_path,
                condition,
                default_settings,
                np.asarray([[3.0, 3.2, 3.1, 3.2], [3.1, 3.0, 3.2, 3.2]]),
                replicate=2,
            ),
        ]
        two_replicate_condition = Condition(
            label=condition.label,
            config_path=condition.config_path,
            replicates=(1, 2),
            sim_config=condition.sim_config,
        )
        ctx = AggregateContext(
            condition=two_replicate_condition,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="10ns",
            settings=default_settings,
        )

        artifact = triad_analysis.aggregate(ctx, artifacts)
        legacy = condition_artifact_to_legacy_result(artifact)

        assert artifact.artifact_type == "condition"
        assert (tmp_path / "aggregated" / "result.json").exists()
        assert artifact.payload["metrics"][SIMULTANEOUS_CONTACT_METRIC]["values"] == [50.0, 100.0]
        assert artifact.payload["metric_metadata"][SIMULTANEOUS_CONTACT_METRIC] == {
            "label": "Simultaneous Contact",
            "unit": "%",
            "higher_is_better": True,
            "direction_labels": ("worsening", "unchanged", "improving"),
        }
        assert legacy.overall_simultaneous_contact == pytest.approx(0.75)
        assert legacy.per_replicate_simultaneous == [0.5, 1.0]
        assert legacy.get_pair_labels() == ["Asp133-His156", "His156-Ser77"]

    def test_compare_condition_artifacts_preserves_triad_metadata(
        self,
        triad_analysis: CatalyticTriadAnalysis,
        condition: Condition,
        tmp_path: Path,
        default_settings: CatalyticTriadSettings,
    ) -> None:
        """Real condition artifacts should compare through the MDA artifact engine."""

        control = Condition(
            label="No Polymer",
            config_path=Path("/fake/control.yaml"),
            replicates=(1, 2),
            sim_config=condition.sim_config,
        )
        peg = Condition(
            label="PEG",
            config_path=Path("/fake/peg.yaml"),
            replicates=(1, 2),
            sim_config=condition.sim_config,
        )
        control_base = tmp_path / "control" / "catalytic_triad"
        peg_base = tmp_path / "peg" / "catalytic_triad"
        control_artifact = triad_analysis.aggregate(
            AggregateContext(
                condition=control,
                replicates=(1, 2),
                output_dir=control_base / "aggregated",
                equilibration="10ns",
                settings=default_settings,
            ),
            [
                _replicate_artifact(
                    control_base,
                    control,
                    default_settings,
                    np.asarray([[3.0, 3.8, 3.1, 3.2], [3.1, 3.0, 3.7, 3.2]]),
                    replicate=1,
                ),
                _replicate_artifact(
                    control_base,
                    control,
                    default_settings,
                    np.asarray([[3.0, 3.2, 3.1, 3.8], [3.1, 3.0, 3.2, 3.2]]),
                    replicate=2,
                ),
            ],
        )
        peg_artifact = triad_analysis.aggregate(
            AggregateContext(
                condition=peg,
                replicates=(1, 2),
                output_dir=peg_base / "aggregated",
                equilibration="10ns",
                settings=default_settings,
            ),
            [
                _replicate_artifact(
                    peg_base,
                    peg,
                    default_settings,
                    np.asarray([[3.0, 3.2, 3.1, 3.8], [3.1, 3.0, 3.2, 3.2]]),
                    replicate=1,
                ),
                _replicate_artifact(
                    peg_base,
                    peg,
                    default_settings,
                    np.asarray([[3.0, 3.2, 3.1, 3.2], [3.1, 3.0, 3.2, 3.2]]),
                    replicate=2,
                ),
            ],
        )
        ctx = ComparisonContext(
            name="triad-project",
            conditions=[control, peg],
            excluded_conditions=[],
            control_label="No Polymer",
            analysis_dirs={"No Polymer": control_base, "PEG": peg_base},
            results_dir=tmp_path / "results",
            equilibration="10ns",
            settings=default_settings,
            fdr_alpha=0.05,
            ttest_method="welch",
            posthoc_method="ttest_bh",
            aggregated_results={"No Polymer": control_artifact, "PEG": peg_artifact},
        )

        comparison = triad_analysis.compare(ctx)

        assert isinstance(comparison, ComparisonArtifact)
        assert comparison.payload["ranking"] == ["PEG", "No Polymer"]
        metadata = comparison.payload["metric_metadata"][SIMULTANEOUS_CONTACT_METRIC]
        assert metadata["label"] == "Simultaneous Contact"
        assert metadata["unit"] == "%"
        assert metadata["higher_is_better"] is True
        assert metadata["direction_labels"] == ["worsening", "unchanged", "improving"]
        assert comparison.payload["pairwise_comparisons"][0]["direction"] == "improving"

    def test_format_accepts_comparison_artifact(
        self,
        triad_analysis: CatalyticTriadAnalysis,
    ) -> None:
        """ComparisonArtifact output should use triad-specific scalar formatting."""

        comparison = ComparisonArtifact(
            analysis_name="catalytic_triad",
            conditions=["No Polymer", "PEG"],
            control_label="No Polymer",
            effective_control="No Polymer",
            payload={
                "condition_summaries": [
                    {
                        "label": "No Polymer",
                        "n_replicates": 2,
                        "simultaneous_contact_fraction_mean": 62.5,
                        "simultaneous_contact_fraction_sem": 12.5,
                        "simultaneous_contact_fraction_replicate_values": [50.0, 75.0],
                    },
                    {
                        "label": "PEG",
                        "n_replicates": 2,
                        "simultaneous_contact_fraction_mean": 87.5,
                        "simultaneous_contact_fraction_sem": 12.5,
                        "simultaneous_contact_fraction_replicate_values": [75.0, 100.0],
                    },
                ],
                "pairwise_comparisons": [
                    {
                        "condition_a": "No Polymer",
                        "condition_b": "PEG",
                        "metric": SIMULTANEOUS_CONTACT_METRIC,
                        "t_statistic": 1.0,
                        "p_value": 0.5,
                        "p_value_adjusted": 0.5,
                        "posthoc_method": "ttest_bh",
                        "cohens_d": 1.0,
                        "effect_size_interpretation": "large",
                        "direction": "improving",
                        "significant": False,
                        "percent_change": 40.0,
                        "testable": True,
                    }
                ],
                "anova": [],
                "ranking": ["PEG", "No Polymer"],
                "rankings_by_metric": {SIMULTANEOUS_CONTACT_METRIC: ["PEG", "No Polymer"]},
                "metric_metadata": {
                    SIMULTANEOUS_CONTACT_METRIC: {
                        "label": "Simultaneous Contact",
                        "unit": "%",
                        "higher_is_better": True,
                        "direction_labels": ["worsening", "unchanged", "improving"],
                    }
                },
                "statistical_parameters": {
                    "fdr_alpha": 0.05,
                    "ttest_method": "welch",
                    "posthoc_method": "ttest_bh",
                    "control_label": "No Polymer",
                    "effective_control": "No Polymer",
                    "equilibration": "10ns",
                },
            },
        )

        formatted = triad_analysis.format(comparison, "text")
        json_formatted = triad_analysis.format(comparison, "json")

        assert "Catalytic Triad Comparison" in formatted
        assert "Simultaneous Contact" in formatted
        assert "PEG" in formatted
        assert "artifact_type" not in formatted
        assert '"artifact_type": "comparison"' in json_formatted
        assert '"metric_metadata"' in json_formatted

    def test_extract_metrics_preserves_percent_scaling(
        self, triad_analysis: CatalyticTriadAnalysis
    ) -> None:
        """Legacy metric extraction should still report percentages."""

        summary = MagicMock(
            overall_simultaneous_contact=0.72,
            sem_simultaneous_contact=0.04,
            per_replicate_simultaneous=[0.70, 0.72, 0.74],
        )

        metrics = triad_analysis.extract_metrics(summary)

        metric = metrics[SIMULTANEOUS_CONTACT_METRIC]
        assert isinstance(metric, MetricValue)
        assert metric.mean == 72.0
        assert metric.sem == 4.0
        assert metric.replicate_values == [70.0, 72.0, 74.0]

    def test_deserialize_rejects_legacy_json(
        self, triad_analysis: CatalyticTriadAnalysis, tmp_path: Path
    ) -> None:
        """Legacy triad JSON should not be loaded as a canonical aggregate."""

        legacy_path = tmp_path / "triad_legacy.json"
        legacy_path.write_text('{"analysis_type": "catalytic_triad_aggregated"}', encoding="utf-8")

        with pytest.raises(ValueError, match="canonical MDAnalysis condition artifact"):
            triad_analysis._deserialize_result(legacy_path)


class TestTriadPlot:
    """Tests for plotting helpers and artifact-only loading."""

    def test_threshold_bars_overlay_pair_and_simultaneous_replicates(self) -> None:
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
                [result], ["Control"], colors=["blue"], plot_settings=PlotSettings()
            )

        mock_scatter.assert_called_once()
        assert mock_scatter.call_args.args[2] == [[40.0, 60.0], [20.0, 30.0], [8.0, 12.0]]
        plt.close(fig)

    @patch("polyzymd.analyses.catalytic_triad.plot_triad_threshold_bars_from_data")
    @patch("polyzymd.analyses.catalytic_triad.plot_triad_kde_panel_from_data")
    def test_plot_delegates_to_artifact_helpers(
        self,
        mock_kde_fn: MagicMock,
        mock_bars_fn: MagicMock,
        triad_analysis: CatalyticTriadAnalysis,
        condition: Condition,
        tmp_path: Path,
        default_settings: CatalyticTriadSettings,
    ) -> None:
        """The public plot hook should delegate to artifact-backed helpers."""

        from polyzymd.config.comparison import PlotSettings

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

        plots = triad_analysis.plot(ctx)

        assert len(plots) == 2
        mock_kde_fn.assert_called_once()
        mock_bars_fn.assert_called_once()


class TestTriadLifecycle:
    """Tests for the MDA lifecycle hooks."""

    def test_build_mda_collector_returns_triad_collector(
        self,
        triad_analysis: CatalyticTriadAnalysis,
        condition: Condition,
        tmp_path: Path,
        default_settings: CatalyticTriadSettings,
    ) -> None:
        """The migrated plugin should expose a triad artifact collector."""

        ctx = _collector_context(tmp_path, condition, default_settings)
        assert isinstance(triad_analysis.build_mda_collector(ctx), TriadArtifactCollector)

    def test_build_mda_jobs_delegates_to_triad_job_builder(
        self, triad_analysis: CatalyticTriadAnalysis
    ) -> None:
        """The public hook should delegate to the triad job builder."""

        mda_ctx = MagicMock()
        mda_ctx.settings = MagicMock()
        sentinel_jobs = [object()]
        with patch(
            "polyzymd.analyses.catalytic_triad.build_triad_jobs", return_value=sentinel_jobs
        ):
            assert triad_analysis.build_mda_jobs(mda_ctx) is sentinel_jobs


class TestTriadDirectPayloadCompatibility:
    """Tests for reusable payload structures used by triad adapters."""

    def test_distance_payload_imports_from_mda_module(self) -> None:
        """Distance payload classes should no longer require the deleted runner wrapper."""

        payload = DistanceReplicatePayload(
            n_frames_total=2,
            n_frames_used=2,
            pair_payloads=[
                DistancePairPayload(
                    pair_label="pair",
                    selection1="a",
                    selection2="b",
                    distances=np.asarray([1.0, 2.0]),
                    mean_distance=1.5,
                    std_distance=0.5,
                    median_distance=1.5,
                    min_distance=1.0,
                    max_distance=2.0,
                    sem_distance=None,
                    correlation_time=None,
                    correlation_time_unit=None,
                    n_independent_frames=None,
                    statistical_inefficiency=None,
                    autocorrelation_warning=None,
                    threshold=3.5,
                    fraction_below_threshold=1.0,
                    histogram_edges=np.asarray([1.0, 2.0]),
                    histogram_counts=np.asarray([2]),
                    kde_x=None,
                    kde_y=None,
                    kde_peak=None,
                    kde_bandwidth=None,
                    n_frames_total=2,
                    n_frames_used=2,
                )
            ],
        )

        assert payload.n_frames_used == 2
