"""Tests for the RMSD analysis plugin."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from polyzymd.analyses._analysis_lifecycle import AnalysisLifecycle
from polyzymd.analyses.base import Condition, PlotContext
from polyzymd.analyses.discovery import get_analysis
from polyzymd.analyses.mda import (
    ArtifactStore,
    FrameSelection,
    MDAAnalysisJob,
    MDABackendPolicy,
    MDACollectorContext,
    MDAJobResult,
    MDAReplicateJobContext,
    MDAUniversePolicy,
    ReplicateArtifact,
)
from polyzymd.analyses.rmsd import RMSDAnalysis, RMSDRunSettings, RMSDSettings
from polyzymd.analyses.rmsd._comparison_results import (
    RMSDComparisonResult,
    RMSDConditionSummary,
    RMSDRunANOVA,
    RMSDRunPairwiseComparison,
    RMSDRunSummary,
)
from polyzymd.analyses.rmsd._formatters import format_rmsd_comparison
from polyzymd.analyses.rmsd._mda import (
    _build_rmsd_analysis,
    _sidecar_filename,
    condition_artifact_to_legacy_result,
)
from polyzymd.analyses.rmsd._plot_settings import RMSDPlotSettings
from polyzymd.analyses.rmsd._plotters import _resolve_npz_sidecar_path
from polyzymd.analyses.rmsd._results import (
    RMSDAggregatedResult,
    RMSDResult,
    RMSDRunAggregatedResult,
    RMSDRunResult,
)
from polyzymd.analyses.shared.config_hash import settings_fingerprint
from polyzymd.config.comparison import PlotSettings
from tests._support.analysis_testkit import (
    make_aggregate_context,
    make_comparison_context,
    make_condition,
    make_replicate_context,
)

EXPECTED_STRETCH_RMSD = [0.0, 4.0 / 3.0, 4.0]


def _make_run_settings() -> list[RMSDRunSettings]:
    """Create two valid RMSD run settings for tests."""
    return [
        RMSDRunSettings(label="protein_backbone"),
        RMSDRunSettings(label="polymer_core", selection="segid C and backbone"),
    ]


def _settings_hash(settings: RMSDSettings) -> str:
    """Return the shared settings fingerprint used by RMSD caches."""
    return settings_fingerprint(settings)


def _make_stretch_universe() -> object:
    """Create a synthetic triangle trajectory with hand-computable RMSD values."""

    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(
        3,
        n_residues=3,
        atom_resindex=np.arange(3),
        trajectory=True,
    )
    universe.add_TopologyAttr("names", ["A1", "A2", "A3"])
    universe.add_TopologyAttr("masses", [12.0, 12.0, 12.0])
    universe.add_TopologyAttr("resnames", ["TST", "TST", "TST"])
    universe.add_TopologyAttr("resids", [1, 2, 3])
    coords = np.asarray(
        [
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
            [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [0.0, 4.0, 0.0]],
            [[0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [0.0, 8.0, 0.0]],
        ],
        dtype=np.float32,
    )
    universe.load_new(coords, format=MemoryReader, dt=5.0)
    return universe


def _make_stretch_run_settings(label: str = "stretch_triangle") -> RMSDRunSettings:
    """Create RMSD settings for the synthetic triangle trajectory."""

    return RMSDRunSettings(
        label=label,
        selection="all",
        alignment_selection="all",
        centroid_selection="all",
        reference_mode="frame",
        reference_frame=0,
    )


def _make_run_result(
    replicate: int,
    run_label: str,
    mean_rmsd: float,
    convergence_time_ns: float | None = None,
    convergence_assessable: bool = True,
) -> RMSDRunResult:
    """Create a minimal, valid RMSDRunResult for test fixtures."""
    return RMSDRunResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=replicate,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA",
        run_label=run_label,
        selection="protein and name CA",
        alignment_selection="protein and name CA",
        reference_mode="centroid",
        reference_frame=1,
        mean_rmsd=mean_rmsd,
        std_rmsd=0.2,
        median_rmsd=mean_rmsd,
        min_rmsd=mean_rmsd - 0.3,
        max_rmsd=mean_rmsd + 0.3,
        final_rmsd=mean_rmsd + 0.1,
        sem_rmsd=0.1,
        converged=convergence_time_ns is not None,
        convergence_assessable=convergence_assessable,
        convergence_time_ns=convergence_time_ns,
        convergence_message="test convergence",
        n_frames_total=100,
        n_frames_used=90,
        npz_path=f"/fake/{run_label}_rep{replicate}.npz",
        time_unit="ns",
        timestep_ps=10.0,
    )


def _make_run_payload(
    replicate: int,
    run_label: str,
    mean_rmsd: float,
    convergence_time_ns: float | None = None,
    convergence_assessable: bool = True,
) -> dict[str, object]:
    """Create a run payload stored inside an RMSD replicate artifact."""

    return {
        "run_label": run_label,
        "selection": (
            "protein and name CA" if run_label == "protein_backbone" else "segid C and backbone"
        ),
        "alignment_selection": (
            "protein and name CA" if run_label == "protein_backbone" else "segid C and backbone"
        ),
        "reference_mode": "centroid",
        "reference_frame": 1,
        "reference_file": None,
        "mean_rmsd": mean_rmsd,
        "std_rmsd": 0.2,
        "median_rmsd": mean_rmsd,
        "min_rmsd": mean_rmsd - 0.3,
        "max_rmsd": mean_rmsd + 0.3,
        "final_rmsd": mean_rmsd + 0.1,
        "sem_rmsd": 0.1,
        "converged": convergence_time_ns is not None,
        "convergence_assessable": convergence_assessable,
        "convergence_time_ns": convergence_time_ns,
        "convergence_message": "test convergence",
        "n_frames_total": 100,
        "n_frames_used": 90,
        "npz_path": f"sidecars/{run_label}_rep{replicate}.npz",
        "time_unit": "ns",
        "timestep_ps": 10.0,
        "raw_timestep_ps": 10.0,
        "frame_stride": 1,
    }


def _make_replicate_artifact(
    replicate: int,
    run_payloads: list[dict[str, object]],
    *,
    condition_label: str = "Test Condition",
    settings_hash: str = "hash12345",
) -> ReplicateArtifact:
    """Create an RMSD replicate artifact for aggregation tests."""

    metrics = {
        f"{payload['run_label']}.mean_rmsd": float(payload["mean_rmsd"]) for payload in run_payloads
    }
    return ReplicateArtifact(
        analysis_name="rmsd",
        condition_label=condition_label,
        replicate=replicate,
        payload={
            "runs": run_payloads,
            "metrics": metrics,
            "replicate_metrics": metrics,
            "n_runs": len(run_payloads),
            "n_frames_total": 100,
            "n_frames_used": 90,
        },
        provenance={
            "frame_selection": {
                "start": 10,
                "stop": 100,
                "step": 1,
                "n_frames_total": 100,
                "n_frames_selected": 90,
            }
        },
        metadata={
            "settings_fingerprint": settings_hash,
            "config_hash": "hash123",
            "polyzymd_version": "1.2.1",
            "equilibration_time": 10.0,
            "equilibration_unit": "ns",
        },
    )


def _make_rmsd_collector_context(
    *,
    condition: Condition,
    tmp_path: Path,
    settings: RMSDSettings,
    frame_selection: FrameSelection,
    replicate: int = 1,
) -> tuple[MDACollectorContext, ArtifactStore]:
    """Create an RMSD collector context and artifact store for tests."""

    replicate_ctx = make_replicate_context(
        condition=condition,
        replicate=replicate,
        output_dir=tmp_path / f"run_{replicate}",
        settings=settings,
        equilibration="10ns",
    )
    store = ArtifactStore(replicate_ctx.output_dir)
    return (
        MDACollectorContext(
            analysis_name="rmsd",
            replicate_context=replicate_ctx,
            frame_selection=frame_selection,
            universe_policy=MDAUniversePolicy(condition_label=condition.label, replicate=replicate),
            artifact_store=store,
            settings_fingerprint=_settings_hash(settings),
        ),
        store,
    )


def _make_rmsd_job_result(
    *,
    condition: Condition,
    frame_selection: FrameSelection,
    run: RMSDRunSettings,
    rmsd_table: np.ndarray,
    analysis_obj: object | None = None,
) -> MDAJobResult:
    """Create a completed RMSD MDA job result for collector tests."""

    if analysis_obj is None:
        analysis_obj = SimpleNamespace(
            _polyzymd_rmsd_metadata={"reference_frame": 1, "reference_file": None}
        )

    return MDAJobResult(
        name=run.label,
        analysis=analysis_obj,
        results=SimpleNamespace(rmsd=np.asarray(rmsd_table, dtype=np.float64)),
        run_kwargs=frame_selection.run_kwargs(),
        frame_selection=frame_selection,
        backend_policy=MDABackendPolicy(),
        universe_policy=MDAUniversePolicy(
            condition_label=condition.label,
            replicate=1,
            metadata={"rmsd_run": run.model_dump(mode="json")},
        ),
    )


def _make_aggregated_run(
    run_label: str,
    selection: str,
    per_replicate_means: list[float],
) -> RMSDRunAggregatedResult:
    """Create an aggregated run result with deterministic values."""
    mean_value = sum(per_replicate_means) / len(per_replicate_means)
    return RMSDRunAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string=selection,
        replicates=[1, 2, 3][: len(per_replicate_means)],
        n_replicates=len(per_replicate_means),
        run_label=run_label,
        selection=selection,
        alignment_selection=selection,
        overall_mean=mean_value,
        overall_sem=0.05,
        overall_median=mean_value,
        per_replicate_means=per_replicate_means,
        per_replicate_stds=[0.2 for _ in per_replicate_means],
        per_replicate_medians=per_replicate_means,
        per_replicate_convergence_times_ns=[None for _ in per_replicate_means],
        per_replicate_convergence_assessable=[True for _ in per_replicate_means],
        n_converged_replicates=0,
        n_assessable_replicates=len(per_replicate_means),
        convergence_fraction=0.0,
        all_converged=False,
    )


def _make_condition_summary(
    label: str, backbone_mean: float, polymer_mean: float
) -> RMSDConditionSummary:
    """Create a condition summary with two run summaries."""
    return RMSDConditionSummary(
        label=label,
        config_path=f"/fake/{label}/config.yaml",
        n_replicates=3,
        run_summaries=[
            RMSDRunSummary(
                label="protein_backbone",
                selection="protein and name CA",
                mean_rmsd=backbone_mean,
                sem_rmsd=0.05,
                per_replicate_means=[backbone_mean - 0.05, backbone_mean, backbone_mean + 0.05],
                n_converged_replicates=2,
                n_assessable_replicates=3,
                convergence_fraction=2.0 / 3.0,
                all_converged=False,
                mean_convergence_time_ns=35.0,
                median_convergence_time_ns=35.0,
            ),
            RMSDRunSummary(
                label="polymer_core",
                selection="segid C and backbone",
                mean_rmsd=polymer_mean,
                sem_rmsd=0.04,
                per_replicate_means=[polymer_mean - 0.04, polymer_mean, polymer_mean + 0.04],
                n_converged_replicates=3,
                n_assessable_replicates=3,
                convergence_fraction=1.0,
                all_converged=True,
                mean_convergence_time_ns=30.0,
                median_convergence_time_ns=30.0,
            ),
        ],
    )


def _make_comparison_result() -> RMSDComparisonResult:
    """Create a complete RMSDComparisonResult for formatter/model tests."""
    conditions = [
        _make_condition_summary("Control", 1.20, 2.00),
        _make_condition_summary("Treatment_A", 1.05, 1.80),
        _make_condition_summary("Treatment_B", 1.35, 2.10),
    ]
    pairwise = [
        RMSDRunPairwiseComparison(
            run_label="protein_backbone",
            condition_a="Control",
            condition_b="Treatment_A",
            t_statistic=2.5,
            p_value=0.03,
            cohens_d=1.1,
            effect_interpretation="large",
            direction="stabilizing",
            significant=True,
            percent_change=-12.5,
        ),
        RMSDRunPairwiseComparison(
            run_label="protein_backbone",
            condition_a="Control",
            condition_b="Treatment_B",
            t_statistic=-2.1,
            p_value=0.04,
            cohens_d=-0.9,
            effect_interpretation="large",
            direction="destabilizing",
            significant=True,
            percent_change=12.5,
        ),
        RMSDRunPairwiseComparison(
            run_label="polymer_core",
            condition_a="Control",
            condition_b="Treatment_A",
            t_statistic=3.0,
            p_value=0.02,
            cohens_d=1.4,
            effect_interpretation="large",
            direction="stabilizing",
            significant=True,
            percent_change=-10.0,
        ),
    ]
    anova = [
        RMSDRunANOVA(run_label="protein_backbone", f_statistic=8.5, p_value=0.01, significant=True),
        RMSDRunANOVA(run_label="polymer_core", f_statistic=9.2, p_value=0.01, significant=True),
    ]
    return RMSDComparisonResult(
        metric="mean_rmsd",
        name="test_project",
        n_runs=2,
        run_labels=["protein_backbone", "polymer_core"],
        control_label="Control",
        conditions=conditions,
        pairwise_comparisons=pairwise,
        anova_by_run=anova,
        ranking_by_run={
            "protein_backbone": ["Treatment_A", "Control", "Treatment_B"],
            "polymer_core": ["Treatment_A", "Control", "Treatment_B"],
        },
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.2.1",
    )


def test_rmsd_run_settings_defaults() -> None:
    """RMSDRunSettings should provide expected defaults."""
    settings = RMSDRunSettings(label="protein_backbone")
    assert settings.selection == "protein and name CA"
    assert settings.alignment_selection == "protein and name CA"
    assert settings.reference_mode == "centroid"
    assert settings.convergence_window_size_ns == pytest.approx(15.0)
    assert settings.convergence_step_size_ns == pytest.approx(5.0)
    assert settings.convergence_slope_threshold == pytest.approx(0.0005)
    assert settings.convergence_sustained_for_ns == pytest.approx(15.0)


def test_rmsd_run_settings_reject_invalid_convergence_step() -> None:
    """RMSDRunSettings should reject convergence step larger than window."""
    with pytest.raises(ValueError, match="convergence_step_size_ns"):
        RMSDRunSettings(
            label="protein_backbone",
            convergence_window_size_ns=10.0,
            convergence_step_size_ns=15.0,
        )


def test_rmsd_settings_valid() -> None:
    """RMSDSettings should accept a valid multi-run configuration."""
    settings = RMSDSettings(runs=_make_run_settings())
    assert len(settings.runs) == 2


def test_rmsd_settings_empty_runs_rejected() -> None:
    """RMSDSettings should reject an empty runs list."""
    with pytest.raises(ValueError, match="At least one RMSD run"):
        RMSDSettings(runs=[])


def test_rmsd_settings_duplicate_labels_rejected() -> None:
    """RMSDSettings should reject duplicate run labels."""
    with pytest.raises(ValueError, match="labels must be unique"):
        RMSDSettings(
            runs=[
                RMSDRunSettings(label="dup"),
                RMSDRunSettings(label="dup", selection="backbone"),
            ]
        )


def test_rmsd_settings_single_run() -> None:
    """RMSDSettings should allow a single-run configuration."""
    settings = RMSDSettings(runs=[RMSDRunSettings(label="single")])
    assert len(settings.runs) == 1
    assert settings.runs[0].label == "single"


def test_settings_cache_tag_changes_with_settings() -> None:
    """Settings cache tag should change when run settings change."""
    analysis = RMSDAnalysis()
    settings_a = RMSDSettings(
        runs=[RMSDRunSettings(label="run_a", selection="protein and name CA")]
    )
    settings_b = RMSDSettings(runs=[RMSDRunSettings(label="run_a", selection="protein")])

    tag_a = analysis._make_settings_cache_tag(settings_a)
    tag_b = analysis._make_settings_cache_tag(settings_b)

    assert tag_a != tag_b


def test_build_mda_jobs_uses_one_job_per_run(condition: Condition, tmp_path: Path) -> None:
    """RMSD should expose MDAnalysis-native jobs instead of legacy runners."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    replicate_ctx = make_replicate_context(
        condition=condition,
        replicate=1,
        output_dir=tmp_path / "run_1",
        settings=settings,
        equilibration="10ns",
    )
    frame_selection = FrameSelection(
        start=10, stop=100, step=2, n_frames_total=120, timestep_ps=10.0
    )
    mda_ctx = MDAReplicateJobContext(
        replicate_context=replicate_ctx,
        universe=object(),
        frame_selection=frame_selection,
        universe_policy=MDAUniversePolicy(condition_label=condition.label, replicate=1),
        artifact_store=ArtifactStore(tmp_path / "run_1"),
    )

    jobs = analysis.build_mda_jobs(mda_ctx)

    assert jobs is not None
    assert [job.name for job in jobs] == ["protein_backbone", "polymer_core"]
    assert all(job.analysis_factory is not None for job in jobs)
    assert all(job.universe_policy.metadata["fresh_universe_per_run"] is True for job in jobs)


def test_mda_collector_writes_replicate_artifact_sidecar(
    condition: Condition, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """RMSD collector should preserve MDAnalysis time values in NPZ sidecars."""
    monkeypatch.setattr("polyzymd.analyses.rmsd._mda.compute_config_hash", lambda _cfg: "hash123")
    monkeypatch.setattr(
        "polyzymd.analyses.shared.convergence.find_convergence_time",
        lambda *args, **kwargs: SimpleNamespace(
            window_start_times_ns=[0.0],
            window_mean_values=[1.1],
            slope_times_ns=[1.0],
            slopes=[0.0],
            converged=True,
            assessable=True,
            convergence_time_ns=1.0,
            message="Converged at 1 ns",
        ),
    )
    settings = RMSDSettings(runs=[RMSDRunSettings(label="protein_backbone")])
    frame_selection = FrameSelection(
        start=1000, stop=1003, step=1, n_frames_total=2000, timestep_ps=10.0
    )
    collector_ctx, store = _make_rmsd_collector_context(
        condition=condition,
        tmp_path=tmp_path,
        settings=settings,
        frame_selection=frame_selection,
    )
    job = _make_rmsd_job_result(
        condition=condition,
        frame_selection=frame_selection,
        run=settings.runs[0],
        rmsd_table=np.asarray(
            [[1000.0, 5.0, 1.0], [1001.0, 25.0, 1.2], [1002.0, 50.0, 1.4]],
            dtype=np.float64,
        ),
    )

    artifact = RMSDAnalysis().build_mda_collector(collector_ctx)(collector_ctx, [job])

    assert artifact.payload["metrics"]["protein_backbone.mean_rmsd"] == pytest.approx(1.2)
    assert artifact.metadata["settings_fingerprint"] == _settings_hash(settings)
    assert artifact.sidecars[0].path == (f"sidecars/{_sidecar_filename('protein_backbone', 0)}")
    with np.load(store.resolve_sidecar(artifact.sidecars[0])) as payload:
        assert payload["rmsd_values"].tolist() == pytest.approx([1.0, 1.2, 1.4])
        assert payload["time_ns"].tolist() == pytest.approx([0.005, 0.025, 0.05])
        assert float(payload["raw_timestep_ps"]) == pytest.approx(10.0)
        assert float(payload["effective_timestep_ps"]) == pytest.approx(22.5)


def test_mda_collector_uses_collision_proof_sidecar_names(
    condition: Condition, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """RMSD sidecar names should not collide after label sanitization."""
    monkeypatch.setattr("polyzymd.analyses.rmsd._mda.compute_config_hash", lambda _cfg: "hash123")
    monkeypatch.setattr(
        "polyzymd.analyses.shared.convergence.find_convergence_time",
        lambda *args, **kwargs: SimpleNamespace(
            window_start_times_ns=[],
            window_mean_values=[],
            slope_times_ns=[],
            slopes=[],
            converged=False,
            assessable=False,
            convergence_time_ns=None,
            message="Too few frames",
        ),
    )
    settings = RMSDSettings(
        runs=[
            RMSDRunSettings(label="a-b"),
            RMSDRunSettings(label="a_b"),
        ]
    )
    frame_selection = FrameSelection(start=0, stop=2, step=1, n_frames_total=2, timestep_ps=5.0)
    collector_ctx, store = _make_rmsd_collector_context(
        condition=condition,
        tmp_path=tmp_path,
        settings=settings,
        frame_selection=frame_selection,
    )
    jobs = [
        _make_rmsd_job_result(
            condition=condition,
            frame_selection=frame_selection,
            run=run,
            rmsd_table=np.asarray([[0.0, 0.0, 0.0], [1.0, 5.0, 1.0]], dtype=np.float64),
        )
        for run in settings.runs
    ]

    artifact = RMSDAnalysis().build_mda_collector(collector_ctx)(collector_ctx, jobs)

    sidecar_paths = [sidecar.path for sidecar in artifact.sidecars]
    assert sidecar_paths == [
        f"sidecars/{_sidecar_filename('a-b', 0)}",
        f"sidecars/{_sidecar_filename('a_b', 1)}",
    ]
    assert len(set(sidecar_paths)) == 2
    for sidecar in artifact.sidecars:
        assert store.resolve_sidecar(sidecar).exists()


def test_mda_rmsd_job_synthetic_universe_matches_expected_artifacts(
    condition: Condition, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Real MDAnalysis RMSD jobs should feed aggregation and comparison artifacts."""
    condition = make_condition(label=condition.label, replicates=(1,))
    monkeypatch.setattr("polyzymd.analyses.rmsd._mda.compute_config_hash", lambda _cfg: "hash123")
    monkeypatch.setattr(
        "polyzymd.analyses.shared.convergence.find_convergence_time",
        lambda *args, **kwargs: SimpleNamespace(
            window_start_times_ns=[],
            window_mean_values=[],
            slope_times_ns=[],
            slopes=[],
            converged=False,
            assessable=False,
            convergence_time_ns=None,
            message="Too few frames",
        ),
    )
    universe = _make_stretch_universe()
    frame_selection = FrameSelection(start=0, stop=3, step=1, n_frames_total=3, timestep_ps=10.0)
    run = _make_stretch_run_settings()
    settings = RMSDSettings(runs=[run])
    mda_analysis = _build_rmsd_analysis(
        universe=universe,
        run=run,
        frame_selection=frame_selection,
    )
    job = MDAAnalysisJob(
        name=run.label,
        analysis=mda_analysis,
        frame_selection=frame_selection,
        universe_policy=MDAUniversePolicy(
            condition_label=condition.label,
            replicate=1,
            metadata={"rmsd_run": run.model_dump(mode="json")},
        ),
    )

    completed = job.run()

    assert completed.results.rmsd[:, 1].tolist() == pytest.approx([0.0, 5.0, 10.0])
    assert completed.results.rmsd[:, 2].tolist() == pytest.approx(EXPECTED_STRETCH_RMSD)
    expected_rmsd_values = np.asarray(EXPECTED_STRETCH_RMSD, dtype=np.float64)
    expected_mean_rmsd = float(np.mean(expected_rmsd_values))

    collector_ctx, store = _make_rmsd_collector_context(
        condition=condition,
        tmp_path=tmp_path,
        settings=settings,
        frame_selection=frame_selection,
    )
    artifact = RMSDAnalysis().build_mda_collector(collector_ctx)(collector_ctx, [completed])
    run_payload = artifact.payload["runs"][0]
    assert run_payload["mean_rmsd"] == pytest.approx(expected_mean_rmsd)
    assert run_payload["npz_path"] == f"sidecars/{_sidecar_filename(run.label, 0)}"
    with np.load(store.resolve_sidecar(artifact.sidecars[0])) as payload:
        assert payload["rmsd_values"].tolist() == pytest.approx(expected_rmsd_values.tolist())
        assert payload["time_ns"].tolist() == pytest.approx([0.0, 0.005, 0.01])
        assert float(payload["effective_timestep_ps"]) == pytest.approx(5.0)

    analysis = RMSDAnalysis()
    settings_hash = _settings_hash(settings)
    aggregate_ctx = make_aggregate_context(
        condition=condition,
        replicates=(1,),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    aggregated = condition_artifact_to_legacy_result(analysis.aggregate(aggregate_ctx, [artifact]))
    assert aggregated.run_results[0].overall_mean == pytest.approx(expected_mean_rmsd)
    assert aggregated.settings_fingerprint == settings_hash

    comparison_ctx = make_comparison_context(
        name="synthetic_rmsd",
        conditions=[condition],
        analysis_dirs={condition.label: tmp_path / "rmsd"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label=condition.label,
        equilibration="10ns",
        recompute=False,
        aggregated_results={condition.label: aggregated},
    )
    comparison = analysis.compare(comparison_ctx)

    assert comparison is not None
    assert comparison.n_runs == 1
    assert comparison.conditions[0].run_summaries[0].mean_rmsd == pytest.approx(expected_mean_rmsd)


def test_rmsd_full_mda_replicate_lifecycle_persists_artifact(
    condition: Condition, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """RMSD replicate lifecycle should execute MDA jobs and persist artifacts."""

    class FakeLoader:
        """Synthetic trajectory loader for full-lifecycle RMSD tests."""

        def __init__(self, sim_config: object, *args: object, **kwargs: object) -> None:
            self.sim_config = sim_config

        def load_universe(self, replicate: int, cache: bool = True) -> object:
            del replicate, cache
            return _make_stretch_universe()

        def get_timestep(self, replicate: int, unit: str = "ps") -> float:
            del replicate
            assert unit == "ps"
            return 5.0

    class FakeProvider:
        """Universe provider that delegates to the synthetic loader."""

        def __init__(self, loader: FakeLoader) -> None:
            self.loader = loader

        @classmethod
        def from_config(cls, config: object, **kwargs: object) -> "FakeProvider":
            del config
            return cls(kwargs["loader"])

        def load_universe(self, replicate: int) -> object:
            return self.loader.load_universe(replicate)

        def provenance_for(self, replicate: int) -> dict[str, object]:
            return {"replicate": replicate, "warnings": []}

    monkeypatch.setattr("polyzymd.analyses.rmsd._mda.compute_config_hash", lambda _cfg: "hash123")
    monkeypatch.setattr("polyzymd.analyses.shared.loader.TrajectoryLoader", FakeLoader)
    monkeypatch.setattr(
        "polyzymd.analyses.shared.convergence.find_convergence_time",
        lambda *args, **kwargs: SimpleNamespace(
            window_start_times_ns=[],
            window_mean_values=[],
            slope_times_ns=[],
            slopes=[],
            converged=False,
            assessable=False,
            convergence_time_ns=None,
            message="Too few frames",
        ),
    )
    condition = make_condition(label=condition.label, replicates=(1,))
    analysis = RMSDAnalysis()
    monkeypatch.setattr(analysis, "_trajectory_loader_factory", lambda: FakeLoader)
    monkeypatch.setattr(analysis, "_mda_universe_provider_factory", lambda: FakeProvider)
    run = _make_stretch_run_settings(label="lifecycle_stretch")
    settings = RMSDSettings(runs=[run])
    output_dir = tmp_path / "rmsd" / "run_1"

    artifact = AnalysisLifecycle(analysis).run_replicate_once(
        condition=condition,
        settings=settings,
        equilibration="0ns",
        output_dir=output_dir,
        replicate=1,
        recompute=True,
    )

    assert isinstance(artifact, ReplicateArtifact)
    assert analysis.replicate_result_path(output_dir).exists()
    store = ArtifactStore(output_dir)
    persisted = store.read_replicate_result()
    assert persisted.analysis_name == "rmsd"
    assert persisted.condition_label == condition.label
    assert persisted.replicate == 1
    run_payload = persisted.payload["runs"][0]
    assert run_payload["mean_rmsd"] == pytest.approx(np.mean(EXPECTED_STRETCH_RMSD))
    assert run_payload["npz_path"] == f"sidecars/{_sidecar_filename(run.label, 0)}"
    assert persisted.sidecars[0].path == run_payload["npz_path"]
    with np.load(store.resolve_sidecar(persisted.sidecars[0])) as payload:
        assert payload["rmsd_values"].tolist() == pytest.approx(EXPECTED_STRETCH_RMSD)
        assert payload["time_ns"].tolist() == pytest.approx([0.0, 0.005, 0.01])


def test_plotter_resolves_npz_from_replicate_artifact(tmp_path: Path) -> None:
    """Plotter should resolve NPZ paths from canonical MDA artifacts."""
    run_dir = tmp_path / "rmsd" / "run_1"
    run_dir.mkdir(parents=True, exist_ok=True)

    store = ArtifactStore(run_dir)
    sidecar = store.write_npz_sidecar(
        "sidecars/rmsd_protein_backbone_timeseries.npz",
        rmsd_values=np.asarray([1.0, 1.2], dtype=np.float64),
        time_ns=np.asarray([0.0, 0.01], dtype=np.float64),
    )
    artifact = _make_replicate_artifact(
        1,
        [
            _make_run_payload(1, "protein_backbone", 1.1)
            | {"npz_path": sidecar.path, "sidecar": sidecar.model_dump(mode="json")}
        ],
        condition_label="Control",
    )
    store.write_replicate_result(artifact)

    resolved = _resolve_npz_sidecar_path(
        tmp_path / "rmsd",
        "protein_backbone",
        1,
    )

    assert resolved == store.resolve_sidecar(sidecar)


def test_rmsd_external_reference_requires_file() -> None:
    """RMSDRunSettings should require reference_file for external mode."""
    with pytest.raises(ValueError, match="reference_file is required"):
        RMSDRunSettings(label="protein_backbone", reference_mode="external")


def test_rmsd_run_result_creation() -> None:
    """RMSDRunResult should be constructible with all required fields."""
    result = _make_run_result(replicate=1, run_label="protein_backbone", mean_rmsd=1.2)
    assert result.run_label == "protein_backbone"
    assert result.mean_rmsd == pytest.approx(1.2)
    assert result.convergence_assessable is True


def test_rmsd_run_result_summary() -> None:
    """RMSDRunResult.summary should produce readable output."""
    result = _make_run_result(replicate=1, run_label="protein_backbone", mean_rmsd=1.2)
    summary = result.summary()
    assert "RMSD Run: protein_backbone" in summary
    assert "Mean:" in summary


def test_rmsd_result_n_runs() -> None:
    """RMSDResult.n_runs should match number of run entries."""
    result = RMSDResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=1,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        run_results=[
            _make_run_result(replicate=1, run_label="protein_backbone", mean_rmsd=1.2),
            _make_run_result(replicate=1, run_label="polymer_core", mean_rmsd=2.0),
        ],
        n_frames_total=100,
        n_frames_used=90,
        trajectory_files=["/fake/traj.dcd"],
    )
    assert result.n_runs == 2


def test_rmsd_aggregated_result_creation() -> None:
    """RMSDAggregatedResult should be constructible with run aggregates."""
    aggregated = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [1.1, 1.2, 1.3]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [1.9, 2.0, 2.1]),
        ],
        source_result_files=[],
    )
    assert aggregated.n_runs == 2
    assert aggregated.n_replicates == 3


def test_rmsd_aggregated_result_summary() -> None:
    """RMSDAggregatedResult.summary should include key metadata."""
    aggregated = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA",
        replicates=[1, 2],
        n_replicates=2,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [1.1, 1.3]),
        ],
        source_result_files=[],
    )
    summary = aggregated.summary()
    assert "RMSD Aggregated Analysis" in summary
    assert "Runs analyzed: 1" in summary


def test_rmsd_comparison_result_creation() -> None:
    """RMSDComparisonResult should be constructible with full fields."""
    result = _make_comparison_result()
    assert result.n_runs == 2
    assert len(result.conditions) == 3


def test_rmsd_comparison_get_run() -> None:
    """RMSDComparisonResult.get_comparisons_for_run should filter by run label."""
    result = _make_comparison_result()
    comps = result.get_comparisons_for_run("protein_backbone")
    assert len(comps) == 2
    assert all(comp.run_label == "protein_backbone" for comp in comps)


def test_rmsd_comparison_ranking() -> None:
    """RMSDComparisonResult.get_ranking should return run-specific ranking."""
    result = _make_comparison_result()
    ranking = result.get_ranking("protein_backbone")
    assert ranking[0] == "Treatment_A"


def test_rmsd_condition_summary_get_run() -> None:
    """RMSDConditionSummary.get_run should return matching run summary."""
    condition = _make_condition_summary("Control", 1.2, 2.0)
    run = condition.get_run("protein_backbone")
    assert run.label == "protein_backbone"
    assert run.mean_rmsd == pytest.approx(1.2)


def test_rmsd_plugin_discovered() -> None:
    """RMSD plugin should be auto-discovered by analyses discovery."""
    from polyzymd.analyses.discovery import clear_cache, list_analyses

    clear_cache()
    analyses = list_analyses()
    assert "rmsd" in analyses
    assert analyses["rmsd"] is RMSDAnalysis


def test_rmsd_plugin_name() -> None:
    """RMSD plugin should expose correct analysis name."""
    assert RMSDAnalysis.name == "rmsd"


def test_rmsd_min_replicates() -> None:
    """RMSD plugin should allow one-replicate smoke-test runs."""
    assert RMSDAnalysis.min_replicates == 1


def test_format_table() -> None:
    """Table formatting should produce human-readable RMSD output."""
    text = format_rmsd_comparison(_make_comparison_result(), "table")
    assert "RMSD Comparison: protein_backbone" in text
    assert "Condition" in text
    assert "Convergence:" in text


def test_format_markdown() -> None:
    """Markdown formatting should include markdown sections and table headers."""
    text = format_rmsd_comparison(_make_comparison_result(), "markdown")
    assert "## RMSD Comparison: protein_backbone" in text
    assert "| Condition | Mean RMSD (Å) | SEM | Rank |" in text


def test_format_markdown_regression_includes_non_significant_and_non_testable() -> None:
    """Markdown formatter should keep stable output for edge comparison states."""
    result = RMSDComparisonResult(
        metric="mean_rmsd",
        name="regression_case",
        n_runs=1,
        run_labels=["run_1"],
        control_label="Control",
        conditions=[
            RMSDConditionSummary(
                label="Control",
                config_path="/fake/control.yaml",
                n_replicates=2,
                run_summaries=[
                    RMSDRunSummary(
                        label="run_1",
                        selection="protein and name CA",
                        mean_rmsd=1.00,
                        sem_rmsd=0.10,
                        per_replicate_means=[0.95, 1.05],
                        n_converged_replicates=2,
                        n_assessable_replicates=2,
                        convergence_fraction=1.0,
                        all_converged=True,
                        mean_convergence_time_ns=10.0,
                        median_convergence_time_ns=10.0,
                    )
                ],
            ),
            RMSDConditionSummary(
                label="Treatment_A",
                config_path="/fake/treatment_a.yaml",
                n_replicates=2,
                run_summaries=[
                    RMSDRunSummary(
                        label="run_1",
                        selection="protein and name CA",
                        mean_rmsd=1.05,
                        sem_rmsd=0.09,
                        per_replicate_means=[1.00, 1.10],
                        n_converged_replicates=1,
                        n_assessable_replicates=2,
                        convergence_fraction=0.5,
                        all_converged=False,
                        mean_convergence_time_ns=12.5,
                        median_convergence_time_ns=12.5,
                    )
                ],
            ),
            RMSDConditionSummary(
                label="Treatment_B",
                config_path="/fake/treatment_b.yaml",
                n_replicates=1,
                run_summaries=[
                    RMSDRunSummary(
                        label="run_1",
                        selection="protein and name CA",
                        mean_rmsd=1.08,
                        sem_rmsd=0.00,
                        per_replicate_means=[1.08],
                        n_converged_replicates=0,
                        n_assessable_replicates=1,
                        convergence_fraction=0.0,
                        all_converged=False,
                        mean_convergence_time_ns=None,
                        median_convergence_time_ns=None,
                    )
                ],
            ),
        ],
        pairwise_comparisons=[
            RMSDRunPairwiseComparison(
                run_label="run_1",
                condition_a="Control",
                condition_b="Treatment_A",
                t_statistic=1.5,
                p_value=0.120,
                cohens_d=0.6,
                effect_interpretation="medium",
                direction="destabilizing",
                significant=False,
                percent_change=5.0,
                testable=True,
            ),
            RMSDRunPairwiseComparison(
                run_label="run_1",
                condition_a="Control",
                condition_b="Treatment_B",
                effect_interpretation="negligible",
                direction="destabilizing",
                significant=False,
                percent_change=8.0,
                testable=False,
                note="Insufficient replicate count",
            ),
        ],
        anova_by_run=[
            RMSDRunANOVA(
                run_label="run_1",
                significant=False,
                testable=False,
                note="Insufficient data for ANOVA",
            )
        ],
        ranking_by_run={"run_1": ["Control", "Treatment_A", "Treatment_B"]},
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.2.1",
    )

    text = format_rmsd_comparison(result, "markdown")

    assert text == (
        "## RMSD Comparison: run_1\n"
        "\n"
        "| Condition | Mean RMSD (Å) | SEM | Rank |\n"
        "|-----------|---------------|-----|------|\n"
        "| Control | 1.00 | 0.10 | 1 |\n"
        "  - Convergence: 2/2 replicates converged (2 assessable), median t_conv = 10.0 ns\n"
        "| Treatment_A | 1.05 | 0.09 | 2 |\n"
        "  - Convergence: 1/2 replicates converged (2 assessable), median t_conv = 12.5 ns\n"
        "| Treatment_B | 1.08 | n/a | 3 |\n"
        "  - Convergence: 0/1 replicates converged (1 assessable), median t_conv = n/a\n"
        "\n"
        "*SEM: n/a (single replicate; not estimable).*\n"
        "\n"
        "- Pairwise: Treatment_A vs Control — Δ=+5.0%, p=0.120 , d=0.60 (medium), destabilizing\n"
        "- Pairwise: Treatment_B vs Control — Δ=+8.0%, destabilizing; "
        "Insufficient replicate count\n"
        "\n"
        "- ANOVA: Insufficient data for ANOVA\n"
    )


def test_format_json() -> None:
    """JSON formatting should return valid JSON with expected keys."""
    text = format_rmsd_comparison(_make_comparison_result(), "json")
    parsed = json.loads(text)
    assert parsed["metric"] == "mean_rmsd"
    assert "conditions" in parsed


def test_rmsd_plot_settings_model_attribute() -> None:
    """RMSD analysis should expose its plot settings model attribute."""
    cls = get_analysis("rmsd")
    assert cls.PlotSettingsModel is RMSDPlotSettings


def test_rmsd_plot_settings_defaults() -> None:
    """RMSD plot settings defaults should match plugin defaults."""
    settings = RMSDPlotSettings()
    assert settings.show_per_replicate is False
    assert settings.figsize == (10, 6)
    assert settings.timeseries_figsize == (12, 5)
    assert settings.show_convergence_plots is False
    assert settings.convergence_figsize == (12.0, 5.0)


def test_aggregate_single_replicate(
    condition: Condition, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """aggregate should handle a single artifact and log SEM warning."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    settings_hash = _settings_hash(settings)
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1,),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    artifact = _make_replicate_artifact(
        1,
        [_make_run_payload(1, "protein_backbone", 1.2), _make_run_payload(1, "polymer_core", 2.0)],
        condition_label=condition.label,
        settings_hash=settings_hash,
    )

    with caplog.at_level("WARNING"):
        aggregated = analysis.aggregate(ctx, [artifact])

    legacy = condition_artifact_to_legacy_result(aggregated)
    assert aggregated.replicates == [1]
    assert legacy.n_replicates == 1
    assert len(legacy.run_results) == 2
    assert "Only one replicate available for RMSD aggregation" in caplog.text


def test_aggregate_multiple_replicates(condition: Condition, tmp_path: Path) -> None:
    """aggregate should compute expected means and preserve run structure."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    settings_hash = _settings_hash(settings)
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2, 3),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    results = [
        _make_replicate_artifact(
            rep,
            [
                _make_run_payload(
                    rep,
                    "protein_backbone",
                    1.0 + 0.1 * rep,
                    convergence_time_ns=(20.0 + 10.0 * rep) if rep < 3 else None,
                    convergence_assessable=True,
                ),
                _make_run_payload(
                    rep,
                    "polymer_core",
                    2.0 + 0.05 * rep,
                    convergence_time_ns=15.0 + 5.0 * rep,
                    convergence_assessable=True,
                ),
            ],
            condition_label=condition.label,
            settings_hash=settings_hash,
        )
        for rep in (1, 2, 3)
    ]

    aggregated = condition_artifact_to_legacy_result(analysis.aggregate(ctx, results))

    assert aggregated.n_replicates == 3
    backbone = next(run for run in aggregated.run_results if run.run_label == "protein_backbone")
    assert backbone.overall_mean == pytest.approx(1.2)
    assert backbone.per_replicate_means == pytest.approx([1.1, 1.2, 1.3])
    assert backbone.n_assessable_replicates == 3
    assert backbone.n_converged_replicates == 2
    assert backbone.convergence_fraction == pytest.approx(2.0 / 3.0)
    assert backbone.all_converged is False
    assert backbone.median_convergence_time_ns == pytest.approx(35.0)


def test_aggregate_orders_complete_out_of_order_inputs(
    condition: Condition, tmp_path: Path
) -> None:
    """aggregate should align per-replicate arrays to declared replicate order."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=[RMSDRunSettings(label="protein_backbone")])
    settings_hash = _settings_hash(settings)
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2, 3),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    results = [
        _make_replicate_artifact(
            rep,
            [_make_run_payload(rep, "protein_backbone", mean_rmsd=1.0 + 0.1 * rep)],
            condition_label=condition.label,
            settings_hash=settings_hash,
        )
        for rep in (3, 1, 2)
    ]

    aggregated = condition_artifact_to_legacy_result(analysis.aggregate(ctx, results))

    backbone = aggregated.run_results[0]
    assert backbone.replicates == [1, 2, 3]
    assert backbone.per_replicate_means == pytest.approx([1.1, 1.2, 1.3])
    assert backbone.per_replicate_medians == pytest.approx([1.1, 1.2, 1.3])


def test_aggregate_empty_results_raises(condition: Condition, tmp_path: Path) -> None:
    """aggregate should fail clearly when replicate inputs are missing."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2, 3),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )

    with pytest.raises(
        ValueError,
        match=r"RMSD aggregation for condition '.+' requires at least one replicate artifact",
    ):
        analysis.aggregate(ctx, [])


def test_aggregate_missing_configured_run_raises(condition: Condition, tmp_path: Path) -> None:
    """aggregate should fail when a configured run is missing for any replicate."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    settings_hash = _settings_hash(settings)
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2, 3),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    results = [
        _make_replicate_artifact(
            1,
            [
                _make_run_payload(1, "protein_backbone", 1.1),
                _make_run_payload(1, "polymer_core", 2.1),
            ],
            condition_label=condition.label,
            settings_hash=settings_hash,
        ),
        _make_replicate_artifact(
            2,
            [_make_run_payload(2, "protein_backbone", 1.2)],
            condition_label=condition.label,
            settings_hash=settings_hash,
        ),
        _make_replicate_artifact(
            3,
            [
                _make_run_payload(3, "protein_backbone", 1.3),
                _make_run_payload(3, "polymer_core", 2.3),
            ],
            condition_label=condition.label,
            settings_hash=settings_hash,
        ),
    ]

    with pytest.raises(ValueError, match="missing configured run"):
        analysis.aggregate(ctx, results)


def test_aggregate_overall_median_uses_median(condition: Condition, tmp_path: Path) -> None:
    """aggregate should compute overall_median using np.median."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=[RMSDRunSettings(label="protein_backbone")])
    settings_hash = _settings_hash(settings)
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2, 3),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )

    results = [
        _make_replicate_artifact(
            rep,
            [
                _make_run_result(
                    replicate=rep,
                    run_label="protein_backbone",
                    mean_rmsd=(2.0 if rep < 3 else 101.0),
                ).model_dump()
                | {"median_rmsd": {1: 1.0, 2: 2.0, 3: 100.0}[rep]}
            ],
            condition_label=condition.label,
            settings_hash=settings_hash,
        )
        for rep in (1, 2, 3)
    ]

    aggregated = condition_artifact_to_legacy_result(analysis.aggregate(ctx, results))
    backbone = next(run for run in aggregated.run_results if run.run_label == "protein_backbone")
    assert backbone.overall_median == pytest.approx(2.0)


def test_compare_two_conditions(tmp_path: Path) -> None:
    """compare with two conditions should produce pairwise results without ANOVA."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    settings_hash = _settings_hash(settings)
    control = make_condition(label="Control")
    treated = make_condition(label="Treated")

    control_agg = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [1.20, 1.25, 1.15]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [2.00, 2.05, 1.95]),
        ],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )
    treated_agg = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [1.00, 1.05, 0.95]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [1.80, 1.85, 1.75]),
        ],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )

    ctx = make_comparison_context(
        name="rmsd_compare",
        conditions=[control, treated],
        analysis_dirs={"Control": tmp_path / "control", "Treated": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="10ns",
        recompute=False,
        aggregated_results={"Control": control_agg, "Treated": treated_agg},
    )

    comparison = analysis.compare(ctx)

    assert comparison is not None
    assert comparison.anova_by_run is None
    assert len(comparison.pairwise_comparisons) == 2


def test_compare_three_conditions(tmp_path: Path) -> None:
    """compare with three conditions should include ANOVA per run."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    settings_hash = _settings_hash(settings)
    conditions = [
        make_condition(label="Control"),
        make_condition(label="A"),
        make_condition(label="B"),
    ]

    aggregated_results = {
        "Control": RMSDAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA; segid C and backbone",
            replicates=[1, 2, 3],
            n_replicates=3,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [1.20, 1.25, 1.15]),
                _make_aggregated_run("polymer_core", "segid C and backbone", [2.00, 2.05, 1.95]),
            ],
            settings_fingerprint=settings_hash,
            source_result_files=[],
        ),
        "A": RMSDAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA; segid C and backbone",
            replicates=[1, 2, 3],
            n_replicates=3,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [1.00, 1.05, 0.95]),
                _make_aggregated_run("polymer_core", "segid C and backbone", [1.80, 1.85, 1.75]),
            ],
            settings_fingerprint=settings_hash,
            source_result_files=[],
        ),
        "B": RMSDAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA; segid C and backbone",
            replicates=[1, 2, 3],
            n_replicates=3,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [1.35, 1.40, 1.30]),
                _make_aggregated_run("polymer_core", "segid C and backbone", [2.15, 2.20, 2.10]),
            ],
            settings_fingerprint=settings_hash,
            source_result_files=[],
        ),
    }

    ctx = make_comparison_context(
        name="rmsd_compare",
        conditions=conditions,
        analysis_dirs={c.label: tmp_path / c.label for c in conditions},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="10ns",
        recompute=False,
        aggregated_results=aggregated_results,
    )

    comparison = analysis.compare(ctx)

    assert comparison is not None
    assert comparison.anova_by_run is not None
    assert len(comparison.anova_by_run) == 2


def test_compare_single_replicate_not_testable(tmp_path: Path) -> None:
    """Pairwise and ANOVA results should be marked not testable for n < 2."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=[RMSDRunSettings(label="protein_backbone")])
    settings_hash = _settings_hash(settings)
    conditions = [
        make_condition(label="Control", replicates=(1,)),
        make_condition(label="A", replicates=(1,)),
        make_condition(label="B", replicates=(1,)),
    ]

    aggregated_results = {
        "Control": RMSDAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA",
            replicates=[1],
            n_replicates=1,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [1.20]),
            ],
            settings_fingerprint=settings_hash,
            source_result_files=[],
        ),
        "A": RMSDAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA",
            replicates=[1],
            n_replicates=1,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [1.00]),
            ],
            settings_fingerprint=settings_hash,
            source_result_files=[],
        ),
        "B": RMSDAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA",
            replicates=[1],
            n_replicates=1,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [1.40]),
            ],
            settings_fingerprint=settings_hash,
            source_result_files=[],
        ),
    }

    ctx = make_comparison_context(
        name="rmsd_compare",
        conditions=conditions,
        analysis_dirs={c.label: tmp_path / c.label for c in conditions},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="10ns",
        recompute=False,
        aggregated_results=aggregated_results,
    )

    comparison = analysis.compare(ctx)

    assert comparison is not None
    assert comparison.pairwise_comparisons
    assert all(not comp.testable for comp in comparison.pairwise_comparisons)
    assert comparison.anova_by_run is not None
    assert comparison.anova_by_run[0].testable is False


def test_compare_missing_configured_run_raises(tmp_path: Path) -> None:
    """compare should fail when an aggregated condition omits a configured run."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    settings_hash = _settings_hash(settings)
    control = make_condition(label="Control")
    treated = make_condition(label="Treated")

    control_agg = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [1.2, 1.25, 1.15]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [2.0, 2.05, 1.95]),
        ],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )
    treated_agg = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [1.0, 1.05, 0.95]),
        ],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )

    ctx = make_comparison_context(
        name="rmsd_compare",
        conditions=[control, treated],
        analysis_dirs={"Control": tmp_path / "control", "Treated": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="10ns",
        recompute=False,
        aggregated_results={"Control": control_agg, "Treated": treated_agg},
    )

    with pytest.raises(
        ValueError, match="Aggregated RMSD result for condition 'Treated' is incomplete"
    ):
        analysis.compare(ctx)


def test_compare_missing_condition_aggregated_result_raises(tmp_path: Path) -> None:
    """compare should fail when any configured condition lacks an aggregated result."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    settings_hash = _settings_hash(settings)
    control = make_condition(label="Control")
    treated = make_condition(label="Treated")

    control_agg = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [1.2, 1.25, 1.15]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [2.0, 2.05, 1.95]),
        ],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )

    ctx = make_comparison_context(
        name="rmsd_compare",
        conditions=[control, treated],
        analysis_dirs={"Control": tmp_path / "control", "Treated": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="10ns",
        recompute=False,
        aggregated_results={"Control": control_agg},
    )

    with pytest.raises(
        ValueError,
        match="Missing aggregated RMSD result for condition 'Treated'",
    ):
        analysis.compare(ctx)


def test_compare_incomplete_run_replicate_values_raise(tmp_path: Path) -> None:
    """compare should fail when a configured run has partial replicate values."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=[RMSDRunSettings(label="protein_backbone")])
    settings_hash = _settings_hash(settings)
    control = make_condition(label="Control")
    treated = make_condition(label="Treated")

    control_agg = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [1.2, 1.25, 1.15]),
        ],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )
    treated_run = _make_aggregated_run("protein_backbone", "protein and name CA", [1.0, 1.05])
    treated_agg = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[treated_run],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )

    ctx = make_comparison_context(
        name="rmsd_compare",
        conditions=[control, treated],
        analysis_dirs={"Control": tmp_path / "control", "Treated": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="10ns",
        recompute=False,
        aggregated_results={"Control": control_agg, "Treated": treated_agg},
    )

    with pytest.raises(
        ValueError,
        match="Aggregated RMSD run 'protein_backbone' for condition 'Treated' has incomplete replicate metadata",
    ):
        analysis.compare(ctx)


@pytest.mark.parametrize(
    ("result_kind", "error_match"),
    [
        pytest.param("stale", "current settings require", id="stale-fingerprint"),
        pytest.param(
            "legacy",
            "missing a settings fingerprint",
            id="missing-fingerprint",
        ),
    ],
)
def test_compare_rejects_preloaded_invalid_aggregated_results(
    result_kind: str,
    error_match: str,
    tmp_path: Path,
) -> None:
    """compare should validate preloaded aggregated RMSD results."""
    analysis = RMSDAnalysis()
    current_settings = RMSDSettings(runs=[RMSDRunSettings(label="protein_backbone")])
    stale_settings = RMSDSettings(
        runs=[RMSDRunSettings(label="protein_backbone", selection="protein")]
    )
    condition = make_condition(label="CondA", replicates=(1, 2))

    selection = "protein" if result_kind == "stale" else "protein and name CA"
    settings_fingerprint = _settings_hash(stale_settings) if result_kind == "stale" else None
    preloaded_result = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string=selection,
        replicates=[1, 2],
        n_replicates=2,
        run_results=[_make_aggregated_run("protein_backbone", selection, [1.1, 1.2])],
        settings_fingerprint=settings_fingerprint,
        source_result_files=[],
    ).model_dump()

    ctx = make_comparison_context(
        name="rmsd_compare",
        conditions=[condition],
        analysis_dirs={"CondA": tmp_path / "analysis" / "conda" / "rmsd"},
        results_dir=tmp_path / "comparison",
        settings=current_settings,
        control_label="CondA",
        equilibration="10ns",
        recompute=False,
        aggregated_results={"CondA": preloaded_result},
    )

    with pytest.raises(ValueError, match=error_match):
        analysis.compare(ctx)


def test_compare_rejects_stale_aggregated_result_from_disk(tmp_path: Path) -> None:
    """compare should fail loudly when aggregated RMSD cache settings are stale."""
    analysis = RMSDAnalysis()
    current_settings = RMSDSettings(runs=[RMSDRunSettings(label="protein_backbone")])
    stale_settings = RMSDSettings(
        runs=[RMSDRunSettings(label="protein_backbone", selection="protein")]
    )
    condition = make_condition(label="CondA", replicates=(1, 2))
    analysis_dir = tmp_path / "analysis" / "conda" / "rmsd"
    aggregated_dir = analysis_dir / "aggregated"
    aggregated_dir.mkdir(parents=True)

    RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein",
        replicates=[1, 2],
        n_replicates=2,
        run_results=[_make_aggregated_run("protein_backbone", "protein", [1.1, 1.2])],
        settings_fingerprint=_settings_hash(stale_settings),
        source_result_files=[],
    ).save(aggregated_dir / "result.json")

    ctx = make_comparison_context(
        name="rmsd_compare",
        conditions=[condition],
        analysis_dirs={"CondA": analysis_dir},
        results_dir=tmp_path / "comparison",
        settings=current_settings,
        control_label="CondA",
        equilibration="10ns",
        recompute=False,
    )

    with pytest.raises(ValueError, match="current settings require"):
        analysis.compare(ctx)


def test_plot_rejects_legacy_aggregated_result_from_disk(tmp_path: Path) -> None:
    """plot should fail loudly when aggregated RMSD cache lacks settings identity."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=[RMSDRunSettings(label="protein_backbone")])
    condition = make_condition(label="CondA", replicates=(1, 2))
    analysis_dir = tmp_path / "analysis" / "conda" / "rmsd"
    aggregated_dir = analysis_dir / "aggregated"
    aggregated_dir.mkdir(parents=True)

    RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA",
        replicates=[1, 2],
        n_replicates=2,
        run_results=[_make_aggregated_run("protein_backbone", "protein and name CA", [1.1, 1.2])],
        source_result_files=[],
    ).save(aggregated_dir / "result.json")

    results_dir = tmp_path / "comparison"
    results_dir.mkdir(parents=True)
    _make_comparison_result().save(results_dir / "result.json")

    ctx = PlotContext(
        conditions=[condition],
        analysis_dirs={"CondA": analysis_dir},
        results_dir=results_dir,
        output_dir=tmp_path / "figures",
        settings=settings,
        plot_settings=PlotSettings(),
        equilibration="10ns",
    )

    with pytest.raises(ValueError, match="missing a settings fingerprint"):
        analysis.plot(ctx)
