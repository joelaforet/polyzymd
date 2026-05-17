"""Tests for the SASA analysis plugin."""

from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest

from polyzymd.analyses.base import Condition, PlotContext
from polyzymd.analyses.mda import (
    ArtifactStore,
    ArtifactStoreError,
    ConditionArtifact,
    FrameSelection,
    MDABackendPolicy,
    MDACollectorContext,
    MDAJobResult,
    MDAUniversePolicy,
    ReplicateArtifact,
)
from polyzymd.analyses.sasa import SASAAnalysis, SASARunSettings, SASASettings
from polyzymd.analyses.sasa._comparison_results import (
    SASAComparisonResult,
    SASAConditionSummary,
    SASARunPairwiseComparison,
    SASARunSummary,
)
from polyzymd.analyses.sasa._formatters import format_sasa_comparison
from polyzymd.analyses.sasa._mda import (
    SASAArtifactCollector,
    build_sasa_jobs,
    condition_artifact_to_legacy_result,
    load_replicate_sasa_sidecar,
)
from polyzymd.analyses.sasa._plot_settings import SASAPlotSettings
from polyzymd.analyses.sasa._plotters import (
    _build_sasa_normalized_control_rows,
    _build_sasa_normalized_replicate_deltas,
    _load_condition_result_payloads,
    _load_replicate_timeseries_from_results,
    _propagate_sasa_normalized_sem,
    _sanitize_run_label,
    plot_sasa_normalized_control_bars,
)
from polyzymd.analyses.shared.sasa import SASAComputationResult
from tests._support.analysis_testkit import (
    make_aggregate_context,
    make_comparison_context,
    make_condition,
    make_replicate_context,
)


def _make_condition(label: str, replicates: tuple[int, ...] = (1, 2)) -> Condition:
    """Create a lightweight comparison condition."""

    return make_condition(
        label=label,
        config_path=Path(f"/fake/{label}.yaml"),
        replicates=replicates,
        sim_config=MagicMock(),
    )


def _raw_sasa(
    *,
    total: list[float],
    residue_offset: float = 0.0,
    residue_key: str = "A:1:ALA",
    zero: bool = False,
) -> SASAComputationResult:
    """Create a synthetic raw SASA result."""

    frames = np.arange(len(total), dtype=np.int64)
    if zero:
        return SASAComputationResult(
            atom_sasa_a2=np.empty((len(total), 0), dtype=np.float64),
            residue_sasa_a2=np.empty((len(total), 0), dtype=np.float64),
            total_sasa_a2=np.zeros(len(total), dtype=np.float64),
            frames=frames,
            time_ns=frames.astype(np.float64) * 0.01,
            target_atom_indices=np.asarray([], dtype=np.int64),
            context_atom_indices=np.asarray([], dtype=np.int64),
            residue_keys=[],
            residue_chainids=[],
            residue_resids=[],
            residue_resnames=[],
        )

    residue_resid = int(residue_key.split(":")[1])
    residue_resname = residue_key.split(":")[2]
    residue_values = np.asarray(total, dtype=np.float64).reshape(-1, 1) + residue_offset
    return SASAComputationResult(
        atom_sasa_a2=residue_values.copy(),
        residue_sasa_a2=residue_values,
        total_sasa_a2=np.asarray(total, dtype=np.float64),
        frames=frames,
        time_ns=frames.astype(np.float64) * 0.01,
        target_atom_indices=np.asarray([0], dtype=np.int64),
        context_atom_indices=np.asarray([0], dtype=np.int64),
        residue_keys=[residue_key],
        residue_chainids=[residue_key.split(":")[0]],
        residue_resids=[residue_resid],
        residue_resnames=[residue_resname],
    )


def _job_result(
    raw: SASAComputationResult,
    *,
    run_label: str = "protein",
    frame_selection: FrameSelection | None = None,
) -> MDAJobResult:
    """Create a completed MDA job result around raw SASA arrays."""

    results = SimpleNamespace(
        run_label=run_label,
        target_selection="chainid A",
        context_selection="chainid A",
        probe_radius_nm=0.14,
        n_sphere_points=960,
        chunk_size=100,
        atom_sasa_a2=raw.atom_sasa_a2,
        residue_sasa_a2=raw.residue_sasa_a2,
        total_sasa_a2=raw.total_sasa_a2,
        frames=raw.frames,
        time_ns=raw.time_ns,
        target_atom_indices=raw.target_atom_indices,
        context_atom_indices=raw.context_atom_indices,
        residue_keys=raw.residue_keys,
        residue_chainids=raw.residue_chainids,
        residue_resids=raw.residue_resids,
        residue_resnames=raw.residue_resnames,
    )
    selection = frame_selection or FrameSelection(start=0, stop=len(raw.total_sasa_a2), step=1)
    return MDAJobResult(
        name=f"sasa:{run_label}",
        analysis=object(),
        results=results,
        run_kwargs=selection.run_kwargs(),
        frame_selection=selection,
        backend_policy=MDABackendPolicy(),
        universe_policy=MDAUniversePolicy(),
    )


def _collect_artifact(
    tmp_path: Path,
    *,
    condition_label: str,
    replicate: int,
    raw: SASAComputationResult,
    settings: SASASettings | None = None,
    run_label: str = "protein",
) -> ReplicateArtifact:
    """Collect and persist a canonical SASA replicate artifact."""

    settings = settings or SASASettings(
        runs=[SASARunSettings(label=run_label, target_selection="chainid A")]
    )
    condition = _make_condition(condition_label)
    ctx = make_replicate_context(
        condition=condition,
        replicate=replicate,
        output_dir=tmp_path / f"run_{replicate}",
        settings=settings,
        equilibration="10ns",
        recompute=True,
    )
    collector_ctx = MDACollectorContext(
        analysis_name="sasa",
        replicate_context=ctx,
        frame_selection=FrameSelection(
            start=0,
            stop=len(raw.total_sasa_a2),
            step=1,
            timestep_ps=10.0,
            n_frames_total=len(raw.total_sasa_a2),
            equilibration="10ns",
        ),
        universe_policy=MDAUniversePolicy(condition_label=condition_label, replicate=replicate),
        artifact_store=ArtifactStore(ctx.output_dir),
        settings_fingerprint=SASAAnalysis._make_settings_cache_tag(settings),
    )
    artifact = SASAArtifactCollector()(collector_ctx, [_job_result(raw, run_label=run_label)])
    ArtifactStore(ctx.output_dir).write_replicate_result(artifact, "result.json")
    return artifact


def _collect_multi_run_artifact(
    tmp_path: Path,
    *,
    condition_label: str,
    replicate: int,
    raw_by_run: dict[str, SASAComputationResult],
    settings: SASASettings,
) -> ReplicateArtifact:
    """Collect and persist a canonical multi-run SASA replicate artifact."""

    condition = _make_condition(condition_label)
    ctx = make_replicate_context(
        condition=condition,
        replicate=replicate,
        output_dir=tmp_path / f"run_{replicate}",
        settings=settings,
        equilibration="10ns",
        recompute=True,
    )
    frame_count = max(len(raw.total_sasa_a2) for raw in raw_by_run.values())
    collector_ctx = MDACollectorContext(
        analysis_name="sasa",
        replicate_context=ctx,
        frame_selection=FrameSelection(
            start=0,
            stop=frame_count,
            step=1,
            timestep_ps=10.0,
            n_frames_total=frame_count,
            equilibration="10ns",
        ),
        universe_policy=MDAUniversePolicy(condition_label=condition_label, replicate=replicate),
        artifact_store=ArtifactStore(ctx.output_dir),
        settings_fingerprint=SASAAnalysis._make_settings_cache_tag(settings),
    )
    completed_jobs = [
        _job_result(raw, run_label=run_label) for run_label, raw in raw_by_run.items()
    ]
    artifact = SASAArtifactCollector()(collector_ctx, completed_jobs)
    ArtifactStore(ctx.output_dir).write_replicate_result(artifact, "result.json")
    return artifact


def _aggregate_artifacts(
    tmp_path: Path,
    artifacts: list[ReplicateArtifact],
    *,
    condition_label: str = "control",
    settings: SASASettings | None = None,
) -> ConditionArtifact:
    """Aggregate canonical replicate artifacts."""

    settings = settings or SASASettings(
        runs=[SASARunSettings(label="protein", target_selection="chainid A")]
    )
    analysis = SASAAnalysis()
    condition = _make_condition(
        condition_label, replicates=tuple(artifact.replicate for artifact in artifacts)
    )
    ctx = make_aggregate_context(
        condition=condition,
        replicates=tuple(artifact.replicate for artifact in artifacts),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    return analysis.aggregate(ctx, artifacts)


def _make_sasa_comparison_result(
    condition_runs: dict[str, dict[str, tuple[float, float]]],
    *,
    control_label: str | None = "control",
) -> SASAComparisonResult:
    """Create an SASA comparison result with per-run means and SEMs."""

    run_labels = list(
        dict.fromkeys(run_label for runs in condition_runs.values() for run_label in runs)
    )
    conditions = []
    for condition_label, runs in condition_runs.items():
        conditions.append(
            SASAConditionSummary(
                label=condition_label,
                config_path=f"/fake/{condition_label}.yaml",
                n_replicates=3,
                run_summaries=[
                    SASARunSummary(
                        label=run_label,
                        target_selection="chainid A",
                        context_selection="all",
                        mean_sasa=mean_sasa,
                        sem_sasa=sem_sasa,
                        per_replicate_means=[mean_sasa],
                    )
                    for run_label, (mean_sasa, sem_sasa) in runs.items()
                ],
            )
        )

    return SASAComparisonResult(
        metric="mean_sasa",
        name="sasa_compare",
        n_runs=len(run_labels),
        run_labels=run_labels,
        control_label=control_label,
        conditions=conditions,
        pairwise_comparisons=[],
        anova_by_run=None,
        ranking_by_run={run_label: list(condition_runs) for run_label in run_labels},
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.0.0",
    )


def _make_sasa_plot_context(
    tmp_path: Path, *, control_label: str | None = "control"
) -> PlotContext:
    """Create a SASA plot context for normalized-control plot tests."""

    from polyzymd.config.comparison import PlotSettings

    return PlotContext(
        conditions=[_make_condition("control"), _make_condition("treated")],
        analysis_dirs={},
        results_dir=tmp_path / "results",
        output_dir=tmp_path / "figures",
        settings=SASASettings(
            runs=[SASARunSettings(label="protein", target_selection="chainid A")]
        ),
        plot_settings=PlotSettings(),
        control_label=control_label,
    )


def test_sasa_plugin_discovered() -> None:
    """SASA plugin should be auto-discovered."""

    from polyzymd.analyses.discovery import clear_cache, list_analyses

    clear_cache()
    analyses = list_analyses()
    assert "sasa" in analyses
    assert analyses["sasa"] is SASAAnalysis


def test_sasa_settings_validation() -> None:
    """SASA settings should enforce run and scalar constraints."""

    settings = SASASettings(runs=[SASARunSettings(label="protein", target_selection="chainid A")])
    assert settings.runs[0].context_selection == "chainid A"
    assert settings.runs[0].stride == 1
    assert settings.chunk_size == 100

    with pytest.raises(ValueError, match="At least one SASA run"):
        SASASettings(runs=[])
    with pytest.raises(ValueError, match="must be unique"):
        SASASettings(
            runs=[
                SASARunSettings(label="dup", target_selection="chainid A"),
                SASARunSettings(label="dup", target_selection="chainid A"),
            ]
        )
    with pytest.raises(ValueError, match="must not contain"):
        SASARunSettings(label="bad/name", target_selection="chainid A")
    with pytest.raises(ValueError, match="probe_radius_nm must be > 0"):
        SASASettings(
            runs=[SASARunSettings(label="protein", target_selection="chainid A")],
            probe_radius_nm=0.0,
        )
    with pytest.raises(ValueError, match="n_sphere_points must be >= 100"):
        SASASettings(
            runs=[SASARunSettings(label="protein", target_selection="chainid A")],
            n_sphere_points=10,
        )
    with pytest.raises(ValueError, match="chunk_size must be >= 1"):
        SASASettings(
            runs=[SASARunSettings(label="protein", target_selection="chainid A")],
            chunk_size=0,
        )
    with pytest.raises(ValueError, match="stride must be >= 1"):
        SASARunSettings(label="protein", target_selection="chainid A", stride=0)


def test_build_mda_jobs_applies_run_stride(monkeypatch: pytest.MonkeyPatch) -> None:
    """MDA job construction should fold run stride into FrameSelection."""

    import polyzymd.analyses.sasa._mda as sasa_mda

    monkeypatch.setattr(sasa_mda, "build_sasa_surface_area_analysis", lambda **_kwargs: MagicMock())
    ctx = SimpleNamespace(
        universe=object(),
        frame_selection=FrameSelection(
            start=5, stop=25, step=2, timestep_ps=10.0, n_frames_total=30
        ),
        replicate_context=SimpleNamespace(condition=_make_condition("control")),
        replicate=1,
        universe_policy=MDAUniversePolicy(),
    )
    settings = SASASettings(
        runs=[SASARunSettings(label="protein", target_selection="chainid A", stride=3)]
    )

    jobs = build_sasa_jobs(ctx, settings)

    assert len(jobs) == 1
    assert jobs[0].frame_selection.run_kwargs() == {"start": 5, "stop": 25, "step": 6}


def test_build_mda_jobs_subsamples_explicit_frames(monkeypatch: pytest.MonkeyPatch) -> None:
    """Explicit frame selectors should preserve non-contiguous stride semantics."""

    import polyzymd.analyses.sasa._mda as sasa_mda

    monkeypatch.setattr(sasa_mda, "build_sasa_surface_area_analysis", lambda **_kwargs: MagicMock())
    ctx = SimpleNamespace(
        universe=object(),
        frame_selection=FrameSelection(frames=(0, 3, 7, 11), timestep_ps=10.0, n_frames_total=20),
        replicate_context=SimpleNamespace(condition=_make_condition("control")),
        replicate=1,
        universe_policy=MDAUniversePolicy(),
    )
    settings = SASASettings(
        runs=[SASARunSettings(label="protein", target_selection="chainid A", stride=2)]
    )

    jobs = build_sasa_jobs(ctx, settings)

    assert jobs[0].frame_selection.run_kwargs() == {"frames": (0, 7)}


def test_sasa_analysisbase_wrapper_bounds_coordinate_buffer(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """SASA wrapper should flush coordinates at ``chunk_size`` during iteration."""

    import MDAnalysis as mda

    import polyzymd.analyses.sasa._mda as sasa_mda

    universe = mda.Universe.empty(
        n_atoms=1,
        n_residues=1,
        atom_resindex=[0],
        trajectory=True,
    )
    universe.add_TopologyAttr("name", ["C"])
    universe.add_TopologyAttr("type", ["C"])
    universe.add_TopologyAttr("element", ["C"])
    universe.add_TopologyAttr("resname", ["ALA"])
    universe.add_TopologyAttr("resid", [1])
    universe.add_TopologyAttr("chainID", ["A"])
    coordinates = np.zeros((5, 1, 3), dtype=np.float32)
    universe.load_new(coordinates, order="fac")
    chunk_sizes: list[int] = []

    def fake_compute_chunk(**kwargs):
        positions = kwargs["positions_angstrom"]
        chunk_sizes.append(int(positions.shape[0]))
        assert positions.shape[0] <= 2
        return {
            "atom_sasa_a2": np.ones((positions.shape[0], 1), dtype=np.float64),
            "residue_sasa_a2": np.ones((positions.shape[0], 1), dtype=np.float64),
            "total_sasa_a2": np.ones(positions.shape[0], dtype=np.float64),
        }

    monkeypatch.setattr(sasa_mda, "_mdtraj_topology_from_atom_group", lambda _atoms: object())
    monkeypatch.setattr(sasa_mda, "_compute_sasa_chunk_from_positions", fake_compute_chunk)

    analysis = sasa_mda.build_sasa_surface_area_analysis(
        universe=universe,
        run_label="protein",
        target_selection="all",
        context_selection="all",
        probe_radius_nm=0.14,
        n_sphere_points=100,
        chunk_size=2,
        timestep_ps=10.0,
    )

    analysis.run()

    assert chunk_sizes == [2, 2, 1]
    assert analysis.results.max_buffered_coordinate_frames == 2
    assert analysis.results.total_sasa_a2.tolist() == pytest.approx([1.0] * 5)


def test_sasa_collector_writes_npz_sidecar_only_for_arrays(tmp_path: Path) -> None:
    """Collector should keep large arrays in registered NPZ sidecars."""

    artifact = _collect_artifact(
        tmp_path,
        condition_label="control",
        replicate=1,
        raw=_raw_sasa(total=[10.0, 12.0, 14.0]),
    )

    assert isinstance(artifact, ReplicateArtifact)
    assert artifact.sidecars
    assert "atom_sasa_a2" not in artifact.payload
    assert "residue_sasa_a2" not in artifact.payload
    assert "total_sasa_a2" not in artifact.payload
    run_result = artifact.payload["run_results"][0]
    assert run_result["mean_sasa"] == pytest.approx(12.0)
    sidecar_path = ArtifactStore(tmp_path / "run_1").validate_sidecar(artifact.sidecars[0])
    with np.load(sidecar_path) as payload:
        assert payload["atom_sasa_a2"].shape == (3, 1)
        assert payload["residue_sasa_a2"].shape == (3, 1)
        assert payload["total_sasa_a2"].tolist() == pytest.approx([10.0, 12.0, 14.0])


def test_sasa_collector_preserves_zero_selection_behavior(tmp_path: Path) -> None:
    """Zero-selection runs should produce artifact-safe missing summaries."""

    artifact = _collect_artifact(
        tmp_path,
        condition_label="control",
        replicate=1,
        raw=_raw_sasa(total=[0.0, 0.0], zero=True),
    )

    run_result = artifact.payload["run_results"][0]
    assert run_result["zero_atom_selection"] is True
    assert run_result["mean_sasa"] is None
    assert run_result["missing_sasa_reason"] == "zero_atom_selection"
    assert run_result["n_target_atoms"] == 0
    assert ArtifactStore(tmp_path / "run_1").read_replicate_result("result.json")


def test_sasa_zero_selection_artifact_roundtrip_aggregate_compare(tmp_path: Path) -> None:
    """Zero-selection artifacts should round-trip and be skipped during compare."""

    settings = SASASettings(runs=[SASARunSettings(label="protein", target_selection="chainid Z")])
    analysis = SASAAnalysis()
    condition_artifacts = {}
    for condition_label in ("control", "treated"):
        condition_root = tmp_path / condition_label
        replicate_artifacts = [
            _collect_artifact(
                condition_root,
                condition_label=condition_label,
                replicate=replicate,
                raw=_raw_sasa(total=[0.0, 0.0], zero=True),
                settings=settings,
            )
            for replicate in (1, 2)
        ]
        reloaded_replicates = [
            ArtifactStore(condition_root / f"run_{replicate}").read_replicate_result("result.json")
            for replicate in (1, 2)
        ]
        aggregate = _aggregate_artifacts(
            condition_root,
            reloaded_replicates,
            condition_label=condition_label,
            settings=settings,
        )
        loaded_aggregate = ArtifactStore(condition_root / "aggregated").read_condition_result(
            "result.json"
        )
        condition_artifacts[condition_label] = loaded_aggregate

        replicate_json = (condition_root / "run_1" / "result.json").read_text()
        aggregate_json = (condition_root / "aggregated" / "result.json").read_text()
        assert "NaN" not in replicate_json
        assert "NaN" not in aggregate_json
        assert aggregate.payload["run_results"][0]["overall_mean"] is None
        assert aggregate.payload["run_results"][0]["per_replicate_means"] == [None, None]
        assert replicate_artifacts[0].payload["metrics"]["mean_sasa:protein"] is None

    ctx = make_comparison_context(
        name="sasa_compare",
        conditions=[_make_condition("control"), _make_condition("treated")],
        analysis_dirs={"control": tmp_path / "control", "treated": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="control",
        equilibration="10ns",
        aggregated_results=condition_artifacts,
    )

    assert analysis.compare(ctx) is None


def test_sasa_compare_mixed_finite_and_zero_selection_has_no_nan_leakage(
    tmp_path: Path,
) -> None:
    """Comparison output should omit zero-selection runs without serializing NaN."""

    settings = SASASettings(
        runs=[
            SASARunSettings(label="finite", target_selection="chainid A"),
            SASARunSettings(label="empty", target_selection="chainid Z"),
        ]
    )
    analysis = SASAAnalysis()
    condition_artifacts = {}
    for condition_label, offset in (("control", 0.0), ("treated", 10.0)):
        condition_root = tmp_path / condition_label
        for replicate in (1, 2):
            _collect_multi_run_artifact(
                condition_root,
                condition_label=condition_label,
                replicate=replicate,
                raw_by_run={
                    "finite": _raw_sasa(total=[10.0 + offset + replicate]),
                    "empty": _raw_sasa(total=[0.0], zero=True),
                },
                settings=settings,
            )
        reloaded_replicates = [
            ArtifactStore(condition_root / f"run_{replicate}").read_replicate_result("result.json")
            for replicate in (1, 2)
        ]
        condition_artifacts[condition_label] = _aggregate_artifacts(
            condition_root,
            reloaded_replicates,
            condition_label=condition_label,
            settings=settings,
        )

    ctx = make_comparison_context(
        name="sasa_compare",
        conditions=[_make_condition("control"), _make_condition("treated")],
        analysis_dirs={"control": tmp_path / "control", "treated": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="control",
        equilibration="10ns",
        aggregated_results=condition_artifacts,
    )

    comparison = analysis.compare(ctx)
    comparison_path = tmp_path / "comparison" / "result.json"
    comparison.save(comparison_path)
    comparison_json = comparison_path.read_text()

    assert isinstance(comparison, SASAComparisonResult)
    assert comparison.run_labels == ["finite"]
    assert all(
        [run.label for run in condition.run_summaries] == ["finite"]
        for condition in comparison.conditions
    )
    assert "NaN" not in comparison_json
    assert '"empty"' not in comparison_json


def test_sasa_sidecar_hash_validation_rejects_tampering(tmp_path: Path) -> None:
    """ArtifactStore validation should reject stale SASA sidecars."""

    artifact = _collect_artifact(
        tmp_path,
        condition_label="control",
        replicate=1,
        raw=_raw_sasa(total=[10.0, 12.0]),
    )
    sidecar_path = ArtifactStore(tmp_path / "run_1").resolve_sidecar(artifact.sidecars[0])
    sidecar_path.write_bytes(b"stale")

    with pytest.raises(ArtifactStoreError, match="SHA-256 mismatch|size mismatch"):
        ArtifactStore(tmp_path / "run_1").validate_sidecar(artifact.sidecars[0])


def test_sasa_aggregate_per_residue_profiles(tmp_path: Path) -> None:
    """Aggregation should average per-residue SASA from sidecars."""

    artifacts = [
        _collect_artifact(
            tmp_path,
            condition_label="control",
            replicate=1,
            raw=_raw_sasa(total=[10.0, 12.0], residue_offset=0.0),
        ),
        _collect_artifact(
            tmp_path,
            condition_label="control",
            replicate=2,
            raw=_raw_sasa(total=[14.0, 16.0], residue_offset=0.0),
        ),
    ]

    aggregated = _aggregate_artifacts(tmp_path, artifacts)

    assert isinstance(aggregated, ConditionArtifact)
    legacy = condition_artifact_to_legacy_result(aggregated)
    run = legacy.run_results[0]
    assert run.overall_mean == pytest.approx(13.0)
    assert run.per_replicate_means == pytest.approx([11.0, 15.0])
    assert run.per_residue_mean_sasa == pytest.approx([13.0])
    assert run.residue_keys == ["A:1:ALA"]


def test_sasa_aggregate_rejects_missing_sidecar(tmp_path: Path) -> None:
    """Aggregation should reject artifacts whose NPZ sidecar is missing."""

    artifacts = [
        _collect_artifact(
            tmp_path, condition_label="control", replicate=1, raw=_raw_sasa(total=[10.0])
        ),
        _collect_artifact(
            tmp_path, condition_label="control", replicate=2, raw=_raw_sasa(total=[11.0])
        ),
    ]
    ArtifactStore(tmp_path / "run_2").resolve_sidecar(artifacts[1].sidecars[0]).unlink()

    with pytest.raises(ValueError, match="sidecar"):
        _aggregate_artifacts(tmp_path, artifacts)


def test_sasa_aggregate_rejects_residue_identity_mismatch(tmp_path: Path) -> None:
    """Aggregation should reject mismatched residue identity/order."""

    artifacts = [
        _collect_artifact(
            tmp_path,
            condition_label="control",
            replicate=1,
            raw=_raw_sasa(total=[10.0], residue_key="A:1:ALA"),
        ),
        _collect_artifact(
            tmp_path,
            condition_label="control",
            replicate=2,
            raw=_raw_sasa(total=[11.0], residue_key="A:2:GLY"),
        ),
    ]

    with pytest.raises(ValueError, match="residue metadata mismatch"):
        _aggregate_artifacts(tmp_path, artifacts)


def test_sasa_aggregate_rejects_legacy_inputs(tmp_path: Path) -> None:
    """Aggregation should reject non-artifact legacy replicate results."""

    analysis = SASAAnalysis()
    settings = SASASettings(runs=[SASARunSettings(label="protein", target_selection="chainid A")])
    ctx = make_aggregate_context(
        condition=_make_condition("control"),
        replicates=(1,),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )

    with pytest.raises(TypeError, match="ReplicateArtifact"):
        analysis.aggregate(ctx, [MagicMock()])


def test_sasa_compare_and_format_with_condition_artifacts(tmp_path: Path) -> None:
    """Comparison and formatting should preserve SASA run-wise output."""

    settings = SASASettings(runs=[SASARunSettings(label="protein", target_selection="chainid A")])
    analysis = SASAAnalysis()
    control = _aggregate_artifacts(
        tmp_path / "control",
        [
            _collect_artifact(
                tmp_path / "control",
                condition_label="control",
                replicate=1,
                raw=_raw_sasa(total=[10.0]),
            ),
            _collect_artifact(
                tmp_path / "control",
                condition_label="control",
                replicate=2,
                raw=_raw_sasa(total=[11.0]),
            ),
        ],
        condition_label="control",
        settings=settings,
    )
    treated = _aggregate_artifacts(
        tmp_path / "treated",
        [
            _collect_artifact(
                tmp_path / "treated",
                condition_label="treated",
                replicate=1,
                raw=_raw_sasa(total=[20.0]),
            ),
            _collect_artifact(
                tmp_path / "treated",
                condition_label="treated",
                replicate=2,
                raw=_raw_sasa(total=[21.0]),
            ),
        ],
        condition_label="treated",
        settings=settings,
    )
    ctx = make_comparison_context(
        name="sasa_compare",
        conditions=[_make_condition("control"), _make_condition("treated")],
        analysis_dirs={"control": tmp_path / "control", "treated": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="control",
        equilibration="10ns",
        aggregated_results={"control": control, "treated": treated},
    )

    comparison = analysis.compare(ctx)

    assert isinstance(comparison, SASAComparisonResult)
    assert comparison.run_labels == ["protein"]
    assert "SASA Comparison" in analysis.format(comparison, "table")


def test_sasa_compare_rejects_legacy_aggregate(tmp_path: Path) -> None:
    """Comparison should fail loudly for legacy aggregate inputs."""

    analysis = SASAAnalysis()
    settings = SASASettings(runs=[SASARunSettings(label="protein", target_selection="chainid A")])
    ctx = make_comparison_context(
        name="sasa_compare",
        conditions=[_make_condition("control")],
        analysis_dirs={"control": tmp_path / "control"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="control",
        equilibration="10ns",
        aggregated_results={"control": MagicMock()},
    )

    with pytest.raises(TypeError, match="canonical MDAnalysis condition artifacts"):
        analysis.compare(ctx)


def test_sasa_timeseries_loader_uses_canonical_artifact_sidecars(tmp_path: Path) -> None:
    """Timeseries plot data should load canonical artifact sidecars only."""

    artifact = _collect_artifact(
        tmp_path / "condition",
        condition_label="condition",
        replicate=1,
        raw=_raw_sasa(total=[10.0, 12.0, 11.0]),
    )

    time_ns, matrix = _load_replicate_timeseries_from_results(
        tmp_path / "condition", "protein", [1]
    )
    loaded = load_replicate_sasa_sidecar(artifact, "protein", tmp_path / "condition" / "run_1")

    assert loaded.total_sasa_a2.tolist() == pytest.approx([10.0, 12.0, 11.0])
    assert time_ns.size == 3
    assert matrix.shape == (1, 3)


def test_plot_loader_ignores_legacy_sasa_json_without_result(tmp_path: Path) -> None:
    """Plot payload loader should not scan legacy SASA JSON files."""

    run_dir = tmp_path / "condition" / "run_1"
    run_dir.mkdir(parents=True)
    (run_dir / "sasa_legacy.json").write_text('{"run_results": [{"run_label": "protein"}]}')

    assert _load_condition_result_payloads(tmp_path / "condition") == []


def test_normalized_control_rows_compute_percent_delta(tmp_path: Path) -> None:
    """Normalized-control rows should report percent changes from control."""

    comparison = _make_sasa_comparison_result(
        {
            "control": {"protein": (10.0, 0.0)},
            "condition_5": {"protein": (5.0, 0.0)},
            "condition_8": {"protein": (8.0, 0.0)},
        },
        control_label=None,
    )
    ctx = _make_sasa_plot_context(tmp_path)

    rows = _build_sasa_normalized_control_rows(ctx, comparison, "protein")

    assert [row.condition_label for row in rows] == ["condition_5", "condition_8"]
    assert [row.percent_delta for row in rows] == pytest.approx([-50.0, -20.0])


def test_normalized_control_rows_run_replicate_deltas(tmp_path: Path) -> None:
    """Normalized-control rows should preserve per-replicate percent changes."""

    comparison = _make_sasa_comparison_result(
        {"control": {"protein": (10.0, 0.0)}, "condition_5": {"protein": (5.0, 0.0)}}
    )
    comparison.conditions[1].run_summaries[0].per_replicate_means = [4.0, 5.0, 6.0]
    ctx = _make_sasa_plot_context(tmp_path)

    rows = _build_sasa_normalized_control_rows(ctx, comparison, "protein")

    assert rows[0].replicate_percent_deltas == pytest.approx([-60.0, -50.0, -40.0])


def test_normalized_replicate_deltas_use_control_mean_baseline() -> None:
    """Control replicates should normalize against the aggregate control mean."""

    deltas = _build_sasa_normalized_replicate_deltas([9.0, 10.0, 11.0], 10.0)

    assert deltas == pytest.approx([-10.0, 0.0, 10.0])


def test_normalized_control_sem_propagation() -> None:
    """Normalized-control SEM should follow first-order error propagation."""

    sem_delta = _propagate_sasa_normalized_sem(
        condition_mean=5.0,
        condition_sem=0.5,
        control_mean=10.0,
        control_sem=1.0,
    )

    expected = 100.0 * np.sqrt((0.5 / 10.0) ** 2 + (5.0 * 1.0 / 10.0**2) ** 2)
    assert sem_delta == pytest.approx(expected)


def test_normalized_control_plot_skips_missing_control(tmp_path: Path) -> None:
    """Normalized-control plotter should return no paths when control is unavailable."""

    comparison = _make_sasa_comparison_result(
        {"control": {"protein": (10.0, 0.1)}, "treated": {"protein": (5.0, 0.1)}},
        control_label="missing",
    )
    ctx = _make_sasa_plot_context(tmp_path, control_label=None)

    assert plot_sasa_normalized_control_bars(ctx, comparison) == []


@pytest.mark.parametrize("control_mean", [0.0, float("nan"), float("inf")])
def test_normalized_control_rows_skip_invalid_control_mean(
    tmp_path: Path,
    control_mean: float,
) -> None:
    """Normalized-control rows should skip runs with invalid control means."""

    comparison = _make_sasa_comparison_result(
        {"control": {"protein": (control_mean, 0.1)}, "treated": {"protein": (5.0, 0.1)}}
    )
    ctx = _make_sasa_plot_context(tmp_path)

    assert _build_sasa_normalized_control_rows(ctx, comparison, "protein") == []


def test_normalized_control_plot_generates_one_path_per_valid_run(tmp_path: Path) -> None:
    """Normalized-control plotter should generate one plot for each valid run label."""

    comparison = _make_sasa_comparison_result(
        {
            "control": {"protein": (10.0, 0.1), "active site": (20.0, 0.2)},
            "treated": {"protein": (5.0, 0.1), "active site": (10.0, 0.2)},
        }
    )
    ctx = _make_sasa_plot_context(tmp_path)

    paths = plot_sasa_normalized_control_bars(ctx, comparison)

    assert [path.name for path in paths] == [
        "sasa_normalized_comparison_protein.png",
        "sasa_normalized_comparison_active_site.png",
    ]
    assert all(path.exists() for path in paths)


def test_sasa_analysis_plot_includes_normalized_control_plotter(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """SASAAnalysis.plot should append normalized-control paths to existing plots."""

    import polyzymd.analyses.sasa._plotters as plotters

    analysis = SASAAnalysis()
    comparison = _make_sasa_comparison_result(
        {"control": {"protein": (10.0, 0.1)}, "treated": {"protein": (5.0, 0.1)}}
    )
    comparison_path = tmp_path / "results" / "result.json"
    comparison.save(comparison_path)
    ctx = _make_sasa_plot_context(tmp_path)
    ctx = PlotContext(
        conditions=ctx.conditions,
        analysis_dirs=ctx.analysis_dirs,
        results_dir=ctx.results_dir,
        output_dir=ctx.output_dir,
        settings=ctx.settings,
        plot_settings=ctx.plot_settings,
        comparison_path=comparison_path,
        control_label=ctx.control_label,
        equilibration=ctx.equilibration,
    )
    expected_paths = [
        tmp_path / "figures" / "bars.png",
        tmp_path / "figures" / "timeseries.png",
        tmp_path / "figures" / "profiles.png",
        tmp_path / "figures" / "normalized.png",
    ]

    monkeypatch.setattr(
        plotters, "plot_sasa_comparison_bars", lambda _ctx, _result: [expected_paths[0]]
    )
    monkeypatch.setattr(plotters, "plot_sasa_timeseries", lambda _ctx, _result: [expected_paths[1]])
    monkeypatch.setattr(
        plotters, "plot_sasa_residue_profiles", lambda _ctx, _result: [expected_paths[2]]
    )
    monkeypatch.setattr(
        plotters,
        "plot_sasa_normalized_control_bars",
        lambda _ctx, _result: [expected_paths[3]],
    )

    assert analysis.plot(ctx) == expected_paths


def test_plot_helper_sanitize_label() -> None:
    """Plot helper should normalize run labels for file names."""

    assert _sanitize_run_label("Protein/Core Run") == "protein_core_run"


def test_sasa_plot_settings_model_attribute() -> None:
    """SASA analysis should expose its plot settings model."""

    from polyzymd.analyses.discovery import get_analysis

    analysis_cls = get_analysis("sasa")
    assert analysis_cls.PlotSettingsModel is SASAPlotSettings


def test_sasa_singleton_pairwise_not_testable() -> None:
    """SASA singleton pairwise results should carry not-testable metadata."""

    from polyzymd.analyses.shared.inferential_statistics import (
        cohens_d,
        independent_ttest,
        percent_change,
    )

    run_a = SASARunSummary(
        label="protein",
        target_selection="chainid A",
        context_selection="chainid A",
        mean_sasa=100.0,
        sem_sasa=0.0,
        per_replicate_means=[100.0],
    )
    run_b = SASARunSummary(
        label="protein",
        target_selection="chainid A",
        context_selection="chainid A",
        mean_sasa=120.0,
        sem_sasa=0.0,
        per_replicate_means=[120.0],
    )

    comparison = SASAAnalysis._compare_run(
        run_label="protein",
        condition_a="A",
        condition_b="B",
        run_a=run_a,
        run_b=run_b,
        independent_ttest=independent_ttest,
        cohens_d=cohens_d,
        percent_change=percent_change,
    )

    assert comparison.testable is False
    assert comparison.significant is False
    assert comparison.p_value is None
    assert comparison.p_value_adjusted is None
    assert comparison.note is not None


def test_sasa_formatter_singleton_sem_rendered_as_not_available() -> None:
    """SASA formatter should not show singleton SEM as a numeric value."""

    result = SASAComparisonResult(
        metric="mean_sasa",
        name="singleton_sasa",
        n_runs=1,
        run_labels=["protein"],
        conditions=[
            SASAConditionSummary(
                label="control",
                config_path="/fake/control.yaml",
                n_replicates=1,
                run_summaries=[
                    SASARunSummary(
                        label="protein",
                        target_selection="chainid A",
                        context_selection="all",
                        mean_sasa=100.0,
                        sem_sasa=0.0,
                        per_replicate_means=[100.0],
                    )
                ],
            )
        ],
        pairwise_comparisons=[],
        ranking_by_run={"protein": ["control"]},
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.0.0",
    )

    text_output = format_sasa_comparison(result, "table")
    markdown_output = format_sasa_comparison(result, "markdown")

    assert "n/a" in text_output
    assert "SEM: n/a (single replicate; not estimable)" in text_output
    assert "| control | 100.00 | n/a | 1 |" in markdown_output


def test_comparison_result_models() -> None:
    """SASA comparison models should instantiate and expose helpers."""

    condition = SASAConditionSummary(
        label="Control",
        config_path="/fake/control.yaml",
        n_replicates=3,
        run_summaries=[
            SASARunSummary(
                label="protein",
                target_selection="chainid A",
                context_selection="all",
                mean_sasa=100.0,
                sem_sasa=1.0,
                per_replicate_means=[99.0, 100.0, 101.0],
            )
        ],
    )
    comp = SASAComparisonResult(
        metric="mean_sasa",
        name="proj",
        n_runs=1,
        run_labels=["protein"],
        control_label="Control",
        conditions=[condition],
        pairwise_comparisons=[
            SASARunPairwiseComparison(
                run_label="protein",
                condition_a="Control",
                condition_b="Treatment",
                t_statistic=2.0,
                p_value=0.04,
                cohens_d=1.0,
                effect_interpretation="large",
                direction="exposure",
                significant=True,
                percent_change=10.0,
            )
        ],
        anova_by_run=None,
        ranking_by_run={"protein": ["Control"]},
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.0.0",
    )
    assert condition.get_run("protein").mean_sasa == pytest.approx(100.0)
    assert comp.get_comparisons_for_run("protein")
    assert "SASA Comparison" in format_sasa_comparison(comp, "table")


def test_run_cache_token_changes_when_stride_changes() -> None:
    """Run cache token should differ when stride changes."""

    token_stride_1 = SASAAnalysis._run_cache_token(
        label="protein",
        target_selection="chainid A",
        context_selection="chainid A",
        probe_radius_nm=0.14,
        n_sphere_points=960,
        stride=1,
        equilibration="10ns",
    )
    token_stride_2 = SASAAnalysis._run_cache_token(
        label="protein",
        target_selection="chainid A",
        context_selection="chainid A",
        probe_radius_nm=0.14,
        n_sphere_points=960,
        stride=2,
        equilibration="10ns",
    )

    assert token_stride_1 != token_stride_2


def test_plot_loader_fallback_tie_breaks_are_removed(tmp_path: Path) -> None:
    """Legacy fallback files should be ignored regardless of mtimes."""

    run_dir = tmp_path / "condition" / "run_1"
    run_dir.mkdir(parents=True)
    older_name_path = run_dir / "sasa_aaa.json"
    newer_name_path = run_dir / "sasa_bbb.json"
    older_name_path.write_text('{"run_results": [{"source": "aaa"}]}', encoding="utf-8")
    newer_name_path.write_text('{"run_results": [{"source": "bbb"}]}', encoding="utf-8")
    same_mtime = 1_700_000_000
    os.utime(older_name_path, (same_mtime, same_mtime))
    os.utime(newer_name_path, (same_mtime, same_mtime))

    assert _load_condition_result_payloads(tmp_path / "condition") == []
