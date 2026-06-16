"""Tests for MDAnalysis artifact collector interfaces."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import Condition, ReplicateContext
from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.mda import (
    ArtifactStore,
    FrameSelection,
    MDACollectorContext,
    MDAJobResult,
    MDAUniversePolicy,
    ReplicateArtifact,
    StrictJSONMDAResultCollector,
    frame_selection_payload,
    strict_json_payload,
)
from polyzymd.analyses.mda.job import MDABackendPolicy


class _Settings(BaseModel):
    """Minimal settings model for collector context tests."""

    scale: float = 1.0


class FakeMDAnalysisResults(dict):
    """Import-light stand-in for ``MDAnalysis.analysis.results.Results``."""

    __module__ = "MDAnalysis.analysis.results"


def _collector_context(
    tmp_path: Path,
    frame_selection: FrameSelection | None = None,
) -> MDACollectorContext:
    """Build a minimal collector context.

    Parameters
    ----------
    tmp_path : Path
        Temporary artifact root.
    frame_selection : FrameSelection or None, optional
        Frame selection to attach to the collector context.

    Returns
    -------
    MDACollectorContext
        Context suitable for collector unit tests.
    """

    condition = Condition("Cond", tmp_path / "cond.yaml", (1,), SimpleNamespace())
    settings = _Settings(scale=2.0)
    replicate_context = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="10ns",
        recompute=False,
        settings=settings,
        result_path=tmp_path / "run_1" / "result.json",
    )
    resolved_frame_selection = frame_selection or FrameSelection(
        start=1,
        stop=5,
        step=2,
        n_frames_total=8,
        first_frame_time_ps=198_400.0,
        selected_start_time_ps=198_800.0,
        equilibration_time_reference="trajectory_timestamp",
    )
    return MDACollectorContext(
        analysis_name="example",
        replicate_context=replicate_context,
        frame_selection=resolved_frame_selection,
        universe_policy=MDAUniversePolicy(condition_label="Cond", replicate=1),
        artifact_store=ArtifactStore(tmp_path / "run_1"),
        settings_fingerprint="abc123",
        warnings=["first warning"],
    )


def _job_result(results: object) -> MDAJobResult:
    """Build a completed job reference.

    Parameters
    ----------
    results : object
        Results payload to expose from the fake job.

    Returns
    -------
    MDAJobResult
        Completed job result reference.
    """

    frame_selection = FrameSelection(start=1, stop=5, step=2, n_frames_total=8)
    return MDAJobResult(
        name="fake_job",
        analysis=SimpleNamespace(results=results),
        results=results,
        run_kwargs={"start": 1, "stop": 5, "step": 2},
        frame_selection=frame_selection,
        backend_policy=MDABackendPolicy(),
        universe_policy=MDAUniversePolicy(condition_label="Cond", replicate=1),
    )


def test_collector_context_exposes_framework_shortcuts(tmp_path: Path) -> None:
    """Collector context should expose common framework identity fields."""

    ctx = _collector_context(tmp_path)

    assert ctx.condition_label == "Cond"
    assert ctx.replicate == 1
    assert ctx.output_dir == tmp_path / "run_1"
    assert ctx.result_path == tmp_path / "run_1" / "result.json"
    assert ctx.settings.scale == 2.0
    assert ctx.warnings == ("first warning",)


def test_strict_json_collector_emits_replicate_artifact(tmp_path: Path) -> None:
    """Default collector should preserve JSON-safe simple job results."""

    ctx = _collector_context(tmp_path)
    collector = StrictJSONMDAResultCollector()

    artifact = collector(ctx, [_job_result({"value": 1.25})])

    assert isinstance(artifact, ReplicateArtifact)
    assert artifact.analysis_name == "example"
    assert artifact.condition_label == "Cond"
    assert artifact.replicate == 1
    assert artifact.payload["jobs"][0]["results"] == {"value": 1.25}
    assert artifact.metadata["settings_fingerprint"] == "abc123"
    assert artifact.warnings == ["first warning"]


def test_frame_selection_payload_is_json_safe(tmp_path: Path) -> None:
    """Frame-selection payload should preserve run-window provenance."""

    ctx = _collector_context(tmp_path)

    payload = frame_selection_payload(ctx.frame_selection)

    assert payload["start"] == 1
    assert payload["stop"] == 5
    assert payload["step"] == 2
    assert payload["n_frames_selected"] == 2
    assert payload["first_frame_time_ps"] == 198_400.0
    assert payload["selected_start_time_ps"] == 198_800.0
    assert payload["equilibration_time_reference"] == "trajectory_timestamp"


def test_frame_selection_payload_serializes_numpy_bool_masks() -> None:
    """Frame-selection payload should serialize NumPy boolean masks as booleans."""

    selection = FrameSelection(
        frames=np.asarray([True, False, True], dtype=np.bool_),
        n_frames_total=3,
    )

    payload = frame_selection_payload(selection)

    assert payload["frames"] == [True, False, True]
    assert all(isinstance(value, bool) for value in payload["frames"])


def test_frame_selection_payload_serializes_numpy_slice_scalars() -> None:
    """Frame-selection payload should serialize NumPy slice scalars as Python values."""

    selection = FrameSelection(
        start=np.int64(1),
        stop=np.int64(6),
        step=np.int64(2),
        equilibration_start=np.int64(1),
        equilibration_ps=np.float32(100.0),
        timestep_ps=np.float64(2.0),
        first_frame_time_ps=np.float64(10.0),
        selected_start_time_ps=np.float32(12.0),
        n_frames_total=np.int64(8),
    )

    payload = frame_selection_payload(selection)

    assert payload["start"] == 1
    assert payload["stop"] == 6
    assert payload["step"] == 2
    assert payload["equilibration_start"] == 1
    assert payload["equilibration_ps"] == pytest.approx(100.0)
    assert payload["timestep_ps"] == pytest.approx(2.0)
    assert payload["first_frame_time_ps"] == pytest.approx(10.0)
    assert payload["selected_start_time_ps"] == pytest.approx(12.0)
    assert payload["n_frames_total"] == 8
    assert payload["n_frames_selected"] == 3
    numeric_values = [
        payload["start"],
        payload["stop"],
        payload["step"],
        payload["equilibration_start"],
        payload["equilibration_ps"],
        payload["timestep_ps"],
        payload["first_frame_time_ps"],
        payload["selected_start_time_ps"],
        payload["n_frames_total"],
        payload["n_frames_selected"],
    ]
    assert all(not hasattr(value, "item") for value in numeric_values)


def test_strict_json_collector_serializes_numpy_bool_run_kwargs(tmp_path: Path) -> None:
    """Default collector should emit JSON-safe run kwargs for NumPy boolean masks."""

    ctx = _collector_context(tmp_path)
    frame_selection = FrameSelection(
        frames=np.asarray([True, False, True], dtype=np.bool_),
        n_frames_total=3,
    )
    job = MDAJobResult(
        name="mask_job",
        analysis=SimpleNamespace(results={"value": 1.0}),
        results={"value": 1.0},
        run_kwargs=frame_selection.run_kwargs(),
        frame_selection=frame_selection,
        backend_policy=MDABackendPolicy(),
        universe_policy=MDAUniversePolicy(condition_label="Cond", replicate=1),
    )

    artifact = StrictJSONMDAResultCollector()(ctx, [job])
    run_frames = artifact.payload["jobs"][0]["run_kwargs"]["frames"]

    json.dumps(artifact.model_dump(mode="json"))
    assert run_frames == [True, False, True]
    assert all(isinstance(value, bool) for value in run_frames)


def test_strict_json_collector_serializes_numpy_slice_run_kwargs(tmp_path: Path) -> None:
    """Default collector should emit JSON-safe slice kwargs from NumPy scalars."""

    frame_selection = FrameSelection(
        start=np.int64(1),
        stop=np.int64(6),
        step=np.int64(2),
        equilibration_start=np.int64(1),
        equilibration_ps=np.float32(100.0),
        timestep_ps=np.float64(2.0),
        first_frame_time_ps=np.float64(10.0),
        selected_start_time_ps=np.float32(12.0),
        n_frames_total=np.int64(8),
    )
    ctx = _collector_context(tmp_path, frame_selection=frame_selection)
    job = MDAJobResult(
        name="slice_job",
        analysis=SimpleNamespace(results={"value": 1.0}),
        results={"value": 1.0},
        run_kwargs=frame_selection.run_kwargs(),
        frame_selection=frame_selection,
        backend_policy=MDABackendPolicy(),
        universe_policy=MDAUniversePolicy(condition_label="Cond", replicate=1),
    )

    artifact = StrictJSONMDAResultCollector()(ctx, [job])
    run_kwargs = artifact.payload["jobs"][0]["run_kwargs"]
    job_selection = artifact.payload["jobs"][0]["frame_selection"]
    provenance_selection = artifact.provenance["frame_selection"]

    json.dumps(artifact.model_dump(mode="json"))
    assert run_kwargs == {"start": 1, "stop": 6, "step": 2}
    assert job_selection["equilibration_start"] == 1
    assert job_selection["equilibration_ps"] == pytest.approx(100.0)
    assert provenance_selection["n_frames_total"] == 8
    assert provenance_selection["n_frames_selected"] == 3
    numeric_values = [
        *run_kwargs.values(),
        job_selection["start"],
        job_selection["stop"],
        job_selection["step"],
        job_selection["equilibration_start"],
        job_selection["equilibration_ps"],
        provenance_selection["first_frame_time_ps"],
        provenance_selection["selected_start_time_ps"],
        provenance_selection["n_frames_total"],
        provenance_selection["n_frames_selected"],
    ]
    assert all(not hasattr(value, "item") for value in numeric_values)


def test_strict_json_payload_rejects_raw_mdanalysis_results() -> None:
    """Strict JSON conversion should require custom collectors for raw Results."""

    with pytest.raises(PluginContractError, match="raw MDAnalysis Results"):
        strict_json_payload({"raw": FakeMDAnalysisResults(value=1.0)}, analysis_name="example")


def test_strict_json_collector_rejects_raw_mdanalysis_results(tmp_path: Path) -> None:
    """Default collector should not serialize raw MDAnalysis Results."""

    ctx = _collector_context(tmp_path)
    collector = StrictJSONMDAResultCollector()

    with pytest.raises(PluginContractError, match="raw MDAnalysis Results"):
        collector(ctx, [_job_result(FakeMDAnalysisResults(value=1.0))])
