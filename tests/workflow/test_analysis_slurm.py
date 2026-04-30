"""Tests for analysis SLURM DAG helpers."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar, cast
from unittest.mock import MagicMock

import pytest
from pydantic import BaseModel, ValidationError

from polyzymd.analyses.base import Analysis, Condition
from polyzymd.workflow.analysis_slurm import (
    AnalysisJobManifest,
    AnalysisSlurmResources,
    ConditionTaskSpec,
    ReplicateTaskSpec,
    SubmittedJobGraph,
    _cancel_jobs,
    _map_slurm_state,
    _query_sacct,
    _sanitize_path_for_script,
    _sanitize_slurm_value,
    _submit_sbatch,
    build_manifest,
    generate_aggregate_script,
    generate_array_script,
    generate_finalize_script,
    generate_replicate_script,
    read_analysis_status,
    reconcile_status_with_slurm,
    submit_analysis_graph,
    submit_analysis_graph_with_arrays,
    update_task_status,
)


class _Settings(BaseModel):
    threshold: float = 1.0


class _ToyRunner:
    """Minimal runner used to satisfy the analysis compute-stage contract."""

    def __init__(self, replicate: int) -> None:
        """Initialize an empty result container for one replicate.

        Parameters
        ----------
        replicate : int
            Replicate identifier associated with this runner.
        """

        self.replicate = replicate
        self.results: dict[str, Any] = {}

    def run(self, start: int, stop: int, step: int = 1) -> "_ToyRunner":
        """Record the requested frame window and return the executed runner.

        Parameters
        ----------
        start : int
            Inclusive start frame.
        stop : int
            Exclusive stop frame.
        step : int, optional
            Frame stride, by default 1.

        Returns
        -------
        _ToyRunner
            Runner with populated ``results``.
        """

        self.results = {
            "replicate": self.replicate,
            "start": start,
            "stop": stop,
            "step": step,
        }
        return self


class _ToyAnalysis(Analysis):
    """Runner-backed toy analysis used by SLURM manifest tests."""

    name: ClassVar[str] = "toy_slurm"
    Settings: ClassVar[type] = _Settings

    def build_runner(self, ctx: Any, replicate: int, universe: Any, window: Any) -> _ToyRunner:
        """Build a minimal runner for the requested replicate.

        Parameters
        ----------
        ctx : Any
            Framework-provided replicate context.
        replicate : int
            Replicate identifier.
        universe : Any
            Loaded trajectory universe.
        window : Any
            Resolved trajectory window.

        Returns
        -------
        _ToyRunner
            Minimal runner exposing ``run()`` and ``results``.
        """

        del ctx, universe, window
        return _ToyRunner(replicate)

    def summarize_replicate(
        self,
        ctx: Any,
        replicate: int,
        runner: _ToyRunner,
        window: Any,
    ) -> dict[str, Any]:
        """Convert runner output into a serializable replicate result.

        Parameters
        ----------
        ctx : Any
            Framework-provided replicate context.
        replicate : int
            Replicate identifier.
        runner : _ToyRunner
            Executed runner with populated results.
        window : Any
            Resolved trajectory window.

        Returns
        -------
        dict[str, Any]
            Serialized replicate summary.
        """

        del window
        return {
            "replicate": replicate,
            "threshold": ctx.settings.threshold,
            "runner_results": runner.results,
        }

    def aggregate(self, ctx: Any, results: list[Any]) -> dict[str, Any]:
        """Aggregate toy replicate results for full-stage workflow tests.

        Parameters
        ----------
        ctx : Any
            Framework-provided aggregate context.
        results : list[Any]
            Replicate summaries.

        Returns
        -------
        dict[str, Any]
            Minimal aggregate summary.
        """

        del ctx
        return {
            "n_replicates": len(results),
            "replicates": [result["replicate"] for result in results],
        }


class _CompareOnlyAnalysis(Analysis):
    name: ClassVar[str] = "compare_only"
    Settings: ClassVar[type] = _Settings
    has_compute_stage: ClassVar[bool] = False
    has_aggregate_stage: ClassVar[bool] = False


def _make_manifest(tmp_path: Path) -> AnalysisJobManifest:
    resources = AnalysisSlurmResources()
    return AnalysisJobManifest(
        analysis_name="toy_slurm",
        comparison_yaml=str(tmp_path / "comparison.yaml"),
        condition_specs=[
            ConditionTaskSpec(
                condition_index=0,
                condition_label="Cond A",
                condition_slug="cond_a",
                replicate_specs=[
                    ReplicateTaskSpec(
                        condition_index=0,
                        replicate=1,
                        condition_label="Cond A",
                        condition_slug="cond_a",
                    ),
                    ReplicateTaskSpec(
                        condition_index=0,
                        replicate=2,
                        condition_label="Cond A",
                        condition_slug="cond_a",
                    ),
                ],
            )
        ],
        settings_snapshot={"threshold": 1.0},
        equilibration="10ns",
        recompute=False,
        resources=resources,
        created_at="2026-01-01T00:00:00+00:00",
    )


def _make_manifest_two_conditions(tmp_path: Path) -> AnalysisJobManifest:
    resources = AnalysisSlurmResources()
    return AnalysisJobManifest(
        analysis_name="toy_slurm",
        comparison_yaml=str(tmp_path / "comparison.yaml"),
        condition_specs=[
            ConditionTaskSpec(
                condition_index=0,
                condition_label="Cond A",
                condition_slug="cond_a",
                replicate_specs=[
                    ReplicateTaskSpec(
                        condition_index=0,
                        replicate=1,
                        condition_label="Cond A",
                        condition_slug="cond_a",
                    )
                ],
            ),
            ConditionTaskSpec(
                condition_index=1,
                condition_label="Cond B",
                condition_slug="cond_b",
                replicate_specs=[
                    ReplicateTaskSpec(
                        condition_index=1,
                        replicate=1,
                        condition_label="Cond B",
                        condition_slug="cond_b",
                    )
                ],
            ),
        ],
        settings_snapshot={"threshold": 1.0},
        equilibration="10ns",
        recompute=False,
        resources=resources,
        created_at="2026-01-01T00:00:00+00:00",
    )


def test_build_manifest_serialization(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Manifest should be built and round-trip serialized."""
    analysis = _ToyAnalysis()
    resources = AnalysisSlurmResources(partition="gpu")
    valid_conditions = [
        Condition("A", tmp_path / "a.yaml", (1, 2), cast(Any, SimpleNamespace())),
        Condition("B", tmp_path / "b.yaml", (3,), cast(Any, SimpleNamespace())),
    ]

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.prepare_comparison_run",
        lambda analysis, config, equilibration: {
            "all_conditions": valid_conditions,
            "valid_conditions": valid_conditions,
            "excluded_conditions": [],
            "condition_by_label": {condition.label: condition for condition in valid_conditions},
            "settings": _Settings(),
            "equilibration": "5ns",
            "analysis_root": tmp_path / "analysis",
        },
    )
    config = cast(Any, SimpleNamespace(source_path=tmp_path / "comparison.yaml"))

    manifest = build_manifest(analysis, config, resources, recompute=True, equilibration=None)
    assert manifest.analysis_name == "toy_slurm"
    assert manifest.condition_specs[0].condition_slug == "A"
    assert manifest.condition_specs[1].replicate_specs[0].replicate == 3

    out = tmp_path / "manifest.json"
    manifest.save(out)
    loaded = AnalysisJobManifest.load(out)
    assert loaded.analysis_name == manifest.analysis_name
    assert loaded.recompute is True
    assert loaded.partial_policy == "strict"


def test_script_generation_contains_worker_commands(tmp_path: Path) -> None:
    """Generated scripts should call worker CLI commands."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources(max_retries=3)
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"

    rep_script = generate_replicate_script(
        manifest,
        manifest.condition_specs[0].replicate_specs[0],
        resources,
        hpc_dir,
    )
    agg_script = generate_aggregate_script(
        manifest, manifest.condition_specs[0], resources, hpc_dir
    )
    fin_script = generate_finalize_script(manifest, resources, hpc_dir)

    rep_text = rep_script.read_text()
    assert "#SBATCH --requeue" in rep_text
    assert "worker-replicate" in rep_text
    assert "scontrol requeue" in rep_text
    assert "run -e build python -c" in rep_text
    assert "python3 -c" not in rep_text

    assert "worker-aggregate" in agg_script.read_text()
    assert "worker-finalize" in fin_script.read_text()


def test_slurm_header_omits_partition_when_unset(tmp_path: Path) -> None:
    """Generated scripts should omit partition SBATCH line when unset."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources(partition=None)
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"

    script = generate_finalize_script(manifest, resources, hpc_dir)
    text = script.read_text()
    assert "#SBATCH --partition=" not in text


def test_generate_array_script_contains_array_directive(tmp_path: Path) -> None:
    """Array scripts should include explicit array task specification."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources(max_retries=3)
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    cond_spec = manifest.condition_specs[0]

    script = generate_array_script(cond_spec, manifest, resources, [1, 2, 3, 4, 5], hpc_dir)
    text = script.read_text()

    assert "#SBATCH --array=1,2,3,4,5" in text


def test_generate_array_script_uses_slurm_array_task_id(tmp_path: Path) -> None:
    """Array scripts should dispatch with SLURM_ARRAY_TASK_ID."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources(max_retries=3)
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    cond_spec = manifest.condition_specs[0]

    script = generate_array_script(cond_spec, manifest, resources, [1, 2], hpc_dir)
    text = script.read_text()

    assert "$SLURM_ARRAY_TASK_ID" in text


def test_generate_array_script_requeues_only_failing_task(tmp_path: Path) -> None:
    """Array retries should requeue only the failing array element."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources(max_retries=3)
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    cond_spec = manifest.condition_specs[0]

    script = generate_array_script(cond_spec, manifest, resources, [1, 2], hpc_dir)
    text = script.read_text()

    assert "scontrol requeue ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" in text
    assert "scontrol requeue $SLURM_ARRAY_JOB_ID" not in text


def test_generate_array_script_uses_condition_index_dispatch(tmp_path: Path) -> None:
    """Array workers should dispatch by manifest condition index."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources(max_retries=3)
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    cond_spec = manifest.condition_specs[0]

    script = generate_array_script(cond_spec, manifest, resources, [1, 2], hpc_dir)
    text = script.read_text()

    assert "--condition-index 0" in text
    assert '--condition "' not in text


def test_generate_array_script_non_contiguous_replicates(tmp_path: Path) -> None:
    """Array scripts should preserve non-contiguous replicate IDs."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources(max_retries=3)
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    cond_spec = manifest.condition_specs[0]

    script = generate_array_script(cond_spec, manifest, resources, [1, 5, 3], hpc_dir)
    text = script.read_text()

    assert "#SBATCH --array=1,3,5" in text


def test_build_manifest_allow_partial_policy(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Manifest should store allow_partial policy from submission."""
    analysis = _ToyAnalysis()
    resources = AnalysisSlurmResources(partition="gpu")
    valid_conditions = [
        Condition("A", tmp_path / "a.yaml", (1, 2), cast(Any, SimpleNamespace())),
        Condition("B", tmp_path / "b.yaml", (3,), cast(Any, SimpleNamespace())),
    ]

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.prepare_comparison_run",
        lambda analysis, config, equilibration: {
            "all_conditions": valid_conditions,
            "valid_conditions": valid_conditions,
            "excluded_conditions": [],
            "condition_by_label": {condition.label: condition for condition in valid_conditions},
            "settings": _Settings(),
            "equilibration": "5ns",
            "analysis_root": tmp_path / "analysis",
        },
    )
    config = cast(Any, SimpleNamespace(source_path=tmp_path / "comparison.yaml"))

    manifest = build_manifest(
        analysis,
        config,
        resources,
        recompute=True,
        equilibration=None,
        allow_partial=True,
    )
    assert manifest.partial_policy == "allow_partial"


def test_build_manifest_sets_finalize_only_pipeline_mode(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Comparator-only plugins should emit finalize-only manifests."""
    analysis = _CompareOnlyAnalysis()
    resources = AnalysisSlurmResources(partition="gpu")
    valid_conditions = [
        Condition("A", tmp_path / "a.yaml", (1, 2), cast(Any, SimpleNamespace())),
    ]

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.prepare_comparison_run",
        lambda analysis, config, equilibration: {
            "all_conditions": valid_conditions,
            "valid_conditions": valid_conditions,
            "excluded_conditions": [],
            "condition_by_label": {condition.label: condition for condition in valid_conditions},
            "settings": _Settings(),
            "equilibration": "5ns",
            "analysis_root": tmp_path / "analysis",
        },
    )
    config = cast(Any, SimpleNamespace(source_path=tmp_path / "comparison.yaml"))

    manifest = build_manifest(analysis, config, resources, recompute=True, equilibration=None)
    assert manifest.pipeline_mode == "finalize_only"


def test_submit_analysis_graph_builds_dependencies(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Submission should wire replicate->aggregate->finalize dependencies."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    calls: list[str | None] = []

    def _fake_submit(script_path: Path, dependency: str | None = None) -> str:
        del script_path
        calls.append(dependency)
        return str(1000 + len(calls))

    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._submit_sbatch", _fake_submit)
    graph = submit_analysis_graph(manifest, resources, hpc_dir)

    assert len(graph.replicate_jobs) == 2
    assert len(graph.aggregator_jobs) == 1
    assert graph.finalizer_job_id.isdigit()
    assert calls[0] is None
    assert calls[1] is None
    assert calls[2] == "afterany:1001:1002"
    assert calls[3] == "afterany:1003"

    loaded_graph = SubmittedJobGraph.load(hpc_dir / "job_graph.json")
    assert loaded_graph.aggregator_jobs[0] == graph.aggregator_jobs[0]


def test_submit_analysis_graph_finalize_only_submits_single_job(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Finalize-only mode should submit only the finalizer job."""
    manifest = _make_manifest(tmp_path)
    manifest.pipeline_mode = "finalize_only"
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    calls: list[str | None] = []

    def _fake_submit(script_path: Path, dependency: str | None = None) -> str:
        del script_path
        calls.append(dependency)
        return "9001"

    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._submit_sbatch", _fake_submit)
    graph = submit_analysis_graph(manifest, resources, hpc_dir)

    assert graph.replicate_jobs == {}
    assert graph.aggregator_jobs == {}
    assert graph.finalizer_job_id == "9001"
    assert calls == [None]


def test_submit_analysis_graph_finalize_only_uses_root_dependencies(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Finalize-only mode should attach optional root dependencies."""
    manifest = _make_manifest(tmp_path)
    manifest.pipeline_mode = "finalize_only"
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    calls: list[str | None] = []

    def _fake_submit(script_path: Path, dependency: str | None = None) -> str:
        del script_path
        calls.append(dependency)
        return "9002"

    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._submit_sbatch", _fake_submit)
    submit_analysis_graph(manifest, resources, hpc_dir, root_dependencies=("101", "102"))

    assert calls == ["afterok:101:102"]


def test_manifest_load_defaults_pipeline_mode_for_legacy_json(tmp_path: Path) -> None:
    """Legacy manifests without pipeline_mode should default to full mode."""
    legacy_payload = {
        "analysis_name": "toy_slurm",
        "comparison_yaml": str(tmp_path / "comparison.yaml"),
        "condition_specs": [
            {
                "condition_index": 0,
                "condition_label": "Cond A",
                "condition_slug": "cond_a",
                "replicate_specs": [
                    {
                        "condition_index": 0,
                        "replicate": 1,
                        "condition_label": "Cond A",
                        "condition_slug": "cond_a",
                    }
                ],
            }
        ],
        "settings_snapshot": {"threshold": 1.0},
        "snapshot_hash": "",
        "partial_policy": "strict",
        "equilibration": "10ns",
        "recompute": False,
        "resources": AnalysisSlurmResources().model_dump(mode="json"),
        "created_at": "2026-01-01T00:00:00+00:00",
    }
    legacy_path = tmp_path / "legacy_manifest.json"
    legacy_path.write_text(json.dumps(legacy_payload))

    loaded = AnalysisJobManifest.load(legacy_path)
    assert loaded.pipeline_mode == "full"


def test_submit_analysis_graph_with_arrays_builds_dependencies(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Array submission should wire array->aggregate->finalize dependencies."""
    manifest = AnalysisJobManifest(
        analysis_name="toy_slurm",
        comparison_yaml=str(tmp_path / "comparison.yaml"),
        condition_specs=[
            ConditionTaskSpec(
                condition_index=0,
                condition_label="Cond A",
                condition_slug="cond_a",
                replicate_specs=[
                    ReplicateTaskSpec(
                        condition_index=0,
                        replicate=1,
                        condition_label="Cond A",
                        condition_slug="cond_a",
                    ),
                    ReplicateTaskSpec(
                        condition_index=0,
                        replicate=2,
                        condition_label="Cond A",
                        condition_slug="cond_a",
                    ),
                ],
            ),
            ConditionTaskSpec(
                condition_index=1,
                condition_label="Cond B",
                condition_slug="cond_b",
                replicate_specs=[
                    ReplicateTaskSpec(
                        condition_index=1,
                        replicate=1,
                        condition_label="Cond B",
                        condition_slug="cond_b",
                    ),
                    ReplicateTaskSpec(
                        condition_index=1,
                        replicate=3,
                        condition_label="Cond B",
                        condition_slug="cond_b",
                    ),
                ],
            ),
        ],
        settings_snapshot={"threshold": 1.0},
        equilibration="10ns",
        recompute=False,
        resources=AnalysisSlurmResources(),
        created_at="2026-01-01T00:00:00+00:00",
    )
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    calls: list[str | None] = []

    def _fake_submit(script_path: Path, dependency: str | None = None) -> str:
        del script_path
        calls.append(dependency)
        return str(2000 + len(calls))

    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._submit_sbatch", _fake_submit)
    graph = submit_analysis_graph_with_arrays(manifest, resources, hpc_dir)

    assert graph.replicate_jobs == {}
    assert graph.array_jobs is not None
    assert len(graph.array_jobs) == 2
    assert len(graph.aggregator_jobs) == 2
    assert graph.finalizer_job_id.isdigit()
    assert len(calls) == 5
    assert calls[0] is None
    assert calls[1] is None
    assert calls[2] == "afterany:2001"
    assert calls[3] == "afterany:2002"
    assert calls[4] == "afterany:2003:2004"


def test_submit_analysis_graph_with_arrays_aggregate_depends_on_array_job_id(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Aggregate dependency should reference parent array job ID."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    calls: list[str | None] = []

    def _fake_submit(script_path: Path, dependency: str | None = None) -> str:
        del script_path
        calls.append(dependency)
        return str(3000 + len(calls))

    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._submit_sbatch", _fake_submit)
    submit_analysis_graph_with_arrays(manifest, resources, hpc_dir)

    assert calls[1] == "afterany:3001"


def test_submit_analysis_graph_with_arrays_rolls_back_on_aggregate_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Array submission rollback should cancel array and aggregate jobs."""
    manifest = _make_manifest_two_conditions(tmp_path)
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    submit_mock = MagicMock(side_effect=["2001", "2002", "2003", RuntimeError("aggregate failure")])
    cancel_mock = MagicMock(
        return_value={
            "2001": {"attempted": True, "cancelled": True, "attempts": 1, "error": None},
            "2002": {"attempted": True, "cancelled": True, "attempts": 1, "error": None},
            "2003": {"attempted": True, "cancelled": True, "attempts": 1, "error": None},
        }
    )

    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._submit_sbatch", submit_mock)
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._cancel_jobs", cancel_mock)

    with pytest.raises(RuntimeError, match="aggregate failure"):
        submit_analysis_graph_with_arrays(manifest, resources, hpc_dir)

    cancel_mock.assert_called_once_with(["2001", "2002", "2003"])


def test_submit_analysis_graph_rolls_back_on_replicate_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Replicate submission failure should cancel already-submitted jobs."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    submit_mock = MagicMock(side_effect=["1001", RuntimeError("replicate failure")])
    cancel_mock = MagicMock(
        return_value={"1001": {"attempted": True, "cancelled": True, "attempts": 1, "error": None}}
    )

    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._submit_sbatch", submit_mock)
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._cancel_jobs", cancel_mock)

    with pytest.raises(RuntimeError, match="replicate failure"):
        submit_analysis_graph(manifest, resources, hpc_dir)

    cancel_mock.assert_called_once_with(["1001"])


def test_submit_analysis_graph_rolls_back_on_aggregate_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Aggregate submission failure should cancel replicate and aggregate jobs."""
    manifest = _make_manifest_two_conditions(tmp_path)
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    submit_mock = MagicMock(side_effect=["1001", "1002", "1003", RuntimeError("aggregate failure")])
    cancel_mock = MagicMock(
        return_value={
            "1001": {"attempted": True, "cancelled": True, "attempts": 1, "error": None},
            "1002": {"attempted": True, "cancelled": True, "attempts": 1, "error": None},
            "1003": {"attempted": True, "cancelled": True, "attempts": 1, "error": None},
        }
    )

    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._submit_sbatch", submit_mock)
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._cancel_jobs", cancel_mock)

    with pytest.raises(RuntimeError, match="aggregate failure"):
        submit_analysis_graph(manifest, resources, hpc_dir)

    cancel_mock.assert_called_once_with(["1001", "1002", "1003"])


def test_submit_analysis_graph_rolls_back_on_finalize_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Finalize submission failure should cancel all upstream jobs."""
    manifest = _make_manifest_two_conditions(tmp_path)
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    submit_mock = MagicMock(
        side_effect=["1001", "1002", "1003", "1004", RuntimeError("finalize failure")]
    )
    cancel_mock = MagicMock(
        return_value={
            "1001": {"attempted": True, "cancelled": True, "attempts": 1, "error": None},
            "1002": {"attempted": True, "cancelled": True, "attempts": 1, "error": None},
            "1003": {"attempted": True, "cancelled": True, "attempts": 1, "error": None},
            "1004": {"attempted": True, "cancelled": True, "attempts": 1, "error": None},
        }
    )

    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._submit_sbatch", submit_mock)
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._cancel_jobs", cancel_mock)

    with pytest.raises(RuntimeError, match="finalize failure"):
        submit_analysis_graph(manifest, resources, hpc_dir)

    cancel_mock.assert_called_once_with(["1001", "1002", "1003", "1004"])


def test_cancel_jobs_success(monkeypatch: pytest.MonkeyPatch) -> None:
    """Cancel helper should submit scancel command and return outcomes."""
    run_mock = MagicMock()
    monkeypatch.setattr("subprocess.run", run_mock)

    result = _cancel_jobs(["1001", "1002"])

    assert result["1001"]["cancelled"] is True
    assert result["1002"]["cancelled"] is True
    assert run_mock.call_count == 2
    run_mock.assert_any_call(
        ["scancel", "1001"],
        capture_output=True,
        text=True,
        check=True,
    )
    run_mock.assert_any_call(
        ["scancel", "1002"],
        capture_output=True,
        text=True,
        check=True,
    )


def test_cancel_jobs_failure_is_suppressed(monkeypatch: pytest.MonkeyPatch) -> None:
    """Cancel helper should not raise when scancel fails."""

    def _raise_called_process_error(*args: Any, **kwargs: Any) -> None:
        del args, kwargs
        raise subprocess.CalledProcessError(returncode=1, cmd=["scancel"], stderr="cannot cancel")

    monkeypatch.setattr("subprocess.run", _raise_called_process_error)

    result = _cancel_jobs(["1001"])
    assert result["1001"]["cancelled"] is False
    assert result["1001"]["attempts"] == 2


def test_submit_analysis_graph_writes_submission_error_sidecar(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Submission failure should persist rollback metadata to sidecar JSON."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    submit_mock = MagicMock(side_effect=["1001", RuntimeError("replicate failure")])
    cancel_mock = MagicMock(
        return_value={"1001": {"attempted": True, "cancelled": True, "attempts": 1, "error": None}}
    )

    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._submit_sbatch", submit_mock)
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._cancel_jobs", cancel_mock)

    with pytest.raises(RuntimeError, match="replicate failure"):
        submit_analysis_graph(manifest, resources, hpc_dir)

    error_path = hpc_dir / "submission_error.json"
    assert error_path.exists()
    payload = json.loads(error_path.read_text())
    assert payload["error"] == "replicate failure"
    assert payload["cancelled_job_ids"] == ["1001"]
    assert payload["cancellation_results"]["1001"]["cancelled"] is True
    assert isinstance(payload["timestamp"], str)


def test_rollback_on_first_submission_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """First submission failure should not attempt cancellation."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm._submit_sbatch",
        MagicMock(side_effect=RuntimeError("first failure")),
    )
    cancel_mock = MagicMock(return_value={})
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._cancel_jobs", cancel_mock)

    with pytest.raises(RuntimeError, match="first failure"):
        submit_analysis_graph(manifest, resources, hpc_dir)

    cancel_mock.assert_not_called()


def test_cancel_jobs_validates_job_id_format(monkeypatch: pytest.MonkeyPatch) -> None:
    """Cancellation should skip malformed job IDs and allow array format IDs."""
    run_mock = MagicMock()
    monkeypatch.setattr("subprocess.run", run_mock)

    result = _cancel_jobs(["1001", "1001_1", "--name=evil", "abc"])

    assert result["1001"]["cancelled"] is True
    assert result["1001_1"]["cancelled"] is True
    assert result["--name=evil"]["attempted"] is False
    assert result["abc"]["attempted"] is False
    assert run_mock.call_count == 2


def test_stale_submission_error_cleaned_on_new_submission(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Old submission_error sidecar should be removed before a new submission."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    hpc_dir.mkdir(parents=True, exist_ok=True)
    stale_sidecar = hpc_dir / "submission_error.json"
    stale_sidecar.write_text(json.dumps({"error": "old"}))

    def _fake_submit(script_path: Path, dependency: str | None = None) -> str:
        del script_path, dependency
        return "1001"

    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._submit_sbatch", _fake_submit)
    graph = submit_analysis_graph(manifest, resources, hpc_dir)

    assert graph.finalizer_job_id == "1001"
    assert not stale_sidecar.exists()


def test_sidecar_write_failure_preserves_original_exception(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Sidecar write errors should not mask submission failures."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm._submit_sbatch",
        MagicMock(side_effect=RuntimeError("replicate failure")),
    )

    original_write_text = Path.write_text

    def _selective_write_text(path: Path, data: str, *args: Any, **kwargs: Any) -> int:
        if path.name == "submission_error.json":
            raise OSError("disk full")
        return original_write_text(path, data, *args, **kwargs)

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.Path.write_text",
        _selective_write_text,
    )

    with pytest.raises(RuntimeError, match="replicate failure"):
        submit_analysis_graph(manifest, resources, hpc_dir)


def test_cancel_jobs_retries_on_transient_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    """Cancellation should retry once after a transient scancel failure."""
    run_mock = MagicMock(
        side_effect=[
            subprocess.CalledProcessError(returncode=1, cmd=["scancel"], stderr="not found"),
            subprocess.CompletedProcess(args=["scancel"], returncode=0, stdout="", stderr=""),
        ]
    )
    sleep_mock = MagicMock()
    monkeypatch.setattr("subprocess.run", run_mock)
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm.time.sleep", sleep_mock)

    result = _cancel_jobs(["1001"])

    assert result["1001"]["cancelled"] is True
    assert result["1001"]["attempts"] == 2
    sleep_mock.assert_called_once_with(2)


def test_submit_analysis_graph_parse_failure_logs_raw_stdout_in_sidecar(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Parse failures should persist raw sbatch stdout for manual recovery."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    parse_error = RuntimeError(
        "SBATCH_PARSE_FAILURE: Could not parse job id from sbatch output. "
        "Raw stdout: 'Submitted batch with unknown format'"
    )

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm._submit_sbatch",
        MagicMock(side_effect=parse_error),
    )

    with pytest.raises(RuntimeError, match="SBATCH_PARSE_FAILURE"):
        submit_analysis_graph(manifest, resources, hpc_dir)

    payload = json.loads((hpc_dir / "submission_error.json").read_text())
    assert payload["raw_sbatch_stdout"] == "'Submitted batch with unknown format'"


def test_status_update_and_read_summary(tmp_path: Path) -> None:
    """Status utilities should write atomically and summarize correctly."""
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    rep_status = hpc_dir / "status" / "replicates" / "cond_a" / "rep_1.json"
    cond_status = hpc_dir / "status" / "conditions" / "cond_a.json"
    fin_status = hpc_dir / "status" / "finalize.json"

    update_task_status(rep_status, "running", 1)
    update_task_status(cond_status, "pending", 0)
    update_task_status(fin_status, "failed", 2, "boom")

    summary = read_analysis_status(hpc_dir)
    assert summary["replicates"]["cond_a"]["rep_1"]["state"] == "running"
    assert summary["conditions"]["cond_a"]["state"] == "pending"
    assert summary["finalize"]["state"] == "failed"
    assert summary["counts"]["failed"] == 1


def test_resource_validation_rejects_invalid_cpu_count() -> None:
    """Resources model should validate CPU count bounds."""
    with pytest.raises(ValidationError):
        AnalysisSlurmResources(cpus_per_task=0)


def test_sanitize_slurm_rejects_newlines() -> None:
    """SLURM sanitizer should reject newline injection attempts."""
    with pytest.raises(ValueError, match="Unsafe SLURM value"):
        _sanitize_slurm_value("normal\n#SBATCH --wrap=bad", "partition")


def test_sanitize_slurm_rejects_shell_metacharacters() -> None:
    """SLURM sanitizer should reject shell metacharacters."""
    for bad_value in (
        "`whoami`",
        "$(whoami)",
        "aa100|cat",
        "normal; rm -rf /",
        "${HOME}",
        "normal < /etc/passwd",
        "normal > /tmp/evil",
    ):
        with pytest.raises(ValueError, match="Unsafe SLURM value"):
            _sanitize_slurm_value(bad_value, "partition")


def test_sanitize_path_for_script_accepts_normal_path() -> None:
    """Path sanitizer should accept standard absolute paths."""
    path = Path("/home/user/project/analysis")
    assert _sanitize_path_for_script(path) == str(path.resolve())


def test_sanitize_path_for_script_accepts_spaces() -> None:
    """Path sanitizer should allow spaces in paths."""
    path = Path("/home/user/my project/analysis")
    assert _sanitize_path_for_script(path) == str(path.resolve())


def test_sanitize_path_for_script_accepts_common_symbols() -> None:
    """Path sanitizer should allow hyphens, underscores, and dots."""
    path = Path("/home/user/my-project_name.v1/analysis")
    assert _sanitize_path_for_script(path) == str(path.resolve())


@pytest.mark.parametrize(
    "bad_path",
    [
        Path("/home/user/bad\tpath/analysis"),
        Path("/home/user/bad\vpath/analysis"),
        Path("/home/user/bad\fpath/analysis"),
    ],
)
def test_sanitize_path_for_script_rejects_control_characters(bad_path: Path) -> None:
    """Path sanitizer should reject ASCII control characters."""
    with pytest.raises(ValueError, match="Unsafe path for generated scripts"):
        _sanitize_path_for_script(bad_path)


def test_sanitize_path_for_script_rejects_null_byte() -> None:
    """Path sanitizer should reject null bytes."""
    bad_path = Path("/home/user/bad\x00path/analysis")

    with pytest.raises(ValueError, match="Unsafe path for generated scripts"):
        _sanitize_path_for_script(bad_path)


@pytest.mark.parametrize(
    "bad_path",
    [
        Path("/home/user/bad'path/analysis"),
        Path('/home/user/bad"path/analysis'),
        Path("/home/user/bad\npath/analysis"),
        Path("/home/user/bad`path/analysis"),
        Path("/home/user/bad$(path)/analysis"),
        Path("/home/user/bad|path/analysis"),
    ],
)
def test_sanitize_path_for_script_rejects_unsafe_tokens(bad_path: Path) -> None:
    """Path sanitizer should reject shell-unsafe and quote-breaking tokens."""
    with pytest.raises(ValueError, match="Unsafe path for generated scripts"):
        _sanitize_path_for_script(bad_path)


def test_script_generation_rejects_unsafe_hpc_path(tmp_path: Path) -> None:
    """Script generation should fail fast when interpolated paths are unsafe."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources()
    bad_hpc_dir = Path("/tmp/polyzymd_bad'path")

    with pytest.raises(ValueError, match="Unsafe path for generated scripts"):
        generate_replicate_script(
            manifest,
            manifest.condition_specs[0].replicate_specs[0],
            resources,
            bad_hpc_dir,
        )


def test_script_generation_quotes_sbatch_output_path_with_spaces(tmp_path: Path) -> None:
    """Generated scripts should quote SBATCH output paths containing spaces."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources()
    hpc_dir = tmp_path / "comparison root" / "toy slurm" / "_hpc"

    script = generate_finalize_script(manifest, resources, hpc_dir)
    text = script.read_text()
    expected_log_path = (hpc_dir / "logs" / "finalize.%j.out").resolve()

    assert f'#SBATCH --output="{expected_log_path}"' in text


def test_slurm_resources_validates_mem_format() -> None:
    """Resources model should reject invalid memory strings."""
    with pytest.raises(ValidationError, match="Invalid mem format"):
        AnalysisSlurmResources(mem="four-gigabytes")


def test_slurm_resources_validates_time_format() -> None:
    """Resources model should reject invalid time strings."""
    with pytest.raises(ValidationError, match="Invalid time format"):
        AnalysisSlurmResources(time="99:99")


def test_slurm_resources_validates_task_bounds() -> None:
    """Resources model should enforce ntasks and cpu bounds."""
    with pytest.raises(ValidationError):
        AnalysisSlurmResources(ntasks=0)
    with pytest.raises(ValidationError):
        AnalysisSlurmResources(cpus_per_task=257)


def test_submission_without_sbatch(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Missing sbatch should raise a clear runtime error."""

    def _raise_not_found(*args, **kwargs):
        del args, kwargs
        raise FileNotFoundError("sbatch")

    monkeypatch.setattr("subprocess.run", _raise_not_found)
    with pytest.raises(RuntimeError, match="sbatch' not found on PATH"):
        _submit_sbatch(tmp_path / "script.sh")


def test_corrupted_status_json_reported(tmp_path: Path) -> None:
    """Corrupted status files should be surfaced in summary warnings."""
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    rep_status = hpc_dir / "status" / "replicates" / "cond_a" / "rep_1.json"
    cond_status = hpc_dir / "status" / "conditions" / "cond_a.json"
    fin_status = hpc_dir / "status" / "finalize.json"
    rep_status.parent.mkdir(parents=True, exist_ok=True)
    cond_status.parent.mkdir(parents=True, exist_ok=True)
    fin_status.parent.mkdir(parents=True, exist_ok=True)
    rep_status.write_text("{bad json")
    cond_status.write_text("{bad json")
    fin_status.write_text("{bad json")

    summary = read_analysis_status(hpc_dir)
    assert summary["replicates"]["cond_a"]["rep_1"]["state"] == "unknown"
    assert summary["conditions"]["cond_a"]["state"] == "unknown"
    assert summary["finalize"]["state"] == "unknown"
    assert len(summary["warnings"]) == 3


def test_pixi_path_configurable(tmp_path: Path) -> None:
    """Generated scripts should use the configured pixi executable path."""
    manifest = _make_manifest(tmp_path)
    resources = AnalysisSlurmResources(pixi_path="/opt/pixi/bin/pixi")
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    script = generate_finalize_script(manifest, resources, hpc_dir)
    assert '"/opt/pixi/bin/pixi" run -e build' in script.read_text()


def test_query_sacct_parses_output(monkeypatch: pytest.MonkeyPatch) -> None:
    """sacct parser should map job IDs to state strings."""

    def _fake_run(*args: Any, **kwargs: Any) -> subprocess.CompletedProcess[str]:
        del args, kwargs
        return subprocess.CompletedProcess(
            args=["sacct"],
            returncode=0,
            stdout="1001|FAILED\n1002|COMPLETED\n",
            stderr="",
        )

    monkeypatch.setattr("subprocess.run", _fake_run)
    parsed = _query_sacct(["1001", "1002"])
    assert parsed == {"1001": "FAILED", "1002": "COMPLETED"}


def test_query_sacct_duplicate_rows_prefer_terminal_state(monkeypatch: pytest.MonkeyPatch) -> None:
    """sacct parser should prefer terminal states for duplicate job rows."""

    def _fake_run(*args: Any, **kwargs: Any) -> subprocess.CompletedProcess[str]:
        del args, kwargs
        return subprocess.CompletedProcess(
            args=["sacct"],
            returncode=0,
            stdout=(
                "1001|RUNNING\n"
                "1001|PENDING\n"
                "1001|FAILED\n"
                "1001|COMPLETED\n"
                "1002|RUNNING\n"
                "1002|PENDING\n"
            ),
            stderr="",
        )

    monkeypatch.setattr("subprocess.run", _fake_run)
    parsed = _query_sacct(["1001", "1002"])
    assert parsed == {"1001": "COMPLETED", "1002": "PENDING"}


def test_query_sacct_empty_output_returns_empty(monkeypatch: pytest.MonkeyPatch) -> None:
    """Empty sacct output should return an empty state map."""

    def _fake_run(*args: Any, **kwargs: Any) -> subprocess.CompletedProcess[str]:
        del args, kwargs
        return subprocess.CompletedProcess(args=["sacct"], returncode=0, stdout="", stderr="")

    monkeypatch.setattr("subprocess.run", _fake_run)
    assert _query_sacct(["1001"]) == {}


def test_query_sacct_malformed_output_rows_ignored(monkeypatch: pytest.MonkeyPatch) -> None:
    """Malformed sacct output rows should be ignored safely."""

    def _fake_run(*args: Any, **kwargs: Any) -> subprocess.CompletedProcess[str]:
        del args, kwargs
        return subprocess.CompletedProcess(
            args=["sacct"],
            returncode=0,
            stdout=(
                "1001 FAILED\n"  # missing delimiter
                "1002|RUNNING|EXTRA\n"  # extra columns
                "1003|FAILED\n"  # valid
            ),
            stderr="",
        )

    monkeypatch.setattr("subprocess.run", _fake_run)
    assert _query_sacct(["1001", "1002", "1003"]) == {"1003": "FAILED"}


def test_query_sacct_missing_command_returns_empty(monkeypatch: pytest.MonkeyPatch) -> None:
    """Missing sacct command should return an empty state map."""

    def _raise_not_found(*args: Any, **kwargs: Any) -> None:
        del args, kwargs
        raise FileNotFoundError("sacct")

    monkeypatch.setattr("subprocess.run", _raise_not_found)
    assert _query_sacct(["1001"]) == {}


def test_query_sacct_error_returns_empty(monkeypatch: pytest.MonkeyPatch) -> None:
    """sacct non-zero exit should return an empty state map."""

    def _raise_called_process_error(*args: Any, **kwargs: Any) -> None:
        del args, kwargs
        raise subprocess.CalledProcessError(returncode=1, cmd=["sacct"], stderr="accounting down")

    monkeypatch.setattr("subprocess.run", _raise_called_process_error)
    assert _query_sacct(["1001"]) == {}


def test_query_sacct_filters_substep_ids(monkeypatch: pytest.MonkeyPatch) -> None:
    """sacct parser should ignore sub-step job IDs like .batch."""

    def _fake_run(*args: Any, **kwargs: Any) -> subprocess.CompletedProcess[str]:
        del args, kwargs
        return subprocess.CompletedProcess(
            args=["sacct"],
            returncode=0,
            stdout="1001|FAILED\n1001.batch|COMPLETED\n1002.extern|FAILED\n",
            stderr="",
        )

    monkeypatch.setattr("subprocess.run", _fake_run)
    assert _query_sacct(["1001", "1002"]) == {"1001": "FAILED"}


@pytest.mark.parametrize(
    ("slurm_state", "expected"),
    [
        ("COMPLETED", "succeeded"),
        ("FAILED", "failed"),
        ("OUT_OF_MEMORY", "failed"),
        ("TIMEOUT", "failed"),
        ("NODE_FAIL", "failed"),
        ("BOOT_FAIL", "failed"),
        ("DEADLINE", "failed"),
        ("REVOKED", "failed"),
        ("CANCELLED", "failed"),
        ("CANCELLED+", "failed"),
        ("PREEMPTED", "failed"),
        ("RUNNING", None),
        ("PENDING", None),
        ("SUSPENDED", None),
        ("REQUEUED", None),
    ],
)
def test_map_slurm_state_documented_cases(slurm_state: str, expected: str | None) -> None:
    """SLURM state mapping should follow reconciliation rules."""
    assert _map_slurm_state(slurm_state) == expected


def test_reconcile_status_updates_running_failed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Reconciliation should update running tasks when sacct reports failure."""
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    rep_status = hpc_dir / "status" / "replicates" / "cond_a" / "rep_1.json"
    update_task_status(rep_status, "running", 1)

    payload = json.loads(rep_status.read_text())
    payload["slurm_job_id"] = "1001"
    rep_status.write_text(json.dumps(payload, indent=2))

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm._query_sacct",
        lambda job_ids: {"1001": "FAILED"},
    )
    summary = reconcile_status_with_slurm(hpc_dir)

    updated_payload = json.loads(rep_status.read_text())
    assert summary["checked"] == 1
    assert summary["updated"] == 1
    assert updated_payload["state"] == "failed"
    assert updated_payload["reconciled_from"] == "FAILED"
    assert "reconciled_at" in updated_payload


def test_reconcile_retrying_task_is_actionable(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Reconciliation should include retrying tasks with SLURM job IDs."""
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    rep_status = hpc_dir / "status" / "replicates" / "cond_a" / "rep_1.json"
    update_task_status(rep_status, "retrying", 1)

    payload = json.loads(rep_status.read_text())
    payload["slurm_job_id"] = "3001"
    rep_status.write_text(json.dumps(payload, indent=2))

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm._query_sacct",
        lambda job_ids: {"3001": "FAILED"},
    )
    summary = reconcile_status_with_slurm(hpc_dir)

    updated_payload = json.loads(rep_status.read_text())
    assert summary["checked"] == 1
    assert summary["updated"] == 1
    assert updated_payload["state"] == "failed"
    assert any(
        change["job_id"] == "3001"
        and change["old_state"] == "retrying"
        and change["new_state"] == "failed"
        for change in summary["changes"]
    )


def test_reconcile_status_skips_completed_files(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Reconciliation should not query or update terminal completed status files."""
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    rep_status = hpc_dir / "status" / "replicates" / "cond_a" / "rep_1.json"
    rep_status.parent.mkdir(parents=True, exist_ok=True)
    rep_status.write_text(
        json.dumps(
            {
                "state": "completed",
                "attempt_count": 1,
                "error_message": None,
                "last_updated": "2026-01-01T00:00:00+00:00",
                "slurm_job_id": "1001",
            },
            indent=2,
        )
    )

    query_mock = MagicMock(return_value={"1001": "FAILED"})
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm._query_sacct", query_mock)
    summary = reconcile_status_with_slurm(hpc_dir)

    query_mock.assert_not_called()
    assert summary == {"checked": 0, "updated": 0, "changes": []}
    assert json.loads(rep_status.read_text())["state"] == "completed"


def test_reconcile_status_empty_when_no_actionable_jobs(tmp_path: Path) -> None:
    """Reconciliation should return empty summary when no pending/running jobs exist."""
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    rep_status = hpc_dir / "status" / "replicates" / "cond_a" / "rep_1.json"
    update_task_status(rep_status, "succeeded", 1)

    summary = reconcile_status_with_slurm(hpc_dir)
    assert summary == {"checked": 0, "updated": 0, "changes": []}


def test_reconcile_status_updates_finalize_status(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Reconciliation should include and update finalize status file."""
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    finalize_status = hpc_dir / "status" / "finalize.json"
    update_task_status(finalize_status, "running", 1)

    payload = json.loads(finalize_status.read_text())
    payload["slurm_job_id"] = "2001"
    finalize_status.write_text(json.dumps(payload, indent=2))

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm._query_sacct",
        lambda job_ids: {"2001": "COMPLETED"},
    )
    summary = reconcile_status_with_slurm(hpc_dir)

    updated_payload = json.loads(finalize_status.read_text())
    assert summary["checked"] == 1
    assert summary["updated"] == 1
    assert updated_payload["state"] == "succeeded"
    assert updated_payload["reconciled_from"] == "COMPLETED"


def test_reconcile_status_skips_when_state_changes_during_reconcile(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Reconciliation should skip files changed between read and write."""
    hpc_dir = tmp_path / "comparison" / "toy_slurm" / "_hpc"
    rep_status = hpc_dir / "status" / "replicates" / "cond_a" / "rep_1.json"
    update_task_status(rep_status, "running", 1)

    initial_payload = json.loads(rep_status.read_text())
    initial_payload["slurm_job_id"] = "1001"
    rep_status.write_text(json.dumps(initial_payload, indent=2))

    changed_payload = dict(initial_payload)
    changed_payload["state"] = "succeeded"

    original_read_text = Path.read_text
    read_count = {"rep_status": 0}

    def _racey_read_text(path: Path, *args: Any, **kwargs: Any) -> str:
        if path == rep_status:
            read_count["rep_status"] += 1
            if read_count["rep_status"] == 2:
                return json.dumps(changed_payload, indent=2)
        return original_read_text(path, *args, **kwargs)

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.Path.read_text",
        _racey_read_text,
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm._query_sacct",
        lambda job_ids: {"1001": "FAILED"},
    )

    summary = reconcile_status_with_slurm(hpc_dir)
    final_payload = json.loads(rep_status.read_text())

    assert summary["checked"] == 1
    assert summary["updated"] == 0
    assert final_payload["state"] == "running"
    assert "reconciled_from" not in final_payload
