"""Tests for analysis SLURM DAG helpers."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar, cast

import pytest
from pydantic import BaseModel, ValidationError

from polyzymd.analyses.base import Analysis, Condition
from polyzymd.workflow.analysis_slurm import (
    AnalysisJobManifest,
    AnalysisSlurmResources,
    ConditionTaskSpec,
    ReplicateTaskSpec,
    SubmittedJobGraph,
    build_manifest,
    generate_aggregate_script,
    generate_finalize_script,
    generate_replicate_script,
    read_analysis_status,
    submit_analysis_graph,
    update_task_status,
    _sanitize_slurm_value,
    _submit_sbatch,
)


class _Settings(BaseModel):
    threshold: float = 1.0


class _ToyAnalysis(Analysis):
    name: ClassVar[str] = "toy_slurm"
    Settings: ClassVar[type] = _Settings


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
        lambda analysis, config, equilibration: (
            valid_conditions,
            _Settings(),
            "5ns",
            tmp_path / "analysis",
        ),
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
        lambda analysis, config, equilibration: (
            valid_conditions,
            _Settings(),
            "5ns",
            tmp_path / "analysis",
        ),
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
