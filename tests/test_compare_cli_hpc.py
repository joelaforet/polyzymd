"""CLI tests for HPC compare commands."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

from click.testing import CliRunner
from pydantic import BaseModel
import pytest

from polyzymd.compare.cli import compare
from polyzymd.workflow.analysis_slurm import (
    ConditionTaskSpec,
    ReplicateTaskSpec,
    compute_manifest_snapshot_hash,
    validate_manifest_snapshot,
)


class _Settings(BaseModel):
    threshold: float = 1.0


def test_submit_dry_run_prints_summary(monkeypatch, tmp_path: Path) -> None:
    """compare submit --dry-run should print submission summary."""
    runner = CliRunner()

    class _FakeAnalysis:
        name = "toy"

    class _Cond:
        def __init__(self, reps):
            self.replicate_specs = [SimpleNamespace(replicate=r) for r in reps]

    manifest = SimpleNamespace(
        condition_specs=[_Cond([1, 2]), _Cond([1])],
        save=lambda path: Path(path).write_text("{}"),
    )

    monkeypatch.setattr(
        "polyzymd.compare.config.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: lambda: _FakeAnalysis()
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.build_manifest", lambda *args, **kwargs: manifest
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.generate_replicate_script",
        lambda *args, **kwargs: tmp_path / "rep.sh",
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.generate_aggregate_script",
        lambda *args, **kwargs: tmp_path / "agg.sh",
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.generate_finalize_script",
        lambda *args, **kwargs: tmp_path / "fin.sh",
    )

    result = runner.invoke(
        compare,
        [
            "submit",
            "toy",
            "--comparison-yaml",
            str(tmp_path / "comparison.yaml"),
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert "Submitted 6 jobs (3 replicate + 2 aggregate + 1 finalize)" in result.output


def test_submit_allow_partial_sets_manifest_policy(monkeypatch, tmp_path: Path) -> None:
    """compare submit should pass allow-partial policy to manifest builder."""
    runner = CliRunner()

    class _FakeAnalysis:
        name = "toy"

    class _Cond:
        def __init__(self, reps):
            self.replicate_specs = [SimpleNamespace(replicate=r) for r in reps]

    manifest = SimpleNamespace(
        condition_specs=[_Cond([1, 2]), _Cond([1])],
        save=lambda path: Path(path).write_text("{}"),
    )
    captured: dict[str, Any] = {}

    monkeypatch.setattr(
        "polyzymd.compare.config.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: lambda: _FakeAnalysis()
    )

    def _capture_build_manifest(*args, **kwargs):
        captured.update(kwargs)
        return manifest

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.build_manifest",
        _capture_build_manifest,
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.generate_replicate_script",
        lambda *args, **kwargs: tmp_path / "rep.sh",
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.generate_aggregate_script",
        lambda *args, **kwargs: tmp_path / "agg.sh",
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.generate_finalize_script",
        lambda *args, **kwargs: tmp_path / "fin.sh",
    )

    result = runner.invoke(
        compare,
        [
            "submit",
            "toy",
            "--comparison-yaml",
            str(tmp_path / "comparison.yaml"),
            "--allow-partial",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert captured["allow_partial"] is True


def test_status_json_output(monkeypatch, tmp_path: Path) -> None:
    """compare status --json should emit JSON summary."""
    runner = CliRunner()

    monkeypatch.setattr(
        "polyzymd.compare.config.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: SimpleNamespace(name="toy")
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.read_analysis_status",
        lambda hpc_dir: {
            "counts": {"pending": 1, "running": 0, "retrying": 0, "succeeded": 2, "failed": 0}
        },
    )

    result = runner.invoke(
        compare,
        ["status", "toy", "--comparison-yaml", str(tmp_path / "comparison.yaml"), "--json"],
    )
    assert result.exit_code == 0
    payload = json.loads(result.output)
    assert payload["counts"]["succeeded"] == 2


def test_status_human_output_includes_unknown_and_warnings(monkeypatch, tmp_path: Path) -> None:
    """compare status should show unknown counts and warnings."""
    runner = CliRunner()

    monkeypatch.setattr(
        "polyzymd.compare.config.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: SimpleNamespace(name="toy")
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.read_analysis_status",
        lambda hpc_dir: {
            "counts": {
                "pending": 1,
                "running": 0,
                "retrying": 0,
                "succeeded": 2,
                "failed": 0,
                "unknown": 2,
            },
            "warnings": ["corrupted status file: /tmp/bad.json"],
        },
    )

    result = runner.invoke(
        compare,
        ["status", "toy", "--comparison-yaml", str(tmp_path / "comparison.yaml")],
    )
    assert result.exit_code == 0
    assert "unknown=2" in result.output
    assert "⚠ Warnings:" in result.output
    assert "- corrupted status file: /tmp/bad.json" in result.output


def test_finalize_command_loads_aggregated_and_runs(monkeypatch, tmp_path: Path) -> None:
    """compare finalize should call finalizer helper and print path."""
    runner = CliRunner()

    class _Analysis:
        name = "toy"

        def _load_aggregated_result(self, path):
            return {"ok": True}

        def figures_output_dir(self, base):
            return Path(base) / "toy"

    condition = SimpleNamespace(label="A")
    config = SimpleNamespace(
        source_path=tmp_path / "comparison.yaml",
        control=None,
        plot_settings=SimpleNamespace(output_dir=tmp_path / "figures"),
        defaults=SimpleNamespace(equilibration_time="10ns"),
        model_copy=lambda deep=True: SimpleNamespace(
            source_path=tmp_path / "comparison.yaml",
            control=None,
            plot_settings=SimpleNamespace(output_dir=tmp_path / "figures"),
            defaults=SimpleNamespace(equilibration_time="10ns"),
        ),
    )

    monkeypatch.setattr("polyzymd.compare.config.ComparisonConfig.from_yaml", lambda path: config)
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: lambda: _Analysis()
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.prepare_comparison_run",
        lambda plugin, config, equilibration: (
            [condition],
            SimpleNamespace(),
            "10ns",
            tmp_path / "analysis",
        ),
    )
    monkeypatch.setattr("polyzymd.compare.io.paths.sanitize_label", lambda label: label)
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator._resolve_settings",
        lambda plugin, config: SimpleNamespace(),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.finalize_comparison_from_disk",
        lambda **kwargs: {"comparison_path": tmp_path / "comparison" / "toy" / "result.json"},
    )

    result = runner.invoke(
        compare,
        ["finalize", "toy", "--comparison-yaml", str(tmp_path / "comparison.yaml")],
    )
    assert result.exit_code == 0
    assert "Saved result:" in result.output


def test_worker_commands_invoke_helpers(monkeypatch, tmp_path: Path) -> None:
    """Worker commands should invoke orchestrator helper functions."""
    runner = CliRunner()

    class _Analysis:
        name = "toy"
        recompute = False
        Settings = _Settings

        def _load_aggregated_result(self, path):
            return {"mean": 1.0}

        def figures_output_dir(self, base):
            return Path(base) / "toy"

    manifest = SimpleNamespace(
        analysis_name="toy",
        comparison_yaml=str(tmp_path / "comparison.yaml"),
        equilibration="10ns",
        recompute=False,
        settings_snapshot={"threshold": 1.0},
        condition_specs=[
            SimpleNamespace(
                condition_index=0,
                condition_label="A",
                replicate_specs=[SimpleNamespace(replicate=1)],
            )
        ],
    )
    condition = SimpleNamespace(label="A", replicates=(1,))
    called = {"rep": 0, "agg": 0, "fin": 0}

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.AnalysisJobManifest.load", lambda path: manifest
    )
    monkeypatch.setattr(
        "polyzymd.compare.config.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(
            source_path=tmp_path / "comparison.yaml",
            control=None,
            plot_settings=SimpleNamespace(output_dir=tmp_path / "figures"),
            defaults=SimpleNamespace(equilibration_time="10ns"),
            model_copy=lambda deep=True: SimpleNamespace(
                source_path=tmp_path / "comparison.yaml",
                control=None,
                plot_settings=SimpleNamespace(output_dir=tmp_path / "figures"),
                defaults=SimpleNamespace(equilibration_time="10ns"),
            ),
        ),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: lambda: _Analysis()
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.validate_manifest_snapshot",
        lambda manifest, plugin, config: ([condition], "10ns", tmp_path / "analysis"),
    )
    monkeypatch.setattr("polyzymd.compare.io.paths.sanitize_label", lambda label: label)
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.run_replicate_once",
        lambda *args, **kwargs: called.__setitem__("rep", called["rep"] + 1) or {"ok": True},
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.aggregate_condition_from_disk",
        lambda *args, **kwargs: called.__setitem__("agg", called["agg"] + 1) or {"ok": True},
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.finalize_comparison_from_disk",
        lambda **kwargs: called.__setitem__("fin", called["fin"] + 1) or {},
    )

    rep_res = runner.invoke(
        compare,
        [
            "worker-replicate",
            "--manifest",
            str(tmp_path / "manifest.json"),
            "--condition-index",
            "0",
            "--replicate",
            "1",
        ],
    )
    agg_res = runner.invoke(
        compare,
        [
            "worker-aggregate",
            "--manifest",
            str(tmp_path / "manifest.json"),
            "--condition-index",
            "0",
        ],
    )
    fin_res = runner.invoke(
        compare,
        ["worker-finalize", "--manifest", str(tmp_path / "manifest.json")],
    )

    assert rep_res.exit_code == 0
    assert agg_res.exit_code == 0
    assert fin_res.exit_code == 0
    assert called == {"rep": 1, "agg": 1, "fin": 1}


def test_finalize_with_missing_conditions(monkeypatch, tmp_path: Path) -> None:
    """Finalize should fail unless allow_partial is requested."""
    runner = CliRunner()

    class _Analysis:
        name = "toy"

        def _load_aggregated_result(self, path):
            del path
            return None

        def figures_output_dir(self, base):
            return Path(base) / "toy"

    condition = SimpleNamespace(label="A")
    config = SimpleNamespace(
        source_path=tmp_path / "comparison.yaml",
        control=None,
        plot_settings=SimpleNamespace(output_dir=tmp_path / "figures"),
        defaults=SimpleNamespace(equilibration_time="10ns"),
        model_copy=lambda deep=True: SimpleNamespace(
            source_path=tmp_path / "comparison.yaml",
            control=None,
            plot_settings=SimpleNamespace(output_dir=tmp_path / "figures"),
            defaults=SimpleNamespace(equilibration_time="10ns"),
        ),
    )

    monkeypatch.setattr("polyzymd.compare.config.ComparisonConfig.from_yaml", lambda path: config)
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: lambda: _Analysis()
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.prepare_comparison_run",
        lambda plugin, config, equilibration: (
            [condition],
            SimpleNamespace(),
            "10ns",
            tmp_path / "analysis",
        ),
    )
    monkeypatch.setattr("polyzymd.compare.io.paths.sanitize_label", lambda label: label)

    blocked = runner.invoke(
        compare,
        ["finalize", "toy", "--comparison-yaml", str(tmp_path / "comparison.yaml")],
    )
    assert blocked.exit_code != 0
    assert "Finalize aborted" in blocked.output

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.finalize_comparison_from_disk",
        lambda **kwargs: {"comparison_path": tmp_path / "comparison" / "toy" / "result.json"},
    )
    allowed = runner.invoke(
        compare,
        [
            "finalize",
            "toy",
            "--comparison-yaml",
            str(tmp_path / "comparison.yaml"),
            "--allow-partial",
        ],
    )
    assert allowed.exit_code == 0
    assert "Warning: missing aggregated results" in allowed.output


def test_manifest_config_drift_detection(monkeypatch, tmp_path: Path) -> None:
    """Workers should fail when live config drifts from manifest snapshot."""
    runner = CliRunner()

    class _Analysis:
        name = "toy"
        Settings = _Settings

    manifest = SimpleNamespace(
        analysis_name="toy",
        comparison_yaml=str(tmp_path / "comparison.yaml"),
        equilibration="10ns",
        recompute=False,
        settings_snapshot={"threshold": 1.0},
        condition_specs=[
            SimpleNamespace(
                condition_index=0,
                condition_label="A",
                replicate_specs=[SimpleNamespace(replicate=1)],
            )
        ],
    )
    condition = SimpleNamespace(label="A", replicates=(1,))

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.AnalysisJobManifest.load", lambda path: manifest
    )
    monkeypatch.setattr(
        "polyzymd.compare.config.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: lambda: _Analysis()
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.validate_manifest_snapshot",
        lambda manifest, plugin, config: (_ for _ in ()).throw(RuntimeError("drift")),
    )
    monkeypatch.setattr("polyzymd.compare.io.paths.sanitize_label", lambda label: label)

    result = runner.invoke(
        compare,
        [
            "worker-replicate",
            "--manifest",
            str(tmp_path / "manifest.json"),
            "--condition-index",
            "0",
            "--replicate",
            "1",
        ],
    )
    assert result.exit_code != 0
    assert result.exception is not None
    assert "drift" in str(result.exception)


def test_validate_manifest_snapshot_detects_real_drift(monkeypatch, tmp_path: Path) -> None:
    """Snapshot validation should fail when settings drift after submission."""

    submitted_conditions = [
        SimpleNamespace(label="A", replicates=(1, 2)),
        SimpleNamespace(label="B", replicates=(1, 2)),
    ]
    submitted_settings = _Settings(threshold=1.0)
    submitted_specs = [
        ConditionTaskSpec(
            condition_index=0,
            condition_label="A",
            condition_slug="A",
            replicate_specs=[
                ReplicateTaskSpec(
                    condition_index=0,
                    replicate=1,
                    condition_label="A",
                    condition_slug="A",
                ),
                ReplicateTaskSpec(
                    condition_index=0,
                    replicate=2,
                    condition_label="A",
                    condition_slug="A",
                ),
            ],
        ),
        ConditionTaskSpec(
            condition_index=1,
            condition_label="B",
            condition_slug="B",
            replicate_specs=[
                ReplicateTaskSpec(
                    condition_index=1,
                    replicate=1,
                    condition_label="B",
                    condition_slug="B",
                ),
                ReplicateTaskSpec(
                    condition_index=1,
                    replicate=2,
                    condition_label="B",
                    condition_slug="B",
                ),
            ],
        ),
    ]
    snapshot_hash = compute_manifest_snapshot_hash(
        analysis_name="toy",
        settings_snapshot=submitted_settings.model_dump(mode="json"),
        condition_specs=submitted_specs,
        equilibration="10ns",
    )
    manifest = SimpleNamespace(
        analysis_name="toy",
        settings_snapshot=submitted_settings.model_dump(mode="json"),
        condition_specs=submitted_specs,
        snapshot_hash=snapshot_hash,
        equilibration="10ns",
    )

    drifted_settings = _Settings(threshold=2.5)
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.prepare_comparison_run",
        lambda analysis, config, equilibration: (
            submitted_conditions,
            drifted_settings,
            "10ns",
            tmp_path / "analysis",
        ),
    )

    with pytest.raises(RuntimeError, match="Manifest/config drift detected"):
        validate_manifest_snapshot(
            cast(Any, manifest),
            cast(Any, SimpleNamespace(name="toy")),
            cast(Any, SimpleNamespace()),
        )


def test_submit_without_sbatch(monkeypatch, tmp_path: Path) -> None:
    """CLI submit should fail cleanly when sbatch is unavailable."""
    runner = CliRunner()

    monkeypatch.setattr("shutil.which", lambda name: None if name == "sbatch" else "/usr/bin/other")
    result = runner.invoke(
        compare,
        ["submit", "toy", "--comparison-yaml", str(tmp_path / "comparison.yaml")],
    )
    assert result.exit_code != 0
    assert "sbatch' not found on PATH" in result.output
