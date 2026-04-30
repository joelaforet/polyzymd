"""CLI tests for HPC compare commands."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import pytest
from click.testing import CliRunner
from pydantic import BaseModel

from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.cli.compare import compare
from polyzymd.config.comparison import PlotSettings
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
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
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
            "-f",
            str(tmp_path / "comparison.yaml"),
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert "Would submit 6 jobs (3 replicate + 2 aggregate + 1 finalize)" in result.output


def test_submit_dry_run_with_job_arrays_prints_array_summary(monkeypatch, tmp_path: Path) -> None:
    """compare submit --job-arrays --dry-run should print array mode summary."""
    runner = CliRunner()

    class _FakeAnalysis:
        name = "toy"

    class _Cond:
        def __init__(self, reps):
            self.replicate_specs = [SimpleNamespace(replicate=r) for r in reps]
            self.condition_slug = "cond"

    manifest = SimpleNamespace(
        condition_specs=[_Cond([1, 2]), _Cond([1])],
        save=lambda path: Path(path).write_text("{}"),
    )

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: lambda: _FakeAnalysis()
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.build_manifest", lambda *args, **kwargs: manifest
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.generate_array_script",
        lambda *args, **kwargs: tmp_path / "arr.sh",
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
            "-f",
            str(tmp_path / "comparison.yaml"),
            "--job-arrays",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert "Would submit 2 array jobs + 2 aggregate + 1 finalize = 5 total" in result.output
    assert "Mode: job arrays" in result.output


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
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
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
            "-f",
            str(tmp_path / "comparison.yaml"),
            "--allow-partial",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert captured["allow_partial"] is True


def test_submit_dry_run_finalize_only_prints_summary(monkeypatch, tmp_path: Path) -> None:
    """compare submit --dry-run should describe finalize-only submissions."""
    runner = CliRunner()

    class _FakeAnalysis:
        name = "toy"

    class _Cond:
        def __init__(self, reps):
            self.replicate_specs = [SimpleNamespace(replicate=r) for r in reps]

    manifest = SimpleNamespace(
        condition_specs=[_Cond([1, 2])],
        pipeline_mode="finalize_only",
        save=lambda path: Path(path).write_text("{}"),
    )

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: lambda: _FakeAnalysis()
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.build_manifest", lambda *args, **kwargs: manifest
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
            "-f",
            str(tmp_path / "comparison.yaml"),
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert "Would submit 1 finalize job (compare-only plugin)" in result.output


def test_status_json_output(monkeypatch, tmp_path: Path) -> None:
    """compare status --json should emit JSON summary."""
    runner = CliRunner()

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
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
        ["status", "toy", "-f", str(tmp_path / "comparison.yaml"), "--json"],
    )
    assert result.exit_code == 0
    payload = json.loads(result.output)
    assert payload["counts"]["succeeded"] == 2


def test_status_human_output_includes_unknown_and_warnings(monkeypatch, tmp_path: Path) -> None:
    """compare status should show unknown counts and warnings."""
    runner = CliRunner()

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
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
        ["status", "toy", "-f", str(tmp_path / "comparison.yaml")],
    )
    assert result.exit_code == 0
    assert "unknown=2" in result.output
    assert "⚠ Warnings:" in result.output
    assert "- corrupted status file: /tmp/bad.json" in result.output


def test_status_reconcile_prints_summary(monkeypatch, tmp_path: Path) -> None:
    """compare status --reconcile should print reconciliation summary."""
    runner = CliRunner()

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: SimpleNamespace(name="toy")
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.reconcile_status_with_slurm",
        lambda hpc_dir: {
            "checked": 3,
            "updated": 3,
            "changes": [
                {
                    "job_id": "1001",
                    "path": "/tmp/rep_1.json",
                    "old_state": "running",
                    "new_state": "failed",
                    "slurm_state": "OUT_OF_MEMORY",
                },
                {
                    "job_id": "1002",
                    "path": "/tmp/rep_2.json",
                    "old_state": "running",
                    "new_state": "failed",
                    "slurm_state": "OUT_OF_MEMORY",
                },
                {
                    "job_id": "1003",
                    "path": "/tmp/cond_a.json",
                    "old_state": "pending",
                    "new_state": "failed",
                    "slurm_state": "CANCELLED",
                },
            ],
        },
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.read_analysis_status",
        lambda hpc_dir: {
            "counts": {
                "pending": 0,
                "running": 0,
                "retrying": 0,
                "succeeded": 0,
                "failed": 3,
                "unknown": 0,
            },
            "warnings": [],
        },
    )

    result = runner.invoke(
        compare,
        ["status", "toy", "-f", str(tmp_path / "comparison.yaml"), "--reconcile"],
    )
    assert result.exit_code == 0
    assert "Reconciled 3 jobs:" in result.output
    assert "2 marked failed (SLURM: OUT_OF_MEMORY)" in result.output
    assert "1 marked failed (SLURM: CANCELLED)" in result.output


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
        plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
        defaults=SimpleNamespace(equilibration_time="10ns"),
        model_copy=lambda deep=True: SimpleNamespace(
            source_path=tmp_path / "comparison.yaml",
            control=None,
            plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
            defaults=SimpleNamespace(equilibration_time="10ns"),
        ),
    )

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml", lambda path: config
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: lambda: _Analysis()
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.prepare_comparison_run",
        lambda plugin, config, equilibration: {
            "all_conditions": [condition],
            "valid_conditions": [condition],
            "excluded_conditions": [],
            "condition_by_label": {condition.label: condition},
            "settings": SimpleNamespace(),
            "equilibration": "10ns",
            "analysis_root": tmp_path / "analysis",
        },
    )
    monkeypatch.setattr("polyzymd.analyses.shared.paths.sanitize_label", lambda label: label)
    captured: dict[str, Any] = {}
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator._resolve_settings",
        lambda plugin, config: SimpleNamespace(),
    )

    def _capture_finalize(**kwargs):
        captured.update(kwargs)
        return {"comparison_path": tmp_path / "comparison" / "toy" / "result.json"}

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.finalize_comparison_from_disk",
        _capture_finalize,
    )

    result = runner.invoke(
        compare,
        ["finalize", "toy", "-f", str(tmp_path / "comparison.yaml"), "--recompute"],
    )
    assert result.exit_code == 0
    assert "Saved result:" in result.output
    assert captured["recompute"] is True


def test_run_reports_plugin_contract_error(monkeypatch, tmp_path: Path) -> None:
    """compare run should label plugin contract errors as plugin bugs."""
    runner = CliRunner()

    class _Analysis:
        name = "toy"

        def format(self, result, output_format="text"):
            return "formatted"

    config = SimpleNamespace(
        name="contract_project",
        conditions=[SimpleNamespace(label="A")],
        defaults=SimpleNamespace(equilibration_time="10ns"),
    )

    monkeypatch.setattr("polyzymd.cli.compare.load_comparison_config", lambda path: config)
    monkeypatch.setattr("polyzymd.cli.compare.validate_and_report", lambda config: None)
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis",
        lambda name: _Analysis,
    )

    def _raise_contract(*args, **kwargs):
        raise PluginContractError("toy.run_replicate() returned None")

    monkeypatch.setattr("polyzymd.analyses.orchestrator.run_comparison", _raise_contract)

    result = runner.invoke(compare, ["run", "toy", "-f", str(tmp_path / "comparison.yaml")])

    assert result.exit_code != 0
    assert "Plugin contract error in analysis 'toy'" in result.output
    assert "likely a PolyzyMD/plugin bug, not missing trajectory data" in result.output


def test_finalize_compare_only_plugin_skips_aggregated_precheck(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Public finalize should not require aggregates for compare-only plugins."""
    runner = CliRunner()

    class _CompareOnlyAnalysis:
        name = "exposure"
        has_compute_stage = False
        has_aggregate_stage = False

        def _load_aggregated_result(self, path):
            del path
            raise AssertionError("Compare-only finalize should not load aggregated results")

        def figures_output_dir(self, base):
            return Path(base) / "exposure"

    condition_a = SimpleNamespace(label="A")
    condition_b = SimpleNamespace(label="B")
    config = SimpleNamespace(
        source_path=tmp_path / "comparison.yaml",
        control=None,
        plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
        defaults=SimpleNamespace(equilibration_time="10ns"),
        model_copy=lambda deep=True: SimpleNamespace(
            source_path=tmp_path / "comparison.yaml",
            control=None,
            plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
            defaults=SimpleNamespace(equilibration_time="10ns"),
        ),
    )

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: config,
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis",
        lambda name: lambda: _CompareOnlyAnalysis(),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.prepare_comparison_run",
        lambda plugin, config, equilibration: {
            "all_conditions": [condition_a, condition_b],
            "valid_conditions": [condition_a, condition_b],
            "excluded_conditions": [],
            "condition_by_label": {"A": condition_a, "B": condition_b},
            "settings": SimpleNamespace(),
            "equilibration": "10ns",
            "analysis_root": tmp_path / "analysis",
        },
    )
    monkeypatch.setattr("polyzymd.analyses.shared.paths.sanitize_label", lambda label: label)
    captured: dict[str, Any] = {}

    def _capture_finalize(**kwargs):
        captured.update(kwargs)
        return {"comparison_path": tmp_path / "comparison" / "exposure" / "result.json"}

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.finalize_comparison_from_disk",
        _capture_finalize,
    )

    result = runner.invoke(
        compare,
        ["finalize", "exposure", "-f", str(tmp_path / "comparison.yaml")],
    )

    assert result.exit_code == 0
    assert captured["aggregated_results"] == {}
    assert captured["analysis_dirs"] == {
        "A": tmp_path / "analysis" / "A" / "exposure",
        "B": tmp_path / "analysis" / "B" / "exposure",
    }
    assert "missing aggregated results" not in result.output


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
        recompute=True,
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
    captured_recompute: dict[str, bool] = {}

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.AnalysisJobManifest.load", lambda path: manifest
    )
    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(
            source_path=tmp_path / "comparison.yaml",
            control=None,
            plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
            defaults=SimpleNamespace(equilibration_time="10ns"),
            model_copy=lambda deep=True: SimpleNamespace(
                source_path=tmp_path / "comparison.yaml",
                control=None,
                plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
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
    monkeypatch.setattr("polyzymd.analyses.shared.paths.sanitize_label", lambda label: label)
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.run_replicate_once",
        lambda *args, **kwargs: called.__setitem__("rep", called["rep"] + 1)
        or captured_recompute.__setitem__("rep", args[-1])
        or {"ok": True},
    )

    def _capture_aggregate(*args, **kwargs):
        called["agg"] += 1
        captured_recompute["agg"] = kwargs["recompute"]
        return {"ok": True}

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.aggregate_condition_from_disk",
        _capture_aggregate,
    )

    def _capture_finalize(**kwargs):
        called["fin"] += 1
        captured_recompute["fin"] = kwargs["recompute"]
        return {}

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.finalize_comparison_from_disk",
        _capture_finalize,
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
    assert captured_recompute == {"rep": True, "agg": True, "fin": True}


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
        plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
        defaults=SimpleNamespace(equilibration_time="10ns"),
        model_copy=lambda deep=True: SimpleNamespace(
            source_path=tmp_path / "comparison.yaml",
            control=None,
            plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
            defaults=SimpleNamespace(equilibration_time="10ns"),
        ),
    )

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml", lambda path: config
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: lambda: _Analysis()
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.prepare_comparison_run",
        lambda plugin, config, equilibration: {
            "all_conditions": [condition],
            "valid_conditions": [condition],
            "excluded_conditions": [],
            "condition_by_label": {condition.label: condition},
            "settings": SimpleNamespace(),
            "equilibration": "10ns",
            "analysis_root": tmp_path / "analysis",
        },
    )
    monkeypatch.setattr("polyzymd.analyses.shared.paths.sanitize_label", lambda label: label)

    blocked = runner.invoke(
        compare,
        ["finalize", "toy", "-f", str(tmp_path / "comparison.yaml")],
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
            "-f",
            str(tmp_path / "comparison.yaml"),
            "--allow-partial",
        ],
    )
    assert allowed.exit_code == 0
    assert "Warning: missing aggregated results" in allowed.output


def test_worker_finalize_finalize_only_skips_aggregated_loading(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """worker-finalize should not require aggregated condition results in finalize-only mode."""
    runner = CliRunner()

    class _CompareOnlyAnalysis:
        name = "exposure"
        Settings = _Settings
        has_compute_stage = False
        has_aggregate_stage = False

        def _load_aggregated_result(self, path):
            del path
            raise AssertionError(
                "Comparator-only finalize should not load aggregated condition results"
            )

        def figures_output_dir(self, base):
            return Path(base) / "exposure"

    condition = SimpleNamespace(label="Cond A")
    cond_spec = ConditionTaskSpec(
        condition_index=0,
        condition_label="Cond A",
        condition_slug="cond_a",
        replicate_specs=[],
    )
    manifest = SimpleNamespace(
        analysis_name="exposure",
        comparison_yaml=str(tmp_path / "comparison.yaml"),
        settings_snapshot={"threshold": 1.0},
        condition_specs=[cond_spec],
        pipeline_mode="finalize_only",
        partial_policy="strict",
    )
    config = SimpleNamespace(
        source_path=tmp_path / "comparison.yaml",
        control=None,
        plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
        defaults=SimpleNamespace(equilibration_time="10ns"),
        model_copy=lambda deep=True: SimpleNamespace(
            source_path=tmp_path / "comparison.yaml",
            control=None,
            plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
            defaults=SimpleNamespace(equilibration_time="10ns"),
        ),
    )
    captured: dict[str, Any] = {}

    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.AnalysisJobManifest.load", lambda path: manifest
    )
    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml", lambda path: config
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: lambda: _CompareOnlyAnalysis()
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.validate_manifest_snapshot",
        lambda manifest, plugin, config: ([condition], "10ns", tmp_path / "analysis"),
    )
    monkeypatch.setattr("polyzymd.analyses.shared.paths.sanitize_label", lambda label: label)

    def _capture_finalize(**kwargs):
        captured.update(kwargs)
        return {"comparison_path": tmp_path / "comparison" / "exposure" / "result.json"}

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.finalize_comparison_from_disk",
        _capture_finalize,
    )

    result = runner.invoke(
        compare,
        ["worker-finalize", "--manifest", str(tmp_path / "manifest.json")],
    )

    assert result.exit_code == 0
    assert captured["aggregated_results"] == {}
    assert captured["analysis_dirs"] == {
        "Cond A": tmp_path / "analysis" / "Cond A" / "exposure",
    }
    assert captured["recompute"] is False


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
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: lambda: _Analysis()
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.validate_manifest_snapshot",
        lambda manifest, plugin, config: (_ for _ in ()).throw(RuntimeError("drift")),
    )
    monkeypatch.setattr("polyzymd.analyses.shared.paths.sanitize_label", lambda label: label)

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
        lambda analysis, config, equilibration: {
            "all_conditions": submitted_conditions,
            "valid_conditions": submitted_conditions,
            "excluded_conditions": [],
            "condition_by_label": {c.label: c for c in submitted_conditions},
            "settings": drifted_settings,
            "equilibration": "10ns",
            "analysis_root": tmp_path / "analysis",
        },
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
        ["submit", "toy", "-f", str(tmp_path / "comparison.yaml")],
    )
    assert result.exit_code != 0
    assert "sbatch' not found on PATH" in result.output


def test_submit_with_satisfied_dependencies_succeeds(monkeypatch, tmp_path: Path) -> None:
    """compare submit should continue when dependency result files exist."""
    runner = CliRunner()

    class _FakeAnalysis:
        name = "dep_child"

    class _DepChildClass:
        name = "dep_child"
        dependencies = ("contacts",)

        def __call__(self):
            return _FakeAnalysis()

    class _Cond:
        def __init__(self, reps):
            self.replicate_specs = [SimpleNamespace(replicate=r) for r in reps]

    manifest = SimpleNamespace(
        condition_specs=[_Cond([1])],
        save=lambda path: Path(path).write_text("{}"),
    )

    dep_result = tmp_path / "comparison" / "contacts" / "result.json"
    dep_result.parent.mkdir(parents=True, exist_ok=True)
    dep_result.write_text("{}")

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr("polyzymd.analyses.discovery.get_analysis", lambda name: _DepChildClass())
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
        ["submit", "dep_child", "-f", str(tmp_path / "comparison.yaml"), "--dry-run"],
    )
    assert result.exit_code == 0


def test_submit_with_missing_dependencies_fails(monkeypatch, tmp_path: Path) -> None:
    """compare submit should fail with clear error when dependencies are missing."""
    runner = CliRunner()

    class _FakeAnalysis:
        name = "dep_child"

    class _DepChildClass:
        name = "dep_child"
        dependencies = ("contacts",)

        def __call__(self):
            return _FakeAnalysis()

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr("polyzymd.analyses.discovery.get_analysis", lambda name: _DepChildClass())

    result = runner.invoke(
        compare,
        ["submit", "dep_child", "-f", str(tmp_path / "comparison.yaml"), "--dry-run"],
    )
    assert result.exit_code != 0
    assert "depends on 'contacts'" in result.output
    assert "polyzymd compare submit contacts" in result.output


def test_submit_with_no_dependencies_skips_preflight(monkeypatch, tmp_path: Path) -> None:
    """compare submit should skip dependency preflight for independent plugins."""
    runner = CliRunner()

    class _FakeAnalysis:
        name = "independent"

    class _IndependentClass:
        name = "independent"
        dependencies = ()

        def __call__(self):
            return _FakeAnalysis()

    class _Cond:
        def __init__(self, reps):
            self.replicate_specs = [SimpleNamespace(replicate=r) for r in reps]

    manifest = SimpleNamespace(
        condition_specs=[_Cond([1])],
        save=lambda path: Path(path).write_text("{}"),
    )

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis", lambda name: _IndependentClass()
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
        ["submit", "independent", "-f", str(tmp_path / "comparison.yaml"), "--dry-run"],
    )
    assert result.exit_code == 0


def test_submit_all_dry_run_orders_and_submits_enabled_analyses(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """compare submit-all should dry-run all enabled analyses in dependency order."""
    runner = CliRunner()

    class _Plugins:
        def get_enabled_plugins(self):
            return ["exposure", "contacts"]

    config = SimpleNamespace(source_path=tmp_path / "comparison.yaml", plugins=_Plugins())

    class _Contacts:
        name = "contacts"
        dependencies = ()

        def __call__(self):
            return SimpleNamespace(name="contacts")

    class _Exposure:
        name = "exposure"
        dependencies = ("contacts",)

        def __call__(self):
            return SimpleNamespace(name="exposure")

    class _Cond:
        def __init__(self, reps):
            self.replicate_specs = [SimpleNamespace(replicate=r) for r in reps]

    build_calls: list[str] = []

    def _build_manifest(plugin, *args, **kwargs):
        del args, kwargs
        build_calls.append(plugin.name)
        return SimpleNamespace(
            condition_specs=[_Cond([1, 2])],
            pipeline_mode="full",
            save=lambda path: Path(path).write_text("{}"),
        )

    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/sbatch")
    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml", lambda path: config
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.order_analyses_for_execution",
        lambda names, satisfied=None: ["contacts", "exposure"],
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis",
        lambda name: _Contacts() if name == "contacts" else _Exposure(),
    )
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm.build_manifest", _build_manifest)
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
        ["submit-all", "-f", str(tmp_path / "comparison.yaml"), "--dry-run"],
    )

    assert result.exit_code == 0
    assert build_calls == ["contacts", "exposure"]
    assert "Submitted analyses summary:" in result.output
    assert "Dry run only: no jobs were submitted" in result.output
    assert result.output.index("contacts") < result.output.index("exposure")


def test_submit_all_exclude_filters_requested_plugins(monkeypatch, tmp_path: Path) -> None:
    """compare submit-all should skip analyses listed with --exclude."""
    runner = CliRunner()

    class _Plugins:
        def get_enabled_plugins(self):
            return ["contacts", "exposure", "rmsf"]

    config = SimpleNamespace(source_path=tmp_path / "comparison.yaml", plugins=_Plugins())

    captured: dict[str, Any] = {}

    def _order(analysis_names, satisfied=None):
        del satisfied
        captured["names"] = list(analysis_names)
        return ["contacts", "rmsf"]

    class _AnalysisClass:
        def __init__(self, name):
            self.name = name
            self.dependencies = ()

        def __call__(self):
            return SimpleNamespace(name=self.name)

    class _Cond:
        def __init__(self):
            self.replicate_specs = [SimpleNamespace(replicate=1)]

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml", lambda path: config
    )
    monkeypatch.setattr("polyzymd.analyses.orchestrator.order_analyses_for_execution", _order)
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis",
        lambda name: _AnalysisClass(name),
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.build_manifest",
        lambda *args, **kwargs: SimpleNamespace(
            condition_specs=[_Cond()],
            pipeline_mode="full",
            save=lambda path: Path(path).write_text("{}"),
        ),
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
            "submit-all",
            "-f",
            str(tmp_path / "comparison.yaml"),
            "--exclude",
            "exposure",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert captured["names"] == ["contacts", "rmsf"]
    assert "exposure" not in result.output


def test_submit_all_excluded_dependency_with_result_is_satisfied(
    monkeypatch, tmp_path: Path
) -> None:
    """submit-all should satisfy excluded dependencies from existing results."""
    runner = CliRunner()

    class _Plugins:
        def get_enabled_plugins(self):
            return ["contacts", "exposure"]

    config = SimpleNamespace(source_path=tmp_path / "comparison.yaml", plugins=_Plugins())

    contacts_result = tmp_path / "comparison" / "contacts" / "result.json"
    contacts_result.parent.mkdir(parents=True, exist_ok=True)
    contacts_result.write_text("{}")

    captured: dict[str, Any] = {}

    def _order(analysis_names, satisfied=None):
        captured["names"] = list(analysis_names)
        captured["satisfied"] = set() if satisfied is None else set(satisfied)
        return ["exposure"]

    class _AnalysisClass:
        def __init__(self, name):
            self.name = name
            self.dependencies = ()

        def __call__(self):
            return SimpleNamespace(name=self.name)

    class _Cond:
        def __init__(self):
            self.replicate_specs = [SimpleNamespace(replicate=1)]

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml", lambda path: config
    )
    monkeypatch.setattr("polyzymd.analyses.orchestrator.order_analyses_for_execution", _order)
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis",
        lambda name: _AnalysisClass(name),
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.build_manifest",
        lambda *args, **kwargs: SimpleNamespace(
            condition_specs=[_Cond()],
            pipeline_mode="full",
            save=lambda path: Path(path).write_text("{}"),
        ),
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
            "submit-all",
            "-f",
            str(tmp_path / "comparison.yaml"),
            "--exclude",
            "contacts",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert captured["names"] == ["exposure"]
    assert captured["satisfied"] == {"contacts"}


def test_submit_all_finalize_only_plugin_dry_run(monkeypatch, tmp_path: Path) -> None:
    """compare submit-all should account for finalize-only plugin manifests."""
    runner = CliRunner()

    class _Plugins:
        def get_enabled_plugins(self):
            return ["compare_only"]

    config = SimpleNamespace(source_path=tmp_path / "comparison.yaml", plugins=_Plugins())

    class _CompareOnly:
        name = "compare_only"
        dependencies = ()

        def __call__(self):
            return SimpleNamespace(name="compare_only")

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml", lambda path: config
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.order_analyses_for_execution",
        lambda names, satisfied=None: ["compare_only"],
    )
    monkeypatch.setattr("polyzymd.analyses.discovery.get_analysis", lambda name: _CompareOnly())
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.build_manifest",
        lambda *args, **kwargs: SimpleNamespace(
            condition_specs=[],
            pipeline_mode="finalize_only",
            save=lambda path: Path(path).write_text("{}"),
        ),
    )
    monkeypatch.setattr(
        "polyzymd.workflow.analysis_slurm.generate_finalize_script",
        lambda *args, **kwargs: tmp_path / "fin.sh",
    )

    result = runner.invoke(
        compare,
        ["submit-all", "-f", str(tmp_path / "comparison.yaml"), "--dry-run"],
    )

    assert result.exit_code == 0
    assert "compare_only\tfinalize_only\t1\tdry-run:compare_only:finalize" in result.output


def test_submit_all_passes_cross_plugin_root_dependencies(monkeypatch, tmp_path: Path) -> None:
    """compare submit-all should wire finalize job dependencies across plugins."""
    runner = CliRunner()

    class _Plugins:
        def get_enabled_plugins(self):
            return ["contacts", "exposure"]

    config = SimpleNamespace(source_path=tmp_path / "comparison.yaml", plugins=_Plugins())

    class _Contacts:
        name = "contacts"
        dependencies = ()

        def __call__(self):
            return SimpleNamespace(name="contacts")

    class _Exposure:
        name = "exposure"
        dependencies = ("contacts",)

        def __call__(self):
            return SimpleNamespace(name="exposure")

    class _Cond:
        def __init__(self):
            self.replicate_specs = [SimpleNamespace(replicate=1)]

    calls: list[tuple[str, tuple[str, ...]]] = []

    def _submit_graph(manifest, resources, hpc_dir, root_dependencies=()):
        del resources, hpc_dir
        analysis_name = (
            Path(manifest.comparison_yaml).name if hasattr(manifest, "comparison_yaml") else ""
        )
        calls.append((analysis_name, tuple(root_dependencies)))
        job_id = "101" if len(calls) == 1 else "202"
        return SimpleNamespace(
            finalizer_job_id=job_id, save=lambda path: Path(path).write_text("{}")
        )

    def _build_manifest(plugin, *args, **kwargs):
        del args, kwargs
        return SimpleNamespace(
            comparison_yaml=f"{plugin.name}.yaml",
            condition_specs=[_Cond()],
            pipeline_mode="full",
            save=lambda path: Path(path).write_text("{}"),
        )

    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/sbatch")
    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml", lambda path: config
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.order_analyses_for_execution",
        lambda names, satisfied=None: ["contacts", "exposure"],
    )
    monkeypatch.setattr(
        "polyzymd.analyses.discovery.get_analysis",
        lambda name: _Contacts() if name == "contacts" else _Exposure(),
    )
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm.build_manifest", _build_manifest)
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm.submit_analysis_graph", _submit_graph)

    result = runner.invoke(
        compare,
        ["submit-all", "-f", str(tmp_path / "comparison.yaml")],
    )

    assert result.exit_code == 0
    assert calls[0][1] == ()
    assert calls[1][1] == ("101",)


def test_submit_uses_plugin_memory_hint_when_cli_default(monkeypatch, tmp_path: Path) -> None:
    """compare submit should apply plugin memory hint when --mem is not passed."""
    runner = CliRunner()

    class _Hint:
        mem = "16G"
        time = None
        cpus_per_task = None

    class _FakeAnalysisClass:
        name = "toy"
        dependencies = ()
        slurm_resource_hint = _Hint()

        def __init__(self):
            self.name = "toy"

    class _Cond:
        def __init__(self, reps):
            self.replicate_specs = [SimpleNamespace(replicate=r) for r in reps]

    captured: dict[str, Any] = {}
    manifest = SimpleNamespace(
        condition_specs=[_Cond([1])],
        save=lambda path: Path(path).write_text("{}"),
    )

    def _capture_manifest(plugin, config, resources, *args, **kwargs):
        del plugin, config, args, kwargs
        captured["mem"] = resources.mem
        return manifest

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr("polyzymd.analyses.discovery.get_analysis", lambda name: _FakeAnalysisClass)
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm.build_manifest", _capture_manifest)
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
        ["submit", "toy", "-f", str(tmp_path / "comparison.yaml"), "--dry-run"],
    )

    assert result.exit_code == 0
    assert captured["mem"] == "16G"


def test_submit_cli_memory_overrides_plugin_hint(monkeypatch, tmp_path: Path) -> None:
    """Explicit --mem should override plugin memory hints."""
    runner = CliRunner()

    class _Hint:
        mem = "16G"
        time = None
        cpus_per_task = None

    class _FakeAnalysisClass:
        name = "toy"
        dependencies = ()
        slurm_resource_hint = _Hint()

        def __init__(self):
            self.name = "toy"

    class _Cond:
        def __init__(self, reps):
            self.replicate_specs = [SimpleNamespace(replicate=r) for r in reps]

    captured: dict[str, Any] = {}
    manifest = SimpleNamespace(
        condition_specs=[_Cond([1])],
        save=lambda path: Path(path).write_text("{}"),
    )

    def _capture_manifest(plugin, config, resources, *args, **kwargs):
        del plugin, config, args, kwargs
        captured["mem"] = resources.mem
        return manifest

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr("polyzymd.analyses.discovery.get_analysis", lambda name: _FakeAnalysisClass)
    monkeypatch.setattr("polyzymd.workflow.analysis_slurm.build_manifest", _capture_manifest)
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
            "-f",
            str(tmp_path / "comparison.yaml"),
            "--mem",
            "8G",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert captured["mem"] == "8G"


def test_submit_sets_noisy_loggers_to_warning(monkeypatch, tmp_path: Path) -> None:
    """compare submit should suppress noisy third-party loggers at planning time."""
    runner = CliRunner()

    class _FakeAnalysis:
        name = "toy"
        dependencies = ()

        def __call__(self):
            return SimpleNamespace(name="toy")

    class _Cond:
        def __init__(self, reps):
            self.replicate_specs = [SimpleNamespace(replicate=r) for r in reps]

    manifest = SimpleNamespace(
        condition_specs=[_Cond([1])],
        save=lambda path: Path(path).write_text("{}"),
    )

    mda_logger = logging.getLogger("MDAnalysis")
    numexpr_logger = logging.getLogger("numexpr")
    old_mda = mda_logger.level
    old_numexpr = numexpr_logger.level
    mda_logger.setLevel(logging.INFO)
    numexpr_logger.setLevel(logging.INFO)

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr("polyzymd.analyses.discovery.get_analysis", lambda name: _FakeAnalysis())
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

    try:
        result = runner.invoke(
            compare,
            ["submit", "toy", "-f", str(tmp_path / "comparison.yaml"), "--dry-run"],
        )
        assert result.exit_code == 0
        assert logging.getLogger("MDAnalysis").level == logging.WARNING
        assert logging.getLogger("numexpr").level == logging.WARNING
    finally:
        mda_logger.setLevel(old_mda)
        numexpr_logger.setLevel(old_numexpr)


def test_submit_qos_tip_shown_when_partition_set_without_qos(monkeypatch, tmp_path: Path) -> None:
    """compare submit should show QoS hint when partition is set and qos is omitted."""
    runner = CliRunner()

    class _FakeAnalysis:
        name = "toy"
        dependencies = ()

        def __call__(self):
            return SimpleNamespace(name="toy")

    class _Cond:
        def __init__(self, reps):
            self.replicate_specs = [SimpleNamespace(replicate=r) for r in reps]

    manifest = SimpleNamespace(
        condition_specs=[_Cond([1])],
        save=lambda path: Path(path).write_text("{}"),
    )

    monkeypatch.setattr(
        "polyzymd.config.comparison.ComparisonConfig.from_yaml",
        lambda path: SimpleNamespace(source_path=tmp_path / "comparison.yaml"),
    )
    monkeypatch.setattr("polyzymd.analyses.discovery.get_analysis", lambda name: _FakeAnalysis())
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
            "-f",
            str(tmp_path / "comparison.yaml"),
            "--partition",
            "gpu",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert "Tip: Many HPC clusters require --qos to be set." in result.output
