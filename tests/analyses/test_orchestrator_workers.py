"""Tests for orchestrator worker helper functions."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar, cast

import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import Analysis, Condition, MetricValue
from polyzymd.analyses.exceptions import ReplicateError
from polyzymd.analyses.orchestrator import (
    aggregate_condition_from_disk,
    finalize_comparison_from_disk,
    run_analysis,
    run_replicate_once,
)
from polyzymd.config.comparison import PlotSettings


class _WorkerSettings(BaseModel):
    scale: float = 1.0


class _WorkerAnalysis(Analysis):
    name: ClassVar[str] = "worker_toy"
    Settings: ClassVar[type] = _WorkerSettings
    min_replicates: ClassVar[int] = 1

    def compute_replicate(self, ctx: Any, replicate: int) -> dict[str, Any]:
        return {"value": float(replicate) * float(ctx.settings.scale), "replicate": replicate}

    def aggregate(self, ctx, results) -> dict[str, Any]:
        values = [float(result["value"]) for result in results]
        mean = sum(values) / len(values)
        return {
            "mean_value": mean,
            "sem_value": 0.0,
            "replicate_values": values,
            "n_replicates": len(values),
        }

    def extract_metrics(self, summary: dict[str, Any]) -> dict[str, MetricValue]:
        return {
            "worker_metric": MetricValue(
                name="worker_metric",
                mean=float(summary["mean_value"]),
                sem=float(summary["sem_value"]),
                replicate_values=[float(v) for v in summary["replicate_values"]],
                higher_is_better=True,
            )
        }

    def plot(self, ctx):
        out = ctx.output_dir / "worker_toy_plot.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("plot")
        return [out]


class _FailingWorkerAnalysis(_WorkerAnalysis):
    def compute_replicate(self, ctx: Any, replicate: int) -> dict[str, Any]:
        raise RuntimeError("test error")


def test_run_replicate_once_saves_canonical_result(tmp_path: Path) -> None:
    """run_replicate_once should save result.json in the run directory."""
    analysis = _WorkerAnalysis()
    condition = Condition("Cond", tmp_path / "cfg.yaml", (1,), cast(Any, SimpleNamespace()))
    settings = _WorkerSettings(scale=2.0)
    run_dir = tmp_path / "analysis" / "cond" / "worker_toy" / "run_1"

    result = run_replicate_once(
        analysis,
        condition,
        settings,
        "10ns",
        run_dir,
        replicate=1,
        recompute=False,
    )
    assert result["value"] == 2.0
    assert (run_dir / "result.json").exists()


def test_run_analysis_raises_structured_replicate_error_on_worker_exception(tmp_path: Path) -> None:
    """run_analysis should raise ReplicateError on unexpected worker failures."""
    analysis = _FailingWorkerAnalysis()
    condition = Condition("Cond", tmp_path / "cfg.yaml", (1,), cast(Any, SimpleNamespace()))
    settings = _WorkerSettings(scale=1.0)

    with pytest.raises(ReplicateError, match="compute_replicate failed"):
        run_analysis(
            analysis=analysis,
            condition=condition,
            settings=settings,
            equilibration="10ns",
            output_dir=tmp_path / "analysis" / "cond" / "worker_toy",
            recompute=False,
        )


def test_aggregate_condition_from_disk_loads_replicates(tmp_path: Path) -> None:
    """aggregate_condition_from_disk should load per-run result.json files."""
    analysis = _WorkerAnalysis()
    condition = Condition("Cond", tmp_path / "cfg.yaml", (1, 2), cast(Any, SimpleNamespace()))
    settings = _WorkerSettings()
    cond_dir = tmp_path / "analysis" / "cond" / "worker_toy"

    run_replicate_once(analysis, condition, settings, "10ns", cond_dir / "run_1", 1, False)
    run_replicate_once(analysis, condition, settings, "10ns", cond_dir / "run_2", 2, False)

    aggregated = aggregate_condition_from_disk(
        analysis,
        condition,
        settings,
        "10ns",
        cond_dir,
        replicates=(1, 2),
    )
    assert aggregated["n_replicates"] == 2
    assert (cond_dir / "aggregated" / "result.json").exists()


def test_aggregate_condition_from_disk_partial_missing(tmp_path: Path) -> None:
    """Aggregator should tolerate missing replicate files when min is met."""
    analysis = _WorkerAnalysis()
    condition = Condition("Cond", tmp_path / "cfg.yaml", (1, 2), cast(Any, SimpleNamespace()))
    settings = _WorkerSettings()
    cond_dir = tmp_path / "analysis" / "cond" / "worker_toy"

    run_replicate_once(analysis, condition, settings, "10ns", cond_dir / "run_1", 1, False)
    aggregated = aggregate_condition_from_disk(
        analysis,
        condition,
        settings,
        "10ns",
        cond_dir,
        replicates=(1, 2),
    )
    assert aggregated["n_replicates"] == 1


def test_finalize_comparison_from_disk_runs_compare_and_plot(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """Finalizer helper should save comparison result and figures."""
    analysis = _WorkerAnalysis()
    cond_a = Condition("A", tmp_path / "a.yaml", (1,), cast(Any, SimpleNamespace()))
    cond_b = Condition("B", tmp_path / "b.yaml", (1,), cast(Any, SimpleNamespace()))

    class _CondCfg:
        def __init__(self, label: str):
            self.label = label

    config = SimpleNamespace(
        name="proj",
        control=None,
        conditions=[_CondCfg("A"), _CondCfg("B")],
        defaults=SimpleNamespace(equilibration_time="10ns"),
        plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
    )

    def _fake_from_cond(cond_cfg):
        return cond_a if cond_cfg.label == "A" else cond_b

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.Condition.from_condition_config", _fake_from_cond
    )

    analysis_dirs = {
        "A": tmp_path / "analysis" / "A" / "worker_toy",
        "B": tmp_path / "analysis" / "B" / "worker_toy",
    }
    aggregated = {
        "A": {"mean_value": 1.0, "sem_value": 0.0, "replicate_values": [1.0], "n_replicates": 1},
        "B": {"mean_value": 2.0, "sem_value": 0.0, "replicate_values": [2.0], "n_replicates": 1},
    }
    out = finalize_comparison_from_disk(
        analysis=analysis,
        config=cast(Any, config),
        analysis_dirs=analysis_dirs,
        aggregated_results=aggregated,
        results_dir=tmp_path / "comparison" / "worker_toy",
        figures_dir=tmp_path / "figures" / "worker_toy",
        settings=_WorkerSettings(),
        effective_control="A",
    )
    assert out["comparison_path"].exists()
    assert out["plots"]


def test_finalize_partial_recomputes_control_when_original_missing(
    monkeypatch,
    caplog,
    tmp_path: Path,
) -> None:
    """Partial finalize should clear control when configured control is dropped."""
    analysis = _WorkerAnalysis()
    cond_a = Condition("A", tmp_path / "a.yaml", (1,), cast(Any, SimpleNamespace()))
    cond_b = Condition("B", tmp_path / "b.yaml", (1,), cast(Any, SimpleNamespace()))
    cond_c = Condition("C", tmp_path / "c.yaml", (1,), cast(Any, SimpleNamespace()))

    class _CondCfg:
        def __init__(self, label: str):
            self.label = label

    config = SimpleNamespace(
        name="proj",
        control="A",
        conditions=[_CondCfg("A"), _CondCfg("B"), _CondCfg("C")],
        defaults=SimpleNamespace(equilibration_time="10ns"),
        plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
    )

    def _fake_from_cond(cond_cfg):
        mapping = {"A": cond_a, "B": cond_b, "C": cond_c}
        return mapping[cond_cfg.label]

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.Condition.from_condition_config", _fake_from_cond
    )

    analysis_dirs = {
        "B": tmp_path / "analysis" / "B" / "worker_toy",
        "C": tmp_path / "analysis" / "C" / "worker_toy",
    }
    aggregated = {
        "B": {"mean_value": 2.0, "sem_value": 0.0, "replicate_values": [2.0], "n_replicates": 1},
        "C": {"mean_value": 3.0, "sem_value": 0.0, "replicate_values": [3.0], "n_replicates": 1},
    }

    caplog.set_level("WARNING")
    out = finalize_comparison_from_disk(
        analysis=analysis,
        config=cast(Any, config),
        analysis_dirs=analysis_dirs,
        aggregated_results=aggregated,
        results_dir=tmp_path / "comparison" / "worker_toy",
        figures_dir=tmp_path / "figures" / "worker_toy",
        settings=_WorkerSettings(),
        effective_control="A",
        allow_partial=True,
    )
    assert out["comparison"] is not None
    assert out["comparison"].control_label is None
    assert "comparison will proceed without a designated control (all-vs-all)" in caplog.text


def test_finalize_partial_with_one_condition_succeeds(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """Partial finalize should succeed when one condition survives."""
    analysis = _WorkerAnalysis()
    cond_a = Condition("A", tmp_path / "a.yaml", (1,), cast(Any, SimpleNamespace()))
    cond_b = Condition("B", tmp_path / "b.yaml", (1,), cast(Any, SimpleNamespace()))

    class _CondCfg:
        def __init__(self, label: str):
            self.label = label

    config = SimpleNamespace(
        name="proj",
        control="A",
        conditions=[_CondCfg("A"), _CondCfg("B")],
        defaults=SimpleNamespace(equilibration_time="10ns"),
        plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
    )

    def _fake_from_cond(cond_cfg):
        return cond_a if cond_cfg.label == "A" else cond_b

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.Condition.from_condition_config", _fake_from_cond
    )

    out = finalize_comparison_from_disk(
        analysis=analysis,
        config=cast(Any, config),
        analysis_dirs={"B": tmp_path / "analysis" / "B" / "worker_toy"},
        aggregated_results={
            "B": {
                "mean_value": 2.0,
                "sem_value": 0.0,
                "replicate_values": [2.0],
                "n_replicates": 1,
            }
        },
        results_dir=tmp_path / "comparison" / "worker_toy",
        figures_dir=tmp_path / "figures" / "worker_toy",
        settings=_WorkerSettings(),
        effective_control="A",
        allow_partial=True,
    )
    assert out["comparison"] is not None
    assert out["comparison"].conditions[0].label == "B"


def test_finalize_partial_raises_with_zero_successful_conditions(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """Partial finalize should fail when no conditions survive."""
    analysis = _WorkerAnalysis()
    cond_a = Condition("A", tmp_path / "a.yaml", (1,), cast(Any, SimpleNamespace()))
    cond_b = Condition("B", tmp_path / "b.yaml", (1,), cast(Any, SimpleNamespace()))

    class _CondCfg:
        def __init__(self, label: str):
            self.label = label

    config = SimpleNamespace(
        name="proj",
        control="A",
        conditions=[_CondCfg("A"), _CondCfg("B")],
        defaults=SimpleNamespace(equilibration_time="10ns"),
        plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
    )

    def _fake_from_cond(cond_cfg):
        return cond_a if cond_cfg.label == "A" else cond_b

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.Condition.from_condition_config", _fake_from_cond
    )

    with pytest.raises(ValueError, match="no successful conditions remain"):
        finalize_comparison_from_disk(
            analysis=analysis,
            config=cast(Any, config),
            analysis_dirs={},
            aggregated_results={},
            results_dir=tmp_path / "comparison" / "worker_toy",
            figures_dir=tmp_path / "figures" / "worker_toy",
            settings=_WorkerSettings(),
            effective_control="A",
            allow_partial=True,
        )
