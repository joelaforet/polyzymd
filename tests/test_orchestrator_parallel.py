"""Tests for extracted orchestrator parallel helpers."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar

from pydantic import BaseModel

from polyzymd.analyses.base import Analysis, Condition, MetricValue
from polyzymd.analyses.orchestrator import (
    aggregate_condition_from_disk,
    prepare_comparison_run,
    run_analysis,
    run_comparison,
    run_replicate_once,
)


class _ParallelSettings(BaseModel):
    factor: float = 1.0


class _ParallelAnalysis(Analysis):
    name: ClassVar[str] = "parallel_toy"
    Settings: ClassVar[type] = _ParallelSettings
    min_replicates: ClassVar[int] = 1

    def compute_replicate(self, ctx, replicate: int) -> dict[str, Any]:
        return {"value": float(replicate) * float(ctx.settings.factor), "replicate": replicate}

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
            "metric": MetricValue(
                name="metric",
                mean=float(summary["mean_value"]),
                sem=float(summary["sem_value"]),
                replicate_values=[float(v) for v in summary["replicate_values"]],
                higher_is_better=True,
            )
        }

    def plot(self, ctx):
        out = ctx.output_dir / "parallel_toy.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("plot")
        return [out]


class _CondCfg:
    def __init__(self, label: str, config: Path, replicates: tuple[int, ...]):
        self.label = label
        self.config = config
        self.replicates = list(replicates)


def _make_config(tmp_path: Path):
    return SimpleNamespace(
        name="parallel_project",
        source_path=tmp_path / "comparison.yaml",
        defaults=SimpleNamespace(equilibration_time="10ns"),
        control="A",
        conditions=[
            _CondCfg("A", tmp_path / "a.yaml", (1, 2)),
            _CondCfg("B", tmp_path / "b.yaml", (1, 2)),
        ],
        plugins=SimpleNamespace(get=lambda name: None),
        plot_settings=SimpleNamespace(output_dir=tmp_path / "figures"),
        model_copy=lambda deep=True: SimpleNamespace(
            name="parallel_project",
            source_path=tmp_path / "comparison.yaml",
            defaults=SimpleNamespace(equilibration_time="10ns"),
            control="A",
            conditions=[
                _CondCfg("A", tmp_path / "a.yaml", (1, 2)),
                _CondCfg("B", tmp_path / "b.yaml", (1, 2)),
            ],
            plugins=SimpleNamespace(get=lambda name: None),
            plot_settings=SimpleNamespace(output_dir=tmp_path / "figures"),
        ),
    )


def test_prepare_comparison_run_returns_conditions_settings(monkeypatch, tmp_path: Path) -> None:
    """prepare_comparison_run should resolve conditions and settings."""
    analysis = _ParallelAnalysis()
    config = _make_config(tmp_path)

    def _from_cond(cond_cfg):
        return Condition(
            cond_cfg.label, cond_cfg.config, tuple(cond_cfg.replicates), SimpleNamespace()
        )

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.Condition.from_condition_config", _from_cond
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator._resolve_settings",
        lambda analysis, config: _ParallelSettings(factor=2.0),
    )

    conditions, settings, equilibration, analysis_root = prepare_comparison_run(
        analysis,
        config,
        "5ns",
    )
    assert [condition.label for condition in conditions] == ["A", "B"]
    assert settings.factor == 2.0
    assert equilibration == "5ns"
    assert analysis_root == tmp_path / "analysis"


def test_run_analysis_and_helper_path_equivalence(tmp_path: Path) -> None:
    """run_analysis and helper pair should produce equivalent aggregates."""
    analysis = _ParallelAnalysis()
    condition = Condition("A", tmp_path / "a.yaml", (1, 2), SimpleNamespace())
    settings = _ParallelSettings(factor=1.0)

    seq_dir = tmp_path / "seq"
    helper_dir = tmp_path / "helper"

    seq = run_analysis(
        analysis,
        condition,
        settings,
        equilibration="10ns",
        output_dir=seq_dir,
        recompute=False,
    )

    run_replicate_once(analysis, condition, settings, "10ns", helper_dir / "run_1", 1, False)
    run_replicate_once(analysis, condition, settings, "10ns", helper_dir / "run_2", 2, False)
    helper = aggregate_condition_from_disk(
        analysis,
        condition,
        settings,
        "10ns",
        helper_dir,
        replicates=(1, 2),
    )

    assert seq["mean_value"] == helper["mean_value"]
    assert seq["n_replicates"] == helper["n_replicates"]


def test_run_comparison_still_works_after_refactor(monkeypatch, tmp_path: Path) -> None:
    """run_comparison should still produce aggregated and comparison output."""
    analysis = _ParallelAnalysis()
    config = _make_config(tmp_path)

    def _from_cond(cond_cfg):
        return Condition(
            cond_cfg.label, cond_cfg.config, tuple(cond_cfg.replicates), SimpleNamespace()
        )

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.Condition.from_condition_config", _from_cond
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator._resolve_settings",
        lambda analysis, config: _ParallelSettings(),
    )

    result = run_comparison(analysis, config, recompute=False, equilibration="10ns")
    assert "A" in result["aggregated"]
    assert "B" in result["aggregated"]
    assert result["comparison"] is not None
    assert result["comparison_path"].exists()
