"""Tests for execution cost hints and local-run warnings."""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar

import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import Analysis, Condition, MetricValue
from polyzymd.analyses.orchestrator import _print_execution_summary, run_comparison


class _HintSettings(BaseModel):
    value: float = 1.0


class _DefaultHintAnalysis(Analysis):
    name: ClassVar[str] = "hint_default"
    Settings: ClassVar[type] = _HintSettings
    min_replicates: ClassVar[int] = 1

    def compute_replicate(self, ctx, replicate: int) -> dict[str, Any]:
        return {"value": float(replicate)}

    def aggregate(self, ctx, results) -> dict[str, Any]:
        values = [float(result["value"]) for result in results]
        return {
            "mean_value": sum(values) / len(values),
            "sem_value": 0.0,
            "replicate_values": values,
        }

    def extract_metrics(self, summary: dict[str, Any]) -> dict[str, MetricValue]:
        return {
            "metric": MetricValue(
                name="metric",
                mean=float(summary["mean_value"]),
                sem=float(summary["sem_value"]),
                replicate_values=[float(v) for v in summary["replicate_values"]],
            )
        }


class _HighHintAnalysis(_DefaultHintAnalysis):
    name: ClassVar[str] = "hint_high"
    execution_cost_hint: ClassVar[str] = "high"


def test_execution_cost_hint_default() -> None:
    """Analysis subclasses should default to medium execution cost."""
    assert _DefaultHintAnalysis.execution_cost_hint == "medium"


@pytest.mark.parametrize(
    "analysis_name",
    ["sasa", "contacts", "binding_free_energy", "exposure", "polymer_affinity"],
)
def test_execution_cost_hint_high(analysis_name: str) -> None:
    """Expensive plugins should declare high execution cost hints."""
    from polyzymd.analyses.discovery import get_analysis

    cls = get_analysis(analysis_name)
    assert cls.execution_cost_hint == "high"


def test_execution_summary_printed(monkeypatch: pytest.MonkeyPatch, caplog, tmp_path: Path) -> None:
    """run_comparison should log an execution summary before compute loop."""
    analysis = _DefaultHintAnalysis()
    conditions = [
        Condition("A", tmp_path / "a.yaml", (1,), SimpleNamespace()),
        Condition("B", tmp_path / "b.yaml", (1,), SimpleNamespace()),
    ]

    config = SimpleNamespace(
        control="A",
        plot_settings=SimpleNamespace(output_dir=tmp_path / "figures"),
        defaults=SimpleNamespace(equilibration_time="10ns"),
        model_copy=lambda deep=True: SimpleNamespace(
            control="A",
            plot_settings=SimpleNamespace(output_dir=tmp_path / "figures"),
            defaults=SimpleNamespace(equilibration_time="10ns"),
        ),
    )

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.prepare_comparison_run",
        lambda analysis, config, equilibration: (
            conditions,
            _HintSettings(),
            "10ns",
            tmp_path / "analysis",
        ),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator._prepare_conditions_with_filter",
        lambda analysis, config, settings: (conditions, conditions, []),
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.run_analysis",
        lambda *args, **kwargs: {
            "mean_value": 1.0,
            "sem_value": 0.0,
            "replicate_values": [1.0],
        },
    )
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.finalize_comparison_from_disk",
        lambda **kwargs: {
            "comparison": {"ok": True},
            "comparison_path": tmp_path / "comparison" / "result.json",
            "plots": [],
        },
    )

    caplog.set_level(logging.INFO)
    run_comparison(analysis, config, recompute=False, equilibration="10ns")
    assert "Mode: sequential (local)" in caplog.text
    assert "total replicate tasks" in caplog.text


def test_slurm_suggestion_when_available(
    monkeypatch: pytest.MonkeyPatch, caplog, tmp_path: Path
) -> None:
    """Warning should suggest compare submit when sbatch is available."""
    analysis = _HighHintAnalysis()
    conditions = [Condition("A", tmp_path / "a.yaml", (1, 2, 3), SimpleNamespace())]
    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.shutil.which", lambda name: "/usr/bin/sbatch"
    )
    caplog.set_level(logging.WARNING)

    _print_execution_summary(analysis, conditions, _HintSettings(), "10ns")
    assert "Consider submitting to SLURM" in caplog.text
    assert "polyzymd compare submit hint_high" in caplog.text


def test_no_slurm_suggestion_when_unavailable(
    monkeypatch: pytest.MonkeyPatch,
    caplog,
    tmp_path: Path,
) -> None:
    """Warning should still suggest HPC path when sbatch is unavailable."""
    analysis = _HighHintAnalysis()
    conditions = [Condition("A", tmp_path / "a.yaml", (1, 2, 3), SimpleNamespace())]
    monkeypatch.setattr("polyzymd.analyses.orchestrator.shutil.which", lambda name: None)
    caplog.set_level(logging.WARNING)

    _print_execution_summary(analysis, conditions, _HintSettings(), "10ns")
    assert "If you have access to an HPC cluster with SLURM" in caplog.text
    assert "polyzymd compare submit hint_high" in caplog.text
