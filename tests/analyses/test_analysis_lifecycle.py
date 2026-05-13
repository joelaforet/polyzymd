"""Tests for the private one-analysis lifecycle engine."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar, Sequence

import pytest
from pydantic import BaseModel

from polyzymd.analyses._analysis_lifecycle import (
    ONE_ANALYSIS_STAGE_ORDER,
    AnalysisCompatibilityAdapter,
    AnalysisLifecycle,
)
from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.orchestrator import run_replicate_once
from polyzymd.config.comparison import PlotSettings


class _LifecycleSettings(BaseModel):
    """Minimal settings model for lifecycle tests."""

    scale: float = 1.0


class _OrderAnalysis(Analysis):
    """Analysis that records compute and aggregate lifecycle calls."""

    name: ClassVar[str] = "lifecycle_order"
    Settings: ClassVar[type] = _LifecycleSettings
    min_replicates: ClassVar[int] = 1

    def __init__(self) -> None:
        self.events: list[str] = []

    def run_replicate(self, ctx: ReplicateContext, replicate: int) -> dict[str, Any]:
        """Record replicate execution and return a scalar payload."""

        self.events.append(f"run:{replicate}")
        return {"value": float(replicate) * ctx.settings.scale, "replicate": replicate}

    def aggregate(self, ctx: AggregateContext, results: Sequence[dict[str, Any]]) -> dict[str, Any]:
        """Record aggregation and return an aggregate payload."""

        self.events.append("aggregate")
        values = [float(result["value"]) for result in results]
        return {
            "mean_value": sum(values) / len(values),
            "sem_value": 0.0,
            "replicate_values": values,
        }


class _DelegatingAnalysis(_OrderAnalysis):
    """Analysis that records every compatibility adapter hook."""

    name: ClassVar[str] = "lifecycle_delegate"

    def filter_conditions(
        self,
        conditions: list[Condition],
        settings: BaseModel | None = None,
    ) -> list[Condition]:
        """Record filter delegation and return original conditions."""

        del settings
        self.events.append("filter")
        return list(conditions)

    def compare(self, ctx: ComparisonContext) -> dict[str, Any]:
        """Record compare delegation and return a payload."""

        self.events.append(f"compare:{ctx.name}")
        return {"compared": True}

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Record plot delegation and return a synthetic path."""

        self.events.append(f"plot:{ctx.output_dir.name}")
        return [ctx.output_dir / "plot.png"]

    def format(self, result: Any, output_format: str = "text") -> str:
        """Record format delegation and return formatted text."""

        self.events.append(f"format:{output_format}")
        return f"{output_format}:{result}"


class _InvalidAggregateAnalysis(_OrderAnalysis):
    """Analysis that violates the aggregate return contract."""

    name: ClassVar[str] = "invalid_lifecycle_aggregate"

    min_replicates: ClassVar[int] = 1

    Settings: ClassVar[type] = _LifecycleSettings

    def aggregate(self, ctx: AggregateContext, results: Sequence[dict[str, Any]]) -> list[str]:
        """Return an invalid aggregate type for validation tests."""

        del ctx, results
        return ["invalid"]


def _condition(tmp_path: Path, replicates: tuple[int, ...] = (1, 2)) -> Condition:
    """Build a lightweight condition for lifecycle tests.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory for synthetic config paths.
    replicates : tuple[int, ...], optional
        Replicate IDs for the condition, by default ``(1, 2)``.

    Returns
    -------
    Condition
        Synthetic condition with a simple config object.
    """

    return Condition("Cond", tmp_path / "cond.yaml", replicates, SimpleNamespace())


def test_stage_order_constant_documents_one_analysis_order() -> None:
    """Stage-order constants should expose the Template Method sequence."""

    assert ONE_ANALYSIS_STAGE_ORDER == ("prepare", "compute", "aggregate", "compare", "plot")


def test_adapter_delegates_existing_analysis_hooks(tmp_path: Path) -> None:
    """Compatibility adapter should delegate to existing Analysis hooks."""

    analysis = _DelegatingAnalysis()
    adapter = AnalysisCompatibilityAdapter(analysis)
    settings = _LifecycleSettings(scale=2.0)
    condition = _condition(tmp_path, replicates=(1,))
    rep_ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="10ns",
        recompute=False,
        settings=settings,
    )
    agg_ctx = AggregateContext(
        condition=condition,
        replicates=(1,),
        output_dir=tmp_path / "aggregated",
        equilibration="10ns",
        settings=settings,
    )
    comp_ctx = ComparisonContext(
        name="project",
        conditions=[condition],
        excluded_conditions=[],
        control_label=None,
        analysis_dirs={"Cond": tmp_path},
        results_dir=tmp_path / "comparison",
        equilibration="10ns",
        settings=settings,
    )
    plot_ctx = PlotContext(
        conditions=[condition],
        analysis_dirs={"Cond": tmp_path},
        results_dir=tmp_path / "comparison",
        output_dir=tmp_path / "figures",
        settings=settings,
        plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
    )

    assert adapter.filter_conditions([condition], settings=settings) == [condition]
    rep_result = adapter.run_replicate(rep_ctx, 1)
    assert rep_result["value"] == 2.0
    assert adapter.aggregate(agg_ctx, [rep_result])["mean_value"] == 2.0
    assert adapter.compare(comp_ctx) == {"compared": True}
    assert adapter.plot(plot_ctx) == [tmp_path / "figures" / "plot.png"]
    assert adapter.format({"ok": True}, fmt="json") == "json:{'ok': True}"
    assert analysis.events == [
        "filter",
        "run:1",
        "aggregate",
        "compare:project",
        "plot:figures",
        "format:json",
    ]


def test_lifecycle_run_analysis_order_and_canonical_save(tmp_path: Path) -> None:
    """Lifecycle run_analysis should preserve order, paths, and identity metadata."""

    analysis = _OrderAnalysis()
    lifecycle = AnalysisLifecycle(analysis)
    condition = _condition(tmp_path)
    settings = _LifecycleSettings(scale=1.5)
    output_dir = tmp_path / "analysis" / analysis.name

    result = lifecycle.run_analysis(
        condition,
        settings,
        equilibration="10ns",
        output_dir=output_dir,
        recompute=False,
    )

    assert analysis.events == ["run:1", "run:2", "aggregate"]
    assert result["replicate_values"] == [1.5, 3.0]
    assert result["replicates"] == [1, 2]
    assert result["n_replicates"] == 2
    assert result["settings_fingerprint"] == analysis.aggregate_settings_fingerprint(settings)
    assert (output_dir / "run_1" / "result.json").exists()
    assert (output_dir / "run_2" / "result.json").exists()
    assert (output_dir / "aggregated" / "result.json").exists()


def test_lifecycle_validates_and_rejects_invalid_aggregate(tmp_path: Path) -> None:
    """Lifecycle should preserve aggregate contract validation behavior."""

    analysis = _InvalidAggregateAnalysis()
    lifecycle = AnalysisLifecycle(analysis)
    condition = _condition(tmp_path, replicates=(1,))

    with pytest.raises(PluginContractError, match="invalid_lifecycle_aggregate.aggregate"):
        lifecycle.run_analysis(
            condition,
            _LifecycleSettings(),
            equilibration="0ns",
            output_dir=tmp_path / "analysis" / analysis.name,
        )


def test_orchestrator_wrapper_delegates_to_lifecycle(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Public wrappers should remain thin delegations to the private lifecycle."""

    calls: list[tuple[str, int, bool]] = []

    def _fake_run_replicate_once(
        self: AnalysisLifecycle,
        condition: Condition,
        settings: BaseModel,
        equilibration: str,
        output_dir: Path,
        replicate: int,
        recompute: bool,
    ) -> dict[str, Any]:
        """Record public wrapper delegation arguments."""

        del self, condition, settings, equilibration, output_dir
        calls.append(("run_replicate_once", replicate, recompute))
        return {"replicate": replicate}

    monkeypatch.setattr(AnalysisLifecycle, "run_replicate_once", _fake_run_replicate_once)

    result = run_replicate_once(
        _OrderAnalysis(),
        _condition(tmp_path, replicates=(3,)),
        _LifecycleSettings(),
        "5ns",
        tmp_path / "run_3",
        replicate=3,
        recompute=True,
    )

    assert result == {"replicate": 3}
    assert calls == [("run_replicate_once", 3, True)]
