"""Tests for the private one-analysis lifecycle engine."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar, Sequence

import pytest
from pydantic import BaseModel

from polyzymd.analyses._analysis_lifecycle import AnalysisLifecycle
from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.orchestrator import (
    finalize_comparison_from_disk,
    prepare_comparison_run,
    run_analysis,
    run_comparison,
    run_plot_only,
    run_replicate_once,
)
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
        path = ctx.output_dir / "plot.png"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("plot")
        return [path]

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


class _RejectingAnalysis(_DelegatingAnalysis):
    """Analysis that rejects all conditions during filtering."""

    name: ClassVar[str] = "rejecting_lifecycle"

    def filter_conditions(
        self,
        conditions: list[Condition],
        settings: BaseModel | None = None,
    ) -> list[Condition]:
        """Reject every condition to exercise validation errors."""

        del conditions, settings
        self.events.append("filter")
        return []


class _ContractReplicateAnalysis(_OrderAnalysis):
    """Analysis that raises a contract error from the replicate hook."""

    name: ClassVar[str] = "contract_replicate_lifecycle"

    def run_replicate(self, ctx: ReplicateContext, replicate: int) -> dict[str, Any]:
        """Raise a plugin contract error without lifecycle wrapping."""

        del ctx, replicate
        raise PluginContractError("contract boom")


class _CondCfg:
    """Small condition config stand-in for public wrapper tests."""

    def __init__(self, label: str, config: Path, replicates: tuple[int, ...]) -> None:
        self.label = label
        self.config = config
        self.replicates = list(replicates)


class _ComparisonConfig:
    """Small comparison config stand-in for public wrapper tests."""

    def __init__(self, tmp_path: Path, labels: tuple[str, ...] = ("Cond",)) -> None:
        self.name = "project"
        self.source_path = tmp_path / "comparison.yaml"
        self.defaults = SimpleNamespace(equilibration_time="10ns")
        self.control = None
        self.conditions = [
            _CondCfg(label, tmp_path / f"{label.lower()}.yaml", (1, 2))
            for label in labels
        ]
        self.plugins = SimpleNamespace(get=lambda name: None)
        self.plot_settings = PlotSettings(output_dir=tmp_path / "figures")

    def model_copy(self, deep: bool = True) -> "_ComparisonConfig":
        """Return a shallow behavioral copy for lifecycle finalization."""

        del deep
        copied = object.__new__(type(self))
        copied.name = self.name
        copied.source_path = self.source_path
        copied.defaults = SimpleNamespace(
            equilibration_time=self.defaults.equilibration_time,
        )
        copied.control = self.control
        copied.conditions = self.conditions
        copied.plugins = self.plugins
        copied.plot_settings = self.plot_settings
        return copied


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


def _patch_condition_loader(monkeypatch: pytest.MonkeyPatch) -> None:
    """Patch config loading while keeping public lifecycle execution intact."""

    def _from_condition_config(cond_cfg: _CondCfg) -> Condition:
        return Condition(
            cond_cfg.label,
            cond_cfg.config,
            tuple(cond_cfg.replicates),
            SimpleNamespace(),
        )

    monkeypatch.setattr(Condition, "from_condition_config", staticmethod(_from_condition_config))


def test_public_run_analysis_order_and_canonical_save(tmp_path: Path) -> None:
    """Public run_analysis should preserve order, paths, and identity metadata."""

    analysis = _OrderAnalysis()
    condition = _condition(tmp_path)
    settings = _LifecycleSettings(scale=1.5)
    output_dir = tmp_path / "analysis" / analysis.name

    result = run_analysis(
        analysis,
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


def test_public_run_analysis_recompute_removes_aggregate_sidecars(tmp_path: Path) -> None:
    """Public run_analysis should clean analysis-owned aggregate outputs on recompute."""

    class _CleanupAnalysis(_OrderAnalysis):
        name: ClassVar[str] = "cleanup_lifecycle"

        def aggregate(
            self,
            ctx: AggregateContext,
            results: Sequence[dict[str, Any]],
        ) -> dict[str, Any]:
            assert ctx.recompute is True
            assert not (ctx.output_dir / "stale_sidecar.txt").exists()
            return super().aggregate(ctx, results)

    analysis = _CleanupAnalysis()
    condition = _condition(tmp_path)
    output_dir = tmp_path / "analysis" / analysis.name
    stale_sidecar = output_dir / "aggregated" / "stale_sidecar.txt"
    stale_sidecar.parent.mkdir(parents=True, exist_ok=True)
    stale_sidecar.write_text("stale")

    run_analysis(
        analysis,
        condition,
        _LifecycleSettings(),
        equilibration="10ns",
        output_dir=output_dir,
        recompute=True,
    )

    assert not stale_sidecar.exists()


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


def test_public_run_replicate_once_writes_canonical_result(tmp_path: Path) -> None:
    """Public run_replicate_once should execute and persist the canonical result."""

    analysis = _OrderAnalysis()
    condition = _condition(tmp_path, replicates=(3,))
    output_dir = tmp_path / "run_3"

    result = run_replicate_once(
        analysis,
        condition,
        _LifecycleSettings(scale=2.0),
        "5ns",
        output_dir,
        replicate=3,
        recompute=True,
    )

    result_path = output_dir / "result.json"
    assert result == {"value": 6.0, "replicate": 3}
    assert result_path.exists()
    assert '"replicate": 3' in result_path.read_text()
    assert analysis.events == ["run:3"]


def test_public_run_comparison_executes_full_lifecycle(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Public run_comparison should produce concrete aggregate, compare, and plot outputs."""

    _patch_condition_loader(monkeypatch)
    analysis = _DelegatingAnalysis()
    config = _ComparisonConfig(tmp_path)

    result = run_comparison(analysis, config, recompute=False, equilibration="5ns")

    comparison_path = tmp_path / "comparison" / analysis.name / "result.json"
    plot_path = tmp_path / "figures" / analysis.name / "plot.png"
    assert result["aggregated"]["Cond"]["mean_value"] == 1.5
    assert result["comparison"] == {"compared": True}
    assert result["comparison_path"] == comparison_path
    assert result["plots"] == [plot_path]
    assert comparison_path.exists()
    assert plot_path.exists()
    assert analysis.events == [
        "filter",
        "run:1",
        "run:2",
        "aggregate",
        "compare:project",
        f"plot:{analysis.name}",
    ]


def test_public_finalize_loads_aggregate_from_disk_and_writes_outputs(tmp_path: Path) -> None:
    """Finalize should load existing aggregates, then save comparison and plot results."""

    analysis = _DelegatingAnalysis()
    condition = _condition(tmp_path)
    settings = _LifecycleSettings()
    analysis_root = tmp_path / "analysis"
    condition_dir = analysis_root / condition.label / analysis.name
    run_analysis(
        analysis,
        condition,
        settings,
        equilibration="10ns",
        output_dir=condition_dir,
    )
    analysis.events.clear()
    config = _ComparisonConfig(tmp_path)
    prepared_state = {
        "all_conditions": [condition],
        "valid_conditions": [condition],
        "excluded_conditions": [],
        "condition_by_label": {condition.label: condition},
        "settings": settings,
        "equilibration": "10ns",
        "analysis_root": analysis_root,
    }

    result = finalize_comparison_from_disk(
        analysis=analysis,
        config=config,
        analysis_dirs={condition.label: condition_dir},
        aggregated_results={},
        results_dir=tmp_path / "comparison" / analysis.name,
        figures_dir=tmp_path / "figures" / analysis.name,
        settings=settings,
        effective_control=None,
        prepared_state=prepared_state,
    )

    assert result["comparison"] == {"compared": True}
    assert result["comparison_path"].exists()
    assert result["plots"] == [tmp_path / "figures" / analysis.name / "plot.png"]
    assert result["plots"][0].exists()
    assert analysis.events == ["compare:project", f"plot:{analysis.name}"]


def test_public_plot_only_returns_paths_and_failures(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Plot-only wrapper should return generated paths or captured failures."""

    class _FailingPlotAnalysis(_DelegatingAnalysis):
        name: ClassVar[str] = "failing_plot_lifecycle"

        def plot(self, ctx: PlotContext) -> list[Path]:
            """Raise a plotting failure for plot-only failure reporting."""

            del ctx
            raise RuntimeError("plot boom")

    _patch_condition_loader(monkeypatch)
    config = _ComparisonConfig(tmp_path)
    analysis = _DelegatingAnalysis()

    paths, failures = run_plot_only(analysis, config, equilibration="1ns")

    assert failures == []
    assert paths == [tmp_path / "figures" / analysis.name / "plot.png"]
    assert paths[0].exists()
    assert analysis.events == ["filter", f"plot:{analysis.name}"]

    failing_paths, failing_failures = run_plot_only(
        _FailingPlotAnalysis(),
        config,
        equilibration="1ns",
    )
    assert failing_paths == []
    assert failing_failures == [("failing_plot_lifecycle", "plot boom")]


def test_public_prepare_validation_error_uses_lifecycle_filter(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Prepare should raise a concrete error when filtering removes every condition."""

    _patch_condition_loader(monkeypatch)

    with pytest.raises(ValueError, match="no valid conditions remain"):
        prepare_comparison_run(_RejectingAnalysis(), _ComparisonConfig(tmp_path), "10ns")


def test_public_wrapper_preserves_contract_exceptions(tmp_path: Path) -> None:
    """Public lifecycle wrappers should propagate plugin contract errors."""

    with pytest.raises(PluginContractError, match="contract boom"):
        run_analysis(
            _ContractReplicateAnalysis(),
            _condition(tmp_path, replicates=(1,)),
            _LifecycleSettings(),
            output_dir=tmp_path / "analysis" / "contract",
        )
