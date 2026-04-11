"""Tests for extracted orchestrator parallel helpers."""

from __future__ import annotations

import statistics
from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar, Sequence
from unittest.mock import MagicMock

import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    Condition,
    MetricValue,
    ReplicateContext,
)
from polyzymd.analyses.exceptions import (
    PluginContractError,
    ReplicateError,
    ReplicateSkippedError,
)
from polyzymd.analyses.orchestrator import (
    _run_replicates,
    aggregate_condition_from_disk,
    prepare_comparison_run,
    run_analysis,
    run_comparison,
    run_replicate_once,
)
from polyzymd.config.comparison import PlotSettings


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


class ToySettings(BaseModel):
    """Minimal settings for orchestrator unit tests."""

    threshold: float = 1.0


class ToyResult(BaseModel):
    """Minimal replicate result for orchestrator unit tests."""

    value: float
    replicate: int


class ToyAggregatedResult(BaseModel):
    """Minimal aggregated result for orchestrator unit tests."""

    mean_value: float
    sem_value: float
    replicate_values: list[float]
    n_replicates: int


class ToyAnalysis(Analysis):
    """Concrete analysis for orchestrator tests."""

    name: ClassVar[str] = "toy"
    Settings: ClassVar[type] = ToySettings
    AggregatedResultClass: ClassVar[type] = ToyAggregatedResult
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> ToyResult:
        return ToyResult(value=replicate * 1.5, replicate=replicate)

    def aggregate(self, ctx: AggregateContext, results: Sequence[ToyResult]) -> ToyAggregatedResult:
        values = [result.value for result in results]
        mean_val = statistics.mean(values)
        sem_val = statistics.stdev(values) / len(values) ** 0.5 if len(values) > 1 else 0.0
        return ToyAggregatedResult(
            mean_value=mean_val,
            sem_value=sem_val,
            replicate_values=values,
            n_replicates=len(values),
        )


class ToyDependentAnalysis(Analysis):
    """Analysis that depends on ToyAnalysis."""

    name: ClassVar[str] = "toy_dependent"
    Settings: ClassVar[type] = ToySettings
    dependencies: ClassVar[tuple[str, ...]] = ("toy",)

    def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> ToyResult:
        return ToyResult(value=replicate * 2.0, replicate=replicate)

    def aggregate(self, ctx: AggregateContext, results: Sequence[ToyResult]) -> ToyAggregatedResult:
        values = [result.value for result in results]
        return ToyAggregatedResult(
            mean_value=sum(values) / len(values),
            sem_value=0.0,
            replicate_values=values,
            n_replicates=len(values),
        )


@pytest.fixture
def toy_analysis():
    return ToyAnalysis()


@pytest.fixture
def toy_condition(tmp_path):
    return Condition(
        label="Test Condition",
        config_path=tmp_path / "config.yaml",
        replicates=(1, 2, 3),
        sim_config=MagicMock(),
    )


@pytest.fixture
def toy_settings():
    return ToySettings(threshold=2.0)


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
        plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
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
            plot_settings=PlotSettings(output_dir=tmp_path / "figures"),
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


class TestOrchestrator:
    """Test the orchestrator's replicate running and dependency sorting."""

    def test_run_replicates_success(self, toy_analysis, toy_condition, toy_settings, tmp_path):
        results, successful, failed = _run_replicates(
            toy_analysis,
            toy_condition,
            toy_settings,
            equilibration="10ns",
            output_dir=tmp_path,
        )
        assert len(results) == 3
        assert successful == [1, 2, 3]
        assert failed == []
        assert (tmp_path / "run_1").exists()
        assert (tmp_path / "run_2").exists()
        assert (tmp_path / "run_3").exists()

    def test_run_replicates_partial_failure(self, toy_condition, toy_settings, tmp_path):
        """If some replicates fail but min_replicates is met, continue."""

        class FailingAnalysis(Analysis):
            name: ClassVar[str] = "failing"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 1

            def compute_replicate(self, ctx, replicate):
                if replicate == 2:
                    raise FileNotFoundError("Trajectory not found")
                return ToyResult(value=float(replicate), replicate=replicate)

            def aggregate(self, ctx, results):
                return None

        analysis = FailingAnalysis()
        results, successful, failed = _run_replicates(
            analysis,
            toy_condition,
            toy_settings,
            equilibration="0ns",
            output_dir=tmp_path,
        )
        assert len(results) == 2
        assert 2 not in successful
        assert 2 in failed

    def test_run_replicates_skipped_error_is_recoverable(
        self,
        toy_condition,
        toy_settings,
        tmp_path,
        caplog,
    ):
        """ReplicateSkippedError should skip replicate without aborting condition."""

        class SkipReplicateAnalysis(Analysis):
            name: ClassVar[str] = "skip_replicate"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 1

            def compute_replicate(self, ctx, replicate):
                if replicate == 2:
                    raise ReplicateSkippedError("No trajectory data found for replicate 2")
                return ToyResult(value=float(replicate), replicate=replicate)

            def aggregate(self, ctx, results):
                return {"ok": True}

        caplog.set_level("WARNING")
        results, successful, failed = _run_replicates(
            SkipReplicateAnalysis(),
            toy_condition,
            toy_settings,
            equilibration="0ns",
            output_dir=tmp_path,
        )
        assert len(results) == 2
        assert successful == [1, 3]
        assert failed == [2]
        assert "No trajectory data found for replicate 2" in caplog.text

    def test_run_replicates_below_minimum_raises(self, toy_condition, toy_settings, tmp_path):
        """If fewer than min_replicates succeed, raise ValueError."""

        class AlwaysFailAnalysis(Analysis):
            name: ClassVar[str] = "always_fail"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 2

            def compute_replicate(self, ctx, replicate):
                raise FileNotFoundError("Missing trajectory")

            def aggregate(self, ctx, results):
                return None

        with pytest.raises(ValueError, match="need at least 2"):
            _run_replicates(
                AlwaysFailAnalysis(),
                toy_condition,
                toy_settings,
                equilibration="0ns",
                output_dir=tmp_path,
            )

    def test_run_replicates_all_skipped_fails_minimum(self, toy_condition, toy_settings, tmp_path):
        """When all replicates skip, minimum replicate validation should fail."""

        class AllSkippedAnalysis(Analysis):
            name: ClassVar[str] = "all_skipped"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 1

            def compute_replicate(self, ctx, replicate):
                raise ReplicateSkippedError(f"No trajectory data found for replicate {replicate}")

            def aggregate(self, ctx, results):
                return {"ok": True}

        with pytest.raises(ValueError, match="need at least 1"):
            _run_replicates(
                AllSkippedAnalysis(),
                toy_condition,
                toy_settings,
                equilibration="0ns",
                output_dir=tmp_path,
            )

    def test_run_replicates_unexpected_failure_raises_replicate_error(
        self, toy_condition, toy_settings, tmp_path
    ):
        """Unexpected compute failures should raise structured ReplicateError."""

        class ExplodingAnalysis(Analysis):
            name: ClassVar[str] = "exploding"
            Settings: ClassVar[type] = ToySettings

            def compute_replicate(self, ctx, replicate):
                raise RuntimeError("boom")

            def aggregate(self, ctx, results):
                return None

        with pytest.raises(ReplicateError, match="condition='Test Condition' replicate=1"):
            _run_replicates(
                ExplodingAnalysis(),
                toy_condition,
                toy_settings,
                equilibration="0ns",
                output_dir=tmp_path,
            )

    def test_run_analysis_rejects_invalid_compute_return_type(self, toy_condition, tmp_path):
        """Invalid compute return types should fail plugin contract validation."""

        class InvalidComputeAnalysis(Analysis):
            name: ClassVar[str] = "invalid_compute"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 1

            def compute_replicate(self, ctx, replicate):
                return ["not", "valid"]

            def aggregate(self, ctx, results):
                return {"ok": True}

        with pytest.raises(PluginContractError, match="invalid_compute.compute_replicate"):
            run_analysis(
                InvalidComputeAnalysis(),
                toy_condition,
                ToySettings(),
                equilibration="0ns",
                output_dir=tmp_path,
            )

    def test_run_analysis_rejects_invalid_aggregate_return_type(self, toy_condition, tmp_path):
        """Invalid aggregate return types should fail plugin contract validation."""

        class InvalidAggregateAnalysis(Analysis):
            name: ClassVar[str] = "invalid_aggregate"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 1

            def compute_replicate(self, ctx, replicate):
                return {"replicate": replicate}

            def aggregate(self, ctx, results):
                return "not-a-dict-or-basemodel"

        with pytest.raises(PluginContractError, match="invalid_aggregate.aggregate"):
            run_analysis(
                InvalidAggregateAnalysis(),
                toy_condition,
                ToySettings(),
                equilibration="0ns",
                output_dir=tmp_path,
            )

    def test_topological_sort_no_deps(self):
        from polyzymd.analyses.orchestrator import _topological_sort

        analysis = ToyAnalysis()
        result = _topological_sort([analysis])
        assert result == [analysis]

    def test_topological_sort_with_deps(self):
        from polyzymd.analyses.orchestrator import _topological_sort

        toy = ToyAnalysis()
        dependent = ToyDependentAnalysis()
        result = _topological_sort([dependent, toy])
        assert result.index(toy) < result.index(dependent)

    def test_topological_sort_circular_raises(self):
        from polyzymd.analyses.orchestrator import _topological_sort

        class CircA(Analysis):
            name: ClassVar[str] = "circ_a"
            Settings: ClassVar[type] = ToySettings
            dependencies: ClassVar[tuple[str, ...]] = ("circ_b",)

            def compute_replicate(self, ctx, replicate):
                pass

            def aggregate(self, ctx, results):
                pass

        class CircB(Analysis):
            name: ClassVar[str] = "circ_b"
            Settings: ClassVar[type] = ToySettings
            dependencies: ClassVar[tuple[str, ...]] = ("circ_a",)

            def compute_replicate(self, ctx, replicate):
                pass

            def aggregate(self, ctx, results):
                pass

        with pytest.raises(ValueError, match="Circular dependency"):
            _topological_sort([CircA(), CircB()])

    def test_run_analysis_full(self, toy_analysis, toy_condition, toy_settings, tmp_path):
        """Test the full run_analysis path (compute + aggregate)."""
        result = run_analysis(
            toy_analysis,
            toy_condition,
            toy_settings,
            equilibration="10ns",
            output_dir=tmp_path,
        )
        assert isinstance(result, ToyAggregatedResult)
        assert result.n_replicates == 3
        assert len(result.replicate_values) == 3
        assert (tmp_path / "aggregated").exists()
        assert (tmp_path / "run_1" / "result.json").exists()
        assert (tmp_path / "aggregated" / "result.json").exists()

    def test_compare_only_analysis_skips_compute_and_aggregate(
        self, toy_condition, toy_settings, tmp_path
    ):
        class CompareOnlyAnalysis(Analysis):
            name: ClassVar[str] = "compare_only"
            Settings: ClassVar[type] = ToySettings
            has_compute_stage: ClassVar[bool] = False
            has_aggregate_stage: ClassVar[bool] = False

        analysis = CompareOnlyAnalysis()
        result = run_analysis(
            analysis,
            toy_condition,
            toy_settings,
            equilibration="0ns",
            output_dir=tmp_path,
        )

        assert result is None
        assert not (tmp_path / "aggregated").exists()


class TestContractEnforcement:
    """Tests for plugin contract violation detection."""

    def test_compute_replicate_none_raises_contract_error(self, tmp_path: Path) -> None:
        """compute_replicate() returning None should raise PluginContractError."""

        class NoneComputeAnalysis(Analysis):
            name: ClassVar[str] = "none_compute"
            Settings: ClassVar[type] = ToySettings

            def compute_replicate(self, ctx, replicate):
                return None

            def aggregate(self, ctx, results):
                return {"ok": True}

        analysis = NoneComputeAnalysis()
        condition = Condition(
            label="Cond",
            config_path=tmp_path / "cfg.yaml",
            replicates=(1,),
            sim_config=SimpleNamespace(),
        )

        with pytest.raises(
            PluginContractError,
            match=r"none_compute.compute_replicate\(\) returned None",
        ):
            run_replicate_once(
                analysis,
                condition,
                ToySettings(),
                equilibration="0ns",
                output_dir=tmp_path / "run_1",
                replicate=1,
                recompute=False,
            )

    def test_contract_error_not_wrapped_as_replicate_error(self, tmp_path: Path) -> None:
        """PluginContractError should propagate, not be wrapped as ReplicateError."""

        class RaisesContractAnalysis(Analysis):
            name: ClassVar[str] = "raises_contract"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 1

            def compute_replicate(self, ctx, replicate):
                raise PluginContractError("contract boom")

            def aggregate(self, ctx, results):
                return {"ok": True}

        condition = Condition(
            label="Cond",
            config_path=tmp_path / "cfg.yaml",
            replicates=(1,),
            sim_config=SimpleNamespace(),
        )

        with pytest.raises(PluginContractError, match="contract boom"):
            _run_replicates(
                RaisesContractAnalysis(),
                condition,
                ToySettings(),
                equilibration="0ns",
                output_dir=tmp_path / "analysis",
            )

    def test_run_comparison_fails_fast_on_contract_error(self, monkeypatch, tmp_path: Path) -> None:
        """Contract violations should abort the comparison, not drop a condition."""

        class NoneComputeAnalysis(Analysis):
            name: ClassVar[str] = "none_compute_comparison"
            Settings: ClassVar[type] = _ParallelSettings
            min_replicates: ClassVar[int] = 1

            def compute_replicate(self, ctx, replicate):
                return None

            def aggregate(self, ctx, results):
                return {"mean_value": 0.0, "sem_value": 0.0, "replicate_values": [0.0]}

        analysis = NoneComputeAnalysis()
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

        with pytest.raises(PluginContractError, match="returned None"):
            run_comparison(analysis, config, recompute=False, equilibration="10ns")

    def test_run_comparison_contract_error_skips_later_conditions(
        self, monkeypatch, tmp_path: Path
    ) -> None:
        """A contract error in one condition should stop processing subsequent conditions."""
        processed_labels: list[str] = []

        class FailsOnFirstCondition(Analysis):
            name: ClassVar[str] = "contract_stop"
            Settings: ClassVar[type] = _ParallelSettings
            min_replicates: ClassVar[int] = 1

            def compute_replicate(self, ctx, replicate):
                del replicate
                processed_labels.append(ctx.condition.label)
                if ctx.condition.label == "A":
                    raise PluginContractError("contract failure on A")
                return {"value": 1.0, "replicate": ctx.replicate}

            def aggregate(self, ctx, results):
                del ctx, results
                return {
                    "mean_value": 1.0,
                    "sem_value": 0.0,
                    "replicate_values": [1.0],
                    "n_replicates": 1,
                }

        analysis = FailsOnFirstCondition()
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

        with pytest.raises(PluginContractError, match="contract failure on A"):
            run_comparison(analysis, config, recompute=False, equilibration="10ns")

        assert processed_labels == ["A"]
