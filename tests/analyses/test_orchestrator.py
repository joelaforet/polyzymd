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
    DependencyError,
    PluginContractError,
    ReplicateError,
    ReplicateSkippedError,
)
from polyzymd.analyses.orchestrator import (
    _validate_dependencies,
    aggregate_condition_from_disk,
    finalize_comparison_from_disk,
    order_analyses_for_execution,
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

    def run_replicate(self, ctx, replicate: int) -> dict[str, Any]:
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

    def run_replicate(self, ctx: ReplicateContext, replicate: int) -> ToyResult:
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

    def run_replicate(self, ctx: ReplicateContext, replicate: int) -> ToyResult:
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

    prepared = prepare_comparison_run(
        analysis,
        config,
        "5ns",
    )
    conditions = prepared["valid_conditions"]
    settings = prepared["settings"]
    equilibration = prepared["equilibration"]
    analysis_root = prepared["analysis_root"]
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


def test_aggregate_from_disk_recompute_removes_stale_aggregate_dir(tmp_path: Path) -> None:
    """aggregate_condition_from_disk should clean aggregate sidecars on recompute."""

    class _AggregateCleanupAnalysis(_ParallelAnalysis):
        name: ClassVar[str] = "aggregate_cleanup_toy"

        def aggregate(self, ctx, results) -> dict[str, Any]:
            assert ctx.recompute is True
            assert not (ctx.output_dir / "stale_sidecar.txt").exists()
            return super().aggregate(ctx, results)

    analysis = _AggregateCleanupAnalysis()
    condition = Condition("A", tmp_path / "a.yaml", (1, 2), SimpleNamespace())
    settings = _ParallelSettings(factor=1.0)
    output_dir = tmp_path / "analysis" / "A" / analysis.name

    run_replicate_once(analysis, condition, settings, "10ns", output_dir / "run_1", 1, False)
    run_replicate_once(analysis, condition, settings, "10ns", output_dir / "run_2", 2, False)
    stale_sidecar = output_dir / "aggregated" / "stale_sidecar.txt"
    stale_sidecar.parent.mkdir(parents=True, exist_ok=True)
    stale_sidecar.write_text("stale")

    aggregate_condition_from_disk(
        analysis,
        condition,
        settings,
        "10ns",
        output_dir,
        replicates=(1, 2),
        recompute=True,
    )

    assert not stale_sidecar.exists()


def test_aggregate_from_disk_missing_replicates_reports_expected_paths(tmp_path: Path) -> None:
    """Missing replicate outputs should report expected run_N/result.json paths."""

    class _NeedsTwoReplicates(_ParallelAnalysis):
        name: ClassVar[str] = "needs_two_replicates"
        min_replicates: ClassVar[int] = 2

    analysis = _NeedsTwoReplicates()
    condition = Condition("A", tmp_path / "a.yaml", (1, 2), SimpleNamespace())
    settings = _ParallelSettings(factor=1.0)
    output_dir = tmp_path / "analysis" / "A" / analysis.name

    run_replicate_once(analysis, condition, settings, "10ns", output_dir / "run_1", 1, False)

    with pytest.raises(ValueError, match=r"replicate result\(s\) on disk") as exc_info:
        aggregate_condition_from_disk(
            analysis,
            condition,
            settings,
            "10ns",
            output_dir,
            replicates=(1, 2),
        )

    message = str(exc_info.value)
    assert str(output_dir / "run_2" / "result.json") in message
    assert "Expected missing replicate output path" in message


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


def test_run_comparison_recompute_contexts_and_cleanup(monkeypatch, tmp_path: Path) -> None:
    """run_comparison should propagate recompute and clean stale owned outputs."""

    class _RecomputeCleanupAnalysis(_ParallelAnalysis):
        name: ClassVar[str] = "recompute_cleanup_toy"

        def __init__(self) -> None:
            self.aggregate_recompute: list[bool] = []
            self.compare_recompute: bool | None = None
            self.plot_recompute: bool | None = None

        def aggregate(self, ctx, results) -> dict[str, Any]:
            self.aggregate_recompute.append(ctx.recompute)
            assert not (ctx.output_dir / "stale_sidecar.txt").exists()
            return super().aggregate(ctx, results)

        def compare(self, ctx) -> dict[str, Any]:
            self.compare_recompute = ctx.recompute
            assert ctx.result_path is not None
            assert not ctx.result_path.exists()
            return {"ok": True}

        def plot(self, ctx) -> list[Path]:
            self.plot_recompute = ctx.recompute
            assert not (ctx.output_dir / "stale_plot.png").exists()
            assert (ctx.output_dir.parent / "unrelated" / "keep.txt").exists()
            out = ctx.output_dir / "new_plot.png"
            out.write_text("new")
            return [out]

    analysis = _RecomputeCleanupAnalysis()
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

    for label in ("A", "B"):
        stale_sidecar = (
            tmp_path / "analysis" / label / analysis.name / "aggregated" / "stale_sidecar.txt"
        )
        stale_sidecar.parent.mkdir(parents=True, exist_ok=True)
        stale_sidecar.write_text("stale")
    stale_result = tmp_path / "comparison" / analysis.name / "result.json"
    stale_result.parent.mkdir(parents=True, exist_ok=True)
    stale_result.write_text('{"stale": true}')
    stale_plot = tmp_path / "figures" / analysis.name / "stale_plot.png"
    stale_plot.parent.mkdir(parents=True, exist_ok=True)
    stale_plot.write_text("stale")
    unrelated = tmp_path / "figures" / "unrelated" / "keep.txt"
    unrelated.parent.mkdir(parents=True, exist_ok=True)
    unrelated.write_text("keep")

    result = run_comparison(analysis, config, recompute=True, equilibration="10ns")

    assert result["comparison"] == {"ok": True}
    assert analysis.aggregate_recompute == [True, True]
    assert analysis.compare_recompute is True
    assert analysis.plot_recompute is True
    assert not stale_plot.exists()
    assert unrelated.exists()


def test_finalize_missing_aggregates_reports_expected_paths(tmp_path: Path) -> None:
    """Finalize errors should include aggregate paths and partial-finalize hints."""
    analysis = _ParallelAnalysis()
    condition_a = Condition("A", tmp_path / "a.yaml", (1,), SimpleNamespace())
    condition_b = Condition("B", tmp_path / "b.yaml", (1,), SimpleNamespace())
    config = _make_config(tmp_path)
    analysis_root = tmp_path / "analysis"
    analysis_dirs = {"A": analysis_root / "A" / analysis.name}
    prepared_state = {
        "all_conditions": [condition_a, condition_b],
        "valid_conditions": [condition_a, condition_b],
        "excluded_conditions": [],
        "condition_by_label": {"A": condition_a, "B": condition_b},
        "settings": _ParallelSettings(),
        "equilibration": "10ns",
        "analysis_root": analysis_root,
    }

    with pytest.raises(ValueError, match="missing aggregated results") as exc_info:
        finalize_comparison_from_disk(
            analysis=analysis,
            config=config,
            analysis_dirs=analysis_dirs,
            aggregated_results={"A": {"ok": True}},
            results_dir=tmp_path / "comparison" / analysis.name,
            figures_dir=tmp_path / "figures" / analysis.name,
            settings=_ParallelSettings(),
            effective_control=None,
            prepared_state=prepared_state,
            allow_partial=False,
        )

    message = str(exc_info.value)
    assert str(analysis_root / "B" / analysis.name / "aggregated" / "result.json") in message
    assert "--allow-partial (CLI) / allow_partial=True (API)" in message


def test_finalize_all_conditions_dropped_mentions_aggregate_jobs(tmp_path: Path) -> None:
    """Dropping every condition should explain that aggregate files were absent."""
    analysis = _ParallelAnalysis()
    condition_a = Condition("A", tmp_path / "a.yaml", (1,), SimpleNamespace())
    config = _make_config(tmp_path)
    prepared_state = {
        "all_conditions": [condition_a],
        "valid_conditions": [condition_a],
        "excluded_conditions": [],
        "condition_by_label": {"A": condition_a},
        "settings": _ParallelSettings(),
        "equilibration": "10ns",
        "analysis_root": tmp_path / "analysis",
    }

    with pytest.raises(ValueError, match="no aggregate files were found") as exc_info:
        finalize_comparison_from_disk(
            analysis=analysis,
            config=config,
            analysis_dirs={},
            aggregated_results={},
            results_dir=tmp_path / "comparison" / analysis.name,
            figures_dir=tmp_path / "figures" / analysis.name,
            settings=_ParallelSettings(),
            effective_control=None,
            prepared_state=prepared_state,
            allow_partial=True,
        )

    assert "Re-run aggregate jobs" in str(exc_info.value)


class TestOrchestrator:
    """Test the orchestrator's replicate running and dependency sorting."""

    def test_run_replicate_once_uses_canonical_entry_point(
        self,
        toy_condition,
        toy_settings,
        tmp_path,
    ):
        """run_replicate_once should call the canonical entry point."""

        class CanonicalRunAnalysis(Analysis):
            """Analysis that tracks canonical entry point calls."""

            name: ClassVar[str] = "canonical_run"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 1

            def __init__(self) -> None:
                """Initialize call tracking for the canonical hook."""

                self.called_replicates: list[int] = []

            def run_replicate(self, ctx: ReplicateContext, replicate: int) -> ToyResult:
                """Return a result through the canonical hook."""

                del ctx
                self.called_replicates.append(replicate)
                return ToyResult(value=float(replicate), replicate=replicate)

            def aggregate(self, ctx, results):
                """Return a simple aggregate result."""

                del ctx, results
                return {"ok": True}

        analysis = CanonicalRunAnalysis()
        result = run_replicate_once(
            analysis,
            toy_condition,
            toy_settings,
            "10ns",
            tmp_path / "run_1",
            1,
            False,
        )

        assert result == ToyResult(value=1.0, replicate=1)
        assert analysis.called_replicates == [1]
        assert (tmp_path / "run_1" / "result.json").exists()

    def test_run_analysis_success(self, toy_analysis, toy_condition, toy_settings, tmp_path):
        result = run_analysis(
            toy_analysis,
            toy_condition,
            toy_settings,
            equilibration="10ns",
            output_dir=tmp_path,
        )
        assert isinstance(result, ToyAggregatedResult)
        assert result.n_replicates == 3
        assert (tmp_path / "run_1").exists()
        assert (tmp_path / "run_2").exists()
        assert (tmp_path / "run_3").exists()

    def test_run_analysis_partial_failure(self, toy_condition, toy_settings, tmp_path):
        """If some replicates fail but min_replicates is met, continue."""

        class FailingAnalysis(Analysis):
            name: ClassVar[str] = "failing"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 1

            def run_replicate(self, ctx, replicate):
                if replicate == 2:
                    raise FileNotFoundError("Trajectory not found")
                return ToyResult(value=float(replicate), replicate=replicate)

            def aggregate(self, ctx, results):
                return {
                    "mean_value": float(sum(result.value for result in results) / len(results)),
                    "sem_value": 0.0,
                    "replicate_values": [float(result.value) for result in results],
                    "n_replicates": len(results),
                }

        analysis = FailingAnalysis()
        result = run_analysis(
            analysis,
            toy_condition,
            toy_settings,
            equilibration="0ns",
            output_dir=tmp_path,
        )
        assert result["n_replicates"] == 2
        assert sorted(result["replicate_values"]) == [1.0, 3.0]

    def test_run_analysis_skipped_error_is_recoverable(
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

            def run_replicate(self, ctx, replicate):
                if replicate == 2:
                    raise ReplicateSkippedError("No trajectory data found for replicate 2")
                return ToyResult(value=float(replicate), replicate=replicate)

            def aggregate(self, ctx, results):
                return {
                    "mean_value": float(sum(result.value for result in results) / len(results)),
                    "sem_value": 0.0,
                    "replicate_values": [float(result.value) for result in results],
                    "n_replicates": len(results),
                }

        caplog.set_level("WARNING")
        result = run_analysis(
            SkipReplicateAnalysis(),
            toy_condition,
            toy_settings,
            equilibration="0ns",
            output_dir=tmp_path,
        )
        assert result["n_replicates"] == 2
        assert sorted(result["replicate_values"]) == [1.0, 3.0]
        assert "No trajectory data found for replicate 2" in caplog.text

    def test_run_analysis_below_minimum_raises(self, toy_condition, toy_settings, tmp_path):
        """If fewer than min_replicates succeed, raise ValueError."""

        class AlwaysFailAnalysis(Analysis):
            name: ClassVar[str] = "always_fail"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 2

            def run_replicate(self, ctx, replicate):
                raise FileNotFoundError("Missing trajectory")

            def aggregate(self, ctx, results):
                return None

        with pytest.raises(ValueError, match="need at least 2"):
            run_analysis(
                AlwaysFailAnalysis(),
                toy_condition,
                toy_settings,
                equilibration="0ns",
                output_dir=tmp_path,
            )

    def test_run_analysis_all_skipped_fails_minimum(self, toy_condition, toy_settings, tmp_path):
        """When all replicates skip, minimum replicate validation should fail."""

        class AllSkippedAnalysis(Analysis):
            name: ClassVar[str] = "all_skipped"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 1

            def run_replicate(self, ctx, replicate):
                raise ReplicateSkippedError(f"No trajectory data found for replicate {replicate}")

            def aggregate(self, ctx, results):
                return {"ok": True}

        with pytest.raises(ValueError, match="need at least 1"):
            run_analysis(
                AllSkippedAnalysis(),
                toy_condition,
                toy_settings,
                equilibration="0ns",
                output_dir=tmp_path,
            )

    def test_run_analysis_unexpected_failure_raises_replicate_error(
        self, toy_condition, toy_settings, tmp_path
    ):
        """Unexpected compute failures should raise structured ReplicateError."""

        class ExplodingAnalysis(Analysis):
            name: ClassVar[str] = "exploding"
            Settings: ClassVar[type] = ToySettings

            def run_replicate(self, ctx, replicate):
                raise RuntimeError("boom")

            def aggregate(self, ctx, results):
                return None

        with pytest.raises(ReplicateError, match="condition='Test Condition' replicate=1"):
            run_analysis(
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

            def run_replicate(self, ctx, replicate):
                return ["not", "valid"]

            def aggregate(self, ctx, results):
                return {"ok": True}

        with pytest.raises(PluginContractError, match="invalid_compute.run_replicate"):
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

            def run_replicate(self, ctx, replicate):
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

            def run_replicate(self, ctx, replicate):
                pass

            def aggregate(self, ctx, results):
                pass

        class CircB(Analysis):
            name: ClassVar[str] = "circ_b"
            Settings: ClassVar[type] = ToySettings
            dependencies: ClassVar[tuple[str, ...]] = ("circ_a",)

            def run_replicate(self, ctx, replicate):
                pass

            def aggregate(self, ctx, results):
                pass

        with pytest.raises(ValueError, match="Circular dependency"):
            _topological_sort([CircA(), CircB()])

    def test_order_analyses_for_execution_returns_dependency_order(self, monkeypatch) -> None:
        """Public ordering helper should return canonical dependency order."""

        class _A(Analysis):
            name: ClassVar[str] = "a"
            aliases: ClassVar[tuple[str, ...]] = ("alias_a",)
            Settings: ClassVar[type] = ToySettings
            dependencies: ClassVar[tuple[str, ...]] = ()

            def run_replicate(self, ctx, replicate):
                return {"replicate": replicate}

            def aggregate(self, ctx, results):
                return {"n": len(results)}

        class _B(Analysis):
            name: ClassVar[str] = "b"
            Settings: ClassVar[type] = ToySettings
            dependencies: ClassVar[tuple[str, ...]] = ("a",)

            def run_replicate(self, ctx, replicate):
                return {"replicate": replicate}

            def aggregate(self, ctx, results):
                return {"n": len(results)}

        monkeypatch.setattr(
            "polyzymd.analyses.discovery.get_analysis",
            lambda name: _A if name in {"a", "alias_a"} else _B,
        )
        monkeypatch.setattr(
            "polyzymd.analyses.discovery.list_all_names",
            lambda: ["a", "alias_a", "b"],
        )

        ordered = order_analyses_for_execution(["b", "alias_a"])
        assert ordered == ["a", "b"]

    def test_order_analyses_allows_satisfied_external_dependencies(self, monkeypatch) -> None:
        """Ordering should allow dependencies satisfied outside the run list."""

        class _Contacts(Analysis):
            name: ClassVar[str] = "contacts"
            Settings: ClassVar[type] = ToySettings
            dependencies: ClassVar[tuple[str, ...]] = ()

            def run_replicate(self, ctx, replicate):
                return {"replicate": replicate}

            def aggregate(self, ctx, results):
                return {"n": len(results)}

        class _DependentAnalysis(Analysis):
            name: ClassVar[str] = "dependent_analysis"
            Settings: ClassVar[type] = ToySettings
            dependencies: ClassVar[tuple[str, ...]] = ("contacts",)

            def run_replicate(self, ctx, replicate):
                return {"replicate": replicate}

            def aggregate(self, ctx, results):
                return {"n": len(results)}

        monkeypatch.setattr(
            "polyzymd.analyses.discovery.get_analysis",
            lambda name: _DependentAnalysis if name == "dependent_analysis" else _Contacts,
        )
        monkeypatch.setattr(
            "polyzymd.analyses.discovery.list_all_names",
            lambda: ["contacts", "dependent_analysis"],
        )

        with pytest.raises(DependencyError, match="not in the current run list"):
            order_analyses_for_execution(["dependent_analysis"])

        ordered = order_analyses_for_execution(["dependent_analysis"], satisfied={"contacts"})
        assert ordered == ["dependent_analysis"]

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

    def test_run_replicate_none_raises_contract_error(self, tmp_path: Path) -> None:
        """run_replicate() returning None should raise PluginContractError."""

        class NoneComputeAnalysis(Analysis):
            name: ClassVar[str] = "none_compute"
            Settings: ClassVar[type] = ToySettings

            def run_replicate(self, ctx, replicate):
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
            match=r"none_compute.run_replicate\(\) returned None",
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

            def run_replicate(self, ctx, replicate):
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
            run_analysis(
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

            def run_replicate(self, ctx, replicate):
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

            def run_replicate(self, ctx, replicate):
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

    def test_validate_dependencies_raises_when_dependency_missing(self, monkeypatch) -> None:
        """Declared dependencies missing from scheduled analyses should fail."""
        monkeypatch.setattr(
            "polyzymd.analyses.discovery.list_all_names",
            lambda: ["toy", "toy_dependent"],
        )
        with pytest.raises(DependencyError, match="not in the current run list"):
            _validate_dependencies([ToyDependentAnalysis()])
