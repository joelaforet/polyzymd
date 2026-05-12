"""Tests for orchestrator worker helper functions."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar, cast

import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import Analysis, Condition, MetricValue
from polyzymd.analyses.exceptions import PluginContractError, ReplicateError
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

    def run_replicate(self, ctx: Any, replicate: int) -> dict[str, Any]:
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


def _worker_aggregate_payload(
    analysis: Analysis,
    settings: BaseModel,
    value: float,
    *,
    replicate: int = 1,
) -> dict[str, Any]:
    """Build a validation-ready worker aggregate payload.

    Parameters
    ----------
    analysis : Analysis
        Analysis instance used to compute the expected settings identity.
    settings : BaseModel
        Active worker settings.
    value : float
        Scalar value represented by the synthetic aggregate.
    replicate : int, optional
        Replicate ID represented by the aggregate, by default 1.

    Returns
    -------
    dict[str, Any]
        Aggregate payload with settings and replicate identity metadata.
    """

    return {
        "mean_value": value,
        "sem_value": 0.0,
        "replicate_values": [value],
        "n_replicates": 1,
        "replicates": [replicate],
        "settings_fingerprint": analysis.aggregate_settings_fingerprint(settings),
    }


class _FailingWorkerAnalysis(_WorkerAnalysis):
    def run_replicate(self, ctx: Any, replicate: int) -> dict[str, Any]:
        raise RuntimeError("test error")


class _RunnerSettings(BaseModel):
    step: int = 2


class _FakeTrajectory:
    def __len__(self) -> int:
        return 10


class _FakeUniverse:
    def __init__(self) -> None:
        self.trajectory = _FakeTrajectory()


class _FakeLoader:
    def __init__(self, sim_config: Any) -> None:
        self.sim_config = sim_config

    def load_universe(self, replicate: int) -> _FakeUniverse:
        self.replicate = replicate
        return _FakeUniverse()

    def get_timestep(self, replicate: int, unit: str = "ps") -> float:
        assert replicate == 1
        assert unit == "ps"
        return 100.0


class _FakeRunner:
    def __init__(self) -> None:
        self.results = SimpleNamespace(values=[])

    def run(self, *, start: int, stop: int, step: int) -> "_FakeRunner":
        self.results.values = list(range(start, stop, step))
        self.results.window = {"start": start, "stop": stop, "step": step}
        return self


class _RunnerWorkerAnalysis(Analysis):
    name: ClassVar[str] = "runner_worker"
    Settings: ClassVar[type] = _RunnerSettings
    min_replicates: ClassVar[int] = 1

    def build_runner(self, ctx: Any, replicate: int, universe: Any, window: Any) -> _FakeRunner:
        assert replicate == 1
        assert len(universe.trajectory) == 10
        assert window.step == ctx.settings.step
        return _FakeRunner()

    def get_trajectory_window(self, ctx: Any, replicate: int, loader: Any, universe: Any) -> Any:
        from polyzymd.analyses.shared.window import resolve_replicate_trajectory_window

        return resolve_replicate_trajectory_window(
            loader=loader,
            replicate=replicate,
            equilibration=ctx.equilibration,
            n_frames_total=len(universe.trajectory),
            step=ctx.settings.step,
        )

    def summarize_replicate(
        self, ctx: Any, replicate: int, runner: Any, window: Any
    ) -> dict[str, Any]:
        return {
            "replicate": replicate,
            "selected_frames": runner.results.values,
            "window": runner.results.window,
            "n_frames_selected": window.n_frames_selected,
        }

    def aggregate(self, ctx: Any, results: list[dict[str, Any]]) -> dict[str, Any]:
        return {
            "n_replicates": len(results),
            "replicate_values": [len(result["selected_frames"]) for result in results],
        }


class _RunnerWithoutRunAnalysis(Analysis):
    name: ClassVar[str] = "runner_without_run"
    Settings: ClassVar[type] = _RunnerSettings
    min_replicates: ClassVar[int] = 1

    def build_runner(self, ctx: Any, replicate: int, universe: Any, window: Any) -> Any:
        del ctx, replicate, universe, window
        return SimpleNamespace()

    def summarize_replicate(
        self, ctx: Any, replicate: int, runner: Any, window: Any
    ) -> dict[str, Any]:
        del ctx, replicate, runner, window
        return {"ok": True}

    def aggregate(self, ctx: Any, results: list[dict[str, Any]]) -> dict[str, Any]:
        del ctx
        return {"n_replicates": len(results)}


class _RunnerWithoutResults:
    def run(self, *, start: int, stop: int, step: int) -> "_RunnerWithoutResults":
        del start, stop, step
        return self


class _RunnerWithoutResultsAnalysis(Analysis):
    name: ClassVar[str] = "runner_without_results"
    Settings: ClassVar[type] = _RunnerSettings
    min_replicates: ClassVar[int] = 1

    def build_runner(self, ctx: Any, replicate: int, universe: Any, window: Any) -> Any:
        del ctx, replicate, universe, window
        return _RunnerWithoutResults()

    def summarize_replicate(
        self, ctx: Any, replicate: int, runner: Any, window: Any
    ) -> dict[str, Any]:
        del ctx, replicate, runner, window
        return {"ok": True}

    def aggregate(self, ctx: Any, results: list[dict[str, Any]]) -> dict[str, Any]:
        del ctx
        return {"n_replicates": len(results)}


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


def test_run_replicate_once_supports_runner_based_compute(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """run_replicate_once should support the runner-based seam."""

    monkeypatch.setattr("polyzymd.analyses.shared.loader.TrajectoryLoader", _FakeLoader)

    analysis = _RunnerWorkerAnalysis()
    condition = Condition("Cond", tmp_path / "cfg.yaml", (1,), cast(Any, SimpleNamespace()))
    settings = _RunnerSettings(step=2)
    run_dir = tmp_path / "analysis" / "cond" / "runner_worker" / "run_1"

    result = run_replicate_once(
        analysis,
        condition,
        settings,
        "200ps",
        run_dir,
        replicate=1,
        recompute=False,
    )

    assert result["replicate"] == 1
    assert result["selected_frames"] == [2, 4, 6, 8]
    assert result["window"] == {"start": 2, "stop": 10, "step": 2}
    assert result["n_frames_selected"] == 4
    assert (run_dir / "result.json").exists()


def test_run_replicate_once_rejects_runner_without_callable_run(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """Runner seam should reject objects without a callable run()."""

    monkeypatch.setattr("polyzymd.analyses.shared.loader.TrajectoryLoader", _FakeLoader)

    analysis = _RunnerWithoutRunAnalysis()
    condition = Condition("Cond", tmp_path / "cfg.yaml", (1,), cast(Any, SimpleNamespace()))
    settings = _RunnerSettings(step=2)
    run_dir = tmp_path / "analysis" / "cond" / "runner_without_run" / "run_1"

    with pytest.raises(PluginContractError, match="must return an object with callable run"):
        run_replicate_once(
            analysis,
            condition,
            settings,
            "200ps",
            run_dir,
            replicate=1,
            recompute=False,
        )


def test_run_replicate_once_rejects_runner_without_results(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """Runner seam should reject executed runners without results."""

    monkeypatch.setattr("polyzymd.analyses.shared.loader.TrajectoryLoader", _FakeLoader)

    analysis = _RunnerWithoutResultsAnalysis()
    condition = Condition("Cond", tmp_path / "cfg.yaml", (1,), cast(Any, SimpleNamespace()))
    settings = _RunnerSettings(step=2)
    run_dir = tmp_path / "analysis" / "cond" / "runner_without_results" / "run_1"

    with pytest.raises(PluginContractError, match="must expose results after run"):
        run_replicate_once(
            analysis,
            condition,
            settings,
            "200ps",
            run_dir,
            replicate=1,
            recompute=False,
        )


def test_run_analysis_raises_structured_replicate_error_on_worker_exception(tmp_path: Path) -> None:
    """run_analysis should raise ReplicateError on unexpected worker failures."""
    analysis = _FailingWorkerAnalysis()
    condition = Condition("Cond", tmp_path / "cfg.yaml", (1,), cast(Any, SimpleNamespace()))
    settings = _WorkerSettings(scale=1.0)

    with pytest.raises(ReplicateError, match="run_replicate failed"):
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
    settings = _WorkerSettings()
    aggregated = {
        "A": _worker_aggregate_payload(analysis, settings, 1.0),
        "B": _worker_aggregate_payload(analysis, settings, 2.0),
    }
    out = finalize_comparison_from_disk(
        analysis=analysis,
        config=cast(Any, config),
        analysis_dirs=analysis_dirs,
        aggregated_results=aggregated,
        results_dir=tmp_path / "comparison" / "worker_toy",
        figures_dir=tmp_path / "figures" / "worker_toy",
        settings=settings,
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
    settings = _WorkerSettings()
    aggregated = {
        "B": _worker_aggregate_payload(analysis, settings, 2.0),
        "C": _worker_aggregate_payload(analysis, settings, 3.0),
    }

    caplog.set_level("WARNING")
    out = finalize_comparison_from_disk(
        analysis=analysis,
        config=cast(Any, config),
        analysis_dirs=analysis_dirs,
        aggregated_results=aggregated,
        results_dir=tmp_path / "comparison" / "worker_toy",
        figures_dir=tmp_path / "figures" / "worker_toy",
        settings=settings,
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
    settings = _WorkerSettings()

    out = finalize_comparison_from_disk(
        analysis=analysis,
        config=cast(Any, config),
        analysis_dirs={"B": tmp_path / "analysis" / "B" / "worker_toy"},
        aggregated_results={"B": _worker_aggregate_payload(analysis, settings, 2.0)},
        results_dir=tmp_path / "comparison" / "worker_toy",
        figures_dir=tmp_path / "figures" / "worker_toy",
        settings=settings,
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


def test_finalize_filtered_control_proceeds_without_allow_partial(
    monkeypatch,
    caplog,
    tmp_path: Path,
) -> None:
    """Filtered controls should auto-switch to all-vs-all without allow_partial."""
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

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.Condition.from_condition_config",
        lambda cond_cfg: cond_a if cond_cfg.label == "A" else cond_b,
    )

    prepared_state = {
        "all_conditions": [cond_a, cond_b],
        "valid_conditions": [cond_b],
        "excluded_conditions": [cond_a],
        "condition_by_label": {"A": cond_a, "B": cond_b},
        "settings": _WorkerSettings(),
        "equilibration": "10ns",
        "analysis_root": tmp_path / "analysis",
    }
    settings = _WorkerSettings()

    caplog.set_level("WARNING")
    out = finalize_comparison_from_disk(
        analysis=analysis,
        config=cast(Any, config),
        analysis_dirs={"B": tmp_path / "analysis" / "B" / "worker_toy"},
        aggregated_results={"B": _worker_aggregate_payload(analysis, settings, 2.0)},
        results_dir=tmp_path / "comparison" / "worker_toy",
        figures_dir=tmp_path / "figures" / "worker_toy",
        settings=settings,
        effective_control="A",
        prepared_state=prepared_state,
        allow_partial=False,
    )

    assert out["comparison"] is not None
    assert out["comparison"].control_label is None
    assert "was excluded by worker_toy.filter_conditions()" in caplog.text


def test_finalize_missing_non_filtered_control_still_requires_allow_partial(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """Missing controls from runtime failures should still fail in strict mode."""
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

    monkeypatch.setattr(
        "polyzymd.analyses.orchestrator.Condition.from_condition_config",
        lambda cond_cfg: cond_a if cond_cfg.label == "A" else cond_b,
    )

    prepared_state = {
        "all_conditions": [cond_a, cond_b],
        "valid_conditions": [cond_a, cond_b],
        "excluded_conditions": [],
        "condition_by_label": {"A": cond_a, "B": cond_b},
        "settings": _WorkerSettings(),
        "equilibration": "10ns",
        "analysis_root": tmp_path / "analysis",
    }

    with pytest.raises(ValueError, match="missing aggregated results"):
        finalize_comparison_from_disk(
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
            prepared_state=prepared_state,
            allow_partial=False,
        )


# === Typed replicate result round-trip tests ===


class _TypedReplicateResult(BaseModel):
    """Pydantic replicate result for typed round-trip tests."""

    value: float
    replicate: int


class _TypedAggregatedResult(BaseModel):
    """Pydantic aggregated result for typed round-trip tests."""

    mean_value: float
    replicate_values: list[float]
    n_replicates: int
    replicates: list[int] | None = None
    settings_fingerprint: str | None = None


class _TypedWorkerAnalysis(Analysis):
    """Analysis plugin that uses Pydantic models for replicate results."""

    name: ClassVar[str] = "typed_worker_toy"
    Settings: ClassVar[type] = _WorkerSettings
    ReplicateResultClass: ClassVar[type | None] = _TypedReplicateResult
    AggregatedResultClass: ClassVar[type | None] = _TypedAggregatedResult
    min_replicates: ClassVar[int] = 1

    def run_replicate(self, ctx: Any, replicate: int) -> _TypedReplicateResult:
        return _TypedReplicateResult(
            value=float(replicate) * float(ctx.settings.scale),
            replicate=replicate,
        )

    def aggregate(self, ctx, results) -> _TypedAggregatedResult:
        # This accesses .value and would fail with dicts
        values = [result.value for result in results]
        return _TypedAggregatedResult(
            mean_value=sum(values) / len(values),
            replicate_values=values,
            n_replicates=len(values),
        )

    def extract_metrics(self, summary):
        return {}


class TestTypedReplicateRoundTrip:
    """Verify Pydantic model round-trip through disk serialization.

    This is the regression test for the HPC aggregate deserialization bug
    where replicate results loaded as dicts instead of Pydantic models.
    """

    def test_aggregate_from_disk_returns_typed_results(self, tmp_path: Path) -> None:
        """Replicate results loaded from disk should be Pydantic models, not dicts."""
        analysis = _TypedWorkerAnalysis()
        condition = Condition("Typed", tmp_path / "cfg.yaml", (1, 2), cast(Any, SimpleNamespace()))
        settings = _WorkerSettings(scale=3.0)
        cond_dir = tmp_path / "analysis" / "cond" / "typed_worker_toy"

        # Phase 1: compute and save replicate results to disk
        run_replicate_once(analysis, condition, settings, "10ns", cond_dir / "run_1", 1, False)
        run_replicate_once(analysis, condition, settings, "10ns", cond_dir / "run_2", 2, False)

        # Verify JSON files exist
        assert (cond_dir / "run_1" / "result.json").exists()
        assert (cond_dir / "run_2" / "result.json").exists()

        # Phase 2: aggregate from disk using the worker aggregate path
        aggregated = aggregate_condition_from_disk(
            analysis,
            condition,
            settings,
            "10ns",
            cond_dir,
            replicates=(1, 2),
        )

        # Aggregate should succeed with typed replicate models
        assert isinstance(aggregated, _TypedAggregatedResult)
        assert aggregated.n_replicates == 2
        assert aggregated.mean_value == pytest.approx(4.5)
        assert aggregated.replicate_values == [3.0, 6.0]

    def test_deserialize_replicate_result_returns_model(self, tmp_path: Path) -> None:
        """Deserialize should return a model when ReplicateResultClass is set."""
        analysis = _TypedWorkerAnalysis()
        result = _TypedReplicateResult(value=42.0, replicate=1)

        result_path = tmp_path / "result.json"
        result_path.write_text(result.model_dump_json())

        loaded = analysis._deserialize_replicate_result(result_path)
        assert isinstance(loaded, _TypedReplicateResult)
        assert loaded.value == 42.0
        assert loaded.replicate == 1

    def test_dict_analysis_still_works(self, tmp_path: Path) -> None:
        """Dict-based analyses with no ReplicateResultClass should still work."""
        analysis = _WorkerAnalysis()
        condition = Condition("Dict", tmp_path / "cfg.yaml", (1,), cast(Any, SimpleNamespace()))
        settings = _WorkerSettings(scale=2.0)
        cond_dir = tmp_path / "analysis" / "cond" / "worker_toy"

        run_replicate_once(analysis, condition, settings, "10ns", cond_dir / "run_1", 1, False)

        aggregated = aggregate_condition_from_disk(
            analysis,
            condition,
            settings,
            "10ns",
            cond_dir,
            replicates=(1,),
        )
        assert isinstance(aggregated, dict)
        assert aggregated["n_replicates"] == 1
