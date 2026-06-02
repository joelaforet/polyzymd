"""Tests for the private one-analysis lifecycle engine."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar, Sequence

import pytest
from pydantic import BaseModel

from polyzymd.analyses._framework.lifecycle import AnalysisLifecycle
from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.exceptions import AggregationError, PluginContractError
from polyzymd.analyses.mda import (
    ArtifactStore,
    MDAAnalysisJob,
    MDABackendPolicy,
    MDAReplicateJobContext,
    MDAUniversePolicy,
)
from polyzymd.analyses.mda.artifacts import ConditionArtifact, ReplicateArtifact
from polyzymd.analyses.mda.plugin import MDACollectorContext
from polyzymd.analyses.mda.store import ArtifactStoreError
from polyzymd.analyses.orchestrator import (
    finalize_comparison_from_disk,
    prepare_comparison_run,
    run_all_plots,
    run_analysis,
    run_comparison,
    run_plot_only,
    run_replicate_once,
)
from polyzymd.config.comparison import PlotSettings


class _LifecycleSettings(BaseModel):
    """Minimal settings model for lifecycle tests."""

    scale: float = 1.0


class _MDAContractMixin:
    """Provide the required MDA lifecycle seam for direct compute fakes."""

    def build_mda_jobs(self, ctx):
        """Return no jobs for tests that override the internal dispatcher."""

        del ctx
        return []


class _OrderAnalysis(_MDAContractMixin, Analysis):
    """Analysis that records compute and aggregate lifecycle calls."""

    name: ClassVar[str] = "lifecycle_order"
    Settings: ClassVar[type] = _LifecycleSettings
    min_replicates: ClassVar[int] = 1

    def __init__(self) -> None:
        self.events: list[str] = []

    def _run_compute_stage(self, ctx: ReplicateContext, replicate: int) -> dict[str, Any]:
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

    def _run_compute_stage(self, ctx: ReplicateContext, replicate: int) -> dict[str, Any]:
        """Raise a plugin contract error without lifecycle wrapping."""

        del ctx, replicate
        raise PluginContractError("contract boom")


class _FakeMDATrajectory:
    """Minimal trajectory exposing a frame count."""

    def __len__(self) -> int:
        """Return the fake frame count."""

        return 8


class _FakeMDAUniverse:
    """Minimal universe-like object for MDA lifecycle tests."""

    trajectory = _FakeMDATrajectory()


class _FakeMDAWindow:
    """Trajectory window compatible with ``FrameSelection.from_trajectory_window``."""

    start = 2
    stop = 8
    step = 2
    equilibration_start = 2
    equilibration_ps = 10.0
    timestep_ps = 5.0
    first_frame_time_ps = None
    selected_start_time_ps = 10.0
    equilibration_time_reference = "loaded_frame_zero"
    n_frames_total = 8
    n_frames_selected = 3
    warning_message = None

    def run_kwargs(self) -> dict[str, int]:
        """Return MDAnalysis run keyword arguments."""

        return {"start": self.start, "stop": self.stop, "step": self.step}


class _FakeMDALoader:
    """Loader seam used by the MDA lifecycle bridge."""

    def __init__(self, sim_config: object) -> None:
        """Store the simulation config."""

        self.sim_config = sim_config

    def load_universe(self, replicate: int) -> _FakeMDAUniverse:
        """Return a fake universe for the requested replicate."""

        assert replicate == 1
        return _FakeMDAUniverse()


class _FakeMDAProvenance:
    """Tiny provenance object exposing the production ``as_dict`` shape."""

    def as_dict(self) -> dict[str, Any]:
        """Return primitive provenance metadata."""

        return {"warnings": ["fake provenance warning"], "loader_class": "_FakeMDALoader"}


class _FakeMDAUniverseProvider:
    """Universe provider that delegates loading to the injected loader."""

    def __init__(self, config: object, *, loader: _FakeMDALoader) -> None:
        """Store provider collaborators."""

        self.config = config
        self.loader = loader

    @classmethod
    def from_config(
        cls,
        config: object,
        *,
        loader: _FakeMDALoader,
    ) -> "_FakeMDAUniverseProvider":
        """Create the fake provider from a simulation config."""

        return cls(config, loader=loader)

    def load_universe(self, replicate: int) -> _FakeMDAUniverse:
        """Load the universe through the fake loader."""

        return self.loader.load_universe(replicate)

    def provenance_for(self, replicate: int) -> _FakeMDAProvenance:
        """Return fake provenance for the requested replicate."""

        assert replicate == 1
        return _FakeMDAProvenance()


class _FakeMDAAnalysisBase:
    """AnalysisBase-like object that records run kwargs."""

    def __init__(self) -> None:
        """Initialize empty results."""

        self.results: dict[str, Any] = {}

    def run(self, **kwargs: Any) -> "_FakeMDAAnalysisBase":
        """Store deterministic results and return self."""

        self.results = {"value": 5.0, "run_kwargs": dict(kwargs)}
        return self


class _FakeMDAnalysisResults(dict):
    """Import-light stand-in for ``MDAnalysis.analysis.results.Results``."""

    __module__ = "MDAnalysis.analysis.results"


class _FakeRawResultsMDAAnalysisBase:
    """AnalysisBase-like object that exposes raw MDAnalysis-style Results."""

    def __init__(self) -> None:
        """Initialize empty results."""

        self.results: _FakeMDAnalysisResults = _FakeMDAnalysisResults()

    def run(self, **kwargs: Any) -> "_FakeRawResultsMDAAnalysisBase":
        """Store deterministic raw results and return self."""

        self.results = _FakeMDAnalysisResults(value=7.0, run_kwargs=dict(kwargs))
        return self


class _RawResultsCollector:
    """Collector that maps fake raw Results to a replicate artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[Any],
    ) -> ReplicateArtifact:
        """Map fake raw results to JSON-safe artifact payloads."""

        job = completed_jobs[0]
        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={
                "jobs": [
                    {
                        "name": job.name,
                        "results": {
                            "value": job.results["value"],
                            "run_kwargs": dict(job.results["run_kwargs"]),
                        },
                    }
                ],
                "n_jobs": 1,
            },
            provenance={"source": "custom_test_collector"},
            metadata={"result_kind": "mapped_raw_results"},
            warnings=list(ctx.warnings),
        )


class _BypassingRawExtraCollector:
    """Collector that bypasses Pydantic validation with a raw extra field."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[Any],
    ) -> ReplicateArtifact:
        """Return an artifact with raw Results stored in an extra field."""

        del completed_jobs
        return ReplicateArtifact.model_construct(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={"jobs": [], "n_jobs": 0},
            sidecars=[],
            provenance={},
            metadata={},
            warnings=[],
            raw=_FakeMDAnalysisResults(value=99.0),
        )


class _MDAJobOnlyAnalysis(Analysis):
    """Analysis that uses only the MDA job lifecycle hook."""

    name: ClassVar[str] = "mda_job_lifecycle"
    Settings: ClassVar[type] = _LifecycleSettings
    min_replicates: ClassVar[int] = 1

    def __init__(self) -> None:
        self.events: list[str] = []

    def _trajectory_loader_factory(self) -> type[_FakeMDALoader]:
        """Return the fake trajectory loader."""

        return _FakeMDALoader

    def _mda_universe_provider_factory(self) -> type[_FakeMDAUniverseProvider]:
        """Return the fake MDA universe provider."""

        return _FakeMDAUniverseProvider

    def get_trajectory_window(self, ctx, replicate, loader, universe) -> _FakeMDAWindow:
        """Return a deterministic frame window."""

        del ctx, replicate, loader, universe
        return _FakeMDAWindow()

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> list[MDAAnalysisJob]:
        """Build one fake MDAnalysis-compatible job."""

        self.events.append(f"build_mda:{ctx.replicate}")
        assert ctx.universe is not None
        assert ctx.frame_selection.run_kwargs() == {"start": 2, "stop": 8, "step": 2}
        assert isinstance(ctx.universe_policy, MDAUniversePolicy)
        return [
            MDAAnalysisJob(
                name="fake_job",
                analysis=_FakeMDAAnalysisBase(),
                frame_selection=ctx.frame_selection,
                backend_policy=ctx.backend_policy,
                universe_policy=ctx.universe_policy,
            )
        ]

    def aggregate(
        self, ctx: AggregateContext, results: Sequence[ReplicateArtifact]
    ) -> dict[str, Any]:
        """Aggregate the saved MDA replicate artifact payload."""

        self.events.append("aggregate")
        assert ctx.replicates == (1,)
        artifact = results[0]
        job = artifact.payload["jobs"][0]
        return {
            "mean_value": job["results"]["value"],
            "run_kwargs": job["results"]["run_kwargs"],
            "warnings": artifact.warnings,
        }


class _RawResultsMDAAnalysis(_MDAJobOnlyAnalysis):
    """MDA analysis whose job exposes raw MDAnalysis-style Results."""

    name: ClassVar[str] = "raw_mda_job_lifecycle"

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> list[MDAAnalysisJob]:
        """Build one fake job with raw Results output."""

        self.events.append(f"build_mda:{ctx.replicate}")
        return [
            MDAAnalysisJob(
                name="raw_job",
                analysis=_FakeRawResultsMDAAnalysisBase(),
                frame_selection=ctx.frame_selection,
                universe_policy=ctx.universe_policy,
            )
        ]


class _CustomCollectorMDAAnalysis(_RawResultsMDAAnalysis):
    """MDA analysis that maps raw Results through a custom collector."""

    name: ClassVar[str] = "custom_collector_mda_lifecycle"

    min_replicates: ClassVar[int] = 1

    def build_mda_collector(self, ctx: MDACollectorContext) -> _RawResultsCollector:
        """Return the custom raw-results collector."""

        assert ctx.analysis_name == self.name
        return _RawResultsCollector()


class _BypassingRawExtraCollectorMDAAnalysis(_MDAJobOnlyAnalysis):
    """MDA analysis whose collector returns raw Results in model extras."""

    name: ClassVar[str] = "bypassing_raw_extra_collector_mda_lifecycle"

    def build_mda_collector(self, ctx: MDACollectorContext) -> _BypassingRawExtraCollector:
        """Return the bypassing collector."""

        assert ctx.analysis_name == self.name
        return _BypassingRawExtraCollector()


class _MetricArtifactAnalysis(Analysis):
    """Analysis that relies on default MDA artifact aggregation."""

    name: ClassVar[str] = "metric_artifact_lifecycle"
    Settings: ClassVar[type] = _LifecycleSettings
    min_replicates: ClassVar[int] = 2

    def __init__(self) -> None:
        self.loader_requests = 0

    def _trajectory_loader_factory(self) -> type[Any]:
        """Return a loader that must not be touched during disk aggregation."""

        self.loader_requests += 1
        raise AssertionError("trajectory loader should not be used for artifact aggregation")

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> list[MDAAnalysisJob]:
        """Satisfy the MDA compute contract without being used in this test."""

        del ctx
        raise AssertionError("MDA jobs should not run during disk aggregation")


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
            _CondCfg(label, tmp_path / f"{label.lower()}.yaml", (1, 2)) for label in labels
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


def test_public_run_analysis_overwrites_stale_aggregate(tmp_path: Path) -> None:
    """Fresh aggregation should always replace stale canonical aggregate JSON."""

    analysis = _OrderAnalysis()
    condition = _condition(tmp_path)
    settings = _LifecycleSettings(scale=2.0)
    output_dir = tmp_path / "analysis" / analysis.name
    aggregate_path = output_dir / "aggregated" / "result.json"
    aggregate_path.parent.mkdir(parents=True, exist_ok=True)
    aggregate_path.write_text(json.dumps({"mean_value": -1.0, "stale": True}))

    result = run_analysis(
        analysis,
        condition,
        settings,
        equilibration="10ns",
        output_dir=output_dir,
        recompute=False,
    )

    saved = json.loads(aggregate_path.read_text())
    assert result["mean_value"] == 3.0
    assert saved["mean_value"] == 3.0
    assert saved["replicate_values"] == [2.0, 4.0]
    assert "stale" not in saved


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


def test_public_lifecycle_runs_mda_jobs_and_saves_artifact(tmp_path: Path) -> None:
    """Public lifecycle should run MDA jobs and save artifacts."""

    analysis = _MDAJobOnlyAnalysis()
    condition = _condition(tmp_path, replicates=(1,))
    output_dir = tmp_path / "analysis" / analysis.name

    result = run_analysis(
        analysis,
        condition,
        _LifecycleSettings(),
        equilibration="10ns",
        output_dir=output_dir,
        recompute=False,
    )

    result_path = output_dir / "run_1" / "result.json"
    saved = ReplicateArtifact.model_validate_json(result_path.read_text())
    assert analysis.events == ["build_mda:1", "aggregate"]
    assert result["mean_value"] == 5.0
    assert result["run_kwargs"] == {"start": 2, "stop": 8, "step": 2}
    assert saved.artifact_type == "replicate"
    assert saved.payload["jobs"][0]["name"] == "fake_job"
    assert saved.payload["jobs"][0]["results"]["value"] == 5.0
    assert "fake provenance warning" in saved.warnings
    assert (output_dir / "aggregated" / "result.json").exists()


def test_public_lifecycle_propagates_backend_policy_to_mda_job(tmp_path: Path) -> None:
    """Replicate contexts should pass configured MDA backend policy to jobs."""

    analysis = _MDAJobOnlyAnalysis()
    condition = _condition(tmp_path, replicates=(1,))
    output_dir = tmp_path / "analysis" / analysis.name

    result = run_analysis(
        analysis,
        condition,
        _LifecycleSettings(),
        equilibration="10ns",
        output_dir=output_dir,
        recompute=False,
        backend_policy=MDABackendPolicy(
            backend="multiprocessing",
            n_workers=2,
            n_parts=4,
        ),
    )

    assert result["run_kwargs"] == {
        "start": 2,
        "stop": 8,
        "step": 2,
        "backend": "multiprocessing",
        "n_workers": 2,
        "n_parts": 4,
    }


def test_public_lifecycle_custom_collector_maps_raw_results(tmp_path: Path) -> None:
    """Custom collectors should map raw MDAnalysis Results to artifacts."""

    analysis = _CustomCollectorMDAAnalysis()
    condition = _condition(tmp_path, replicates=(1,))
    output_dir = tmp_path / "analysis" / analysis.name

    result = run_analysis(
        analysis,
        condition,
        _LifecycleSettings(),
        equilibration="10ns",
        output_dir=output_dir,
        recompute=False,
    )

    saved = ReplicateArtifact.model_validate_json(
        (output_dir / "run_1" / "result.json").read_text()
    )
    assert result["mean_value"] == 7.0
    assert saved.payload["jobs"][0]["name"] == "raw_job"
    assert saved.payload["jobs"][0]["results"]["value"] == 7.0
    assert saved.metadata["result_kind"] == "mapped_raw_results"


def test_public_lifecycle_default_collector_rejects_raw_results(tmp_path: Path) -> None:
    """Default MDA collector should require explicit mapping for raw Results."""

    analysis = _RawResultsMDAAnalysis()
    condition = _condition(tmp_path, replicates=(1,))

    with pytest.raises(PluginContractError, match="raw MDAnalysis Results"):
        run_analysis(
            analysis,
            condition,
            _LifecycleSettings(),
            equilibration="10ns",
            output_dir=tmp_path / "analysis" / analysis.name,
            recompute=False,
        )


def test_public_lifecycle_rejects_collector_raw_extra_before_save(tmp_path: Path) -> None:
    """Lifecycle validation should catch collector artifacts with raw Results extras."""

    analysis = _BypassingRawExtraCollectorMDAAnalysis()
    condition = _condition(tmp_path, replicates=(1,))
    output_dir = tmp_path / "analysis" / analysis.name

    with pytest.raises(PluginContractError, match="raw MDAnalysis Results"):
        run_analysis(
            analysis,
            condition,
            _LifecycleSettings(),
            equilibration="10ns",
            output_dir=output_dir,
            recompute=False,
        )

    assert not (output_dir / "run_1" / "result.json").exists()


def test_aggregate_from_disk_reports_malformed_mda_artifact_context(tmp_path: Path) -> None:
    """MDA artifact loading failures should include analysis, path, and replicate context."""

    analysis = _MDAJobOnlyAnalysis()
    condition = _condition(tmp_path, replicates=(1,))
    output_dir = tmp_path / "analysis" / analysis.name
    result_path = output_dir / "run_1" / "result.json"
    result_path.parent.mkdir(parents=True)
    result_path.write_text('{"artifact_type": "replicate", "analysis_name": "mda_job_lifecycle"}')

    with pytest.raises(ArtifactStoreError) as exc_info:
        AnalysisLifecycle(analysis).aggregate_condition_from_disk(
            condition,
            _LifecycleSettings(),
            "10ns",
            output_dir,
            (1,),
        )

    message = str(exc_info.value)
    assert "mda_job_lifecycle" in message
    assert "condition='Cond'" in message
    assert "replicate=1" in message
    assert str(result_path) in message
    assert "Failed to validate replicate artifact" in message


def test_default_mda_aggregation_from_disk_uses_artifacts_only(tmp_path: Path) -> None:
    """Default MDA aggregation should not load trajectories or universes."""

    analysis = _MetricArtifactAnalysis()
    condition = _condition(tmp_path, replicates=(1, 2))
    settings = _LifecycleSettings()
    output_dir = tmp_path / "analysis" / analysis.name
    settings_fp = analysis.aggregate_settings_fingerprint(settings)
    for replicate, value in ((1, 1.0), (2, 1.4)):
        artifact = ReplicateArtifact(
            analysis_name=analysis.name,
            condition_label=condition.label,
            replicate=replicate,
            payload={"metrics": {"mean_value": value}},
            provenance={"frame_selection": {"start": 0, "stop": 4, "step": 1}},
            metadata={"settings_fingerprint": settings_fp},
        )
        ArtifactStore(output_dir / f"run_{replicate}").write_replicate_result(artifact)

    result = AnalysisLifecycle(analysis).aggregate_condition_from_disk(
        condition,
        settings,
        "10ns",
        output_dir,
        (1, 2),
    )

    assert isinstance(result, ConditionArtifact)
    assert analysis.loader_requests == 0
    assert result.payload["metrics"]["mean_value"]["mean"] == pytest.approx(1.2)
    assert (output_dir / "aggregated" / "result.json").exists()
    loaded = analysis._load_aggregated_result(output_dir / "aggregated")
    assert isinstance(loaded, ConditionArtifact)
    assert loaded.payload == result.payload


def test_default_mda_aggregation_from_disk_records_partial_success(tmp_path: Path) -> None:
    """Default MDA aggregation should preserve skipped replicate provenance."""

    analysis = _MetricArtifactAnalysis()
    condition = _condition(tmp_path, replicates=(1, 2, 3))
    settings = _LifecycleSettings()
    output_dir = tmp_path / "analysis" / analysis.name
    settings_fp = analysis.aggregate_settings_fingerprint(settings)
    for replicate, value in ((1, 1.0), (3, 1.4)):
        artifact = ReplicateArtifact(
            analysis_name=analysis.name,
            condition_label=condition.label,
            replicate=replicate,
            payload={"metrics": {"mean_value": value}},
            provenance={"frame_selection": {"start": 0, "stop": 4, "step": 1}},
            metadata={"settings_fingerprint": settings_fp},
        )
        ArtifactStore(output_dir / f"run_{replicate}").write_replicate_result(artifact)

    result = AnalysisLifecycle(analysis).aggregate_condition_from_disk(
        condition,
        settings,
        "10ns",
        output_dir,
        (1, 2, 3),
    )

    assert isinstance(result, ConditionArtifact)
    assert result.replicates == [1, 3]
    assert result.skipped_replicates == [
        {
            "replicate": 2,
            "reason": "missing artifact",
            "path": str(output_dir / "run_2" / "result.json"),
        }
    ]
    assert result.payload["metrics"]["mean_value"]["mean"] == pytest.approx(1.2)


def test_default_mda_aggregation_rejects_unexpected_disk_artifacts(tmp_path: Path) -> None:
    """Default lifecycle aggregation should use canonical disk discovery."""

    analysis = _MetricArtifactAnalysis()
    condition = _condition(tmp_path, replicates=(1, 2))
    settings = _LifecycleSettings()
    output_dir = tmp_path / "analysis" / analysis.name
    settings_fp = analysis.aggregate_settings_fingerprint(settings)
    for replicate, value in ((1, 1.0), (2, 1.4), (3, 1.8)):
        artifact = ReplicateArtifact(
            analysis_name=analysis.name,
            condition_label=condition.label,
            replicate=replicate,
            payload={"metrics": {"mean_value": value}},
            provenance={"frame_selection": {"start": 0, "stop": 4, "step": 1}},
            metadata={"settings_fingerprint": settings_fp},
        )
        ArtifactStore(output_dir / f"run_{replicate}").write_replicate_result(artifact)

    with pytest.raises(AggregationError, match="unexpected replicate artifact"):
        AnalysisLifecycle(analysis).aggregate_condition_from_disk(
            condition,
            settings,
            "10ns",
            output_dir,
            (1, 2),
        )


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
            raise ValueError("plot boom")

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
    assert failing_failures == [
        (
            "failing_plot_lifecycle",
            "failing_plot_lifecycle: plot failed for comparison='project': ValueError: plot boom",
        )
    ]


def test_public_plot_only_propagates_unexpected_plot_errors(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Plot-only execution should not hide unexpected plugin bugs."""

    class _RuntimePlotAnalysis(_DelegatingAnalysis):
        name: ClassVar[str] = "runtime_plot_lifecycle"

        def plot(self, ctx: PlotContext) -> list[Path]:
            """Raise an unexpected runtime failure from plotting."""

            del ctx
            raise RuntimeError("unexpected plot boom")

    _patch_condition_loader(monkeypatch)
    config = _ComparisonConfig(tmp_path)

    with pytest.raises(RuntimeError, match="unexpected plot boom"):
        run_plot_only(_RuntimePlotAnalysis(), config, equilibration="1ns")


def test_public_run_all_plots_passes_equilibration_override(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Plot-all wrapper should pass equilibration overrides to plot-only execution."""

    captured: dict[str, Any] = {}

    def _get_analysis(name: str):
        captured["analysis_name"] = name
        return _DelegatingAnalysis

    def _run_plot_only(
        analysis: Analysis,
        config: _ComparisonConfig,
        equilibration: str | None = None,
    ) -> tuple[list[Path], list[tuple[str, str]]]:
        captured["analysis"] = analysis.name
        captured["config"] = config
        captured["equilibration"] = equilibration
        return [tmp_path / "plot.png"], []

    monkeypatch.setattr("polyzymd.analyses.discovery.get_analysis", _get_analysis)
    monkeypatch.setattr("polyzymd.analyses.orchestrator.run_plot_only", _run_plot_only)

    paths, failures = run_all_plots(_ComparisonConfig(tmp_path), ["lifecycle_delegate"], "2ns")

    assert paths == [tmp_path / "plot.png"]
    assert failures == []
    assert captured["analysis_name"] == "lifecycle_delegate"
    assert captured["analysis"] == "lifecycle_delegate"
    assert captured["equilibration"] == "2ns"


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
