"""Tests for MDAnalysis-compatible job execution."""

from __future__ import annotations

from typing import Any

import pytest

from polyzymd.analyses.mda import (
    FrameSelection,
    MDAAnalysisJob,
    MDAAnalysisJobError,
    MDABackendPolicy,
    MDAFunctionAdapter,
    MDAJobResult,
    MDAUniversePolicy,
)


class FakeAnalysis:
    """Minimal AnalysisBase-like object used by job tests."""

    def __init__(self, results: Any | None = None) -> None:
        """Initialize fake analysis state."""

        self.results = {} if results is None else results
        self.run_kwargs: dict[str, Any] | None = None
        self.run_count = 0

    def run(self, **kwargs: Any) -> FakeAnalysis:
        """Store run keyword arguments and return this analysis."""

        self.run_kwargs = dict(kwargs)
        self.run_count += 1
        return self


class FakeReturnedAnalysis:
    """Analysis object returned by ``run()`` with its own results."""

    def __init__(self) -> None:
        """Initialize returned analysis results."""

        self.results = {"source": "returned"}


class FakeAnalysisReturningOther(FakeAnalysis):
    """Fake analysis whose ``run()`` returns a result-bearing object."""

    def run(self, **kwargs: Any) -> FakeReturnedAnalysis:
        """Store kwargs and return a separate completed analysis object."""

        self.run_kwargs = dict(kwargs)
        self.run_count += 1
        return FakeReturnedAnalysis()


class FakeAnalysisReturningScalar(FakeAnalysis):
    """Fake analysis whose ``run()`` returns a non-result value."""

    def run(self, **kwargs: Any) -> str:
        """Store kwargs and return a value without ``results``."""

        self.run_kwargs = dict(kwargs)
        self.run_count += 1
        return "ignored"


class FailingAnalysis(FakeAnalysis):
    """Fake analysis that raises during execution."""

    def run(self, **kwargs: Any) -> FailingAnalysis:
        """Raise a runtime error to exercise job wrapping."""

        raise RuntimeError("boom")


class NoResultsAnalysis:
    """Fake analysis with ``run()`` but no ``results`` attribute."""

    def run(self, **kwargs: Any) -> NoResultsAnalysis:
        """Return this invalid completed analysis."""

        return self


def test_job_runs_analysis_with_frame_and_backend_kwargs() -> None:
    """Jobs should merge frame selection and explicit backend policy kwargs."""

    analysis = FakeAnalysis(results={"rmsd": [1.0, 2.0]})
    backend_policy = MDABackendPolicy(
        backend="multiprocessing",
        n_workers=2,
        n_parts=4,
        unsupported_backend=True,
        verbose=False,
        progressbar_kwargs={"desc": "job"},
    )
    universe_policy = MDAUniversePolicy(condition_label="PEG", replicate=1)
    job = MDAAnalysisJob(
        name="fake",
        analysis=analysis,
        frame_selection=FrameSelection(start=2, stop=10, step=2, n_frames_total=12),
        backend_policy=backend_policy,
        universe_policy=universe_policy,
    )

    result = job.run()

    assert isinstance(result, MDAJobResult)
    assert analysis.run_count == 1
    assert analysis.run_kwargs == {
        "start": 2,
        "stop": 10,
        "step": 2,
        "backend": "multiprocessing",
        "n_workers": 2,
        "n_parts": 4,
        "unsupported_backend": True,
        "verbose": False,
        "progressbar_kwargs": {"desc": "job"},
    }
    assert result.name == "fake"
    assert result.analysis is analysis
    assert result.results == {"rmsd": [1.0, 2.0]}
    assert result.run_kwargs == analysis.run_kwargs
    assert result.backend_policy is backend_policy
    assert result.universe_policy is universe_policy
    assert job.result is result
    assert job.results == {"rmsd": [1.0, 2.0]}


def test_default_backend_policy_forwards_no_kwargs() -> None:
    """Default backend policy should keep per-replicate execution serial."""

    analysis = FakeAnalysis()
    job = MDAAnalysisJob(name="serial", analysis=analysis, frame_selection=FrameSelection())

    job.run()

    assert analysis.run_kwargs == {}
    assert job.result is not None
    assert job.result.run_kwargs == {}


def test_job_uses_returned_object_when_it_exposes_results() -> None:
    """MDAnalysis analyses may return a completed result-bearing object."""

    analysis = FakeAnalysisReturningOther(results={"source": "original"})
    job = MDAAnalysisJob(name="returned", analysis=analysis)

    result = job.run()

    assert isinstance(result.analysis, FakeReturnedAnalysis)
    assert result.results == {"source": "returned"}


def test_job_falls_back_to_original_analysis_when_run_return_has_no_results() -> None:
    """Jobs should use the original analysis when ``run()`` returns a scalar."""

    analysis = FakeAnalysisReturningScalar(results={"source": "original"})
    job = MDAAnalysisJob(name="fallback", analysis=analysis)

    result = job.run()

    assert result.analysis is analysis
    assert result.results == {"source": "original"}


def test_job_constructs_analysis_from_factory_for_each_run() -> None:
    """Factories should build a fresh AnalysisBase-like object per execution."""

    analyses: list[FakeAnalysis] = []

    def factory() -> FakeAnalysis:
        """Create a fake analysis and record construction."""

        analysis = FakeAnalysis(results={"value": len(analyses) + 1})
        analyses.append(analysis)
        return analysis

    job = MDAAnalysisJob(name="factory", analysis_factory=factory)

    first = job.execute()
    second = job.run()

    assert len(analyses) == 2
    assert first.analysis is analyses[0]
    assert second.analysis is analyses[1]
    assert first.results == {"value": 1}
    assert second.results == {"value": 2}
    assert job.analysis is None
    assert job.result is second


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"analysis": FakeAnalysis(), "analysis_factory": FakeAnalysis},
    ],
)
def test_job_requires_exactly_one_analysis_source(kwargs: dict[str, Any]) -> None:
    """Constructor misuse should fail before runtime execution."""

    with pytest.raises(ValueError, match="exactly one"):
        MDAAnalysisJob(name="bad", **kwargs)


def test_job_rejects_analysis_without_run() -> None:
    """Analysis objects must provide the MDAnalysis-compatible run method."""

    with pytest.raises(TypeError, match="callable run"):
        MDAAnalysisJob(name="bad", analysis=object())


def test_job_rejects_empty_name() -> None:
    """Job names should be non-empty for useful errors and provenance."""

    with pytest.raises(ValueError, match="non-empty"):
        MDAAnalysisJob(name="", analysis=FakeAnalysis())


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("frame_selection", object(), "FrameSelection"),
        ("backend_policy", object(), "MDABackendPolicy"),
        ("universe_policy", object(), "MDAUniversePolicy"),
    ],
)
def test_job_rejects_invalid_policy_collaborators(
    field: str,
    value: object,
    message: str,
) -> None:
    """Invalid policy collaborators should fail at construction time."""

    with pytest.raises(TypeError, match=message):
        MDAAnalysisJob(name="bad-policy", analysis=FakeAnalysis(), **{field: value})


def test_factory_runtime_contract_failure_is_job_error() -> None:
    """Invalid factory outputs should be reported as job contract failures."""

    job = MDAAnalysisJob(name="factory", analysis_factory=lambda: object())

    with pytest.raises(MDAAnalysisJobError, match="factory returned"):
        job.run()


def test_analysis_runtime_failure_is_job_error() -> None:
    """Analysis execution errors should include job context."""

    job = MDAAnalysisJob(name="failing", analysis=FailingAnalysis())
    with pytest.raises(MDAAnalysisJobError, match="failing.*analysis.run"):
        job.run()


def test_missing_results_after_run_is_job_error() -> None:
    """Completed analyses must expose a results attribute."""

    job = MDAAnalysisJob(name="no-results", analysis=NoResultsAnalysis())
    with pytest.raises(MDAAnalysisJobError, match="without a results"):
        job.run()


def test_results_before_run_is_job_error() -> None:
    """The convenience results property should not mask missing execution."""

    job = MDAAnalysisJob(name="not-run", analysis=FakeAnalysis())
    with pytest.raises(MDAAnalysisJobError, match="has not run"):
        _ = job.results


@pytest.mark.parametrize("field", ["n_workers", "n_parts"])
@pytest.mark.parametrize("value", [0, -1, True, 1.5])
def test_backend_policy_rejects_non_positive_counts(field: str, value: Any) -> None:
    """Worker and part counts must be positive integers."""

    with pytest.raises(ValueError, match=field):
        MDABackendPolicy(backend="multiprocessing", **{field: value})


@pytest.mark.parametrize(
    "kwargs",
    [
        {"n_workers": 2},
        {"n_parts": 2},
        {"unsupported_backend": True},
    ],
)
def test_backend_policy_rejects_backend_ambiguity(kwargs: dict[str, Any]) -> None:
    """Backend-specific options require explicit backend opt-in."""

    with pytest.raises(ValueError, match="explicit backend"):
        MDABackendPolicy(**kwargs)


def test_backend_policy_allows_progress_options_without_backend() -> None:
    """Progress options are run controls and do not imply internal backend use."""

    policy = MDABackendPolicy(verbose=True, progressbar_kwargs={"leave": False})

    assert policy.run_kwargs() == {
        "verbose": True,
        "progressbar_kwargs": {"leave": False},
    }


def test_universe_policy_serializes_provenance_and_metadata() -> None:
    """Universe policy should be provenance-only and JSON-friendly."""

    class Provenance:
        """Small provenance object with an ``as_dict`` method."""

        def as_dict(self) -> dict[str, str]:
            """Return serialized provenance."""

            return {"topology": "top.pdb"}

    policy = MDAUniversePolicy(
        condition_label="control",
        replicate=2,
        provenance=Provenance(),
        metadata={"loader": "fake"},
    )

    assert policy.as_dict() == {
        "condition_label": "control",
        "replicate": 2,
        "provenance": {"topology": "top.pdb"},
        "metadata": {"loader": "fake"},
    }


def test_function_adapter_normalizes_mapping_result() -> None:
    """Function adapters should call functions once with frame and function kwargs."""

    calls: list[dict[str, Any]] = []

    def function(universe: str, **kwargs: Any) -> dict[str, Any]:
        """Return a mapping and record call context."""

        calls.append({"universe": universe, "kwargs": dict(kwargs)})
        return {"mean": 1.5}

    adapter = MDAFunctionAdapter(function, "universe", function_kwargs={"selection": "protein"})

    result = adapter.run(start=1, stop=5)

    assert result is adapter
    assert calls == [
        {
            "universe": "universe",
            "kwargs": {"start": 1, "stop": 5, "selection": "protein"},
        }
    ]
    assert adapter.results == {"mean": 1.5}


def test_function_adapter_normalizes_scalar_result() -> None:
    """Scalar function outputs should be exposed under a stable value key."""

    adapter = MDAFunctionAdapter(lambda universe, **kwargs: 3.25, object())

    adapter.run()
    assert adapter.results == {"value": 3.25}


def test_from_function_creates_executable_job() -> None:
    """Function jobs should use the same job executor contract."""

    def function(universe: str, **kwargs: Any) -> dict[str, Any]:
        """Return universe and frame kwargs for assertion."""

        return {"universe": universe, "kwargs": dict(kwargs)}

    job = MDAAnalysisJob.from_function(
        "function-job",
        function,
        "fake-universe",
        frame_selection=FrameSelection(frames=[0, 2, 4], n_frames_total=5),
        function_kwargs={"metric": "rg"},
    )
    result = job.run()

    assert result.results == {
        "universe": "fake-universe",
        "kwargs": {"frames": (0, 2, 4), "metric": "rg"},
    }
    assert job.results == result.results


def test_function_job_rejects_backend_policy_kwargs() -> None:
    """Function adapters should receive frame kwargs only, not backend kwargs."""

    with pytest.raises(ValueError, match="Function-adapter jobs"):
        MDAAnalysisJob.from_function(
            "function-backend",
            lambda universe, **kwargs: {"ok": True},
            object(),
            backend_policy=MDABackendPolicy(backend="multiprocessing", n_workers=2),
        )


def test_direct_function_adapter_job_rejects_backend_policy_kwargs() -> None:
    """Direct function-adapter jobs should reject non-default backend policy."""

    adapter = MDAFunctionAdapter(lambda universe, **kwargs: {"ok": True}, object())

    with pytest.raises(ValueError, match="Function-adapter jobs"):
        MDAAnalysisJob(
            name="direct-function-backend",
            analysis=adapter,
            backend_policy=MDABackendPolicy(backend="multiprocessing"),
        )


def test_factory_function_adapter_job_rejects_backend_policy_at_runtime() -> None:
    """Factory-backed function adapters should fail before backend kwargs are forwarded."""

    calls: list[dict[str, Any]] = []

    def function(universe: object, **kwargs: Any) -> dict[str, bool]:
        """Record unexpected function calls."""

        calls.append(dict(kwargs))
        return {"ok": True}

    job = MDAAnalysisJob(
        name="factory-function-backend",
        analysis_factory=lambda: MDAFunctionAdapter(function, object()),
        backend_policy=MDABackendPolicy(backend="multiprocessing", n_workers=2),
    )

    with pytest.raises(MDAAnalysisJobError, match="Function-adapter jobs"):
        job.run()
    assert calls == []
