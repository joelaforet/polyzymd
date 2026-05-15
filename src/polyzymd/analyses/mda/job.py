"""Job execution primitives for MDAnalysis-compatible analyses."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Any

from polyzymd.analyses.mda.base import MDAnalysisExtensionError, MDARunKwargs
from polyzymd.analyses.mda.frame_selection import FrameSelection


class MDAAnalysisJobError(MDAnalysisExtensionError):
    """Runtime failure raised when an MDAnalysis job violates its contract."""


@dataclass(frozen=True)
class MDABackendPolicy:
    """Validated MDAnalysis internal backend options for one job.

    The default policy forwards no backend-related keyword arguments so that
    PolyzyMD-level parallelism remains the default. Supplying worker, part, or
    unsupported-backend options requires an explicit backend to avoid ambiguous
    nested-parallelism requests.
    """

    backend: Any = None
    n_workers: int | None = None
    n_parts: int | None = None
    unsupported_backend: bool | None = None
    verbose: bool | None = None
    progressbar_kwargs: Mapping[str, Any] | None = None

    def __post_init__(self) -> None:
        """Validate backend opt-in and worker/part counts.

        Raises
        ------
        ValueError
            Raised when worker/part counts are not positive or backend-specific
            options are supplied without an explicit backend.
        """

        self._validate_positive_count(self.n_workers, field_name="n_workers")
        self._validate_positive_count(self.n_parts, field_name="n_parts")
        if self.backend is None and (
            self.n_workers is not None
            or self.n_parts is not None
            or self.unsupported_backend is not None
        ):
            raise ValueError(
                "n_workers, n_parts, and unsupported_backend require an explicit backend"
            )
        if self.progressbar_kwargs is not None:
            object.__setattr__(self, "progressbar_kwargs", dict(self.progressbar_kwargs))

    def run_kwargs(self) -> MDARunKwargs:
        """Return keyword arguments for ``AnalysisBase.run``.

        Returns
        -------
        MDARunKwargs
            Non-``None`` backend, worker, progress, and verbosity settings.
        """

        kwargs = MDARunKwargs()
        if self.backend is not None:
            kwargs["backend"] = self.backend
        if self.n_workers is not None:
            kwargs["n_workers"] = self.n_workers
        if self.n_parts is not None:
            kwargs["n_parts"] = self.n_parts
        if self.unsupported_backend is not None:
            kwargs["unsupported_backend"] = self.unsupported_backend
        if self.verbose is not None:
            kwargs["verbose"] = self.verbose
        if self.progressbar_kwargs is not None:
            kwargs["progressbar_kwargs"] = dict(self.progressbar_kwargs)
        return kwargs

    def is_default(self) -> bool:
        """Return whether the policy forwards no ``run()`` keyword arguments.

        Returns
        -------
        bool
            ``True`` when no backend, progress, or verbosity options are set.
        """

        return len(self.run_kwargs()) == 0

    @staticmethod
    def _validate_positive_count(value: int | None, *, field_name: str) -> None:
        """Validate an optional positive integer count.

        Parameters
        ----------
        value : int or None
            Candidate count value.
        field_name : str
            Name used in validation errors.

        Raises
        ------
        ValueError
            Raised when ``value`` is not a positive integer.
        """

        if value is None:
            return
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ValueError(f"{field_name} must be a positive integer")


@dataclass(frozen=True)
class MDAUniversePolicy:
    """Lightweight universe provenance and execution policy for one job.

    This policy intentionally does not load universes. It carries only metadata
    and provenance from the layer that already resolved or supplied a universe.
    """

    condition_label: str | None = None
    replicate: int | None = None
    provenance: Any = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Freeze metadata to avoid accidental mutation after job execution."""

        object.__setattr__(self, "metadata", dict(self.metadata))

    def as_dict(self) -> dict[str, Any]:
        """Serialize lightweight policy metadata to primitive values.

        Returns
        -------
        dict[str, Any]
            Dictionary containing condition, replicate, provenance, and metadata.
        """

        provenance = self.provenance
        if hasattr(provenance, "as_dict"):
            provenance = provenance.as_dict()
        return {
            "condition_label": self.condition_label,
            "replicate": self.replicate,
            "provenance": provenance,
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True)
class MDAJobResult:
    """Completed result reference for one MDAnalysis-compatible job."""

    name: str
    analysis: Any
    results: Any
    run_kwargs: Mapping[str, Any]
    frame_selection: FrameSelection
    backend_policy: MDABackendPolicy
    universe_policy: MDAUniversePolicy


class MDAFunctionAdapter:
    """Adapt a simple function to the MDAnalysis ``run()``/``results`` shape.

    The function is called exactly once during ``run()`` as
    ``function(universe, **frame_kwargs, **function_kwargs)``. The adapter does
    not implement trajectory loops; functions that need frame iteration should
    use MDAnalysis primitives internally.
    """

    def __init__(
        self,
        function: Callable[..., Any],
        universe: Any,
        *,
        function_kwargs: Mapping[str, Any] | None = None,
    ) -> None:
        """Initialize the function adapter.

        Parameters
        ----------
        function : Callable[..., Any]
            Function receiving the universe, frame kwargs, and function kwargs.
        universe : Any
            Already-loaded universe or universe-like object supplied by the
            caller.
        function_kwargs : Mapping[str, Any] or None, optional
            Additional keyword arguments forwarded after frame kwargs.

        Raises
        ------
        TypeError
            Raised when ``function`` is not callable.
        """

        if not callable(function):
            raise TypeError("function must be callable")
        self.function = function
        self.universe = universe
        self.function_kwargs = dict(function_kwargs or {})
        self.results: dict[str, Any] | Any = {}

    def run(self, **frame_kwargs: Any) -> MDAFunctionAdapter:
        """Run the adapted function once and store normalized results.

        Parameters
        ----------
        **frame_kwargs : Any
            Frame-selection keyword arguments from ``FrameSelection``.

        Returns
        -------
        MDAFunctionAdapter
            This adapter with ``results`` populated.
        """

        value = self.function(self.universe, **frame_kwargs, **self.function_kwargs)
        self.results = self._normalize_result(value)
        return self

    @staticmethod
    def _normalize_result(value: Any) -> dict[str, Any]:
        """Normalize mapping or scalar function output.

        Parameters
        ----------
        value : Any
            Return value from the adapted function.

        Returns
        -------
        dict[str, Any]
            Mapping output converted to ``dict`` or scalar output under the
            ``"value"`` key.
        """

        if isinstance(value, Mapping):
            return dict(value)
        return {"value": value}


@dataclass
class MDAAnalysisJob:
    """Execute exactly one MDAnalysis ``AnalysisBase``-compatible job.

    A job receives either a ready analysis object or a zero-argument factory
    that constructs one. Execution merges ``FrameSelection.run_kwargs()`` with
    validated backend policy kwargs and calls ``analysis.run(**kwargs)`` once.
    """

    name: str
    frame_selection: FrameSelection = field(default_factory=FrameSelection)
    analysis: Any = None
    analysis_factory: Callable[[], Any] | None = None
    backend_policy: MDABackendPolicy = field(default_factory=MDABackendPolicy)
    universe_policy: MDAUniversePolicy = field(default_factory=MDAUniversePolicy)
    result: MDAJobResult | None = field(default=None, init=False)

    def __post_init__(self) -> None:
        """Validate constructor misuse before execution.

        Raises
        ------
        ValueError
            Raised when both or neither analysis inputs are supplied.
        TypeError
            Raised when supplied inputs do not satisfy the job construction
            contract.
        """

        if not isinstance(self.name, str) or not self.name:
            raise ValueError("name must be a non-empty string")
        self._validate_collaborators()
        has_analysis = self.analysis is not None
        has_factory = self.analysis_factory is not None
        if has_analysis == has_factory:
            raise ValueError("Provide exactly one of analysis or analysis_factory")
        if has_analysis and not callable(getattr(self.analysis, "run", None)):
            raise TypeError("analysis must provide a callable run(**kwargs) method")
        if has_factory and not callable(self.analysis_factory):
            raise TypeError("analysis_factory must be callable")
        if isinstance(self.analysis, MDAFunctionAdapter):
            self._validate_function_backend_policy()

    @classmethod
    def from_function(
        cls,
        name: str,
        function: Callable[..., Any],
        universe: Any,
        *,
        frame_selection: FrameSelection | None = None,
        backend_policy: MDABackendPolicy | None = None,
        universe_policy: MDAUniversePolicy | None = None,
        function_kwargs: Mapping[str, Any] | None = None,
    ) -> MDAAnalysisJob:
        """Create a job from a function and already-loaded universe.

        Parameters
        ----------
        name : str
            Job name used in result metadata and error messages.
        function : Callable[..., Any]
            Function called once by ``MDAFunctionAdapter.run``.
        universe : Any
            Already-loaded universe or universe-like object.
        frame_selection : FrameSelection or None, optional
            Frame selection to forward to the function adapter.
        backend_policy : MDABackendPolicy or None, optional
            Backend policy for this job.
        universe_policy : MDAUniversePolicy or None, optional
            Lightweight universe provenance policy.
        function_kwargs : Mapping[str, Any] or None, optional
            Additional keyword arguments forwarded to the function.

        Returns
        -------
        MDAAnalysisJob
            Job wrapping the function adapter.
        """

        resolved_backend_policy = backend_policy or MDABackendPolicy()
        if not resolved_backend_policy.is_default():
            raise ValueError("Function-adapter jobs do not accept backend policy kwargs")

        return cls(
            name=name,
            analysis=MDAFunctionAdapter(
                function,
                universe,
                function_kwargs=function_kwargs,
            ),
            frame_selection=frame_selection or FrameSelection(),
            backend_policy=resolved_backend_policy,
            universe_policy=universe_policy or MDAUniversePolicy(),
        )

    @property
    def results(self) -> Any:
        """Return completed job results.

        Returns
        -------
        Any
            Results object from the completed analysis.

        Raises
        ------
        MDAAnalysisJobError
            Raised when the job has not run yet.
        """

        if self.result is None:
            raise MDAAnalysisJobError(f"MDAAnalysisJob '{self.name}' has not run yet")
        return self.result.results

    def run(self) -> MDAJobResult:
        """Execute the wrapped analysis and store the completed job result.

        Returns
        -------
        MDAJobResult
            Completed job result reference.

        Raises
        ------
        MDAAnalysisJobError
            Raised when factory construction, analysis execution, or result
            collection violates the runtime job contract.
        """

        analysis = self._resolve_analysis()
        run_kwargs = self._build_run_kwargs(analysis)
        try:
            run_return = analysis.run(**run_kwargs)
        except Exception as exc:
            raise MDAAnalysisJobError(
                f"MDAAnalysisJob '{self.name}' failed during analysis.run(): {exc}"
            ) from exc

        completed_analysis = run_return if hasattr(run_return, "results") else analysis
        if not hasattr(completed_analysis, "results"):
            raise MDAAnalysisJobError(
                f"MDAAnalysisJob '{self.name}' completed without a results attribute"
            )

        self.result = MDAJobResult(
            name=self.name,
            analysis=completed_analysis,
            results=completed_analysis.results,
            run_kwargs=run_kwargs,
            frame_selection=self.frame_selection,
            backend_policy=self.backend_policy,
            universe_policy=self.universe_policy,
        )
        return self.result

    def execute(self) -> MDAJobResult:
        """Execute the job.

        Returns
        -------
        MDAJobResult
            Completed job result reference.
        """

        return self.run()

    def _resolve_analysis(self) -> Any:
        """Return the analysis object for this run.

        Returns
        -------
        Any
            Analysis object with a callable ``run`` method.

        Raises
        ------
        MDAAnalysisJobError
            Raised when the factory fails or returns an invalid object.
        """

        if self.analysis is not None:
            return self.analysis
        try:
            analysis = self.analysis_factory()
        except Exception as exc:
            raise MDAAnalysisJobError(
                f"MDAAnalysisJob '{self.name}' failed to construct analysis: {exc}"
            ) from exc
        if not callable(getattr(analysis, "run", None)):
            raise MDAAnalysisJobError(
                f"MDAAnalysisJob '{self.name}' factory returned an object without run(**kwargs)"
            )
        return analysis

    def _build_run_kwargs(self, analysis: Any) -> dict[str, Any]:
        """Merge frame selection and backend policy keyword arguments.

        Parameters
        ----------
        analysis : Any
            Analysis object selected for this execution.

        Returns
        -------
        dict[str, Any]
            Keyword arguments forwarded to ``analysis.run``.
        """

        if isinstance(analysis, MDAFunctionAdapter):
            self._validate_function_backend_policy()

        return {
            **self.frame_selection.run_kwargs(),
            **self.backend_policy.run_kwargs(),
        }

    def _validate_collaborators(self) -> None:
        """Validate policy collaborators supplied to the job.

        Raises
        ------
        TypeError
            Raised when a collaborator does not implement the expected concrete
            extension-layer policy type.
        """

        if not isinstance(self.frame_selection, FrameSelection):
            raise TypeError("frame_selection must be a FrameSelection instance")
        if not isinstance(self.backend_policy, MDABackendPolicy):
            raise TypeError("backend_policy must be an MDABackendPolicy instance")
        if not isinstance(self.universe_policy, MDAUniversePolicy):
            raise TypeError("universe_policy must be an MDAUniversePolicy instance")

    def _validate_function_backend_policy(self) -> None:
        """Reject backend policy kwargs for function-adapter jobs.

        Raises
        ------
        MDAAnalysisJobError
            Raised when a factory-backed function adapter is paired with a
            non-default backend policy at execution time.
        ValueError
            Raised when a direct function-adapter job is constructed with a
            non-default backend policy.
        """

        if self.backend_policy.is_default():
            return
        message = "Function-adapter jobs do not accept backend policy kwargs"
        if self.analysis is None:
            raise MDAAnalysisJobError(message)
        raise ValueError(message)
