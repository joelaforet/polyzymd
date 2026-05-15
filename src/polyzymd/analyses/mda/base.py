"""Import-light primitives for the MDAnalysis extension layer."""

from __future__ import annotations

from typing import Any, Final, Protocol, TypedDict, runtime_checkable

MDA_EXTENSION_API_VERSION: Final[str] = "1"


class MDARunKwargs(TypedDict, total=False):
    """Keyword arguments accepted by ``MDAnalysis.analysis.base.AnalysisBase.run``.

    The value types stay intentionally broad where MDAnalysis accepts multiple
    backend or frame-selector forms. This keeps the public extension-layer type
    import-light and independent of optional MDAnalysis runtime objects.
    """

    start: int | None
    stop: int | None
    step: int | None
    frames: Any
    verbose: bool | None
    progressbar_kwargs: dict[str, Any] | None
    backend: Any
    n_workers: int | None
    n_parts: int | None
    unsupported_backend: bool | None


@runtime_checkable
class AnalysisBaseLike(Protocol):
    """Structural protocol for MDAnalysis ``AnalysisBase``-style objects."""

    results: Any

    def run(self, **kwargs: Any) -> AnalysisBaseLike:
        """Run the analysis and return the analysis object.

        Parameters
        ----------
        **kwargs : Any
            Keyword arguments forwarded to the wrapped analysis object's
            ``run()`` method.

        Returns
        -------
        AnalysisBaseLike
            The analysis object after execution.
        """


class MDAnalysisExtensionError(RuntimeError):
    """Base runtime error for MDAnalysis extension-layer failures."""
