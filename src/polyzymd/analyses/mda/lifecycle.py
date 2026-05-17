"""Lifecycle bridge between PolyzyMD replicates and MDAnalysis jobs."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from pydantic import BaseModel

from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.mda.artifacts import ReplicateArtifact, raw_mdanalysis_results_path
from polyzymd.analyses.mda.frame_selection import FrameSelection
from polyzymd.analyses.mda.job import MDAAnalysisJob, MDAJobResult, MDAUniversePolicy
from polyzymd.analyses.mda.plugin import MDACollectorContext
from polyzymd.analyses.mda.store import ArtifactStore, ArtifactStoreError

if TYPE_CHECKING:
    from polyzymd.analyses._contexts import ReplicateContext
    from polyzymd.analyses.mda.job import MDABackendPolicy

logger = logging.getLogger("polyzymd.analyses")


@dataclass(frozen=True)
class MDAReplicateJobContext:
    """Context passed to ``Analysis.build_mda_jobs()`` for one replicate."""

    replicate_context: ReplicateContext
    universe: Any
    frame_selection: FrameSelection
    universe_policy: MDAUniversePolicy
    artifact_store: ArtifactStore

    @property
    def output_dir(self) -> Path:
        """Return the replicate output directory.

        Returns
        -------
        Path
            Directory owned by this replicate analysis run.
        """

        return self.replicate_context.output_dir

    @property
    def replicate(self) -> int:
        """Return the one-indexed replicate ID.

        Returns
        -------
        int
            Replicate ID from the framework context.
        """

        return self.replicate_context.replicate

    @property
    def settings(self) -> BaseModel:
        """Return resolved plugin settings.

        Returns
        -------
        BaseModel
            Settings model supplied by the public lifecycle.
        """

        return self.replicate_context.settings

    @property
    def backend_policy(self) -> MDABackendPolicy:
        """Return the MDAnalysis backend policy for job construction.

        Returns
        -------
        MDABackendPolicy
            Policy resolved from comparison configuration, or the serial default.
        """

        return self.replicate_context.backend_policy


def build_mda_replicate_job_context(
    analysis: Any, ctx: ReplicateContext, replicate: int
) -> MDAReplicateJobContext:
    """Build the MDAnalysis job context for one replicate.

    Parameters
    ----------
    analysis : Any
        Analysis instance requesting an MDAnalysis job context.
    ctx : ReplicateContext
        Framework-provided replicate context.
    replicate : int
        One-indexed replicate ID.

    Returns
    -------
    MDAReplicateJobContext
        Context containing the loaded universe, resolved frame selection,
        universe policy, and artifact store.
    """

    loader = analysis._trajectory_loader_factory()(ctx.sim_config)
    provider = _build_universe_provider(analysis, ctx, loader)
    universe = provider.load_universe(replicate)
    window = analysis.get_trajectory_window(ctx, replicate, loader, universe)
    if getattr(window, "warning_message", None):
        logger.warning(
            "%s: %s [condition=%s, replicate=%d]",
            analysis.name,
            window.warning_message,
            ctx.condition.label,
            replicate,
        )
    frame_selection = FrameSelection.from_trajectory_window(window)
    provenance = _provenance_for(provider, replicate)
    universe_policy = MDAUniversePolicy(
        condition_label=ctx.condition.label,
        replicate=replicate,
        provenance=provenance,
        metadata={"equilibration": ctx.equilibration},
    )
    artifact_store = analysis._mda_artifact_store_factory()(ctx.output_dir)
    return MDAReplicateJobContext(
        replicate_context=ctx,
        universe=universe,
        frame_selection=frame_selection,
        universe_policy=universe_policy,
        artifact_store=artifact_store,
    )


def run_mda_replicate_jobs(
    analysis: Any, ctx: ReplicateContext, replicate: int
) -> ReplicateArtifact | None:
    """Run MDAnalysis jobs for one replicate and collect a strict artifact.

    Parameters
    ----------
    analysis : Any
        Analysis instance with a ``build_mda_jobs()`` hook.
    ctx : ReplicateContext
        Framework-provided replicate context.
    replicate : int
        One-indexed replicate ID.

    Returns
    -------
    ReplicateArtifact or None
        Collected replicate artifact, or ``None`` when the hook declines the MDA
        path.
    """

    mda_ctx = build_mda_replicate_job_context(analysis, ctx, replicate)
    jobs = analysis.build_mda_jobs(mda_ctx)
    if jobs is None:
        return None
    normalized_jobs = _validate_jobs(jobs, analysis_name=analysis.name)
    completed_jobs = [job.run() for job in normalized_jobs]
    return _artifact_from_completed_jobs(analysis, mda_ctx, completed_jobs)


def _build_universe_provider(analysis: Any, ctx: ReplicateContext, loader: Any) -> Any:
    """Create the universe provider using the analysis injection hook.

    Parameters
    ----------
    analysis : Any
        Analysis instance supplying the provider factory.
    ctx : ReplicateContext
        Framework-provided replicate context.
    loader : Any
        Shared trajectory loader already allocated for this replicate.

    Returns
    -------
    Any
        Universe provider compatible with ``UniverseProvider``.
    """

    provider_factory = analysis._mda_universe_provider_factory()
    if hasattr(provider_factory, "from_config"):
        return provider_factory.from_config(ctx.sim_config, loader=loader)
    return provider_factory(ctx.sim_config, loader=loader)


def _provenance_for(provider: Any, replicate: int) -> Any:
    """Return provider provenance when the provider exposes it.

    Parameters
    ----------
    provider : Any
        Universe provider used to load the replicate universe.
    replicate : int
        One-indexed replicate ID.

    Returns
    -------
    Any
        Provider-specific provenance object, or ``None`` when unavailable.
    """

    if hasattr(provider, "provenance_for"):
        return provider.provenance_for(replicate)
    if hasattr(provider, "get_provenance"):
        return provider.get_provenance(replicate)
    return None


def _validate_jobs(
    jobs: Sequence[MDAAnalysisJob], *, analysis_name: str
) -> tuple[MDAAnalysisJob, ...]:
    """Validate job-builder output.

    Parameters
    ----------
    jobs : sequence of MDAAnalysisJob
        Jobs returned by ``build_mda_jobs()``.
    analysis_name : str
        Analysis name for diagnostics.

    Returns
    -------
    tuple[MDAAnalysisJob, ...]
        Concrete job tuple ready for execution.
    """

    if isinstance(jobs, (str, bytes)) or not isinstance(jobs, Sequence):
        raise PluginContractError(
            f"{analysis_name}.build_mda_jobs() must return a sequence of MDAAnalysisJob objects"
        )
    normalized_jobs = tuple(jobs)
    if not normalized_jobs:
        raise PluginContractError(f"{analysis_name}.build_mda_jobs() returned no jobs")
    invalid = [job for job in normalized_jobs if not isinstance(job, MDAAnalysisJob)]
    if invalid:
        raise PluginContractError(
            f"{analysis_name}.build_mda_jobs() returned {type(invalid[0]).__name__}; "
            "expected MDAAnalysisJob"
        )
    return normalized_jobs


def _artifact_from_completed_jobs(
    analysis: Any,
    ctx: MDAReplicateJobContext,
    completed_jobs: Sequence[MDAJobResult],
) -> ReplicateArtifact:
    """Collect completed MDA job results through the analysis collector.

    Parameters
    ----------
    analysis : Any
        Analysis instance that owns the jobs.
    ctx : MDAReplicateJobContext
        MDA replicate context used for execution.
    completed_jobs : sequence of MDAJobResult
        Completed job result references.

    Returns
    -------
    ReplicateArtifact
        Strict JSON-compatible artifact envelope.
    """

    collector_ctx = _build_collector_context(analysis, ctx)
    collector = analysis.build_mda_collector(collector_ctx)
    if not callable(collector):
        raise PluginContractError(
            f"{analysis.name}.build_mda_collector() must return a callable collector, "
            f"got {type(collector).__name__}"
        )
    artifact = collector(collector_ctx, tuple(completed_jobs))
    _validate_collected_artifact(artifact, collector_ctx)
    return artifact


def _build_collector_context(analysis: Any, ctx: MDAReplicateJobContext) -> MDACollectorContext:
    """Build collector context from a completed MDA replicate context.

    Parameters
    ----------
    analysis : Any
        Analysis instance that owns the jobs.
    ctx : MDAReplicateJobContext
        MDA replicate context used for execution.

    Returns
    -------
    MDACollectorContext
        Context passed to the artifact collector.
    """

    settings_fingerprint = analysis.aggregate_settings_fingerprint(ctx.settings)
    return MDACollectorContext(
        analysis_name=analysis.name,
        replicate_context=ctx.replicate_context,
        frame_selection=ctx.frame_selection,
        universe_policy=ctx.universe_policy,
        artifact_store=ctx.artifact_store,
        settings_fingerprint=settings_fingerprint,
        warnings=_collector_warnings(ctx),
    )


def _collector_warnings(ctx: MDAReplicateJobContext) -> tuple[str, ...]:
    """Return warning messages known before collection.

    Parameters
    ----------
    ctx : MDAReplicateJobContext
        MDA replicate context used for execution.

    Returns
    -------
    tuple of str
        Warning messages from frame selection and universe provenance.
    """

    warnings: list[str] = []
    if ctx.frame_selection.warning_message:
        warnings.append(ctx.frame_selection.warning_message)
    policy = ctx.universe_policy.as_dict()
    provenance = policy.get("provenance") if isinstance(policy, dict) else None
    provider_warnings = provenance.get("warnings", []) if isinstance(provenance, dict) else []
    if isinstance(provider_warnings, list):
        warnings.extend(str(warning) for warning in provider_warnings)
    return tuple(warnings)


def _validate_collected_artifact(
    artifact: Any,
    ctx: MDACollectorContext,
) -> None:
    """Validate collector output before lifecycle persistence.

    Parameters
    ----------
    artifact : Any
        Collector output to validate.
    ctx : MDACollectorContext
        Collector context with expected identity and artifact store.
    """

    if not isinstance(artifact, ReplicateArtifact):
        raise PluginContractError(
            f"{ctx.analysis_name}.build_mda_collector() returned {type(artifact).__name__}; "
            "expected ReplicateArtifact"
        )
    expected_identity = {
        "analysis_name": ctx.analysis_name,
        "condition_label": ctx.condition_label,
        "replicate": ctx.replicate,
    }
    actual_identity = {
        "analysis_name": artifact.analysis_name,
        "condition_label": artifact.condition_label,
        "replicate": artifact.replicate,
    }
    if actual_identity != expected_identity:
        raise PluginContractError(
            f"{ctx.analysis_name}.build_mda_collector() returned artifact identity "
            f"{actual_identity!r}; expected {expected_identity!r}"
        )
    raw_path = raw_mdanalysis_results_path(artifact)
    if raw_path is not None:
        raise PluginContractError(
            f"{ctx.analysis_name}.build_mda_collector() returned raw MDAnalysis Results "
            f"at {raw_path}; map Results to JSON primitives or sidecars before persistence"
        )
    for sidecar in artifact.sidecars:
        try:
            ctx.artifact_store.validate_sidecar(sidecar)
        except ArtifactStoreError as exc:
            raise PluginContractError(
                f"{ctx.analysis_name}.build_mda_collector() returned invalid sidecar "
                f"{sidecar.path!r}: {exc}"
            ) from exc
