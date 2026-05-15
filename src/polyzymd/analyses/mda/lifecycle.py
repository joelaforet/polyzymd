"""Lifecycle bridge between PolyzyMD replicates and MDAnalysis jobs."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from pydantic import BaseModel

from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.mda.artifacts import ReplicateArtifact
from polyzymd.analyses.mda.frame_selection import FrameSelection
from polyzymd.analyses.mda.job import MDAAnalysisJob, MDAJobResult, MDAUniversePolicy
from polyzymd.analyses.mda.store import ArtifactStore

if TYPE_CHECKING:
    from polyzymd.analyses._contexts import ReplicateContext

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
    """Collect completed MDA job results into a replicate artifact.

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

    job_payloads = [_job_result_payload(job, analysis_name=analysis.name) for job in completed_jobs]
    warnings = []
    if ctx.frame_selection.warning_message:
        warnings.append(ctx.frame_selection.warning_message)
    provenance = _json_payload(ctx.universe_policy.as_dict(), analysis_name=analysis.name)
    provider_warnings = (
        provenance.get("provenance", {}).get("warnings", []) if isinstance(provenance, dict) else []
    )
    if isinstance(provider_warnings, list):
        warnings.extend(str(warning) for warning in provider_warnings)
    return ReplicateArtifact(
        analysis_name=analysis.name,
        condition_label=ctx.replicate_context.condition.label,
        replicate=ctx.replicate,
        payload={"jobs": job_payloads, "n_jobs": len(job_payloads)},
        provenance={
            "source": "mda_job_lifecycle",
            "frame_selection": _frame_selection_payload(ctx.frame_selection),
            "universe_policy": provenance,
        },
        metadata={"result_kind": "mda_replicate_jobs"},
        warnings=warnings,
    )


def _job_result_payload(job: MDAJobResult, *, analysis_name: str) -> dict[str, Any]:
    """Serialize one completed job result to JSON-compatible primitives.

    Parameters
    ----------
    job : MDAJobResult
        Completed job result reference.
    analysis_name : str
        Analysis name for diagnostics.

    Returns
    -------
    dict[str, Any]
        JSON-compatible job payload.
    """

    return {
        "name": job.name,
        "results": _json_payload(job.results, analysis_name=f"{analysis_name}.{job.name}"),
        "run_kwargs": _json_payload(
            dict(job.run_kwargs), analysis_name=f"{analysis_name}.{job.name}"
        ),
        "frame_selection": _frame_selection_payload(job.frame_selection),
        "backend_policy": _json_payload(
            job.backend_policy.run_kwargs(), analysis_name=f"{analysis_name}.{job.name}"
        ),
        "universe_policy": _json_payload(
            job.universe_policy.as_dict(), analysis_name=f"{analysis_name}.{job.name}"
        ),
    }


def _frame_selection_payload(frame_selection: FrameSelection) -> dict[str, Any]:
    """Serialize frame-selection provenance to primitive values.

    Parameters
    ----------
    frame_selection : FrameSelection
        Frame selection used for a job or replicate context.

    Returns
    -------
    dict[str, Any]
        JSON-compatible frame-selection metadata.
    """

    return {
        "start": frame_selection.start,
        "stop": frame_selection.stop,
        "step": frame_selection.step,
        "frames": _json_payload(frame_selection.frames, analysis_name="frame_selection"),
        "equilibration": frame_selection.equilibration,
        "equilibration_start": frame_selection.equilibration_start,
        "equilibration_ps": frame_selection.equilibration_ps,
        "timestep_ps": frame_selection.timestep_ps,
        "n_frames_total": frame_selection.n_frames_total,
        "n_frames_selected": frame_selection.n_frames_selected,
        "warning_message": frame_selection.warning_message,
    }


def _json_payload(value: Any, *, analysis_name: str) -> Any:
    """Convert supported values to strict JSON-compatible primitives.

    Parameters
    ----------
    value : Any
        Candidate payload returned by an MDA job.
    analysis_name : str
        Analysis name for diagnostics.

    Returns
    -------
    Any
        JSON-compatible primitive, list, or dictionary.
    """

    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise PluginContractError(
                f"{analysis_name}.build_mda_jobs() produced non-finite float result {value!r}; "
                "add a collector in a future MDA lifecycle task"
            )
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, BaseModel):
        return _json_payload(value.model_dump(mode="json"), analysis_name=analysis_name)
    if isinstance(value, Mapping):
        payload = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise PluginContractError(
                    f"{analysis_name}.build_mda_jobs() produced non-string mapping key "
                    f"{key!r}; add a collector in a future MDA lifecycle task"
                )
            payload[key] = _json_payload(item, analysis_name=analysis_name)
        return payload
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_json_payload(item, analysis_name=analysis_name) for item in value]
    raise PluginContractError(
        f"{analysis_name}.build_mda_jobs() produced non-JSON-serializable "
        f"{type(value).__name__} results; add a collector in a future MDA lifecycle task"
    )
