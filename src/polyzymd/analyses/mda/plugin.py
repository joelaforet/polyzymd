"""Collector interfaces for MDAnalysis job outputs."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from operator import index as operator_index
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol

from pydantic import BaseModel

from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.mda.artifacts import (
    ReplicateArtifact,
    raw_mdanalysis_results_path,
)
from polyzymd.analyses.mda.frame_selection import FrameSelection
from polyzymd.analyses.mda.job import MDAJobResult, MDAUniversePolicy
from polyzymd.analyses.mda.store import ArtifactStore

if TYPE_CHECKING:
    from polyzymd.analyses._framework.contexts import ReplicateContext


@dataclass(frozen=True)
class MDACollectorContext:
    """Context supplied to an MDAnalysis artifact collector.

    The context contains framework identity and provenance for one replicate so
    collectors can map raw job results into PolyzyMD-owned artifacts without
    reaching back into the orchestrator.
    """

    analysis_name: str
    replicate_context: ReplicateContext
    frame_selection: FrameSelection
    universe_policy: MDAUniversePolicy
    artifact_store: ArtifactStore
    settings_fingerprint: str | None = None
    warnings: Sequence[str] = ()

    def __post_init__(self) -> None:
        """Freeze warning messages as strings for artifact reuse."""

        object.__setattr__(self, "warnings", tuple(str(warning) for warning in self.warnings))

    @property
    def condition_label(self) -> str:
        """Return the simulation condition label.

        Returns
        -------
        str
            Condition label from the framework context.
        """

        return self.replicate_context.condition.label

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
    def output_dir(self) -> Path:
        """Return the replicate output directory.

        Returns
        -------
        Path
            Directory owned by this replicate analysis run.
        """

        return self.replicate_context.output_dir

    @property
    def result_path(self) -> Path:
        """Return the canonical replicate result path.

        Returns
        -------
        Path
            Canonical artifact JSON path for this replicate.
        """

        return self.replicate_context.result_path

    @property
    def settings(self) -> BaseModel:
        """Return resolved plugin settings.

        Returns
        -------
        BaseModel
            Settings model supplied by the public lifecycle.
        """

        return self.replicate_context.settings


class MDAArtifactCollector(Protocol):
    """Protocol for converting completed MDAnalysis jobs to an artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[MDAJobResult],
    ) -> ReplicateArtifact:
        """Collect completed jobs into one replicate artifact.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context for one replicate.
        completed_jobs : sequence of MDAJobResult
            Completed MDAnalysis-compatible jobs.

        Returns
        -------
        ReplicateArtifact
            PolyzyMD-owned artifact for this replicate.
        """


class StrictJSONMDAResultCollector:
    """Default collector for jobs that already return strict JSON values.

    This collector preserves the P2-001 simple job behavior while rejecting raw
    MDAnalysis ``Results`` containers and other non-JSON values. Analyses with
    rich ``Results`` objects should implement a custom collector that maps data
    to primitive payloads and sidecars.
    """

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[MDAJobResult],
    ) -> ReplicateArtifact:
        """Collect completed jobs into a strict JSON artifact.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context for one replicate.
        completed_jobs : sequence of MDAJobResult
            Completed MDAnalysis-compatible jobs.

        Returns
        -------
        ReplicateArtifact
            JSON-compatible replicate artifact.
        """

        job_payloads = [
            _job_result_payload(job, analysis_name=ctx.analysis_name) for job in completed_jobs
        ]
        provenance = strict_json_payload(
            ctx.universe_policy.as_dict(), analysis_name=ctx.analysis_name
        )
        metadata: dict[str, Any] = {"result_kind": "mda_replicate_jobs"}
        if ctx.settings_fingerprint is not None:
            metadata["settings_fingerprint"] = ctx.settings_fingerprint
        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={"jobs": job_payloads, "n_jobs": len(job_payloads)},
            provenance={
                "source": "mda_job_lifecycle",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": provenance,
            },
            metadata=metadata,
            warnings=list(ctx.warnings),
        )


def frame_selection_payload(frame_selection: FrameSelection) -> dict[str, Any]:
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
        "frames": _frame_selector_payload(frame_selection.frames),
        "equilibration": frame_selection.equilibration,
        "equilibration_start": frame_selection.equilibration_start,
        "equilibration_ps": frame_selection.equilibration_ps,
        "timestep_ps": frame_selection.timestep_ps,
        "first_frame_time_ps": frame_selection.first_frame_time_ps,
        "selected_start_time_ps": frame_selection.selected_start_time_ps,
        "equilibration_time_reference": frame_selection.equilibration_time_reference,
        "n_frames_total": frame_selection.n_frames_total,
        "n_frames_selected": frame_selection.n_frames_selected,
        "warning_message": frame_selection.warning_message,
    }


def _frame_selector_payload(frames: Any) -> list[int | bool] | None:
    """Serialize explicit frame selectors to JSON-safe Python scalars."""

    if frames is None:
        return None
    payload: list[int | bool] = []
    for frame in frames:
        if _is_boolean_frame_value(frame):
            payload.append(bool(frame))
        else:
            try:
                payload.append(operator_index(frame))
            except TypeError as exc:
                raise PluginContractError(
                    "frame_selection.build_mda_jobs() produced non-integer explicit frame "
                    f"selector {frame!r}; use integer indices or a boolean mask"
                ) from exc
    return payload


def _is_boolean_frame_value(frame: Any) -> bool:
    """Return whether a frame selector value is boolean-like."""

    if isinstance(frame, bool):
        return True
    frame_type = type(frame)
    return (
        frame_type.__name__ == "bool_"
        and frame_type.__module__.split(".", maxsplit=1)[0] == "numpy"
    )


def strict_json_payload(value: Any, *, analysis_name: str) -> Any:
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

    raw_path = raw_mdanalysis_results_path(value)
    if raw_path is not None:
        raise PluginContractError(
            f"{analysis_name}.build_mda_jobs() produced raw MDAnalysis Results at {raw_path}; "
            "implement build_mda_collector() to map Results to JSON primitives or sidecars"
        )
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise PluginContractError(
                f"{analysis_name}.build_mda_jobs() produced non-finite float result {value!r}; "
                "implement build_mda_collector() to map Results to JSON primitives or sidecars"
            )
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, BaseModel):
        return strict_json_payload(value.model_dump(mode="json"), analysis_name=analysis_name)
    if isinstance(value, Mapping):
        payload = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise PluginContractError(
                    f"{analysis_name}.build_mda_jobs() produced non-string mapping key "
                    f"{key!r}; implement build_mda_collector() to map Results to JSON "
                    "primitives or sidecars"
                )
            payload[key] = strict_json_payload(item, analysis_name=analysis_name)
        return payload
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [strict_json_payload(item, analysis_name=analysis_name) for item in value]
    raise PluginContractError(
        f"{analysis_name}.build_mda_jobs() produced non-JSON-serializable "
        f"{type(value).__name__} results; implement build_mda_collector() to map Results "
        "to JSON primitives or sidecars"
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
        "results": strict_json_payload(job.results, analysis_name=f"{analysis_name}.{job.name}"),
        "run_kwargs": strict_json_payload(
            dict(job.run_kwargs), analysis_name=f"{analysis_name}.{job.name}"
        ),
        "frame_selection": frame_selection_payload(job.frame_selection),
        "backend_policy": strict_json_payload(
            job.backend_policy.run_kwargs(), analysis_name=f"{analysis_name}.{job.name}"
        ),
        "universe_policy": strict_json_payload(
            job.universe_policy.as_dict(), analysis_name=f"{analysis_name}.{job.name}"
        ),
    }
