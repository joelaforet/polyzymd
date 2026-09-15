"""Run one replicate of a contract plugin and collect its artifact.

The lifecycle owns everything around ``plugin.compute()``: it loads the
universe, resolves the production window, builds the artifact store, calls the
analysis, and checks the artifact that comes back before the framework writes
it.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from pydantic import BaseModel

from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.mda.artifacts import ReplicateArtifact, raw_mdanalysis_results_path
from polyzymd.analyses.mda.frame_selection import (
    FrameSelection,
    _normalize_frame_selector_values,
    _normalize_scalar_value,
)
from polyzymd.analyses.mda.store import ArtifactStore, ArtifactStoreError
from polyzymd.analyses.mda.universe import UniverseProvider
from polyzymd.analyses.shared.loader import TrajectoryLoader

if TYPE_CHECKING:
    from polyzymd.analyses._framework.contexts import ReplicateContext

logger = logging.getLogger("polyzymd.analyses")


@dataclass(frozen=True)
class MDAUniversePolicy:
    """Identity of the files one replicate was measured from.

    It carries no universe, only the provenance the loader already recorded, so
    an artifact can say which topology and trajectories produced its numbers.
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
    """What one call to a plugin's ``compute()`` produced for one replicate."""

    name: str
    results: Any
    frame_selection: FrameSelection
    universe_policy: MDAUniversePolicy


@dataclass(frozen=True)
class MDAReplicateJobContext:
    """What the lifecycle hands a plugin for one replicate.

    It carries the loaded universe, the production window, the identity of the
    files behind them, and the store the plugin writes sidecars to. The
    warnings are the ones the framework already knows about, from frame
    selection and from the universe provider; the plugin's own warnings reach
    the artifact through the metadata of its observables.
    """

    replicate_context: ReplicateContext
    universe: Any
    frame_selection: FrameSelection
    universe_policy: MDAUniversePolicy
    artifact_store: ArtifactStore
    warnings: tuple[str, ...] = ()

    @property
    def output_dir(self) -> Path:
        """Directory owned by this replicate analysis run."""
        return self.replicate_context.output_dir

    @property
    def replicate(self) -> int:
        """One-indexed replicate ID."""
        return self.replicate_context.replicate

    @property
    def condition_label(self) -> str:
        """Label of the condition this replicate belongs to."""
        return self.replicate_context.condition.label

    @property
    def settings(self) -> BaseModel:
        """Resolved plugin settings."""
        return self.replicate_context.settings


def build_trajectory_loader(sim_config: Any) -> TrajectoryLoader:
    """Create the trajectory loader used for one condition."""
    return TrajectoryLoader(sim_config)


def run_replicate(analysis: Any, ctx: ReplicateContext, replicate: int) -> ReplicateArtifact:
    """Compute one replicate and return the artifact the framework will write.

    Parameters
    ----------
    analysis : Analysis
        Analysis running the plugin.
    ctx : ReplicateContext
        Framework-provided replicate context.
    replicate : int
        One-indexed replicate ID.

    Returns
    -------
    ReplicateArtifact
        Artifact holding the reduced observables, the sidecars and the identity
        block.
    """
    job_ctx = _job_context(analysis, ctx, replicate)
    measured = analysis.measure_replicate(job_ctx)
    artifact = analysis.collect_replicate(job_ctx, measured)
    _validate_artifact(artifact, job_ctx, analysis_name=analysis.name)
    return artifact


def _job_context(analysis: Any, ctx: ReplicateContext, replicate: int) -> MDAReplicateJobContext:
    """Load the universe and resolve the production window for one replicate."""
    loader = build_trajectory_loader(ctx.sim_config)
    provider = _build_universe_provider(ctx, loader)
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
    universe_policy = MDAUniversePolicy(
        condition_label=ctx.condition.label,
        replicate=replicate,
        provenance=_provenance_for(provider, replicate),
        metadata={"equilibration": ctx.equilibration},
    )
    return MDAReplicateJobContext(
        replicate_context=ctx,
        universe=universe,
        frame_selection=frame_selection,
        universe_policy=universe_policy,
        artifact_store=ArtifactStore(ctx.output_dir),
        warnings=_known_warnings(frame_selection, universe_policy),
    )


def _build_universe_provider(ctx: ReplicateContext, loader: Any) -> Any:
    """Create the universe provider for one replicate."""
    return UniverseProvider.from_config(ctx.sim_config, loader=loader)


def _provenance_for(provider: Any, replicate: int) -> Any:
    """Return provider provenance when the provider exposes it."""
    if hasattr(provider, "provenance_for"):
        return provider.provenance_for(replicate)
    if hasattr(provider, "get_provenance"):
        return provider.get_provenance(replicate)
    return None


def _known_warnings(
    frame_selection: FrameSelection, universe_policy: MDAUniversePolicy
) -> tuple[str, ...]:
    """Warnings the framework already holds before the plugin runs."""
    warnings: list[str] = []
    if frame_selection.warning_message:
        warnings.append(frame_selection.warning_message)
    policy = universe_policy.as_dict()
    provenance = policy.get("provenance") if isinstance(policy, dict) else None
    provider_warnings = provenance.get("warnings", []) if isinstance(provenance, dict) else []
    if isinstance(provider_warnings, list):
        warnings.extend(str(warning) for warning in provider_warnings)
    return tuple(warnings)


def _validate_artifact(artifact: Any, ctx: MDAReplicateJobContext, *, analysis_name: str) -> None:
    """Check the artifact's identity, its payload and its sidecars.

    A raw MDAnalysis ``Results`` object anywhere in the artifact is rejected
    rather than serialized, because it pickles a whole universe into the JSON.
    """
    if not isinstance(artifact, ReplicateArtifact):
        raise PluginContractError(
            f"{analysis_name} produced {type(artifact).__name__}; expected ReplicateArtifact"
        )
    expected = {
        "analysis_name": analysis_name,
        "condition_label": ctx.condition_label,
        "replicate": ctx.replicate,
    }
    actual = {
        "analysis_name": artifact.analysis_name,
        "condition_label": artifact.condition_label,
        "replicate": artifact.replicate,
    }
    if actual != expected:
        raise PluginContractError(
            f"{analysis_name} produced artifact identity {actual!r}; expected {expected!r}"
        )
    raw_path = raw_mdanalysis_results_path(artifact)
    if raw_path is not None:
        raise PluginContractError(
            f"{analysis_name} produced raw MDAnalysis Results at {raw_path}; map them to JSON "
            "primitives or sidecars before persistence"
        )
    for sidecar in artifact.sidecars:
        try:
            ctx.artifact_store.validate_sidecar(sidecar)
        except ArtifactStoreError as exc:
            raise PluginContractError(
                f"{analysis_name} produced an invalid sidecar {sidecar.path!r}: {exc}"
            ) from exc


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
        "start": _normalize_scalar_value(frame_selection.start),
        "stop": _normalize_scalar_value(frame_selection.stop),
        "step": _normalize_scalar_value(frame_selection.step),
        "frames": _frame_selector_payload(frame_selection.frames),
        "equilibration": frame_selection.equilibration,
        "equilibration_start": _normalize_scalar_value(frame_selection.equilibration_start),
        "equilibration_ps": _normalize_scalar_value(frame_selection.equilibration_ps),
        "timestep_ps": _normalize_scalar_value(frame_selection.timestep_ps),
        "first_frame_time_ps": _normalize_scalar_value(frame_selection.first_frame_time_ps),
        "selected_start_time_ps": _normalize_scalar_value(frame_selection.selected_start_time_ps),
        "equilibration_time_reference": frame_selection.equilibration_time_reference,
        "n_frames_total": _normalize_scalar_value(frame_selection.n_frames_total),
        "n_frames_selected": _normalize_scalar_value(frame_selection.n_frames_selected),
        "warning_message": frame_selection.warning_message,
    }


def _frame_selector_payload(frames: Any) -> list[int | bool] | None:
    """Serialize explicit frame selectors to JSON-safe Python scalars."""

    if frames is None:
        return None
    try:
        return _normalize_frame_selector_values(frames)
    except ValueError as exc:
        raise PluginContractError(
            "The plugin produced a non-integer explicit frame selector; "
            "use integer indices or a boolean mask"
        ) from exc
