"""Lifecycle helpers for the contacts analysis facade."""

from __future__ import annotations

import logging
from typing import Any, Sequence

from polyzymd.analyses.base import AggregateContext, ReplicateContext
from polyzymd.analyses.contacts._aggregator import (
    CONTACTS_LEGACY_RECOMPUTE_GUIDANCE,
    aggregate_contact_artifacts,
)

logger = logging.getLogger("polyzymd.analyses.contacts")


def align_replicate_results(results: Sequence[Any], replicates: Sequence[int]) -> list[Any]:
    """Return replicate results ordered to match the requested IDs.

    Parameters
    ----------
    results : Sequence[Any]
        Per-replicate results returned by the orchestrator.
    replicates : Sequence[int]
        Requested replicate IDs from the aggregate context.

    Returns
    -------
    list[Any]
        Results ordered by ``replicates``.

    Raises
    ------
    ValueError
        If any result is missing a replicate ID, duplicates an ID, or does not
        match the requested replicate set.
    """

    requested = tuple(int(rep) for rep in replicates)
    by_replicate: dict[int, Any] = {}
    for result in results:
        stored_replicate = getattr(result, "replicate", None)
        try:
            replicate = int(stored_replicate)
        except (TypeError, ValueError) as exc:
            raise ValueError("contacts aggregate input lacks a valid replicate ID") from exc
        if replicate in by_replicate:
            raise ValueError(f"contacts aggregate input duplicates replicate {replicate}")
        by_replicate[replicate] = result

    requested_set = set(requested)
    stored_set = set(by_replicate)
    if stored_set != requested_set:
        raise ValueError(
            "contacts aggregate input replicate IDs do not match requested replicates: "
            f"stored={sorted(stored_set)}, requested={sorted(requested_set)}"
        )
    return [by_replicate[replicate] for replicate in requested]


def get_trajectory_window(
    ctx: ReplicateContext,
    replicate: int,
    loader: Any,
    universe: Any,
) -> Any:
    """Resolve the contacts analysis window with a dt-aware fallback.

    Parameters
    ----------
    ctx : ReplicateContext
        Framework-provided context.
    replicate : int
        Replicate ID accepted for lifecycle signature compatibility.
    loader : Any
        Trajectory loader used for timestep fallback.
    universe : Any
        Loaded trajectory universe.

    Returns
    -------
    Any
        Resolved trajectory window.
    """

    from polyzymd.analyses.shared.window import resolve_trajectory_window

    del replicate
    try:
        timestep_ps = float(universe.trajectory.dt)
        if timestep_ps <= 0:
            raise ValueError
    except (AttributeError, TypeError, ValueError):
        try:
            timestep_ps = float(loader.get_timestep(ctx.replicate, unit="ps"))
        except (AttributeError, TypeError, ValueError):
            timestep_ps = 1.0

    return resolve_trajectory_window(
        equilibration=ctx.equilibration,
        n_frames_total=len(universe.trajectory),
        timestep_ps=timestep_ps,
    )


def aggregate(analysis: Any, ctx: AggregateContext, results: Sequence[Any]) -> Any:
    """Aggregate contact artifacts across replicates for one condition.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for cache, path, and save hooks.
    ctx : AggregateContext
        Framework-provided context.
    results : Sequence[Any]
        Per-replicate contact results.

    Returns
    -------
    Any
        Aggregated result with per-residue statistics.
    """

    from polyzymd.analyses.mda import ReplicateArtifact

    del analysis
    if not results or not all(isinstance(result, ReplicateArtifact) for result in results):
        raise ValueError(CONTACTS_LEGACY_RECOMPUTE_GUIDANCE)
    logger.info("  Aggregating %d contacts replicate artifacts...", len(results))
    return aggregate_contact_artifacts(results, ctx)
