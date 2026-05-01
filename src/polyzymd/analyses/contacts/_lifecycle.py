"""Lifecycle helpers for the contacts analysis facade."""

from __future__ import annotations

import logging
from typing import Any, Sequence

from polyzymd.analyses.base import AggregateContext, Analysis, ReplicateContext
from polyzymd.analyses.contacts._aggregator import AggregatedContactResult
from polyzymd.analyses.contacts._identity import (
    CONTACTS_SETTINGS_FINGERPRINT_DOMAIN,
    contacts_settings_fingerprint,
    normalize_polymer_types,
)
from polyzymd.analyses.contacts._results import ContactResult
from polyzymd.analyses.contacts._runner import ContactsReplicateRunner, build_contact_grouping
from polyzymd.analyses.exceptions import ReplicateSkippedError

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


def run_replicate(analysis: Any, ctx: ReplicateContext, replicate: int) -> Any:
    """Run contacts for a single replicate.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for cache, path, and save hooks.
    ctx : ReplicateContext
        Framework-provided context.
    replicate : int
        1-indexed replicate number.

    Returns
    -------
    Any
        Per-replicate contact result.
    """

    result_file = analysis._replicate_sidecar_path(
        ctx.output_dir,
        ctx.settings,
        ctx.equilibration,
        replicate,
    )
    canonical_file = ctx.result_path or ctx.output_dir / "result.json"
    result_candidates = ()
    if not ctx.recompute:
        result_candidates = analysis._replicate_cache_candidates(
            ctx.output_dir,
            ctx.settings,
            ctx.equilibration,
            replicate,
            result_path=ctx.result_path,
        )

    for candidate in result_candidates:
        cached = analysis._load_cache_candidate(
            ContactResult,
            candidate,
            recompute=ctx.recompute,
            sim_config=ctx.sim_config,
            settings=None,
        )
        if (
            cached is not None
            and analysis._cache_matches_window(cached, ctx.equilibration)
            and analysis._cache_matches_replicate_id(cached, replicate, source=candidate)
            and analysis._cache_matches_contacts_settings(cached, ctx.settings, source=candidate)
        ):
            if (
                ctx.result_path is not None
                and candidate != ctx.result_path
                and not ctx.result_path.exists()
            ):
                analysis._attach_contacts_identity_metadata(cached, ctx.settings)
                analysis.save_result(cached, ctx.result_path)
            if candidate == canonical_file and not result_file.exists():
                analysis._attach_contacts_identity_metadata(cached, ctx.settings)
                cached.save(result_file)
            return cached

    try:
        result = Analysis.run_replicate(analysis, ctx, replicate)
    except FileNotFoundError as e:
        logger.warning(f"  Skipping replicate {replicate}: trajectory data not found. {e}")
        raise ReplicateSkippedError(
            f"No trajectory data found for replicate {replicate}: {e}"
        ) from e
    except (ValueError, OSError, IndexError) as e:
        logger.debug("Full traceback:", exc_info=True)
        logger.warning(f"  Skipping replicate {replicate}: analysis failed with error: {e}")
        raise ReplicateSkippedError(
            f"Contacts analysis failed for replicate {replicate}: {type(e).__name__}: {e}"
        ) from e

    ctx.output_dir.mkdir(parents=True, exist_ok=True)
    analysis._attach_contacts_identity_metadata(result, ctx.settings)
    result.save(result_file)
    if ctx.result_path is not None:
        analysis.save_result(result, ctx.result_path)
    logger.info(f"  Saved: {result_file}")
    return result


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


def build_runner(
    analysis: Any,
    ctx: ReplicateContext,
    replicate: int,
    universe: Any,
    window: Any,
) -> Any:
    """Build the runner-backed contacts execution object.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for selection helpers.
    ctx : ReplicateContext
        Framework-provided context.
    replicate : int
        Replicate ID accepted for lifecycle signature compatibility.
    universe : Any
        Loaded trajectory universe.
    window : Any
        Resolved trajectory window accepted for signature compatibility.

    Returns
    -------
    Any
        Runner-backed contacts execution object.
    """

    import polyzymd.analyses.contacts as contacts_facade
    from polyzymd.analyses.shared.selectors import MDAnalysisSelector

    del replicate, window
    settings = ctx.settings
    polymer_selection = analysis._effective_polymer_selection(settings)
    return ContactsReplicateRunner(
        universe=universe,
        target_selector=MDAnalysisSelector(settings.protein_selection),
        query_selector=MDAnalysisSelector(polymer_selection),
        cutoff=settings.cutoff,
        grouping=build_contact_grouping(settings.grouping),
        analyzer_cls=contacts_facade.ParallelContactAnalyzer,
    )


def summarize_replicate(
    analysis: Any,
    ctx: ReplicateContext,
    replicate: int,
    runner: Any,
    window: Any,
) -> Any:
    """Attach framework metadata to the runner-produced contact result.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for selection helpers.
    ctx : ReplicateContext
        Framework-provided context.
    replicate : int
        Replicate ID.
    runner : Any
        Runner-backed contacts execution object.
    window : Any
        Resolved trajectory window.

    Returns
    -------
    Any
        Contact result with framework metadata attached.
    """

    from polyzymd.analyses._results_base import get_polyzymd_version
    from polyzymd.analyses.shared.config_hash import compute_config_hash, settings_fingerprint
    from polyzymd.analyses.shared.loader import parse_time_string

    result = runner.results
    eq_value, eq_unit = parse_time_string(ctx.equilibration)

    result.config_hash = compute_config_hash(ctx.sim_config)
    result.polyzymd_version = get_polyzymd_version()
    result.replicate = replicate
    result.equilibration_time = eq_value
    result.equilibration_unit = eq_unit
    polymer_selection = analysis._effective_polymer_selection(ctx.settings)
    result.selection_string = f"{ctx.settings.protein_selection} : {polymer_selection}"
    result.start_frame = window.start
    result.timestep_ps = float(window.timestep_ps) * window.step

    metadata = dict(getattr(result, "metadata", {}) or {})
    contact_settings_fp = contacts_settings_fingerprint(ctx.settings)
    compute_residence_times = bool(ctx.settings.compute_residence_times)
    metadata["settings_fingerprint"] = contact_settings_fp
    metadata["contacts_settings_fingerprint"] = contact_settings_fp
    metadata["legacy_full_settings_fingerprint"] = settings_fingerprint(ctx.settings)
    metadata["settings_fingerprint_domain"] = CONTACTS_SETTINGS_FINGERPRINT_DOMAIN
    metadata["compute_residence_times"] = compute_residence_times
    metadata["residence_times_computed"] = compute_residence_times
    metadata["protein_selection"] = ctx.settings.protein_selection
    metadata["polymer_selection"] = ctx.settings.polymer_selection
    metadata["effective_polymer_selection"] = polymer_selection
    metadata["polymer_types_filter"] = normalize_polymer_types(ctx.settings.polymer_types)
    metadata["grouping"] = ctx.settings.grouping
    metadata["cutoff"] = float(ctx.settings.cutoff)
    metadata["criteria_cutoff"] = float(ctx.settings.cutoff)
    result.metadata = metadata

    has_stats = getattr(result, "has_per_residue_statistics", None)
    compute_stats = getattr(result, "compute_per_residue_statistics", None)
    if callable(has_stats) and callable(compute_stats) and not has_stats():
        compute_stats()

    return result


def aggregate(analysis: Any, ctx: AggregateContext, results: Sequence[Any]) -> Any:
    """Aggregate contact results across replicates for one condition.

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

    from polyzymd.analyses._results_base import get_polyzymd_version
    from polyzymd.analyses.contacts._aggregator import aggregate_contact_results
    from polyzymd.analyses.shared.config_hash import compute_config_hash, settings_fingerprint
    from polyzymd.analyses.shared.loader import parse_time_string

    agg_file = analysis._aggregated_sidecar_path(
        ctx.output_dir,
        ctx.settings,
        ctx.equilibration,
        ctx.replicates,
    )
    canonical_file = ctx.result_path or ctx.output_dir / "result.json"
    recompute = getattr(ctx, "recompute", False)
    aggregate_candidates = ()
    if not recompute:
        aggregate_candidates = analysis._aggregated_cache_candidates(
            ctx.output_dir,
            ctx.settings,
            ctx.equilibration,
            ctx.replicates,
            result_path=ctx.result_path,
        )

    for candidate in aggregate_candidates:
        cached = analysis._load_cache_candidate(
            AggregatedContactResult,
            candidate,
            recompute=recompute,
            sim_config=ctx.condition.sim_config,
            settings=None,
        )
        if cached is not None and analysis._cache_matches_context(
            cached,
            settings=ctx.settings,
            equilibration=ctx.equilibration,
            sim_config=ctx.condition.sim_config,
            replicates=ctx.replicates,
            source=candidate,
        ):
            if (
                ctx.result_path is not None
                and candidate != ctx.result_path
                and not ctx.result_path.exists()
            ):
                analysis._attach_contacts_identity_metadata(cached, ctx.settings)
                analysis.save_result(cached, ctx.result_path)
            if candidate == canonical_file and not agg_file.exists():
                analysis._attach_contacts_identity_metadata(cached, ctx.settings)
                cached.save(agg_file)
            return cached

    aligned_results = analysis._align_replicate_results(results, ctx.replicates)
    logger.info(f"  Aggregating {len(aligned_results)} replicates...")
    agg_result = aggregate_contact_results(
        aligned_results,
        compute_residence_times=ctx.settings.compute_residence_times,
    )

    eq_value, eq_unit = parse_time_string(ctx.equilibration)
    agg_result.config_hash = compute_config_hash(ctx.condition.sim_config)
    agg_result.polyzymd_version = get_polyzymd_version()
    agg_result.replicates = list(ctx.replicates)
    agg_result.equilibration_time = eq_value
    agg_result.equilibration_unit = eq_unit
    polymer_selection = analysis._effective_polymer_selection(ctx.settings)
    agg_result.selection_string = f"{ctx.settings.protein_selection} : {polymer_selection}"
    metadata = dict(getattr(agg_result, "metadata", {}) or {})
    contact_settings_fp = contacts_settings_fingerprint(ctx.settings)
    compute_residence_times = bool(ctx.settings.compute_residence_times)
    metadata["settings_fingerprint"] = contact_settings_fp
    metadata["contacts_settings_fingerprint"] = contact_settings_fp
    metadata["legacy_full_settings_fingerprint"] = settings_fingerprint(ctx.settings)
    metadata["settings_fingerprint_domain"] = CONTACTS_SETTINGS_FINGERPRINT_DOMAIN
    metadata["compute_residence_times"] = compute_residence_times
    metadata["residence_times_computed"] = compute_residence_times
    metadata["replicates"] = list(ctx.replicates)
    metadata["protein_selection"] = ctx.settings.protein_selection
    metadata["polymer_selection"] = ctx.settings.polymer_selection
    metadata["effective_polymer_selection"] = polymer_selection
    metadata["polymer_types_filter"] = normalize_polymer_types(ctx.settings.polymer_types)
    metadata["grouping"] = ctx.settings.grouping
    metadata["cutoff"] = float(ctx.settings.cutoff)
    metadata["criteria_cutoff"] = float(ctx.settings.cutoff)
    agg_result.metadata = metadata

    # Save both sidecar and canonical outputs for downstream consumers
    ctx.output_dir.mkdir(parents=True, exist_ok=True)
    agg_result.save(agg_file)
    if ctx.result_path is not None:
        analysis.save_result(agg_result, ctx.result_path)
    logger.info(f"  Saved aggregated contacts: {agg_file}")

    return agg_result
