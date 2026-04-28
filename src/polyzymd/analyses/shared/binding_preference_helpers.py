"""Shared helpers for computing binding preference from contacts data.

Used by the contacts, binding_free_energy, and polymer_affinity analysis
plugins to compute SASA-based binding preference enrichment.

Public functions
----------------
- find_enzyme_pdb()
- resolve_enzyme_pdb()
- try_load_cached_binding_preference()
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol, Sequence

if TYPE_CHECKING:
    from polyzymd.analyses.shared.binding_preference import (
        AggregatedBindingPreferenceResult,
        BindingPreferenceResult,
    )

logger = logging.getLogger("polyzymd.analyses")

BINDING_PREFERENCE_SETTINGS_METADATA_KEY = "binding_preference_settings_fingerprint"
CONTACTS_SETTINGS_METADATA_KEY = "contacts_settings_fingerprint"


class ConditionLike(Protocol):
    """Minimal condition protocol required by BP helper functions.

    Attributes
    ----------
    label : str
        Condition label
    replicates : Sequence[int]
        Replicate IDs associated with this condition
    """

    label: str
    replicates: Sequence[int]


def find_enzyme_pdb(sim_config: Any) -> Path | None:
    """Find enzyme PDB file from simulation config.

    Searches common locations relative to the project directory.

    Parameters
    ----------
    sim_config : SimulationConfig
        Simulation configuration (must have ``output.projects_directory``).

    Returns
    -------
    Path or None
        Path to enzyme PDB, or None if not found.

    Raises
    ------
    ValueError
        If a glob pattern matches multiple candidate enzyme PDB files.
    """
    import glob as glob_module

    project_dir = sim_config.output.projects_directory
    possible_paths = [
        project_dir / "structures" / "enzyme.pdb",
        project_dir / "input" / "enzyme.pdb",
        project_dir.parent / "structures" / "enzyme.pdb",
        project_dir.parent / "enzyme.pdb",
    ]

    for path in possible_paths:
        if path.exists():
            return path

    # Try glob for any PDB with "enzyme" in name
    patterns = [
        str(project_dir / "**" / "*enzyme*.pdb"),
        str(project_dir.parent / "*enzyme*.pdb"),
    ]
    for pattern in patterns:
        matches = sorted(glob_module.glob(pattern, recursive=True))
        if matches:
            if len(matches) > 1:
                raise ValueError(
                    f"Ambiguous enzyme PDB auto-discovery for pattern '{pattern}': "
                    f"{len(matches)} matches found: " + ", ".join(matches)
                )
            return Path(matches[0])

    return None


def resolve_enzyme_pdb(
    enzyme_pdb_setting: str | None,
    source_path: Path | None,
    sim_config: Any,
) -> Path | None:
    """Resolve the enzyme PDB path from settings or auto-discovery.

    Parameters
    ----------
    enzyme_pdb_setting : str or None
        Explicit enzyme PDB path from analysis settings (e.g.,
        ``enzyme_pdb_for_sasa``). If relative, resolved against
        *source_path*'s parent directory.
    source_path : Path or None
        Path to the comparison.yaml file (used to resolve relative paths).
    sim_config : Any
        Simulation configuration for auto-discovery fallback.

    Returns
    -------
    Path or None
        Resolved enzyme PDB path, or None if not found.
    """
    if enzyme_pdb_setting:
        if source_path:
            enzyme_pdb = source_path.parent / enzyme_pdb_setting
        else:
            enzyme_pdb = Path(enzyme_pdb_setting)
        if enzyme_pdb.exists():
            return enzyme_pdb
        logger.warning(f"Explicit enzyme_pdb_for_sasa not found at {enzyme_pdb}")
        return None

    return find_enzyme_pdb(sim_config)


def binding_preference_settings_fingerprint(settings: Any) -> str:
    """Compute the binding-preference settings fingerprint for a settings object.

    Parameters
    ----------
    settings : Any
        Analysis settings object containing binding-preference fields.

    Returns
    -------
    str
        First eight hexadecimal characters of the canonical settings digest.
    """

    return compute_binding_preference_settings_fingerprint(
        surface_exposure_threshold=getattr(settings, "surface_exposure_threshold", 0.2),
        enzyme_pdb_for_sasa=getattr(settings, "enzyme_pdb_for_sasa", None),
        include_default_aa_groups=getattr(settings, "include_default_aa_groups", True),
        protein_groups=getattr(settings, "protein_groups", None),
        protein_partitions=getattr(settings, "protein_partitions", None),
        polymer_type_selections=getattr(settings, "polymer_type_selections", None),
        polymer_chain=getattr(settings, "polymer_chain", "C"),
    )


def compute_binding_preference_settings_fingerprint(
    *,
    surface_exposure_threshold: float,
    enzyme_pdb_for_sasa: str | None = None,
    include_default_aa_groups: bool,
    protein_groups: dict[str, list[int]] | None,
    protein_partitions: dict[str, list[str]] | None,
    polymer_type_selections: dict[str, str] | None,
    polymer_chain: str,
) -> str:
    """Compute a fingerprint from parameters that define BP outputs.

    Parameters
    ----------
    surface_exposure_threshold : float
        Relative SASA threshold used for surface filtering.
    enzyme_pdb_for_sasa : str or None, optional
        Configured enzyme PDB path, when explicitly provided.
    include_default_aa_groups : bool
        Whether default amino-acid groups are included.
    protein_groups : dict[str, list[int]] or None
        Custom protein groups.
    protein_partitions : dict[str, list[str]] or None
        Partition definitions for mutually exclusive group comparisons.
    polymer_type_selections : dict[str, str] or None
        Polymer type selection map.
    polymer_chain : str
        Polymer chain ID used for auto-detection.

    Returns
    -------
    str
        First eight hexadecimal characters of the canonical settings digest.
    """

    payload = {
        "surface_exposure_threshold": float(surface_exposure_threshold),
        "enzyme_pdb_for_sasa": enzyme_pdb_for_sasa,
        "include_default_aa_groups": bool(include_default_aa_groups),
        "protein_groups": protein_groups,
        "protein_partitions": protein_partitions,
        "polymer_type_selections": polymer_type_selections,
        "polymer_chain": polymer_chain,
    }
    serialized = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()[:8]


def try_load_cached_binding_preference(
    cond: ConditionLike,
    analysis_dir: Path,
    *,
    settings_fp: str | None = None,
    contact_settings_fp: str | None = None,
    equilibration: str | None = None,
    successful_replicates: Sequence[int] | None = None,
) -> "AggregatedBindingPreferenceResult | BindingPreferenceResult | None":
    """Try to load cached binding preference results for a condition.

    Searches for binding preference files in order of preference:

    1. binding_preference_aggregated.json
    2. binding_preference_aggregated_reps*.json (glob pattern)
    3. binding_preference.json (single replicate)
    4. Per-replicate files (binding_preference_rep{N}.json)

    Parameters
    ----------
    cond : ConditionLike
        Condition to load.
    analysis_dir : Path
        Analysis directory for this condition.
    settings_fp : str or None, optional
        Binding-preference settings fingerprint for cache lookup. When
        provided, fingerprinted cache files are searched first, and legacy
        filenames are used only when metadata proves compatibility.
    contact_settings_fp : str or None, optional
        Contacts settings fingerprint that upstream contacts artifacts must
        have used. When provided, binding-preference caches must record this
        contact identity in metadata.
    equilibration : str or None, optional
        Requested contacts window. When provided, cached binding-preference
        artifacts must record the same window in metadata.
    successful_replicates : Sequence[int] or None, optional
        Replicate IDs that should be represented in an aggregate cache. When
        omitted, ``cond.replicates`` is used as the requested replicate
        contract.

    Returns
    -------
    AggregatedBindingPreferenceResult | BindingPreferenceResult | None
        Loaded result, or None if not found.
    """
    import glob as glob_module

    from polyzymd.analyses.shared.binding_preference import (
        AggregatedBindingPreferenceResult,
        BindingPreferenceResult,
        aggregate_binding_preference,
    )

    expected_replicates = _expected_binding_preference_replicates(cond, successful_replicates)

    if settings_fp is not None:
        fp_agg_path = analysis_dir / f"binding_preference_aggregated_s{settings_fp}.json"
        if fp_agg_path.exists():
            result = AggregatedBindingPreferenceResult.load(fp_agg_path)
            if _binding_preference_aggregate_cache_is_compatible(
                result,
                fp_agg_path,
                settings_fp=None,
                contact_settings_fp=contact_settings_fp,
                equilibration=equilibration,
                expected_replicates=expected_replicates,
                condition_label=cond.label,
            ):
                logger.debug(f"Loaded aggregated binding preference for {cond.label}")
                return result

        fp_agg_pattern = str(
            analysis_dir / f"binding_preference_aggregated_s{settings_fp}_reps*.json"
        )
        fp_agg_matches = sorted(glob_module.glob(fp_agg_pattern))
        if len(fp_agg_matches) == 1:
            result = AggregatedBindingPreferenceResult.load(fp_agg_matches[0])
            if _binding_preference_aggregate_cache_is_compatible(
                result,
                Path(fp_agg_matches[0]),
                settings_fp=None,
                contact_settings_fp=contact_settings_fp,
                equilibration=equilibration,
                expected_replicates=expected_replicates,
                condition_label=cond.label,
            ):
                logger.debug(f"Loaded aggregated binding preference for {cond.label}")
                return result
        if len(fp_agg_matches) > 1:
            raise ValueError(
                f"Ambiguous binding preference cache for {cond.label}: "
                f"found {len(fp_agg_matches)} files matching '{fp_agg_pattern}': "
                + ", ".join(fp_agg_matches)
            )

        fp_single_path = analysis_dir / f"binding_preference_s{settings_fp}.json"
        if fp_single_path.exists():
            result = BindingPreferenceResult.load(fp_single_path)
            if _binding_preference_cache_matches_contact_settings(
                result,
                fp_single_path,
                contact_settings_fp,
            ) and _binding_preference_cache_matches_window(result, fp_single_path, equilibration):
                logger.debug(f"Loaded single binding preference for {cond.label}")
                return result

        fp_rep_results = _collect_compatible_binding_preference_replicates(
            cond,
            path_for_replicate=lambda rep: analysis_dir
            / f"binding_preference_s{settings_fp}_rep{rep}.json",
            load_result=BindingPreferenceResult.load,
            settings_fp=settings_fp,
            contact_settings_fp=contact_settings_fp,
            equilibration=equilibration,
            expected_replicates=expected_replicates,
        )

        if fp_rep_results is None:
            return None
        if fp_rep_results:
            agg_result = aggregate_binding_preference(fp_rep_results)
            logger.debug(
                f"Aggregated {len(fp_rep_results)} replicate binding preference results "
                f"for {cond.label}"
            )
            return agg_result

    # Try aggregated result first (multi-replicate)
    agg_path = analysis_dir / "binding_preference_aggregated.json"
    if agg_path.exists():
        result = AggregatedBindingPreferenceResult.load(agg_path)
        if _binding_preference_aggregate_cache_is_compatible(
            result,
            agg_path,
            settings_fp=settings_fp,
            contact_settings_fp=contact_settings_fp,
            equilibration=equilibration,
            expected_replicates=expected_replicates,
            condition_label=cond.label,
        ):
            logger.debug(f"Loaded aggregated binding preference for {cond.label}")
            return result

    # Try aggregated result with rep range in name (e.g., _reps1-3.json)
    agg_pattern = str(analysis_dir / "binding_preference_aggregated_reps*.json")
    agg_matches = sorted(glob_module.glob(agg_pattern))
    if len(agg_matches) == 1:
        agg_match_path = Path(agg_matches[0])
        result = AggregatedBindingPreferenceResult.load(agg_match_path)
        if _binding_preference_aggregate_cache_is_compatible(
            result,
            agg_match_path,
            settings_fp=settings_fp,
            contact_settings_fp=contact_settings_fp,
            equilibration=equilibration,
            expected_replicates=expected_replicates,
            condition_label=cond.label,
        ):
            logger.debug(f"Loaded aggregated binding preference for {cond.label}")
            return result
    if len(agg_matches) > 1:
        raise ValueError(
            f"Ambiguous binding preference cache for {cond.label}: "
            f"found {len(agg_matches)} files: " + ", ".join(agg_matches)
        )

    # Try single replicate result
    single_path = analysis_dir / "binding_preference.json"
    if single_path.exists():
        result = BindingPreferenceResult.load(single_path)
        if not _binding_preference_cache_matches_settings(result, single_path, settings_fp):
            logger.warning(
                "Ignoring binding preference cache for %s without matching settings fingerprint: %s",
                cond.label,
                single_path,
            )
        elif not _binding_preference_cache_matches_contact_settings(
            result,
            single_path,
            contact_settings_fp,
        ):
            logger.warning(
                "Ignoring binding preference cache for %s without matching contacts settings "
                "fingerprint: %s",
                cond.label,
                single_path,
            )
        elif not _binding_preference_cache_matches_window(result, single_path, equilibration):
            logger.warning(
                "Ignoring binding preference cache for %s without matching window: %s",
                cond.label,
                single_path,
            )
        else:
            logger.debug(f"Loaded single binding preference for {cond.label}")
            return result

    # Try per-replicate results and aggregate them
    rep_results = _collect_compatible_binding_preference_replicates(
        cond,
        path_for_replicate=lambda rep: analysis_dir / f"binding_preference_rep{rep}.json",
        load_result=BindingPreferenceResult.load,
        settings_fp=settings_fp,
        contact_settings_fp=contact_settings_fp,
        equilibration=equilibration,
        expected_replicates=expected_replicates,
    )

    if rep_results:
        agg_result = aggregate_binding_preference(rep_results)
        logger.debug(
            f"Aggregated {len(rep_results)} replicate binding preference results for {cond.label}"
        )
        return agg_result

    return None


def _collect_compatible_binding_preference_replicates(
    cond: ConditionLike,
    *,
    path_for_replicate: Any,
    load_result: Any,
    settings_fp: str | None,
    contact_settings_fp: str | None,
    equilibration: str | None,
    expected_replicates: tuple[int, ...] | None,
) -> list[Any] | None:
    """Load per-replicate binding-preference caches only for an exact identity.

    Parameters
    ----------
    cond : ConditionLike
        Condition whose replicate IDs define candidate file locations.
    path_for_replicate : callable
        Function returning the expected path for a replicate ID.
    load_result : callable
        Function used to deserialize a binding-preference result.
    settings_fp : str or None
        Required binding-preference settings fingerprint, if known.
    contact_settings_fp : str or None
        Required contacts settings fingerprint, if known.
    equilibration : str or None
        Requested contacts window, if known.
    expected_replicates : tuple[int, ...] or None
        Exact replicate identity required for aggregation.

    Returns
    -------
    list[Any] | None
        Loaded per-replicate results, an empty list when no cache files exist,
        or ``None`` when the cache files are incomplete or mismatched.
    """

    rep_results: list[Any] = []
    loaded_replicates: list[int] = []
    for rep in cond.replicates:
        rep_path = path_for_replicate(rep)
        if not rep_path.exists():
            continue
        result = load_result(rep_path)
        if not _binding_preference_cache_matches_settings(result, rep_path, settings_fp):
            logger.warning(
                "Ignoring binding preference cache for %s without matching settings "
                "fingerprint: %s",
                cond.label,
                rep_path,
            )
            continue
        if not _binding_preference_cache_matches_contact_settings(
            result,
            rep_path,
            contact_settings_fp,
        ):
            logger.warning(
                "Ignoring binding preference cache for %s without matching contacts settings "
                "fingerprint: %s",
                cond.label,
                rep_path,
            )
            continue
        if not _binding_preference_cache_matches_window(result, rep_path, equilibration):
            logger.warning(
                "Ignoring binding preference cache for %s without matching window: %s",
                cond.label,
                rep_path,
            )
            continue
        if not _binding_preference_replicate_cache_matches_identity(result, rep_path, rep):
            logger.warning(
                "Ignoring binding preference cache for %s with mismatched replicate metadata: %s",
                cond.label,
                rep_path,
            )
            return None
        rep_results.append(result)
        loaded_replicates.append(int(rep))

    if not rep_results:
        return []

    loaded_identity = _normalize_replicate_ids(loaded_replicates)
    if expected_replicates is not None and loaded_identity != expected_replicates:
        logger.debug(
            "Rejecting per-replicate binding preference caches for %s: loaded=%s expected=%s",
            cond.label,
            loaded_identity,
            expected_replicates,
        )
        return None
    return rep_results


def _binding_preference_replicate_cache_matches_identity(
    result: Any,
    path: Path,
    expected_replicate: int,
) -> bool:
    """Return whether per-replicate metadata matches the filename identity.

    Parameters
    ----------
    result : Any
        Loaded per-replicate binding-preference result.
    path : Path
        Cache path used in diagnostics.
    expected_replicate : int
        Replicate ID implied by the lookup path.

    Returns
    -------
    bool
        ``True`` when metadata is absent or declares the expected replicate ID.
    """

    metadata = getattr(result, "metadata", None)
    if not isinstance(metadata, dict) or "replicate" not in metadata:
        return True
    try:
        stored_replicate = int(metadata["replicate"])
    except (TypeError, ValueError):
        logger.debug("Binding preference replicate cache has invalid metadata: %s", path)
        return False
    if stored_replicate != int(expected_replicate):
        logger.debug(
            "Binding preference replicate cache mismatch for %s: stored=%d expected=%d",
            path,
            stored_replicate,
            expected_replicate,
        )
        return False
    return True


def _binding_preference_cache_matches_settings(
    result: Any,
    path: Path,
    settings_fp: str | None,
) -> bool:
    """Return whether a binding-preference cache proves settings identity.

    Parameters
    ----------
    result : Any
        Loaded binding-preference result.
    path : Path
        Cache path used for filename fingerprint inspection.
    settings_fp : str or None
        Active settings fingerprint, if known.

    Returns
    -------
    bool
        ``True`` when no fingerprint is required, or when the path or result
        metadata declares the active fingerprint.
    """

    if settings_fp is None:
        return True

    from polyzymd.analyses.shared.config_hash import extract_settings_fingerprint_from_path

    path_fp = extract_settings_fingerprint_from_path(path)
    if path_fp == settings_fp:
        return True

    result_fp = getattr(result, BINDING_PREFERENCE_SETTINGS_METADATA_KEY, None) or getattr(
        result,
        "settings_fingerprint",
        None,
    )
    if result_fp is None:
        result_fp = getattr(result, "settings_fp", None)
    metadata = getattr(result, "metadata", None)
    if result_fp is None and isinstance(metadata, dict):
        result_fp = (
            metadata.get(BINDING_PREFERENCE_SETTINGS_METADATA_KEY)
            or metadata.get("settings_fingerprint")
            or metadata.get("settings_fp")
        )
    if result_fp is None:
        result_fp = _read_binding_preference_settings_fingerprint(path)

    return result_fp == settings_fp


def _binding_preference_cache_matches_contact_settings(
    result: Any,
    path: Path,
    contact_settings_fp: str | None,
) -> bool:
    """Return whether a BP cache records the required contacts identity.

    Parameters
    ----------
    result : Any
        Loaded binding-preference result.
    path : Path
        Cache path used for JSON metadata fallback.
    contact_settings_fp : str or None
        Required contacts settings fingerprint, if known.

    Returns
    -------
    bool
        ``True`` when no contacts fingerprint is required, or when metadata
        records the active upstream contacts fingerprint.
    """

    if contact_settings_fp is None:
        return True

    metadata = getattr(result, "metadata", None)
    result_fp = None
    if isinstance(metadata, dict):
        result_fp = metadata.get(CONTACTS_SETTINGS_METADATA_KEY) or metadata.get(
            "contact_settings_fingerprint"
        )
    if result_fp is None:
        result_fp = _read_binding_preference_contact_settings_fingerprint(path)

    return result_fp == contact_settings_fp


def _binding_preference_aggregate_cache_is_compatible(
    result: Any,
    path: Path,
    *,
    settings_fp: str | None,
    contact_settings_fp: str | None,
    equilibration: str | None,
    expected_replicates: tuple[int, ...] | None,
    condition_label: str | None = None,
) -> bool:
    """Return whether an aggregate cache matches all active cache contracts.

    Parameters
    ----------
    result : Any
        Loaded aggregate binding-preference result.
    path : Path
        Cache path used for metadata and filename inspection.
    settings_fp : str or None
        Required binding-preference settings fingerprint, if known.
    contact_settings_fp : str or None
        Required contacts settings fingerprint, if known.
    equilibration : str or None
        Requested equilibration window, if known.
    expected_replicates : tuple[int, ...] or None
        Replicate IDs that the aggregate must prove it represents.
    condition_label : str or None, optional
        Condition label used in warning messages.

    Returns
    -------
    bool
        ``True`` when the aggregate is safe to reuse.
    """

    label = condition_label or "condition"
    if not _binding_preference_cache_matches_settings(result, path, settings_fp):
        logger.warning(
            "Ignoring binding preference cache for %s without matching settings fingerprint: %s",
            label,
            path,
        )
        return False
    if not _binding_preference_cache_matches_contact_settings(result, path, contact_settings_fp):
        logger.warning(
            "Ignoring binding preference cache for %s without matching contacts settings "
            "fingerprint: %s",
            label,
            path,
        )
        return False
    if not _binding_preference_cache_matches_window(result, path, equilibration):
        logger.warning(
            "Ignoring binding preference cache for %s without matching window: %s",
            label,
            path,
        )
        return False
    if not _binding_preference_aggregate_matches_replicates(result, path, expected_replicates):
        logger.warning(
            "Ignoring binding preference aggregate cache for %s without matching replicates: %s",
            label,
            path,
        )
        return False
    return True


def _binding_preference_aggregate_matches_replicates(
    result: Any,
    path: Path,
    expected_replicates: tuple[int, ...] | None,
) -> bool:
    """Return whether an aggregate cache proves replicate identity.

    Parameters
    ----------
    result : Any
        Loaded aggregate binding-preference result.
    path : Path
        Cache path used in diagnostics.
    expected_replicates : tuple[int, ...] or None
        Requested or successful replicate contract. ``None`` means no
        replicate-aware contract is available.

    Returns
    -------
    bool
        ``True`` when the aggregate metadata and replicate count match the
        expected replicate IDs.
    """

    metadata = getattr(result, "metadata", None)
    stored_replicates = None
    if isinstance(metadata, dict):
        stored_replicates = _normalize_replicate_ids(metadata.get("replicates"))

    if expected_replicates is None:
        if stored_replicates is None:
            return True
        return _binding_preference_n_replicates_matches(result, path, len(stored_replicates))

    if stored_replicates is None:
        logger.debug("Binding preference aggregate lacks replicate metadata: %s", path)
        return False
    if stored_replicates != expected_replicates:
        logger.debug(
            "Binding preference aggregate replicate mismatch for %s: stored=%s expected=%s",
            path,
            stored_replicates,
            expected_replicates,
        )
        return False
    return _binding_preference_n_replicates_matches(result, path, len(expected_replicates))


def _binding_preference_n_replicates_matches(result: Any, path: Path, expected_count: int) -> bool:
    """Return whether ``n_replicates`` matches an expected aggregate count.

    Parameters
    ----------
    result : Any
        Loaded aggregate binding-preference result.
    path : Path
        Cache path used in diagnostics.
    expected_count : int
        Expected number of replicate IDs represented by the aggregate.

    Returns
    -------
    bool
        ``True`` when the top-level replicate count matches.
    """

    try:
        n_replicates = int(getattr(result, "n_replicates", None))
    except (TypeError, ValueError):
        logger.debug("Binding preference aggregate has invalid n_replicates: %s", path)
        return False
    if n_replicates != expected_count:
        logger.debug(
            "Binding preference aggregate n_replicates mismatch for %s: stored=%d expected=%d",
            path,
            n_replicates,
            expected_count,
        )
        return False
    return True


def _expected_binding_preference_replicates(
    cond: ConditionLike,
    successful_replicates: Sequence[int] | None,
) -> tuple[int, ...] | None:
    """Return the replicate contract expected for an aggregate cache.

    Parameters
    ----------
    cond : ConditionLike
        Condition supplying requested replicate IDs.
    successful_replicates : Sequence[int] or None
        Explicit successful replicate IDs from resolved upstream artifacts.

    Returns
    -------
    tuple[int, ...] or None
        Normalized replicate IDs, or ``None`` when no contract is available.
    """

    if successful_replicates is not None:
        return _normalize_replicate_ids(successful_replicates) or ()
    return _normalize_replicate_ids(getattr(cond, "replicates", None))


def _normalize_replicate_ids(value: Any) -> tuple[int, ...] | None:
    """Normalize replicate IDs for order-independent cache comparisons.

    Parameters
    ----------
    value : Any
        Replicate ID collection from metadata or a condition.

    Returns
    -------
    tuple[int, ...] or None
        Sorted unique replicate IDs, or ``None`` when the value is invalid or
        missing.
    """

    if value is None or isinstance(value, str):
        return None
    try:
        replicates = tuple(sorted(int(rep) for rep in value))
    except (TypeError, ValueError):
        return None
    if len(set(replicates)) != len(replicates):
        return None
    return replicates


def _binding_preference_cache_matches_window(
    result: Any,
    path: Path,
    equilibration: str | None,
) -> bool:
    """Return whether a binding-preference cache matches the requested window.

    Parameters
    ----------
    result : Any
        Loaded binding-preference result.
    path : Path
        Cache path used for JSON metadata fallback.
    equilibration : str or None
        Requested equilibration window.

    Returns
    -------
    bool
        ``True`` when no window is requested, or metadata proves the same
        window.
    """

    if equilibration is None:
        return True

    from polyzymd.analyses.contacts._paths import binding_preference_artifact_matches_window

    metadata = getattr(result, "metadata", None)
    if isinstance(metadata, dict) and binding_preference_artifact_matches_window(
        path,
        equilibration,
    ):
        return True
    return binding_preference_artifact_matches_window(path, equilibration)


def _read_binding_preference_settings_fingerprint(path: Path) -> str | None:
    """Read a settings fingerprint directly from a cache JSON file.

    Parameters
    ----------
    path : Path
        Binding-preference cache path.

    Returns
    -------
    str or None
        Settings fingerprint from top-level or metadata keys, if present.
    """

    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    direct = (
        payload.get(BINDING_PREFERENCE_SETTINGS_METADATA_KEY)
        or payload.get("settings_fingerprint")
        or payload.get("settings_fp")
    )
    if direct is not None:
        return str(direct)
    metadata = payload.get("metadata")
    if isinstance(metadata, dict):
        stored = (
            metadata.get(BINDING_PREFERENCE_SETTINGS_METADATA_KEY)
            or metadata.get("settings_fingerprint")
            or metadata.get("settings_fp")
        )
        if stored is not None:
            return str(stored)
    return None


def _read_binding_preference_contact_settings_fingerprint(path: Path) -> str | None:
    """Read the upstream contacts settings fingerprint from a BP cache file.

    Parameters
    ----------
    path : Path
        Binding-preference cache path.

    Returns
    -------
    str or None
        Contacts settings fingerprint from top-level or metadata keys, if
        present.
    """

    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    direct = payload.get(CONTACTS_SETTINGS_METADATA_KEY) or payload.get(
        "contact_settings_fingerprint"
    )
    if direct is not None:
        return str(direct)
    metadata = payload.get("metadata")
    if isinstance(metadata, dict):
        stored = metadata.get(CONTACTS_SETTINGS_METADATA_KEY) or metadata.get(
            "contact_settings_fingerprint"
        )
        if stored is not None:
            return str(stored)
    return None
