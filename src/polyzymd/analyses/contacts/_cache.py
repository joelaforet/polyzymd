"""Cache helpers for the contacts analysis facade."""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Sequence

from pydantic import BaseModel

from polyzymd.analyses.contacts._aggregator import AggregatedContactResult
from polyzymd.analyses.contacts._identity import (
    CONTACTS_SETTINGS_FINGERPRINT_DOMAIN,
    contacts_settings_fingerprint,
    contacts_settings_fingerprint_candidates,
    normalize_polymer_types,
)
from polyzymd.analyses.contacts._identity import (
    effective_polymer_selection as build_effective_polymer_selection,
)

logger = logging.getLogger("polyzymd.analyses.contacts")


def replicate_sidecar_path(
    output_dir: Path,
    settings: BaseModel,
    equilibration: str,
    replicate: int,
) -> Path:
    """Build the legacy per-replicate contacts side-output path.

    Parameters
    ----------
    output_dir : Path
        Directory containing per-replicate contacts outputs.
    settings : BaseModel
        Active contacts settings.
    equilibration : str
        Requested analysis equilibration window.
    replicate : int
        Requested replicate ID.

    Returns
    -------
    Path
        Fingerprinted per-replicate sidecar path.
    """

    from polyzymd.analyses.shared.loader import parse_time_string

    eq_value, eq_unit = parse_time_string(equilibration)
    eq_str = f"eq{eq_value:g}{eq_unit}"
    cut_str = f"cut{settings.cutoff:.1f}"
    settings_fp = contacts_settings_fingerprint(settings)
    return output_dir / f"contacts_{eq_str}_{cut_str}_s{settings_fp}_rep{replicate}.json"


def replicate_sidecar_candidates(
    output_dir: Path,
    settings: BaseModel,
    equilibration: str,
    replicate: int,
) -> tuple[Path, ...]:
    """Build accepted per-replicate contacts side-output paths.

    Parameters
    ----------
    output_dir : Path
        Directory containing per-replicate contacts outputs.
    settings : BaseModel
        Active contacts settings.
    equilibration : str
        Requested analysis equilibration window.
    replicate : int
        Requested replicate ID.

    Returns
    -------
    tuple[Path, ...]
        Current sidecar path followed by legacy compatibility candidates.
    """

    from polyzymd.analyses.shared.loader import parse_time_string

    eq_value, eq_unit = parse_time_string(equilibration)
    eq_str = f"eq{eq_value:g}{eq_unit}"
    cut_str = f"cut{settings.cutoff:.1f}"
    return tuple(
        output_dir / f"contacts_{eq_str}_{cut_str}_s{settings_fp}_rep{replicate}.json"
        for settings_fp in contacts_settings_fingerprint_candidates(settings)
    )


def replicate_cache_candidates(
    output_dir: Path,
    settings: BaseModel,
    equilibration: str,
    replicate: int,
    *,
    result_path: Path | None = None,
) -> tuple[Path, ...]:
    """Build validated lookup candidates for a per-replicate contacts cache.

    Parameters
    ----------
    output_dir : Path
        Per-replicate output directory, usually ``run_<N>``.
    settings : BaseModel
        Active contacts settings.
    equilibration : str
        Requested analysis equilibration window.
    replicate : int
        Requested replicate ID.
    result_path : Path or None, optional
        Canonical cache path supplied by the orchestrator.

    Returns
    -------
    tuple[Path, ...]
        Sidecar candidates followed by the canonical result path. All returned
        paths are still subject to content validation by the lifecycle.
    """

    analysis_dir = _analysis_dir_for_replicate_output(output_dir, replicate)
    sidecars = [*replicate_sidecar_candidates(output_dir, settings, equilibration, replicate)]
    sidecars.extend(
        _replicate_sidecar_discovery_candidates(
            analysis_dir,
            output_dir,
            settings,
            equilibration,
            replicate,
        )
    )

    canonical = result_path or output_dir / "result.json"
    return _dedupe_paths([*sidecars, canonical])


def _analysis_dir_for_replicate_output(output_dir: Path, replicate: int) -> Path:
    """Return the contacts analysis directory for a replicate output path.

    Parameters
    ----------
    output_dir : Path
        Per-replicate output directory.
    replicate : int
        Requested replicate ID.

    Returns
    -------
    Path
        Parent contacts analysis directory when ``output_dir`` is ``run_<N>``,
        otherwise ``output_dir``.
    """

    if output_dir.name == f"run_{int(replicate)}":
        return output_dir.parent
    return output_dir


def _replicate_sidecar_discovery_candidates(
    analysis_dir: Path,
    output_dir: Path,
    settings: BaseModel,
    equilibration: str,
    replicate: int,
) -> tuple[Path, ...]:
    """Return sidecar candidates mirrored from contacts path discovery.

    Parameters
    ----------
    analysis_dir : Path
        Contacts analysis directory for the condition.
    output_dir : Path
        Per-replicate output directory.
    settings : BaseModel
        Active contacts settings.
    equilibration : str
        Requested analysis equilibration window.
    replicate : int
        Requested replicate ID.

    Returns
    -------
    tuple[Path, ...]
        Fingerprinted and legacy sidecar paths in deterministic precedence
        order, excluding canonical ``result.json``.
    """

    from polyzymd.analyses.contacts._paths import contact_equilibration_filename_token

    eq_token = contact_equilibration_filename_token(equilibration)
    eq_glob = eq_token if eq_token is not None else "eq*"
    run_dir = analysis_dir / f"run_{int(replicate)}"
    candidate_paths: list[Path] = []

    for settings_fp in contacts_settings_fingerprint_candidates(settings):
        candidate_paths.extend(
            [
                analysis_dir / f"contacts_{eq_glob}_cut*_s{settings_fp}_rep{replicate}.json",
                output_dir / f"contacts_{eq_glob}_cut*_s{settings_fp}_rep{replicate}.json",
                run_dir / f"contacts_{eq_glob}_cut*_s{settings_fp}_rep{replicate}.json",
            ]
        )

    candidate_paths.extend(
        [
            analysis_dir / f"contacts_{eq_glob}_cut*_rep{replicate}.json",
            output_dir / f"contacts_{eq_glob}_cut*_rep{replicate}.json",
            run_dir / f"contacts_{eq_glob}_cut*_rep{replicate}.json",
            analysis_dir / f"contacts_rep{replicate}.json",
            output_dir / f"contacts_rep{replicate}.json",
            run_dir / f"contacts_rep{replicate}.json",
        ]
    )
    return _expand_glob_candidates(candidate_paths)


def aggregated_sidecar_path(
    output_dir: Path,
    settings: BaseModel,
    equilibration: str,
    replicates: Sequence[int],
) -> Path:
    """Build the legacy aggregated contacts side-output path.

    Parameters
    ----------
    output_dir : Path
        Directory containing aggregated contacts outputs.
    settings : BaseModel
        Active contacts settings.
    equilibration : str
        Requested analysis equilibration window.
    replicates : Sequence[int]
        Requested replicate IDs.

    Returns
    -------
    Path
        Fingerprinted aggregated sidecar path.
    """

    from polyzymd.analyses.shared.loader import parse_time_string
    from polyzymd.analyses.shared.paths import format_replicate_cache_token

    rep_token = format_replicate_cache_token(replicates)
    eq_value, eq_unit = parse_time_string(equilibration)
    eq_str = f"eq{eq_value:g}{eq_unit}"
    cut_str = f"cut{settings.cutoff:.1f}"
    settings_fp = contacts_settings_fingerprint(settings)
    return output_dir / f"contacts_aggregated_{eq_str}_{cut_str}_s{settings_fp}_{rep_token}.json"


def aggregated_sidecar_candidates(
    output_dir: Path,
    settings: BaseModel,
    equilibration: str,
    replicates: Sequence[int],
) -> tuple[Path, ...]:
    """Build accepted aggregated contacts side-output paths.

    Parameters
    ----------
    output_dir : Path
        Directory containing aggregated contacts outputs.
    settings : BaseModel
        Active contacts settings.
    equilibration : str
        Requested analysis equilibration window.
    replicates : Sequence[int]
        Requested replicate IDs.

    Returns
    -------
    tuple[Path, ...]
        Current sidecar path followed by legacy compatibility candidates.
    """

    from polyzymd.analyses.shared.loader import parse_time_string
    from polyzymd.analyses.shared.paths import format_replicate_cache_token

    rep_token = format_replicate_cache_token(replicates)
    eq_value, eq_unit = parse_time_string(equilibration)
    eq_str = f"eq{eq_value:g}{eq_unit}"
    cut_str = f"cut{settings.cutoff:.1f}"
    return tuple(
        output_dir / f"contacts_aggregated_{eq_str}_{cut_str}_s{settings_fp}_{rep_token}.json"
        for settings_fp in contacts_settings_fingerprint_candidates(settings)
    )


def aggregated_cache_candidates(
    output_dir: Path,
    settings: BaseModel,
    equilibration: str,
    replicates: Sequence[int],
    *,
    result_path: Path | None = None,
) -> tuple[Path, ...]:
    """Build validated lookup candidates for an aggregated contacts cache.

    Parameters
    ----------
    output_dir : Path
        Directory containing aggregated contacts outputs.
    settings : BaseModel
        Active contacts settings.
    equilibration : str
        Requested analysis equilibration window.
    replicates : Sequence[int]
        Requested replicate IDs.
    result_path : Path or None, optional
        Canonical aggregate cache path supplied by the orchestrator.

    Returns
    -------
    tuple[Path, ...]
        Fingerprinted sidecars, legacy aggregate JSONs, and the canonical
        result path in validation order.
    """

    canonical = result_path or output_dir / "result.json"
    return _dedupe_paths(
        [
            *aggregated_sidecar_candidates(output_dir, settings, equilibration, replicates),
            *_legacy_aggregated_sidecar_candidates(output_dir, equilibration),
            canonical,
        ]
    )


def _legacy_aggregated_sidecar_candidates(output_dir: Path, equilibration: str) -> tuple[Path, ...]:
    """Return legacy aggregate sidecar JSON paths in deterministic order.

    Parameters
    ----------
    output_dir : Path
        Directory containing aggregated contacts outputs.
    equilibration : str
        Requested analysis equilibration window.

    Returns
    -------
    tuple[Path, ...]
        Non-canonical aggregate sidecars whose contents must still pass full
        context validation before use.
    """

    from polyzymd.analyses.contacts._paths import contact_equilibration_filename_token

    eq_token = contact_equilibration_filename_token(equilibration)
    eq_glob = eq_token if eq_token is not None else "eq*"
    patterns = (
        f"contacts_aggregated_{eq_glob}_cut*_reps*.json",
        f"contacts_aggregated_{eq_glob}_cut*.json",
        f"contacts_aggregated_{eq_glob}_*.json",
        f"contacts_aggregated_{eq_glob}.json",
        "contacts_aggregated_reps*.json",
        "contacts_aggregated.json",
        "aggregated_contacts*.json",
    )
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(sorted(output_dir.glob(pattern)))
    paths.extend(sorted(output_dir.glob("*.json")))
    return tuple(path for path in _dedupe_paths(paths) if path.name != "result.json")


def _expand_glob_candidates(candidates: Sequence[Path]) -> tuple[Path, ...]:
    """Expand literal and glob-style path candidates.

    Parameters
    ----------
    candidates : Sequence[Path]
        Candidate paths that may contain glob metacharacters in the file name.

    Returns
    -------
    tuple[Path, ...]
        Existing paths in deterministic order.
    """

    paths: list[Path] = []
    for candidate in candidates:
        name = candidate.name
        if any(char in name for char in "*?["):
            paths.extend(sorted(candidate.parent.glob(name)))
        elif candidate.exists():
            paths.append(candidate)
    return _dedupe_paths(paths)


def _dedupe_paths(paths: Sequence[Path]) -> tuple[Path, ...]:
    """Return paths with duplicates removed while preserving order.

    Parameters
    ----------
    paths : Sequence[Path]
        Paths to deduplicate.

    Returns
    -------
    tuple[Path, ...]
        Deduplicated paths.
    """

    return tuple(dict.fromkeys(paths))


def load_cache_candidate(
    analysis: Any,
    result_cls: type,
    candidate: Path,
    *,
    recompute: bool,
    sim_config: Any | None = None,
    settings: BaseModel | None = None,
) -> Any | None:
    """Load one cache candidate without letting invalid legacy JSON abort lookup.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for the underlying cache hook.
    result_cls : type
        Result model class with a ``load()`` method.
    candidate : Path
        Candidate JSON path to load.
    recompute : bool
        Whether cache loading should be skipped.
    sim_config : Any or None, optional
        Simulation configuration used by the underlying cache hook.
    settings : BaseModel or None, optional
        Settings model used by the underlying cache hook.

    Returns
    -------
    Any or None
        Loaded cache result, or ``None`` when the candidate is absent, invalid,
        or incompatible with the supplied cache hook arguments.
    """

    try:
        return analysis._check_cache(
            result_cls,
            candidate,
            recompute=recompute,
            sim_config=sim_config,
            settings=settings,
        )
    except (OSError, ValueError) as exc:
        logger.warning("contacts: ignoring unreadable cache candidate at %s: %s", candidate, exc)
        return None


def effective_polymer_selection(settings: BaseModel) -> str:
    """Return the polymer selection constrained by polymer type filters.

    Parameters
    ----------
    settings : BaseModel
        Active contacts settings.

    Returns
    -------
    str
        MDAnalysis selection used for contact query atoms.
    """

    polymer_types = [
        str(polymer_type).strip()
        for polymer_type in (settings.polymer_types or [])
        if str(polymer_type).strip()
    ]
    return build_effective_polymer_selection(settings.polymer_selection, polymer_types)


def cache_matches_window(result: Any, equilibration: str) -> bool:
    """Return whether a cached result matches the requested analysis window.

    Parameters
    ----------
    result : Any
        Loaded contacts result.
    equilibration : str
        Requested equilibration window from the comparison context.

    Returns
    -------
    bool
        ``True`` when stored equilibration metadata matches the requested
        window, otherwise ``False``.
    """

    from polyzymd.analyses.shared.loader import convert_time, parse_time_string

    expected_time, expected_unit = parse_time_string(equilibration)
    stored_time = getattr(result, "equilibration_time", None)
    stored_unit = getattr(result, "equilibration_unit", None)
    if stored_time is None or stored_unit is None:
        return False
    try:
        stored_time_ps = convert_time(float(stored_time), str(stored_unit), "ps")
        expected_time_ps = convert_time(float(expected_time), str(expected_unit), "ps")
    except (TypeError, ValueError):
        return False
    return math.isclose(stored_time_ps, expected_time_ps, rel_tol=0.0, abs_tol=1.0e-9)


def result_window_string(result: Any) -> str | None:
    """Return the equilibration window encoded on a result.

    Parameters
    ----------
    result : Any
        Loaded contacts result with equilibration metadata.

    Returns
    -------
    str or None
        Window string such as ``"10ns"``, or ``None`` when metadata is missing
        or malformed.
    """

    stored_time = getattr(result, "equilibration_time", None)
    stored_unit = getattr(result, "equilibration_unit", None)
    if stored_time is None or stored_unit is None:
        return None
    try:
        stored_time_value = float(stored_time)
    except (TypeError, ValueError):
        return None
    return f"{stored_time_value:g}{stored_unit}"


def cache_matches_replicate_id(
    result: Any,
    replicate: int,
    *,
    source: Path | None = None,
) -> bool:
    """Return whether a cached per-replicate result matches the request.

    Parameters
    ----------
    result : Any
        Loaded per-replicate contacts result.
    replicate : int
        Requested replicate ID.
    source : Path or None, optional
        Cache source path used for diagnostics.

    Returns
    -------
    bool
        ``True`` when the stored replicate ID equals ``replicate``.
    """

    stored_replicate = getattr(result, "replicate", None)
    try:
        stored = int(stored_replicate)
    except (TypeError, ValueError):
        logger.warning(
            "contacts: ignoring replicate cache without valid replicate identity%s",
            f" at {source}" if source is not None else "",
        )
        return False
    if stored != int(replicate):
        logger.warning(
            "contacts: ignoring replicate cache for replicate %d; requested %d%s",
            stored,
            int(replicate),
            f" at {source}" if source is not None else "",
        )
        return False
    return True


def embedded_contacts_settings_fingerprint(result: Any) -> str | None:
    """Return contacts settings fingerprint embedded in cache content.

    Parameters
    ----------
    result : Any
        Loaded contacts artifact.

    Returns
    -------
    str or None
        Embedded fingerprint from result fields or metadata. Filename-only
        fingerprints are intentionally excluded.
    """

    stored_fingerprint = getattr(result, "settings_fingerprint", None)
    if stored_fingerprint is None:
        stored_fingerprint = getattr(result, "settings_fp", None)
    if stored_fingerprint is None:
        metadata = getattr(result, "metadata", None)
        if isinstance(metadata, dict):
            stored_fingerprint = (
                metadata.get("contacts_settings_fingerprint")
                or metadata.get("contact_settings_fingerprint")
                or metadata.get("settings_fingerprint")
                or metadata.get("settings_fp")
            )
    return str(stored_fingerprint) if stored_fingerprint is not None else None


def cache_has_contacts_identity_proof(result: Any, settings: BaseModel) -> bool:
    """Return whether cache content fully proves contacts settings identity.

    Parameters
    ----------
    result : Any
        Loaded contacts artifact.
    settings : BaseModel
        Active contacts settings.

    Returns
    -------
    bool
        ``True`` when content metadata/schema fields match the active cutoff,
        selections, and grouping without relying on the filename.
    """

    metadata = getattr(result, "metadata", None)
    metadata = metadata if isinstance(metadata, dict) else {}
    cutoff = getattr(result, "criteria_cutoff", None)
    if cutoff is None:
        cutoff = metadata.get("criteria_cutoff", metadata.get("cutoff"))
    try:
        cutoff_matches = abs(float(cutoff) - float(settings.cutoff)) <= 1e-6
    except (TypeError, ValueError):
        return False
    if not cutoff_matches:
        return False

    selection_string = str(getattr(result, "selection_string", "") or "")
    protein_selection = metadata.get("protein_selection")
    effective_selection = metadata.get("effective_polymer_selection")
    if (protein_selection is None or effective_selection is None) and " : " in selection_string:
        protein_text, polymer_text = selection_string.split(" : ", 1)
        protein_selection = protein_selection or protein_text
        effective_selection = effective_selection or polymer_text

    grouping = metadata.get("grouping", metadata.get("contact_grouping"))
    expected_effective = effective_polymer_selection(settings)
    return (
        str(protein_selection).strip() == str(settings.protein_selection).strip()
        and str(effective_selection).strip() == expected_effective.strip()
        and str(grouping).strip() == str(settings.grouping).strip()
    )


def coerce_optional_bool(value: Any) -> bool | None:
    """Coerce common metadata values to booleans.

    Parameters
    ----------
    value : Any
        Metadata value to normalize.

    Returns
    -------
    bool or None
        Boolean value, or ``None`` when the input is missing or ambiguous.
    """

    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1", "yes", "y"}:
            return True
        if normalized in {"false", "0", "no", "n"}:
            return False
    if isinstance(value, int):
        if value == 1:
            return True
        if value == 0:
            return False
    return None


def cache_matches_residence_time_setting(
    result: Any,
    settings: BaseModel,
    *,
    allow_missing: bool,
    source: Path | None = None,
) -> bool:
    """Return whether a cache proves the active residence-time setting.

    Parameters
    ----------
    result : Any
        Loaded contacts artifact.
    settings : BaseModel
        Active contacts settings.
    allow_missing : bool
        Whether legacy artifacts without RT identity markers may be accepted.
    source : Path or None, optional
        Cache source path used for diagnostics.

    Returns
    -------
    bool
        ``True`` when RT metadata and aggregate content match the request.
    """

    expected = bool(getattr(settings, "compute_residence_times", True))
    metadata = getattr(result, "metadata", None)
    metadata = metadata if isinstance(metadata, dict) else {}
    source_suffix = f" at {source}" if source is not None else ""

    compute_marker = coerce_optional_bool(
        getattr(result, "compute_residence_times", metadata.get("compute_residence_times"))
    )
    computed_marker = coerce_optional_bool(
        getattr(result, "residence_times_computed", metadata.get("residence_times_computed"))
    )
    markers = [marker for marker in (compute_marker, computed_marker) if marker is not None]
    if markers:
        if any(marker != expected for marker in markers):
            logger.warning(
                "contacts: ignoring cache with residence-time setting mismatch%s",
                source_suffix,
            )
            return False
    elif not allow_missing:
        logger.warning(
            "contacts: ignoring cache without residence-time setting proof%s",
            source_suffix,
        )
        return False

    if not expected and cache_has_residence_time_summaries(result):
        logger.warning(
            "contacts: ignoring cache with residence-time summaries while disabled%s",
            source_suffix,
        )
        return False

    return True


def cache_has_residence_time_summaries(result: Any) -> bool:
    """Return whether an aggregate artifact contains RT summary maps.

    Parameters
    ----------
    result : Any
        Loaded contacts artifact.

    Returns
    -------
    bool
        ``True`` when global or per-residue residence-time maps are non-empty.
    """

    if getattr(result, "residence_time_by_polymer_type", None):
        return True
    if getattr(result, "residence_time_by_polymer_type_replicates", None):
        return True
    for residue in getattr(result, "residue_stats", []) or []:
        if getattr(residue, "residence_time_by_polymer_type", None):
            return True
        if getattr(residue, "residence_time_by_polymer_type_per_replicate", None):
            return True
        if getattr(residue, "residence_time_by_polymer_type_replicates", None):
            return True
    return False


def attach_contacts_identity_metadata(result: Any, settings: BaseModel) -> None:
    """Attach current contacts identity metadata to a cache result.

    Parameters
    ----------
    result : Any
        Contacts result object to update in place.
    settings : BaseModel
        Active contacts settings.
    """

    from polyzymd.analyses.shared.config_hash import settings_fingerprint

    metadata = dict(getattr(result, "metadata", {}) or {})
    contact_settings_fp = contacts_settings_fingerprint(settings)
    compute_residence_times = bool(getattr(settings, "compute_residence_times", True))
    effective_selection = effective_polymer_selection(settings)
    metadata["settings_fingerprint"] = contact_settings_fp
    metadata["contacts_settings_fingerprint"] = contact_settings_fp
    metadata["legacy_full_settings_fingerprint"] = settings_fingerprint(settings)
    metadata["settings_fingerprint_domain"] = CONTACTS_SETTINGS_FINGERPRINT_DOMAIN
    metadata["compute_residence_times"] = compute_residence_times
    metadata["residence_times_computed"] = compute_residence_times
    metadata["protein_selection"] = str(settings.protein_selection).strip()
    metadata["polymer_selection"] = str(settings.polymer_selection).strip()
    metadata["effective_polymer_selection"] = effective_selection
    metadata["polymer_types_filter"] = normalize_polymer_types(
        getattr(settings, "polymer_types", None)
    )
    metadata["grouping"] = str(settings.grouping).strip()
    metadata["cutoff"] = float(settings.cutoff)
    metadata["criteria_cutoff"] = float(settings.cutoff)
    result.metadata = metadata


def cache_matches_contacts_settings(
    result: Any,
    settings: BaseModel,
    *,
    source: Path | None = None,
) -> bool:
    """Return whether a cached contacts artifact matches active settings.

    Parameters
    ----------
    result : Any
        Loaded contacts artifact.
    settings : BaseModel
        Active contacts settings.
    source : Path or None, optional
        Cache source path used for diagnostics only.

    Returns
    -------
    bool
        ``True`` when embedded identity matches a compatible contacts
        fingerprint, or cache content provides complete explicit settings proof.
        Filename-only fingerprints are not trusted.
    """

    source_suffix = f" at {source}" if source is not None else ""
    stored_fingerprint = embedded_contacts_settings_fingerprint(result)
    if stored_fingerprint is None:
        if cache_has_contacts_identity_proof(result, settings):
            if not cache_matches_residence_time_setting(
                result,
                settings,
                allow_missing=bool(getattr(settings, "compute_residence_times", True)),
                source=source,
            ):
                return False
            attach_contacts_identity_metadata(result, settings)
            return True
        logger.warning(
            "contacts: ignoring cache without embedded settings fingerprint or "
            "complete settings proof%s",
            source_suffix,
        )
        return False

    candidates = contacts_settings_fingerprint_candidates(settings)
    if str(stored_fingerprint) in candidates:
        allow_missing = bool(getattr(settings, "compute_residence_times", True)) or str(
            stored_fingerprint
        ) == contacts_settings_fingerprint(settings)
        if not cache_matches_residence_time_setting(
            result,
            settings,
            allow_missing=allow_missing,
            source=source,
        ):
            return False
        return True

    logger.warning(
        "contacts: ignoring cache with settings fingerprint mismatch%s: stored=%s, current=%s",
        source_suffix,
        stored_fingerprint,
        candidates[0],
    )
    return False


def cache_matches_context(
    analysis: Any,
    result: Any,
    *,
    settings: BaseModel,
    equilibration: str,
    sim_config: Any,
    replicates: Sequence[int] | None = None,
    allow_replicate_subset: bool = False,
    source: Path | None = None,
) -> bool:
    """Return whether an aggregate result matches the active context.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for compatibility hooks.
    result : Any
        Loaded aggregated contacts result.
    settings : BaseModel
        Active contacts settings.
    equilibration : str
        Active equilibration/window request.
    sim_config : Any
        Active condition simulation configuration.
    replicates : Sequence[int] or None, optional
        Requested replicate tuple for aggregated cache validation.
    allow_replicate_subset : bool, optional
        Whether a stored subset of the requested replicates is accepted.
    source : Path or None, optional
        Cache source path used for settings-fingerprint fallback and diagnostics.

    Returns
    -------
    bool
        ``True`` when settings, window, and known config identity match.
    """

    from polyzymd.analyses.shared.config_hash import validate_config_hash

    if not analysis._cache_matches_window(result, equilibration):
        logger.warning(
            "contacts: ignoring aggregated cache with mismatched equilibration window%s",
            f" at {source}" if source is not None else "",
        )
        return False

    if not analysis._cache_matches_contacts_settings(result, settings, source=source):
        return False

    if replicates is not None and not analysis._cache_matches_replicates(
        result,
        replicates,
        allow_subset=allow_replicate_subset,
        source=source,
    ):
        return False

    stored_hash = getattr(result, "config_hash", None)
    if stored_hash not in (None, "", "unknown") and not validate_config_hash(
        str(stored_hash),
        sim_config,
    ):
        return False

    return True


def cache_matches_replicates(
    analysis: Any,
    result: Any,
    replicates: Sequence[int],
    *,
    allow_subset: bool = False,
    source: Path | None = None,
) -> bool:
    """Return whether cached aggregate replicate identity matches the request.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for compatibility hooks.
    result : Any
        Loaded aggregated contacts result.
    replicates : Sequence[int]
        Requested replicate numbers for the active aggregate.
    allow_subset : bool, optional
        Whether stored replicate identity may be a non-empty subset of the
        requested set.
    source : Path or None, optional
        Cache source path used for diagnostics.

    Returns
    -------
    bool
        ``True`` when the cached result declares the requested replicate set,
        or a subset when explicitly allowed. Missing metadata is treated as
        incompatible because legacy aggregates cannot prove their identity.
    """

    requested = tuple(sorted(int(rep) for rep in replicates))
    stored_replicates = getattr(result, "replicates", None)
    if stored_replicates is None:
        metadata = getattr(result, "metadata", None)
        if isinstance(metadata, dict):
            stored_replicates = metadata.get("replicates")
    if not stored_replicates:
        logger.warning(
            "contacts: ignoring aggregated cache without replicate identity%s",
            f" at {source}" if source is not None else "",
        )
        return False

    try:
        stored = tuple(sorted(int(rep) for rep in stored_replicates))
    except (TypeError, ValueError):
        logger.warning(
            "contacts: ignoring aggregated cache with invalid replicate identity%s",
            f" at {source}" if source is not None else "",
        )
        return False

    if len(set(stored)) != len(stored):
        logger.warning(
            "contacts: ignoring aggregated cache with duplicate replicate identity%s",
            f" at {source}" if source is not None else "",
        )
        return False

    is_allowed_subset = allow_subset and set(stored).issubset(set(requested))
    if stored != requested and not is_allowed_subset:
        logger.warning(
            "contacts: ignoring aggregated cache for replicates %s; requested %s%s",
            stored,
            requested,
            f" at {source}" if source is not None else "",
        )
        return False

    stored_count = len(stored)
    if allow_subset and stored_count < analysis.min_replicates:
        logger.warning(
            "contacts: ignoring aggregated cache with %d replicates below minimum %d%s",
            stored_count,
            analysis.min_replicates,
            f" at {source}" if source is not None else "",
        )
        return False

    declared_count = getattr(result, "n_replicates", None)
    try:
        declared_count = int(declared_count)
    except (TypeError, ValueError):
        logger.warning(
            "contacts: ignoring aggregated cache with invalid n_replicates%s",
            f" at {source}" if source is not None else "",
        )
        return False
    if declared_count != stored_count:
        logger.warning(
            "contacts: ignoring aggregated cache with n_replicates=%d but %d stored replicate IDs%s",
            declared_count,
            stored_count,
            f" at {source}" if source is not None else "",
        )
        return False

    mismatched_vectors = analysis._replicate_vector_length_mismatches(result, stored)
    if mismatched_vectors:
        logger.warning(
            "contacts: ignoring aggregated cache with per-replicate vector length mismatch for %s%s",
            ", ".join(mismatched_vectors),
            f" at {source}" if source is not None else "",
        )
        return False

    return True


def replicate_vector_length_mismatches(
    result: Any,
    expected_replicates: Sequence[int],
) -> list[str]:
    """Return per-replicate vector fields whose lengths are inconsistent.

    Parameters
    ----------
    result : Any
        Loaded aggregated contacts result.
    expected_replicates : Sequence[int]
        Stored replicate identities for the aggregate.

    Returns
    -------
    list[str]
        Field labels for vectors that do not match stored replicate identity.
    """

    mismatches: list[str] = []
    expected = tuple(int(rep) for rep in expected_replicates)
    expected_count = len(expected)
    expected_set = set(expected)
    total_frames = getattr(result, "total_frames_per_replicate", None)
    if isinstance(total_frames, (list, tuple)) and len(total_frames) != expected_count:
        mismatches.append("total_frames_per_replicate")

    residue_stats = getattr(result, "residue_stats", []) or []
    for index, residue in enumerate(residue_stats):
        contact_fractions = getattr(residue, "contact_fraction_per_replicate", None)
        if (
            isinstance(contact_fractions, (list, tuple))
            and len(contact_fractions) != expected_count
        ):
            mismatches.append(f"residue_stats[{index}].contact_fraction_per_replicate")

        per_type_fractions = getattr(residue, "by_polymer_type_per_replicate", {}) or {}
        if isinstance(per_type_fractions, dict):
            for polymer_type, values in per_type_fractions.items():
                if isinstance(values, (list, tuple)) and len(values) != expected_count:
                    mismatches.append(
                        f"residue_stats[{index}].by_polymer_type_per_replicate[{polymer_type}]"
                    )

        residence_times = getattr(residue, "residence_time_by_polymer_type_per_replicate", {}) or {}
        residence_replicates = (
            getattr(residue, "residence_time_by_polymer_type_replicates", {}) or {}
        )
        if isinstance(residence_times, dict):
            for polymer_type, values in residence_times.items():
                field = (
                    "residue_stats"
                    f"[{index}].residence_time_by_polymer_type_per_replicate"
                    f"[{polymer_type}]"
                )
                if not isinstance(values, (list, tuple)):
                    continue
                sparse_ids = None
                if isinstance(residence_replicates, dict):
                    sparse_ids = residence_replicates.get(polymer_type)
                if sparse_ids is None:
                    # Legacy sparse caches lack identity metadata; accept compact vectors
                    if len(values) > expected_count:
                        mismatches.append(field)
                    continue
                if not isinstance(sparse_ids, (list, tuple)):
                    mismatches.append(f"{field}.replicates")
                    continue
                try:
                    sparse_replicates = [int(rep) for rep in sparse_ids]
                except (TypeError, ValueError):
                    mismatches.append(f"{field}.replicates")
                    continue
                if len(values) != len(sparse_replicates):
                    mismatches.append(field)
                    continue
                if len(set(sparse_replicates)) != len(sparse_replicates):
                    mismatches.append(f"{field}.replicates")
                    continue
                if not set(sparse_replicates).issubset(expected_set):
                    mismatches.append(f"{field}.replicates")

    return mismatches


def load_validated_aggregated_result(
    analysis: Any,
    aggregated_dir: Path,
    *,
    settings: BaseModel,
    equilibration: str,
    replicates: Sequence[int],
    sim_config: Any,
    recompute: bool,
    allow_replicate_subset: bool = False,
) -> Any | None:
    """Load an aggregated contacts result after cache identity validation.

    Fingerprinted sidecars are preferred over canonical ``result.json`` so a
    stale canonical aggregate cannot shadow a matching window-aware contacts
    aggregate.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for cache and path compatibility hooks.
    aggregated_dir : Path
        Directory containing aggregate cache files.
    settings : BaseModel
        Active contacts settings.
    equilibration : str
        Active equilibration/window request.
    replicates : Sequence[int]
        Replicates expected in the aggregate.
    sim_config : Any
        Active condition simulation configuration.
    recompute : bool
        Whether cache loading should be skipped.
    allow_replicate_subset : bool, optional
        Whether finalized aggregates with successful replicate subsets may be
        loaded.

    Returns
    -------
    Any or None
        Valid aggregated contacts result, or ``None`` when no compatible cache
        is available.
    """

    candidates = analysis._aggregated_cache_candidates(
        aggregated_dir,
        settings,
        equilibration,
        replicates,
        result_path=analysis.aggregate_result_path(aggregated_dir),
    )
    saw_json = aggregated_dir.exists() and any(aggregated_dir.glob("*.json"))

    for candidate in candidates:
        cached = analysis._load_cache_candidate(
            AggregatedContactResult,
            candidate,
            recompute=recompute,
            sim_config=None,
            settings=None,
        )
        if cached is None:
            continue
        if analysis._cache_matches_context(
            cached,
            settings=settings,
            equilibration=equilibration,
            sim_config=sim_config,
            replicates=replicates,
            allow_replicate_subset=allow_replicate_subset,
            source=candidate,
        ):
            return cached

    if saw_json:
        return None

    return analysis._load_aggregated_result(aggregated_dir)
