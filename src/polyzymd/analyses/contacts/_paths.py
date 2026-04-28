"""Contacts-specific path discovery helpers.

This module centralizes contacts file naming conventions so shared analysis
helpers remain storage-agnostic.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Sequence


def find_contact_result_for_replicate(
    analysis_dir: Path,
    replicate: int,
    *,
    settings_fp: str | None = None,
    equilibration: str | None = None,
    strict_identity: bool = False,
) -> Path | None:
    """Locate a contacts JSON result for one replicate.

    Search order is deterministic and returns the first matching path.
    Fingerprinted sidecars are preferred before canonical
    ``run_{replicate}/result.json`` only when a settings fingerprint is
    requested because the canonical path may be stale after a settings change:

    - ``contacts_eq*_cut*_s{settings_fp}_rep{replicate}.json`` (when ``settings_fp`` is set)
    - ``run_{replicate}/contacts_eq*_cut*_s{settings_fp}_rep{replicate}.json``
      (when ``settings_fp`` is set)
    - ``run_{replicate}/result.json`` (when it declares ``settings_fp`` for
      explicit settings lookups, otherwise after fingerprinted sidecars)
    - Any ``_s*`` fingerprinted sidecar blocks explicit legacy fallback when
      it does not embed matching contacts identity metadata
    - ``contacts_eq*_cut*_rep{replicate}.json`` with matching metadata
    - ``run_{replicate}/contacts_eq*_cut*_rep{replicate}.json`` with matching metadata
    - ``contacts_rep{replicate}.json`` with matching metadata
    - ``run_{replicate}/contacts_rep{replicate}.json`` with matching metadata
    - Fingerprinted and legacy paths without explicit settings lookups follow
      the compatibility order documented by the code below

    Parameters
    ----------
    analysis_dir : Path
        Contacts analysis directory for one condition.
    replicate : int
        Replicate number.
    settings_fp : str or None, optional
        Settings fingerprint. When provided, fingerprinted cache filenames are
        searched before canonical and legacy patterns. Non-matching
        fingerprinted sidecars block fallback unless canonical metadata proves
        the requested settings identity.
    equilibration : str or None, optional
        Requested contacts analysis window. When provided, only artifacts whose
        filename or metadata proves the same equilibration window are returned.
    strict_identity : bool, optional
        Require fingerprinted contacts sidecars to carry embedded contacts
        identity metadata when no explicit settings fingerprint is known.

    Returns
    -------
    Path | None
        Contact result path if found, otherwise ``None``.
    """
    run_dir = analysis_dir / f"run_{replicate}"
    eq_token = contact_equilibration_filename_token(equilibration)
    eq_glob = eq_token if eq_token is not None else "eq*"
    fingerprinted_patterns = (
        (analysis_dir, f"contacts_{eq_glob}_cut*_s*_rep{replicate}.json"),
        (run_dir, f"contacts_{eq_glob}_cut*_s*_rep{replicate}.json"),
    )

    canonical = analysis_dir / f"run_{replicate}" / "result.json"

    if settings_fp is not None:
        fingerprinted = _first_unique_match(
            (
                (
                    analysis_dir,
                    f"contacts_{eq_glob}_cut*_s{settings_fp}_rep{replicate}.json",
                ),
                (run_dir, f"contacts_{eq_glob}_cut*_s{settings_fp}_rep{replicate}.json"),
            ),
            replicate,
            equilibration=equilibration,
        )
        if fingerprinted is not None and contact_artifact_matches_settings_fingerprint(
            fingerprinted,
            settings_fp,
        ):
            return fingerprinted

        if (
            canonical.exists()
            and contact_artifact_matches_window(canonical, equilibration)
            and contact_artifact_matches_settings_fingerprint(canonical, settings_fp)
        ):
            return canonical

        if _has_any_fingerprinted_match(fingerprinted_patterns):
            return None

        return _first_unique_metadata_proven_legacy_match(
            analysis_dir,
            run_dir,
            replicate,
            settings_fp,
            eq_glob=eq_glob,
            equilibration=equilibration,
        )

    if settings_fp is None:
        if strict_identity:
            discovered_fingerprinted = _first_unique_proven_contact_identity_match(
                fingerprinted_patterns,
                replicate,
                equilibration=equilibration,
            )
            if discovered_fingerprinted is not None:
                return discovered_fingerprinted
            if _has_any_fingerprinted_match(fingerprinted_patterns):
                return None
        elif canonical.exists() and contact_artifact_matches_window(canonical, equilibration):
            return canonical
        else:
            discovered_fingerprinted = _first_unique_match(
                fingerprinted_patterns,
                replicate,
                equilibration=equilibration,
            )
            if discovered_fingerprinted is not None:
                return discovered_fingerprinted

    if canonical.exists() and contact_artifact_matches_window(canonical, equilibration):
        return canonical

    search_patterns: list[tuple[Path, str]] = [
        (analysis_dir, f"contacts_{eq_glob}_cut*_rep{replicate}.json"),
        (run_dir, f"contacts_{eq_glob}_cut*_rep{replicate}.json"),
    ]

    matched = _first_unique_match(search_patterns, replicate, equilibration=equilibration)
    if matched is not None:
        return matched

    legacy = analysis_dir / f"contacts_rep{replicate}.json"
    if legacy.exists() and contact_artifact_matches_window(legacy, equilibration):
        return legacy

    legacy_run = analysis_dir / f"run_{replicate}" / f"contacts_rep{replicate}.json"
    if legacy_run.exists() and contact_artifact_matches_window(legacy_run, equilibration):
        return legacy_run

    return None


def _first_unique_metadata_proven_legacy_match(
    analysis_dir: Path,
    run_dir: Path,
    replicate: int,
    settings_fp: str,
    *,
    eq_glob: str,
    equilibration: str | None = None,
) -> Path | None:
    """Return a non-fingerprinted legacy contacts artifact with matching metadata.

    Parameters
    ----------
    analysis_dir : Path
        Contacts analysis directory for one condition.
    run_dir : Path
        Per-replicate contacts output directory.
    replicate : int
        Replicate number used in ambiguity diagnostics.
    settings_fp : str
        Required contacts settings fingerprint from embedded JSON metadata.
    eq_glob : str
        Equilibration filename token or glob used for windowed legacy names.
    equilibration : str or None, optional
        Requested window used to discard mismatched candidates.

    Returns
    -------
    Path or None
        The unique metadata-proven legacy contacts path, or ``None`` when no
        legacy candidate proves the requested contacts identity.
    """

    search_patterns: tuple[tuple[Path, str], ...] = (
        (analysis_dir, f"contacts_{eq_glob}_cut*_rep{replicate}.json"),
        (run_dir, f"contacts_{eq_glob}_cut*_rep{replicate}.json"),
        (analysis_dir, f"contacts_rep{replicate}.json"),
        (run_dir, f"contacts_rep{replicate}.json"),
    )
    for search_dir, pattern in search_patterns:
        matches = sorted(
            path
            for path in search_dir.glob(pattern)
            if contact_artifact_matches_window(path, equilibration)
            and _contact_artifact_matches_explicit_contacts_fingerprint(path, settings_fp)
        )
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValueError(
                f"Ambiguous contacts cache for replicate {replicate}: "
                f"found {len(matches)} files matching '{pattern}' in {search_dir}: "
                + ", ".join(str(m.name) for m in matches)
            )

    return None


def _contact_artifact_matches_explicit_contacts_fingerprint(path: Path, settings_fp: str) -> bool:
    """Return whether a contacts artifact embeds the requested contacts fingerprint.

    Parameters
    ----------
    path : Path
        Contacts artifact path to inspect.
    settings_fp : str
        Required contacts settings fingerprint.

    Returns
    -------
    bool
        ``True`` only when the artifact records ``contacts_settings_fingerprint``
        or the legacy singular spelling. Generic settings fingerprints are not
        accepted for legacy sidecar fallback.
    """

    try:
        data = json.loads(path.read_text())
    except (FileNotFoundError, OSError, json.JSONDecodeError):
        return False
    if not isinstance(data, dict):
        return False

    direct = data.get("contacts_settings_fingerprint") or data.get("contact_settings_fingerprint")
    if isinstance(direct, str) and direct == settings_fp:
        return True

    metadata = data.get("metadata")
    if not isinstance(metadata, dict):
        return False
    stored = metadata.get("contacts_settings_fingerprint") or metadata.get(
        "contact_settings_fingerprint"
    )
    return isinstance(stored, str) and stored == settings_fp


def _has_any_fingerprinted_match(search_patterns: Sequence[tuple[Path, str]]) -> bool:
    """Return whether any fingerprinted contacts sidecar exists.

    Parameters
    ----------
    search_patterns : Sequence[tuple[Path, str]]
        Ordered ``(directory, glob_pattern)`` pairs to inspect.

    Returns
    -------
    bool
        ``True`` when at least one path matches any pattern.
    """

    return any(any(search_dir.glob(pattern)) for search_dir, pattern in search_patterns)


def _first_unique_proven_contact_identity_match(
    search_patterns: Sequence[tuple[Path, str]],
    replicate: int,
    *,
    equilibration: str | None = None,
) -> Path | None:
    """Return one fingerprinted contacts match with embedded identity proof.

    Parameters
    ----------
    search_patterns : Sequence[tuple[Path, str]]
        Ordered ``(directory, glob_pattern)`` pairs to inspect.
    replicate : int
        Replicate number used in ambiguity diagnostics.
    equilibration : str or None, optional
        Requested window used to discard mismatched candidates.

    Returns
    -------
    Path or None
        The unique proven matched path, or ``None`` when no candidate proves
        contacts identity.
    """

    for search_dir, pattern in search_patterns:
        matches = sorted(
            path
            for path in search_dir.glob(pattern)
            if contact_artifact_matches_window(path, equilibration)
            and contact_artifact_has_contacts_identity_proof(path)
        )
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValueError(
                f"Ambiguous contacts cache for replicate {replicate}: "
                f"found {len(matches)} files matching '{pattern}' in {search_dir}: "
                + ", ".join(str(m.name) for m in matches)
            )

    return None


def _first_unique_match(
    search_patterns: Sequence[tuple[Path, str]],
    replicate: int,
    *,
    equilibration: str | None = None,
) -> Path | None:
    """Return the first unique glob match from ordered contacts patterns.

    Parameters
    ----------
    search_patterns : Sequence[tuple[Path, str]]
        Ordered ``(directory, glob_pattern)`` pairs to inspect.
    replicate : int
        Replicate number used in ambiguity diagnostics.
    equilibration : str or None, optional
        Requested window used to discard mismatched candidates.

    Returns
    -------
    Path or None
        The unique matched path, or ``None`` when no pattern matches.
    """

    for search_dir, pattern in search_patterns:
        matches = sorted(
            path
            for path in search_dir.glob(pattern)
            if contact_artifact_matches_window(path, equilibration)
        )
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValueError(
                f"Ambiguous contacts cache for replicate {replicate}: "
                f"found {len(matches)} files matching '{pattern}' in {search_dir}: "
                + ", ".join(str(m.name) for m in matches)
            )

    return None


def find_contact_results_for_replicates(
    analysis_dir: Path,
    replicates: Sequence[int],
    *,
    settings_fp: str | None = None,
    equilibration: str | None = None,
    strict_identity: bool = False,
) -> dict[int, Path]:
    """Resolve contacts result paths for all requested replicates.

    Parameters
    ----------
    analysis_dir : Path
        Contacts analysis directory for one condition.
    replicates : Sequence[int]
        Replicate IDs to resolve.
    settings_fp : str or None, optional
        Settings fingerprint. When provided, fingerprinted cache filenames are
        searched before legacy patterns.
    equilibration : str or None, optional
        Requested contacts analysis window.
    strict_identity : bool, optional
        Require fingerprinted sidecars to carry embedded contacts identity when
        no explicit settings fingerprint is known.

    Returns
    -------
    dict[int, Path]
        Mapping of replicate ID to discovered contacts result path.
    """
    found: dict[int, Path] = {}
    for replicate in replicates:
        path = find_contact_result_for_replicate(
            analysis_dir,
            replicate,
            settings_fp=settings_fp,
            equilibration=equilibration,
            strict_identity=strict_identity,
        )
        if path is not None:
            found[replicate] = path
    return found


def extract_contact_settings_fingerprint(path: Path) -> str | None:
    """Extract the contacts settings fingerprint recorded by a contacts artifact.

    Parameters
    ----------
    path : Path
        Contacts result path.

    Returns
    -------
    str or None
        Settings fingerprint from result metadata, if available.
    """

    try:
        data = json.loads(path.read_text())
    except (FileNotFoundError, OSError, json.JSONDecodeError):
        return None

    if isinstance(data, dict):
        if path.name.startswith("binding_preference"):
            return _extract_contacts_fingerprint_from_binding_preference(data)
        direct = (
            data.get("contacts_settings_fingerprint")
            or data.get("contact_settings_fingerprint")
            or data.get("settings_fingerprint")
            or data.get("settings_fp")
        )
        if isinstance(direct, str):
            return direct
        metadata = data.get("metadata")
        if isinstance(metadata, dict):
            stored = (
                metadata.get("contacts_settings_fingerprint")
                or metadata.get("contact_settings_fingerprint")
                or metadata.get("settings_fingerprint")
                or metadata.get("settings_fp")
            )
            if isinstance(stored, str):
                return stored

    return None


def contact_artifact_has_contacts_identity_proof(path: Path) -> bool:
    """Return whether a contacts artifact embeds contacts identity metadata.

    Filename ``_s`` tokens are intentionally ignored because they are only
    candidate locators and do not prove the artifact contents were computed with
    that contacts settings domain.

    Parameters
    ----------
    path : Path
        Contacts artifact path to inspect.

    Returns
    -------
    bool
        ``True`` when JSON metadata records a contacts settings fingerprint.
    """

    return extract_contact_settings_fingerprint(path) is not None


def has_unproven_fingerprinted_contact_artifacts(
    analysis_dir: Path,
    replicates: Sequence[int],
    *,
    equilibration: str | None = None,
) -> bool:
    """Return whether fingerprinted contacts sidecars block strict fallback.

    A replicate blocks strict downstream use when matching ``_s*`` contacts
    sidecars exist for the requested window but none embed contacts identity
    metadata. Legacy artifacts remain usable when no fingerprinted sidecar is
    present for that replicate and window.

    Parameters
    ----------
    analysis_dir : Path
        Contacts analysis directory for one condition.
    replicates : Sequence[int]
        Replicate IDs to inspect.
    equilibration : str or None, optional
        Requested contacts analysis window.

    Returns
    -------
    bool
        ``True`` when at least one requested replicate has only unproven
        fingerprinted sidecars.
    """

    for replicate in replicates:
        artifacts = _iter_fingerprinted_contact_artifacts(
            analysis_dir,
            (replicate,),
            equilibration=equilibration,
        )
        if artifacts and not any(
            contact_artifact_has_contacts_identity_proof(path) for path in artifacts
        ):
            return True
    return False


def _extract_contacts_fingerprint_from_binding_preference(
    data: dict[str, Any],
) -> str | None:
    """Extract upstream contacts identity from a BP artifact.

    Parameters
    ----------
    data : dict[str, Any]
        Parsed JSON payload.

    Returns
    -------
    str or None
        Contacts settings fingerprint from explicitly contact-domain BP metadata.
    """

    direct = data.get("contacts_settings_fingerprint") or data.get("contact_settings_fingerprint")
    if isinstance(direct, str):
        return direct

    metadata = data.get("metadata")
    if isinstance(metadata, dict):
        stored = metadata.get("contacts_settings_fingerprint") or metadata.get(
            "contact_settings_fingerprint"
        )
        if isinstance(stored, str):
            return stored
    return None


def infer_contact_results_settings_fingerprint(paths_by_replicate: dict[int, Path]) -> str | None:
    """Infer the contacts settings fingerprint shared by resolved artifacts.

    Parameters
    ----------
    paths_by_replicate : dict[int, Path]
        Resolved contacts result paths keyed by replicate.

    Returns
    -------
    str or None
        Unique settings fingerprint found across contacts artifacts, or ``None``
        when no artifact records one.

    Raises
    ------
    ValueError
        If resolved artifacts record conflicting settings fingerprints.
    """

    fingerprints = {
        fingerprint
        for path in paths_by_replicate.values()
        if (fingerprint := extract_contact_settings_fingerprint(path)) is not None
    }
    if len(fingerprints) > 1:
        raise ValueError(
            "Resolved contacts artifacts have inconsistent settings fingerprints: "
            + ", ".join(sorted(fingerprints))
        )
    return next(iter(fingerprints), None)


def infer_contacts_artifact_settings_fingerprint(
    analysis_dir: Path,
    replicates: Sequence[int],
    *,
    equilibration: str | None = None,
) -> str | None:
    """Infer the contacts settings fingerprint from available artifacts.

    Contacts-derived downstream analyses may not know the exact contacts
    settings used to produce upstream artifacts. This helper treats
    fingerprinted contacts sidecars as candidate locators only, then requires
    embedded contacts metadata or content proof before inferring identity.

    Parameters
    ----------
    analysis_dir : Path
        Contacts analysis directory for one condition.
    replicates : Sequence[int]
        Replicate IDs associated with the condition.
    equilibration : str or None, optional
        Requested contacts analysis window. Artifacts from other windows are
        ignored while inferring the upstream settings identity.

    Returns
    -------
    str or None
        Unique discovered contacts settings fingerprint, or ``None`` when no
        artifact records one.

    Raises
    ------
    ValueError
        If artifacts record conflicting settings fingerprints.
    """

    fingerprinted_contact_artifacts = _iter_fingerprinted_contact_artifacts(
        analysis_dir,
        replicates,
        equilibration=equilibration,
    )
    fingerprinted_contact_fingerprints = {
        fingerprint
        for path in fingerprinted_contact_artifacts
        if (fingerprint := extract_contact_settings_fingerprint(path)) is not None
    }
    if fingerprinted_contact_fingerprints:
        return _unique_contacts_settings_fingerprint(fingerprinted_contact_fingerprints)
    if fingerprinted_contact_artifacts:
        return None

    resolved_contact_results = find_contact_results_for_replicates(
        analysis_dir,
        replicates,
        equilibration=equilibration,
        strict_identity=True,
    )
    canonical_contact_fingerprints = {
        fingerprint
        for path in resolved_contact_results.values()
        if _contact_artifact_defines_current_identity(path)
        if (fingerprint := extract_contact_settings_fingerprint(path)) is not None
    }
    if canonical_contact_fingerprints:
        return _unique_contacts_settings_fingerprint(canonical_contact_fingerprints)

    legacy_contact_fingerprints = {
        fingerprint
        for path in _iter_non_fingerprinted_contact_artifacts(
            analysis_dir,
            replicates,
            equilibration=equilibration,
        )
        if (fingerprint := extract_contact_settings_fingerprint(path)) is not None
    }
    if legacy_contact_fingerprints:
        return _unique_contacts_settings_fingerprint(legacy_contact_fingerprints)

    resolved_contact_fingerprints = {
        fingerprint
        for path in resolved_contact_results.values()
        if (fingerprint := extract_contact_settings_fingerprint(path)) is not None
    }
    if resolved_contact_fingerprints:
        return _unique_contacts_settings_fingerprint(resolved_contact_fingerprints)

    fingerprints = set()
    for path in _iter_fingerprinted_binding_preference_artifacts(
        analysis_dir,
        replicates,
        equilibration=equilibration,
    ):
        fingerprint = extract_contact_settings_fingerprint(path)
        if fingerprint is not None:
            fingerprints.add(fingerprint)

    return _unique_contacts_settings_fingerprint(fingerprints)


def _iter_non_fingerprinted_contact_artifacts(
    analysis_dir: Path,
    replicates: Sequence[int],
    *,
    equilibration: str | None = None,
) -> list[Path]:
    """Return legacy contacts artifacts in deterministic order.

    Filename ``_s`` tokens are excluded so unproven fingerprinted sidecars keep
    their barrier semantics in callers that already checked them first.

    Parameters
    ----------
    analysis_dir : Path
        Contacts analysis directory for one condition.
    replicates : Sequence[int]
        Replicate IDs associated with the condition.
    equilibration : str or None, optional
        Requested contacts analysis window.

    Returns
    -------
    list[Path]
        Unique non-fingerprinted contacts artifact paths.
    """

    eq_token = contact_equilibration_filename_token(equilibration)
    eq_glob = eq_token if eq_token is not None else "eq*"
    paths: list[Path] = []
    for replicate in replicates:
        run_dir = analysis_dir / f"run_{replicate}"
        paths.extend(sorted(analysis_dir.glob(f"contacts_{eq_glob}_cut*_rep{replicate}.json")))
        paths.extend(sorted(run_dir.glob(f"contacts_{eq_glob}_cut*_rep{replicate}.json")))
        paths.append(analysis_dir / f"contacts_rep{replicate}.json")
        paths.append(run_dir / f"contacts_rep{replicate}.json")

    return [
        path
        for path in dict.fromkeys(paths)
        if path.exists()
        and "_s" not in path.name
        and contact_artifact_matches_window(path, equilibration)
    ]


def _unique_contacts_settings_fingerprint(fingerprints: set[str]) -> str | None:
    """Return the unique contacts settings fingerprint from a candidate set.

    Parameters
    ----------
    fingerprints : set[str]
        Candidate settings fingerprints discovered from compatible artifacts.

    Returns
    -------
    str or None
        Unique settings fingerprint, or ``None`` when no candidates were found.

    Raises
    ------
    ValueError
        If artifacts record conflicting settings fingerprints.
    """

    if len(fingerprints) > 1:
        raise ValueError(
            "Resolved contacts artifacts have inconsistent settings fingerprints: "
            + ", ".join(sorted(fingerprints))
        )
    return next(iter(fingerprints), None)


def _iter_fingerprinted_contact_artifacts(
    analysis_dir: Path,
    replicates: Sequence[int],
    *,
    equilibration: str | None = None,
) -> list[Path]:
    """Return fingerprinted contacts artifacts in deterministic order.

    Parameters
    ----------
    analysis_dir : Path
        Contacts analysis directory for one condition.
    replicates : Sequence[int]
        Replicate IDs associated with the condition.
    equilibration : str or None, optional
        Requested contacts analysis window.

    Returns
    -------
    list[Path]
        Unique fingerprinted contacts artifact paths.
    """

    eq_token = contact_equilibration_filename_token(equilibration)
    eq_glob = eq_token if eq_token is not None else "eq*"
    paths: list[Path] = []
    for replicate in replicates:
        paths.extend(sorted(analysis_dir.glob(f"contacts_{eq_glob}_cut*_s*_rep{replicate}.json")))
        paths.extend(
            sorted(
                (analysis_dir / f"run_{replicate}").glob(
                    f"contacts_{eq_glob}_cut*_s*_rep{replicate}.json"
                )
            )
        )
    return [
        path
        for path in dict.fromkeys(paths)
        if contact_artifact_matches_window(path, equilibration)
    ]


def _iter_fingerprinted_binding_preference_artifacts(
    analysis_dir: Path,
    replicates: Sequence[int],
    *,
    equilibration: str | None = None,
) -> list[Path]:
    """Return fingerprinted binding-preference artifacts in deterministic order.

    Parameters
    ----------
    analysis_dir : Path
        Contacts analysis directory for one condition.
    replicates : Sequence[int]
        Replicate IDs associated with the condition.
    equilibration : str or None, optional
        Requested contacts analysis window.

    Returns
    -------
    list[Path]
        Unique fingerprinted binding-preference artifact paths.
    """

    patterns = [
        "binding_preference_aggregated_s*.json",
        "binding_preference_aggregated_s*_reps*.json",
        "binding_preference_s*.json",
    ]
    paths: list[Path] = []
    for pattern in patterns:
        paths.extend(sorted(analysis_dir.glob(pattern)))
    for replicate in replicates:
        paths.extend(sorted(analysis_dir.glob(f"binding_preference_s*_rep{replicate}.json")))
    return [
        path
        for path in dict.fromkeys(paths)
        if binding_preference_artifact_matches_window(path, equilibration)
    ]


def contact_equilibration_filename_token(equilibration: str | None) -> str | None:
    """Return the contacts filename token for an equilibration window.

    Parameters
    ----------
    equilibration : str or None
        Equilibration string such as ``"10ns"``.

    Returns
    -------
    str or None
        Token such as ``"eq10ns"``, or ``None`` when no window was requested.
    """

    if equilibration is None:
        return None

    from polyzymd.analyses.shared.loader import parse_time_string

    eq_value, eq_unit = parse_time_string(equilibration)
    return f"eq{eq_value:g}{eq_unit}"


def contact_artifact_matches_window(path: Path, equilibration: str | None) -> bool:
    """Return whether a contacts artifact proves the requested window.

    Parameters
    ----------
    path : Path
        Contacts result path.
    equilibration : str or None
        Requested equilibration window.

    Returns
    -------
    bool
        ``True`` when no window was requested, or when filename or JSON metadata
        declares the requested equilibration window.
    """

    if equilibration is None:
        return True

    if _filename_has_equilibration_token(path, equilibration):
        return True
    return _json_artifact_matches_window(path, equilibration)


def contact_artifact_matches_settings_fingerprint(path: Path, settings_fp: str) -> bool:
    """Return whether a contacts artifact declares the requested settings identity.

    Parameters
    ----------
    path : Path
        Contacts result path.
    settings_fp : str
        Requested contacts settings fingerprint.

    Returns
    -------
    bool
        ``True`` when JSON metadata declares ``settings_fp``.
    """

    return extract_contact_settings_fingerprint(path) == settings_fp


def binding_preference_artifact_matches_window(path: Path, equilibration: str | None) -> bool:
    """Return whether a binding-preference artifact proves the requested window.

    Parameters
    ----------
    path : Path
        Binding-preference result path.
    equilibration : str or None
        Requested equilibration window.

    Returns
    -------
    bool
        ``True`` when no window was requested, or when artifact metadata records
        the requested contacts window.
    """

    if equilibration is None:
        return True
    return _json_artifact_matches_window(path, equilibration)


def _filename_has_equilibration_token(path: Path, equilibration: str) -> bool:
    """Return whether a cache filename includes the requested equilibration token."""

    token = contact_equilibration_filename_token(equilibration)
    return token is not None and f"_{token}_" in f"_{path.name}_"


def _contact_artifact_defines_current_identity(path: Path) -> bool:
    """Return whether a canonical artifact is a complete contacts result."""

    if path.name != "result.json" or not path.parent.name.startswith("run_"):
        return False

    try:
        data = json.loads(path.read_text())
    except (FileNotFoundError, OSError, json.JSONDecodeError):
        return False
    if not isinstance(data, dict):
        return False

    required_fields = {"criteria_cutoff", "criteria_label", "n_frames", "residue_contacts"}
    return required_fields.issubset(data)


def _json_artifact_matches_window(path: Path, equilibration: str) -> bool:
    """Read JSON metadata and compare it to the requested window."""

    try:
        data = json.loads(path.read_text())
    except (FileNotFoundError, OSError, json.JSONDecodeError):
        return False
    if not isinstance(data, dict):
        return False

    candidates = [data]
    metadata = data.get("metadata")
    if isinstance(metadata, dict):
        candidates.append(metadata)

    for candidate in candidates:
        if _window_mapping_matches(candidate, equilibration):
            return True
    return False


def _window_mapping_matches(mapping: dict[str, Any], equilibration: str) -> bool:
    """Return whether a JSON mapping declares the requested window."""

    stored_label = mapping.get("equilibration") or mapping.get("equilibration_time_label")
    if isinstance(stored_label, str) and _window_strings_match(stored_label, equilibration):
        return True

    stored_time = mapping.get("equilibration_time")
    stored_unit = mapping.get("equilibration_unit")
    if stored_time is None or stored_unit is None:
        return False
    return _window_strings_match(f"{stored_time}{stored_unit}", equilibration)


def _window_strings_match(stored: str, requested: str) -> bool:
    """Compare two time strings after conversion to picoseconds."""

    from polyzymd.analyses.shared.loader import convert_time, parse_time_string

    try:
        stored_time, stored_unit = parse_time_string(stored)
        requested_time, requested_unit = parse_time_string(requested)
        stored_ps = convert_time(float(stored_time), stored_unit, "ps")
        requested_ps = convert_time(float(requested_time), requested_unit, "ps")
    except (TypeError, ValueError):
        return False
    return math.isclose(stored_ps, requested_ps, rel_tol=0.0, abs_tol=1.0e-9)
