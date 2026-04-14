"""Contacts-specific path discovery helpers.

This module centralizes contacts file naming conventions so shared analysis
helpers remain storage-agnostic.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence


def find_contact_result_for_replicate(
    analysis_dir: Path,
    replicate: int,
    *,
    settings_fp: str | None = None,
) -> Path | None:
    """Locate a contacts JSON result for one replicate.

    Search order is deterministic and returns the first matching path:

    - ``run_{replicate}/result.json``
    - ``contacts_eq*_cut*_s{settings_fp}_rep{replicate}.json`` (when ``settings_fp`` is set)
    - ``run_{replicate}/contacts_eq*_cut*_s{settings_fp}_rep{replicate}.json``
      (when ``settings_fp`` is set)
    - ``contacts_eq*_cut*_rep{replicate}.json``
    - ``run_{replicate}/contacts_eq*_cut*_rep{replicate}.json``
    - ``contacts_rep{replicate}.json``
    - ``run_{replicate}/contacts_rep{replicate}.json``

    Parameters
    ----------
    analysis_dir : Path
        Contacts analysis directory for one condition.
    replicate : int
        Replicate number.
    settings_fp : str or None, optional
        Settings fingerprint. When provided, fingerprinted cache filenames are
        searched before legacy patterns.

    Returns
    -------
    Path | None
        Contact result path if found, otherwise ``None``.
    """
    canonical = analysis_dir / f"run_{replicate}" / "result.json"
    if canonical.exists():
        return canonical

    run_dir = analysis_dir / f"run_{replicate}"
    search_patterns: list[tuple[Path, str]] = []
    if settings_fp is not None:
        search_patterns.extend(
            [
                (analysis_dir, f"contacts_eq*_cut*_s{settings_fp}_rep{replicate}.json"),
                (run_dir, f"contacts_eq*_cut*_s{settings_fp}_rep{replicate}.json"),
            ]
        )
    search_patterns.extend(
        [
            (analysis_dir, f"contacts_eq*_cut*_rep{replicate}.json"),
            (run_dir, f"contacts_eq*_cut*_rep{replicate}.json"),
        ]
    )

    for search_dir, pattern in search_patterns:
        matches = sorted(search_dir.glob(pattern))
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValueError(
                f"Ambiguous contacts cache for replicate {replicate}: "
                f"found {len(matches)} files matching '{pattern}' in {search_dir}: "
                + ", ".join(str(m.name) for m in matches)
            )

    legacy = analysis_dir / f"contacts_rep{replicate}.json"
    if legacy.exists():
        return legacy

    legacy_run = analysis_dir / f"run_{replicate}" / f"contacts_rep{replicate}.json"
    if legacy_run.exists():
        return legacy_run

    return None


def find_contact_results_for_replicates(
    analysis_dir: Path,
    replicates: Sequence[int],
    *,
    settings_fp: str | None = None,
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
        )
        if path is not None:
            found[replicate] = path
    return found
