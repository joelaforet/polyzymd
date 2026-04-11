"""Contacts-specific path discovery helpers.

This module centralizes contacts file naming conventions so shared analysis
helpers remain storage-agnostic.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence


def find_contact_result_for_replicate(analysis_dir: Path, replicate: int) -> Path | None:
    """Locate a contacts JSON result for one replicate.

    Search order is deterministic and returns the first matching path:

    - ``run_{replicate}/result.json``
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

    Returns
    -------
    Path | None
        Contact result path if found, otherwise ``None``.
    """
    canonical = analysis_dir / f"run_{replicate}" / "result.json"
    if canonical.exists():
        return canonical

    matches = sorted(analysis_dir.glob(f"contacts_eq*_cut*_rep{replicate}.json"))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(
            f"Ambiguous contacts cache for replicate {replicate}: "
            f"found {len(matches)} files in {analysis_dir}: "
            + ", ".join(str(m.name) for m in matches)
        )

    matches = sorted(
        (analysis_dir / f"run_{replicate}").glob(f"contacts_eq*_cut*_rep{replicate}.json")
    )
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(
            f"Ambiguous contacts cache for replicate {replicate}: "
            f"found {len(matches)} files in {analysis_dir / f'run_{replicate}'}: "
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
    analysis_dir: Path, replicates: Sequence[int]
) -> dict[int, Path]:
    """Resolve contacts result paths for all requested replicates.

    Parameters
    ----------
    analysis_dir : Path
        Contacts analysis directory for one condition.
    replicates : Sequence[int]
        Replicate IDs to resolve.

    Returns
    -------
    dict[int, Path]
        Mapping of replicate ID to discovered contacts result path.
    """
    found: dict[int, Path] = {}
    for replicate in replicates:
        path = find_contact_result_for_replicate(analysis_dir, replicate)
        if path is not None:
            found[replicate] = path
    return found
