"""Diagnostics for archived PolyzyMD features."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ArchivedAnalysisPlugin:
    """Archive metadata for a removed analysis plugin.

    Parameters
    ----------
    name : str
        Canonical analysis plugin name.
    tag : str
        Git tag containing the archived implementation.
    branch : str
        Branch where the implementation last shipped as active code.
    """

    name: str
    tag: str
    branch: str


ARCHIVED_ANALYSIS_PLUGINS: dict[str, ArchivedAnalysisPlugin] = {
    "binding_free_energy": ArchivedAnalysisPlugin(
        name="binding_free_energy",
        tag="archive_experimental_analysis",
        branch="feature/mda-analysis-migration",
    ),
    "bfe": ArchivedAnalysisPlugin(
        name="binding_free_energy",
        tag="archive_experimental_analysis",
        branch="feature/mda-analysis-migration",
    ),
    "exposure": ArchivedAnalysisPlugin(
        name="exposure",
        tag="archive_experimental_analysis",
        branch="feature/mda-analysis-migration",
    ),
    "pa": ArchivedAnalysisPlugin(
        name="polymer_affinity",
        tag="archive_experimental_analysis",
        branch="feature/mda-analysis-migration",
    ),
    "polymer_affinity": ArchivedAnalysisPlugin(
        name="polymer_affinity",
        tag="archive_experimental_analysis",
        branch="feature/mda-analysis-migration",
    ),
    "polymer_bridging": ArchivedAnalysisPlugin(
        name="polymer_bridging",
        tag="archive_experimental_analysis",
        branch="feature/mda-analysis-migration",
    ),
    "bridging": ArchivedAnalysisPlugin(
        name="polymer_bridging",
        tag="archive_experimental_analysis",
        branch="feature/mda-analysis-migration",
    ),
}


def get_archived_analysis_plugin(name: str) -> ArchivedAnalysisPlugin | None:
    """Return archive metadata for an analysis name, if available.

    Parameters
    ----------
    name : str
        User-provided analysis plugin name.

    Returns
    -------
    ArchivedAnalysisPlugin | None
        Archive metadata when *name* is a removed plugin, otherwise ``None``.
    """
    return ARCHIVED_ANALYSIS_PLUGINS.get(name.lower())


def format_archived_analysis_message(name: str, *, context: str) -> str:
    """Build a reusable archived-plugin diagnostic.

    Parameters
    ----------
    name : str
        User-provided analysis plugin name.
    context : str
        Configuration or CLI context where the removed plugin was requested.

    Returns
    -------
    str
        User-facing message explaining where to recover the archived code.
    """
    archived = get_archived_analysis_plugin(name)
    if archived is None:
        return f"Analysis plugin '{name}' is not archived."

    return (
        f"Analysis plugin '{name}' in {context} has been archived and is no longer "
        "shipped as an active PolyzyMD analysis. To recover the archived implementation, "
        f"use git tag '{archived.tag}' from branch '{archived.branch}'. "
        "Remove this entry or choose a currently available analysis to continue."
    )
