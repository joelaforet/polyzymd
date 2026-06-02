"""Condition filtering helpers for the contacts analysis facade."""

from __future__ import annotations

import logging
from typing import Any

from pydantic import BaseModel

from polyzymd.analyses.base import Condition

logger = logging.getLogger("polyzymd.analyses.contacts")


class PolymerTopologyUnavailableError(OSError):
    """Raised when no replicate topology can be inspected for polymer atoms."""


def filter_conditions(
    analysis: Any, conditions: list[Condition], settings: BaseModel | None
) -> list[Condition]:
    """Filter conditions to only those with polymer atoms.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance that provides settings and selection helpers.
    conditions : list[Condition]
        Conditions from the comparison configuration.
    settings : BaseModel or None
        Resolved plugin settings from the orchestrator.

    Returns
    -------
    list[Condition]
        Conditions with polymer atoms, plus conditions that could not be checked.
    """

    resolved = settings if isinstance(settings, analysis.Settings) else analysis.Settings()
    valid: list[Condition] = []

    for cond in conditions:
        try:
            if analysis._condition_has_polymer(cond, resolved):
                valid.append(cond)
            else:
                logger.info(
                    f"  Excluding '{cond.label}': no polymer atoms found with selection "
                    f"'{analysis._effective_polymer_selection(resolved)}'"
                )
        except (PolymerTopologyUnavailableError, AttributeError, KeyError, OSError) as e:
            logger.warning(f"  Error checking condition '{cond.label}': {e} — including anyway")
            valid.append(cond)

    return valid


def condition_has_polymer(analysis: Any, cond: Condition, settings: BaseModel) -> bool:
    """Check whether a condition has polymer atoms.

    Detection uses topology inspection via MDAnalysis with the effective polymer
    selection. Stale contacts caches are not used as polymer evidence.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance that provides selection helpers.
    cond : Condition
        Condition to check.
    settings : BaseModel
        Resolved contacts settings with polymer query filters.

    Returns
    -------
    bool
        ``True`` when polymer atoms are detected in any replicate topology.
    """

    inspected_any = False

    for rep in cond.replicates:
        run_dir = cond.sim_config.get_working_directory(rep)

        try:
            from polyzymd.analyses.shared.loader import TrajectoryLoader

            loader = TrajectoryLoader(cond.sim_config)
            topo_path = loader.find_topology(run_dir)
        except (FileNotFoundError, ImportError) as e:
            logger.debug(f"  Could not locate topology for {cond.label} rep {rep}: {e}")
            continue

        try:
            import MDAnalysis as mda

            universe = mda.Universe(str(topo_path))
            polymer_atoms = universe.select_atoms(analysis._effective_polymer_selection(settings))
            inspected_any = True
            if len(polymer_atoms) > 0:
                logger.debug(f"  {cond.label} rep {rep}: {len(polymer_atoms)} polymer atoms")
                return True
            logger.debug(f"  {cond.label} rep {rep}: 0 polymer atoms")
        except (AttributeError, ValueError, KeyError, OSError, ImportError) as e:
            logger.warning(f"  Error checking {cond.label} rep {rep}: {e}")
            continue

    if not inspected_any:
        raise PolymerTopologyUnavailableError(
            f"no replicate topology could be inspected for condition '{cond.label}'"
        )

    return False
