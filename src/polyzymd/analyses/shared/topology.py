"""Topology bond checks shared by fragment-based observables.

Fragment-based observables (Rg in fragment mode, contacts polymer chain
identity) are only meaningful when the atoms they measure are connected by
bonds. Two failure modes matter. A topology can carry no bonds at all, which is
what MDAnalysis produces for a PDB whose atom serials run above 99999, and a
topology can carry bonds for some molecules and not others, which is what
happens when only the standard residues have usable CONECT records. The second
case does not raise inside MDAnalysis; it silently yields one singleton
fragment per atom. Both are checked here against the selection actually being
measured, not against the universe as a whole.

References
----------
Michaud-Agrawal, N., Denning, E. J., Woolf, T. B., and Beckstein, O. (2011).
MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
Journal of Computational Chemistry, 32(10), 2319-2327. doi:10.1002/jcc.21787
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup
    from MDAnalysis.core.universe import Universe

__all__ = [
    "SINGLETON_ATOM_FRACTION_LIMIT",
    "require_topology_bonds",
    "topology_bond_source",
]

#: Largest share of a selection that may sit in single-atom fragments before the
#: selection counts as unbonded. A real polymer selection has none. A protein
#: whose CONECT records MDAnalysis refused has almost all of its atoms here, and
#: a mixed protein-and-polymer selection in that state still reports a quarter of
#: its atoms as singletons while producing fragment Rg values near zero, so an
#: all-or-nothing test would let it through.
SINGLETON_ATOM_FRACTION_LIMIT = 0.05


def _no_data_error() -> tuple[type[BaseException], ...]:
    """Return the exception classes MDAnalysis raises for absent bond topology.

    Returns
    -------
    tuple[type[BaseException], ...]
        ``NoDataError`` when MDAnalysis is importable, plus ``AttributeError``
        for objects that expose no bond attributes at all.
    """

    try:
        from MDAnalysis.exceptions import NoDataError
    except ImportError:  # pragma: no cover - MDAnalysis is a hard dependency here
        return (AttributeError,)
    return (NoDataError, AttributeError)


def _missing_bond_reason(atom_group: AtomGroup) -> str | None:
    """Explain why a selection has no usable bond topology.

    Parameters
    ----------
    atom_group : AtomGroup
        Selection that a fragment-based observable is about to measure.

    Returns
    -------
    str or None
        A sentence naming the problem, or ``None`` when the selection has
        usable bonds.
    """

    errors = _no_data_error()
    try:
        n_bonds = len(atom_group.bonds)
    except errors:
        return "The topology carries no bond information."
    if n_bonds == 0:
        return "No bond in the topology touches the selected atoms."
    try:
        fragments = list(atom_group.fragments)
    except errors:
        return "The topology carries no bond information."
    if not fragments:
        return "The selection resolved to no topology fragments."
    n_atoms = len(atom_group)
    singleton_atoms = sum(len(fragment) for fragment in fragments if len(fragment) == 1)
    if n_atoms > 1 and singleton_atoms > n_atoms * SINGLETON_ATOM_FRACTION_LIMIT:
        percent = 100.0 * singleton_atoms / n_atoms
        return (
            f"{singleton_atoms} of the {n_atoms} selected atoms are single-atom fragments "
            f"({percent:.0f} percent), so most of the selection carries no bonds."
        )
    return None


def require_topology_bonds(
    atom_group: AtomGroup,
    *,
    context: str,
    topology_path: Path | str | None = None,
    allow_fallback: bool = False,
) -> tuple[list[Any], str | None]:
    """Return the bonded fragments of a selection, or fail loudly.

    Parameters
    ----------
    atom_group : AtomGroup
        Selection whose fragments are needed.
    context : str
        What needed the fragments, for example ``"Rg run 'polymer' in fragment
        mode"``. It opens the error and warning text.
    topology_path : Path or str or None, optional
        Topology file named in the error message.
    allow_fallback : bool, optional
        Return the whole selection as one fragment instead of raising, by
        default False.

    Returns
    -------
    tuple[list[Any], str or None]
        The fragments, and a warning message when the fallback was used.

    Raises
    ------
    TopologyBondsMissingError
        If the selection has no usable bonds and ``allow_fallback`` is False.
    """

    from polyzymd.analyses.exceptions import TopologyBondsMissingError

    reason = _missing_bond_reason(atom_group)
    if reason is None:
        return list(atom_group.fragments), None
    if not allow_fallback:
        raise TopologyBondsMissingError(
            context=context,
            n_atoms=len(atom_group),
            topology=topology_path,
            detail=reason,
        )
    return [atom_group], f"{context}: {reason} Treating the whole selection as one fragment."


def topology_bond_source(universe: Universe) -> tuple[bool, str]:
    """Report whether a universe carries bonds and where they came from.

    Parameters
    ----------
    universe : Universe
        Loaded MDAnalysis universe. Objects without a ``bonds`` attribute are
        reported as carrying no bonds.

    Returns
    -------
    tuple[bool, str]
        Whether bonds are present, and the bond source as ``"conect"`` (read
        from the topology file), ``"guessed"`` (inferred by MDAnalysis), or
        ``"none"``.

    Notes
    -----
    This describes the topology as a whole. It says nothing about whether the
    atoms a given analysis selects are bonded; use
    :func:`require_topology_bonds` for that.
    """

    bonds = getattr(universe, "bonds", None)
    try:
        n_bonds = len(bonds) if bonds is not None else 0
    except TypeError:
        n_bonds = 0
    if n_bonds == 0:
        return False, "none"
    # MDAnalysis exposes the per-bond guessed flags only as TopologyGroup._guessed.
    # There is no public accessor on the group, so reach it through the public
    # AtomGroup.bonds rather than through universe._topology.
    guessed = getattr(bonds, "_guessed", None)
    if guessed is None:
        return True, "conect"
    try:
        all_guessed = bool(np.all(np.asarray(guessed, dtype=bool)))
    except (TypeError, ValueError):
        all_guessed = False
    return True, "guessed" if all_guessed else "conect"
