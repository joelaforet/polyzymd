"""Per-frame distances between labelled pairs of positions.

Both distance-based plugins measure the same thing, so the measurement lives
here and neither plugin owns a copy. A pair names two selections written in the
extended syntax of :mod:`polyzymd.analyses.shared.selections`, so an endpoint
can be one atom, the midpoint of several atoms, or the centre of mass of a
group. Distances are taken with the minimum-image convention when the timestep
carries a valid box.

Coordinates are read as the trajectory stores them. Nothing is aligned first,
because a distance is invariant under rigid-body motion and rotating the
coordinates while keeping the original box vectors corrupts the minimum image.

References
----------
Michaud-Agrawal, N., Denning, E. J., Woolf, T. B. & Beckstein, O. (2011).
MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
*Journal of Computational Chemistry*, 32(10), 2319-2327. doi:10.1002/jcc.21787
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np
from numpy.typing import NDArray
from pydantic import BaseModel, Field

from polyzymd.analyses.contract import iter_frames
from polyzymd.analyses.exceptions import SelectionError

if TYPE_CHECKING:
    from polyzymd.analyses.mda.frame_selection import FrameSelection

LOGGER = logging.getLogger(__name__)

__all__ = ["PairSelection", "pair_distance_matrix"]


class PairSelection(BaseModel):
    """One labelled pair of positions to measure between.

    Parameters
    ----------
    label : str
        Name the pair is reported under, for example ``"Ser77-His156"``.
    selection_a, selection_b : str
        Selections for the two endpoints, in the extended syntax of
        :mod:`polyzymd.analyses.shared.selections`.
    """

    label: str = Field(min_length=1, description="Name the pair is reported under")
    selection_a: str = Field(min_length=1, description="First endpoint selection")
    selection_b: str = Field(min_length=1, description="Second endpoint selection")


def pair_distance_matrix(
    universe: Any,
    frames: FrameSelection,
    pairs: Sequence[PairSelection],
    *,
    use_pbc: bool = True,
    notes: list[str] | None = None,
) -> NDArray[np.float64]:
    """Measure every pair on every production frame.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe loaded by the framework.
    frames : FrameSelection
        Production window resolved by the framework.
    pairs : sequence of PairSelection
        Pairs to measure, in report order.
    use_pbc : bool, optional
        Take minimum-image distances when the timestep has a valid box, by
        default True. A frame without one is measured without periodicity.
    notes : list of str or None, optional
        List that messages about the measurement are appended to, such as an
        endpoint that spans several chains or a frame with no usable box. A
        plugin passes one in and puts it under the ``warnings`` key of its
        observables' metadata, so the message reaches the replicate artifact
        instead of only the log.

    Returns
    -------
    numpy.ndarray
        Distances in angstrom with shape ``(n_pairs, n_frames)``.

    Raises
    ------
    SelectionError
        If a selection matches no atoms, or matches several atoms without
        saying which point of the group is meant.
    """
    from MDAnalysis.lib.distances import calc_bonds

    from polyzymd.analyses.shared.selections import get_position

    collected: list[str] = [] if notes is None else notes
    resolved = [
        (
            _resolve(universe, pair.selection_a, pair.label, collected),
            _resolve(universe, pair.selection_b, pair.label, collected),
        )
        for pair in pairs
    ]
    rows: list[NDArray[np.float64]] = []
    warned = False
    for timestep in iter_frames(universe, frames):
        box = _box(timestep) if use_pbc else None
        if use_pbc and box is None and not warned:
            warned = True
            _note(
                collected,
                "pair distances: the timestep carries no usable box, so affected frames are "
                "measured without the minimum-image convention",
            )
        positions_a = np.asarray(
            [get_position(atoms, mode) for (atoms, mode), _ in resolved], dtype=np.float64
        )
        positions_b = np.asarray(
            [get_position(atoms, mode) for _, (atoms, mode) in resolved], dtype=np.float64
        )
        rows.append(calc_bonds(positions_a, positions_b, box=box).astype(np.float64))
    if not rows:
        raise SelectionError("pair distances: the production window selected no frames")
    return np.asarray(rows, dtype=np.float64).T


def _resolve(universe: Any, selection: str, label: str, notes: list[str]) -> tuple[Any, Any]:
    """Atom group and position mode for one endpoint, or a diagnostic error.

    A selection that spans several chains is measured, not refused, because
    residue numbers restart per chain and a midpoint or centre of mass then
    silently averages over the copies. The note says so.
    """
    from polyzymd.analyses.shared.diagnostics import (
        get_selection_diagnostics,
        warn_if_multi_chain_selection,
    )
    from polyzymd.analyses.shared.selections import SelectionMode, parse_selection_string

    parsed = parse_selection_string(selection)
    atoms = universe.select_atoms(parsed.selection)
    if len(atoms) == 0:
        raise SelectionError(
            f"selection {selection!r} matched no atoms.\n\n"
            f"{get_selection_diagnostics(universe, selection)}"
        )
    if parsed.mode == SelectionMode.SINGLE and len(atoms) > 1:
        raise SelectionError(
            f"selection {selection!r} matched {len(atoms)} atoms, and a pair endpoint is one "
            "point. Wrap it in midpoint(...) or com(...) to reduce the group to one position."
        )
    if warn_if_multi_chain_selection(atoms, selection, f"for distance pair {label!r}"):
        notes.append(
            f"selection {selection!r} of pair {label!r} matched atoms from several chains, so "
            "the measured point averages over the copies. Restrict it with 'protein and' or "
            "'chainid A and'."
        )
    return atoms, parsed.mode


def _note(notes: list[str], message: str) -> None:
    """Record a measurement message once, and log it."""
    if message not in notes:
        notes.append(message)
    LOGGER.warning("%s", message)


def _box(timestep: Any) -> NDArray[np.float32] | None:
    """Unit cell of one timestep, or None when it is missing or degenerate."""
    dimensions = getattr(timestep, "dimensions", None)
    if dimensions is None:
        return None
    box = np.asarray(dimensions, dtype=np.float32)
    if box.shape[0] < 6 or not np.all(np.isfinite(box[:6])) or np.any(box[:3] <= 0):
        return None
    return box[:6]
