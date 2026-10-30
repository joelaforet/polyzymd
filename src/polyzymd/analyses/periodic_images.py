"""How close a solute comes to its own periodic image over a trajectory.

A cube or a rhombic dodecahedron is sized so the solute stays clear of its
images in any orientation. A rectangular box is sized for the solute's starting
orientation only, so a long protein that turns during the run can come within
the nonbonded cutoff of its own image. :func:`solute_image_distances` measures
that distance frame by frame; :func:`image_distance_findings` turns it into
plain warnings.
"""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import MDAnalysis as mda


def _lattice_vectors(dimensions: np.ndarray) -> list[np.ndarray]:
    """The 26 non-zero lattice translations of a cell, in Angstrom."""
    from MDAnalysis.lib.mdamath import triclinic_vectors

    box = triclinic_vectors(dimensions)
    return [
        i * box[0] + j * box[1] + k * box[2]
        for i, j, k in itertools.product((-1, 0, 1), repeat=3)
        if (i, j, k) != (0, 0, 0)
    ]


def solute_image_distances(
    universe: "mda.Universe",
    selection: str = "segid PROT and not name H*",
    *,
    cutoff_nm: float = 2.0,
    step: int = 1,
) -> list[tuple[int, float, float | None]]:
    """Closest distance from the selected atoms to their own periodic images, per frame.

    The selection is treated as one piece. Its coordinates in the first frame
    are taken as whole, which they are in a PolyzyMD build (the solutes are
    placed whole in the centre of the cell), and later frames are unwrapped
    with MDAnalysis' :class:`~MDAnalysis.transformations.nojump.NoJump`, so
    no atom of the selection jumps across the box between frames. Chains in
    contact therefore stay in contact rather than being counted as
    self-image contacts, and no bonds or masses are needed. Every frame is
    read in order for the unwrapping; ``step`` only thins the frames measured.
    The universe's own trajectory is not changed.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Topology and trajectory, with box dimensions on every frame; the
        selection must be whole in the first frame.
    selection : str
        The solute; by default the protein segment's heavy atoms.
    cutoff_nm : float
        Distances are resolved only below this value (it keeps the search fast);
        a frame with no image closer is reported as None.
    step : int
        Analyse every ``step``-th frame.

    Returns
    -------
    list of (frame, time_ps, distance_nm or None)
    """
    from MDAnalysis.coordinates.timestep import Timestep
    from MDAnalysis.lib.distances import capped_distance
    from MDAnalysis.transformations.nojump import NoJump

    atoms = universe.select_atoms(selection)
    if len(atoms) == 0:
        raise ValueError(f"the selection {selection!r} matches no atoms")
    nojump = NoJump()
    results: list[tuple[int, float, float | None]] = []
    for index, ts in enumerate(universe.trajectory):
        # NoJump on a copy of the selection only: the universe is left as it was.
        piece = Timestep(len(atoms))
        piece.frame = index
        piece.dimensions = ts.dimensions
        piece.positions = atoms.positions
        positions = np.asarray(nojump(piece).positions, dtype=np.float32)
        if index % step:
            continue
        best: float | None = None
        for vector in _lattice_vectors(ts.dimensions):
            _, distances = capped_distance(
                positions,
                (positions + vector).astype(np.float32),
                cutoff_nm * 10.0,
                return_distances=True,
            )
            if len(distances):
                low = float(distances.min()) / 10.0
                best = low if best is None else min(best, low)
        results.append((int(ts.frame), float(ts.time), best))
    return results


def image_distance_findings(
    distances: list[tuple[int, float, float | None]], nonbonded_cutoff_nm: float
) -> list[str]:
    """Warnings for the frames in which the solute comes within the cutoff of its image."""
    close = [(frame, time, d) for frame, time, d in distances if d is not None]
    below = [item for item in close if item[2] < nonbonded_cutoff_nm]
    if not below:
        return []
    frame, time, distance = min(below, key=lambda item: item[2])
    return [
        f"the solute comes within {distance:.2f} nm of its own periodic image (frame {frame}, "
        f"{time:.0f} ps), inside the {nonbonded_cutoff_nm:.2f} nm nonbonded cutoff, in "
        f"{len(below)} of {len(distances)} analysed frames: it interacts with its own image. "
        "Use a larger padding, or a rhombic dodecahedron, which is sized for any orientation."
    ]
