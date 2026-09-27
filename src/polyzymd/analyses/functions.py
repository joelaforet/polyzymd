"""Per-frame measurements shipped with PolyzyMD, as plain functions.

Each function takes MDAnalysis ``AtomGroup`` arguments positioned at one frame
and returns one number, so it runs through
:meth:`polyzymd.analyses.study.Study.timeseries` like any function you write.
"""

from __future__ import annotations

import warnings
from typing import Any


def radius_of_gyration(atoms: Any) -> float:
    """Return the mass-weighted radius of gyration of ``atoms`` at the current frame.

    This calls ``AtomGroup.radius_of_gyration()``, which weights each atom by
    its mass, uses the coordinates as loaded without unwrapping molecules
    split across periodic boundaries, and returns Å.

    Parameters
    ----------
    atoms : MDAnalysis.core.groups.AtomGroup
        Atoms to measure.

    Returns
    -------
    float
        Radius of gyration in Å.
    """
    return float(atoms.radius_of_gyration())


def rmsd(atoms: Any, reference: Any) -> float:
    """Return the RMSD of ``atoms`` from ``reference`` after optimal superposition.

    This calls ``MDAnalysis.analysis.rms.rmsd`` with ``center=True`` and
    ``superposition=True``, which moves both sets of coordinates to their
    centres of geometry, rotates ``atoms`` onto ``reference`` and returns the
    square root of the mean squared distance per atom, in Å, with every atom
    weighted equally. This is what ``MDAnalysis.analysis.rms.RMSD`` computes
    for one selection, as the legacy rmsd plugin did.

    Parameters
    ----------
    atoms : MDAnalysis.core.groups.AtomGroup
        Atoms to measure.
    reference : MDAnalysis.core.groups.AtomGroup
        The same number of atoms at their reference positions, usually from
        :func:`polyzymd.analyses.reference.reference`.

    Returns
    -------
    float
        RMSD in Å.
    """
    from MDAnalysis.analysis.rms import rmsd as _rmsd

    return float(_rmsd(atoms.positions, reference.positions, center=True, superposition=True))


def pair_distance(
    atoms_a: Any, atoms_b: Any, mode_a: str = "single", mode_b: str = "single", pbc: bool = True
) -> float:
    """Return the distance in Å between one point of ``atoms_a`` and one of ``atoms_b``.

    Each point is the position of the only atom for mode ``"single"``, the
    center of geometry for ``"midpoint"`` or ``"centroid"``, and the center
    of mass for ``"com"``, from
    :func:`polyzymd.analyses.shared.selections.get_position`, as the legacy
    distances plugin wrote them as ``midpoint(...)`` and ``com(...)``. The
    distance comes from ``MDAnalysis.lib.distances.calc_bonds``, with the
    minimum image of the current box when ``pbc`` is true and the box has
    finite, positive lengths, and without periodic images otherwise. A
    frame with ``pbc`` true and no valid box raises a warning, as the legacy
    plugin did.

    Parameters
    ----------
    atoms_a, atoms_b : MDAnalysis.core.groups.AtomGroup
        Atoms of the two points.
    mode_a, mode_b : {"single", "midpoint", "centroid", "com"}, optional
        How each point is taken from its atoms.
    pbc : bool, optional
        Use the minimum image of the current box.

    Returns
    -------
    float
        Distance in Å.

    Raises
    ------
    ValueError
        If mode ``"single"`` is given more than one atom.
    """
    import numpy as np
    from MDAnalysis.lib.distances import calc_bonds

    from polyzymd.analyses.shared.selections import SelectionMode, get_position

    box = atoms_a.dimensions if pbc else None
    if box is not None:
        box = np.asarray(box, dtype=np.float32)[:6]
        if not np.all(np.isfinite(box)) or np.any(box[:3] <= 0):
            box = None
    if pbc and box is None:
        warnings.warn(
            "pair_distance: the frame has no valid box, so its distance uses no periodic image.",
            stacklevel=2,
        )
    a = get_position(atoms_a, SelectionMode(mode_a))[np.newaxis]
    b = get_position(atoms_b, SelectionMode(mode_b))[np.newaxis]
    return float(calc_bonds(a, b, box=box)[0])


def all_below(*distances: Any, thresholds: Any) -> Any:
    """Return 1 for each frame in which every distance is below its threshold, else 0.

    Parameters
    ----------
    *distances : numpy.ndarray
        One series per pair, with one distance per frame.
    thresholds : sequence of float
        One threshold per series, in the same unit. A distance counts as
        below only when it is strictly less than its threshold.

    Returns
    -------
    numpy.ndarray
        1.0 or 0.0 per frame, for :meth:`~polyzymd.analyses.timeseries.Timeseries.transform`.
    """
    import numpy as np

    below = [np.asarray(d) < t for d, t in zip(distances, thresholds, strict=True)]
    return np.logical_and.reduce(below).astype(np.float64)


def rmsf(atoms: Any, fit: Any, reference: Any, frames: Any, about_reference: bool = False) -> Any:
    """Return the RMSF in Å of each residue of ``atoms`` over ``frames``.

    Every frame is superposed on the reference by its ``fit`` atoms with
    :func:`~polyzymd.analyses.reference.superpose`, which rotates as
    ``MDAnalysis.analysis.align.AlignTraj`` does, without moving the
    trajectory. ``MDAnalysis.analysis.rms.RMSF`` then gives each atom's
    root mean square fluctuation about its mean position over the frames.
    With ``about_reference`` true, each atom's root mean square deviation
    from its reference position is taken instead, as the legacy rmsf plugin
    did in external mode. Each residue's value is the mean over its atoms in
    ``atoms``, in the order of ``atoms.residues``.

    Parameters
    ----------
    atoms : MDAnalysis.core.groups.AtomGroup
        Atoms to measure.
    fit : MDAnalysis.core.groups.AtomGroup
        Atoms that are superposed on the reference.
    reference : MDAnalysis.core.groups.AtomGroup
        Reference positions of ``atoms | fit`` in index order, from
        :func:`polyzymd.analyses.reference.reference` with the selection
        ``"(<atoms>) or (<fit>)"``.
    frames : numpy.ndarray
        Trajectory frame indices to use.
    about_reference : bool, optional
        Measure deviations from the reference positions instead of
        fluctuations about the mean.

    Returns
    -------
    numpy.ndarray
        One value per residue of ``atoms``.

    Raises
    ------
    ProtocolError
        If ``reference`` does not hold one position per atom of ``atoms | fit``.
    """
    import MDAnalysis as mda
    import numpy as np
    from MDAnalysis.analysis.rms import RMSF
    from MDAnalysis.coordinates.memory import MemoryReader

    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.reference import superpose

    group = atoms | fit
    if len(reference) != len(group):
        raise ProtocolError(
            f"rmsf: the reference has {len(reference)} atoms and the selection and fit atoms "
            f"together have {len(group)}.",
            hint="Build the reference from the selection '(<selection>) or (<alignment>)'.",
        )
    where_fit = np.searchsorted(group.indices, fit.indices)
    where = np.searchsorted(group.indices, atoms.indices)
    coordinates = np.array([group.positions for _ in atoms.universe.trajectory[frames]], float)
    target = reference.positions.astype(float)
    moved = superpose(coordinates, where_fit, target[where_fit])[:, where]
    if about_reference:
        per_atom = np.sqrt(np.mean(np.sum((moved - target[where]) ** 2, axis=2), axis=0))
    else:
        copy = mda.Merge(atoms)
        copy.load_new(moved.astype(np.float32), format=MemoryReader)
        per_atom = RMSF(copy.atoms).run().results.rmsf
    _, residue = np.unique(atoms.resindices, return_inverse=True)
    return np.bincount(residue, weights=per_atom) / np.bincount(residue)
