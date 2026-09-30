"""Measurements shipped with PolyzyMD, as plain functions.

The per-frame functions take MDAnalysis ``AtomGroup`` arguments positioned at
one frame and return one number, so they run through
:meth:`polyzymd.analyses.study.Study.timeseries` like any function you write.
The per-replicate functions (:func:`rmsf`, :func:`rms_deviation`,
:func:`rms_decomposition`, :func:`residue_sasa`, :func:`dssp_occupancy` and
:func:`residue_contacts` and :func:`residue_occlusion`) also take the production
frame indices and return one value per residue, and run through
:meth:`polyzymd.analyses.study.Study.per_replicate`.
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


def _superposed_deviations(atoms: Any, fit: Any, reference: Any, frames: Any) -> tuple:
    """Superpose every frame on the reference and return the per-atom deviation, RMSF and offset.

    Every frame is superposed on the reference by its ``fit`` atoms with
    :func:`~polyzymd.analyses.reference.superpose`, which rotates as
    ``MDAnalysis.analysis.align.AlignTraj`` does, without moving the
    trajectory. The superposed positions of ``atoms`` are then read once to
    give, per atom, the root mean square deviation from the reference
    position, ``MDAnalysis.analysis.rms.RMSF`` (the fluctuation about the
    mean position), and the distance of the mean position from the
    reference position. The squared deviation equals the squared RMSF plus
    the squared offset for every atom.
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
    deviation = np.sqrt(np.mean(np.sum((moved - target[where]) ** 2, axis=2), axis=0))
    copy = mda.Merge(atoms)
    copy.load_new(moved.astype(np.float32), format=MemoryReader)
    fluctuation = RMSF(copy.atoms).run().results.rmsf
    offset = np.linalg.norm(moved.mean(axis=0) - target[where], axis=1)
    return deviation, fluctuation, offset


def _per_residue(atoms: Any, per_atom: Any) -> Any:
    """Average per-atom values over the atoms of each residue, in the order of ``atoms.residues``."""
    import numpy as np

    _, residue = np.unique(atoms.resindices, return_inverse=True)
    return np.bincount(residue, weights=per_atom) / np.bincount(residue)


def rmsf(atoms: Any, fit: Any, reference: Any, frames: Any) -> Any:
    """Return the RMSF in Å of each residue of ``atoms`` over ``frames``.

    Every frame is superposed on the reference by its ``fit`` atoms, and
    ``MDAnalysis.analysis.rms.RMSF`` gives each atom's root mean square
    fluctuation about its mean position over the frames, as ``gmx rmsf -o``
    does. The reference only decides what the frames are superposed on.
    Each residue's value is the mean over its atoms in ``atoms``, in the
    order of ``atoms.residues``.

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

    Returns
    -------
    numpy.ndarray
        One value per residue of ``atoms``.

    Raises
    ------
    ProtocolError
        If ``reference`` does not hold one position per atom of ``atoms | fit``.
    """
    return _per_residue(atoms, _superposed_deviations(atoms, fit, reference, frames)[1])


def rms_deviation(atoms: Any, fit: Any, reference: Any, frames: Any) -> Any:
    """Return each residue's root mean square deviation in Å from the reference over ``frames``.

    After the superposition of :func:`rmsf`, each atom's value is
    ``sqrt(<|x(t) - x_ref|^2>)``, as ``gmx rmsf -od`` gives, and the legacy
    rmsf plugin gave in external mode. Each residue's value is the mean over
    its atoms. The arguments are those of :func:`rmsf`.
    """
    return _per_residue(atoms, _superposed_deviations(atoms, fit, reference, frames)[0])


#: Names of the per-residue means that :func:`rms_decomposition` returns first.
RMS_PARTS = ("rms_deviation", "rmsf", "offset")

#: Names of the per-residue mean squares that :func:`rms_decomposition` returns after them.
MS_PARTS = ("ms_deviation", "msf", "ms_offset")


def rms_decomposition(atoms: Any, fit: Any, reference: Any, frames: Any) -> Any:
    """Return each residue's RMS deviation, RMSF and offset, and their mean squares, in one pass.

    The first three rows, named in :data:`RMS_PARTS`, are the values of
    :func:`rms_deviation` and :func:`rmsf`, and the offset, the distance in
    Å of each atom's mean position from its reference position, each
    averaged over the residue's atoms. The last three, named in
    :data:`MS_PARTS`, are the means over the residue's atoms of the squares
    of the same per-atom values, in Å². For every atom the squared deviation
    is the squared RMSF plus the squared offset, so ``ms_deviation`` equals
    ``msf + ms_offset`` for every residue. The arguments are those of
    :func:`rmsf`.

    Returns
    -------
    numpy.ndarray
        Shape ``(6, n_residues)``.
    """
    import numpy as np

    parts = _superposed_deviations(atoms, fit, reference, frames)
    means = [_per_residue(atoms, values) for values in parts]
    return np.vstack(means + [_per_residue(atoms, values**2) for values in parts])


#: Probe radius in nm and sphere point count of :func:`sasa` and :func:`residue_sasa`.
SASA_PROBE_RADIUS_NM = 0.14
SASA_SPHERE_POINTS = 960

#: MDTraj topologies of SASA contexts and target positions, keyed by universe and atom indices.
_SASA_TOPOLOGIES: dict[tuple[int, bytes, bytes], tuple[Any, Any, Any]] = {}


def _sasa_topology(target: Any, context: Any) -> tuple[Any, Any]:
    """Return the MDTraj topology of ``context`` and the positions of ``target`` in it.

    The topology holds one MDTraj atom per atom of ``context``, in index
    order, with its name, residue and element from the MDAnalysis universe;
    MDTraj's Shrake-Rupley code takes each atom's radius from its element.
    Residues are split by MDAnalysis residue index. It is built once per
    universe, context and target and reused for every frame.
    """
    import mdtraj as md
    import numpy as np

    from polyzymd.analyses.exceptions import ProtocolError

    key = (id(context.universe), context.indices.tobytes(), target.indices.tobytes())
    cached = _SASA_TOPOLOGIES.get(key)
    if cached is not None and cached[0] is context.universe:
        return cached[1], cached[2]
    ordered = bool(np.all(np.diff(context.indices) > 0))
    where = np.searchsorted(context.indices, target.indices)
    inside = (
        ordered
        and len(target) > 0
        and bool(np.all(where < len(context)))
        and np.array_equal(context.indices[np.minimum(where, len(context) - 1)], target.indices)
    )
    if not inside:
        raise ProtocolError(
            f"sasa: the target has {len(target)} atoms and must be a non-empty part of the "
            f"context, which has {len(context)} atoms in index order.",
            hint="Choose a context selection that contains every target atom, such as "
            "'protein or resname SBM EGM' for the target 'protein'.",
        )
    if not hasattr(context, "elements"):
        raise ProtocolError(
            "sasa: the universe has no element for its atoms, and the atomic radii come from them.",
            hint="Load the replicate with PolyzyMD, which fills in elements from atom types or names.",
        )
    topology = md.Topology()
    chain = topology.add_chain()
    residues: dict[int, Any] = {}
    for atom in context:
        residue = residues.get(atom.resindex)
        if residue is None:
            residue = residues[atom.resindex] = topology.add_residue(
                str(atom.resname), chain, resSeq=int(atom.resid)
            )
        symbol = str(atom.element).strip().capitalize()
        try:
            element = md.element.get_by_symbol(symbol)
        except KeyError as exc:
            raise ProtocolError(
                f"sasa: atom {atom.index} ({atom.name}) has element {atom.element!r}, which "
                "MDTraj does not know.",
                hint="Check the topology's element column or the atom types.",
            ) from exc
        topology.add_atom(str(atom.name), element, residue)
    if len(_SASA_TOPOLOGIES) > 32:
        _SASA_TOPOLOGIES.clear()
    _SASA_TOPOLOGIES[key] = (context.universe, topology, where)
    return topology, where


def _atom_sasa(
    target: Any, context: Any, positions: Any, probe_radius_nm: float, n_sphere_points: int
) -> Any:
    """Return the SASA in Å² of each ``target`` atom in each frame of ``positions`` (Å, of ``context``)."""
    import mdtraj as md
    import numpy as np
    from scipy.spatial import cKDTree

    from polyzymd.analyses.exceptions import ProtocolError

    topology, where = _sasa_topology(target, context)
    xyz = np.asarray(positions, dtype=np.float32) / 10.0
    # MDTraj's Shrake-Rupley code ends the process, without raising, when two
    # atoms are closer than 1e-5 nm, so such frames are refused here first.
    for frame in xyz:
        close = cKDTree(frame).query_pairs(1e-5, output_type="ndarray")
        if len(close):
            first, second = context.indices[close[0]]
            raise ProtocolError(
                f"sasa: atoms {first} and {second} are at the same position, and MDTraj's "
                "Shrake-Rupley code cannot compute SASA with two atoms on top of each other.",
                hint="Check the coordinates of those atoms, or leave one of them out of the "
                "selection.",
            )
    trajectory = md.Trajectory(xyz=xyz, topology=topology)
    atom_nm2 = md.shrake_rupley(
        trajectory, mode="atom", probe_radius=probe_radius_nm, n_sphere_points=n_sphere_points
    )
    return np.asarray(atom_nm2, dtype=np.float64)[:, where] * 100.0


def sasa(
    target: Any,
    context: Any,
    probe_radius_nm: float = SASA_PROBE_RADIUS_NM,
    n_sphere_points: int = SASA_SPHERE_POINTS,
) -> float:
    """Return the solvent-accessible surface area in Å² of ``target`` at the current frame.

    ``mdtraj.shrake_rupley`` computes the SASA of every atom of ``context``
    by the Shrake-Rupley method, with a probe of ``probe_radius_nm`` and
    ``n_sphere_points`` points per atom and MDTraj's atomic radius for each
    atom's element, and the values of the ``target`` atoms are summed. Atoms
    of ``context`` outside ``target``, such as polymer, occlude the target
    without being counted. Periodic images are not considered, so atoms only
    occlude each other within the coordinates as loaded.

    The frame is computed in a call of its own, because MDTraj 1.11.1 gives
    every frame after the first that each OpenMP thread computes in one call
    about 0.1 percent too much area; see
    :func:`residue_sasa`.

    Parameters
    ----------
    target : MDAnalysis.core.groups.AtomGroup
        Atoms whose SASA is returned; they must all be in ``context``.
    context : MDAnalysis.core.groups.AtomGroup
        Atoms present in the calculation.
    probe_radius_nm : float, optional
        Probe radius in nm, 0.14 by default.
    n_sphere_points : int, optional
        Points on each atom's sphere, 960 by default.

    Returns
    -------
    float
        SASA of ``target`` in Å².
    """
    positions = context.positions[None]
    return float(_atom_sasa(target, context, positions, probe_radius_nm, n_sphere_points).sum())


def residue_sasa(
    target: Any,
    context: Any,
    frames: Any,
    probe_radius_nm: float = SASA_PROBE_RADIUS_NM,
    n_sphere_points: int = SASA_SPHERE_POINTS,
) -> Any:
    """Return each ``target`` residue's SASA in Å², averaged over ``frames``.

    Each frame's atom SASA is computed as in :func:`sasa` and summed over
    the atoms of each residue of ``target``. The result has one value per
    residue, in the order of ``target.residues``.

    Every frame goes to ``mdtraj.shrake_rupley`` in a call of its own. In
    MDTraj 1.11.1, a frame that follows another on the same OpenMP thread in
    one call comes out about 0.1 percent larger than the same coordinates
    alone, and than an
    independent Shrake-Rupley calculation, so frames are never batched.
    """
    import numpy as np

    residue = np.unique(target.resindices, return_inverse=True)[1]
    total = np.zeros(len(target.residues))
    for _ in context.universe.trajectory[frames]:
        atom = _atom_sasa(
            target, context, context.positions[None], probe_radius_nm, n_sphere_points
        )
        total += np.bincount(residue, weights=atom[0], minlength=len(total))
    return total / len(frames)


#: Names of the eight DSSP classes, with their codes in
#: ``mdtraj.compute_dssp(simplified=False)``. ``unassigned`` is MDTraj's ``"NA"``,
#: given to a residue it cannot assign.
DSSP_CLASSES = {
    "alpha_helix": "H",
    "3_10_helix": "G",
    "pi_helix": "I",
    "extended_strand": "E",
    "isolated_bridge": "B",
    "turn": "T",
    "bend": "S",
    "loop": " ",
    "unassigned": "NA",
}

#: Names of the three classes of ``mdtraj.compute_dssp(simplified=True)``, with
#: their codes, and ``unassigned``.
DSSP_SIMPLIFIED = {"helix": "H", "strand": "E", "coil": "C", "unassigned": "NA"}

#: The eight classes each simplified class joins, as MDTraj translates them.
DSSP_GROUPS = {
    "helix": ("alpha_helix", "3_10_helix", "pi_helix"),
    "strand": ("extended_strand", "isolated_bridge"),
    "coil": ("turn", "bend", "loop"),
}

#: MDTraj topologies of DSSP selections, keyed by universe and atom indices.
_DSSP_TOPOLOGIES: dict[tuple[int, bytes], tuple[Any, Any]] = {}


def _dssp_topology(atoms: Any) -> Any:
    """Return an MDTraj topology of ``atoms``, which must hold whole residues.

    One MDTraj chain is made per chain ID (or segment, when the topology has
    no chain IDs; one chain when it has neither), so DSSP never pairs
    residues of different chains; each
    residue keeps its name and number, and each atom its name and its element
    from the MDAnalysis universe. It is built once per universe and selection.
    """
    import mdtraj as md

    from polyzymd.analyses.exceptions import ProtocolError

    key = (id(atoms.universe), atoms.indices.tobytes())
    cached = _DSSP_TOPOLOGIES.get(key)
    if cached is not None and cached[0] is atoms.universe:
        return cached[1]
    if len(atoms) == 0 or len(atoms.residues.atoms) != len(atoms):
        raise ProtocolError(
            f"dssp: the selection has {len(atoms)} atoms and must hold whole residues "
            f"({len(atoms.residues.atoms)} atoms in its residues).",
            hint="Select whole residues, such as 'protein'; DSSP needs every backbone atom.",
        )
    if not hasattr(atoms, "elements"):
        raise ProtocolError(
            "dssp: the universe has no element for its atoms.",
            hint="Load the replicate with PolyzyMD, which fills in elements from atom types or names.",
        )
    topology = md.Topology()
    chains: dict[str, Any] = {}
    residues: dict[int, Any] = {}
    for atom in atoms:
        chain_key = str(getattr(atom, "chainID", "") or getattr(atom, "segid", ""))
        chain = chains.get(chain_key)
        if chain is None:
            chain = chains[chain_key] = topology.add_chain()
        residue = residues.get(atom.resindex)
        if residue is None:
            residue = residues[atom.resindex] = topology.add_residue(
                str(atom.resname), chain, resSeq=int(atom.resid)
            )
        try:
            element = md.element.get_by_symbol(str(atom.element).strip().capitalize())
        except KeyError as exc:
            raise ProtocolError(
                f"dssp: atom {atom.index} ({atom.name}) has element {atom.element!r}, which "
                "MDTraj does not know.",
                hint="Check the topology's element column or the atom types.",
            ) from exc
        topology.add_atom(str(atom.name), element, residue)
    if len(_DSSP_TOPOLOGIES) > 32:
        _DSSP_TOPOLOGIES.clear()
    _DSSP_TOPOLOGIES[key] = (atoms.universe, topology)
    return topology


def dssp_occupancy(atoms: Any, frames: Any, simplified: bool = True, chunk: int = 200) -> Any:
    """Return, for each residue of ``atoms``, the fraction of ``frames`` in each DSSP class.

    ``mdtraj.compute_dssp(simplified=simplified)`` assigns every residue of
    ``atoms`` a DSSP code, or ``"NA"`` when it cannot, on every frame,
    ``chunk`` frames per call. With ``simplified``, the rows are the classes of
    :data:`DSSP_SIMPLIFIED`: helix (DSSP H, G and I), strand (E and B), coil
    (T, S and loop) and unassigned. Otherwise they are the eight classes and
    unassigned of :data:`DSSP_CLASSES`. Columns follow ``atoms.residues``, and
    each residue's row values sum to 1. The coordinates are used as loaded, so
    a protein split across a periodic boundary should be made whole first.

    Returns
    -------
    numpy.ndarray
        Shape ``(len(DSSP_SIMPLIFIED), n_residues)`` with ``simplified``,
        otherwise ``(len(DSSP_CLASSES), n_residues)``.
    """
    import mdtraj as md
    import numpy as np

    topology = _dssp_topology(atoms)
    codes = list((DSSP_SIMPLIFIED if simplified else DSSP_CLASSES).values())
    counts = np.zeros((len(codes), len(atoms.residues)))
    trajectory = atoms.universe.trajectory
    for start in range(0, len(frames), chunk):
        xyz = np.array([atoms.positions for _ in trajectory[frames[start : start + chunk]]])
        assigned = md.compute_dssp(
            md.Trajectory(xyz=xyz.astype(np.float32) / 10.0, topology=topology),
            simplified=simplified,
        )
        for row, code in enumerate(codes):
            counts[row] += (assigned == code).sum(axis=0)
    return counts / len(frames)


#: Cutoff in Å within which a protein residue and a polymer are in contact.
#: Distance in Å within which :func:`residue_contacts` counts a contact.
CONTACT_CUTOFF = 4.0


def residue_contacts(
    protein: Any,
    polymer: Any,
    frames: Any,
    cutoff: float = CONTACT_CUTOFF,
    types: Any = (),
    pbc: bool = True,
) -> Any:
    """Return, for each residue of ``protein``, the fraction of ``frames`` it touches ``polymer``.

    On every frame, ``MDAnalysis.lib.distances.capped_distance`` finds every
    pair of a ``polymer`` atom and a ``protein`` atom closer than ``cutoff``
    Å, with the minimum image of the frame's box when ``pbc`` is true and the
    box is known. A protein residue is in contact on a frame when any of its
    atoms is in such a pair. The first row is each residue's fraction of
    frames in contact with any polymer atom; then one row per residue name in
    ``types``, the fraction of frames in contact with polymer atoms of that
    residue name, such as one monomer type. The type rows do not sum to the
    first, since a residue can touch several types on one frame. Columns
    follow ``protein.residues``.

    The atoms compared are the atoms given; select heavy atoms, such as
    ``chainid A and not element H``, to leave hydrogens out.

    Returns
    -------
    numpy.ndarray
        Shape ``(1 + len(types), n_residues)``.
    """
    import numpy as np

    counts = np.zeros((1 + len(types), len(protein.residues)))
    for _, touched in _distance_frames(protein, polymer, frames, cutoff, types, pbc):
        counts += touched
    return counts / len(frames)


def _distance_frames(
    protein: Any, polymer: Any, frames: Any, cutoff: float, types: Any, pbc: bool
) -> Any:
    """Yield each frame's time in ps and contact rows of :func:`residue_contacts`.

    The rows are booleans of shape ``(1 + len(types), n_residues)``: in
    contact with any ``polymer`` atom, then with the atoms of each residue
    name in ``types``.
    """
    import numpy as np
    from MDAnalysis.lib.distances import capped_distance

    residue = np.unique(protein.resindices, return_inverse=True)[1]
    n_residues = len(protein.residues)
    type_of_atom = np.full(len(polymer), -1)
    for row, name in enumerate(types):
        type_of_atom[polymer.resnames == name] = row
    for ts in protein.universe.trajectory[frames]:
        touched = np.zeros((1 + len(types), n_residues), dtype=bool)
        box = ts.dimensions if pbc else None
        pairs = capped_distance(
            polymer.positions, protein.positions, max_cutoff=cutoff, box=box, return_distances=False
        )
        if len(pairs):
            touched[0, residue[pairs[:, 1]]] = True
            kinds = type_of_atom[pairs[:, 0]]
            known = kinds >= 0
            touched[1 + kinds[known], residue[pairs[known, 1]]] = True
        yield float(ts.time), touched


#: Relative SASA, as a fraction of the residue's maximum ASA, at or above which
#: :func:`residue_occlusion` counts a residue exposed with the protein alone.
EXPOSED_THRESHOLD = 0.2

#: Relative SASA below which :func:`residue_occlusion` counts an exposed residue
#: buried by the occluder.
BURIED_THRESHOLD = 0.2

#: Rows of :func:`residue_occlusion` before the per-type rows.
OCCLUSION_PARTS = ("contact_fraction", "exposed_fraction", "occluded_area", "exposed_area")


def _nearest_images(anchor: Any, atoms: Any, box: Any) -> Any:
    """Return the positions of ``atoms`` with each molecule moved to its image nearest ``anchor``.

    Atoms are grouped into molecules by their bonded fragment. Each molecule
    is translated as a whole by the box vector that brings its centroid to
    the minimum image of its offset from the centroid of ``anchor``, with
    ``MDAnalysis.lib.distances.minimize_vectors``. Molecules must be whole in
    the coordinates as loaded. Without a box the positions are returned as
    they are.
    """
    import numpy as np
    from MDAnalysis.exceptions import NoDataError
    from MDAnalysis.lib.distances import minimize_vectors

    from polyzymd.analyses.exceptions import ProtocolError

    positions = atoms.positions
    if box is None:
        return positions
    try:
        fragments = atoms.fragindices
    except NoDataError as exc:
        raise ProtocolError(
            "occlusion: the topology has no bonds, so polymer molecules cannot be moved whole "
            "to the image nearest the protein.",
            hint="Load a topology with bonds, or pass pbc=False (use_pbc=false) to use the "
            "coordinates as loaded.",
        ) from exc
    inverse = np.unique(fragments, return_inverse=True)[1]
    counts = np.bincount(inverse)
    centroids = (
        np.stack([np.bincount(inverse, weights=positions[:, axis]) for axis in range(3)], axis=1)
        / counts[:, None]
    )
    offsets = (centroids - anchor.positions.mean(axis=0)).astype(np.float32)
    shifts = minimize_vectors(offsets, box) - offsets
    return positions + shifts[inverse]


def residue_occlusion(
    protein: Any,
    occluder: Any,
    frames: Any,
    exposed_threshold: float = EXPOSED_THRESHOLD,
    buried_threshold: float = BURIED_THRESHOLD,
    types: Any = (),
    max_asa: str = "theoretical",
    pbc: bool = True,
    probe_radius_nm: float = SASA_PROBE_RADIUS_NM,
    n_sphere_points: int = SASA_SPHERE_POINTS,
) -> Any:
    """Return, for each residue of ``protein``, how much ``occluder`` covers its surface.

    On every frame, the SASA of each protein residue is computed twice as in
    :func:`residue_sasa`, one frame per MDTraj call: with the protein alone,
    and with the protein and ``occluder``, whose atoms cover the protein
    without being counted. When ``pbc`` is true and the frame has a box, each
    occluder molecule is first moved whole to its periodic image nearest the
    protein, since the SASA calculation ignores periodic images. A residue's
    relative SASA is its SASA over its maximum ASA from Tien et al. 2013 (the
    ``theoretical`` or ``empirical`` column, see
    :func:`~polyzymd.analyses.shared.aa_classification.get_max_asa`). It is
    exposed on a frame when its relative SASA with the protein alone is at
    least ``exposed_threshold``, and in contact when it is exposed, its
    relative SASA with the occluder is below ``buried_threshold``, and the
    occluder lowers its SASA. ``exposed_threshold=0`` counts every residue
    exposed, so any residue the occluder brings below ``buried_threshold``
    is in contact; ``buried_threshold`` must be above 0.

    The rows, named in :data:`OCCLUSION_PARTS`, are each residue's fraction
    of frames in contact; its fraction of frames exposed; its mean occluded
    area, the SASA lost to the occluder, ``max(0, alone - with)``, in Å²; and
    its mean SASA with the protein alone, in Å². Then one row per residue
    name in ``types``: the fraction of frames in contact when only the
    occluder atoms of that residue name are present. Columns follow the
    residues of ``protein`` that have a maximum ASA; residues without one,
    such as terminal caps, still cover their neighbours but are not measured.

    Returns
    -------
    numpy.ndarray
        Shape ``(4 + len(types), n_measured_residues)``.
    """
    import numpy as np

    n_measured = len(_measured_residues(protein, max_asa))
    sums = np.zeros((len(OCCLUSION_PARTS) + len(types), n_measured))
    for _, contact, exposed, alone, covered in _occlusion_frames(
        protein,
        occluder,
        frames,
        exposed_threshold,
        buried_threshold,
        types,
        max_asa,
        pbc,
        probe_radius_nm,
        n_sphere_points,
    ):
        sums[0] += contact[0]
        sums[1] += exposed
        sums[2] += np.maximum(0.0, alone - covered)
        sums[3] += alone
        sums[len(OCCLUSION_PARTS) :] += contact[1:]
    return sums / len(frames)


def _measured_residues(protein: Any, max_asa: str) -> Any:
    """Return the positions in ``protein.residues`` of the residues with a maximum ASA."""
    import numpy as np

    from polyzymd.analyses.shared.aa_classification import get_max_asa

    return np.array(
        [
            i
            for i, residue in enumerate(protein.residues)
            if get_max_asa(str(residue.resname), max_asa) is not None
        ],
        dtype=int,
    )


def _occlusion_frames(
    protein: Any,
    occluder: Any,
    frames: Any,
    exposed_threshold: float,
    buried_threshold: float,
    types: Any,
    max_asa: str,
    pbc: bool,
    probe_radius_nm: float,
    n_sphere_points: int,
) -> Any:
    """Yield, for each frame, what :func:`residue_occlusion` measures on it.

    Each item is the frame's time in ps; the contact rows, booleans of shape
    ``(1 + len(types), n_measured)``, with all ``occluder`` atoms and then
    with each residue name's; whether each measured residue is exposed; and
    its SASA in Å² alone and with all occluder atoms.
    """
    import numpy as np

    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.shared.aa_classification import get_max_asa

    if not (exposed_threshold >= 0 and buried_threshold > 0):
        raise ProtocolError(
            f"occlusion: exposed_threshold must be at least 0 and buried_threshold above 0, "
            f"got {exposed_threshold} and {buried_threshold}.",
            hint="Pass fractions of the maximum ASA, such as exposed_threshold=0.2 and "
            "buried_threshold=0.2.",
        )
    measured = _measured_residues(protein, max_asa)
    if len(measured) == 0:
        raise ProtocolError(
            "occlusion: no residue of the protein selection has a maximum ASA.",
            hint="Select standard amino acids, such as 'protein'.",
        )
    maximum = np.array([get_max_asa(str(protein.residues[i].resname), max_asa) for i in measured])
    exposed_limit, buried_limit = exposed_threshold * maximum, buried_threshold * maximum
    residue = np.unique(protein.resindices, return_inverse=True)[1]
    n_residues = len(protein.residues)
    groups = [occluder, *(occluder[occluder.resnames == name] for name in types)]
    contexts = [protein | group for group in groups]
    rows = [
        np.searchsorted(context.indices, group.indices) for context, group in zip(contexts, groups)
    ]
    picks = [np.searchsorted(occluder.indices, group.indices) for group in groups]

    def per_residue(context: Any, positions: Any) -> Any:
        atom = _atom_sasa(protein, context, positions[None], probe_radius_nm, n_sphere_points)
        return np.bincount(residue, weights=atom[0], minlength=n_residues)[measured]

    for ts in protein.universe.trajectory[frames]:
        alone = per_residue(protein, protein.positions)
        exposed = alone >= exposed_limit
        imaged = _nearest_images(protein, occluder, ts.dimensions if pbc else None)
        contact = np.zeros((len(groups), len(measured)), dtype=bool)
        covered_all = alone
        for k, (context, group) in enumerate(zip(contexts, groups)):
            if len(group) == 0:
                covered = alone
            else:
                positions = context.positions
                positions[rows[k]] = imaged[picks[k]]
                covered = per_residue(context, positions)
            contact[k] = exposed & (covered < buried_limit) & (covered < alone)
            if k == 0:
                covered_all = covered
        yield float(ts.time), contact, exposed, alone, covered_all


#: Rows of :func:`contact_lifetimes`, one column per group.
LIFETIME_PARTS = ("mean_lifetime", "n_events", "censored_fraction")


def contact_events(mask: Any, gap: int = 0) -> tuple[Any, Any]:
    """Return the lengths in frames of the runs of ``True`` in each column of ``mask``, and which are censored.

    ``mask`` has one row per frame and one column per series, such as one
    residue. An event is a run of consecutive frames in which a column is
    true. With ``gap``, absences of at most ``gap`` consecutive frames between
    two presences are filled first, with
    ``MDAnalysis.lib.correlations.correct_intermittency``, so they do not end
    the event. An event that includes the first or the last frame is
    censored: it may have begun before or lasted beyond the frames, so its
    length is only a lower bound.

    Returns
    -------
    tuple of numpy.ndarray
        The lengths in frames of every event, column by column, and a boolean
        array that is true for the censored ones.
    """
    import numpy as np

    mask = np.asarray(mask, dtype=bool)
    if gap > 0:
        from MDAnalysis.lib.correlations import correct_intermittency

        sets = [set(np.flatnonzero(row).tolist()) for row in mask]
        filled = correct_intermittency(sets, gap)
        mask = np.zeros_like(mask)
        for i, present in enumerate(filled):
            mask[i, list(present)] = True
    padded = np.zeros((mask.shape[0] + 2, mask.shape[1]), dtype=np.int8)
    padded[1:-1] = mask
    change = np.diff(padded, axis=0)
    lengths, censored = [], []
    for column in range(mask.shape[1]):
        starts = np.flatnonzero(change[:, column] == 1)
        ends = np.flatnonzero(change[:, column] == -1)
        lengths.append(ends - starts)
        censored.append((starts == 0) | (ends == mask.shape[0]))
    return np.concatenate(lengths).astype(int), np.concatenate(censored).astype(bool)


def restricted_mean_lifetime(durations: Any, censored: Any, horizon: float) -> float:
    """Return the Kaplan-Meier restricted mean of ``durations`` up to ``horizon``.

    ``scipy.stats.ecdf`` estimates the survival function, the probability that
    an event lasts longer than ``t``, from the durations, the ``censored``
    ones being only lower bounds (Kaplan and Meier 1958). The restricted mean
    is its area from 0 to ``horizon`` (Royston and Parmar 2013): the mean
    duration with every duration cut at ``horizon``. It is ``nan`` without
    durations.
    """
    import numpy as np
    from scipy.stats import CensoredData, ecdf

    durations = np.asarray(durations, dtype=float)
    censored = np.asarray(censored, dtype=bool)
    if len(durations) == 0:
        return float("nan")
    survival = ecdf(CensoredData(uncensored=durations[~censored], right=durations[censored])).sf
    times = np.concatenate([[0.0], survival.quantiles, [np.inf]])
    levels = np.concatenate([[1.0], survival.probabilities])
    ends = np.minimum(times[1:], horizon)
    starts = np.minimum(times[:-1], horizon)
    return float(np.sum(levels * (ends - starts)))


def contact_lifetimes(
    protein: Any,
    polymer: Any,
    frames: Any,
    method: str = "occlusion",
    types: Any = (),
    tolerance_ps: float = 0.0,
    **options: Any,
) -> Any:
    """Return how long contacts of ``protein`` residues with ``polymer`` last.

    Each frame's contacts come from :func:`residue_occlusion` for ``method``
    ``occlusion`` or :func:`residue_contacts` for ``distance``, with
    ``options`` passed on (``exposed_threshold``, ``buried_threshold``, ``max_asa``, ``pbc``,
    ``probe_radius_nm`` and ``n_sphere_points``, or ``cutoff`` and ``pbc``).
    An event is a run of consecutive frames in which one residue is in
    contact, found by :func:`contact_events`; absences of at most
    ``tolerance_ps`` ps are filled first. A run of ``k`` frames lasts ``k``
    times the frame spacing. Events of all residues are pooled.

    The columns are the polymer, then each residue name in ``types``, whose
    contacts are those of its atoms alone. The rows, named in
    :data:`LIFETIME_PARTS`, are the Kaplan-Meier restricted mean lifetime in
    ns, up to the time the frames span (see :func:`restricted_mean_lifetime`),
    ``nan`` without events; the number of events; and the fraction of them
    censored by the first or last frame.

    Returns
    -------
    numpy.ndarray
        Shape ``(3, 1 + len(types))``.
    """
    import numpy as np

    from polyzymd.analyses.exceptions import ProtocolError

    frames = list(frames)
    if tolerance_ps < 0:
        raise ProtocolError(
            f"contact_lifetimes: tolerance_ps must be at least 0, got {tolerance_ps}.",
            hint="Pass tolerance_ps=0 for events that end at the first absent frame.",
        )
    if method == "occlusion":
        known = {
            "exposed_threshold",
            "buried_threshold",
            "max_asa",
            "pbc",
            "probe_radius_nm",
            "n_sphere_points",
        }
        settings = {
            "exposed_threshold": EXPOSED_THRESHOLD,
            "buried_threshold": BURIED_THRESHOLD,
            "max_asa": "theoretical",
            "pbc": True,
            "probe_radius_nm": SASA_PROBE_RADIUS_NM,
            "n_sphere_points": SASA_SPHERE_POINTS,
        }
    elif method == "distance":
        known = {"cutoff", "pbc"}
        settings = {"cutoff": CONTACT_CUTOFF, "pbc": True}
    else:
        raise ProtocolError(
            f"contact_lifetimes: method must be 'occlusion' or 'distance', got {method!r}.",
            hint="Pass method='occlusion' or method='distance'.",
        )
    if set(options) - known:
        raise ProtocolError(
            f"contact_lifetimes: method {method} takes no option "
            f"{', '.join(sorted(set(options) - known))}.",
            hint=f"Pass only {', '.join(sorted(known))}.",
        )
    settings.update(options)
    if method == "occlusion":
        stream = (
            (time, contact)
            for time, contact, *_ in _occlusion_frames(
                protein,
                polymer,
                frames,
                settings["exposed_threshold"],
                settings["buried_threshold"],
                types,
                settings["max_asa"],
                settings["pbc"],
                settings["probe_radius_nm"],
                settings["n_sphere_points"],
            )
        )
    else:
        stream = _distance_frames(
            protein, polymer, frames, settings["cutoff"], types, settings["pbc"]
        )
    times, masks = [], []
    for time, contact in stream:
        times.append(time)
        masks.append(contact)
    if len(times) < 2:
        raise ProtocolError(
            f"contact_lifetimes: {len(times)} frame, and a lifetime needs at least two.",
            hint="Use more production frames, or a smaller stride.",
        )
    spacing = np.diff(times)
    step = float(np.median(spacing))
    if step <= 0 or np.max(np.abs(spacing - step)) > 1e-3 * step:
        raise ProtocolError(
            f"contact_lifetimes: the frames are not evenly spaced in time (from "
            f"{spacing.min():.6g} to {spacing.max():.6g} ps apart), so a run of frames has no "
            "single duration.",
            hint="Check the trajectory's time axis for missing or repeated frames.",
        )
    # Frame times carry float rounding, so 40 ps at a 40.0000001 ps spacing is one frame.
    gap = int(np.floor(tolerance_ps / step * (1 + 1e-6))) if tolerance_ps > 0 else 0
    masks = np.asarray(masks)  # frames, groups, residues
    horizon = len(times) * step / 1000.0
    result = np.zeros((len(LIFETIME_PARTS), masks.shape[1]))
    for group in range(masks.shape[1]):
        lengths, censored = contact_events(masks[:, group, :], gap)
        durations = lengths * step / 1000.0
        result[0, group] = restricted_mean_lifetime(durations, censored, horizon)
        result[1, group] = len(lengths)
        result[2, group] = float(np.mean(censored)) if len(lengths) else float("nan")
    return result
