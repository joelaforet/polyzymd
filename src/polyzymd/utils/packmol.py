"""
Custom Packmol input generation and execution utilities.

This module provides a thin replacement for the OpenFF Interchange
``pack_box()`` polymer-packing path, adding support for Packmol keywords
that the OpenFF wrapper does not expose (currently: ``movebadrandom``).

All heavy imports (openff.interchange, openff.units, numpy) are lazy so
this module can be imported in environments without the full simulation
stack.

Typical usage
-------------
>>> from polyzymd.utils.packmol import build_packmol_input, run_packmol
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from numpy.typing import NDArray

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_PACKMOL_INPUT_FILE = "packmol_input.txt"
_PACKMOL_OUTPUT_FILE = "packmol_output.pdb"
_PACKMOL_SOLUTE_FILE = "_PACKING_SOLUTE.pdb"
_PACKMOL_MOLECULE_PREFIX = "_PACKING_MOLECULE"


# ---------------------------------------------------------------------------
# Post-assembly geometry assertion
# ---------------------------------------------------------------------------

#: Number of packed atoms closer than ``0.5 * tolerance`` to the solute above which
#: the assembled system is treated as a frame mismatch rather than imperfect packing.
#:
#: Packmol exit code 173 ("imperfect packing") legitimately leaves a handful of
#: residual contacts (observed: 1 atom at 0.6 A and 3-5 atoms below the tolerance in
#: dense 30-chain polymer shells); energy minimisation resolves those.  A solute/solvent
#: frame mismatch (d96b1fcd) instead puts hundreds to thousands of solvent atoms inside
#: the solute (observed: 1182 atoms below 1 A in a defective CALB control build).  The
#: limit sits two orders of magnitude from both observations.
SOLVATION_CLASH_ATOM_LIMIT = 20


class SolvationClashError(ValueError):
    """Raised when packed solvent/polymer atoms overlap the fixed solute.

    Packmol guarantees that every placed atom is at least ``tolerance`` away
    from the fixed solute *in the Packmol frame*.  If the assembled topology
    combines solute and solvent coordinates expressed in different frames
    (the solute/solvent frame offset fixed in d96b1fcd), hundreds of solvent
    molecules end up inside the solute while Packmol reports success.  This
    error turns that silent defect into a hard build failure.
    """


class PeriodicImageClashError(ValueError):
    """Raised when an atom overlaps one of its own periodic images.

    PolyzyMD packs into the rectangular *brick* that represents a triclinic
    (rhombic-dodecahedron) cell.  If molecules are packed in a region larger
    than that brick — as they were when polymers were packed in a separate
    rectangular box before the final cell was known — atoms protrude through
    the brick faces and land on top of their images across a lattice vector.
    Packmol never sees those contacts (it is run without periodicity), so the
    defect is silent until minimisation blows up or the run dies with NaN.
    This error turns it into a hard build failure.
    """


def periodic_image_statistics(
    positions_angstrom: "NDArray",
    box_vectors_angstrom: "NDArray",
    *,
    tolerance_angstrom: float,
) -> dict[str, object]:
    """Closest-periodic-image statistics for a set of coordinates.

    Every atom is translated by each of the 26 non-zero lattice vectors
    ``i*a + j*b + k*c`` (``i, j, k`` in ``{-1, 0, 1}``) and queried against a
    KD-tree of the untranslated coordinates.  The zero translation is skipped,
    so an atom is never counted against itself; an atom *is* counted against
    its own image across a lattice vector, which is exactly the defect being
    looked for.

    Parameters
    ----------
    positions_angstrom : NDArray
        Coordinates in Angstrom, shape ``(N, 3)``.
    box_vectors_angstrom : NDArray
        Row-major box vectors in Angstrom, shape ``(3, 3)``.
    tolerance_angstrom : float
        Packmol tolerance for the run; distances are only resolved below it
        (``distance_upper_bound``), which keeps the 26 queries fast for
        100k-atom systems.

    Returns
    -------
    dict
        ``n_atoms``, ``n_atoms_below_tolerance``,
        ``n_atoms_below_half_tolerance`` (distinct atoms involved in at least
        one such contact, both partners counted), ``min_distance_angstrom``
        (``inf`` when no image lies within the tolerance),
        ``worst_pair`` (``(shifted_atom, image_partner)`` indices or ``None``)
        and ``worst_lattice_vector`` (``(i, j, k)`` or ``None``).
    """
    import itertools

    import numpy as np
    from scipy.spatial import cKDTree

    positions = np.asarray(positions_angstrom, dtype=float).reshape(-1, 3)
    box = np.asarray(box_vectors_angstrom, dtype=float).reshape(3, 3)

    n_atoms = int(positions.shape[0])
    stats: dict[str, object] = {
        "n_atoms": n_atoms,
        "n_atoms_below_tolerance": 0,
        "n_atoms_below_half_tolerance": 0,
        "min_distance_angstrom": float("inf"),
        "worst_pair": None,
        "worst_lattice_vector": None,
    }
    if n_atoms == 0:
        return stats

    tree = cKDTree(positions)
    below_tolerance = np.zeros(n_atoms, dtype=bool)
    below_half = np.zeros(n_atoms, dtype=bool)
    half_tolerance = 0.5 * tolerance_angstrom

    for shift in itertools.product((-1, 0, 1), repeat=3):
        if shift == (0, 0, 0):
            continue
        translation = np.asarray(shift, dtype=float) @ box
        distances, neighbours = tree.query(
            positions + translation,
            k=1,
            distance_upper_bound=tolerance_angstrom,
        )
        hits = distances < tolerance_angstrom
        if not hits.any():
            continue
        below_tolerance[hits] = True
        below_tolerance[neighbours[hits]] = True
        close = distances < half_tolerance
        if close.any():
            below_half[close] = True
            below_half[neighbours[close]] = True
        nearest = int(np.argmin(distances))
        if float(distances[nearest]) < float(stats["min_distance_angstrom"]):
            stats["min_distance_angstrom"] = float(distances[nearest])
            stats["worst_pair"] = (nearest, int(neighbours[nearest]))
            stats["worst_lattice_vector"] = tuple(int(v) for v in shift)

    stats["n_atoms_below_tolerance"] = int(np.count_nonzero(below_tolerance))
    stats["n_atoms_below_half_tolerance"] = int(np.count_nonzero(below_half))
    return stats


def _box_vectors_as_angstrom(box_vectors) -> "NDArray":
    """Return row-major box vectors as a plain ``(3, 3)`` Angstrom array."""
    import numpy as np

    if hasattr(box_vectors, "m_as"):
        box_vectors = box_vectors.m_as("angstrom")
    return np.asarray(box_vectors, dtype=float).reshape(3, 3)


def _assert_periodic_image_separation(
    topology,
    box_vectors,
    *,
    tolerance_angstrom: float,
    label: str = "system",
) -> dict[str, object]:
    """Fail loudly if any atom overlaps one of its own periodic images.

    Called after polymer packing (solute + polymers) and after solvation (the
    whole system).  Atoms closer than ``0.5 * tolerance_angstrom`` to an image
    atom raise :class:`PeriodicImageClashError`; contacts between
    ``0.5 * tolerance`` and ``tolerance`` only warn, since Packmol's own
    tolerance is not enforced across the periodic boundary and minimisation
    resolves a marginal contact.

    Parameters
    ----------
    topology : openff.toolkit.Topology
        Assembled topology carrying positions.
    box_vectors : openff.units.Quantity or NDArray
        Row-major periodic box vectors of the assembled system.
    tolerance_angstrom : float
        Packmol tolerance used for the run.
    label : str
        Human-readable name of the stage for log and error messages.

    Returns
    -------
    dict
        The statistics from :func:`periodic_image_statistics`.

    Raises
    ------
    PeriodicImageClashError
        If any atom lies within ``0.5 * tolerance_angstrom`` of an image atom.
    """
    import numpy as np

    positions = np.asarray(topology.get_positions().m_as("angstrom"), dtype=float)
    box = _box_vectors_as_angstrom(box_vectors)
    stats = periodic_image_statistics(positions, box, tolerance_angstrom=tolerance_angstrom)

    if stats["n_atoms_below_half_tolerance"] > 0:
        pair = stats["worst_pair"]
        raise PeriodicImageClashError(
            f"{stats['n_atoms_below_half_tolerance']} atom(s) of the {label} lie within "
            f"{0.5 * tolerance_angstrom:.2f} A of a periodic image "
            f"({stats['n_atoms_below_tolerance']} within the {tolerance_angstrom:.2f} A "
            f"Packmol tolerance; minimum image separation "
            f"{stats['min_distance_angstrom']:.3f} A between atoms {pair} across lattice "
            f"vector {stats['worst_lattice_vector']}; {stats['n_atoms']} atoms checked). "
            "Molecules were packed outside the periodic brick, so they overlap themselves "
            "across the cell boundary; minimisation cannot resolve this and the run would "
            "die with NaN. Refusing to continue the build."
        )

    if stats["n_atoms_below_tolerance"] > 0:
        logger.warning(
            "%d atom(s) of the %s lie between %.2f and %.2f A of a periodic image "
            "(minimum %.3f A); Packmol does not enforce its tolerance across the "
            "periodic boundary and minimisation will resolve this.",
            stats["n_atoms_below_tolerance"],
            label,
            0.5 * tolerance_angstrom,
            tolerance_angstrom,
            stats["min_distance_angstrom"],
        )
    else:
        logger.info(
            "Periodic-image separation check passed for the %s: %d atoms, no image "
            "closer than %.2f A.",
            label,
            stats["n_atoms"],
            tolerance_angstrom,
        )
    return stats


def separation_statistics(
    solute_xyz: "NDArray",
    other_xyz: "NDArray",
    *,
    tolerance_angstrom: float,
) -> dict[str, float | int]:
    """Nearest-solute-atom distance statistics for a set of packed atoms.

    Parameters
    ----------
    solute_xyz, other_xyz : NDArray
        Coordinates in Angstrom, shape ``(N, 3)`` and ``(M, 3)``.  Direct
        (non-periodic) distances are used: a frame mismatch shows up as
        direct overlap, and direct distances never under-report a clash.
    tolerance_angstrom : float
        Packmol tolerance that was requested for the packing run.

    Returns
    -------
    dict
        ``n_other`` (M), ``n_below_tolerance``, ``n_below_half_tolerance``,
        ``min_distance_angstrom`` (``inf`` when either set is empty).
    """
    import numpy as np
    from scipy.spatial import cKDTree

    solute_xyz = np.asarray(solute_xyz, dtype=float).reshape(-1, 3)
    other_xyz = np.asarray(other_xyz, dtype=float).reshape(-1, 3)
    stats: dict[str, float | int] = {
        "n_other": int(other_xyz.shape[0]),
        "n_below_tolerance": 0,
        "n_below_half_tolerance": 0,
        "min_distance_angstrom": float("inf"),
    }
    if solute_xyz.shape[0] == 0 or other_xyz.shape[0] == 0:
        return stats

    distances, _ = cKDTree(solute_xyz).query(other_xyz, k=1)
    stats["n_below_tolerance"] = int(np.count_nonzero(distances < tolerance_angstrom))
    stats["n_below_half_tolerance"] = int(np.count_nonzero(distances < 0.5 * tolerance_angstrom))
    stats["min_distance_angstrom"] = float(distances.min())
    return stats


def _assert_solute_solvent_separation(
    topology,
    n_solute_atoms: int,
    *,
    tolerance_angstrom: float,
    label: str = "solvent",
) -> dict[str, float | int]:
    """Fail loudly if packed atoms sit inside the solute after assembly.

    The first *n_solute_atoms* atoms of *topology* are the fixed solute; the
    remainder are the freshly packed molecules.  Packed atoms closer than
    ``0.5 * tolerance_angstrom`` to a solute atom count as clashes.  More than
    :data:`SOLVATION_CLASH_ATOM_LIMIT` clashing atoms means the solute and the
    packed coordinates are expressed in different frames and the build is
    aborted; a handful of clashes is the residue of an imperfect Packmol run
    (exit code 173) that minimisation resolves, and only produces a warning,
    as do distances between ``0.5 * tolerance`` and ``tolerance``.

    Parameters
    ----------
    topology : openff.toolkit.Topology
        Assembled topology (solute first, packed molecules after).
    n_solute_atoms : int
        Number of leading solute atoms.  ``0`` disables the check.
    tolerance_angstrom : float
        Packmol tolerance used for the run.
    label : str
        Human-readable name for the packed species in messages.

    Returns
    -------
    dict
        The statistics from :func:`separation_statistics`.

    Raises
    ------
    SolvationClashError
        If more than :data:`SOLVATION_CLASH_ATOM_LIMIT` packed atoms are closer
        than ``0.5 * tolerance_angstrom`` to the solute.
    """
    import numpy as np

    if n_solute_atoms <= 0:
        return {
            "n_other": 0,
            "n_below_tolerance": 0,
            "n_below_half_tolerance": 0,
            "min_distance_angstrom": float("inf"),
        }

    positions = np.asarray(topology.get_positions().m_as("angstrom"), dtype=float)
    stats = separation_statistics(
        positions[:n_solute_atoms],
        positions[n_solute_atoms:],
        tolerance_angstrom=tolerance_angstrom,
    )

    if stats["n_below_half_tolerance"] > SOLVATION_CLASH_ATOM_LIMIT:
        raise SolvationClashError(
            f"{stats['n_below_half_tolerance']} {label} atom(s) lie within "
            f"{0.5 * tolerance_angstrom:.2f} A of the solute "
            f"({stats['n_below_tolerance']} within the {tolerance_angstrom:.2f} A Packmol "
            f"tolerance; minimum separation {stats['min_distance_angstrom']:.3f} A; "
            f"{stats['n_other']} {label} atoms checked; limit {SOLVATION_CLASH_ATOM_LIMIT}). "
            "Imperfect Packmol runs leave at most a handful of such contacts, so the "
            "assembled solute and packed coordinates are almost certainly expressed in "
            "different frames (solute/solvent frame mismatch, see d96b1fcd). Refusing to "
            "continue the build."
        )

    if stats["n_below_half_tolerance"] > 0:
        logger.warning(
            "%d %s atom(s) lie within %.2f A of the solute (minimum %.3f A); this is the "
            "residue of an imperfect Packmol run and will be resolved by minimisation.",
            stats["n_below_half_tolerance"],
            label,
            0.5 * tolerance_angstrom,
            stats["min_distance_angstrom"],
        )
    elif stats["n_below_tolerance"] > 0:
        logger.warning(
            "%d %s atom(s) lie between %.2f and %.2f A of the solute "
            "(minimum %.3f A); Packmol tolerance was not fully honoured.",
            stats["n_below_tolerance"],
            label,
            0.5 * tolerance_angstrom,
            tolerance_angstrom,
            stats["min_distance_angstrom"],
        )
    else:
        logger.info(
            "Solute/%s separation check passed: %d atoms, minimum %.3f A (tolerance %.2f A)",
            label,
            stats["n_other"],
            stats["min_distance_angstrom"],
            tolerance_angstrom,
        )
    return stats


# ---------------------------------------------------------------------------
# Input-file builder
# ---------------------------------------------------------------------------


def build_packmol_input(
    molecule_pdb_paths: list[str],
    molecule_counts: list[int],
    box_size_angstrom: "NDArray",
    tolerance_angstrom: float,
    solute_pdb_path: str | None = None,
    use_pbc: bool = False,
    movebadrandom: bool = False,
    ignore_conect: bool = False,
    inner_exclusion_box_angstrom: "NDArray | None" = None,
    inside_sphere_angstrom: "NDArray | None" = None,
    nloop: int | None = None,
    seed: int | None = None,
) -> str:
    """Build the text content of a Packmol input file.

    Parameters
    ----------
    molecule_pdb_paths : list[str]
        Paths to PDB files for each unique molecule type to pack.
    molecule_counts : list[int]
        Number of copies of each molecule type. Entries of zero are skipped.
    box_size_angstrom : NDArray
        1-D array of length 3 giving the box edge lengths in Angstrom.
        For non-PBC runs the effective packing box is shrunk by *tolerance*
        (matching the OpenFF convention).
    tolerance_angstrom : float
        Minimum distance between atoms of different molecules, in Angstrom.
    solute_pdb_path : str or None, optional
        Path to a PDB file for the fixed solute (protein, substrate, …).
        When provided a ``fixed`` section is emitted first.
    use_pbc : bool, optional
        If ``True``, emit a ``pbc`` global keyword instead of per-structure
        ``inside box`` constraints (requires Packmol ≥ 20.15.0).
    movebadrandom : bool, optional
        If ``True``, add the ``movebadrandom`` global keyword, which places
        badly-packed molecules at random positions in the box rather than
        near well-packed neighbours. This improves convergence when the
        restraints are complex (many unique chain types, dense packing).
        Default is ``False``.
    ignore_conect : bool, optional
        If ``True``, add the ``ignore_conect`` Packmol keyword so that
        CONECT records in PDB input files are not parsed.  This prevents
        failures when atom indices exceed the 5-digit PDB fixed-width
        field limit (>99,999 atoms).  Default is ``False``.
    inner_exclusion_box_angstrom : NDArray or None, optional
        When provided, a 1-D array of 6 floats
        ``[xmin, ymin, zmin, xmax, ymax, zmax]`` in Angstrom defining a
        rectangular region that packed molecules must avoid.  An
        ``outside box`` constraint is added to every non-fixed structure
        block, creating a rectangular *shell* between the inner exclusion
        zone and the outer packing box.  Ignored when *use_pbc* is
        ``True``.  Default is ``None`` (no exclusion zone).
    inside_sphere_angstrom : NDArray or None, optional
        When provided, a 1-D array of 4 floats ``[cx, cy, cz, radius]`` in
        Angstrom.  An ``inside sphere`` constraint is added to every
        non-fixed structure block, so packed molecules stay within that
        sphere *in addition to* the ``inside box`` packing box.  Used to
        keep polymer chains in a shell around the solute when the packing
        box is the full periodic brick.  Ignored when *use_pbc* is ``True``.
        Default is ``None`` (no spherical confinement).
    nloop : int or None, optional
        Maximum number of GENCAN optimisation loops *per molecule type*.
        Packmol's default is 50; for dense shell-packing with many
        molecules, higher values (200–500) improve convergence.
        Default is ``None`` (use Packmol's built-in default).
    seed : int or None, optional
        Seed for Packmol's random number generator (``seed`` keyword).
        ``None`` omits the keyword, so Packmol uses its fixed built-in
        default and every run with identical inputs produces identical
        coordinates.  Pass the replicate index to obtain independent
        starting configurations per replicate.

    Returns
    -------
    str
        Complete Packmol input file text, ready to be written to disk.
    """
    import numpy as np

    if use_pbc:
        effective_box = np.asarray(box_size_angstrom, dtype=float)
    else:
        effective_box = np.asarray(box_size_angstrom, dtype=float) - tolerance_angstrom

    # Pre-format the exclusion constraint line (if requested and not PBC)
    _exclusion_line: str | None = None
    if inner_exclusion_box_angstrom is not None and not use_pbc:
        ebox = np.asarray(inner_exclusion_box_angstrom, dtype=float)
        if ebox.shape != (6,):
            raise ValueError(f"inner_exclusion_box_angstrom must have shape (6,), got {ebox.shape}")
        _exclusion_line = (
            f"  outside box"
            f" {ebox[0]:.6f} {ebox[1]:.6f} {ebox[2]:.6f}"
            f" {ebox[3]:.6f} {ebox[4]:.6f} {ebox[5]:.6f}"
        )

    _sphere_line: str | None = None
    if inside_sphere_angstrom is not None and not use_pbc:
        sphere = np.asarray(inside_sphere_angstrom, dtype=float)
        if sphere.shape != (4,):
            raise ValueError(f"inside_sphere_angstrom must have shape (4,), got {sphere.shape}")
        _sphere_line = (
            f"  inside sphere" f" {sphere[0]:.6f} {sphere[1]:.6f} {sphere[2]:.6f} {sphere[3]:.6f}"
        )

    lines: list[str] = [
        f"tolerance {tolerance_angstrom:f}",
        "filetype pdb",
        f"output {_PACKMOL_OUTPUT_FILE}",
        "",
    ]

    if ignore_conect:
        lines.append("ignore_conect")
        lines.append("")

    if movebadrandom:
        lines.append("movebadrandom")
        lines.append("")

    if nloop is not None:
        lines.append(f"nloop {nloop}")
        lines.append("")

    if seed is not None:
        lines.append(f"seed {int(seed)}")
        lines.append("")

    if use_pbc:
        lines.append(
            f"pbc 0. 0. 0. {effective_box[0]:.6f} {effective_box[1]:.6f} {effective_box[2]:.6f}"
        )
        lines.append("")

    # Fixed solute block
    if solute_pdb_path is not None:
        lines.extend(
            [
                f"structure {solute_pdb_path}",
                "  number 1",
                "  fixed 0. 0. 0. 0. 0. 0.",
                "end structure",
                "",
            ]
        )

    # One block per unique molecule type
    for pdb_path, count in zip(molecule_pdb_paths, molecule_counts):
        if count == 0:
            continue
        block = [
            f"structure {pdb_path}",
            f"  number {count}",
        ]
        if not use_pbc:
            block.append(
                f"  inside box 0. 0. 0."
                f" {effective_box[0]:.6f} {effective_box[1]:.6f} {effective_box[2]:.6f}"
            )
        if _sphere_line is not None:
            block.append(_sphere_line)
        if _exclusion_line is not None:
            block.append(_exclusion_line)
        block.append("end structure")
        block.append("")
        lines.extend(block)

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Packmol executor
# ---------------------------------------------------------------------------


def run_packmol(
    input_text: str,
    working_directory: str | Path,
    retain_working_files: bool = True,
) -> Path:
    """Write a Packmol input file and execute Packmol.

    Parameters
    ----------
    input_text : str
        Complete Packmol input file content (from :func:`build_packmol_input`).
    working_directory : str or Path
        Directory in which to write working files and invoke Packmol.
        The directory is created if it does not exist.
    retain_working_files : bool, optional
        When ``True`` (default), all files in *working_directory* are kept
        after the run.  When ``False`` the directory is removed on success
        (mimicking OpenFF behaviour for temporary directories).

    Returns
    -------
    Path
        Absolute path to the Packmol output PDB file.

    Raises
    ------
    OSError
        If the ``packmol`` binary cannot be found on ``PATH``.
    RuntimeError
        If Packmol exits with a non-zero return code or does not print
        ``'Success!'`` in its output.
    """
    packmol_binary = shutil.which("packmol")
    if packmol_binary is None:
        raise OSError(
            "Packmol binary not found on PATH. "
            "Install Packmol and make sure it is accessible as 'packmol'."
        )

    working_directory = Path(working_directory)
    working_directory.mkdir(parents=True, exist_ok=True)

    _temporary = False
    _actual_dir = working_directory

    input_path = working_directory / _PACKMOL_INPUT_FILE
    output_path = working_directory / _PACKMOL_OUTPUT_FILE
    error_log_path = working_directory / "packmol_error.log"

    input_path.write_text(input_text)

    logger.debug("Running Packmol in %s", working_directory)
    logger.debug("Input file:\n%s", input_text)

    original_cwd = Path.cwd()
    try:
        os.chdir(working_directory)
        with input_path.open() as fh:
            result = subprocess.run(
                packmol_binary,
                stdin=fh,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
    finally:
        os.chdir(original_cwd)

    stdout = result.stdout.decode("utf-8", errors="replace")

    # Exit code 173 means "ended without perfect packing".  Packmol still
    # writes its best solution to the output file.  Since systems are
    # energy-minimised before MD, minor steric violations are acceptable.
    # Treat 173 as a warning rather than a fatal error.
    _PACKMOL_EXIT_IMPERFECT = 173

    if result.returncode == _PACKMOL_EXIT_IMPERFECT:
        error_log_path.write_text(stdout)
        logger.warning(
            "Packmol exited with code %d (imperfect packing). "
            "The best solution was written to %s and will be used. "
            "Minor steric clashes will be resolved during energy minimisation. "
            "See %s for details.",
            _PACKMOL_EXIT_IMPERFECT,
            output_path.name,
            error_log_path,
        )
    elif result.returncode != 0:
        error_log_path.write_text(stdout)
        raise RuntimeError(
            f"Packmol exited with return code {result.returncode}. "
            f"See {error_log_path} for details."
        )
    elif "Success!" not in stdout:
        raise RuntimeError(
            "Packmol did not raise an error code but 'Success!' was not found "
            "in its output. The packing may not have converged. "
            f"Working directory: {working_directory}"
        )

    if not retain_working_files and _temporary:
        shutil.rmtree(_actual_dir, ignore_errors=True)

    return output_path.resolve()


# ---------------------------------------------------------------------------
# High-level polymer packing helper
# ---------------------------------------------------------------------------


def pack_polymers(
    molecules: list,
    number_of_copies: list[int],
    solute,
    box_vectors,
    *,
    tolerance_angstrom: float = 2.0,
    movebadrandom: bool = False,
    nloop: int | None = 200,
    seed: int | None = None,
    exclude_solute_bbox: bool = False,
    confine_to_sphere: bool = True,
    sphere_padding_angstrom: float = 20.0,
    working_directory: str | Path | None = None,
    retain_working_files: bool = True,
):
    """Pack polymer chains around a fixed solute using Packmol.

    This is a drop-in replacement for the OpenFF ``pack_box()`` call in
    :meth:`~polyzymd.builders.system_builder.SystemBuilder.pack_polymers`.
    It adds support for the ``movebadrandom`` Packmol keyword.  *box_vectors*
    are the **final periodic box of the simulation**: the solute is centred in
    the rectangular brick of that cell and the chains are packed inside the
    brick (shrunk by the tolerance), so no atom can protrude through a brick
    face and overlap its own periodic image.  By default the chains are also
    confined to a sphere centred on the solute (``confine_to_sphere``) so that
    they stay in a shell around the protein instead of spreading into the
    corners of the brick; the fixed solute and the Packmol tolerance keep them
    off the protein itself.  Setting ``exclude_solute_bbox`` additionally
    restores the older behaviour of an ``outside box`` annulus around the
    solute's bounding box, which over-constrains long chains into a thin
    annulus and makes Packmol converge slowly or not at all.

    Parameters
    ----------
    molecules : list[openff.toolkit.Molecule]
        Unique molecule types to pack.
    number_of_copies : list[int]
        Number of copies of each molecule type (parallel to *molecules*).
    solute : openff.toolkit.Topology
        Fixed topology (protein + substrate) placed at the origin.
    box_vectors : openff.units.Quantity
        **Final** periodic box vectors of the simulation cell, shape (3, 3)
        with length units (e.g. nanometers).  Packing happens inside the
        rectangular brick of this cell.
    tolerance_angstrom : float, optional
        Packmol tolerance in Angstrom (default 2.0).
    movebadrandom : bool, optional
        Pass the ``movebadrandom`` keyword to Packmol (default ``False``).
    nloop : int or None, optional
        Maximum GENCAN optimisation loops per molecule type.  Packmol's
        default is 50; for dense shell-packing with many molecules, higher
        values (200-500) improve convergence.  Default is ``200``.
    seed : int or None, optional
        Packmol random seed (see :func:`build_packmol_input`).
    exclude_solute_bbox : bool, optional
        Add an ``outside box`` constraint equal to the solute's bounding box
        inflated by the tolerance, forcing the chains into a rectangular shell.
        Default ``False`` (chains pack anywhere; the tolerance against the
        fixed solute prevents overlap).
    confine_to_sphere : bool, optional
        Add an ``inside sphere`` constraint centred on the solute with radius
        ``0.5 * |solute bounding-box diagonal| + sphere_padding_angstrom``,
        keeping the chains in a shell around the protein rather than in the
        corners of the brick.  Default ``True``.
    sphere_padding_angstrom : float, optional
        Padding added to the solute's bounding-box circumradius to obtain the
        confinement sphere, in Angstrom (default 20.0 = 2.0 nm, matching the
        ``polymers.packing.padding`` default).
    working_directory : str, Path, or None, optional
        Directory for Packmol input/output files.  A temporary directory is
        created when ``None``.
    retain_working_files : bool, optional
        Keep working files after the run (default ``True``).

    Returns
    -------
    openff.toolkit.Topology
        Packed topology with solute + all polymer chains and box vectors set.

    Raises
    ------
    SolvationClashError
        If any polymer atom ends up closer than ``0.5 * tolerance_angstrom``
        to the solute after assembly (solute/polymer frame mismatch).
    PeriodicImageClashError
        If any atom of the packed system ends up closer than
        ``0.5 * tolerance_angstrom`` to one of its own periodic images.
    """
    import numpy as np
    from openff.packmol._packmol import (
        _center_topology_at,
        _compute_brick_from_box_vectors,
        _create_molecule_pdbs,
        _create_solute_pdb,
        _load_positions,
    )
    from openff.toolkit import Topology
    from openff.units import Quantity

    # --- sort molecule types by count descending ---
    # Higher-count molecules first improves Packmol convergence (fewer
    # restarts when the most-replicated species is placed first).
    paired = list(zip(molecules, number_of_copies))
    paired.sort(key=lambda pair: pair[1], reverse=True)
    molecules = [p[0] for p in paired]
    number_of_copies = [p[1] for p in paired]

    logger.info(
        "Packmol molecule order (count-descending): %s",
        ", ".join(f"{n}x" for n in number_of_copies),
    )

    # --- resolve working directory ---
    _temporary = False
    if working_directory is None:
        working_directory = Path(tempfile.mkdtemp())
        _temporary = True
    else:
        working_directory = Path(working_directory).resolve()
        working_directory.mkdir(parents=True, exist_ok=True)

    # --- compute brick dimensions ---
    brick_size = _compute_brick_from_box_vectors(box_vectors)
    box_size_angstrom = np.asarray(brick_size.m_as("angstrom"), dtype=float)

    # --- center solute in the brick ---
    centered_solute = _center_topology_at("BRICK", solute, box_vectors, brick_size)

    # --- optional inner exclusion box (polymer shell constraint) ---
    # By default the chains may occupy the whole packing box; Packmol keeps
    # every polymer atom at least ``tolerance`` from the fixed solute, which
    # is all the geometry we need.  The legacy ``outside box`` shell equal to
    # the solute's bounding box (inflated by the tolerance) over-constrains
    # long chains into an annulus that is often thinner than the chains
    # themselves, and Packmol then grinds to its loop limit without
    # converging.
    inner_exclusion_box = None
    if exclude_solute_bbox:
        inner_exclusion_box = _solute_bbox_exclusion(
            centered_solute, box_size_angstrom, tolerance_angstrom, molecules
        )
    else:
        logger.info(
            "Polymers pack throughout the box (no solute bounding-box exclusion); "
            "the %.1f A Packmol tolerance against the fixed solute prevents overlap.",
            tolerance_angstrom,
        )

    # --- optional spherical confinement around the solute ---
    # The packing box is the *final* periodic brick, which is much larger than
    # the solute.  Without a further constraint Packmol happily fills the brick
    # corners, so chains end up far from the protein they are meant to shield.
    # The sphere keeps them in a shell around the solute while every atom still
    # lies inside the brick, which is what makes the periodic images safe.
    inside_sphere = None
    if confine_to_sphere and solute is not None:
        inside_sphere = solute_sphere_constraint(
            centered_solute, padding_angstrom=sphere_padding_angstrom
        )
        logger.info(
            "Confining polymers to a sphere at (%.2f, %.2f, %.2f) A with radius %.2f A "
            "inside a %.2f x %.2f x %.2f A brick.",
            *inside_sphere,
            *box_size_angstrom,
        )

    # Force PBC off for polymer packing so per-structure ``inside box`` (and
    # the optional ``outside box``) constraints apply.  PBC mode removes
    # per-structure spatial constraints.  (Solvation can still use PBC
    # since it doesn't need shell constraints.)
    _use_pbc = False

    original_cwd = Path.cwd()
    try:
        os.chdir(working_directory)

        # write PDB files
        solute_pdb = _create_solute_pdb(centered_solute, box_vectors)
        molecule_pdbs = _create_molecule_pdbs(molecules)

        # build and run packmol
        input_text = build_packmol_input(
            molecule_pdb_paths=molecule_pdbs,
            molecule_counts=number_of_copies,
            box_size_angstrom=box_size_angstrom,
            tolerance_angstrom=tolerance_angstrom,
            solute_pdb_path=solute_pdb,
            use_pbc=_use_pbc,
            movebadrandom=movebadrandom,
            inner_exclusion_box_angstrom=inner_exclusion_box,
            inside_sphere_angstrom=inside_sphere,
            nloop=nloop,
            seed=seed,
        )

        output_path = run_packmol(
            input_text=input_text,
            working_directory=working_directory,
            retain_working_files=True,  # always keep; we clean up below
        )

        positions = _load_positions(str(output_path.name))

    finally:
        os.chdir(original_cwd)

    # --- assemble output topology ---
    # Use the *centered* solute so that protein coordinates match the
    # Packmol frame of reference (polymers were packed around the centered
    # copy).  Previous code incorrectly used the original un-centered
    # solute, creating a spatial mismatch between protein and polymers.
    added_molecules = []
    for mol, n in zip(molecules, number_of_copies):
        added_molecules.extend([mol] * n)
    packed_topology = Topology.from_molecules(added_molecules)

    n_solute_atoms = len(positions) - packed_topology.n_atoms
    packed_topology.set_positions(Quantity(positions[n_solute_atoms:], "angstrom"))

    if solute is not None:
        packed_topology = centered_solute + packed_topology

    packed_topology.box_vectors = box_vectors

    if solute is not None:
        _assert_solute_solvent_separation(
            packed_topology,
            centered_solute.n_atoms,
            tolerance_angstrom=tolerance_angstrom,
            label="polymer",
        )

    _assert_periodic_image_separation(
        packed_topology,
        box_vectors,
        tolerance_angstrom=tolerance_angstrom,
        label="packed solute + polymers",
    )

    if _temporary and not retain_working_files:
        shutil.rmtree(working_directory, ignore_errors=True)

    return packed_topology


# ---------------------------------------------------------------------------
# High-level solvation helper
# ---------------------------------------------------------------------------


def solvate_with_packmol(
    molecules: list,
    number_of_copies: list[int],
    solute,
    box_vectors,
    *,
    tolerance_angstrom: float = 2.0,
    movebadrandom: bool = False,
    seed: int | None = None,
    center_solute: bool = True,
    working_directory: str | Path | None = None,
    retain_working_files: bool = True,
):
    """Solvate a system using Packmol, replacing OpenFF's ``pack_box()``.

    This is a drop-in replacement for
    ``openff.interchange.components._packmol.pack_box()`` that handles the
    PDB CONECT-record overflow affecting systems with >99,999 atoms.
    OpenMM's PDB writer hex-encodes atom indices beyond 5 digits, which
    Packmol's fixed-width Fortran parser cannot read.

    The fix has two layers:

    1. **Strip CONECT records** from the solute PDB after writing.
       Packmol only needs atomic coordinates; bond topology is already
       captured in the OpenFF ``Topology`` object.
    2. **``ignore_conect`` keyword** in the Packmol input file as a
       belt-and-suspenders safeguard.

    Parameters
    ----------
    molecules : list[openff.toolkit.Molecule]
        Solvent molecule types (water, ions, co-solvents).
    number_of_copies : list[int]
        Number of copies of each molecule type (parallel to *molecules*).
    solute : openff.toolkit.Topology
        Topology to solvate (protein + polymers + substrate).
    box_vectors : openff.units.Quantity
        Box vectors with shape (3, 3) and length units (e.g. nanometers).
    tolerance_angstrom : float, optional
        Packmol tolerance in Angstrom (default 2.0).
    movebadrandom : bool, optional
        Pass the ``movebadrandom`` keyword to Packmol (default ``False``).
    seed : int or None, optional
        Packmol random seed (see :func:`build_packmol_input`).
    center_solute : bool, optional
        Re-centre the solute at the centre of the periodic brick before
        packing (default ``True``).  Pass ``False`` when the solute is a
        topology that was already framed in *this* brick — for example the
        output of :func:`pack_polymers`, whose chains are packed against the
        brick faces.  Re-centring such a topology by its centre of geometry
        shifts it as a rigid body and pushes atoms back out of the brick.
    working_directory : str, Path, or None, optional
        Directory for Packmol input/output files.  A temporary directory is
        created when ``None``.
    retain_working_files : bool, optional
        Keep working files after the run (default ``True``).

    Returns
    -------
    openff.toolkit.Topology
        Solvated topology with solute + solvent and box vectors set.

    Raises
    ------
    SolvationClashError
        If any solvent atom ends up closer than ``0.5 * tolerance_angstrom``
        to the solute after assembly (solute/solvent frame mismatch).
    PeriodicImageClashError
        If any atom of the solvated system ends up closer than
        ``0.5 * tolerance_angstrom`` to one of its own periodic images.
    """
    import numpy as np
    from openff.packmol._packmol import (
        _center_topology_at,
        _compute_brick_from_box_vectors,
        _create_molecule_pdbs,
        _create_solute_pdb,
        _load_positions,
    )
    from openff.toolkit import Topology
    from openff.units import Quantity

    # --- resolve working directory ---
    _temporary = False
    if working_directory is None:
        working_directory = Path(tempfile.mkdtemp())
        _temporary = True
    else:
        working_directory = Path(working_directory).resolve()
        working_directory.mkdir(parents=True, exist_ok=True)

    # --- compute brick dimensions ---
    brick_size = _compute_brick_from_box_vectors(box_vectors)
    box_size_angstrom = np.asarray(brick_size.m_as("angstrom"), dtype=float)

    # --- center solute in the brick ---
    # Skipped when the caller already framed the topology in this brick
    # (polymer packing does), because a centre-of-geometry shift of an
    # already-framed system pushes atoms back out through the brick faces.
    if center_solute:
        centered_solute = _center_topology_at("BRICK", solute, box_vectors, brick_size)
    else:
        logger.info("Solute is already framed in the periodic brick; skipping re-centring.")
        centered_solute = solute

    # detect whether PBC is usable (rectangular box + packmol >= 20.15.0)
    _use_pbc = _check_pbc_available(box_vectors)

    original_cwd = Path.cwd()
    try:
        os.chdir(working_directory)

        # write PDB files
        solute_pdb = _create_solute_pdb(centered_solute, box_vectors)
        molecule_pdbs = _create_molecule_pdbs(molecules)

        # Strip CONECT records from the solute PDB.
        # OpenMM hex-encodes atom indices > 99999, which Packmol cannot parse.
        if solute_pdb is not None:
            n_stripped = _strip_conect_records(solute_pdb)
            if n_stripped > 0:
                logger.info(
                    "Stripped %d CONECT records from solute PDB (%d atoms)",
                    n_stripped,
                    centered_solute.n_atoms,
                )

        # build and run packmol
        _ignore_conect = _check_ignore_conect_supported()
        input_text = build_packmol_input(
            molecule_pdb_paths=molecule_pdbs,
            molecule_counts=number_of_copies,
            box_size_angstrom=box_size_angstrom,
            tolerance_angstrom=tolerance_angstrom,
            solute_pdb_path=solute_pdb,
            use_pbc=_use_pbc,
            movebadrandom=movebadrandom,
            ignore_conect=_ignore_conect,
            seed=seed,
        )

        output_path = run_packmol(
            input_text=input_text,
            working_directory=working_directory,
            retain_working_files=True,  # always keep; we clean up below
        )

        positions = _load_positions(str(output_path.name))

    finally:
        os.chdir(original_cwd)

    # --- assemble output topology ---
    added_molecules = []
    for mol, n in zip(molecules, number_of_copies):
        added_molecules.extend([mol] * n)
    solvent_topology = Topology.from_molecules(added_molecules)

    n_solute_atoms = len(positions) - solvent_topology.n_atoms
    solvent_topology.set_positions(Quantity(positions[n_solute_atoms:], "angstrom"))

    # Use the *centered* solute so that solute coordinates match the
    # Packmol frame of reference (waters were placed around the centered
    # copy).  Previous code incorrectly used the original un-centered
    # solute, creating a spatial mismatch between solute and solvent.
    if solute is not None:
        solvated_topology = centered_solute + solvent_topology
    else:
        solvated_topology = solvent_topology

    solvated_topology.box_vectors = box_vectors

    if solute is not None:
        _assert_solute_solvent_separation(
            solvated_topology,
            centered_solute.n_atoms,
            tolerance_angstrom=tolerance_angstrom,
            label="solvent",
        )

    _assert_periodic_image_separation(
        solvated_topology,
        box_vectors,
        tolerance_angstrom=tolerance_angstrom,
        label="solvated system",
    )

    if _temporary and not retain_working_files:
        shutil.rmtree(working_directory, ignore_errors=True)

    return solvated_topology


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def solute_sphere_constraint(
    solute,
    *,
    padding_angstrom: float,
) -> "NDArray":
    """Spherical confinement region around a solute, in Angstrom.

    The sphere is centred on the solute's centre of geometry — which is where
    :func:`openff.packmol._packmol._center_topology_at` puts it in the brick —
    and its radius is the circumradius of the solute's bounding box plus
    *padding_angstrom*.  The radius therefore depends only on the solute's
    shape, never on where the chains happen to land, so replicates of one
    condition share it.

    Parameters
    ----------
    solute : openff.toolkit.Topology
        Solute topology carrying positions (typically already brick-centred).
    padding_angstrom : float
        Extra room for the chains beyond the solute's bounding-box
        circumradius, in Angstrom.

    Returns
    -------
    NDArray
        ``[cx, cy, cz, radius]`` in Angstrom, ready for
        :func:`build_packmol_input`.
    """
    import numpy as np

    positions = np.asarray(solute.get_positions().m_as("angstrom"), dtype=float).reshape(-1, 3)
    center = positions.mean(axis=0)
    extent = positions.max(axis=0) - positions.min(axis=0)
    radius = 0.5 * float(np.linalg.norm(extent)) + float(padding_angstrom)
    return np.array([center[0], center[1], center[2], radius], dtype=float)


def _solute_bbox_exclusion(
    centered_solute,
    box_size_angstrom: "NDArray",
    tolerance_angstrom: float,
    molecules: list,
) -> "NDArray":
    """Legacy rectangular shell: solute bounding box inflated by the tolerance.

    Logs the resulting shell thickness and warns when it is thinner than the
    largest polymer diameter, in which case Packmol is unlikely to converge.
    """
    import numpy as np

    from polyzymd.utils.boxvectors import get_topology_bbox_bounds

    min_coords, max_coords = get_topology_bbox_bounds(centered_solute)
    inner_exclusion_box = np.array(
        [
            min_coords[0] - tolerance_angstrom,
            min_coords[1] - tolerance_angstrom,
            min_coords[2] - tolerance_angstrom,
            max_coords[0] + tolerance_angstrom,
            max_coords[1] + tolerance_angstrom,
            max_coords[2] + tolerance_angstrom,
        ]
    )

    # Effective packing box (non-PBC) runs from 0 to ``box_size - tolerance``.
    effective_box = np.asarray(box_size_angstrom, dtype=float) - tolerance_angstrom
    shell_lo = inner_exclusion_box[:3]
    shell_hi = effective_box - inner_exclusion_box[3:]
    min_shell = float(min(np.min(shell_lo), np.min(shell_hi)))

    logger.info(
        "Protein bbox (A): [%.1f, %.1f, %.1f] to [%.1f, %.1f, %.1f]", *min_coords, *max_coords
    )
    logger.info("Exclusion box (A): [%.1f, %.1f, %.1f] to [%.1f, %.1f, %.1f]", *inner_exclusion_box)
    logger.info(
        "Shell thickness lo (A): [%.1f, %.1f, %.1f]  hi: [%.1f, %.1f, %.1f]  min: %.1f",
        *shell_lo,
        *shell_hi,
        min_shell,
    )

    max_diameter = max(_max_molecule_diameter_angstrom(mol) for mol in molecules)
    if min_shell < max_diameter:
        deficit_nm = (max_diameter - min_shell) / 10.0
        logger.warning(
            "Polymer shell thickness (%.1f A) is less than the largest "
            "polymer diameter (%.1f A). Packing may fail or produce poor "
            "results. Consider increasing packing.padding by at least "
            "%.1f nm, or disable exclude_solute_bbox.",
            min_shell,
            max_diameter,
            deficit_nm,
        )
    return inner_exclusion_box


def _max_molecule_diameter_angstrom(mol) -> float:
    """Return the maximum internal distance of a molecule in Angstrom.

    This is the largest pairwise distance between any two atoms in the
    molecule's first conformer — effectively the molecule's "diameter".
    Used to check whether the polymer shell is thick enough to contain
    the molecule.

    Parameters
    ----------
    mol : openff.toolkit.Molecule
        Molecule with at least one conformer.

    Returns
    -------
    float
        Maximum pairwise distance in Angstrom, or 0.0 if the molecule
        has fewer than 2 atoms.
    """
    import numpy as np

    coords = mol.conformers[0].m_as("angstrom")
    if len(coords) <= 1:
        return 0.0
    diffs = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dists_sq = np.sum(diffs**2, axis=-1)
    return float(np.sqrt(np.max(dists_sq)))


def _strip_conect_records(pdb_path: str | Path) -> int:
    """Remove CONECT lines from a PDB file in-place.

    Parameters
    ----------
    pdb_path : str or Path
        Path to the PDB file to modify.

    Returns
    -------
    int
        Number of CONECT lines removed.
    """
    pdb_path = Path(pdb_path)
    lines = pdb_path.read_text().splitlines(keepends=True)
    filtered = [line for line in lines if not line.startswith("CONECT")]
    n_removed = len(lines) - len(filtered)
    if n_removed > 0:
        pdb_path.write_text("".join(filtered))
    return n_removed


def _check_ignore_conect_supported() -> bool:
    """Return True if the installed Packmol supports ``ignore_conect``.

    The ``ignore_conect`` keyword was added after version 21.1.3.
    We check for version > 21.1.3.  If the version cannot be determined,
    return False (safe default — CONECT stripping handles the issue).
    """
    try:
        from openff.packmol._packmol import _get_packmol_version
        from packaging.version import Version

        return _get_packmol_version() > Version("21.1.3")
    except (ImportError, RuntimeError, TypeError, ValueError) as exc:
        logger.warning(
            "Could not determine Packmol ignore_conect support (%s). "
            "Continuing with CONECT stripping fallback; verify the Packmol executable "
            "is installed if packing fails.",
            exc,
        )
        return False


def _check_pbc_available(box_vectors) -> bool:
    """Return True if the box is rectangular and Packmol supports PBC."""
    try:
        import numpy as np
        from openff.packmol._packmol import _get_packmol_version
        from packaging.version import Version

        box_arr = np.asarray(box_vectors.m)
        is_rectangular = bool(np.all(box_arr == np.diag(np.diagonal(box_arr))))
        if not is_rectangular:
            return False
        return _get_packmol_version() >= Version("20.15.0")
    except (AttributeError, ImportError, RuntimeError, TypeError, ValueError) as exc:
        logger.warning(
            "Could not determine Packmol PBC support (%s). "
            "Continuing without Packmol PBC mode; verify box vectors and Packmol "
            "installation if packing fails.",
            exc,
        )
        return False
