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

    lines: list[str] = [
        f"tolerance {tolerance_angstrom:f}",
        "filetype pdb",
        f"output {_PACKMOL_OUTPUT_FILE}",
        "",
    ]

    if movebadrandom:
        lines.append("movebadrandom")
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
            try:
                result = subprocess.check_output(
                    packmol_binary,
                    stdin=fh,
                    stderr=subprocess.STDOUT,
                )
            except subprocess.CalledProcessError as exc:
                error_log_path.write_bytes(exc.output)
                raise RuntimeError(
                    f"Packmol exited with return code {exc.returncode}. "
                    f"See {error_log_path} for details."
                ) from exc
    finally:
        os.chdir(original_cwd)

    stdout = result.decode("utf-8", errors="replace")
    if "Success!" not in stdout:
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
    working_directory: str | Path | None = None,
    retain_working_files: bool = True,
):
    """Pack polymer chains around a fixed solute using Packmol.

    This is a drop-in replacement for the OpenFF ``pack_box()`` call in
    :meth:`~polyzymd.builders.system_builder.SystemBuilder.pack_polymers`.
    It adds support for the ``movebadrandom`` Packmol keyword.

    Parameters
    ----------
    molecules : list[openff.toolkit.Molecule]
        Unique molecule types to pack.
    number_of_copies : list[int]
        Number of copies of each molecule type (parallel to *molecules*).
    solute : openff.toolkit.Topology
        Fixed topology (protein + substrate) placed at the origin.
    box_vectors : openff.units.Quantity
        Box vectors with shape (3, 3) and length units (e.g. nanometers).
    tolerance_angstrom : float, optional
        Packmol tolerance in Angstrom (default 2.0).
    movebadrandom : bool, optional
        Pass the ``movebadrandom`` keyword to Packmol (default ``False``).
    working_directory : str, Path, or None, optional
        Directory for Packmol input/output files.  A temporary directory is
        created when ``None``.
    retain_working_files : bool, optional
        Keep working files after the run (default ``True``).

    Returns
    -------
    openff.toolkit.Topology
        Packed topology with solute + all polymer chains and box vectors set.
    """
    import numpy as np
    from openff.interchange.components._packmol import (
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
        working_directory = Path(working_directory)
        working_directory.mkdir(parents=True, exist_ok=True)

    # --- compute brick dimensions ---
    brick_size = _compute_brick_from_box_vectors(box_vectors)
    box_size_angstrom = np.asarray(brick_size.m_as("angstrom"), dtype=float)

    # --- center solute in the brick ---
    centered_solute = _center_topology_at("BRICK", solute, box_vectors, brick_size)

    # detect whether PBC is usable (rectangular box + packmol >= 20.15.0)
    _use_pbc = _check_pbc_available(box_vectors)

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
    packed_topology = Topology.from_molecules(added_molecules)

    n_solute_atoms = len(positions) - packed_topology.n_atoms
    packed_topology.set_positions(Quantity(positions[n_solute_atoms:], "angstrom"))

    if solute is not None:
        packed_topology = solute + packed_topology

    packed_topology.box_vectors = box_vectors

    if _temporary and not retain_working_files:
        shutil.rmtree(working_directory, ignore_errors=True)

    return packed_topology


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _check_pbc_available(box_vectors) -> bool:
    """Return True if the box is rectangular and Packmol supports PBC."""
    try:
        import numpy as np
        from openff.interchange.components._packmol import _get_packmol_version
        from packaging.version import Version

        box_arr = np.asarray(box_vectors.m)
        is_rectangular = bool(np.all(box_arr == np.diag(np.diagonal(box_arr))))
        if not is_rectangular:
            return False
        return _get_packmol_version() >= Version("20.15.0")
    except Exception:
        return False
