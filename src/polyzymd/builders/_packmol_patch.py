"""
Patched PACKMOL wrapper with support for optimization keywords.

This module extends OpenFF Interchange's PACKMOL wrapper to support
additional keywords that can speed up convergence for complex systems:

- movebadrandom: Move badly-placed molecules to random positions
- discale: Distance tolerance scale factor for local optimization
- maxit: Maximum iterations per GENCAN optimization loop
- nloop: Maximum number of optimization loops
- movefrac: Fraction of molecules to move in bad-molecule heuristic
- seed: Random seed for reproducibility

These keywords are documented in the PACKMOL manual:
https://m3g.github.io/packmol/userguide.shtml

Example usage:
    >>> from polyzymd.builders._packmol_patch import pack_box_with_optimization
    >>>
    >>> topology = pack_box_with_optimization(
    ...     molecules=[water, dmso],
    ...     number_of_copies=[1000, 100],
    ...     solute=protein_topology,
    ...     box_vectors=box_vecs,
    ...     optimization_kwargs={
    ...         "movebadrandom": True,
    ...         "discale": 1.5,
    ...         "nloop": 200,
    ...     }
    ... )
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray
from openff.toolkit import Molecule, Quantity, RDKitToolkitWrapper, Topology
from openff.utilities.utilities import temporary_cd

from openff.interchange.exceptions import PACKMOLRuntimeError, PACKMOLValueError

LOGGER = logging.getLogger(__name__)


@dataclass
class PackmolOptimizationOptions:
    """PACKMOL optimization options for faster convergence.

    Attributes:
        movebadrandom: If True, move badly-placed molecules to random positions
            instead of near well-packed molecules. Helps with complex restraints.
        discale: Distance tolerance scale factor for local optimization.
            Setting to 1.5 can help convergence.
        maxit: Maximum iterations per GENCAN optimization loop (default: 20).
        nloop: Maximum number of optimization loops before giving up.
        movefrac: Fraction of molecules to move in bad-molecule heuristic.
            For large systems, smaller values (e.g., 0.01) may work better.
        seed: Random seed. Use -1 for automatic seed from system time.
    """

    movebadrandom: bool = False
    discale: Optional[float] = None
    maxit: Optional[int] = None
    nloop: Optional[int] = None
    movefrac: Optional[float] = None
    seed: Optional[int] = None

    def to_packmol_lines(self) -> List[str]:
        """Generate PACKMOL input file lines for these options.

        Returns:
            List of lines to insert into PACKMOL input file.
        """
        lines = []

        if self.movebadrandom:
            lines.append("movebadrandom")

        if self.discale is not None:
            lines.append(f"discale {self.discale}")

        if self.maxit is not None:
            lines.append(f"maxit {self.maxit}")

        if self.nloop is not None:
            lines.append(f"nloop {self.nloop}")

        if self.movefrac is not None:
            lines.append(f"movefrac {self.movefrac}")

        if self.seed is not None:
            lines.append(f"seed {self.seed}")

        return lines


def _find_packmol() -> Optional[str]:
    """Find the packmol binary."""
    return shutil.which("packmol")


def _build_input_file_with_optimization(
    molecule_file_names: List[str],
    molecule_counts: List[int],
    structure_to_solvate: Optional[str],
    box_size: Quantity,
    tolerance: Quantity,
    optimization: Optional[PackmolOptimizationOptions] = None,
) -> tuple[str, str]:
    """
    Construct the PACKMOL input file with optimization keywords.

    This is a patched version of OpenFF's _build_input_file() that supports
    additional PACKMOL keywords for faster convergence.

    Parameters
    ----------
    molecule_file_names
        The paths to the molecule PDB files.
    molecule_counts
        The number of each molecule to add.
    structure_to_solvate
        The path to the structure to solvate (solute).
    box_size
        The lengths of each side of the rectangular brick box.
    tolerance
        The PACKMOL convergence tolerance.
    optimization
        Optional optimization settings.

    Returns
    -------
    tuple[str, str]
        Paths to the input file and expected output file.
    """
    box_size_ang = (box_size - tolerance).m_as("angstrom")
    tolerance_ang = tolerance.m_as("angstrom")

    # Add the global header options
    output_file_path = "packmol_output.pdb"
    input_lines = [
        f"tolerance {tolerance_ang:f}",
        "filetype pdb",
        f"output {output_file_path}",
    ]

    # Add optimization keywords right after header (before structures)
    if optimization is not None:
        opt_lines = optimization.to_packmol_lines()
        if opt_lines:
            input_lines.extend(opt_lines)
            LOGGER.info(f"PACKMOL optimization keywords: {', '.join(opt_lines)}")

    input_lines.append("")  # Blank line before structures

    # Add the section for the molecule to solvate if provided
    if structure_to_solvate is not None:
        input_lines.extend(
            [
                f"structure {structure_to_solvate}",
                "  number 1",
                "  fixed 0. 0. 0. 0. 0. 0.",
                "end structure",
                "",
            ]
        )

    # Add a section for each type of molecule to add
    for file_name, count in zip(molecule_file_names, molecule_counts):
        if count == 0:
            continue

        input_lines.extend(
            [
                f"structure {file_name}",
                f"  number {count}",
                f"  inside box 0. 0. 0. {box_size_ang[0]} {box_size_ang[1]} {box_size_ang[2]}",
                "end structure",
                "",
            ]
        )

    packmol_input = "\n".join(input_lines)

    # Write PACKMOL input
    packmol_file_name = "packmol_input.txt"
    with open(packmol_file_name, "w") as f:
        f.write(packmol_input)

    return packmol_file_name, output_file_path


def _create_solute_pdb(
    topology: Optional[Topology],
    box_vectors: Quantity,
) -> Optional[str]:
    """Write out the solute topology to PDB so that PACKMOL can read it."""
    if topology is None:
        return None

    from openff.interchange.components._packmol import _wrap_into_brick

    # Copy the topology so we can change it
    topology = Topology(topology)
    # Wrap the positions into the brick representation
    topology.set_positions(
        _wrap_into_brick(
            topology.get_positions(),
            box_vectors,
        ),
    )
    # Write to PDB
    solute_pdb_filename = "_PACKING_SOLUTE.pdb"
    topology.to_file(
        solute_pdb_filename,
        file_format="PDB",
    )
    return solute_pdb_filename


def _create_molecule_pdbs(molecules: List[Molecule]) -> List[str]:
    """Write out PDBs of the molecules so that PACKMOL can read them."""
    pdb_file_names = []
    for index, molecule in enumerate(molecules):
        # Make a copy of the molecule so we don't change the input
        molecule = Molecule(molecule)

        # Generate conformers if they're missing
        if molecule.n_conformers <= 0:
            molecule.generate_conformers(n_conformers=1)
        else:
            # Possible issues writing multi-conformer PDBs with RDKit
            molecule._conformers = [molecule.conformers[0]]

        pdb_file_name = f"_PACKING_MOLECULE{index}.pdb"
        pdb_file_names.append(pdb_file_name)

        molecule.to_file(
            pdb_file_name,
            file_format="PDB",
            toolkit_registry=RDKitToolkitWrapper(),
        )
    return pdb_file_names


def _load_positions(output_file_path: str) -> NDArray:
    """Load positions from PACKMOL output PDB."""
    try:
        return np.asarray(
            [
                [line[31:39], line[39:46], line[47:54]]
                for line in open(output_file_path).readlines()
                if line.startswith("HETATM") or line.startswith("ATOM")
            ],
            dtype=np.float32,
        )
    except Exception as error:
        raise PACKMOLRuntimeError(
            "PACKMOL output could not be parsed by native coordinate parser, "
            "please raise an issue with code reproducing this error.",
        ) from error


def pack_box_with_optimization(
    molecules: List[Molecule],
    number_of_copies: List[int],
    solute: Optional[Topology] = None,
    tolerance: Quantity = Quantity(2.0, "angstrom"),
    box_vectors: Optional[Quantity] = None,
    target_density: Optional[Quantity] = None,
    box_shape: ArrayLike = None,
    center_solute: Union[bool, Literal["BOX_VECS", "ORIGIN", "BRICK"]] = False,
    working_directory: Optional[str] = None,
    retain_working_files: bool = True,
    optimization_kwargs: Optional[Dict[str, Any]] = None,
) -> Topology:
    """
    Run PACKMOL to generate a box containing a mixture of molecules.

    This is a wrapper around OpenFF's pack_box() that supports additional
    PACKMOL optimization keywords for faster convergence.

    Parameters
    ----------
    molecules
        The molecules in the system.
    number_of_copies
        A list of the number of copies of each molecule type.
    solute
        An OpenFF Topology to include in the box.
    tolerance
        The minimum spacing between molecules during packing.
    box_vectors
        The box vectors to fill. Required if target_density not provided.
    target_density
        Target mass density for final system. Required if box_vectors not provided.
    box_shape
        The shape of the simulation box (default: RHOMBIC_DODECAHEDRON).
    center_solute
        How to center the solute in the simulation box.
    working_directory
        The directory in which to generate temporary working files.
    retain_working_files
        If True, retain all working files.
    optimization_kwargs
        Dictionary of PACKMOL optimization options:
        - movebadrandom (bool): Move badly-placed molecules randomly
        - discale (float): Distance tolerance scale factor
        - maxit (int): Max iterations per GENCAN loop
        - nloop (int): Max optimization loops
        - movefrac (float): Fraction of molecules to move
        - seed (int): Random seed (-1 for auto)

    Returns
    -------
    Topology
        An OpenFF Topology with the packed system.

    Raises
    ------
    PACKMOLRuntimeError
        When PACKMOL fails to execute or converge.
    PACKMOLValueError
        When inputs are invalid or inconsistent.
    """
    from openff.interchange.components._packmol import (
        RHOMBIC_DODECAHEDRON,
        _box_from_density,
        _box_vectors_are_in_reduced_form,
        _center_topology_at,
        _compute_brick_from_box_vectors,
        _validate_inputs,
    )

    # Set default box shape
    if box_shape is None:
        box_shape = RHOMBIC_DODECAHEDRON

    box_shape = np.asarray(box_shape)
    if box_shape.shape == (3,):
        box_shape = box_shape * np.identity(3)

    # Find PACKMOL
    packmol_path = _find_packmol()
    if packmol_path is None:
        raise OSError("PACKMOL not found, cannot run pack_box_with_optimization()")

    # Validate inputs
    _validate_inputs(
        molecules,
        number_of_copies,
        solute,
        box_shape,
        box_vectors,
        target_density,
    )

    # Estimate box vectors from mass density if not provided
    if target_density is not None:
        box_vectors = _box_from_density(
            molecules,
            number_of_copies,
            target_density,
            box_shape,
        )

    # If neither box vectors nor density given, take box vectors from solute
    if box_vectors is None:
        box_vectors = solute.box_vectors

    if not _box_vectors_are_in_reduced_form(box_vectors):
        raise PACKMOLValueError(
            "pack_box_with_optimization requires box vectors to be in OpenMM reduced form.\n"
            "See http://docs.openmm.org/latest/userguide/theory/"
            "05_other_features.html#periodic-boundary-conditions",
        )

    # Compute dimensions of equivalent brick
    brick_size = _compute_brick_from_box_vectors(box_vectors)

    # Center the solute
    if center_solute and solute is not None:
        solute = _center_topology_at(
            center_solute,
            solute,
            box_vectors,
            brick_size,
        )

    # Create optimization options object
    optimization = None
    if optimization_kwargs:
        optimization = PackmolOptimizationOptions(**optimization_kwargs)

    # Set up working directory
    temporary_directory = False
    if working_directory is None:
        working_directory = tempfile.mkdtemp()
        temporary_directory = True

    if len(working_directory) > 0:
        os.makedirs(working_directory, exist_ok=True)

    with temporary_cd(working_directory):
        solute_pdb_filename = _create_solute_pdb(solute, box_vectors)

        # Create PDB files for all molecules
        pdb_file_names = _create_molecule_pdbs(molecules)

        # Generate input file with optimization keywords
        input_file_path, output_file_path = _build_input_file_with_optimization(
            pdb_file_names,
            number_of_copies,
            solute_pdb_filename,
            brick_size,
            tolerance,
            optimization=optimization,
        )

        # Run PACKMOL
        with open(input_file_path) as file_handle:
            try:
                result = subprocess.check_output(
                    packmol_path,
                    stdin=file_handle,
                    stderr=subprocess.STDOUT,
                )
            except subprocess.CalledProcessError as error:
                open("packmol_error.log", "w").write(error.stdout.decode("utf-8"))
                raise PACKMOLRuntimeError(
                    f"PACKMOL failed with error code {error.returncode}. "
                    f"Wrote file packmol_error.log in working directory: {working_directory}",
                ) from error

            packmol_succeeded = result.decode("utf-8").find("Success!") > 0

        if not packmol_succeeded:
            raise PACKMOLRuntimeError(
                "PACKMOL did not raise an error code, but 'Success!' not found in output. "
                "Please raise an issue showing how you arrived at this error.",
            )

        positions = _load_positions(output_file_path)

    # Clean up temporary directory if requested
    if temporary_directory and not retain_working_files:
        shutil.rmtree(working_directory)

    # Construct output topology
    added_molecules = []
    for mol, n in zip(molecules, number_of_copies):
        added_molecules.extend([mol] * n)
    topology = Topology.from_molecules(added_molecules)

    # Set positions, skipping the positions from solute
    n_solute_atoms = len(positions) - topology.n_atoms
    topology.set_positions(Quantity(positions[n_solute_atoms:], "angstrom"))

    # Add solute back in with original, unwrapped positions
    if solute is not None:
        topology = solute + topology

    # Set box vectors
    topology.box_vectors = box_vectors

    return topology
