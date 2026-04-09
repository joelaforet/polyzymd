"""Post-vacuum-equilibration free polymer placement using Packmol.

This module places free polymers around a relaxed protein-polymer conjugate
using the project's existing Packmol shell-packing utilities. It is called
by run_conjugate_simulation.py after the conjugation notebook has produced
a vacuum-equilibrated conjugate.
"""

from __future__ import annotations

import copy
import logging
import subprocess
from pathlib import Path

logger = logging.getLogger("conjugate_sim.free_placement")


def place_free_polymers(
    conjugate_off,
    free_polymer_offs: list,
    *,
    padding_nm_ladder: tuple[float, ...] = (2.0, 3.0, 4.0, 5.0),
    tolerance_angstrom: float = 2.0,
    movebadrandom: bool = True,
    nloop: int = 500,
    working_directory: Path | None = None,
) -> list:
    """Place free polymers around the relaxed conjugate using Packmol.

    Parameters
    ----------
    conjugate_off : openff.toolkit.Molecule
        Relaxed conjugate molecule (protein + conjugated polymers) with
        equilibrated conformer.
    free_polymer_offs : list[openff.toolkit.Molecule]
        Unplaced free polymer molecules (with conformers from generation).
    padding_nm_ladder : tuple[float, ...]
        Padding values to try in order. If packing fails at one padding,
        the next (larger) value is tried.
    tolerance_angstrom : float
        Packmol tolerance in Angstrom.
    movebadrandom : bool
        Pass movebadrandom to Packmol.
    nloop : int
        Max GENCAN optimization loops.
    working_directory : Path or None
        Working directory for Packmol files.

    Returns
    -------
    list[openff.toolkit.Molecule]
        Free polymer molecules with updated conformers (placed coordinates).
    """
    import numpy as np
    from openff.toolkit import Molecule
    from openff.units import Quantity, unit

    from polyzymd.utils.packmol import pack_polymers

    # Build solute topology from conjugate
    solute_topology = conjugate_off.to_topology()

    # Deduplicate free polymers by SMILES for Packmol input
    # Packmol takes unique molecule types + counts
    smiles_to_mol: dict[str, Molecule] = {}
    smiles_order: list[str] = []
    smiles_counts: dict[str, int] = {}
    smiles_indices: dict[str, list[int]] = {}

    for i, mol in enumerate(free_polymer_offs):
        smi = mol.to_smiles()
        if smi not in smiles_to_mol:
            smiles_to_mol[smi] = mol
            smiles_order.append(smi)
            smiles_counts[smi] = 0
            smiles_indices[smi] = []
        smiles_counts[smi] += 1
        smiles_indices[smi].append(i)

    unique_mols = [smiles_to_mol[s] for s in smiles_order]
    counts = [smiles_counts[s] for s in smiles_order]

    logger.info(
        "Free polymer placement: %d total polymers, %d unique types",
        len(free_polymer_offs),
        len(unique_mols),
    )

    # Try padding ladder
    last_error = None
    packed_topology = None

    for padding_nm in padding_nm_ladder:
        # Compute box vectors from solute bbox + padding
        solute_positions = solute_topology.get_positions()
        if solute_positions is None:
            raise RuntimeError("Conjugate topology has no positions")

        pos_nm = solute_positions.m_as(unit.nanometer)
        min_coords = pos_nm.min(axis=0)
        max_coords = pos_nm.max(axis=0)
        box_lengths = max_coords - min_coords + 2 * padding_nm

        box_vectors = Quantity(
            np.diag(box_lengths),
            unit.nanometer,
        )

        logger.info(
            "Trying padding=%.1f nm → box %.1f x %.1f x %.1f nm",
            padding_nm,
            box_lengths[0],
            box_lengths[1],
            box_lengths[2],
        )

        try:
            packed_topology = pack_polymers(
                molecules=unique_mols,
                number_of_copies=counts,
                solute=solute_topology,
                box_vectors=box_vectors,
                tolerance_angstrom=tolerance_angstrom,
                movebadrandom=movebadrandom,
                nloop=nloop,
                working_directory=working_directory,
            )
            logger.info("Packing succeeded with padding=%.1f nm", padding_nm)
            break
        except (subprocess.CalledProcessError, RuntimeError, FileNotFoundError, ValueError) as e:
            last_error = e
            logger.warning("Packing failed with padding=%.1f nm: %s", padding_nm, e)
            continue

    if packed_topology is None:
        raise RuntimeError(
            "Free polymer placement failed after trying padding values "
            f"{padding_nm_ladder}: {last_error}"
        )

    # Extract placed free polymer coordinates from packed topology
    # packed_topology contains solute atoms first, then packed polymer atoms
    # pack_polymers sorts molecules by count descending, so we need to
    # reconstruct the original ordering
    packed_positions = packed_topology.get_positions()
    if packed_positions is None:
        raise RuntimeError("Packed topology has no positions")

    packed_coords_angstrom = packed_positions.m_as(unit.angstrom)
    n_solute = solute_topology.n_atoms
    polymer_coords = packed_coords_angstrom[n_solute:]

    # pack_polymers sorts by count descending. We need to map back to
    # original order. Build the sorted order that pack_polymers used
    paired = list(zip(smiles_order, counts, strict=False))
    paired_sorted = sorted(enumerate(paired), key=lambda x: x[1][1], reverse=True)

    # Build placed molecules in the sorted order first, then reorder
    placed_by_sorted_group: list[tuple[int, list]] = []
    coord_offset = 0
    for _sort_idx, (_orig_group_idx, (smi, count)) in enumerate(paired_sorted):
        group_placed = []
        mol_template = smiles_to_mol[smi]
        n_atoms = mol_template.n_atoms
        for _copy in range(count):
            mol_copy = Molecule(mol_template)
            copy_coords = polymer_coords[coord_offset : coord_offset + n_atoms]
            mol_copy._conformers = [Quantity(copy_coords, unit.angstrom)]
            group_placed.append(mol_copy)
            coord_offset += n_atoms
        placed_by_sorted_group.append((_orig_group_idx, group_placed))

    # Rebuild in original smiles group order
    placed_by_original_group: dict[int, list] = {}
    for orig_group_idx, mols in placed_by_sorted_group:
        placed_by_original_group[orig_group_idx] = mols

    # Map back to original free_polymer_offs order
    placed_free_polymer_offs = [None] * len(free_polymer_offs)
    for group_idx, smi in enumerate(smiles_order):
        orig_indices = smiles_indices[smi]
        group_mols = placed_by_original_group[group_idx]
        for local_idx, orig_idx in enumerate(orig_indices):
            placed_free_polymer_offs[orig_idx] = group_mols[local_idx]

    # Verify all placed
    assert all(m is not None for m in placed_free_polymer_offs), "Some free polymers not placed"

    # Also update conjugate_off conformer with the centered coordinates
    # from packed topology (pack_polymers centers the solute in the box)
    conjugate_coords = packed_coords_angstrom[:n_solute]
    conjugate_off._conformers = [Quantity(conjugate_coords, unit.angstrom)]
    logger.info("Updated conjugate coordinates to Packmol-centered frame")

    return placed_free_polymer_offs


def finalize_component_metadata(
    base_metadata: dict,
    conjugate_off,
    placed_free_polymer_offs: list,
) -> dict:
    """Add free polymer global indices to component metadata.

    Parameters
    ----------
    base_metadata : dict
        Metadata from the notebook (conjugate-only).
    conjugate_off : openff.toolkit.Molecule
        Conjugate molecule.
    placed_free_polymer_offs : list[openff.toolkit.Molecule]
        Placed free polymer molecules.

    Returns
    -------
    dict
        Complete metadata with free polymer indices added.
    """
    metadata = copy.deepcopy(base_metadata)

    free_polymer_ranges = []
    free_offset = conjugate_off.n_atoms
    all_free_heavy = []

    for fi, free_mol in enumerate(placed_free_polymer_offs):
        n_atoms = free_mol.n_atoms
        atom_indices = list(range(free_offset, free_offset + n_atoms))
        heavy_indices = [
            idx
            for idx, local_idx in zip(atom_indices, range(n_atoms), strict=False)
            if free_mol.atom(local_idx).atomic_number > 1
        ]
        free_polymer_ranges.append(
            {
                "polymer_id": f"free_{fi}",
                "atom_indices": atom_indices,
                "heavy_atom_indices": heavy_indices,
            }
        )
        all_free_heavy.extend(heavy_indices)
        free_offset += n_atoms

    all_free_heavy.sort()

    metadata["n_free_polymer_atoms"] = sum(m.n_atoms for m in placed_free_polymer_offs)
    metadata["free_polymers"] = free_polymer_ranges
    metadata["restraint_groups"]["free_polymer_heavy"] = all_free_heavy

    logger.info(
        "Finalized metadata: %d free polymers, %d free heavy atoms",
        len(placed_free_polymer_offs),
        len(all_free_heavy),
    )

    return metadata
