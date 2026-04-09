"""Helpers to patch OpenFF topology metadata to PolyzyMD PDB convention

This module is used after solvation and before Interchange creation so that
downstream OpenMM PDB writes use stable chain, residue, and atom labeling
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

_CONJUGATED_RESNAME_MAP = {
    "SB1": "SBM",
    "SB2": "SBM",
    "EG1": "EGP",
    "EG2": "EGP",
}


def _parse_protein_records(crystal_pdb_path: Path) -> list[tuple[str, str, str]]:
    """Parse canonical atom records from the crystal protein PDB

    Parameters
    ----------
    crystal_pdb_path : Path
        Path to the crystal protein PDB used as canonical naming source

    Returns
    -------
    list[tuple[str, str, str]]
        Ordered ``(atom_name, residue_name, residue_number)`` tuples for all
        ATOM/HETATM records in file order
    """

    records: list[tuple[str, str, str]] = []
    with crystal_pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            residue_number = line[22:26].strip()
            records.append((atom_name, residue_name, residue_number))

    if not records:
        raise ValueError(f"No ATOM/HETATM records found in {crystal_pdb_path}")

    return records


def _canonicalize_conjugated_resname(
    raw_resname: str,
    assembled_resid: int | None,
    crosslink_resids: set[int],
) -> str:
    """Map conjugated polymer residue names to canonical POC names

    Parameters
    ----------
    raw_resname : str
        Original residue name from the assembled topology metadata
    assembled_resid : int | None
        Original assembled-PDB residue number for this residue
    crosslink_resids : set[int]
        Assembled-PDB residue numbers that represent crosslinked NHS residues

    Returns
    -------
    str
        Canonicalized residue name
    """

    if raw_resname in {"NH1", "NH2"}:
        if assembled_resid is not None and assembled_resid in crosslink_resids:
            return "NHX"
        return "NHS"

    return _CONJUGATED_RESNAME_MAP.get(raw_resname, raw_resname)


def _chain_for_solvent_offset(offset: int) -> str:
    """Convert solvent chain offset to single-character chain identifier

    Parameters
    ----------
    offset : int
        Zero-based solvent chain offset where 0 maps to chain D

    Returns
    -------
    str
        Chain identifier
    """

    chain_code = ord("D") + offset
    if chain_code > ord("Z"):
        raise ValueError("Solvent chain overflow beyond chain Z")
    return chr(chain_code)


def apply_chain_convention(
    solvated_topology: Any,
    component_metadata: dict[str, Any],
    crystal_pdb_path: Path,
    crosslink_resids: tuple[int, ...],
) -> dict[str, Any]:
    """Patch an OpenFF solvated topology to PolyzyMD chain convention

    Parameters
    ----------
    solvated_topology : Any
        Solvated OpenFF ``Topology`` object to patch in place
    component_metadata : dict[str, Any]
        Metadata dictionary with component atom indices from the conjugation POC
    crystal_pdb_path : Path
        Path to the crystal protein PDB used to recover canonical protein naming
    crosslink_resids : tuple[int, ...]
        Assembled-PDB residue numbers corresponding to crosslinked NHS residues

    Returns
    -------
    dict[str, Any]
        Summary counts for patched components and used chains
    """

    canonical_records = _parse_protein_records(Path(crystal_pdb_path))
    protein_atom_indices = list(component_metadata["protein"]["atom_indices"])
    if len(canonical_records) != len(protein_atom_indices):
        raise ValueError(
            "Protein atom count mismatch between crystal PDB and metadata: "
            f"{len(canonical_records)} != {len(protein_atom_indices)}"
        )

    molecules = list(solvated_topology.molecules)
    if not molecules:
        raise ValueError("Solvated topology has no molecules")

    conjugate_molecule = molecules[0]
    protein_index_to_record = dict(zip(protein_atom_indices, canonical_records, strict=True))
    protein_index_set = set(protein_atom_indices)
    crosslink_set = set(crosslink_resids)
    chains_used: set[str] = set()

    n_protein_atoms = 0
    n_conjugated_polymer_atoms = 0

    conjugated_resid_map: dict[str, int] = {}
    last_old_resid: str | None = None
    next_polymer_resid = 1

    for atom in conjugate_molecule.atoms:
        atom_idx = atom.molecule_atom_index
        if atom_idx in protein_index_set:
            atom_name, residue_name, residue_number = protein_index_to_record[atom_idx]
            atom.metadata["chain_id"] = "A"
            atom.metadata["residue_name"] = residue_name
            atom.metadata["residue_number"] = residue_number
            atom.name = atom_name
            n_protein_atoms += 1
            chains_used.add("A")
            continue

        old_resid_text = str(atom.metadata.get("residue_number", "")).strip()
        if old_resid_text != last_old_resid:
            if old_resid_text not in conjugated_resid_map:
                conjugated_resid_map[old_resid_text] = next_polymer_resid
                next_polymer_resid += 1
            last_old_resid = old_resid_text

        old_resname = str(atom.metadata.get("residue_name", "")).strip()
        assembled_resid = int(old_resid_text) if old_resid_text.isdigit() else None
        canonical_resname = _canonicalize_conjugated_resname(
            old_resname, assembled_resid, crosslink_set
        )

        atom.metadata["chain_id"] = "C"
        atom.metadata["residue_name"] = canonical_resname
        atom.metadata["residue_number"] = str(conjugated_resid_map[old_resid_text])
        n_conjugated_polymer_atoms += 1
        chains_used.add("C")

    polymer_resid_counter = next_polymer_resid - 1

    n_free_polymer_molecules = len(component_metadata.get("free_polymers", []))
    n_free_polymer_atoms = 0

    for free_molecule in molecules[1 : 1 + n_free_polymer_molecules]:
        local_resid_map: dict[str, int] = {}
        for atom in free_molecule.atoms:
            old_resid_text = str(atom.metadata.get("residue_number", "")).strip()
            if old_resid_text not in local_resid_map:
                polymer_resid_counter += 1
                local_resid_map[old_resid_text] = polymer_resid_counter

            atom.metadata["chain_id"] = "C"
            atom.metadata["residue_number"] = str(local_resid_map[old_resid_text])
            n_free_polymer_atoms += 1
            chains_used.add("C")

    n_water = 0
    n_na = 0
    n_cl = 0
    solvent_chain_offset = 0
    solvent_resid_counter = 0

    for solvent_molecule in molecules[1 + n_free_polymer_molecules :]:
        if solvent_resid_counter >= 9999:
            solvent_chain_offset += 1
            solvent_resid_counter = 1
        else:
            solvent_resid_counter += 1

        chain_id = _chain_for_solvent_offset(solvent_chain_offset)
        residue_number = str(solvent_resid_counter)
        chains_used.add(chain_id)

        atoms = list(solvent_molecule.atoms)
        atomic_numbers = [atom.atomic_number for atom in atoms]

        is_water = len(atoms) == 3 and atomic_numbers.count(8) == 1 and atomic_numbers.count(1) == 2
        is_na = len(atoms) == 1 and atomic_numbers[0] == 11
        is_cl = len(atoms) == 1 and atomic_numbers[0] == 17

        if is_water:
            oxygen_atom = next(atom for atom in atoms if atom.atomic_number == 8)
            hydrogen_atoms = sorted(
                (atom for atom in atoms if atom.atomic_number == 1),
                key=lambda atom: atom.molecule_atom_index,
            )

            oxygen_atom.metadata["chain_id"] = chain_id
            oxygen_atom.metadata["residue_name"] = "HOH"
            oxygen_atom.metadata["residue_number"] = residue_number
            oxygen_atom.name = "O"

            for idx, hydrogen_atom in enumerate(hydrogen_atoms, start=1):
                hydrogen_atom.metadata["chain_id"] = chain_id
                hydrogen_atom.metadata["residue_name"] = "HOH"
                hydrogen_atom.metadata["residue_number"] = residue_number
                hydrogen_atom.name = f"H{idx}"

            n_water += 1
            continue

        if is_na:
            atom = atoms[0]
            atom.metadata["chain_id"] = chain_id
            atom.metadata["residue_name"] = "NA"
            atom.metadata["residue_number"] = residue_number
            atom.name = "NA"
            n_na += 1
            continue

        if is_cl:
            atom = atoms[0]
            atom.metadata["chain_id"] = chain_id
            atom.metadata["residue_name"] = "CL"
            atom.metadata["residue_number"] = residue_number
            atom.name = "CL"
            n_cl += 1
            continue

        for atom in atoms:
            atom.metadata["chain_id"] = chain_id
            atom.metadata["residue_number"] = residue_number

    return {
        "n_protein_atoms": n_protein_atoms,
        "n_conjugated_polymer_atoms": n_conjugated_polymer_atoms,
        "n_free_polymer_atoms": n_free_polymer_atoms,
        "n_water": n_water,
        "n_na": n_na,
        "n_cl": n_cl,
        "chains_used": sorted(chains_used),
    }
