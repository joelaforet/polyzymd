"""Tests for RDKit-first PDB input identity mapping."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.builders.conjugation._linkage import parse_pdb_atom_records
from polyzymd.builders.conjugation.rdkit_inputs import (
    AtomIdentity,
    build_rdkit_input_bundle,
    load_pdb_as_rdkit_input,
)

pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402


def test_load_pdb_as_rdkit_input_keeps_explicit_hydrogens_and_identity(tmp_path: Path):
    """PDB loading should preserve atom order and map explicit hydrogens."""
    pdb_path = tmp_path / "lys_fragment.pdb"
    pdb_path.write_text(
        _pdb_atom(10, "NZ", "LYS", "A", 23, 0.0, 0.0, 0.0, element="N")
        + _pdb_atom(11, "H10", "LYS", "A", 23, 0.0, 1.0, 0.0, element="H")
        + _pdb_atom(12, "H11", "LYS", "A", 23, 1.0, 0.0, 0.0, element="H")
        + "END\n",
        encoding="utf-8",
    )

    bundle = load_pdb_as_rdkit_input(pdb_path)

    assert bundle.mol.GetNumAtoms() == 3
    assert bundle.source_path == pdb_path
    assert bundle.pdb_index_to_rdkit_index == {0: 0, 1: 1, 2: 2}
    assert bundle.rdkit_index_to_pdb_index == {0: 0, 1: 1, 2: 2}
    assert bundle.serial_to_rdkit_index[11] == 1
    assert bundle.atom_identities[1] == AtomIdentity(
        chain_id="A",
        residue_name="LYS",
        residue_number=23,
        atom_name="H10",
        atom_index=1,
        serial=11,
        element="H",
        rdkit_index=1,
    )


def test_build_rdkit_input_bundle_rejects_atom_count_mismatch(tmp_path: Path):
    """Identity mapping should fail clearly if RDKit and PDB counts differ."""
    pdb_path = tmp_path / "two_atoms.pdb"
    pdb_path.write_text(
        _pdb_atom(1, "C1", "MOL", "C", 1, 0.0, 0.0, 0.0, element="C")
        + _pdb_atom(2, "H1", "MOL", "C", 1, 0.0, 1.0, 0.0, element="H")
        + "END\n",
        encoding="utf-8",
    )
    mol = Chem.MolFromSmiles("C")
    assert mol is not None

    with pytest.raises(ValueError, match="atom-count mismatch"):
        build_rdkit_input_bundle(
            mol,
            pdb_atoms=parse_pdb_atom_records(pdb_path),
            source_path=pdb_path,
        )


def _pdb_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x_coord: float,
    y_coord: float,
    z_coord: float,
    *,
    element: str = "C",
    record: str = "ATOM",
) -> str:
    """Format one PDB atom line for RDKit input tests."""
    return (
        f"{record:<6}{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    {x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )
