"""Tests for SMILES-defined one-residue moiety fragment generation."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.builders.conjugation.polymer import (
    GeneratedMoietyFragment,
    build_smiles_moiety_fragment,
)

pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402


def test_builds_neutral_smiles_moiety_fragment(tmp_path: Path):
    """Neutral SMILES should produce one residue with explicit hydrogens and sidecars."""
    fragment = build_smiles_moiety_fragment(
        "CCO",
        "MOE",
        name="ethanol",
        output_dir=tmp_path,
        random_seed=17,
    )

    assert isinstance(fragment, GeneratedMoietyFragment)
    assert fragment.name == "ethanol"
    assert fragment.residue_name == "MOE"
    assert len(fragment.residues) == 1
    assert {atom.residue_name for atom in fragment.atoms} == {"MOE"}
    assert {atom.residue_number for atom in fragment.atoms} == {1}
    assert len(fragment.atoms) == 9
    assert fragment.pdb_path == tmp_path / "ethanol.pdb"
    assert fragment.sdf_path == tmp_path / "ethanol.sdf"
    assert fragment.pdb_path.exists()
    assert fragment.sdf_path.exists()


def test_builds_charged_smiles_moiety_fragment():
    """Formal charges from charged SMILES should be preserved on fragment atoms."""
    fragment = build_smiles_moiety_fragment("[NH4+]", "AMM", random_seed=7)

    charges = {atom.atom_name: atom.formal_charge for atom in fragment.atoms}
    assert sum(charge or 0 for charge in charges.values()) == 1
    assert next(atom for atom in fragment.atoms if atom.element == "N").formal_charge == 1
    assert all(order == 1.0 for *_atoms, order in fragment.bond_orders)


@pytest.mark.parametrize("residue_name", ["MO", "MOEY", "M-E", "M E"])
def test_rejects_invalid_residue_names(residue_name: str):
    """Moiety residue names must be exactly three PDB-safe characters."""
    with pytest.raises(ValueError, match="residue names must be exactly three"):
        build_smiles_moiety_fragment("CO", residue_name)


def test_sidecar_atom_counts_align(tmp_path: Path):
    """PDB and SDF sidecars should carry the same atom count as the fragment."""
    fragment = build_smiles_moiety_fragment(
        "C(=O)[O-]",
        "ACE",
        output_dir=tmp_path,
        random_seed=11,
    )
    assert fragment.pdb_path is not None
    assert fragment.sdf_path is not None

    pdb_atom_count = sum(
        1
        for line in fragment.pdb_path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM", "HETATM"))
    )
    supplier = Chem.SDMolSupplier(str(fragment.sdf_path), removeHs=False, sanitize=False)
    sdf_mols = [mol for mol in supplier if mol is not None]

    assert len(sdf_mols) == 1
    assert pdb_atom_count == len(fragment.atoms)
    assert sdf_mols[0].GetNumAtoms() == len(fragment.atoms)


def test_preserves_bond_orders_and_formal_charges():
    """Fragments should retain RDKit bond orders and non-zero formal charges."""
    fragment = build_smiles_moiety_fragment("C(=O)[O-]", "ACE", random_seed=23)

    assert any(order == 2.0 for *_atoms, order in fragment.bond_orders)
    assert any(atom.formal_charge == -1 for atom in fragment.atoms)
    assert set(fragment.bonds) == {(left, right) for left, right, _order in fragment.bond_orders}


def test_atom_naming_is_deterministic_with_seed():
    """Atom names should be deterministic and unique for repeat generation."""
    first = build_smiles_moiety_fragment("CC(=O)O", "ACY", random_seed=101)
    second = build_smiles_moiety_fragment("CC(=O)O", "ACY", random_seed=101)

    first_names = [atom.atom_name for atom in first.atoms]
    assert first_names == [atom.atom_name for atom in second.atoms]
    assert len(first_names) == len(set(first_names))
    assert all(1 <= len(atom_name) <= 4 for atom_name in first_names)
    assert any(atom_name.startswith("C") and atom_name != "C" for atom_name in first_names)
