"""Tests for NHS-Lys conjugation graph edit primitives."""

from __future__ import annotations

import pytest

from polyzymd.builders.conjugation.nhs_lys import (
    detect_nhs_reactive_group,
    execute_nhs_lys_amide_rdkit_graph_edit,
    extract_lysine_reactive_site,
    plan_nhs_lys_amide,
)
from polyzymd.builders.conjugation.sites import AttachmentSite

pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402


def _mol_from_smiles(smiles: str):
    """Build an explicit-hydrogen RDKit molecule for tests.

    Parameters
    ----------
    smiles : str
        Input SMILES string.

    Returns
    -------
    rdkit.Chem.Mol
        RDKit molecule with explicit hydrogens.
    """
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    return Chem.AddHs(mol)


def test_detect_nhs_reactive_group_on_nhs_ester_like_molecule():
    """NHS detection should find the acyl carbon and traversed leaving group."""
    mol = _mol_from_smiles("CC(=O)ON1C(=O)CCC1=O")

    group = detect_nhs_reactive_group(mol)

    assert mol.GetAtomWithIdx(group.reactive_carbon_index).GetSymbol() == "C"
    assert mol.GetAtomWithIdx(group.bridging_oxygen_index).GetSymbol() == "O"
    assert group.bridging_oxygen_index in group.leaving_group_atom_indices
    assert group.reactive_carbon_index in group.retained_atom_indices
    assert group.reactive_carbon_index not in group.leaving_group_atom_indices


def test_detect_nhs_reactive_group_rejects_missing_or_ambiguous_groups():
    """NHS detection should fail clearly when assignments are absent or ambiguous."""
    methyl_acetate = _mol_from_smiles("CC(=O)OC")
    with pytest.raises(ValueError, match="No unambiguous NHS ester reactive group"):
        detect_nhs_reactive_group(methyl_acetate)

    two_groups = _mol_from_smiles("CC(=O)ON1C(=O)CCC1=O.CC(=O)ON1C(=O)CCC1=O")
    with pytest.raises(ValueError, match="Ambiguous NHS ester reactive group"):
        detect_nhs_reactive_group(two_groups)


def test_execute_nhs_lys_graph_edit_adds_bond_and_removes_leaving_group():
    """RDKit graph surgery should remove supplied atoms and add one NZ-C bond."""
    protein = _mol_from_smiles("CN")
    moiety = _mol_from_smiles("CC(=O)ON1C(=O)CCC1=O")
    AllChem.EmbedMolecule(protein, randomSeed=7)
    AllChem.EmbedMolecule(moiety, randomSeed=11)
    group = detect_nhs_reactive_group(moiety)
    site_atom_index = next(atom.GetIdx() for atom in protein.GetAtoms() if atom.GetSymbol() == "N")
    site_hydrogen_index = next(
        neighbor.GetIdx()
        for neighbor in protein.GetAtomWithIdx(site_atom_index).GetNeighbors()
        if neighbor.GetSymbol() == "H"
    )

    result = execute_nhs_lys_amide_rdkit_graph_edit(
        protein_mol=protein,
        moiety_mol=moiety,
        site_atom_index=site_atom_index,
        reactive_carbon_index=group.reactive_carbon_index,
        leaving_atom_indices=group.leaving_group_atom_indices,
        site_hydrogen_indices=(site_hydrogen_index,),
    )

    product = result.product_mol
    expected_atoms = (
        protein.GetNumAtoms() - 1 + moiety.GetNumAtoms() - len(group.leaving_group_atom_indices)
    )
    assert product.GetNumAtoms() == expected_atoms
    assert product.GetBondBetweenAtoms(
        result.added_bond.begin_atom_index,
        result.added_bond.end_atom_index,
    )
    assert result.removed_protein_atom_indices == (site_hydrogen_index,)
    assert result.removed_moiety_atom_indices == group.leaving_group_atom_indices


def test_execute_nhs_lys_graph_edit_strips_stale_conformers():
    """Graph edit results should strip input conformers and warn about coordinates."""
    protein = _mol_from_smiles("CN")
    moiety = _mol_from_smiles("CC(=O)ON1C(=O)CCC1=O")
    AllChem.EmbedMolecule(protein, randomSeed=17)
    AllChem.EmbedMolecule(moiety, randomSeed=19)
    group = detect_nhs_reactive_group(moiety)
    site_atom_index = next(atom.GetIdx() for atom in protein.GetAtoms() if atom.GetSymbol() == "N")
    site_hydrogen_index = next(
        neighbor.GetIdx()
        for neighbor in protein.GetAtomWithIdx(site_atom_index).GetNeighbors()
        if neighbor.GetSymbol() == "H"
    )

    result = execute_nhs_lys_amide_rdkit_graph_edit(
        protein_mol=protein,
        moiety_mol=moiety,
        site_atom_index=site_atom_index,
        reactive_carbon_index=group.reactive_carbon_index,
        leaving_atom_indices=group.leaving_group_atom_indices,
        site_hydrogen_indices=(site_hydrogen_index,),
    )

    assert result.product_mol.GetNumConformers() == 0
    assert any("conformers" in warning for warning in result.warnings)


def test_extract_lysine_site_from_mocked_topology_atom_metadata():
    """Explicit site matching should resolve NZ, CE, and NZ-bound hydrogens."""
    atoms = [
        {
            "index": 0,
            "name": "CE",
            "atom_name": "CE",
            "atomic_number": 6,
            "metadata": {"chain_id": "A", "residue_name": "LYS", "residue_number": "23"},
        },
        {
            "index": 1,
            "name": "NZ",
            "atom_name": "NZ",
            "atomic_number": 7,
            "metadata": {"chain_id": "A", "residue_name": "LYS", "residue_number": "23"},
        },
        {
            "index": 2,
            "name": "HZ1",
            "atom_name": "HZ1",
            "atomic_number": 1,
            "metadata": {"chain_id": "A", "residue_name": "LYS", "residue_number": "23"},
        },
        {
            "index": 3,
            "name": "HZ2",
            "atom_name": "HZ2",
            "atomic_number": 1,
            "metadata": {"chain_id": "A", "residue_name": "LYS", "residue_number": "23"},
        },
    ]
    site = AttachmentSite(
        chain_id="A",
        residue_name="LYS",
        residue_number=23,
        atom_name="NZ",
    )

    resolved = extract_lysine_reactive_site(
        site,
        atoms,
        bonds=[(0, 1), (1, 2), (1, 3)],
        positions={1: (1.0, 2.0, 3.0), 0: (0.0, 2.0, 3.0)},
    )
    plan = plan_nhs_lys_amide(
        resolved,
        detect_nhs_reactive_group(_mol_from_smiles("CC(=O)ON1C(=O)CCC1=O")),
        site_hydrogen_indices_to_remove=(3,),
    )

    assert resolved.nz_atom_index == 1
    assert resolved.ce_atom_index == 0
    assert resolved.nz_hydrogen_indices == (2, 3)
    assert resolved.nz_position == (1.0, 2.0, 3.0)
    assert plan.mechanism == "nhs_lys_amide"
    assert plan.add_bond[0] == 1
    assert plan.site_hydrogen_indices_to_remove == (3,)
