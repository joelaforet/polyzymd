"""Tests for NHS-Lys conjugation graph edit primitives."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

import polyzymd.builders.conjugation._linkage as linkers_module
from polyzymd.builders.conjugation._linkage import NhsLysModifierLinker, resolve_modifier_nhs_atoms
from polyzymd.builders.conjugation.polymer import (
    GeneratedPolymerFragment,
    PolymerFragmentAtom,
    PreparedFragment,
)
from polyzymd.builders.conjugation.reactions.nhs_lys import (
    NhsLysAttachmentSite,
    NhsLysReaction,
    detect_nhs_reactive_group,
    execute_nhs_lys_amide_rdkit_graph_edit,
    extract_lysine_reactive_site,
    plan_nhs_lys_amide,
)

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
    site = NhsLysAttachmentSite(
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


def test_nhs_lys_linker_resolves_poc_like_hydrogen_names_by_rdkit(tmp_path):
    """NHS-Lys convenience should resolve NZ hydrogens by chemistry, not names."""
    protein_path = _lysine_pdb(
        tmp_path,
        hydrogens=(("H10", 2.0, 0.7, 0.0), ("H11", 2.0, -0.7, 0.0), ("H13", 2.0, 0.0, 0.7)),
    )
    linker = NhsLysModifierLinker(target_residue_number=23)

    site = linker.resolve_site(protein_path)
    attachment = linker.attachment(protein_path)

    assert site.atom.atom_name == "NZ"
    assert tuple(atom.atom_name for atom in site.removable_hydrogens) == ("H11", "H13")
    assert attachment.nz_hydrogen_atom_names_to_remove == ("H11", "H13")
    assert attachment.nz_hydrogen_atom_serials_to_remove == (7, 8)
    assert any("Protonated lysine" in warning for warning in site.warnings)


def test_nhs_lys_reaction_resolves_generic_attachment_plan(tmp_path):
    """The NHS-Lys reaction should own generic plan resolution defaults."""
    protein_path = _lysine_pdb(
        tmp_path,
        hydrogens=(("HZ1", 2.0, 0.7, 0.0), ("HZ2", 2.0, -0.7, 0.0), ("HZ3", 2.0, 0.0, 0.7)),
    )
    site = SimpleNamespace(chain_id="A", residue_name="LYS", residue_number=23, atom_name="NZ")
    raw_modifier = _modifier_with_explicit_selectors(rdkit_mol=None)
    modifier = GeneratedPolymerFragment(
        atoms=raw_modifier.atoms,
        reactive_atom_index=raw_modifier.reactive_atom_index,
        leaving_atom_indices=raw_modifier.leaving_atom_indices,
        name=raw_modifier.name,
    )

    plan = NhsLysReaction().resolve_attachment(protein_path, site, modifier)

    assert plan.protein_product_residue_name == "LYX"
    assert plan.modifier_product_residue_name == "NHX"
    assert plan.pablo_crosslink_requirement.residues == ("LYX", "NHX")
    assert plan.pablo_crosslink_requirement.linking_atoms == ("NZ", "C1")
    assert plan.target_bond_length_angstrom == pytest.approx(1.33)


def test_nhs_lys_reaction_rejects_mismatched_prepared_fragment(tmp_path):
    """Reaction endpoints and the canonical prepared structure must share one atom graph."""
    protein_path = _lysine_pdb(
        tmp_path,
        hydrogens=(("HZ1", 2.0, 0.7, 0.0), ("HZ2", 2.0, -0.7, 0.0), ("HZ3", 2.0, 0.0, 0.7)),
    )
    site = SimpleNamespace(chain_id="A", residue_name="LYS", residue_number=23, atom_name="NZ")
    raw_modifier = _modifier_with_explicit_selectors(rdkit_mol=None)
    modifier = GeneratedPolymerFragment(
        atoms=raw_modifier.atoms,
        reactive_atom_index=raw_modifier.reactive_atom_index,
        leaving_atom_indices=raw_modifier.leaving_atom_indices,
        name=raw_modifier.name,
    )
    mismatched = modifier.model_copy(
        update={
            "atoms": (
                modifier.atoms[0].model_copy(update={"atom_name": "BAD"}),
                *modifier.atoms[1:],
            )
        }
    )
    prepared = PreparedFragment.from_generated_fragment(
        mismatched,
        source_identity="mismatched",
        source_kind="polymer",
    )

    with pytest.raises(ValueError, match="does not represent the reaction fragment"):
        NhsLysReaction().resolve_attachment(
            protein_path,
            site,
            modifier,
            prepared_fragment=prepared,
        )


def test_nhs_lys_linker_resolves_canonical_hz_names_by_rdkit(tmp_path):
    """Canonical HZ names should work because they are N-bound hydrogens."""
    protein_path = _lysine_pdb(
        tmp_path,
        hydrogens=(("HZ1", 2.0, 0.7, 0.0), ("HZ2", 2.0, -0.7, 0.0), ("HZ3", 2.0, 0.0, 0.7)),
    )
    linker = NhsLysModifierLinker(target_residue_number=23)

    site = linker.resolve_site(protein_path)

    assert tuple(atom.atom_name for atom in site.removable_hydrogens) == ("HZ2", "HZ3")


def test_nhs_lys_linker_fails_clearly_without_required_hydrogens(tmp_path):
    """Missing NZ hydrogens should fail until normalization policy is implemented."""
    protein_path = _lysine_pdb(tmp_path, hydrogens=())
    linker = NhsLysModifierLinker(target_residue_number=23)

    with pytest.raises(ValueError, match="Automatic hydrogen addition"):
        linker.resolve_site(protein_path)


def test_resolve_modifier_nhs_atoms_warns_when_rdkit_detection_fallback_succeeds(
    monkeypatch, caplog
):
    """Explicit modifier selectors should warn when RDKit NHS detection fails."""
    modifier = _modifier_with_explicit_selectors(rdkit_mol=object())

    def fail_detection(mol):
        """Raise a synthetic RDKit detection failure."""
        raise ValueError(f"no NHS group in {type(mol).__name__}")

    monkeypatch.setattr(linkers_module, "detect_nhs_reactive_group", fail_detection)

    with caplog.at_level("WARNING"):
        reactive_atom, leaving_atoms = resolve_modifier_nhs_atoms(modifier)

    assert reactive_atom.atom_name == "C1"
    assert tuple(atom.atom_name for atom in leaving_atoms) == ("O1",)
    assert "RDKit NHS detection failed" in caplog.text
    assert "explicit generated-fragment fallback reactive atom C1" in caplog.text
    assert "no NHS group" in caplog.text


def test_resolve_modifier_nhs_atoms_reports_rdkit_and_fallback_failures(monkeypatch):
    """Modifier NHS resolution should report both RDKit and fallback failures."""
    modifier = _modifier_with_explicit_selectors(rdkit_mol=object(), reactive_atom_index=None)

    def fail_detection(mol):
        """Raise a synthetic RDKit detection failure."""
        raise ValueError(f"no NHS group in {type(mol).__name__}")

    monkeypatch.setattr(linkers_module, "detect_nhs_reactive_group", fail_detection)

    with pytest.raises(ValueError, match="RDKit NHS detection failed.*Fallback failure"):
        resolve_modifier_nhs_atoms(modifier)


def _modifier_with_explicit_selectors(*, rdkit_mol, reactive_atom_index=0):
    """Build a modifier-like object with explicit generated-fragment selectors."""
    atoms = (
        PolymerFragmentAtom(
            atom_index=0,
            serial=1,
            atom_name="C1",
            residue_name="NHX",
            residue_number=1,
            chain_id="C",
            x=0.0,
            y=0.0,
            z=0.0,
            element="C",
        ),
        PolymerFragmentAtom(
            atom_index=1,
            serial=2,
            atom_name="O1",
            residue_name="NHX",
            residue_number=1,
            chain_id="C",
            x=1.0,
            y=0.0,
            z=0.0,
            element="O",
        ),
    )
    return SimpleNamespace(
        atoms=atoms,
        name="test-modifier",
        rdkit_mol=rdkit_mol,
        reactive_atom_serial=None,
        reactive_atom_index=reactive_atom_index,
        reactive_atom_name=None,
        leaving_atom_serials=(),
        leaving_atom_indices=(1,),
        leaving_atom_names=(),
    )


def _lysine_pdb(tmp_path, *, hydrogens):
    """Create a lysine PDB with configurable NZ hydrogen names."""
    path = tmp_path / "lysine.pdb"
    lines = [
        _pdb_atom(1, "N", "LYS", "A", 23, 0.0, 0.0, 0.0, element="N"),
        _pdb_atom(2, "CA", "LYS", "A", 23, 1.0, 0.0, 0.0),
        _pdb_atom(3, "CD", "LYS", "A", 23, 1.3, 0.0, 0.0),
        _pdb_atom(4, "CE", "LYS", "A", 23, 1.6, 0.0, 0.0),
        _pdb_atom(5, "NZ", "LYS", "A", 23, 2.0, 0.0, 0.0, element="N"),
    ]
    for offset, (name, x_coord, y_coord, z_coord) in enumerate(hydrogens, start=6):
        lines.append(
            _pdb_atom(
                offset,
                name,
                "LYS",
                "A",
                23,
                x_coord,
                y_coord,
                z_coord,
                element="H",
            )
        )
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


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
    """Format one PDB atom line for NHS-Lys tests."""
    return (
        f"{record:<6}{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    {x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )
