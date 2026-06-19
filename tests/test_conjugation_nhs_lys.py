"""Tests for NHS-Lys conjugation graph edit primitives."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

import polyzymd.builders.conjugation.linkers as linkers_module
from polyzymd.builders.conjugation.builder import CovalentModificationBuilder
from polyzymd.builders.conjugation.diagnostics import DiagnosticCode
from polyzymd.builders.conjugation.exceptions import (
    ConjugationError,
    ConjugationNotImplementedError,
)
from polyzymd.builders.conjugation.linkers import NhsLysModifierLinker, resolve_modifier_nhs_atoms
from polyzymd.builders.conjugation.nhs_lys import (
    detect_nhs_reactive_group,
    execute_nhs_lys_amide_rdkit_graph_edit,
    extract_lysine_reactive_site,
    plan_nhs_lys_amide,
)
from polyzymd.builders.conjugation.polymer_fragment import PolymerFragmentAtom
from polyzymd.builders.conjugation.reactions.nhs_lys import NhsLysReaction
from polyzymd.builders.conjugation.sites import AttachmentSite
from polyzymd.config.schema import ConjugationConfig

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


def _nhs_lys_config(*, mechanism: str = "nhs_lys_amide", enabled_count: int = 1):
    """Build a construct-mode conjugation config for builder execution tests.

    Parameters
    ----------
    mechanism : str, optional
        Mechanism name assigned to the first attachment, by default
        ``"nhs_lys_amide"``.
    enabled_count : int, optional
        Number of enabled attachments to create, by default ``1``.

    Returns
    -------
    ConjugationConfig
        Minimal enabled construct-mode conjugation config.
    """
    attachments = []
    for index in range(enabled_count):
        attachments.append(
            {
                "name": f"lys23-nhs-{index + 1}",
                "site": {
                    "chain_id": "A",
                    "residue_name": "LYS",
                    "residue_number": 23,
                    "atom_name": "NZ",
                },
                "moiety": {
                    "name": "NHS-linker",
                    "role": "moiety",
                    "smiles": "CC(=O)ON1C(=O)CCC1=O",
                },
                "mechanism": {"name": mechanism},
            }
        )
    return ConjugationConfig(enabled=True, mode="construct", attachments=attachments)


def _explicit_execution_context(protein, moiety, *, explicit_nhs_group=True):
    """Build a minimal explicit RDKit execution context.

    Parameters
    ----------
    protein : rdkit.Chem.Mol
        Protein fragment molecule.
    moiety : rdkit.Chem.Mol
        NHS moiety molecule.
    explicit_nhs_group : bool, optional
        Include explicit NHS group indices instead of using autodetection, by
        default ``True``.

    Returns
    -------
    dict
        Builder context payload.
    """
    site_atom_index = next(atom.GetIdx() for atom in protein.GetAtoms() if atom.GetSymbol() == "N")
    site_hydrogen_index = next(
        neighbor.GetIdx()
        for neighbor in protein.GetAtomWithIdx(site_atom_index).GetNeighbors()
        if neighbor.GetSymbol() == "H"
    )
    payload = {
        "protein_mol": protein,
        "moiety_mol": moiety,
        "explicit_site_atom_index": site_atom_index,
        "explicit_site_hydrogen_indices": (site_hydrogen_index,),
    }
    if explicit_nhs_group:
        payload["explicit_nhs_group"] = detect_nhs_reactive_group(moiety).model_dump(mode="json")
    return {"conjugation_rdkit_execution": payload}


def _protein_topology_metadata_from_mol(protein):
    """Build lightweight lysine-like atom and bond metadata from an RDKit fragment.

    Parameters
    ----------
    protein : rdkit.Chem.Mol
        Protein fragment molecule with explicit hydrogens.

    Returns
    -------
    tuple[list[dict], list[tuple[int, int]]]
        Atom records and bond pairs accepted by the lysine site extractor.
    """
    nitrogen_index = next(atom.GetIdx() for atom in protein.GetAtoms() if atom.GetSymbol() == "N")
    carbon_index = next(
        neighbor.GetIdx()
        for neighbor in protein.GetAtomWithIdx(nitrogen_index).GetNeighbors()
        if neighbor.GetSymbol() == "C"
    )
    atoms = []
    for atom in protein.GetAtoms():
        index = atom.GetIdx()
        if index == nitrogen_index:
            atom_name = "NZ"
        elif index == carbon_index:
            atom_name = "CE"
        elif atom.GetSymbol() == "H":
            atom_name = f"HZ{index}"
        else:
            atom_name = atom.GetSymbol()
        atoms.append(
            {
                "index": index,
                "atom_name": atom_name,
                "atomic_number": atom.GetAtomicNum(),
                "metadata": {
                    "chain_id": "A",
                    "residue_name": "LYS",
                    "residue_number": 23,
                },
            }
        )
    bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in protein.GetBonds()]
    return atoms, bonds


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


def test_nhs_lys_reaction_template_creates_current_linker_defaults():
    """The reaction template should produce current NHS-Lys linker defaults."""
    linker = NhsLysReaction.create_linker(target_residue_number=23)

    assert isinstance(linker, NhsLysModifierLinker)
    assert linker.target_chain == "A"
    assert linker.target_residue_name == "LYS"
    assert linker.target_residue_number == 23
    assert linker.target_atom_name == "NZ"
    assert linker.lysine_target_resname == "LYX"
    assert linker.modifier_target_resname == "NHX"
    assert linker.max_nz_hydrogens_to_remove == 2


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


def test_builder_executes_one_nhs_lys_graph_edit_with_explicit_indices():
    """Builder should run one explicit NHS-Lys graph edit and keep topology unchanged."""
    topology = "unchanged-topology"
    protein = _mol_from_smiles("CN")
    moiety = _mol_from_smiles("CC(=O)ON1C(=O)CCC1=O")
    builder = CovalentModificationBuilder(_nhs_lys_config())

    result = builder.build(
        topology,
        context=_explicit_execution_context(protein, moiety, explicit_nhs_group=True),
    )

    assert result.topology == topology
    assert len(result.graph_edit_results) == 1
    assert len(result.graph_edit_summaries) == 1
    graph_result = result.graph_edit_results[0].graph_edit_result
    assert graph_result.product_mol.GetBondBetweenAtoms(
        graph_result.added_bond.begin_atom_index,
        graph_result.added_bond.end_atom_index,
    )
    summary = result.graph_edit_summaries[0]
    assert summary.mechanism == "nhs_lys_amide"
    assert summary.topology_unchanged is True
    assert summary.removed_atoms_count >= 2
    assert any(
        event.code == DiagnosticCode.GRAPH_EDIT_EXECUTION
        for event in result.diagnostics.diagnostics
    )


def test_builder_executes_with_topology_site_extraction_and_nhs_autodetection():
    """Builder should resolve the lysine site from metadata and autodetect NHS atoms."""
    protein = _mol_from_smiles("CN")
    moiety = _mol_from_smiles("CC(=O)ON1C(=O)CCC1=O")
    atoms, bonds = _protein_topology_metadata_from_mol(protein)
    builder = CovalentModificationBuilder(_nhs_lys_config())

    result = builder.build(
        "topology",
        context={
            "conjugation_rdkit_execution": {
                "protein_mol": protein,
                "moiety_mol": moiety,
                "protein_topology_atoms": atoms,
                "protein_topology_bonds": bonds,
            }
        },
    )

    summary = result.graph_edit_summaries[0]
    assert summary.site["residue_name"] == "LYS"
    assert summary.removed_protein_atom_indices
    assert summary.removed_moiety_atom_indices
    assert (
        summary.product_atom_count
        == result.graph_edit_results[0].graph_edit_result.product_mol.GetNumAtoms()
    )


def test_builder_without_explicit_context_preserves_not_implemented_behavior():
    """Construct mode without explicit execution context should still fail as deferred."""
    builder = CovalentModificationBuilder(_nhs_lys_config())

    with pytest.raises(ConjugationNotImplementedError, match="graph surgery"):
        builder.build("topology")


def test_builder_rejects_multiple_attachments_with_explicit_context():
    """The controlled RDKit execution path should reject multi-site requests."""
    protein = _mol_from_smiles("CN")
    moiety = _mol_from_smiles("CC(=O)ON1C(=O)CCC1=O")
    builder = CovalentModificationBuilder(_nhs_lys_config(enabled_count=2))

    with pytest.raises(ConjugationError, match="exactly one enabled attachment"):
        builder.build("topology", context=_explicit_execution_context(protein, moiety))


def test_builder_rejects_non_nhs_mechanism_with_explicit_context():
    """The controlled RDKit execution path should reject non-NHS mechanisms."""
    config = ConjugationConfig(
        enabled=True,
        mode="construct",
        attachments=[
            {
                "name": "asn23-glycan",
                "site": {
                    "chain_id": "A",
                    "residue_name": "ASN",
                    "residue_number": 23,
                    "atom_name": "ND2",
                },
                "moiety": {"name": "glycan", "role": "glycan", "smiles": "CO"},
                "mechanism": {"name": "n_glycosidic_asn"},
            }
        ],
    )
    protein = _mol_from_smiles("CN")
    moiety = _mol_from_smiles("CC(=O)ON1C(=O)CCC1=O")
    builder = CovalentModificationBuilder(config)

    with pytest.raises(ConjugationError, match="nhs_lys_amide"):
        builder.build("topology", context=_explicit_execution_context(protein, moiety))


def test_builder_result_serialization_excludes_rdkit_mols_and_keeps_summary():
    """Serialized builder results should include summaries but exclude RDKit objects."""
    protein = _mol_from_smiles("CN")
    moiety = _mol_from_smiles("CC(=O)ON1C(=O)CCC1=O")
    builder = CovalentModificationBuilder(_nhs_lys_config())

    result = builder.build(
        "topology",
        context=_explicit_execution_context(protein, moiety, explicit_nhs_group=True),
    )
    dumped = result.model_dump(mode="json")

    assert "graph_edit_results" not in dumped
    assert dumped["graph_edit_summaries"][0]["mechanism"] == "nhs_lys_amide"
    assert "product_mol" not in str(dumped)


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
