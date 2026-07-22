"""Tests for the N-glycosylation conjugation reaction template."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation._linkage import PdbAtomSelector
from polyzymd.builders.conjugation.polymer import build_smiles_moiety_fragment
from polyzymd.builders.conjugation.reactions import (
    NGlycosylationReaction,
    OGlycosylationReaction,
    ReactionTemplate,
    get_reaction,
    list_reactions,
)
from polyzymd.builders.conjugation.reactions.n_glycosylation import (
    _resolve_bonded_site_hydrogen,
    detect_glycan_anomeric_group,
)
from polyzymd.config.schema import ConjugationAttachmentConfig

pytest.importorskip("rdkit")

BRANCHED_GLCNAC_SMILES = "CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O"
THREE_BRANCH_GLYCAN_SMILES = (
    "CC(=O)N[C@@H]1[C@H]([C@@H]([C@H](O[C@H]1O)CO)O[C@H]2[C@@H]([C@H]"
    "([C@@H]([C@H](O2)CO)O[C@H]3[C@H]([C@H]([C@@H]([C@H](O3)CO[C@@H]4"
    "[C@H]([C@H]([C@@H]([C@H](O4)CO[C@@H]5[C@H]([C@H]([C@@H]([C@H](O5)CO)O)O)O)O)"
    "O[C@@H]6[C@H]([C@H]([C@@H]([C@H](O6)CO)O)O)O[C@@H]7[C@H]([C@H]([C@@H]([C@H](O7)CO)O)O)O)O)O)"
    "O[C@@H]8[C@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O[C@@H]9[C@H]([C@H]([C@@H]([C@H](O9)CO)O)O)"
    "O[C@@H]1[C@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O)O)NC(=O)C)O"
)


def test_registry_exposes_n_glycosylation_template_and_aliases():
    """The built-in reaction registry should expose N-glycosylation aliases."""
    reactions = list_reactions()

    assert "n_glycosylation" in reactions
    assert get_reaction("n_glycosylation") is NGlycosylationReaction
    assert get_reaction("n_glycan_asn") is NGlycosylationReaction
    assert get_reaction("asn_n_glycosylation") is NGlycosylationReaction
    assert issubclass(get_reaction("n_glycosylation"), ReactionTemplate)


def test_registry_exposes_o_glycosylation_template_and_aliases():
    """The built-in registry should expose O-glycosylation aliases."""
    assert get_reaction("o_glycosylation") is OGlycosylationReaction
    assert get_reaction("o_glycan") is OGlycosylationReaction


@pytest.mark.parametrize(
    ("residue_name", "site_atom", "site_hydrogen", "product_residue"),
    (("SER", "OG", "HG", "OLS"), ("THR", "OG1", "HG1", "OLT")),
)
def test_o_glycosylation_resolves_ser_thr_products(
    tmp_path: Path,
    residue_name: str,
    site_atom: str,
    site_hydrogen: str,
    product_residue: str,
):
    """Ser and Thr sites should use GLYCAM-compatible O-linked product labels."""
    protein_path = _hydroxyl_site_pdb(
        tmp_path,
        residue_name=residue_name,
        site_atom=site_atom,
        site_hydrogen=site_hydrogen,
    )
    glycan = _glycan_fragment(tmp_path)
    site = SimpleNamespace(
        chain_id="A",
        residue_name=residue_name,
        residue_number=42,
        atom_name=None,
        insertion_code="",
        atom_serial=None,
        atom_index=None,
    )
    attachment = SimpleNamespace(
        site=site,
        mechanism=SimpleNamespace(product_residues=None, bond=None),
    )
    settings = OGlycosylationReaction.settings_from_attachment(attachment)

    plan = OGlycosylationReaction.resolve_plan(
        protein_path,
        site,
        glycan,
        settings=settings,
    )

    assert plan.protein_product_residue_name == product_residue
    assert plan.contract.protein_endpoint.selector.atom_name == site_atom
    assert plan.pablo_crosslink_requirement.linking_atoms[0] == site_atom
    assert plan.pablo_crosslink_requirement.leaving_atoms[0] == (site_hydrogen,)
    assert plan.target_bond_length_angstrom == pytest.approx(1.43)


def test_o_glycosylation_rejects_non_ser_thr_site():
    """O-glycosylation should require an explicit Ser or Thr site."""
    attachment = SimpleNamespace(
        site=SimpleNamespace(residue_name="ASN"),
        mechanism=SimpleNamespace(product_residues=None, bond=None),
    )

    with pytest.raises(ValueError, match="requires site.residue_name SER or THR"):
        OGlycosylationReaction.settings_from_attachment(attachment)


def test_site_atom_defaults_to_asn_nd2_when_omitted(tmp_path: Path):
    """Site resolution should default omitted atom_name to ASN ND2."""
    fragment = _glycan_fragment(tmp_path)
    site = SimpleNamespace(chain_id="A", residue_name="ASN", residue_number=42)

    contract = NGlycosylationReaction.build_contract(site, fragment)

    assert contract.protein_endpoint.selector.residue_name == "ASN"
    assert contract.protein_endpoint.selector.atom_name == "ND2"
    assert contract.protein_endpoint.product_residue_name == "NLN"
    assert contract.bond.protein_atom_name == "ND2"
    assert contract.bond.target_bond_length_angstrom == pytest.approx(1.45)


def test_attachment_settings_default_target_atom_and_bond_length_when_omitted():
    """Mechanism settings should fill N-glycosylation defaults from sparse configs."""
    attachment = ConjugationAttachmentConfig.model_validate(
        {
            "name": "glycan",
            "site": {"chain_id": "A", "residue_name": "ASN", "residue_number": 42},
            "moiety": {"name": "glcnac", "smiles": "CO", "residue_name": "NAG"},
            "mechanism": {"name": "n_glycosylation"},
        }
    )

    settings = NGlycosylationReaction.settings_from_attachment(attachment)

    assert settings.target_atom_name == "ND2"
    assert settings.target_bond_length_angstrom == pytest.approx(1.45)


def test_detects_anomeric_carbon_and_hydroxyl_leaving_atoms_in_branched_glycan(
    tmp_path: Path,
):
    """The glycan endpoint should come from SMARTS/connectivity, not atom names."""
    fragment = _glycan_fragment(tmp_path)

    group = NGlycosylationReaction.detect_anomeric_group(fragment)

    atoms_by_index = {atom.atom_index: atom for atom in fragment.atoms}
    assert atoms_by_index[group.reactive_carbon_index].element == "C"
    assert atoms_by_index[group.hydroxyl_oxygen_index].element == "O"
    assert group.hydroxyl_oxygen_index in group.leaving_atom_indices
    assert any(atoms_by_index[index].element == "H" for index in group.leaving_atom_indices)
    assert atoms_by_index[group.reactive_carbon_index].atom_name not in {
        atoms_by_index[index].atom_name for index in group.leaving_atom_indices
    }


def test_detects_validated_three_branch_glycan_checkpoint():
    """The fast E2E glycan must be the validated three-branch input, not mono-NAG."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    heavy_mol = Chem.MolFromSmiles(THREE_BRANCH_GLYCAN_SMILES)
    assert heavy_mol is not None
    mol = Chem.AddHs(heavy_mol)

    group = detect_glycan_anomeric_group(mol)

    assert heavy_mol.GetNumAtoms() == 117
    assert mol.GetNumAtoms() == 225
    assert rdMolDescriptors.CalcMolFormula(mol) == "C64H108N2O51"
    assert Chem.GetFormalCharge(mol) == 0
    assert mol.GetRingInfo().NumRings() == 10
    assert group.reactive_carbon_index == 9
    assert group.hydroxyl_oxygen_index == 10
    assert group.ring_oxygen_index == 8
    assert group.leaving_atom_indices == (10, 126)


def test_resolve_plan_builds_nln_to_user_glycan_residue_linkage(tmp_path: Path):
    """Resolved plans should carry NLN, the user glycan residue, and exact leaving atoms."""
    protein_path = _asn_pdb(tmp_path)
    fragment = _glycan_fragment(tmp_path, residue_name="NAG")
    site = {"chain_id": "A", "residue_name": "ASN", "residue_number": 42}

    plan = NGlycosylationReaction.resolve_plan(protein_path, site, fragment)

    assert plan.contract.mechanism_name == "n_glycosylation"
    assert plan.protein_link_atom.atom_name == "ND2"
    assert plan.protein_product_residue_name == "NLN"
    assert plan.modifier_product_residue_name == "NAG"
    assert tuple(atom.atom_name for atom in plan.protein_leaving_atoms) == ("HD21",)
    assert plan.pablo_crosslink_requirement.residues == ("NLN", "NAG")
    assert plan.pablo_crosslink_requirement.linking_atoms == (
        "ND2",
        plan.modifier_link_atom.atom_name,
    )
    assert plan.pablo_crosslink_requirement.leaving_atoms[0] == ("HD21",)
    assert len(plan.modifier_leaving_atoms) >= 2


def test_asn_nd2_hydrogen_resolution_selects_hd21_from_canonical_pair(
    tmp_path: Path,
) -> None:
    """Canonical HD21/HD22 candidates should remove HD21 by name, not file order."""
    protein_path = _asn_pdb(tmp_path, hydrogens=("HD22", "HD21"))

    hydrogen = NGlycosylationReaction.resolve_protein_leaving_atom(protein_path, _asn_selector())

    assert hydrogen.atom_name == "HD21"


def test_asn_nd2_hydrogen_resolution_returns_single_candidate_unchanged(
    tmp_path: Path,
) -> None:
    """A single Pablo-template ND2 hydrogen should be accepted as supplied."""
    protein_path = _asn_pdb(tmp_path, hydrogens=("HD22",))

    hydrogen = NGlycosylationReaction.resolve_protein_leaving_atom(protein_path, _asn_selector())

    assert hydrogen.atom_name == "HD22"


def test_asn_nd2_hydrogen_resolution_rejects_ambiguous_noncanonical_candidates(
    tmp_path: Path,
) -> None:
    """Multiple noncanonical template candidates should fail without geometry fallback."""
    protein_path = _asn_pdb(tmp_path, hydrogens=("HN1", "HN2"))
    residue_library = _fake_asn_nd2_hydrogen_library(("HN1", "HN2"))

    with pytest.raises(ValueError, match="unambiguous preferred HD21"):
        _resolve_bonded_site_hydrogen(
            protein_path,
            _asn_selector(),
            mechanism_name="n_glycosylation",
            preferred_name="HD21",
            residue_library=residue_library,
        )


def test_rejects_non_asn_target_residue(tmp_path: Path):
    """N-glycosylation should fail clearly when the selected residue is not ASN."""
    fragment = _glycan_fragment(tmp_path)
    site = {"chain_id": "A", "residue_name": "LYS", "residue_number": 42}

    with pytest.raises(ValueError, match="target residue must be ASN"):
        NGlycosylationReaction.build_contract(site, fragment)


def test_rejects_moiety_without_anomeric_motif(tmp_path: Path):
    """Non-glycan moieties should fail with an anomeric-motif diagnostic."""
    fragment = build_smiles_moiety_fragment("CCO", "ETH", output_dir=tmp_path, random_seed=5)
    site = {"chain_id": "A", "residue_name": "ASN", "residue_number": 42}

    with pytest.raises(ValueError, match="No glycan anomeric motif"):
        NGlycosylationReaction.build_contract(site, fragment)


def _glycan_fragment(tmp_path: Path, *, residue_name: str = "NAG"):
    """Build the branched glycan-like one-residue test fragment."""
    return build_smiles_moiety_fragment(
        BRANCHED_GLCNAC_SMILES,
        residue_name,
        name="glcnac",
        output_dir=tmp_path,
        random_seed=17,
    )


def _asn_selector() -> PdbAtomSelector:
    """Return the default ASN ND2 selector for local fixtures."""
    return PdbAtomSelector(
        chain_id="A",
        residue_name="ASN",
        residue_number=42,
        atom_name="ND2",
        insertion_code="",
    )


def _fake_asn_nd2_hydrogen_library(
    hydrogen_names: tuple[str, ...],
) -> dict[str, tuple[object, ...]]:
    """Return a Pablo-like residue library with injectable ASN ND2 hydrogens."""
    atoms = [
        SimpleNamespace(name="ND2", symbol="N"),
        *(SimpleNamespace(name=name, symbol="H") for name in hydrogen_names),
    ]
    bonds = [SimpleNamespace(atom1="ND2", atom2=name) for name in hydrogen_names]
    return {"ASN": (SimpleNamespace(atoms=atoms, bonds=bonds),)}


def _asn_pdb(tmp_path: Path, *, hydrogens: tuple[str, ...] = ("HD21",)) -> Path:
    """Create a small Asn residue with one explicit ND2 leaving hydrogen."""
    path = tmp_path / "asn.pdb"
    lines = [
        _pdb_atom(1, "N", "ASN", "A", 42, 0.0, 0.0, 0.0, element="N"),
        _pdb_atom(2, "CA", "ASN", "A", 42, 1.0, 0.0, 0.0),
        _pdb_atom(3, "CB", "ASN", "A", 42, 1.5, 0.7, 0.0),
        _pdb_atom(4, "CG", "ASN", "A", 42, 2.0, 0.0, 0.0),
        _pdb_atom(5, "OD1", "ASN", "A", 42, 2.0, -1.2, 0.0, element="O"),
        _pdb_atom(6, "ND2", "ASN", "A", 42, 2.0, 1.2, 0.0, element="N"),
    ]
    for offset, hydrogen_name in enumerate(hydrogens, start=7):
        lines.append(
            _pdb_atom(
                offset,
                hydrogen_name,
                "ASN",
                "A",
                42,
                2.0,
                2.1 + offset / 10,
                0.0,
                element="H",
            )
        )
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def _hydroxyl_site_pdb(
    tmp_path: Path,
    *,
    residue_name: str,
    site_atom: str,
    site_hydrogen: str,
) -> Path:
    """Create a small Ser or Thr residue with its explicit hydroxyl hydrogen."""
    path = tmp_path / f"{residue_name.lower()}.pdb"
    lines = [
        _pdb_atom(1, "N", residue_name, "A", 42, 0.0, 0.0, 0.0, element="N"),
        _pdb_atom(2, "CA", residue_name, "A", 42, 1.0, 0.0, 0.0),
        _pdb_atom(3, "CB", residue_name, "A", 42, 1.5, 0.7, 0.0),
        _pdb_atom(4, site_atom, residue_name, "A", 42, 2.0, 1.2, 0.0, element="O"),
        _pdb_atom(5, site_hydrogen, residue_name, "A", 42, 2.0, 2.0, 0.0, element="H"),
    ]
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
) -> str:
    """Format one PDB atom line for N-glycosylation tests."""
    return (
        f"ATOM  {serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    {x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )
