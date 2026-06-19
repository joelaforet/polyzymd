"""Tests for the N-glycosylation conjugation reaction template."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.polymer import build_smiles_moiety_fragment
from polyzymd.builders.conjugation.reactions import (
    NGlycosylationReaction,
    ReactionTemplate,
    get_reaction,
    list_reactions,
)

pytest.importorskip("rdkit")

BRANCHED_GLCNAC_SMILES = "CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O"


def test_registry_exposes_n_glycosylation_template_and_aliases():
    """The built-in reaction registry should expose N-glycosylation aliases."""
    reactions = list_reactions()

    assert "n_glycosylation" in reactions
    assert get_reaction("n_glycosylation") is NGlycosylationReaction
    assert get_reaction("n_glycan_asn") is NGlycosylationReaction
    assert get_reaction("asn_n_glycosylation") is NGlycosylationReaction
    assert issubclass(get_reaction("n_glycosylation"), ReactionTemplate)


def test_site_atom_defaults_to_asn_nd2_when_omitted(tmp_path: Path):
    """Site resolution should default omitted atom_name to ASN ND2."""
    fragment = _glycan_fragment(tmp_path)
    site = SimpleNamespace(chain_id="A", residue_name="ASN", residue_number=42)

    contract = NGlycosylationReaction.build_contract(site, fragment)

    assert contract.protein_endpoint.selector.residue_name == "ASN"
    assert contract.protein_endpoint.selector.atom_name == "ND2"
    assert contract.protein_endpoint.product_residue_name == "ASX"
    assert contract.bond.protein_atom_name == "ND2"
    assert contract.bond.target_bond_length_angstrom == pytest.approx(1.45)


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


def test_resolve_plan_builds_asx_to_user_glycan_residue_linkage(tmp_path: Path):
    """Resolved plans should carry ASX, the user glycan residue, and exact leaving atoms."""
    protein_path = _asn_pdb(tmp_path)
    fragment = _glycan_fragment(tmp_path, residue_name="NAG")
    site = {"chain_id": "A", "residue_name": "ASN", "residue_number": 42}

    plan = NGlycosylationReaction.resolve_plan(protein_path, site, fragment)

    assert plan.contract.mechanism_name == "n_glycosylation"
    assert plan.protein_link_atom.atom_name == "ND2"
    assert plan.protein_product_residue_name == "ASX"
    assert plan.modifier_product_residue_name == "NAG"
    assert tuple(atom.atom_name for atom in plan.protein_leaving_atoms) == ("HD21",)
    assert plan.pablo_crosslink_requirement.residues == ("ASX", "NAG")
    assert plan.pablo_crosslink_requirement.linking_atoms == (
        "ND2",
        plan.modifier_link_atom.atom_name,
    )
    assert plan.pablo_crosslink_requirement.leaving_atoms[0] == ("HD21",)
    assert len(plan.modifier_leaving_atoms) >= 2


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


def _asn_pdb(tmp_path: Path) -> Path:
    """Create a small Asn residue with explicit ND2 hydrogens."""
    path = tmp_path / "asn.pdb"
    lines = [
        _pdb_atom(1, "N", "ASN", "A", 42, 0.0, 0.0, 0.0, element="N"),
        _pdb_atom(2, "CA", "ASN", "A", 42, 1.0, 0.0, 0.0),
        _pdb_atom(3, "CB", "ASN", "A", 42, 1.5, 0.7, 0.0),
        _pdb_atom(4, "CG", "ASN", "A", 42, 2.0, 0.0, 0.0),
        _pdb_atom(5, "OD1", "ASN", "A", 42, 2.0, -1.2, 0.0, element="O"),
        _pdb_atom(6, "ND2", "ASN", "A", 42, 2.0, 1.2, 0.0, element="N"),
        _pdb_atom(7, "HD21", "ASN", "A", 42, 2.0, 2.1, 0.0, element="H"),
        _pdb_atom(8, "HD22", "ASN", "A", 42, 2.8, 1.2, 0.0, element="H"),
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
