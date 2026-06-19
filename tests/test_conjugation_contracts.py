"""Tests for generic PDB linkage contracts."""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from polyzymd.builders.conjugation._linkage import (
    ExplicitLinkageContract,
    LinkageBond,
    NhsLysModifierLinker,
    PdbAtomSelector,
    ReactiveEndpoint,
    explicit_linkage_contract_from_config,
    resolve_explicit_linkage_contract,
)
from polyzymd.builders.conjugation.polymer import GeneratedPolymerFragment
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord
from polyzymd.config.schema import ConjugationAttachmentConfig


def test_explicit_contract_resolves_multi_residue_modifier_atoms(tmp_path: Path):
    """Explicit selectors should resolve the intended modifier residue."""
    protein_path = _protein_pdb(tmp_path)
    modifier_path = _multi_residue_modifier_pdb(tmp_path)
    contract = _generic_contract()

    plan = resolve_explicit_linkage_contract(protein_path, modifier_path, contract)

    assert plan.protein_link_atom.atom_name == "NZ"
    assert plan.modifier_link_atom.residue_number == 2
    assert plan.modifier_link_atom.atom_name == "C1"
    assert [atom.residue_number for atom in plan.modifier_leaving_atoms] == [2]


def test_leaving_atom_names_are_scoped_to_selected_residue(tmp_path: Path):
    """Duplicate leaving atom names in other residues must not be selected."""
    plan = resolve_explicit_linkage_contract(
        _protein_pdb(tmp_path),
        _multi_residue_modifier_pdb(tmp_path),
        _generic_contract(),
    )

    assert tuple(atom.atom_name for atom in plan.modifier_leaving_atoms) == ("LG",)
    assert tuple(atom.serial for atom in plan.modifier_leaving_atoms) == (103,)


def test_missing_and_ambiguous_selectors_fail_clearly(tmp_path: Path):
    """Selectors should fail clearly when they match zero or multiple atoms."""
    protein_path = _protein_pdb(tmp_path)
    missing_contract = _generic_contract(
        modifier_selector=PdbAtomSelector(
            chain_id="C",
            residue_name="NHS",
            residue_number=99,
            atom_name="C1",
        )
    )

    with pytest.raises(ValueError, match="Expected exactly one modifier link atom"):
        resolve_explicit_linkage_contract(
            protein_path,
            _multi_residue_modifier_pdb(tmp_path),
            missing_contract,
        )

    ambiguous_path = tmp_path / "ambiguous_modifier.pdb"
    ambiguous_path.write_text(
        _pdb_atom(101, "C1", "NHS", "C", 2, 0.0, 0.0, 0.0, element="C", record="HETATM")
        + _pdb_atom(
            102,
            "C1",
            "NHS",
            "C",
            2,
            1.0,
            0.0,
            0.0,
            element="C",
            record="HETATM",
        )
        + _pdb_atom(103, "LG", "NHS", "C", 2, 2.0, 0.0, 0.0, element="O", record="HETATM")
        + "END\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="found 2"):
        resolve_explicit_linkage_contract(protein_path, ambiguous_path, _generic_contract())


def test_resolved_plan_produces_pablo_crosslink_requirement(tmp_path: Path):
    """Resolved plans should expose the exact Pablo crosslink requirement."""
    plan = resolve_explicit_linkage_contract(
        _protein_pdb(tmp_path),
        _multi_residue_modifier_pdb(tmp_path),
        _generic_contract(),
    )

    requirement = plan.pablo_crosslink_requirement
    assert requirement.residues == ("LYX", "NHX")
    assert requirement.linking_atoms == ("NZ", "C1")
    assert requirement.leaving_atoms == (("HZ2", "HZ3"), ("LG",))
    assert requirement.bond_order == 1


def test_nhs_lys_linker_emits_equivalent_generic_plan(tmp_path: Path):
    """The NHS-Lys helper should resolve through the generic contract path."""
    linker = NhsLysModifierLinker(target_residue_number=23)
    modifier = _generated_nhs_modifier()

    plan = linker.resolve_plan(_protein_pdb(tmp_path), modifier)

    assert plan.contract.mechanism_name == "nhs_lys"
    assert plan.protein_link_atom.atom_name == "NZ"
    assert plan.modifier_link_atom.atom_name == "RC"
    assert tuple(atom.atom_name for atom in plan.protein_leaving_atoms) == ("HZ2", "HZ3")
    assert tuple(atom.atom_name for atom in plan.modifier_leaving_atoms) == ("LG",)
    assert plan.pablo_crosslink_requirement.residues == ("LYX", "NHX")
    assert plan.pablo_crosslink_requirement.linking_atoms == ("NZ", "RC")


def test_attachment_config_parses_explicit_multi_residue_link_site():
    """YAML-shaped config should parse explicit linkage selector fields."""
    attachment = ConjugationAttachmentConfig.model_validate(
        {
            "name": "lys23-polymer",
            "site": {
                "chain_id": "A",
                "residue_name": "LYS",
                "residue_number": 23,
                "atom_name": "NZ",
                "insertion_code": "B",
                "atom_serial": 4,
                "atom_index": 3,
            },
            "moiety": {
                "name": "polymer",
                "input_path": "modifier.pdb",
                "link_site": {
                    "chain_id": "C",
                    "residue_name": "NHS",
                    "residue_number": 2,
                    "atom_name": "C1",
                    "insertion_code": "A",
                    "atom_serial": 102,
                    "atom_index": 1,
                },
            },
            "mechanism": {
                "name": "explicit_linkage",
                "product_residues": {"site": "LYX", "moiety": "NHX"},
                "bond": {"site_atom": "NZ", "moiety_atom": "C1", "order": 1},
                "leaving_atoms": {"site": ["HZ2", "HZ3"], "moiety": ["LG"]},
            },
        }
    )

    assert attachment.moiety.link_site is not None
    assert attachment.moiety.link_site.residue_number == 2
    assert attachment.site.insertion_code == "B"
    assert attachment.site.atom_serial == 4
    assert attachment.site.atom_index == 3
    assert attachment.moiety.link_site.insertion_code == "A"
    assert attachment.moiety.link_site.atom_serial == 102
    assert attachment.moiety.link_site.atom_index == 1
    assert attachment.mechanism.name == "explicit_linkage"
    assert attachment.mechanism.product_residues.moiety == "NHX"
    assert attachment.mechanism.leaving_atoms.moiety == ["LG"]


def test_moiety_link_site_requires_residue_disambiguation():
    """Moiety linkage selectors should require chain, residue, and atom fields."""
    with pytest.raises(ValidationError):
        ConjugationAttachmentConfig.model_validate(
            {
                "name": "bad-selector",
                "moiety": {
                    "name": "polymer",
                    "link_site": {"residue_name": "NHS", "atom_name": "C1"},
                },
            }
        )


@pytest.mark.parametrize(
    ("patch", "message"),
    [
        ({"site": {"atom_name": None}}, "site.atom_name"),
        ({"moiety": {"link_site": None}}, "moiety.link_site"),
        ({"moiety": {"smiles": "CC"}}, "moiety.smiles"),
        (
            {
                "moiety": {
                    "polymer_recipe": {
                        "monomers": [
                            {
                                "label": "A",
                                "name": "Acrylate",
                                "residue_name": "ACR",
                                "smiles": "C=C",
                                "probability": 1.0,
                            }
                        ],
                        "length": 2,
                    }
                }
            },
            "polymer_recipe",
        ),
        ({"moiety": {"input_path": "modifier.sdf"}}, ".pdb"),
        (
            {"mechanism": {"product_residues": {"site": None, "moiety": "NHX"}}},
            "product_residues.site",
        ),
    ],
)
def test_explicit_linkage_rejects_non_pdb_or_incomplete_config(patch, message):
    """Explicit linkage config should enforce the v1 PDB-only contract."""
    data = _valid_explicit_attachment_data()
    _deep_update(data, patch)

    with pytest.raises(ValidationError, match=message):
        ConjugationAttachmentConfig.model_validate(data)


def test_explicit_linkage_contract_adapter_maps_config_fields():
    """The config adapter should map selectors, residues, leaving atoms, and bond fields."""
    attachment = ConjugationAttachmentConfig.model_validate(_valid_explicit_attachment_data())

    contract = explicit_linkage_contract_from_config(attachment)

    assert contract.protein_endpoint.selector.chain_id == "A"
    assert contract.protein_endpoint.selector.residue_name == "LYS"
    assert contract.protein_endpoint.selector.residue_number == 23
    assert contract.protein_endpoint.selector.insertion_code == "B"
    assert contract.protein_endpoint.selector.atom_serial == 4
    assert contract.protein_endpoint.selector.atom_index == 3
    assert contract.modifier_endpoint.selector.chain_id == "C"
    assert contract.modifier_endpoint.selector.residue_name == "NHS"
    assert contract.modifier_endpoint.selector.insertion_code == "A"
    assert contract.protein_endpoint.product_residue_name == "LYX"
    assert contract.modifier_endpoint.product_residue_name == "NHX"
    assert contract.protein_endpoint.leaving_atom_names == ("HZ2", "HZ3")
    assert contract.modifier_endpoint.leaving_atom_names == ("LG",)
    assert contract.bond.protein_atom_name == "NZ"
    assert contract.bond.modifier_atom_name == "C1"
    assert contract.bond.bond_order == 2
    assert contract.bond.target_bond_length_angstrom == 1.45


def test_nhs_lys_linker_fails_on_target_residue_name_mismatch(tmp_path: Path):
    """The NHS-Lys helper should honor the generic residue-name selector."""
    linker = NhsLysModifierLinker(target_residue_number=23, target_residue_name="ALA")

    with pytest.raises(ValueError, match="protein link atom"):
        linker.resolve_plan(_protein_pdb(tmp_path), _generated_nhs_modifier())


def test_explicit_hz_contract_fails_against_poc_hydrogen_names(tmp_path: Path):
    """Strict explicit contracts must not reinterpret missing atom names chemically."""
    with pytest.raises(ValueError, match="named HZ2"):
        resolve_explicit_linkage_contract(
            _protein_pdb_with_poc_hydrogen_names(tmp_path),
            _multi_residue_modifier_pdb(tmp_path),
            _generic_contract(),
        )


def _generic_contract(
    *, modifier_selector: PdbAtomSelector | None = None
) -> ExplicitLinkageContract:
    """Build a generic explicit linkage contract for tests."""
    protein_selector = PdbAtomSelector(
        chain_id="A",
        residue_name="LYS",
        residue_number=23,
        atom_name="NZ",
    )
    modifier_selector = modifier_selector or PdbAtomSelector(
        chain_id="C",
        residue_name="NHS",
        residue_number=2,
        atom_name="C1",
    )
    return ExplicitLinkageContract(
        protein_endpoint=ReactiveEndpoint(
            participant="protein",
            selector=protein_selector,
            product_residue_name="LYX",
            leaving_atom_names=("HZ2", "HZ3"),
        ),
        modifier_endpoint=ReactiveEndpoint(
            participant="modifier",
            selector=modifier_selector,
            product_residue_name="NHX",
            leaving_atom_names=("LG",),
        ),
        bond=LinkageBond(protein_atom_name="NZ", modifier_atom_name="C1", bond_order=1),
        mechanism_name="explicit_linkage",
    )


def _valid_explicit_attachment_data() -> dict:
    """Return a valid explicit linkage attachment config mapping."""
    return {
        "name": "lys23-polymer",
        "site": {
            "chain_id": "A",
            "residue_name": "LYS",
            "residue_number": 23,
            "insertion_code": "B",
            "atom_name": "NZ",
            "atom_serial": 4,
            "atom_index": 3,
        },
        "moiety": {
            "name": "polymer",
            "input_path": "modifier.pdb",
            "link_site": {
                "chain_id": "C",
                "residue_name": "NHS",
                "residue_number": 2,
                "insertion_code": "A",
                "atom_name": "C1",
                "atom_serial": 102,
                "atom_index": 1,
            },
        },
        "mechanism": {
            "name": "explicit_linkage",
            "product_residues": {"site": "LYX", "moiety": "NHX"},
            "bond": {
                "site_atom": "NZ",
                "moiety_atom": "C1",
                "order": 2,
                "target_bond_length_angstrom": 1.45,
            },
            "leaving_atoms": {"site": ["HZ2", "HZ3"], "moiety": ["LG"]},
        },
    }


def _deep_update(data: dict, patch: dict) -> None:
    """Apply a nested mapping patch in place for tests."""
    for key, value in patch.items():
        if isinstance(value, dict) and isinstance(data.get(key), dict):
            _deep_update(data[key], value)
        else:
            data[key] = value


def _protein_pdb(tmp_path: Path) -> Path:
    """Create a small protein PDB containing a lysine site."""
    path = tmp_path / "protein.pdb"
    path.write_text(
        _pdb_atom(1, "N", "LYS", "A", 23, 0.0, 0.0, 0.0, element="N")
        + _pdb_atom(2, "CA", "LYS", "A", 23, 1.0, 0.0, 0.0)
        + _pdb_atom(3, "CE", "LYS", "A", 23, 1.5, 0.0, 0.0)
        + _pdb_atom(4, "NZ", "LYS", "A", 23, 2.0, 0.0, 0.0, element="N")
        + _pdb_atom(5, "HZ1", "LYS", "A", 23, 2.0, 0.7, 0.0, element="H")
        + _pdb_atom(6, "HZ2", "LYS", "A", 23, 2.0, -0.7, 0.0, element="H")
        + _pdb_atom(7, "HZ3", "LYS", "A", 23, 2.0, 0.0, 0.7, element="H")
        + _pdb_atom(8, "N", "ALA", "A", 24, 4.0, 0.0, 0.0, element="N")
        + "END\n",
        encoding="utf-8",
    )
    return path


def _protein_pdb_with_poc_hydrogen_names(tmp_path: Path) -> Path:
    """Create a lysine PDB using POC-like noncanonical NZ hydrogen names."""
    path = tmp_path / "protein_poc_hydrogens.pdb"
    path.write_text(
        _pdb_atom(1, "N", "LYS", "A", 23, 0.0, 0.0, 0.0, element="N")
        + _pdb_atom(2, "CA", "LYS", "A", 23, 1.0, 0.0, 0.0)
        + _pdb_atom(3, "CE", "LYS", "A", 23, 1.5, 0.0, 0.0)
        + _pdb_atom(4, "NZ", "LYS", "A", 23, 2.0, 0.0, 0.0, element="N")
        + _pdb_atom(5, "H10", "LYS", "A", 23, 2.0, 0.7, 0.0, element="H")
        + _pdb_atom(6, "H11", "LYS", "A", 23, 2.0, -0.7, 0.0, element="H")
        + _pdb_atom(7, "H13", "LYS", "A", 23, 2.0, 0.0, 0.7, element="H")
        + "END\n",
        encoding="utf-8",
    )
    return path


def _multi_residue_modifier_pdb(tmp_path: Path) -> Path:
    """Create a multi-residue modifier PDB with duplicate leaving names."""
    path = tmp_path / "modifier.pdb"
    path.write_text(
        _pdb_atom(101, "P1", "SBM", "C", 1, 0.0, 0.0, 0.0, record="HETATM")
        + _pdb_atom(102, "C1", "NHS", "C", 2, 1.0, 0.0, 0.0, record="HETATM")
        + _pdb_atom(103, "LG", "NHS", "C", 2, 2.0, 0.0, 0.0, element="O", record="HETATM")
        + _pdb_atom(104, "LG", "EGP", "C", 3, 3.0, 0.0, 0.0, element="O", record="HETATM")
        + "END\n",
        encoding="utf-8",
    )
    return path


def _generated_nhs_modifier() -> GeneratedPolymerFragment:
    """Create a generated NHS modifier for convenience-linker tests."""
    atoms = (
        PdbAtomRecord(
            serial=101,
            atom_index=0,
            atom_name="P1",
            residue_name="SBM",
            chain_id="C",
            residue_number=1,
            x=0.0,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=102,
            atom_index=1,
            atom_name="RC",
            residue_name="NHS",
            chain_id="C",
            residue_number=2,
            x=1.0,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=103,
            atom_index=2,
            atom_name="LG",
            residue_name="NHS",
            chain_id="C",
            residue_number=2,
            x=2.0,
            y=0.0,
            z=0.0,
            element="O",
            record_name="HETATM",
        ),
    )
    return GeneratedPolymerFragment.from_atom_records(
        atoms,
        reactive_atom_serial=102,
        reactive_atom_index=1,
        leaving_atom_serials=(103,),
        leaving_atom_indices=(2,),
        leaving_atom_names=("LG",),
        name="nhs_modifier",
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
    """Format one PDB atom line for contract tests."""
    return (
        f"{record:<6}{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    {x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )
