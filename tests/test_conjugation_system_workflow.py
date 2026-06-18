"""Tests for config-driven conjugation system workflow helpers."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.contracts import PabloCrosslinkRequirement
from polyzymd.builders.conjugation.system_workflow import (
    _apply_pdb_atom_names_to_topology,
    _policy_with_resolved_crosslink,
    _require_supported_coordinate_backend,
    _restore_pdb_atom_name_fields,
)
from polyzymd.config.schema import ConjugationCcdPabloPolicyConfig


class _AtomDouble:
    def __init__(self, name: str) -> None:
        self.name = name
        self.metadata: dict[str, str] = {}


class _TopologyDouble:
    def __init__(self, atoms: list[_AtomDouble]) -> None:
        self._atoms = atoms

    def atoms(self) -> tuple[_AtomDouble, ...]:
        return tuple(self._atoms)


def test_restore_pdb_atom_name_fields_updates_only_template_atom_count(tmp_path: Path):
    """Conjugate template names should be restored without touching solvent names."""
    template = tmp_path / "linked.pdb"
    target = tmp_path / "solvated.pdb"
    template.write_text(
        "".join(
            [
                _pdb_line(1, " N  ", "ALA", "A", 1),
                _pdb_line(2, " CA ", "ALA", "A", 1),
            ]
        ),
        encoding="utf-8",
    )
    target.write_text(
        "".join(
            [
                _pdb_line(1, " N1x", "ALA", "A", 1),
                _pdb_line(2, " C1x", "ALA", "A", 1),
                _pdb_line(3, " O  ", "HOH", "D", 1),
            ]
        ),
        encoding="utf-8",
    )

    restored = _restore_pdb_atom_name_fields(target, template)
    lines = target.read_text(encoding="utf-8").splitlines()

    assert restored == 2
    assert lines[0][12:16] == " N  "
    assert lines[1][12:16] == " CA "
    assert lines[2][12:16] == " O  "


def test_apply_pdb_atom_names_to_topology_uses_same_order_template(tmp_path: Path):
    """OpenFF topology atom names should be reset from the linked PDB template."""
    template = tmp_path / "linked.pdb"
    template.write_text(
        "".join(
            [
                _pdb_line(1, " N  ", "ALA", "A", 1),
                _pdb_line(2, " CA ", "ALA", "A", 1),
            ]
        ),
        encoding="utf-8",
    )
    atoms = [_AtomDouble("N1x"), _AtomDouble("C1x")]

    _apply_pdb_atom_names_to_topology(_TopologyDouble(atoms), template)

    assert [atom.name for atom in atoms] == ["N", "CA"]
    assert [atom.metadata["atom_name"] for atom in atoms] == ["N", "CA"]


def test_policy_with_resolved_crosslink_uses_product_state_leaving_atoms():
    """Generated Pablo policies should not remove atoms from product-state PDBs."""
    requirement = PabloCrosslinkRequirement(
        residues=("LYX", "NHX"),
        linking_atoms=("NZ", "C047"),
        leaving_atoms=(("H11", "H13"), ("O020",)),
        bond_order=1,
    )
    policy = ConjugationCcdPabloPolicyConfig(crosslinks=[])

    generated = _policy_with_resolved_crosslink(
        policy,
        SimpleNamespace(pablo_crosslink_requirement=requirement),
    )

    assert generated.crosslinks[0].leaving_atoms == ((), ())


def test_coordinate_backend_gate_allows_explicit_nhs_lys_mechanism():
    """The system workflow should only route the named implemented backend to NHS-Lys."""
    attachment = SimpleNamespace(
        mechanism=SimpleNamespace(name="nhs_lys_amide", reaction_smarts=None)
    )

    _require_supported_coordinate_backend(attachment)


def test_coordinate_backend_gate_preflights_generic_smarts_then_blocks():
    """Generic reaction SMARTS should not silently enter the NHS-Lys coordinate path."""
    attachment = SimpleNamespace(
        mechanism=SimpleNamespace(
            name="generic_amide",
            reaction_smarts="[N:1]([H:2]).[C:3](=[O:4])[O:5]>>[N:1][C:3](=[O:4])",
            atom_roles=[
                {"map_number": 1, "participant": "site", "role": "linking"},
                {"map_number": 2, "participant": "site", "role": "leaving"},
                {"map_number": 3, "participant": "moiety", "role": "linking"},
                {"map_number": 5, "participant": "moiety", "role": "leaving"},
            ],
        )
    )

    with pytest.raises(NotImplementedError, match="generic SMARTS preflight") as excinfo:
        _require_supported_coordinate_backend(attachment)

    message = str(excinfo.value)
    assert "1 added" in message
    assert "coordinate surgery only for mechanism 'nhs_lys_amide'" in message


def test_coordinate_backend_gate_blocks_unspecified_mechanisms_without_smarts():
    """Unsupported mechanisms without SMARTS should fail before polymer generation."""
    attachment = SimpleNamespace(mechanism=SimpleNamespace(name="custom", reaction_smarts=None))

    with pytest.raises(NotImplementedError, match="currently implements coordinate surgery only"):
        _require_supported_coordinate_backend(attachment)


def _pdb_line(
    serial: int,
    atom_name_field: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
) -> str:
    return (
        f"ATOM  {serial:5d} {atom_name_field[:4]} {residue_name:>3} {chain_id}"
        f"{residue_number:4d}       0.000   0.000   0.000  1.00  0.00           C  \n"
    )
