"""Tests for config-driven conjugation system workflow helpers."""

from __future__ import annotations

from pathlib import Path

from polyzymd.builders.conjugation.system_workflow import (
    _apply_pdb_atom_names_to_topology,
    _restore_pdb_atom_name_fields,
)


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
