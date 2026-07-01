"""Tests for the Polymerist PDB to generated fragment adapter."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.builders.conjugation.polymer.polymerist import generated_fragment_from_polymerist_pdb
from polyzymd.builders.conjugation.polymer.recipe import generate_multi_residue_molecule
from polyzymd.builders.conjugation.structure.pdb import NhsLysPdbAttachment, write_crosslinked_pdb
from tests._support.conjugation_polymer_recipes import sbma_egpma_nhs_recipe


def _pdb_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    *,
    element: str,
) -> str:
    """Format one fixed-width PDB atom line for adapter tests."""
    return (
        f"HETATM{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )


def _conect(source: int, *targets: int) -> str:
    """Format one PDB CONECT record."""
    return f"CONECT{source:5d}" + "".join(f"{target:5d}" for target in targets) + "\n"


def _minimal_sidecar_pdb(tmp_path: Path) -> Path:
    """Write a minimal PDB fixture for sidecar bond-order tests."""
    path = tmp_path / "sidecar_fixture.pdb"
    path.write_text(
        _pdb_atom(1, "C1", "MOL", "C", 1, 0.0, 0.0, 0.0, element="C")
        + _pdb_atom(2, "O1", "MOL", "C", 1, 1.2, 0.0, 0.0, element="O")
        + _pdb_atom(3, "C2", "MOL", "C", 1, -1.2, 0.0, 0.0, element="C")
        + _conect(1, 2, 3)
        + _conect(2, 1)
        + _conect(3, 1)
        + "END\n",
        encoding="utf-8",
    )
    return path


def _write_sidecar_sdf(
    path: Path,
    bonds: tuple[tuple[int, int, int], ...],
    *,
    formal_charges: dict[int, int] | None = None,
) -> None:
    """Write an RDKit SDF fixture with caller-specified one-based bond orders."""
    Chem = pytest.importorskip("rdkit.Chem")
    bond_types = {
        0: Chem.BondType.UNSPECIFIED,
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }

    editable = Chem.RWMol()
    formal_charges = formal_charges or {}
    for index, symbol in enumerate(("C", "O", "C"), start=1):
        atom = Chem.Atom(symbol)
        atom.SetFormalCharge(formal_charges.get(index, 0))
        editable.AddAtom(atom)
    for atom_1, atom_2, order in bonds:
        editable.AddBond(atom_1 - 1, atom_2 - 1, bond_types[order])
    mol = editable.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.MolToMolFile(mol, str(path))


def _fragment_bond_order(fragment, atom_1: int, atom_2: int) -> float:
    """Return a generated-fragment bond order by serial pair."""
    serials = tuple(sorted((atom_1, atom_2)))
    matches = [order for left, right, order in fragment.bond_orders if (left, right) == serials]
    assert len(matches) == 1
    return matches[0]


def _nhs_group(*, serial_offset: int = 0, residue_number: int = 1, x_offset: float = 0.0) -> str:
    """Create one small NHS-like ester residue with explicit connectivity."""
    serial = serial_offset
    return (
        _pdb_atom(serial + 1, "C1", "NHS", "C", residue_number, x_offset, 0.0, 0.0, element="C")
        + _pdb_atom(
            serial + 2,
            "O1",
            "NHS",
            "C",
            residue_number,
            x_offset + 1.2,
            0.0,
            0.0,
            element="O",
        )
        + _pdb_atom(
            serial + 3,
            "O2",
            "NHS",
            "C",
            residue_number,
            x_offset,
            1.2,
            0.0,
            element="O",
        )
        + _pdb_atom(
            serial + 4,
            "N1",
            "NHS",
            "C",
            residue_number,
            x_offset,
            2.5,
            0.0,
            element="N",
        )
        + _pdb_atom(
            serial + 5,
            "C2",
            "NHS",
            "C",
            residue_number,
            x_offset,
            3.8,
            0.0,
            element="C",
        )
        + _conect(serial + 1, serial + 2, serial + 3)
        + _conect(serial + 2, serial + 1)
        + _conect(serial + 3, serial + 1, serial + 4)
        + _conect(serial + 4, serial + 3, serial + 5)
        + _conect(serial + 5, serial + 4)
    )


def _polymerist_like_pdb(tmp_path: Path) -> Path:
    """Write a small Polymerist-like PDB with one NHS residue and one side residue."""
    path = tmp_path / "polymerist_like.pdb"
    path.write_text(
        _nhs_group()
        + _pdb_atom(6, "P1", "SBM", "C", 2, 5.0, 0.0, 0.0, element="C")
        + _conect(1, 6)
        + _conect(6, 1)
        + "END\n"
    )
    return path


def _protein_pdb(tmp_path: Path) -> Path:
    """Create a small protein PDB containing a reactive lysine."""
    path = tmp_path / "protein.pdb"
    path.write_text(
        "ATOM      1 N    LYS A  23       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2 CA   LYS A  23       1.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3 CE   LYS A  23       1.500   0.000   0.000  1.00  0.00           C\n"
        "ATOM      4 NZ   LYS A  23       2.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      5 HZ1  LYS A  23       2.000   0.700   0.000  1.00  0.00           H\n"
        "ATOM      6 HZ2  LYS A  23       2.000  -0.700   0.000  1.00  0.00           H\n"
        "ATOM      7 HZ3  LYS A  23       2.000   0.000   0.700  1.00  0.00           H\n"
        "END\n"
    )
    return path


def test_polymerist_like_pdb_converts_to_generated_fragment(tmp_path):
    """Synthetic NHS-containing PDB files should convert to generated fragments."""
    fragment = generated_fragment_from_polymerist_pdb(_polymerist_like_pdb(tmp_path), sequence="CA")

    assert len(fragment.atoms) == 6
    assert fragment.sequence == "CA"
    assert [residue.residue_name for residue in fragment.residues] == ["NHS", "SBM"]
    assert fragment.reactive_atom_serial == 1
    assert fragment.reactive_atom_index == 0
    assert fragment.reactive_atom_name == "C1"


def test_reactive_carbon_leaving_group_and_bonds_are_mapped_to_pdb_records(tmp_path):
    """Reactive data and bonds should be translated from RDKit indices to PDB records."""
    fragment = generated_fragment_from_polymerist_pdb(_polymerist_like_pdb(tmp_path))

    assert fragment.leaving_atom_serials == (3, 4, 5)
    assert fragment.leaving_atom_indices == (2, 3, 4)
    assert fragment.leaving_atom_names == ("O2", "N1", "C2")
    assert (1, 2) in fragment.bonds
    assert (1, 3) in fragment.bonds
    assert (3, 4) in fragment.bonds
    assert (1, 2, 1.0) in fragment.bond_orders
    assert all(
        isinstance(source, int) and isinstance(target, int) for source, target in fragment.bonds
    )


def test_generated_fragment_uses_aligned_sidecar_positive_bond_orders(tmp_path):
    """Generated fragments should carry positive bond orders from aligned SDF sidecars."""
    pdb_path = _minimal_sidecar_pdb(tmp_path)
    _write_sidecar_sdf(pdb_path.with_suffix(".sdf"), ((1, 2, 2), (1, 3, 1)))

    fragment = generated_fragment_from_polymerist_pdb(
        pdb_path,
        reactive_atom_name="C1",
        leaving_atom_names=("O1",),
    )

    assert _fragment_bond_order(fragment, 1, 2) == 2.0
    assert _fragment_bond_order(fragment, 1, 3) == 1.0
    assert all(order > 0.0 for *_atoms, order in fragment.bond_orders)


def test_generated_fragment_preserves_sdf_formal_charges(tmp_path):
    """Generated fragments should carry non-zero formal charges from aligned SDF sidecars."""
    pdb_path = _minimal_sidecar_pdb(tmp_path)
    _write_sidecar_sdf(
        pdb_path.with_suffix(".sdf"),
        ((1, 2, 1), (1, 3, 1)),
        formal_charges={2: -1, 3: 1},
    )

    fragment = generated_fragment_from_polymerist_pdb(
        pdb_path,
        reactive_atom_name="C1",
        leaving_atom_names=("O1",),
    )

    charges = {atom.atom_name: atom.formal_charge for atom in fragment.atoms}
    assert charges["O1"] == -1
    assert charges["C2"] == 1


def test_generated_fragment_rejects_zero_order_sidecar_bonds(tmp_path):
    """Query/zero-order SDF sidecars should fail instead of being chemistry-repaired."""
    pdb_path = _minimal_sidecar_pdb(tmp_path)
    _write_sidecar_sdf(pdb_path.with_suffix(".sdf"), ((1, 2, 0), (1, 3, 1)))

    with pytest.raises(ValueError, match="under-specified zero/unknown bond orders"):
        generated_fragment_from_polymerist_pdb(
            pdb_path,
            reactive_atom_name="C1",
            leaving_atom_names=("O1",),
        )


def test_multiple_nhs_groups_raise_ambiguity_without_selector(tmp_path):
    """Multiple NHS-like groups should require a residue selector."""
    path = tmp_path / "two_nhs_groups.pdb"
    path.write_text(
        _nhs_group(residue_number=1)
        + _nhs_group(serial_offset=10, residue_number=2, x_offset=8.0)
        + "END\n"
    )

    with pytest.raises(ValueError, match="Ambiguous NHS ester reactive group"):
        generated_fragment_from_polymerist_pdb(path)

    selected = generated_fragment_from_polymerist_pdb(path, reactive_residue_number=2)
    assert selected.reactive_atom_serial == 11
    assert selected.leaving_atom_serials == (13, 14, 15)


def test_generated_fragment_feeds_crosslinked_pdb_writer(tmp_path):
    """The adapter output should convert to a placed fragment for PDB assembly."""
    fragment = generated_fragment_from_polymerist_pdb(_polymerist_like_pdb(tmp_path))
    output_path = tmp_path / "assembled.pdb"

    result = write_crosslinked_pdb(
        _protein_pdb(tmp_path),
        fragment.to_placed_fragment(),
        NhsLysPdbAttachment(
            target_chain="A",
            target_residue_number=23,
            nz_hydrogen_atom_names_to_remove=("HZ2", "HZ3"),
        ),
        output_path,
    )

    output_text = output_path.read_text()
    assert result.removed_polymer_atom_count == 3
    assert result.added_conect_pair[1] > result.added_conect_pair[0]
    assert "NHX" in output_text
    assert "LYX" in output_text


def test_real_polymerist_generation_pdb_converts_to_generated_fragment(tmp_path):
    """Real Polymerist PDB generation should feed the generated-fragment adapter."""
    pytest.importorskip("polymerist", exc_type=ImportError)

    recipe = sbma_egpma_nhs_recipe(length=3, seed=5, reactive_monomer_index=1)
    result = generate_multi_residue_molecule(recipe, tmp_path / "polymerist", max_retries=1)
    assert result.pdb_path is not None

    fragment = generated_fragment_from_polymerist_pdb(
        result.pdb_path,
        recipe=recipe,
        sequence=result.sequence,
    )

    assert len(fragment.atoms) > 0
    assert len(fragment.bonds) > 0
    assert len(fragment.bond_orders) > 0
    assert fragment.reactive_atom_serial is not None
    assert fragment.reactive_atom_index is not None
    assert fragment.reactive_atom_name is not None
    assert fragment.leaving_atom_serials
    assert fragment.leaving_atom_indices
