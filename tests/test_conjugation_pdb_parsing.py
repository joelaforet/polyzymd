"""Tests for canonical conjugation PDB parsing helpers."""

from pathlib import Path

from polyzymd.builders.conjugation.structure.parsing import (
    PdbAtomRecord as ParsingPdbAtomRecord,
)
from polyzymd.builders.conjugation.structure.parsing import (
    parse_pdb_atom_line,
    parse_pdb_atom_records,
    parse_pdb_conect_adjacency,
    parse_pdb_conect_pairs,
    parse_pdb_conect_records,
    parse_pdb_conect_target_serials,
    parse_pdb_link_sides,
    pdb_coordinates,
    pdb_has_conect_records,
)
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord

ATOM_LINE = "ATOM      7  CA ASER A  42B     11.125  -2.500   3.250  0.50 12.34" "           C1+"
HETATM_LINE = "HETATM   12  O1  PEG C   2       1.000   2.000   3.000  1.00  0.00" "           O  "


def test_atom_record_parser_preserves_fixed_width_identity_fields() -> None:
    """PDB atom parsing should preserve identity, coordinates, and ordering fields."""
    atom = parse_pdb_atom_line(ATOM_LINE, atom_index=3)

    assert isinstance(atom, PdbAtomRecord)
    assert atom.serial == 7
    assert atom.atom_index == 3
    assert atom.atom_name == "CA"
    assert atom.alt_loc == "A"
    assert atom.residue_name == "SER"
    assert atom.chain_id == "A"
    assert atom.residue_number == 42
    assert atom.insertion_code == "B"
    assert (atom.x, atom.y, atom.z) == (11.125, -2.5, 3.25)
    assert atom.occupancy == 0.5
    assert atom.temp_factor == 12.34
    assert atom.element == "C"
    assert atom.charge == "1+"
    assert atom.record_name == "ATOM"


def test_atom_record_parser_preserves_zero_occupancy() -> None:
    """PDB atom parsing should not treat valid zero occupancy as missing."""
    line = f"{ATOM_LINE[:54]}  0.00{ATOM_LINE[60:]}"

    atom = parse_pdb_atom_line(line)

    assert atom.occupancy == 0.0
    assert atom.temp_factor == 12.34


def test_atom_record_parser_defaults_blank_occupancy_and_temperature() -> None:
    """Blank occupancy and temperature fields should use PDB-safe defaults."""
    line = f"{ATOM_LINE[:54]}            {ATOM_LINE[66:]}"

    atom = parse_pdb_atom_line(line)

    assert atom.occupancy == 1.0
    assert atom.temp_factor == 0.0


def test_pdb_atom_record_import_identity_is_stable() -> None:
    """Public structure import surfaces should expose the canonical PdbAtomRecord."""
    from polyzymd.builders.conjugation import structure

    assert ParsingPdbAtomRecord is PdbAtomRecord
    assert structure.PdbAtomRecord is PdbAtomRecord
    assert PdbAtomRecord.__module__ == "polyzymd.builders.conjugation.structure.pdb"


def test_atom_file_parser_assigns_source_order_indices(tmp_path: Path) -> None:
    """File parsing should skip non-atom records and assign zero-based source order."""
    pdb_path = tmp_path / "atoms.pdb"
    pdb_path.write_text(f"REMARK ignored\n{ATOM_LINE}\n{HETATM_LINE}\nEND\n", encoding="utf-8")

    atoms = parse_pdb_atom_records(pdb_path)

    assert [atom.serial for atom in atoms] == [7, 12]
    assert [atom.atom_index for atom in atoms] == [0, 1]
    assert [atom.record_name for atom in atoms] == ["ATOM", "HETATM"]
    assert pdb_coordinates(pdb_path) == ((11.125, -2.5, 3.25), (1.0, 2.0, 3.0))


def test_conect_parser_normalizes_unique_pairs_and_preserves_record_semantics(
    tmp_path: Path,
) -> None:
    """CONECT parsing should expose pair, target-list, and adjacency semantics."""
    pdb_path = tmp_path / "conect.pdb"
    pdb_path.write_text(
        "CONECT    7   12   13\n" "CONECT   12    7\n" "CONECT   13   13   XX\n" "END\n",
        encoding="utf-8",
    )

    assert pdb_has_conect_records(pdb_path)
    assert parse_pdb_conect_target_serials("CONECT    7   12   13") == (12, 13)
    assert parse_pdb_conect_pairs(pdb_path) == ((7, 12), (7, 13))
    assert parse_pdb_conect_adjacency(pdb_path) == {7: (12, 13), 12: (7,), 13: (7,)}
    assert parse_pdb_conect_records(pdb_path) == {7: (12, 13), 12: (7,), 13: (13,)}


def test_link_parser_preserves_fixed_width_and_whitespace_fallback() -> None:
    """LINK parsing should support strict PDB columns and existing diagnostic fallback."""
    chars = list(" " * 80)
    chars[0:6] = "LINK  "
    chars[12:16] = "NZ  "
    chars[17:20] = "LYS"
    chars[21:22] = "A"
    chars[22:26] = "  10"
    chars[26:27] = "B"
    chars[42:46] = "C1  "
    chars[47:50] = "PEG"
    chars[51:52] = "C"
    chars[52:56] = "   2"
    chars[56:57] = "D"
    fixed = "".join(chars)
    side_1, side_2 = parse_pdb_link_sides(fixed)

    assert side_1.atom_name == "NZ"
    assert side_1.residue_name == "LYS"
    assert side_1.chain_id == "A"
    assert side_1.residue_number == 10
    assert side_1.insertion_code == "B"
    assert side_2.atom_name == "C1"
    assert side_2.residue_name == "PEG"
    assert side_2.chain_id == "C"
    assert side_2.residue_number == 2
    assert side_2.insertion_code == "D"

    fallback_1, fallback_2 = parse_pdb_link_sides("LINK NZ LYS A 10 C1 PEG C 2")
    assert fallback_1.atom_name == "NZ"
    assert fallback_2.atom_name == "C1"
