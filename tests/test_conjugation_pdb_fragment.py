"""Tests for generic residue-resolved PDB-fragment ingestion."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.builders.conjugation._pdb_fragment import load_pdb_fragment
from polyzymd.builders.conjugation.reactions.n_glycosylation import glygen_pdb_profile_from_fragment


def test_pdb_fragment_loader_accepts_multi_residue_fragment_without_roh(tmp_path: Path) -> None:
    """Generic PDB-fragment ingestion should not require GlyGen ROH chemistry."""
    fragment_path = _write_generic_fragment(tmp_path / "fragment.pdb")

    result = load_pdb_fragment(fragment_path)

    assert result.connectivity_provenance == "conect"
    assert [item["residue_name"] for item in result.residue_mapping] == ["AAA", "BBB"]
    assert result.serial_bonds == ((1, 2), (2, 3), (3, 4))


def test_n_glycosylation_profile_rejects_generic_fragment_without_roh(tmp_path: Path) -> None:
    """N-glycosylation should reject generic PDB fragments without its profile."""
    fragment = load_pdb_fragment(_write_generic_fragment(tmp_path / "fragment.pdb"))

    with pytest.raises(ValueError, match="ROH residue"):
        glygen_pdb_profile_from_fragment(fragment)


def _write_generic_fragment(path: Path) -> Path:
    """Write a connected two-residue PDB fragment with no GlyGen/GlyCAM ROH residue."""
    lines = [
        _atom_line(1, "C1", "AAA", 1, 0.0, 0.0, 0.0, "C"),
        _atom_line(2, "O1", "AAA", 1, 1.2, 0.0, 0.0, "O"),
        _atom_line(3, "C2", "BBB", 2, 2.4, 0.0, 0.0, "C"),
        _atom_line(4, "O2", "BBB", 2, 3.6, 0.0, 0.0, "O"),
        "CONECT    1    2\n",
        "CONECT    2    1    3\n",
        "CONECT    3    2    4\n",
        "CONECT    4    3\n",
        "END\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")
    return path


def _atom_line(
    serial: int,
    atom_name: str,
    residue_name: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    """Return one fixed-width PDB HETATM line."""
    return (
        f"HETATM{serial:5d} {atom_name:<4} {residue_name:>3} C{residue_number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n"
    )
