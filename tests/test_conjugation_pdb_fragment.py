"""Tests for generic residue-resolved PDB-fragment ingestion."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation._pdb_fragment import load_pdb_fragment
from polyzymd.builders.conjugation.reactions.n_glycosylation import (
    residue_resolved_glycan_pdb_profile_from_fragment,
)


def test_pdb_fragment_loader_accepts_multi_residue_fragment_without_cap(tmp_path: Path) -> None:
    """Generic PDB-fragment ingestion should not require a hydroxyl cap residue."""
    fragment_path = _write_generic_fragment(tmp_path / "fragment.pdb")

    result = load_pdb_fragment(fragment_path)

    assert result.connectivity_provenance == "conect"
    assert [item["residue_name"] for item in result.residue_mapping] == ["AAA", "BBB"]
    assert result.serial_bonds == ((1, 2), (2, 3), (3, 4))


def test_n_glycosylation_profile_rejects_malformed_generic_fragment(tmp_path: Path) -> None:
    """N-glycosylation should reject fragments without C1 hydroxyl chemistry."""
    fragment = load_pdb_fragment(_write_generic_fragment(tmp_path / "fragment.pdb"))

    with pytest.raises(ValueError, match="explicit hydroxyl O/H group"):
        residue_resolved_glycan_pdb_profile_from_fragment(fragment)


def test_n_glycosylation_profile_accepts_generic_residue_local_hydroxyl(
    tmp_path: Path,
) -> None:
    """N-glycosylation should accept non-source-specific residue-resolved glycans."""
    fragment = load_pdb_fragment(_write_generic_reducing_end_glycan(tmp_path / "external.pdb"))

    result = residue_resolved_glycan_pdb_profile_from_fragment(fragment)

    assert result.reducing_c1_serial == 1
    assert result.leaving_atom_indices == (1, 2)
    assert result.leaving_group_representation == "local_oh"
    assert [item.residue_name for item in result.fragment.residues] == ["QAA", "Z9B"]
    assert result.fragment.reactive_atom_serial == 1
    assert result.fragment.leaving_atom_serials == (2, 3)


def test_n_glycosylation_profile_honors_configured_link_site_and_leaving_atoms(
    tmp_path: Path,
) -> None:
    """N-glycosylation should validate mechanism-owned selectors structurally."""
    fragment = load_pdb_fragment(_write_generic_reducing_end_glycan(tmp_path / "external.pdb"))

    result = residue_resolved_glycan_pdb_profile_from_fragment(
        fragment,
        link_site=SimpleNamespace(
            chain_id="C",
            residue_name="QAA",
            residue_number=1,
            atom_name="C1",
            insertion_code="",
            atom_serial=None,
            atom_index=None,
        ),
        leaving_atom_names=("O1", "HO1"),
    )

    assert result.fragment.reactive_atom_serial == 1
    assert result.fragment.leaving_atom_serials == (2, 3)


def test_n_glycosylation_profile_rejects_configured_leaving_atom_mismatch(
    tmp_path: Path,
) -> None:
    """Configured leaving atom names should be checked against the graph result."""
    fragment = load_pdb_fragment(_write_generic_reducing_end_glycan(tmp_path / "external.pdb"))

    with pytest.raises(ValueError, match="Configured moiety leaving atoms"):
        residue_resolved_glycan_pdb_profile_from_fragment(
            fragment,
            leaving_atom_names=("O9", "H9"),
        )


def _write_generic_fragment(path: Path) -> Path:
    """Write a connected two-residue PDB fragment with no reducing-end hydroxyl."""
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


def _write_generic_reducing_end_glycan(path: Path) -> Path:
    """Write a non-source-specific two-residue glycan with local C1-OH chemistry."""
    lines = [
        _atom_line(1, "C1", "QAA", 1, 0.0, 0.0, 0.0, "C"),
        _atom_line(2, "O1", "QAA", 1, 1.2, 0.0, 0.0, "O"),
        _atom_line(3, "HO1", "QAA", 1, 1.8, 0.0, 0.0, "H"),
        _atom_line(4, "O5", "QAA", 1, -1.2, 0.0, 0.0, "O"),
        _atom_line(5, "C5", "QAA", 1, -1.8, 0.8, 0.0, "C"),
        _atom_line(6, "C2", "QAA", 1, 0.0, 1.2, 0.0, "C"),
        _atom_line(7, "O3", "QAA", 1, 0.0, 2.4, 0.0, "O"),
        _atom_line(8, "C4", "Z9B", 2, 0.0, 3.6, 0.0, "C"),
        _atom_line(9, "O4", "Z9B", 2, 0.0, 4.8, 0.0, "O"),
        "CONECT    1    2    4    6\n",
        "CONECT    2    1    3\n",
        "CONECT    3    2\n",
        "CONECT    4    1    5\n",
        "CONECT    5    4\n",
        "CONECT    6    1    7\n",
        "CONECT    7    6    8\n",
        "CONECT    8    7    9\n",
        "CONECT    9    8\n",
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
