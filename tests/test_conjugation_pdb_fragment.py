"""Tests for generic residue-resolved PDB-fragment ingestion."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation._pdb_fragment import load_pdb_fragment
from polyzymd.builders.conjugation.reactions.n_glycosylation import (
    residue_resolved_glycan_pdb_profile_from_fragment,
)
from polyzymd.builders.conjugation.structure.parsing import (
    parse_pdb_atom_records,
    parse_pdb_conect_pairs,
)

G42666_STYLE_CONECT_PATH = (
    Path(__file__).resolve().parent / "data/conjugation/glygen/synthetic_g42666_style_conect.pdb"
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

    with pytest.raises(ValueError, match="No glycan anomeric motif"):
        residue_resolved_glycan_pdb_profile_from_fragment(fragment)


def test_n_glycosylation_profile_accepts_generic_residue_local_hydroxyl(
    tmp_path: Path,
) -> None:
    """N-glycosylation should accept non-source-specific residue-resolved glycans."""
    fragment = load_pdb_fragment(_write_generic_reducing_end_glycan(tmp_path / "external.pdb"))

    result = residue_resolved_glycan_pdb_profile_from_fragment(fragment)

    assert result.reducing_c1_serial == 4
    assert result.leaving_atom_indices == (0, 1)
    assert result.leaving_group_representation == "separate_residue_hydroxyl_cap"
    assert [item.residue_name for item in result.fragment.residues] == ["ROH", "QAA", "Z9B"]
    assert result.fragment.reactive_atom_serial == 4
    assert result.fragment.leaving_atom_serials == (1, 2)


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
            atom_name="X4",
            insertion_code="",
            atom_serial=None,
            atom_index=None,
        ),
        leaving_atom_names=("O1", "HO1"),
    )

    assert result.fragment.reactive_atom_serial == 4
    assert result.fragment.leaving_atom_serials == (1, 2)


def test_n_glycosylation_profile_rejects_configured_leaving_atom_mismatch(
    tmp_path: Path,
) -> None:
    """Configured leaving atom names should be checked against the graph result."""
    fragment = load_pdb_fragment(_write_generic_reducing_end_glycan(tmp_path / "external.pdb"))

    with pytest.raises(ValueError, match="Configured leaving atoms"):
        residue_resolved_glycan_pdb_profile_from_fragment(
            fragment,
            leaving_atom_names=("O9", "H9"),
        )


def test_n_glycosylation_profile_is_not_source_atom_name_dependent(tmp_path: Path) -> None:
    """Graph-first motif detection should not require ROH, C1, O1, or HO1 names."""
    fragment = load_pdb_fragment(_write_name_independent_reducing_end(tmp_path / "named.pdb"))

    result = residue_resolved_glycan_pdb_profile_from_fragment(fragment)

    assert result.reducing_c1_serial == 4
    assert result.fragment.leaving_atom_serials == (1, 2)
    assert result.leaving_group_representation == "separate_residue_hydroxyl_cap"


def test_pdb_fragment_loader_preserves_explicit_atom_charges(tmp_path: Path) -> None:
    """Nonblank PDB charge fields should remain authoritative per atom."""
    path = _write_nitro_fragment(tmp_path / "explicit_nitro.pdb", explicit_charges=True)

    result = load_pdb_fragment(path)
    fragment = result.to_fragment(reactive_atom_serial=1)

    assert dict(result.serial_formal_charges) == {2: 1, 3: -1}
    assert {atom.serial: atom.formal_charge for atom in fragment.atoms if atom.formal_charge} == {
        2: 1,
        3: -1,
    }


def test_pdb_fragment_loader_propagates_assigned_blank_atom_charges(tmp_path: Path) -> None:
    """Blank PDB charge fields should accept RDKit charges assigned on the exact graph."""
    path = _write_nitro_fragment(tmp_path / "blank_nitro.pdb", explicit_charges=False)

    result = load_pdb_fragment(path)
    fragment = result.to_fragment(reactive_atom_serial=1)

    assigned = dict(result.serial_formal_charges)
    propagated = {atom.serial: atom.formal_charge for atom in fragment.atoms if atom.formal_charge}
    assert assigned[2] == 1
    assert sorted(assigned.values()) == [-1, 1]
    assert propagated == assigned


def test_n_glycosylation_profile_rejects_terminal_retained_oxygen(tmp_path: Path) -> None:
    """Terminal non-hydroxyl oxygens should not be accepted as retained ring O atoms."""
    path = _write_generic_reducing_end_glycan(tmp_path / "terminal_o.pdb")
    path.write_text(
        path.read_text(encoding="utf-8").replace("CONECT   13   14\n", ""),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="No glycan anomeric motif|bond orders could not"):
        residue_resolved_glycan_pdb_profile_from_fragment(load_pdb_fragment(path))


def test_n_glycosylation_profile_rejects_acyclic_retained_oxygen(tmp_path: Path) -> None:
    """Retained O candidates need an alternate heavy path back to anomeric carbon."""
    path = _write_generic_reducing_end_glycan(tmp_path / "acyclic_o.pdb")
    path.write_text(
        path.read_text(encoding="utf-8").replace("CONECT   11   13\n", ""),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="No glycan anomeric motif|disconnected"):
        residue_resolved_glycan_pdb_profile_from_fragment(load_pdb_fragment(path))


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
    residue_map = {"4YB": ("QAA", 1), "0YB": ("Z9B", 2), "ROH": ("ROH", 3)}
    lines = []
    for atom in parse_pdb_atom_records(G42666_STYLE_CONECT_PATH):
        residue_name, residue_number = residue_map[atom.residue_name]
        atom_name = f"X{atom.serial}" if atom.residue_name != "ROH" else atom.atom_name
        lines.append(
            _atom_line(
                atom.serial,
                atom_name,
                residue_name,
                residue_number,
                atom.x,
                atom.y,
                atom.z,
                atom.element,
            )
        )
    lines.extend(
        _conect(left, right) for left, right in parse_pdb_conect_pairs(G42666_STYLE_CONECT_PATH)
    )
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")
    return path


def _write_name_independent_reducing_end(path: Path) -> Path:
    """Write a graph-valid reducing end with arbitrary atom and residue names."""
    return _write_generic_reducing_end_glycan(path)


def _write_nitro_fragment(path: Path, *, explicit_charges: bool) -> Path:
    """Write a neutral nitromethane graph for formal-charge propagation tests."""
    n_charge = "1+" if explicit_charges else ""
    o_charge = "1-" if explicit_charges else ""
    lines = [
        _atom_line(1, "C1", "NTR", 1, 0.0, 0.0, 0.0, "C"),
        _atom_line(2, "N1", "NTR", 1, 1.4, 0.0, 0.0, "N", charge=n_charge),
        _atom_line(3, "O1", "NTR", 1, 2.6, 0.0, 0.0, "O", charge=o_charge),
        _atom_line(4, "O2", "NTR", 1, 1.4, 1.2, 0.0, "O"),
        _atom_line(5, "H1", "NTR", 1, -0.5, 0.9, 0.0, "H"),
        _atom_line(6, "H2", "NTR", 1, -0.5, -0.9, 0.0, "H"),
        _atom_line(7, "H3", "NTR", 1, 0.0, 0.0, 1.0, "H"),
        "CONECT    1    2    5    6    7\n",
        "CONECT    2    3    4\n",
        "END\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")
    return path


def _conect(source: int, *targets: int) -> str:
    """Format one fixed-width CONECT fixture line."""
    return f"CONECT{source:5d}{''.join(f'{target:5d}' for target in targets)}\n"


def _atom_line(
    serial: int,
    atom_name: str,
    residue_name: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    element: str,
    charge: str = "",
) -> str:
    """Return one fixed-width PDB HETATM line."""
    return (
        f"HETATM{serial:5d} {atom_name:<4} {residue_name:>3} C{residue_number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n"
    )[:-1] + f"{charge:>2}\n"
