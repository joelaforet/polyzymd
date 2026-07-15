"""Tests for GlyGen/GlyCAM residue-resolved glycan PDB ingestion."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from polyzymd.builders.conjugation import ConjugatedPolymerSystemSettings as PublicSettings
from polyzymd.builders.conjugation._moiety_provider import validate_moiety_source_config
from polyzymd.builders.conjugation.glygen_pdb import load_glygen_glycan_pdb
from polyzymd.builders.conjugation.system_workflow import (
    ConjugatedPolymerSystemSettings,
    _build_attachment_spec,
    _build_glygen_coordinate_only_result,
)


def test_glygen_loader_uses_conect_and_preserves_residues(tmp_path: Path) -> None:
    """GlyGen loader should validate CONECT graphs and strict ROH leaving atoms."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)

    result = load_glygen_glycan_pdb(glycan_path)

    assert result.connectivity_provenance == "conect"
    assert result.reducing_c1_serial == 1
    assert result.fragment.leaving_atom_serials == (8, 9)
    assert [residue.residue_name for residue in result.fragment.residues] == ["NAG", "MAN", "ROH"]
    assert any(item.plausible_glycosidic for item in result.linkage_diagnostics)


def test_conjugation_facade_exports_workflow_settings() -> None:
    """Public conjugation facade should expose workflow settings."""
    assert PublicSettings is ConjugatedPolymerSystemSettings


def test_glygen_loader_infers_coordinates_without_conect(tmp_path: Path) -> None:
    """GlyGen loader should report coordinate-inferred provenance without CONECT."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_no_conect.pdb", include_conect=False)

    result = load_glygen_glycan_pdb(glycan_path)

    assert result.connectivity_provenance == "coordinate_inferred"
    assert result.reducing_c1_serial == 1


def test_glygen_loader_rejects_missing_roh(tmp_path: Path) -> None:
    """GlyGen loader should fail actionably when ROH labels are absent."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad.pdb", include_conect=True)
    glycan_path.write_text(
        glycan_path.read_text(encoding="utf-8").replace("ROH", "BAD"), encoding="utf-8"
    )

    with pytest.raises(ValueError, match="ROH residue"):
        load_glygen_glycan_pdb(glycan_path)


def test_glygen_loader_rejects_wrong_roh_ho1_element(tmp_path: Path) -> None:
    """GlyGen loader should require ROH:HO1 to be a hydrogen atom."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad_ho1.pdb", include_conect=True)
    _replace_atom_element(glycan_path, serial=9, element="O")

    with pytest.raises(ValueError, match="HO1 element H"):
        load_glygen_glycan_pdb(glycan_path)


def test_glygen_loader_rejects_wrong_roh_o1_element(tmp_path: Path) -> None:
    """GlyGen loader should require ROH:O1 to be an oxygen atom."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad_o1.pdb", include_conect=True)
    _replace_atom_element(glycan_path, serial=8, element="C")

    with pytest.raises(ValueError, match="O1 element O"):
        load_glygen_glycan_pdb(glycan_path)


def test_glygen_loader_rejects_wrong_reducing_c1_element(tmp_path: Path) -> None:
    """GlyGen loader should require the reducing C1 candidate to be carbon."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_bad_c1.pdb", include_conect=True)
    _replace_atom_element(glycan_path, serial=1, element="O")

    with pytest.raises(ValueError, match="unique reducing sugar C1"):
        load_glygen_glycan_pdb(glycan_path)


def test_glygen_loader_rejects_missing_exact_roh_bond(tmp_path: Path) -> None:
    """GlyGen loader should prove the exact ROH:HO1-ROH:O1-C1 path."""
    glycan_path = _write_glygen_fixture(
        tmp_path / "glycan_missing_roh_bond.pdb", include_conect=True
    )
    text = glycan_path.read_text(encoding="utf-8")
    text = text.replace("CONECT    8    9\n", "CONECT    1    9\n")
    glycan_path.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match="exact ROH:HO1-ROH:O1-C1"):
        load_glygen_glycan_pdb(glycan_path)


def test_glygen_loader_rejects_extra_roh_o1_carbon_bond(tmp_path: Path) -> None:
    """GlyGen loader should not accept unrelated C-O bonds around ROH:O1."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_extra_roh_bond.pdb", include_conect=True)
    text = glycan_path.read_text(encoding="utf-8")
    text = text.replace("CONECT    8    9\n", "CONECT    8    9    3\n")
    glycan_path.write_text(text, encoding="utf-8")

    with pytest.raises(ValueError, match="bonded exactly"):
        load_glygen_glycan_pdb(glycan_path)


def test_coordinate_inferred_loader_rejects_false_hydrogen_proximity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Coordinate-inferred validation should reject impossible H-H bonds."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_false_hh.pdb", include_conect=False)
    _patch_coordinate_bonds(
        monkeypatch,
        atom_count=9,
        bonds=((1, 2), (1, 3), (1, 8), (3, 4), (4, 6), (6, 5), (6, 7), (8, 9), (9, 9)),
    )

    with pytest.raises(ValueError, match="unsafe element bond H-H"):
        load_glygen_glycan_pdb(glycan_path)


def test_coordinate_inferred_loader_rejects_overbonded_atoms(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Coordinate-inferred validation should reject unsafe extra valence."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_overbonded.pdb", include_conect=False)
    _patch_coordinate_bonds(
        monkeypatch,
        atom_count=9,
        bonds=((1, 2), (1, 3), (1, 8), (3, 4), (4, 6), (6, 5), (6, 7), (8, 9), (3, 8)),
    )

    with pytest.raises(ValueError, match="overbonded atoms"):
        load_glygen_glycan_pdb(glycan_path)


def test_coordinate_inferred_loader_accepts_glygen_branch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Coordinate-inferred validation should accept a plausible branched glycan graph."""
    glycan_path = _write_glygen_branch_fixture(tmp_path / "glycan_branch.pdb")
    _patch_coordinate_bonds(
        monkeypatch,
        atom_count=13,
        bonds=(
            (1, 2),
            (1, 3),
            (1, 8),
            (3, 4),
            (3, 10),
            (4, 6),
            (6, 5),
            (6, 7),
            (8, 9),
            (10, 11),
            (11, 12),
            (12, 13),
        ),
    )

    result = load_glygen_glycan_pdb(glycan_path)

    assert result.connectivity_provenance == "coordinate_inferred"
    assert [residue.residue_name for residue in result.fragment.residues] == [
        "NAG",
        "MAN",
        "ROH",
        "GAL",
    ]


def test_input_path_is_n_glycosylation_only() -> None:
    """Input-path moieties should be rejected outside the N-glycosylation MVP."""
    moiety = SimpleNamespace(
        input_path=Path("glycan.pdb"), smiles=None, residue_name=None, polymer_recipe=None
    )

    with pytest.raises(ValueError, match="n_glycosylation"):
        validate_moiety_source_config(moiety, mechanism_name="nhs_lys")

    assert validate_moiety_source_config(moiety, mechanism_name="n_glycosylation") == ["input_path"]


def test_coordinate_only_workflow_removes_roh_and_links_asn(tmp_path: Path) -> None:
    """Coordinate-only workflow should write residue-preserved Asn-glycan artifacts."""
    glycan_path = _write_glygen_fixture(tmp_path / "glycan_conect.pdb", include_conect=True)
    protein_path = _write_asn_fixture(tmp_path / "asn.pdb")
    attachment = _attachment(glycan_path)
    settings = ConjugatedPolymerSystemSettings(glygen_pdb_output_mode="coordinate_only")
    spec, *_ = _build_attachment_spec(
        attachment,
        attachment_index=1,
        protein_pdb_path=protein_path,
        artifact_dir=tmp_path,
        workflow_settings=settings,
    )

    result = _build_glygen_coordinate_only_result(
        protein_pdb_path=protein_path,
        specs=(spec,),
        output_dir=tmp_path,
        construction_dir=tmp_path / "construction",
        protein_canonicalization=None,
    )
    output = result.crosslinked_conjugate_pdb_path.read_text(encoding="utf-8")

    assert result.status == "coordinate_only"
    assert " ROH " not in output
    assert " NAG C" in output
    assert " MAN C" in output
    assert "CONECT" in output
    assert " ASX A" in output
    assert "HD21" not in output
    assert Path(result.artifact_paths["glygen_glygen_ingestion"]).exists()


def _attachment(glycan_path: Path) -> SimpleNamespace:
    """Build a minimal n_glycosylation attachment namespace."""
    return SimpleNamespace(
        name="asn_glycan",
        enabled=True,
        site=SimpleNamespace(
            chain_id="A",
            residue_name="ASN",
            residue_number=1,
            atom_name="ND2",
            insertion_code="",
            atom_serial=None,
            atom_index=None,
        ),
        moiety=SimpleNamespace(
            name="glycan",
            input_path=glycan_path,
            smiles=None,
            residue_name=None,
            polymer_recipe=None,
        ),
        mechanism=SimpleNamespace(
            name="n_glycosylation",
            product_residues=SimpleNamespace(site=None, moiety=None),
            bond=SimpleNamespace(order=1, site_atom=None, target_bond_length_angstrom=1.45),
        ),
    )


def _write_asn_fixture(path: Path) -> Path:
    """Write a tiny explicit-hydrogen ASN residue fixture."""
    atoms = [
        (1, "N", "ASN", "A", 1, 0.0, 0.0, 0.0, "N"),
        (2, "CA", "ASN", "A", 1, 1.4, 0.0, 0.0, "C"),
        (3, "C", "ASN", "A", 1, 2.0, 1.2, 0.0, "C"),
        (4, "CB", "ASN", "A", 1, 1.8, -1.2, 0.0, "C"),
        (5, "CG", "ASN", "A", 1, 3.0, -1.2, 0.0, "C"),
        (6, "OD1", "ASN", "A", 1, 3.6, -2.2, 0.0, "O"),
        (7, "ND2", "ASN", "A", 1, 3.5, 0.0, 0.0, "N"),
        (8, "HD21", "ASN", "A", 1, 4.4, 0.0, 0.0, "H"),
        (9, "HD22", "ASN", "A", 1, 3.1, 0.8, 0.0, "H"),
    ]
    path.write_text("".join(_pdb_atom(*atom) for atom in atoms) + "END\n", encoding="utf-8")
    return path


def _write_glygen_fixture(path: Path, *, include_conect: bool) -> Path:
    """Write a small labeled multi-residue GlyGen-style glycan fixture."""
    atoms = [
        (1, "C1", "NAG", "", 1, 0.0, 0.0, 0.0, "C"),
        (2, "O5", "NAG", "", 1, 0.0, 1.4, 0.0, "O"),
        (3, "C2", "NAG", "", 1, 1.5, 0.0, 0.0, "C"),
        (4, "O2", "NAG", "", 1, 2.8, 0.0, 0.0, "O"),
        (5, "O1", "MAN", "", 2, 5.2, 0.0, 0.0, "O"),
        (6, "C1", "MAN", "", 2, 4.1, 0.0, 0.0, "C"),
        (7, "O5", "MAN", "", 2, 4.1, 1.4, 0.0, "O"),
        (8, "O1", "ROH", "", 3, -1.3, 0.0, 0.0, "O"),
        (9, "HO1", "ROH", "", 3, -2.2, 0.0, 0.0, "H"),
    ]
    lines = [_pdb_atom(*atom) for atom in atoms]
    if include_conect:
        lines.extend(
            [
                "CONECT    1    2    3    8\n",
                "CONECT    3    4\n",
                "CONECT    4    6\n",
                "CONECT    6    5\n",
                "CONECT    6    7\n",
                "CONECT    8    9\n",
            ]
        )
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def _write_glygen_branch_fixture(path: Path) -> Path:
    """Write a labeled branched GlyGen-style glycan fixture without CONECT."""
    atoms = [
        (1, "C1", "NAG", "", 1, 0.0, 0.0, 0.0, "C"),
        (2, "O5", "NAG", "", 1, 0.0, 1.4, 0.0, "O"),
        (3, "C2", "NAG", "", 1, 1.5, 0.0, 0.0, "C"),
        (4, "O2", "NAG", "", 1, 2.8, 0.0, 0.0, "O"),
        (5, "O1", "MAN", "", 2, 5.2, 0.0, 0.0, "O"),
        (6, "C1", "MAN", "", 2, 4.1, 0.0, 0.0, "C"),
        (7, "O5", "MAN", "", 2, 4.1, 1.4, 0.0, "O"),
        (8, "O1", "ROH", "", 3, -1.3, 0.0, 0.0, "O"),
        (9, "HO1", "ROH", "", 3, -2.2, 0.0, 0.0, "H"),
        (10, "C3", "NAG", "", 1, 1.5, 1.4, 0.0, "C"),
        (11, "O3", "NAG", "", 1, 2.8, 1.4, 0.0, "O"),
        (12, "C1", "GAL", "", 4, 3.9, 1.4, 0.0, "C"),
        (13, "O5", "GAL", "", 4, 3.9, 2.8, 0.0, "O"),
    ]
    path.write_text("".join(_pdb_atom(*atom) for atom in atoms) + "END\n", encoding="utf-8")
    return path


def _pdb_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    """Format one fixed-width PDB atom line for fixtures."""
    return (
        f"HETATM{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )


def _replace_atom_element(path: Path, *, serial: int, element: str) -> None:
    """Replace one fixture atom element while preserving fixed-width columns."""
    lines = []
    for line in path.read_text(encoding="utf-8").splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")) and int(line[6:11]) == serial:
            line = f"{line[:76]}{element:>2}{line[78:]}"
        lines.append(line)
    path.write_text("".join(lines), encoding="utf-8")


def _patch_coordinate_bonds(
    monkeypatch: pytest.MonkeyPatch,
    *,
    atom_count: int,
    bonds: tuple[tuple[int, int], ...],
) -> None:
    """Patch RDKit proximity bonding to make graph validation deterministic."""
    from rdkit import Chem

    class FakeBond:
        """Minimal RDKit bond double for serial-order fixtures."""

        def __init__(self, begin_serial: int, end_serial: int) -> None:
            self._begin_index = begin_serial - 1
            self._end_index = end_serial - 1

        def GetBeginAtomIdx(self) -> int:
            """Return the zero-based begin atom index."""
            return self._begin_index

        def GetEndAtomIdx(self) -> int:
            """Return the zero-based end atom index."""
            return self._end_index

    class FakeMol:
        """Minimal RDKit molecule double for coordinate inference."""

        def GetNumAtoms(self) -> int:
            """Return the fixture atom count."""
            return atom_count

        def GetBonds(self) -> tuple[FakeBond, ...]:
            """Return deterministic fake bonds."""
            return tuple(FakeBond(left, right) for left, right in bonds)

    def fake_mol_from_pdb_file(*_args: Any, **_kwargs: Any) -> FakeMol:
        """Return the fake molecule for the loader's RDKit call."""
        return FakeMol()

    monkeypatch.setattr(Chem, "MolFromPDBFile", fake_mol_from_pdb_file)
