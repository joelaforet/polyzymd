"""Tests for GlyGen/GlyCAM residue-resolved glycan PDB ingestion."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

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
        (5, "O1", "MAN", "", 2, 3.9, 0.0, 0.0, "O"),
        (6, "C1", "MAN", "", 2, 5.2, 0.0, 0.0, "C"),
        (7, "O5", "MAN", "", 2, 5.2, 1.4, 0.0, "O"),
        (8, "O1", "ROH", "", 3, -1.3, 0.0, 0.0, "O"),
        (9, "HO1", "ROH", "", 3, -2.2, 0.0, 0.0, "H"),
    ]
    lines = [_pdb_atom(*atom) for atom in atoms]
    if include_conect:
        lines.extend(
            [
                "CONECT    1    2    3    8\n",
                "CONECT    3    4\n",
                "CONECT    4    5\n",
                "CONECT    5    6\n",
                "CONECT    6    7\n",
                "CONECT    8    9\n",
            ]
        )
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
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
