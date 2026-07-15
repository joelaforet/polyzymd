"""Tests for Packmol-constrained modifier placement."""

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

from polyzymd.builders.conjugation._linkage import (
    ExplicitLinkageContract,
    LinkageBond,
    NhsLysModifierLinker,
    PdbAtomSelector,
    ReactiveEndpoint,
    resolve_explicit_linkage_contract,
)
from polyzymd.builders.conjugation.placement import (
    _minimum_distance,
    place_modifier_with_packmol,
    place_modifier_with_resolved_plan,
    place_modifiers_with_resolved_plans,
)
from polyzymd.builders.conjugation.polymer import GeneratedPolymerFragment
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord


def test_packmol_input_uses_random_constrained_reactive_placement(monkeypatch, tmp_path: Path):
    """Packmol input should use random constrained placement without deterministic orientation."""
    protein_path = _protein_pdb(tmp_path)
    modifier = _generated_modifier()
    linker = NhsLysModifierLinker(target_residue_number=23)
    monkeypatch.chdir(tmp_path)

    def fake_run_packmol(input_text: str, work_dir: Path) -> Path:
        """Write a simple Packmol-like output preserving input ordering."""
        assert work_dir.is_absolute()
        for line in input_text.splitlines():
            if line.startswith("structure "):
                structure_path = Path(line.split(maxsplit=1)[1])
                assert structure_path.is_absolute()
                assert structure_path.exists()
        output_path = work_dir / "packmol_output.pdb"
        protein_lines = [
            line
            for line in (work_dir / "protein_fixed_sterics.pdb").read_text().splitlines(True)
            if line.startswith(("ATOM", "HETATM"))
        ]
        modifier_lines = [
            line
            for line in (work_dir / "modifier_retained.pdb").read_text().splitlines(True)
            if line.startswith(("ATOM", "HETATM"))
        ]
        output_path.write_text("".join([*protein_lines, *modifier_lines, "END\n"]))
        assert input_text == (work_dir / "packmol.inp").read_text()
        return output_path

    result = place_modifier_with_packmol(
        protein_path,
        modifier,
        linker,
        Path("relative-placement-artifacts"),
        run_packmol_func=fake_run_packmol,
    )

    input_text = result.packmol_input_text
    assert "movebadrandom" in input_text
    assert "nloop 500" in input_text
    assert "fixed 0. 0. 0. 0. 0. 0." in input_text
    assert "structure" in input_text
    assert "atoms 2" in input_text
    assert input_text.count("inside sphere") == 1
    assert "inside sphere" in input_text
    assert "outside sphere" not in input_text
    assert "rotate" not in input_text.lower()
    assert "center" not in input_text.lower()
    assert "NZ" in result.excluded_protein_atom_names
    assert "HZ2" in result.excluded_protein_atom_names
    assert "HZ3" in result.excluded_protein_atom_names
    assert abs(result.placed_bond_length_angstrom - 1.33) < 1.0e-6
    assert result.placed_modifier.reactive_atom_name == "RC"


def test_resolved_plan_placement_uses_resolved_atoms_and_target_length(tmp_path: Path):
    """Placement should use resolved atom identities and contract bond length."""
    protein_path = _protein_pdb(tmp_path)
    modifier = _generated_modifier()
    plan = resolve_explicit_linkage_contract(
        protein_path,
        modifier,
        _explicit_contract(target_bond_length=1.45),
    )

    def fake_run_packmol(input_text: str, work_dir: Path) -> Path:
        """Write a simple Packmol-like output preserving input ordering."""
        output_path = work_dir / "packmol_output.pdb"
        protein_lines = [
            line
            for line in (work_dir / "protein_fixed_sterics.pdb").read_text().splitlines(True)
            if line.startswith(("ATOM", "HETATM"))
        ]
        modifier_lines = [
            line
            for line in (work_dir / "modifier_retained.pdb").read_text().splitlines(True)
            if line.startswith(("ATOM", "HETATM"))
        ]
        output_path.write_text("".join([*protein_lines, *modifier_lines, "END\n"]))
        assert "atoms 2" in input_text
        return output_path

    result = place_modifier_with_resolved_plan(
        protein_path,
        modifier,
        plan,
        tmp_path,
        run_packmol_func=fake_run_packmol,
    )

    assert result.target_bond_length_angstrom == 1.45
    assert abs(result.placed_bond_length_angstrom - 1.45) < 1.0e-6
    retained_lines = [
        line
        for line in result.modifier_pdb_path.read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]
    assert len(retained_lines) == 3


def test_glygen_coordinate_only_packmol_input_constrains_only_c1(tmp_path: Path):
    """GlyGen coordinate-only placement should constrain only reactive C1."""
    protein_path = _protein_pdb(tmp_path)
    modifier = _generated_modifier_with_c1_reactive_atom()
    plan = resolve_explicit_linkage_contract(
        protein_path,
        modifier,
        _explicit_contract(target_bond_length=1.45, modifier_atom_name="C1"),
    )

    def fake_run_packmol(input_text: str, work_dir: Path) -> Path:
        """Write a simple Packmol-like output preserving input ordering."""
        output_path = work_dir / "packmol_output.pdb"
        protein_lines = [
            line
            for line in (work_dir / "protein_fixed_sterics.pdb").read_text().splitlines(True)
            if line.startswith(("ATOM", "HETATM"))
        ]
        modifier_lines = [
            line
            for line in (work_dir / "modifier_retained.pdb").read_text().splitlines(True)
            if line.startswith(("ATOM", "HETATM"))
        ]
        output_path.write_text("".join([*protein_lines, *modifier_lines, "END\n"]))
        return output_path

    result = place_modifier_with_resolved_plan(
        protein_path,
        modifier,
        plan,
        tmp_path,
        run_packmol_func=fake_run_packmol,
    )

    input_text = result.packmol_input_text
    assert "atoms 1\n  inside sphere" in input_text
    assert input_text.count("inside sphere") == 1
    assert "inside sphere" in input_text
    assert " 1.40" not in input_text
    assert " 1.50" not in input_text
    assert "outside sphere" not in input_text
    assert "atoms 2\n  outside sphere" not in input_text
    assert "atoms 3\n  outside sphere" not in input_text


def test_joint_resolved_plan_placement_uses_one_packmol_run_for_two_fragments(
    monkeypatch,
    tmp_path: Path,
):
    """Joint placement should build one Packmol input with all constrained fragments."""
    protein_path = _protein_pdb(tmp_path)
    modifiers = (_generated_modifier(), _generated_modifier())
    plans = tuple(
        resolve_explicit_linkage_contract(
            protein_path,
            modifier,
            _explicit_contract(target_bond_length=1.45),
        )
        for modifier in modifiers
    )
    calls = {"packmol": 0}
    monkeypatch.chdir(tmp_path)

    def fake_run_packmol(input_text: str, work_dir: Path) -> Path:
        """Write a simple Packmol-like output preserving structure ordering."""
        calls["packmol"] += 1
        assert work_dir.is_absolute()
        output_path = work_dir / "packmol_output.pdb"
        atom_lines = []
        for line in input_text.splitlines():
            if not line.startswith("structure "):
                continue
            structure_path = Path(line.split(maxsplit=1)[1])
            assert structure_path.is_absolute()
            assert structure_path.exists()
            atom_lines.extend(
                pdb_line
                for pdb_line in structure_path.read_text().splitlines(True)
                if pdb_line.startswith(("ATOM", "HETATM"))
            )
        output_path.write_text("".join([*atom_lines, "END\n"]))
        assert input_text == (work_dir / "packmol.inp").read_text()
        return output_path

    results = place_modifiers_with_resolved_plans(
        protein_path,
        modifiers,
        plans,
        Path("relative-joint-placement-artifacts"),
        run_packmol_func=fake_run_packmol,
    )

    assert calls["packmol"] == 1
    assert len(results) == 2
    assert results[0].packmol_input_path == results[1].packmol_input_path
    assert results[0].packmol_output_path == results[1].packmol_output_path
    assert results[0].modifier_pdb_path.parent.parent.name == "placement_01"
    assert results[1].modifier_pdb_path.parent.parent.name == "placement_02"
    input_text = results[0].packmol_input_text
    assert input_text.count("structure ") == 3
    assert input_text.count("inside sphere") == 2
    assert "outside sphere" not in input_text
    assert all(abs(result.placed_bond_length_angstrom - 1.45) < 1.0e-6 for result in results)


def test_minimum_distance_uses_mdanalysis_distance_array(monkeypatch):
    """Minimum distance should use MDAnalysis distance utilities when available."""
    points_a = np.asarray([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    points_b = np.asarray([[2.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    calls = {"distance_array": 0}

    def fake_distance_array(observed_a: np.ndarray, observed_b: np.ndarray) -> np.ndarray:
        """Return NumPy-equivalent distances for the fake MDAnalysis module."""
        calls["distance_array"] += 1
        np.testing.assert_array_equal(observed_a, points_a)
        np.testing.assert_array_equal(observed_b, points_b)
        return np.linalg.norm(observed_a[:, np.newaxis, :] - observed_b[np.newaxis, :, :], axis=2)

    _install_fake_mdanalysis_distance_array(monkeypatch, fake_distance_array)

    assert _minimum_distance(points_a, points_b) == 2.0
    assert calls["distance_array"] == 1


def test_minimum_distance_falls_back_when_mdanalysis_distance_array_errors(monkeypatch):
    """Minimum distance should preserve NumPy behavior when MDAnalysis cannot run."""
    points_a = np.asarray([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    points_b = np.asarray([[2.0, 0.0, 0.0], [10.0, 0.0, 0.0]])

    def fake_distance_array(_points_a: np.ndarray, _points_b: np.ndarray) -> np.ndarray:
        """Simulate an expected MDAnalysis runtime failure."""
        raise RuntimeError("MDAnalysis distance backend unavailable")

    _install_fake_mdanalysis_distance_array(monkeypatch, fake_distance_array)

    assert _minimum_distance(points_a, points_b) == 2.0


def test_minimum_distance_propagates_unexpected_mdanalysis_runtime_errors(monkeypatch):
    """Unexpected MDAnalysis runtime errors should not silently use NumPy fallback."""
    points_a = np.asarray([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    points_b = np.asarray([[2.0, 0.0, 0.0], [10.0, 0.0, 0.0]])

    def fake_distance_array(_points_a: np.ndarray, _points_b: np.ndarray) -> np.ndarray:
        """Simulate an unexpected MDAnalysis runtime failure."""
        raise RuntimeError("unexpected coordinate shape regression")

    _install_fake_mdanalysis_distance_array(monkeypatch, fake_distance_array)

    with pytest.raises(RuntimeError, match="unexpected coordinate shape regression"):
        _minimum_distance(points_a, points_b)


@pytest.mark.parametrize(
    "message",
    (
        "failed to import coordinate shape regression",
        "could not import malformed coordinate buffer",
        "no module named coordinate buffer",
    ),
)
def test_minimum_distance_propagates_generic_import_wording_runtime_errors(
    monkeypatch,
    message: str,
):
    """Generic import-wording RuntimeErrors should not trigger backend fallback."""
    points_a = np.asarray([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    points_b = np.asarray([[2.0, 0.0, 0.0], [10.0, 0.0, 0.0]])

    def fake_distance_array(_points_a: np.ndarray, _points_b: np.ndarray) -> np.ndarray:
        """Simulate a non-backend RuntimeError with import-like wording."""
        raise RuntimeError(message)

    _install_fake_mdanalysis_distance_array(monkeypatch, fake_distance_array)

    with pytest.raises(RuntimeError, match=message):
        _minimum_distance(points_a, points_b)


def _install_fake_mdanalysis_distance_array(monkeypatch, distance_array) -> None:
    """Install a fake MDAnalysis distance module for placement helper tests."""
    mdanalysis_module = ModuleType("MDAnalysis")
    lib_module = ModuleType("MDAnalysis.lib")
    distances_module = ModuleType("MDAnalysis.lib.distances")
    distances_module.distance_array = distance_array
    lib_module.distances = distances_module
    mdanalysis_module.lib = lib_module

    monkeypatch.setitem(sys.modules, "MDAnalysis", mdanalysis_module)
    monkeypatch.setitem(sys.modules, "MDAnalysis.lib", lib_module)
    monkeypatch.setitem(sys.modules, "MDAnalysis.lib.distances", distances_module)


def _explicit_contract(
    *,
    target_bond_length: float,
    modifier_atom_name: str = "RC",
) -> ExplicitLinkageContract:
    """Build a generic explicit linkage contract for placement tests."""
    modifier_residue_number = 1 if modifier_atom_name == "C1" else 2
    modifier_leaving_atom_names = () if modifier_atom_name == "C1" else ("LG",)
    return ExplicitLinkageContract(
        protein_endpoint=ReactiveEndpoint(
            participant="protein",
            selector=PdbAtomSelector(
                chain_id="A",
                residue_name="LYS",
                residue_number=23,
                atom_name="NZ",
            ),
            product_residue_name="LYX",
            leaving_atom_names=("HZ2", "HZ3"),
        ),
        modifier_endpoint=ReactiveEndpoint(
            participant="modifier",
            selector=PdbAtomSelector(
                chain_id="Z",
                residue_name="SB1" if modifier_atom_name == "C1" else "NHS",
                residue_number=modifier_residue_number,
                atom_name=modifier_atom_name,
            ),
            product_residue_name="NHX",
            leaving_atom_names=modifier_leaving_atom_names,
        ),
        bond=LinkageBond(
            protein_atom_name="NZ",
            modifier_atom_name=modifier_atom_name,
            bond_order=1,
            target_bond_length_angstrom=target_bond_length,
        ),
        mechanism_name="explicit_linkage",
    )


def _generated_modifier_with_c1_reactive_atom() -> GeneratedPolymerFragment:
    """Create a GlyGen-like coordinate-only modifier with reactive C1."""
    fragment = _generated_modifier()
    return fragment.model_copy(
        update={
            "reactive_atom_serial": 101,
            "reactive_atom_index": 0,
            "reactive_atom_name": "C1",
        }
    )


def _protein_pdb(tmp_path: Path) -> Path:
    """Create a small lysine-containing protein PDB."""
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
        + "END\n"
    )
    return path


def _generated_modifier() -> GeneratedPolymerFragment:
    """Create a small generated modifier with a leaving group."""
    atoms = (
        PdbAtomRecord(
            serial=101,
            atom_index=0,
            atom_name="C1",
            residue_name="SB1",
            chain_id="Z",
            residue_number=1,
            x=5.0,
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
            chain_id="Z",
            residue_number=2,
            x=3.3,
            y=0.0,
            z=0.0,
            element="C",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=103,
            atom_index=2,
            atom_name="O1",
            residue_name="NHS",
            chain_id="Z",
            residue_number=2,
            x=3.8,
            y=0.5,
            z=0.0,
            element="O",
            record_name="HETATM",
        ),
        PdbAtomRecord(
            serial=104,
            atom_index=3,
            atom_name="LG",
            residue_name="NHS",
            chain_id="Z",
            residue_number=2,
            x=4.2,
            y=1.0,
            z=0.0,
            element="O",
            record_name="HETATM",
        ),
    )
    return GeneratedPolymerFragment.from_atom_records(
        atoms,
        bonds=((101, 102), (102, 103), (102, 104)),
        reactive_atom_serial=102,
        reactive_atom_index=1,
        reactive_atom_name="RC",
        leaving_atom_serials=(104,),
        leaving_atom_indices=(3,),
        leaving_atom_names=("LG",),
        name="modifier",
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
) -> str:
    """Format one PDB atom line for tests."""
    return (
        f"ATOM  {serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    {x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )
