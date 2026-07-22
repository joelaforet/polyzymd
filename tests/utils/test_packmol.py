"""
Tests for polyzymd.utils.packmol — input-file generation utilities.

These tests exercise :func:`build_packmol_input` directly and do NOT
require a Packmol binary or any heavy simulation dependencies.
"""

from __future__ import annotations

import json
import subprocess
import sys
import types
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

from polyzymd.utils.packmol import (
    _PACKMOL_OUTPUT_FILE,
    build_packmol_input,
)


class _CompletedPackmol:
    """Minimal subprocess result for mocked Packmol execution."""

    returncode = 0
    stdout = b"Success!\n"


class _TimedOutPackmol:
    """Minimal Popen process that times out, then terminates cleanly."""

    returncode = -15

    def __init__(self, *, require_kill: bool = False) -> None:
        self.communicate_calls = 0
        self.terminated = False
        self.killed = False
        self.require_kill = require_kill

    def communicate(self, timeout=None):
        """Raise immediately on the bounded run and return after terminate."""
        self.communicate_calls += 1
        if self.communicate_calls == 1 or (self.communicate_calls == 2 and self.require_kill):
            raise subprocess.TimeoutExpired("packmol", timeout)
        return b"Packmol partial progress\n", None

    def terminate(self) -> None:
        """Record graceful termination."""
        self.terminated = True

    def kill(self) -> None:
        """Record forced termination."""
        self.killed = True


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

BOX_3A = np.array([30.0, 40.0, 50.0])  # Angstrom box dimensions
PDBS = ["mol0.pdb", "mol1.pdb"]
COUNTS = [5, 10]
TOL = 2.0


# ---------------------------------------------------------------------------
# Header structure tests
# ---------------------------------------------------------------------------


class TestBuildPackmolInputHeader:
    """The generated input must contain the mandatory packmol header lines."""

    def test_tolerance_line_present(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL)
        assert f"tolerance {TOL:f}" in text

    def test_filetype_pdb_line_present(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL)
        assert "filetype pdb" in text

    def test_output_filename_present(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL)
        assert f"output {_PACKMOL_OUTPUT_FILE}" in text


# ---------------------------------------------------------------------------
# movebadrandom keyword tests
# ---------------------------------------------------------------------------


class TestMovebadrandom:
    """The movebadrandom keyword should appear iff the flag is True."""

    def test_movebadrandom_absent_by_default(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL)
        assert "movebadrandom" not in text

    def test_movebadrandom_absent_when_false(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL, movebadrandom=False)
        assert "movebadrandom" not in text

    def test_movebadrandom_present_when_true(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL, movebadrandom=True)
        assert "movebadrandom" in text

    def test_movebadrandom_precedes_structure_blocks(self):
        """movebadrandom must appear before any structure block."""
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL, movebadrandom=True)
        mbr_pos = text.index("movebadrandom")
        struct_pos = text.index("structure")
        assert mbr_pos < struct_pos


# ---------------------------------------------------------------------------
# Solute (fixed) block tests
# ---------------------------------------------------------------------------


class TestSoluteBlock:
    """When a solute PDB path is given a fixed structure block should appear."""

    def test_solute_block_present(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL, solute_pdb_path="protein.pdb")
        assert "structure protein.pdb" in text
        assert "number 1" in text
        assert "fixed 0. 0. 0. 0. 0. 0." in text

    def test_solute_block_absent_when_none(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL, solute_pdb_path=None)
        assert "fixed" not in text

    def test_solute_precedes_molecule_blocks(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL, solute_pdb_path="protein.pdb")
        solute_pos = text.index("structure protein.pdb")
        mol0_pos = text.index("structure mol0.pdb")
        assert solute_pos < mol0_pos


# ---------------------------------------------------------------------------
# Molecule structure block tests
# ---------------------------------------------------------------------------


class TestMoleculeBlocks:
    """One structure block per molecule type, skipping zero-count entries."""

    def test_all_molecule_blocks_present(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL)
        for pdb in PDBS:
            assert f"structure {pdb}" in text

    def test_molecule_counts_in_blocks(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL)
        for count in COUNTS:
            assert f"number {count}" in text

    def test_zero_count_molecule_omitted(self):
        pdbs = ["mol0.pdb", "mol1.pdb", "mol2.pdb"]
        counts = [3, 0, 7]
        text = build_packmol_input(pdbs, counts, BOX_3A, TOL)
        assert "structure mol0.pdb" in text
        assert "structure mol1.pdb" not in text
        assert "structure mol2.pdb" in text

    def test_inside_box_line_in_non_pbc_mode(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL, use_pbc=False)
        assert "inside box" in text

    def test_inside_box_absent_in_pbc_mode(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL, use_pbc=True)
        assert "inside box" not in text

    def test_pbc_keyword_present_in_pbc_mode(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL, use_pbc=True)
        assert text.startswith("tolerance") or "pbc " in text
        assert "pbc " in text


# ---------------------------------------------------------------------------
# Box-size arithmetic tests
# ---------------------------------------------------------------------------


class TestBoxSizeArithmetic:
    """The effective box used in 'inside box' must be shrunk by tolerance."""

    def test_inside_box_shrunk_by_tolerance(self):
        box = np.array([30.0, 40.0, 50.0])
        tol = 2.0
        text = build_packmol_input(["m.pdb"], [1], box, tol, use_pbc=False)
        # Effective max coords = box - tol
        assert f"{28.0:.6f}" in text  # 30 - 2
        assert f"{38.0:.6f}" in text  # 40 - 2
        assert f"{48.0:.6f}" in text  # 50 - 2

    def test_pbc_box_not_shrunk(self):
        box = np.array([30.0, 40.0, 50.0])
        tol = 2.0
        text = build_packmol_input(["m.pdb"], [1], box, tol, use_pbc=True)
        assert f"{30.0:.6f}" in text
        assert f"{40.0:.6f}" in text
        assert f"{50.0:.6f}" in text


# ---------------------------------------------------------------------------
# structure_extra_lines tests
# ---------------------------------------------------------------------------


class TestStructureExtraLines:
    """Per-structure extra lines are inserted correctly and validated."""

    def test_extra_lines_absent_by_default(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL)
        assert "atoms 42" not in text
        assert "inside sphere 10.5 20.3 15.7 5.0" not in text
        assert "end atoms" not in text

    def test_extra_lines_inserted_before_end_structure(self):
        extra_lines = [["atoms 1", "end atoms"], []]
        text = build_packmol_input(
            PDBS,
            COUNTS,
            BOX_3A,
            TOL,
            structure_extra_lines=extra_lines,
        )

        block_start = text.index("structure mol0.pdb")
        block_end = text.index("end structure", block_start)
        atoms_pos = text.index("  atoms 1", block_start)
        end_atoms_pos = text.index("  end atoms", block_start)

        assert block_start < atoms_pos < end_atoms_pos < block_end

    def test_extra_lines_with_atoms_constraint(self):
        extra_lines = [["atoms 42", "inside sphere 10.5 20.3 15.7 5.0", "end atoms"], []]
        text = build_packmol_input(
            PDBS,
            COUNTS,
            BOX_3A,
            TOL,
            structure_extra_lines=extra_lines,
        )

        assert "  atoms 42" in text
        assert "  inside sphere 10.5 20.3 15.7 5.0" in text
        assert "  end atoms" in text

    def test_extra_lines_empty_list_no_effect(self):
        text_default = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL)
        text_empty = build_packmol_input(
            PDBS,
            COUNTS,
            BOX_3A,
            TOL,
            structure_extra_lines=[[], []],
        )
        assert text_empty == text_default

    def test_extra_lines_length_mismatch_raises(self):
        with pytest.raises(ValueError, match="structure_extra_lines"):
            build_packmol_input(
                PDBS,
                COUNTS,
                BOX_3A,
                TOL,
                structure_extra_lines=[[]],
            )

    def test_extra_lines_skips_zero_count(self):
        pdbs = ["mol0.pdb", "mol1.pdb"]
        counts = [0, 2]
        extra_lines = [["atoms 999", "end atoms"], ["atoms 42", "end atoms"]]
        text = build_packmol_input(
            pdbs,
            counts,
            BOX_3A,
            TOL,
            structure_extra_lines=extra_lines,
        )

        assert "structure mol0.pdb" not in text
        assert "  atoms 999" not in text
        assert "structure mol1.pdb" in text
        assert "  atoms 42" in text

    def test_extra_lines_multiple_structures(self):
        extra_lines = [
            ["atoms 11", "end atoms"],
            ["atoms 22", "end atoms"],
        ]
        text = build_packmol_input(
            PDBS,
            COUNTS,
            BOX_3A,
            TOL,
            structure_extra_lines=extra_lines,
        )

        mol0_start = text.index("structure mol0.pdb")
        mol1_start = text.index("structure mol1.pdb")

        mol0_block = text[mol0_start:mol1_start]
        mol1_block = text[mol1_start:]

        assert "  atoms 11" in mol0_block
        assert "  atoms 22" not in mol0_block
        assert "  atoms 22" in mol1_block
        assert "  atoms 11" not in mol1_block


# ---------------------------------------------------------------------------
# PolymerPackingConfig schema test
# ---------------------------------------------------------------------------


class TestPolymerPackingConfigSchema:
    """movebadrandom should be readable from the config schema."""

    def test_default_movebadrandom_is_false(self):
        from polyzymd.config.schema import PolymerPackingConfig

        cfg = PolymerPackingConfig()
        assert cfg.movebadrandom is False

    def test_movebadrandom_can_be_set_true(self):
        from polyzymd.config.schema import PolymerPackingConfig

        cfg = PolymerPackingConfig(movebadrandom=True)
        assert cfg.movebadrandom is True

    def test_default_padding_and_tolerance_unchanged(self):
        from polyzymd.config.schema import PolymerPackingConfig

        cfg = PolymerPackingConfig()
        assert cfg.padding == pytest.approx(2.0)
        assert cfg.tolerance == pytest.approx(2.0)


class TestRunPackmol:
    """Regression tests for Packmol executor path handling."""

    def test_relative_working_directory_is_resolved(self, monkeypatch, tmp_path):
        """run_packmol should reopen input after cwd changes for relative work dirs."""
        from polyzymd.utils import packmol

        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(packmol.shutil, "which", MagicMock(return_value="packmol"))

        def fake_run(_binary: str, *, stdin, stdout, stderr):
            """Validate that Packmol receives an existing absolute input file."""
            del stdout, stderr
            input_path = Path(stdin.name)
            assert input_path.is_absolute()
            assert input_path.exists()
            assert input_path.read_text(encoding="utf-8") == "Success input\n"
            return _CompletedPackmol()

        monkeypatch.setattr(packmol.subprocess, "run", fake_run)

        output_path = packmol.run_packmol("Success input\n", Path("relative_packmol_work"))

        assert output_path.is_absolute()
        assert output_path == (tmp_path / "relative_packmol_work" / _PACKMOL_OUTPUT_FILE).resolve()
        assert (tmp_path / "relative_packmol_work" / "packmol_input.txt").exists()

    @pytest.mark.parametrize("require_kill", [False, True])
    def test_timeout_preserves_candidate_and_writes_explicit_status(
        self, monkeypatch, tmp_path, caplog, require_kill
    ):
        """A timeout should stop Packmol and retain a candidate for validation."""
        from polyzymd.utils import packmol

        monkeypatch.setattr(packmol.shutil, "which", MagicMock(return_value="packmol"))
        process = _TimedOutPackmol(require_kill=require_kill)

        def fake_popen(_binary: str, *, stdin, stdout, stderr):
            """Create a complete best candidate before the simulated timeout."""
            del stdin, stdout, stderr
            Path(_PACKMOL_OUTPUT_FILE).write_text("END\n", encoding="utf-8")
            return process

        monkeypatch.setattr(packmol.subprocess, "Popen", fake_popen)

        output_path = packmol.run_packmol(
            "timeout input\n",
            tmp_path,
            timeout_seconds=0.01,
        )

        status = json.loads((tmp_path / "packmol_run_status.json").read_text())
        assert output_path == (tmp_path / _PACKMOL_OUTPUT_FILE).resolve()
        assert process.terminated is True
        assert process.killed is require_kill
        assert process.communicate_calls == (3 if require_kill else 2)
        assert status == {
            "accepted_after_validation": None,
            "outcome": "timeout_candidate",
            "output_exists": True,
            "output_path": str(output_path),
            "return_code": -15,
            "timed_out": True,
        }
        assert "mandatory structural validation" in caplog.text

    def test_exit_173_writes_explicit_candidate_status(self, monkeypatch, tmp_path, caplog):
        """Exit 173 should be retained as an imperfect candidate, not inferred from text."""
        from polyzymd.utils import packmol

        result = _CompletedPackmol()
        result.returncode = 173
        result.stdout = b"Packmol stopped without a perfect packing\n"
        monkeypatch.setattr(packmol.shutil, "which", MagicMock(return_value="packmol"))
        monkeypatch.setattr(packmol.subprocess, "run", MagicMock(return_value=result))

        output_path = packmol.run_packmol("imperfect input\n", tmp_path)
        status = json.loads((tmp_path / "packmol_run_status.json").read_text())

        assert output_path == (tmp_path / _PACKMOL_OUTPUT_FILE).resolve()
        assert status["timed_out"] is False
        assert status["return_code"] == 173
        assert status["outcome"] == "imperfect_candidate"
        assert status["accepted_after_validation"] is None
        assert "imperfect packing" in caplog.text


class _MockTopology:
    """Lightweight topology stand-in for assembly-frame tests."""

    def __init__(self, positions: np.ndarray):
        self.positions = np.asarray(positions, dtype=float)
        self.n_atoms = int(self.positions.shape[0])
        self.box_vectors = None

    def set_positions(self, positions: np.ndarray) -> None:
        """Store assigned positions for later assertions."""
        self.positions = np.asarray(positions, dtype=float)
        self.n_atoms = int(self.positions.shape[0])

    def __add__(self, other: _MockTopology) -> _MockTopology:
        """Combine two topologies by stacking coordinates."""
        merged = np.vstack([self.positions, other.positions])
        return _MockTopology(merged)


class _MockBrickSize:
    """Simple object exposing ``m_as`` like an OpenFF quantity."""

    def __init__(self, values: np.ndarray):
        self._values = np.asarray(values, dtype=float)

    def m_as(self, _unit: str) -> np.ndarray:
        """Return the stored raw values irrespective of unit string."""
        return self._values


class TestSolvateAssemblyCoordinates:
    """Regression tests for centered-solute topology assembly."""

    @staticmethod
    def _install_fake_openff_modules(
        monkeypatch: pytest.MonkeyPatch,
        *,
        centered_solute: _MockTopology,
        assembled_solvent: _MockTopology,
        loaded_positions: np.ndarray,
    ) -> dict[str, MagicMock]:
        """Install fake OpenFF modules used by lazy imports in packmol helpers."""
        mock_compute_brick = MagicMock(return_value=_MockBrickSize(np.array([30.0, 40.0, 50.0])))
        mock_center = MagicMock(return_value=centered_solute)
        mock_create_solute_pdb = MagicMock(return_value="solute.pdb")
        mock_create_molecule_pdbs = MagicMock(return_value=["water.pdb"])
        mock_load_positions = MagicMock(return_value=np.asarray(loaded_positions, dtype=float))

        packmol_mod = types.ModuleType("openff.packmol._packmol")
        packmol_mod._center_topology_at = mock_center
        packmol_mod._compute_brick_from_box_vectors = mock_compute_brick
        packmol_mod._create_molecule_pdbs = mock_create_molecule_pdbs
        packmol_mod._create_solute_pdb = mock_create_solute_pdb
        packmol_mod._load_positions = mock_load_positions

        openff_pkg = types.ModuleType("openff")
        openff_pkg.__path__ = []

        openff_packmol_pkg = types.ModuleType("openff.packmol")
        openff_packmol_pkg.__path__ = []

        toolkit_mod = types.ModuleType("openff.toolkit")

        class _FakeTopologyFactory:
            @staticmethod
            def from_molecules(_molecules: list[object]) -> _MockTopology:
                return assembled_solvent

        toolkit_mod.Topology = _FakeTopologyFactory

        units_mod = types.ModuleType("openff.units")
        units_mod.Quantity = lambda values, _unit: np.asarray(values, dtype=float)

        monkeypatch.setitem(sys.modules, "openff", openff_pkg)
        monkeypatch.setitem(sys.modules, "openff.packmol", openff_packmol_pkg)
        monkeypatch.setitem(sys.modules, "openff.packmol._packmol", packmol_mod)
        monkeypatch.setitem(sys.modules, "openff.toolkit", toolkit_mod)
        monkeypatch.setitem(sys.modules, "openff.units", units_mod)

        return {
            "center": mock_center,
            "compute_brick": mock_compute_brick,
            "create_solute_pdb": mock_create_solute_pdb,
            "create_molecule_pdbs": mock_create_molecule_pdbs,
            "load_positions": mock_load_positions,
        }

    def test_solvate_with_packmol_assembles_with_centered_solute(self, monkeypatch, tmp_path):
        """solvate_with_packmol should assemble using BRICK-centered solute."""
        from polyzymd.utils import packmol

        original_solute = _MockTopology(
            np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
        )
        centered_solute = _MockTopology(
            np.array([[15.0, 20.0, 25.0], [16.0, 21.0, 26.0], [17.0, 22.0, 27.0]])
        )
        solvent_topology = _MockTopology(np.zeros((3, 3), dtype=float))

        loaded_positions = np.array(
            [
                [15.0, 20.0, 25.0],
                [16.0, 21.0, 26.0],
                [17.0, 22.0, 27.0],
                [30.0, 30.0, 30.0],
                [31.0, 31.0, 31.0],
                [32.0, 32.0, 32.0],
            ]
        )

        openff_mocks = self._install_fake_openff_modules(
            monkeypatch,
            centered_solute=centered_solute,
            assembled_solvent=solvent_topology,
            loaded_positions=loaded_positions,
        )

        monkeypatch.setattr(
            packmol, "run_packmol", MagicMock(return_value=tmp_path / "packmol.pdb")
        )
        monkeypatch.setattr(packmol, "_strip_conect_records", MagicMock(return_value=0))
        monkeypatch.setattr(packmol, "_check_pbc_available", MagicMock(return_value=False))
        monkeypatch.setattr(
            packmol, "_check_ignore_conect_supported", MagicMock(return_value=False)
        )

        result = packmol.solvate_with_packmol(
            molecules=[object()],
            number_of_copies=[1],
            solute=original_solute,
            box_vectors=MagicMock(),
            working_directory=tmp_path,
        )

        openff_mocks["center"].assert_called_once()
        np.testing.assert_allclose(
            result.positions[: original_solute.n_atoms], centered_solute.positions
        )
        assert not np.allclose(
            result.positions[: original_solute.n_atoms], original_solute.positions
        )

    def test_pack_polymers_assembles_with_centered_solute(self, monkeypatch, tmp_path):
        """pack_polymers should assemble using BRICK-centered solute."""
        from polyzymd.utils import packmol

        original_solute = _MockTopology(
            np.array([[0.0, 0.0, 0.0], [2.0, 2.0, 2.0], [4.0, 4.0, 4.0]])
        )
        centered_solute = _MockTopology(
            np.array([[15.0, 20.0, 25.0], [18.0, 23.0, 28.0], [21.0, 26.0, 31.0]])
        )
        polymer_topology = _MockTopology(np.zeros((3, 3), dtype=float))

        loaded_positions = np.array(
            [
                [15.0, 20.0, 25.0],
                [18.0, 23.0, 28.0],
                [21.0, 26.0, 31.0],
                [40.0, 40.0, 40.0],
                [41.0, 41.0, 41.0],
                [42.0, 42.0, 42.0],
            ]
        )

        openff_mocks = self._install_fake_openff_modules(
            monkeypatch,
            centered_solute=centered_solute,
            assembled_solvent=polymer_topology,
            loaded_positions=loaded_positions,
        )

        boxvectors_mod = types.ModuleType("polyzymd.utils.boxvectors")
        boxvectors_mod.get_topology_bbox_bounds = MagicMock(
            return_value=(np.array([10.0, 10.0, 10.0]), np.array([20.0, 20.0, 20.0]))
        )
        monkeypatch.setitem(sys.modules, "polyzymd.utils.boxvectors", boxvectors_mod)

        monkeypatch.setattr(
            packmol, "run_packmol", MagicMock(return_value=tmp_path / "packmol.pdb")
        )
        monkeypatch.setattr(packmol, "_max_molecule_diameter_angstrom", MagicMock(return_value=1.0))

        result = packmol.pack_polymers(
            molecules=[object()],
            number_of_copies=[1],
            solute=original_solute,
            box_vectors=MagicMock(),
            working_directory=tmp_path,
        )

        openff_mocks["center"].assert_called_once()
        np.testing.assert_allclose(
            result.positions[: original_solute.n_atoms], centered_solute.positions
        )
        assert not np.allclose(
            result.positions[: original_solute.n_atoms], original_solute.positions
        )
