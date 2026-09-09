"""
Tests for polyzymd.utils.packmol — input-file generation utilities.

These tests exercise :func:`build_packmol_input` directly and do NOT
require a Packmol binary or any heavy simulation dependencies.
"""

from __future__ import annotations

import sys
import types
from unittest.mock import MagicMock

import numpy as np
import pytest

from polyzymd.utils.packmol import (
    _PACKMOL_OUTPUT_FILE,
    build_packmol_input,
)

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

    def get_positions(self) -> _MockBrickSize:
        """Return positions wrapped in an ``m_as``-capable object."""
        return _MockBrickSize(self.positions)


class _MockBrickSize:
    """Simple object exposing ``m_as`` like an OpenFF quantity."""

    def __init__(self, values: np.ndarray):
        self._values = np.asarray(values, dtype=float)

    def m_as(self, _unit: str) -> np.ndarray:
        """Return the stored raw values irrespective of unit string."""
        return self._values


def _mock_box_vectors(dims=(30.0, 40.0, 50.0)) -> _MockBrickSize:
    """Quantity-like 3x3 box vectors matching the mocked brick size."""
    return _MockBrickSize(np.diag(np.asarray(dims, dtype=float)))


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
            box_vectors=_mock_box_vectors(),
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
            box_vectors=_mock_box_vectors(),
            working_directory=tmp_path,
        )

        openff_mocks["center"].assert_called_once()
        np.testing.assert_allclose(
            result.positions[: original_solute.n_atoms], centered_solute.positions
        )
        assert not np.allclose(
            result.positions[: original_solute.n_atoms], original_solute.positions
        )

    def test_solvate_raises_on_overlapping_solvent(self, monkeypatch, tmp_path):
        """Solvent placed on top of the solute must abort the build."""
        from polyzymd.utils import packmol

        original_solute = _MockTopology(
            np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
        )
        centered_solute = _MockTopology(
            np.array([[15.0, 20.0, 25.0], [16.0, 21.0, 26.0], [17.0, 22.0, 27.0]])
        )
        solvent_topology = _MockTopology(np.zeros((3, 3), dtype=float))

        # Packmol output: solute atoms followed by solvent atoms that land
        # (after a frame mismatch) right on top of the centered solute.
        loaded_positions = np.array(
            [
                [15.0, 20.0, 25.0],
                [16.0, 21.0, 26.0],
                [17.0, 22.0, 27.0],
                [15.1, 20.0, 25.0],
                [16.0, 21.2, 26.0],
                [40.0, 40.0, 40.0],
            ]
        )

        self._install_fake_openff_modules(
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

        monkeypatch.setattr(packmol, "SOLVATION_CLASH_ATOM_LIMIT", 0)
        with pytest.raises(packmol.SolvationClashError) as excinfo:
            packmol.solvate_with_packmol(
                molecules=[object()],
                number_of_copies=[1],
                solute=original_solute,
                box_vectors=_mock_box_vectors(),
                working_directory=tmp_path,
                tolerance_angstrom=2.0,
            )

        message = str(excinfo.value)
        assert "2 solvent atom(s)" in message
        assert "d96b1fcd" in message
        assert "0.100" in message  # minimum separation

    def test_pack_polymers_raises_on_overlapping_polymer(self, monkeypatch, tmp_path):
        """Polymer atoms placed on top of the solute must abort the build."""
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
                [18.0, 23.0, 28.5],
                [41.0, 41.0, 41.0],
                [42.0, 42.0, 42.0],
            ]
        )

        self._install_fake_openff_modules(
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

        monkeypatch.setattr(packmol, "SOLVATION_CLASH_ATOM_LIMIT", 0)
        with pytest.raises(packmol.SolvationClashError, match="1 polymer atom"):
            packmol.pack_polymers(
                molecules=[object()],
                number_of_copies=[1],
                solute=original_solute,
                box_vectors=_mock_box_vectors(),
                working_directory=tmp_path,
            )


class TestSeparationStatistics:
    """Unit tests for the nearest-solute distance statistics."""

    def test_counts_and_minimum(self):
        from polyzymd.utils.packmol import separation_statistics

        solute = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
        other = np.array(
            [
                [0.5, 0.0, 0.0],  # 0.5 A  -> below half tolerance
                [1.5, 0.0, 0.0],  # 1.5 A  -> below tolerance only
                [10.0, 3.0, 0.0],  # 3.0 A -> fine
            ]
        )
        stats = separation_statistics(solute, other, tolerance_angstrom=2.0)
        assert stats["n_other"] == 3
        assert stats["n_below_tolerance"] == 2
        assert stats["n_below_half_tolerance"] == 1
        assert stats["min_distance_angstrom"] == pytest.approx(0.5)

    def test_empty_inputs(self):
        from polyzymd.utils.packmol import separation_statistics

        stats = separation_statistics(np.zeros((0, 3)), np.ones((4, 3)), tolerance_angstrom=2.0)
        assert stats["n_below_half_tolerance"] == 0
        assert stats["min_distance_angstrom"] == float("inf")

    def test_assert_warns_but_passes_between_half_and_full_tolerance(self, caplog):
        from polyzymd.utils.packmol import _assert_solute_solvent_separation

        topo = _MockTopology(np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]))
        with caplog.at_level("WARNING", logger="polyzymd.utils.packmol"):
            stats = _assert_solute_solvent_separation(topo, 1, tolerance_angstrom=2.0)
        assert stats["n_below_tolerance"] == 1
        assert any("not fully honoured" in rec.message for rec in caplog.records)

    def test_assert_skips_when_no_solute(self):
        from polyzymd.utils.packmol import _assert_solute_solvent_separation

        topo = _MockTopology(np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]))
        stats = _assert_solute_solvent_separation(topo, 0, tolerance_angstrom=2.0)
        assert stats["n_other"] == 0


class TestPackmolSeed:
    """Packmol ``seed`` keyword rendering and pass-through."""

    def test_seed_absent_by_default(self):
        from polyzymd.utils.packmol import build_packmol_input

        text = build_packmol_input(["w.pdb"], [10], np.array([30.0, 30.0, 30.0]), 2.0)
        assert "seed" not in text

    def test_seed_rendered_when_given(self):
        from polyzymd.utils.packmol import build_packmol_input

        text = build_packmol_input(["w.pdb"], [10], np.array([30.0, 30.0, 30.0]), 2.0, seed=7)
        lines = text.splitlines()
        assert "seed 7" in lines
        # Global keywords must precede the first structure block.
        assert lines.index("seed 7") < lines.index("structure w.pdb")

    def test_solvate_with_packmol_forwards_seed(self, monkeypatch, tmp_path):
        from polyzymd.utils import packmol

        original_solute = _MockTopology(np.array([[0.0, 0.0, 0.0]]))
        centered_solute = _MockTopology(np.array([[15.0, 20.0, 25.0]]))
        solvent_topology = _MockTopology(np.zeros((1, 3), dtype=float))
        loaded_positions = np.array([[15.0, 20.0, 25.0], [30.0, 30.0, 30.0]])

        TestSolvateAssemblyCoordinates._install_fake_openff_modules(
            monkeypatch,
            centered_solute=centered_solute,
            assembled_solvent=solvent_topology,
            loaded_positions=loaded_positions,
        )
        run_packmol = MagicMock(return_value=tmp_path / "packmol.pdb")
        monkeypatch.setattr(packmol, "run_packmol", run_packmol)
        monkeypatch.setattr(packmol, "_strip_conect_records", MagicMock(return_value=0))
        monkeypatch.setattr(packmol, "_check_pbc_available", MagicMock(return_value=False))
        monkeypatch.setattr(
            packmol, "_check_ignore_conect_supported", MagicMock(return_value=False)
        )

        packmol.solvate_with_packmol(
            molecules=[object()],
            number_of_copies=[1],
            solute=original_solute,
            box_vectors=_mock_box_vectors(),
            working_directory=tmp_path,
            seed=3,
        )

        assert "seed 3" in run_packmol.call_args.kwargs["input_text"].splitlines()

    def test_pack_polymers_forwards_seed(self, monkeypatch, tmp_path):
        from polyzymd.utils import packmol

        original_solute = _MockTopology(np.array([[0.0, 0.0, 0.0]]))
        centered_solute = _MockTopology(np.array([[15.0, 20.0, 25.0]]))
        polymer_topology = _MockTopology(np.zeros((1, 3), dtype=float))
        loaded_positions = np.array([[15.0, 20.0, 25.0], [40.0, 40.0, 40.0]])

        TestSolvateAssemblyCoordinates._install_fake_openff_modules(
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
        run_packmol = MagicMock(return_value=tmp_path / "packmol.pdb")
        monkeypatch.setattr(packmol, "run_packmol", run_packmol)
        monkeypatch.setattr(packmol, "_max_molecule_diameter_angstrom", MagicMock(return_value=1.0))

        packmol.pack_polymers(
            molecules=[object()],
            number_of_copies=[1],
            solute=original_solute,
            box_vectors=_mock_box_vectors(),
            working_directory=tmp_path,
            seed=5,
        )

        assert "seed 5" in run_packmol.call_args.kwargs["input_text"].splitlines()


class TestClashAtomLimit:
    """A few imperfect-packing contacts warn; a frame mismatch raises."""

    @staticmethod
    def _topology(n_clashing: int, n_clean: int = 50) -> _MockTopology:
        solute = np.array([[0.0, 0.0, 0.0]])
        clashing = np.tile([[0.6, 0.0, 0.0]], (n_clashing, 1))
        clean = np.column_stack(
            [np.linspace(10.0, 60.0, n_clean), np.zeros(n_clean), np.zeros(n_clean)]
        )
        return _MockTopology(np.vstack([solute, clashing, clean]))

    def test_residual_contacts_only_warn(self, caplog):
        from polyzymd.utils import packmol

        topo = self._topology(packmol.SOLVATION_CLASH_ATOM_LIMIT)
        with caplog.at_level("WARNING", logger="polyzymd.utils.packmol"):
            stats = packmol._assert_solute_solvent_separation(topo, 1, tolerance_angstrom=2.0)
        assert stats["n_below_half_tolerance"] == packmol.SOLVATION_CLASH_ATOM_LIMIT
        assert any("imperfect Packmol run" in rec.message for rec in caplog.records)

    def test_many_contacts_raise(self):
        from polyzymd.utils import packmol

        topo = self._topology(packmol.SOLVATION_CLASH_ATOM_LIMIT + 1)
        with pytest.raises(packmol.SolvationClashError, match="frame mismatch"):
            packmol._assert_solute_solvent_separation(topo, 1, tolerance_angstrom=2.0)

    def test_limit_is_far_below_defective_builds(self):
        from polyzymd.utils import packmol

        # Defective March-2026 CALB control: 1182 solvent atoms below 1 A.
        assert packmol.SOLVATION_CLASH_ATOM_LIMIT < 1182 / 10


class TestPolymerShellExclusion:
    """The solute bounding-box annulus is opt-in."""

    def _run(self, monkeypatch, tmp_path, **kwargs):
        from polyzymd.utils import packmol

        original_solute = _MockTopology(np.array([[0.0, 0.0, 0.0]]))
        centered_solute = _MockTopology(np.array([[15.0, 20.0, 25.0]]))
        polymer_topology = _MockTopology(np.zeros((1, 3), dtype=float))
        loaded_positions = np.array([[15.0, 20.0, 25.0], [40.0, 40.0, 40.0]])
        TestSolvateAssemblyCoordinates._install_fake_openff_modules(
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
        run_packmol = MagicMock(return_value=tmp_path / "packmol.pdb")
        monkeypatch.setattr(packmol, "run_packmol", run_packmol)
        monkeypatch.setattr(packmol, "_max_molecule_diameter_angstrom", MagicMock(return_value=1.0))
        packmol.pack_polymers(
            molecules=[object()],
            number_of_copies=[1],
            solute=original_solute,
            box_vectors=_mock_box_vectors(),
            working_directory=tmp_path,
            **kwargs,
        )
        return run_packmol.call_args.kwargs["input_text"], boxvectors_mod.get_topology_bbox_bounds

    def test_default_has_no_outside_box(self, monkeypatch, tmp_path):
        text, bbox = self._run(monkeypatch, tmp_path)
        assert "outside box" not in text
        assert "inside box" in text
        bbox.assert_not_called()

    def test_opt_in_restores_annulus(self, monkeypatch, tmp_path):
        text, bbox = self._run(monkeypatch, tmp_path, exclude_solute_bbox=True)
        assert "outside box 8.000000 8.000000 8.000000 22.000000 22.000000 22.000000" in text
        bbox.assert_called_once()

    def test_nloop_is_rendered(self, monkeypatch, tmp_path):
        text, _ = self._run(monkeypatch, tmp_path, nloop=50)
        assert "nloop 50" in text.splitlines()

    def test_packing_config_defaults(self):
        from polyzymd.config.schema import PolymerPackingConfig

        cfg = PolymerPackingConfig()
        assert cfg.exclude_solute_bbox is False
        assert cfg.nloop == 200


# ---------------------------------------------------------------------------
# Spherical confinement constraint
# ---------------------------------------------------------------------------


class TestInsideSphereConstraint:
    """``inside sphere`` is rendered only when a sphere is supplied."""

    SPHERE = np.array([10.0, 20.0, 30.0, 40.0])

    def test_sphere_absent_by_default(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL)
        assert "inside sphere" not in text

    def test_sphere_rendered_when_requested(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL, inside_sphere_angstrom=self.SPHERE)
        assert "  inside sphere 10.000000 20.000000 30.000000 40.000000" in text.splitlines()

    def test_sphere_in_every_molecule_block(self):
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL, inside_sphere_angstrom=self.SPHERE)
        assert text.count("inside sphere") == len(PDBS)

    def test_sphere_accompanies_inside_box(self):
        """The sphere narrows the packing box; it does not replace it."""
        text = build_packmol_input(PDBS, COUNTS, BOX_3A, TOL, inside_sphere_angstrom=self.SPHERE)
        assert text.index("inside box") < text.index("inside sphere")

    def test_sphere_absent_in_pbc_mode(self):
        text = build_packmol_input(
            PDBS, COUNTS, BOX_3A, TOL, use_pbc=True, inside_sphere_angstrom=self.SPHERE
        )
        assert "inside sphere" not in text

    def test_sphere_requires_four_values(self):
        with pytest.raises(ValueError, match=r"shape \(4,\)"):
            build_packmol_input(
                PDBS, COUNTS, BOX_3A, TOL, inside_sphere_angstrom=np.array([1.0, 2.0, 3.0])
            )

    def test_constraint_from_solute_topology(self):
        """Radius = bounding-box circumradius + padding, centre = centre of geometry."""
        from polyzymd.utils.packmol import solute_sphere_constraint

        solute = _MockTopology(np.array([[0.0, 0.0, 0.0], [6.0, 8.0, 0.0], [3.0, 4.0, 0.0]]))
        sphere = solute_sphere_constraint(solute, padding_angstrom=5.0)
        np.testing.assert_allclose(sphere[:3], [3.0, 4.0, 0.0])
        # bbox extent (6, 8, 0) -> diagonal 10 -> circumradius 5, plus 5 A padding
        assert sphere[3] == pytest.approx(10.0)

    def test_radius_is_independent_of_position(self):
        """Two framings of the same solute must give the same radius."""
        from polyzymd.utils.packmol import solute_sphere_constraint

        coords = np.array([[0.0, 0.0, 0.0], [6.0, 8.0, 0.0], [3.0, 4.0, 0.0]])
        here = solute_sphere_constraint(_MockTopology(coords), padding_angstrom=5.0)
        there = solute_sphere_constraint(_MockTopology(coords + 137.0), padding_angstrom=5.0)
        assert here[3] == pytest.approx(there[3])


# ---------------------------------------------------------------------------
# Periodic-image separation assertion
# ---------------------------------------------------------------------------

# A rhombic-dodecahedron cell in reduced form: c has x and y components equal
# to half of a and b, which is where the pre-fix builds put polymers on top of
# their own images.
TRICLINIC_BOX = np.array(
    [
        [40.0, 0.0, 0.0],
        [0.0, 40.0, 0.0],
        [20.0, 20.0, 28.284271],
    ]
)


class TestPeriodicImageSeparation:
    """Atoms must not overlap their own periodic images."""

    def test_clean_system_passes(self):
        from polyzymd.utils import packmol

        # Two atoms in the middle of the cell, far from every image.
        topo = _MockTopology(np.array([[18.0, 18.0, 12.0], [22.0, 22.0, 16.0]]))
        stats = packmol._assert_periodic_image_separation(
            topo, TRICLINIC_BOX, tolerance_angstrom=2.0, label="test"
        )
        assert stats["n_atoms_below_tolerance"] == 0
        assert stats["n_atoms_below_half_tolerance"] == 0

    def test_image_pair_raises(self):
        from polyzymd.utils import packmol

        # Second atom sits one full ``a`` vector away, plus 0.1 A: its image
        # across -a lands 0.1 A from the first atom.
        first = np.array([5.0, 18.0, 12.0])
        second = first + TRICLINIC_BOX[0] + np.array([0.1, 0.0, 0.0])
        topo = _MockTopology(np.vstack([first, second]))

        with pytest.raises(packmol.PeriodicImageClashError) as excinfo:
            packmol._assert_periodic_image_separation(
                topo, TRICLINIC_BOX, tolerance_angstrom=2.0, label="test system"
            )

        message = str(excinfo.value)
        assert "2 atom(s) of the test system" in message
        assert "0.100" in message
        assert "(0, 1)" in message or "(1, 0)" in message

    def test_image_pair_across_the_c_vector_raises(self):
        """The z face of the brick is where the pre-fix builds overlapped."""
        from polyzymd.utils import packmol

        first = np.array([12.0, 14.0, 1.0])
        second = first + TRICLINIC_BOX[2] + np.array([0.0, 0.0, 0.2])
        topo = _MockTopology(np.vstack([first, second]))

        with pytest.raises(packmol.PeriodicImageClashError, match="0.200"):
            packmol._assert_periodic_image_separation(
                topo, TRICLINIC_BOX, tolerance_angstrom=2.0, label="test"
            )

    def test_atom_is_not_counted_against_its_own_image(self):
        """A lone atom is never in contact with itself, only with its images."""
        from polyzymd.utils.packmol import periodic_image_statistics

        stats = periodic_image_statistics(
            np.array([[20.0, 20.0, 14.0]]), TRICLINIC_BOX, tolerance_angstrom=2.0
        )
        assert stats["n_atoms"] == 1
        assert stats["n_atoms_below_tolerance"] == 0
        assert stats["min_distance_angstrom"] == float("inf")

    def test_self_image_contact_is_detected(self):
        """An atom too close to its own image across a short cell must fail."""
        from polyzymd.utils import packmol

        tiny_box = np.array([[0.8, 0.0, 0.0], [0.0, 40.0, 0.0], [0.0, 0.0, 40.0]])
        topo = _MockTopology(np.array([[0.4, 20.0, 20.0]]))
        with pytest.raises(packmol.PeriodicImageClashError, match="0.800"):
            packmol._assert_periodic_image_separation(
                topo, tiny_box, tolerance_angstrom=2.0, label="test"
            )

    def test_marginal_contact_only_warns(self, caplog):
        from polyzymd.utils import packmol

        first = np.array([5.0, 18.0, 12.0])
        second = first + TRICLINIC_BOX[0] + np.array([1.5, 0.0, 0.0])
        topo = _MockTopology(np.vstack([first, second]))

        with caplog.at_level("WARNING", logger="polyzymd.utils.packmol"):
            stats = packmol._assert_periodic_image_separation(
                topo, TRICLINIC_BOX, tolerance_angstrom=2.0, label="test"
            )
        assert stats["n_atoms_below_tolerance"] == 2
        assert stats["n_atoms_below_half_tolerance"] == 0
        assert any("periodic image" in record.message for record in caplog.records)

    def test_accepts_quantity_box_vectors(self):
        from polyzymd.utils import packmol

        topo = _MockTopology(np.array([[18.0, 18.0, 12.0], [22.0, 22.0, 16.0]]))
        stats = packmol._assert_periodic_image_separation(
            topo, _MockBrickSize(TRICLINIC_BOX), tolerance_angstrom=2.0, label="test"
        )
        assert stats["n_atoms_below_tolerance"] == 0

    def test_worst_pair_and_lattice_vector_are_reported(self):
        from polyzymd.utils.packmol import periodic_image_statistics

        first = np.array([5.0, 18.0, 12.0])
        second = first + TRICLINIC_BOX[1] + np.array([0.0, 0.3, 0.0])
        stats = periodic_image_statistics(
            np.vstack([first, second]), TRICLINIC_BOX, tolerance_angstrom=2.0
        )
        assert stats["min_distance_angstrom"] == pytest.approx(0.3)
        assert set(stats["worst_pair"]) == {0, 1}
        assert stats["worst_lattice_vector"] in {(0, 1, 0), (0, -1, 0)}


class TestPackPolymersFinalBox:
    """Polymers are packed inside the supplied final brick, confined to a sphere."""

    def _run(self, monkeypatch, tmp_path, **kwargs):
        from polyzymd.utils import packmol

        original_solute = _MockTopology(np.array([[0.0, 0.0, 0.0]]))
        centered_solute = _MockTopology(np.array([[15.0, 20.0, 25.0], [17.0, 20.0, 25.0]]))
        polymer_topology = _MockTopology(np.zeros((1, 3), dtype=float))
        loaded_positions = np.array([[15.0, 20.0, 25.0], [17.0, 20.0, 25.0], [5.0, 8.0, 10.0]])
        TestSolvateAssemblyCoordinates._install_fake_openff_modules(
            monkeypatch,
            centered_solute=centered_solute,
            assembled_solvent=polymer_topology,
            loaded_positions=loaded_positions,
        )
        run_packmol = MagicMock(return_value=tmp_path / "packmol.pdb")
        monkeypatch.setattr(packmol, "run_packmol", run_packmol)
        monkeypatch.setattr(packmol, "_max_molecule_diameter_angstrom", MagicMock(return_value=1.0))

        packmol.pack_polymers(
            molecules=[object()],
            number_of_copies=[1],
            solute=original_solute,
            box_vectors=_mock_box_vectors(),
            working_directory=tmp_path,
            **kwargs,
        )
        return run_packmol.call_args.kwargs["input_text"]

    def test_packing_box_is_the_supplied_brick(self, monkeypatch, tmp_path):
        """The brick of the final cell (30, 40, 50), shrunk by the 2 A tolerance."""
        text = self._run(monkeypatch, tmp_path)
        assert "  inside box 0. 0. 0. 28.000000 38.000000 48.000000" in text.splitlines()

    def test_sphere_constraint_is_added_by_default(self, monkeypatch, tmp_path):
        text = self._run(monkeypatch, tmp_path, sphere_padding_angstrom=20.0)
        # centred solute spans 2 A in x -> circumradius 1 A, plus 20 A padding
        assert "  inside sphere 16.000000 20.000000 25.000000 21.000000" in text.splitlines()

    def test_sphere_can_be_disabled(self, monkeypatch, tmp_path):
        text = self._run(monkeypatch, tmp_path, confine_to_sphere=False)
        assert "inside sphere" not in text
        assert "inside box" in text


class TestSolvateCenterSolute:
    """``center_solute=False`` keeps an already brick-framed topology in place."""

    def _run(self, monkeypatch, tmp_path, **kwargs):
        from polyzymd.utils import packmol

        original_solute = _MockTopology(np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        centered_solute = _MockTopology(np.array([[15.0, 20.0, 25.0], [16.0, 21.0, 26.0]]))
        solvent_topology = _MockTopology(np.zeros((1, 3), dtype=float))
        loaded_positions = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [25.0, 25.0, 25.0]])

        mocks = TestSolvateAssemblyCoordinates._install_fake_openff_modules(
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
            box_vectors=_mock_box_vectors(),
            working_directory=tmp_path,
            **kwargs,
        )
        return result, mocks, original_solute

    def test_centring_is_skipped(self, monkeypatch, tmp_path):
        result, mocks, original_solute = self._run(monkeypatch, tmp_path, center_solute=False)
        mocks["center"].assert_not_called()
        np.testing.assert_allclose(
            result.positions[: original_solute.n_atoms], original_solute.positions
        )

    def test_centring_is_the_default(self, monkeypatch, tmp_path):
        _, mocks, _ = self._run(monkeypatch, tmp_path)
        mocks["center"].assert_called_once()
