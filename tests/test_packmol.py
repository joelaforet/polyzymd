"""
Tests for polyzymd.utils.packmol — input-file generation utilities.

These tests exercise :func:`build_packmol_input` directly and do NOT
require a Packmol binary or any heavy simulation dependencies.
"""

from __future__ import annotations

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
