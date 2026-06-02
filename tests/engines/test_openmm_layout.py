"""Tests for OpenMM trajectory layout resolution."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from polyzymd.engines.openmm.engine import OpenMMEngine


def _make_engine() -> OpenMMEngine:
    """Create an OpenMMEngine with a minimal mock config.

    Returns
    -------
    OpenMMEngine
        Engine instance with required production settings.
    """
    config = MagicMock()
    config.simulation_phases.production.duration = 100.0
    config.simulation_phases.production.time_step = 2.0
    config.simulation_phases.production.samples = 250
    return OpenMMEngine(config=config)


class TestOpenMMTopologySearch:
    """Tests for topology search priority."""

    def test_prefers_solvated_system_pdb(self, tmp_path: Path) -> None:
        """solvated_system.pdb is found first even with production PDBs."""
        (tmp_path / "solvated_system.pdb").write_text("ATOM")
        prod0 = tmp_path / "production_0"
        prod0.mkdir()
        (prod0 / "production_0_topology.pdb").write_text("ATOM")
        (prod0 / "production_0_trajectory.dcd").write_bytes(b"DCD")

        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.topology_path is not None
        assert layout.topology_path.name == "solvated_system.pdb"

    def test_falls_back_to_production_0_topology(self, tmp_path: Path) -> None:
        """production_0/production_0_topology.pdb is found when no solvated_system."""
        prod0 = tmp_path / "production_0"
        prod0.mkdir()
        (prod0 / "production_0_topology.pdb").write_text("ATOM")
        (prod0 / "production_0_trajectory.dcd").write_bytes(b"DCD")

        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.topology_path is not None
        assert layout.topology_path.name == "production_0_topology.pdb"

    def test_falls_back_to_legacy_production_topology(self, tmp_path: Path) -> None:
        """production/production_topology.pdb is found as legacy fallback."""
        prod = tmp_path / "production"
        prod.mkdir()
        (prod / "production_topology.pdb").write_text("ATOM")

        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.topology_path is not None
        assert layout.topology_path.name == "production_topology.pdb"

    def test_rejects_arbitrary_root_pdb(self, tmp_path: Path) -> None:
        """Arbitrary root PDBs are not treated as OpenMM topology."""
        (tmp_path / "custom_topology.pdb").write_text("ATOM")

        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.topology_path is None

    def test_rejects_non_exact_topology_name(self, tmp_path: Path) -> None:
        """Non-canonical topology names are rejected even in production dirs."""
        prod = tmp_path / "production_0"
        prod.mkdir()
        (prod / "custom_topology.pdb").write_text("ATOM")
        (prod / "production_0_trajectory.dcd").write_bytes(b"DCD")

        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.topology_path is None

    def test_no_files_returns_none(self, tmp_path: Path) -> None:
        """Empty directory returns None topology."""
        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.topology_path is None


class TestOpenMMTrajectorySearch:
    """Tests for trajectory search priority."""

    def test_finds_daisy_chain_segments_in_order(self, tmp_path: Path) -> None:
        """production_N/production_N_trajectory.dcd files are found in order."""
        for i in [0, 1, 2]:
            directory = tmp_path / f"production_{i}"
            directory.mkdir()
            (directory / f"production_{i}_trajectory.dcd").write_bytes(b"DCD")

        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert len(layout.trajectory_paths) == 3
        assert layout.trajectory_paths[0].name == "production_0_trajectory.dcd"
        assert layout.trajectory_paths[2].name == "production_2_trajectory.dcd"

    def test_finds_legacy_single_trajectory(self, tmp_path: Path) -> None:
        """production/production_trajectory.dcd is found when no daisy-chain."""
        prod = tmp_path / "production"
        prod.mkdir()
        (prod / "production_trajectory.dcd").write_bytes(b"DCD")

        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert len(layout.trajectory_paths) == 1
        assert layout.trajectory_paths[0].name == "production_trajectory.dcd"

    def test_rejects_empty_daisy_chain_segment(self, tmp_path: Path) -> None:
        """Empty canonical daisy-chain trajectories are rejected."""
        prod0 = tmp_path / "production_0"
        prod0.mkdir()
        (prod0 / "production_0_trajectory.dcd").write_bytes(b"")

        engine = _make_engine()
        with pytest.raises(ValueError, match="Empty OpenMM trajectory segment"):
            engine.resolve_trajectory_layout(tmp_path, replicate=1)

    def test_rejects_gapped_daisy_chain_segments(self, tmp_path: Path) -> None:
        """Canonical daisy-chain segment indices must be contiguous."""
        for i in [0, 2]:
            directory = tmp_path / f"production_{i}"
            directory.mkdir()
            (directory / f"production_{i}_trajectory.dcd").write_bytes(b"DCD")

        engine = _make_engine()
        with pytest.raises(ValueError, match="production_1"):
            engine.resolve_trajectory_layout(tmp_path, replicate=1)

    def test_rejects_arbitrary_recursive_production_trajectory(self, tmp_path: Path) -> None:
        """Recursive production trajectory globs are intentionally not used."""
        subdir = tmp_path / "subdir"
        subdir.mkdir()
        (subdir / "production_custom_trajectory.dcd").write_bytes(b"DCD")

        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.trajectory_paths == []

    def test_no_trajectories_returns_empty(self, tmp_path: Path) -> None:
        """Empty directory returns empty trajectory list."""
        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.trajectory_paths == []

    def test_layout_format_is_dcd_pdb(self, tmp_path: Path) -> None:
        """Layout always reports DCD/PDB format for OpenMM."""
        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.trajectory_format == "dcd"
        assert layout.topology_format == "pdb"
