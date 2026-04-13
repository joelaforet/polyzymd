"""Tests for GROMACS trajectory layout resolution."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from polyzymd.engines.gromacs.engine import GromacsEngine


def _make_engine(system_name: str = "system") -> GromacsEngine:
    """Create a GromacsEngine with a mock config.

    Parameters
    ----------
    system_name : str, optional
        Name returned by ``generate_system_name``, by default ``"system"``.

    Returns
    -------
    GromacsEngine
        Engine instance with minimal required settings.
    """
    config = MagicMock()
    config.generate_system_name.return_value = system_name
    config.simulation_phases.production.duration = 100.0
    config.simulation_phases.production.time_step = 2.0
    config.simulation_phases.production.samples = 250
    config.gromacs.gmx_binary = None
    config.gromacs.grompp_flags = ""
    config.gromacs.mdrun_flags = ""
    config.gromacs.module_load = None
    return GromacsEngine(config=config, gmx_binary="gmx")


class TestGromacsTopologySearch:
    """Tests for GROMACS topology resolution."""

    def test_prefers_named_pdb(self, tmp_path: Path) -> None:
        """<prefix>.pdb is preferred over generic PDB."""
        (tmp_path / "system.pdb").write_text("ATOM")
        (tmp_path / "other.pdb").write_text("ATOM")

        engine = _make_engine("system")
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.topology_path is not None
        assert layout.topology_path.name == "system.pdb"
        assert layout.topology_format == "pdb"

    def test_falls_back_to_any_pdb(self, tmp_path: Path) -> None:
        """Any PDB is found if named PDB is absent."""
        (tmp_path / "topology.pdb").write_text("ATOM")

        engine = _make_engine("system")
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.topology_path is not None
        assert layout.topology_path.name == "topology.pdb"

    def test_gro_fallback(self, tmp_path: Path) -> None:
        """GRO is used when no PDB exists."""
        (tmp_path / "system.gro").write_text("GRO")

        engine = _make_engine("system")
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.topology_path is not None
        assert layout.topology_path.name == "system.gro"
        assert layout.topology_format == "gro"

    def test_pdb_preferred_over_gro(self, tmp_path: Path) -> None:
        """PDB is always preferred over GRO for chain ID preservation."""
        (tmp_path / "system.pdb").write_text("ATOM")
        (tmp_path / "system.gro").write_text("GRO")

        engine = _make_engine("system")
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.topology_format == "pdb"

    def test_no_topology_returns_none(self, tmp_path: Path) -> None:
        """Empty directory returns None topology."""
        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.topology_path is None


class TestGromacsTrajectorySearch:
    """Tests for GROMACS trajectory resolution."""

    def test_prefers_prod_xtc(self, tmp_path: Path) -> None:
        """prod.xtc is preferred over other XTCs."""
        (tmp_path / "prod.xtc").write_bytes(b"XTC")
        (tmp_path / "other.xtc").write_bytes(b"XTC")

        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert len(layout.trajectory_paths) == 1
        assert layout.trajectory_paths[0].name == "prod.xtc"

    def test_falls_back_to_any_xtc(self, tmp_path: Path) -> None:
        """Any XTC found when prod.xtc absent."""
        (tmp_path / "trajectory.xtc").write_bytes(b"XTC")

        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert len(layout.trajectory_paths) == 1

    def test_no_trajectories(self, tmp_path: Path) -> None:
        """Empty directory returns empty trajectory list."""
        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.trajectory_paths == []

    def test_layout_format_is_xtc(self, tmp_path: Path) -> None:
        """Layout reports XTC format for GROMACS."""
        engine = _make_engine()
        layout = engine.resolve_trajectory_layout(tmp_path, replicate=1)
        assert layout.trajectory_format == "xtc"
