"""Tests for engine base models."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

from polyzymd.engines.base import EngineSubmitRequest, TrajectoryLayout
from polyzymd.engines.gromacs.engine import GromacsEngine
from polyzymd.engines.openmm.engine import OpenMMEngine


def test_trajectory_layout_creation() -> None:
    """TrajectoryLayout should store provided paths and formats."""
    layout = TrajectoryLayout(
        topology_path=Path("system.pdb"),
        trajectory_paths=[Path("prod_0.dcd"), Path("prod_1.dcd")],
        trajectory_format="dcd",
        topology_format="pdb",
    )

    assert layout.topology_path == Path("system.pdb")
    assert layout.trajectory_paths == [Path("prod_0.dcd"), Path("prod_1.dcd")]
    assert layout.trajectory_format == "dcd"
    assert layout.topology_format == "pdb"


def test_engine_submit_request_creation() -> None:
    """EngineSubmitRequest should capture submission metadata."""
    request = EngineSubmitRequest(
        replicate=2,
        config_path=Path("config.yaml"),
        working_dir=Path("scratch/run_2"),
        job_name="run_2",
    )

    assert request.replicate == 2
    assert request.config_path == Path("config.yaml")
    assert request.working_dir == Path("scratch/run_2")
    assert request.job_name == "run_2"


class TestEngineSubdir:
    """Tests for engine_subdir class variable."""

    def test_gromacs_engine_subdir(self):
        """GROMACS engine should have 'gromacs' subdir."""
        assert GromacsEngine.engine_subdir == "gromacs"

    def test_openmm_engine_subdir_is_none(self):
        """OpenMM engine should have no subdir (None)."""
        assert OpenMMEngine.engine_subdir is None


class TestResolveEngineWorkingDirectory:
    """Tests for resolve_engine_working_directory() method."""

    def test_gromacs_appends_subdir(self, tmp_path):
        """GROMACS should append /gromacs/ to replicate root."""
        config = MagicMock()
        config.gromacs.gmx_binary = None
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        result = engine.resolve_engine_working_directory(tmp_path / "system_run1")
        assert result == tmp_path / "system_run1" / "gromacs"

    def test_openmm_returns_root_unchanged(self, tmp_path):
        """OpenMM should return the replicate root as-is."""
        config = MagicMock()
        engine = OpenMMEngine(config=config)
        result = engine.resolve_engine_working_directory(tmp_path / "system_run1")
        assert result == tmp_path / "system_run1"


class TestGetEngineWorkingDirectory:
    """Tests for get_engine_working_directory() method."""

    def test_gromacs_uses_scratch_plus_subdir(self, tmp_path):
        """GROMACS should return scratch/{template}_runN/gromacs/."""
        config = MagicMock()
        config.get_working_directory.return_value = (
            tmp_path / "scratch" / "CALB_PEG_100ns_343K_run1"
        )
        config.gromacs.gmx_binary = None
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        result = engine.get_engine_working_directory(config, 1)
        assert result == tmp_path / "scratch" / "CALB_PEG_100ns_343K_run1" / "gromacs"
        config.get_working_directory.assert_called_once_with(1)

    def test_openmm_uses_scratch_root(self, tmp_path):
        """OpenMM should return scratch/{template}_runN/ (no subdir)."""
        config = MagicMock()
        config.get_working_directory.return_value = (
            tmp_path / "scratch" / "CALB_PEG_100ns_343K_run1"
        )
        engine = OpenMMEngine(config=config)
        result = engine.get_engine_working_directory(config, 1)
        assert result == tmp_path / "scratch" / "CALB_PEG_100ns_343K_run1"
