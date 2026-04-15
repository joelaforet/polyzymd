"""Regression tests for TrajectoryLoader engine-aware refactor.

Tests cover:
- OpenMM daisy-chain and legacy directory layouts
- GROMACS flat production layout
- GRO topology chain-ID warning
- engine_override parameter
- Lazy engine creation
- find_topology() backward compatibility for direct callers
- _infer_replicate() edge cases
- Error paths (missing topology, missing trajectories, missing working dir)
- Fallback when engine name is unrecognisable
"""

from __future__ import annotations

import logging
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from polyzymd.analyses.shared.loader import TrajectoryLoader

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_openmm_config(tmp_path: Path, *, engine: str = "openmm") -> MagicMock:
    """Build a minimal mock SimulationConfig for OpenMM tests."""
    config = MagicMock()
    config.engine = engine
    config.get_working_directory.side_effect = lambda rep: tmp_path / f"run_{rep}"
    config.output.effective_scratch_directory = tmp_path
    config.generate_system_name.return_value = "test_system"
    return config


def _make_gromacs_config(tmp_path: Path) -> MagicMock:
    """Build a minimal mock SimulationConfig for GROMACS tests."""
    config = MagicMock()
    config.engine = "gromacs"
    config.get_working_directory.side_effect = lambda rep: tmp_path / f"run_{rep}"
    config.output.effective_scratch_directory = tmp_path
    config.generate_system_name.return_value = "test_system"
    return config


def _create_openmm_daisy_chain(run_dir: Path, n_segments: int = 3) -> None:
    """Populate *run_dir* with a daisy-chain OpenMM directory layout."""
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "solvated_system.pdb").write_text("ATOM solvated")
    for i in range(n_segments):
        seg_dir = run_dir / f"production_{i}"
        seg_dir.mkdir()
        (seg_dir / f"production_{i}_trajectory.dcd").write_bytes(b"\x00")
        (seg_dir / f"production_{i}_topology.pdb").write_text("ATOM topology")


def _create_openmm_legacy(run_dir: Path) -> None:
    """Populate *run_dir* with a legacy OpenMM single-directory layout."""
    run_dir.mkdir(parents=True, exist_ok=True)
    prod_dir = run_dir / "production"
    prod_dir.mkdir()
    (prod_dir / "production_topology.pdb").write_text("ATOM topology")
    (prod_dir / "production_trajectory.dcd").write_bytes(b"\x00")


def _create_gromacs_flat(run_dir: Path, *, use_gro: bool = False) -> None:
    """Populate *run_dir* with a flat GROMACS production layout."""
    gromacs_dir = run_dir / "gromacs"
    gromacs_dir.mkdir(parents=True, exist_ok=True)
    if use_gro:
        (gromacs_dir / "solvated_system.gro").write_text("GRO content")
    else:
        (gromacs_dir / "solvated_system.pdb").write_text("ATOM topology")
    (gromacs_dir / "prod.xtc").write_bytes(b"\x00")


# ---------------------------------------------------------------------------
# Lazy engine creation
# ---------------------------------------------------------------------------


class TestLazyEngine:
    """Engine should not be created until the first resolution call."""

    def test_engine_not_created_at_init(self, tmp_path):
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)
        assert loader._engine is None

    def test_engine_created_on_first_find_topology(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        loader.find_topology(run_dir)
        assert loader._engine is not None

    def test_engine_cached_across_calls(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        loader.find_topology(run_dir)
        first_engine = loader._engine
        loader.find_topology(run_dir)
        assert loader._engine is first_engine


# ---------------------------------------------------------------------------
# OpenMM layouts
# ---------------------------------------------------------------------------


class TestOpenMMDaisyChain:
    """Loader resolves OpenMM daisy-chain segmented layouts."""

    def test_finds_topology(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=2)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert info.topology_file == run_dir / "solvated_system.pdb"

    def test_finds_all_segments_in_order(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=3)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert info.n_segments == 3
        for i, f in enumerate(info.trajectory_files):
            assert f.name == f"production_{i}_trajectory.dcd"

    def test_replicate_routing(self, tmp_path):
        for rep in (1, 3):
            run_dir = tmp_path / f"run_{rep}"
            _create_openmm_daisy_chain(run_dir, n_segments=1)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        info1 = loader.get_trajectory_info(replicate=1)
        info3 = loader.get_trajectory_info(replicate=3)
        assert info1.working_directory != info3.working_directory
        assert info1.replicate == 1
        assert info3.replicate == 3


class TestOpenMMLegacy:
    """Loader resolves legacy single-directory OpenMM layout."""

    def test_finds_legacy_topology(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_legacy(run_dir)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert info.topology_file.name == "production_topology.pdb"

    def test_finds_legacy_trajectory(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_legacy(run_dir)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert len(info.trajectory_files) == 1
        assert info.trajectory_files[0].name == "production_trajectory.dcd"


# ---------------------------------------------------------------------------
# GROMACS layouts
# ---------------------------------------------------------------------------


class TestGromacsFlat:
    """Loader resolves flat GROMACS production layout."""

    def test_finds_pdb_topology(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert info.topology_file.suffix == ".pdb"

    def test_finds_prod_xtc(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert len(info.trajectory_files) == 1
        assert info.trajectory_files[0].name == "prod.xtc"

    def test_n_segments_is_one(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        assert info.n_segments == 1


# ---------------------------------------------------------------------------
# GRO topology warning
# ---------------------------------------------------------------------------


class TestGROWarning:
    """GRO topology triggers a chain-ID warning exactly once."""

    def test_warns_on_gro_topology(self, tmp_path, caplog):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir, use_gro=True)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        with caplog.at_level(logging.WARNING, logger="polyzymd.analyses.shared.loader"):
            loader.find_topology(run_dir)

        assert any(
            "GRO topology" in record.message and record.name == "polyzymd.analyses.shared.loader"
            for record in caplog.records
        )

    def test_warns_only_once_for_same_path(self, tmp_path, caplog):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir, use_gro=True)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        with caplog.at_level(logging.WARNING, logger="polyzymd.analyses.shared.loader"):
            loader.find_topology(run_dir)
            loader.find_topology(run_dir)

        gro_warnings = [
            record.message
            for record in caplog.records
            if "GRO topology" in record.message and record.name == "polyzymd.analyses.shared.loader"
        ]
        assert len(gro_warnings) == 1

    def test_no_warning_for_pdb_topology(self, tmp_path, caplog):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir, use_gro=False)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config)

        with caplog.at_level(logging.WARNING, logger="polyzymd.analyses.shared.loader"):
            loader.find_topology(run_dir)

        gro_warnings = [m for m in caplog.messages if "GRO topology" in m]
        assert len(gro_warnings) == 0


# ---------------------------------------------------------------------------
# engine_override parameter
# ---------------------------------------------------------------------------


class TestEngineOverride:
    """engine_override forces a specific engine regardless of config.engine."""

    def test_override_gromacs_on_openmm_config(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_gromacs_flat(run_dir)
        # Config says openmm, but override says gromacs
        config = _make_openmm_config(tmp_path, engine="openmm")
        loader = TrajectoryLoader(config, engine_override="gromacs")

        info = loader.get_trajectory_info(replicate=1)
        assert info.trajectory_files[0].name == "prod.xtc"

    def test_override_openmm_on_gromacs_config(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=2)
        config = _make_gromacs_config(tmp_path)
        loader = TrajectoryLoader(config, engine_override="openmm")

        info = loader.get_trajectory_info(replicate=1)
        assert info.topology_file.name == "solvated_system.pdb"
        assert len(info.trajectory_files) == 2


# ---------------------------------------------------------------------------
# find_topology backward compatibility
# ---------------------------------------------------------------------------


class TestFindTopologyCompat:
    """find_topology(working_dir) works for direct plugin callers."""

    def test_arbitrary_working_dir(self, tmp_path):
        """Plugin passes an explicit working_dir, not from config."""
        arbitrary_dir = tmp_path / "some_other_dir"
        arbitrary_dir.mkdir()
        (arbitrary_dir / "solvated_system.pdb").write_text("ATOM")

        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        result = loader.find_topology(arbitrary_dir)
        assert result == arbitrary_dir / "solvated_system.pdb"

    def test_raises_when_no_topology(self, tmp_path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()

        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        with pytest.raises(FileNotFoundError, match="No topology file found"):
            loader.find_topology(empty_dir)


# ---------------------------------------------------------------------------
# _infer_replicate edge cases
# ---------------------------------------------------------------------------


class TestInferReplicate:
    """Best-effort replicate parsing from directory names."""

    def test_parses_run_1(self, tmp_path):
        result = TrajectoryLoader._infer_replicate(tmp_path / "run_1")
        assert result == 1

    def test_parses_run_42(self, tmp_path):
        result = TrajectoryLoader._infer_replicate(tmp_path / "run_42")
        assert result == 42

    def test_fallback_for_unknown_dir_name(self, tmp_path):
        result = TrajectoryLoader._infer_replicate(tmp_path / "my_simulation")
        assert result == 1

    def test_handles_non_path_gracefully(self):
        """MagicMock or other non-Path objects fall back to 1."""
        result = TrajectoryLoader._infer_replicate(MagicMock())
        assert result == 1


# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------


class TestErrorPaths:
    """FileNotFoundError raised for missing files or directories."""

    def test_missing_working_directory(self, tmp_path):
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        with pytest.raises(FileNotFoundError, match="Working directory not found"):
            loader.get_trajectory_info(replicate=99)

    def test_no_topology_in_get_trajectory_info(self, tmp_path):
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        # Trajectory but no topology
        prod_dir = run_dir / "production_0"
        prod_dir.mkdir()
        (prod_dir / "production_0_trajectory.dcd").write_bytes(b"\x00")

        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        with pytest.raises(FileNotFoundError, match="No topology file found"):
            loader.get_trajectory_info(replicate=1)

    def test_no_trajectories_in_get_trajectory_info(self, tmp_path):
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        (run_dir / "solvated_system.pdb").write_text("ATOM")
        # Topology but no trajectories

        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        with pytest.raises(FileNotFoundError, match="No production trajectory files found"):
            loader.get_trajectory_info(replicate=1)

    def test_find_trajectories_raises_on_empty_dir(self, tmp_path):
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()

        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        with pytest.raises(FileNotFoundError, match="No production trajectory files found"):
            loader._find_trajectories(run_dir)


# ---------------------------------------------------------------------------
# Mock config fallback
# ---------------------------------------------------------------------------


class TestMockConfigFallback:
    """Unrecognised engine names fall back to OpenMM resolver."""

    def test_mock_engine_falls_back_to_openmm(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=1)

        config = MagicMock()
        config.engine = "nonexistent_engine"
        config.get_working_directory.side_effect = lambda rep: tmp_path / f"run_{rep}"
        config.output.effective_scratch_directory = tmp_path

        loader = TrajectoryLoader(config)
        topo = loader.find_topology(run_dir)
        assert topo == run_dir / "solvated_system.pdb"

    def test_fully_mocked_config_works(self, tmp_path):
        """Purely mocked config (no engine attr value) still works."""
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=1)

        config = MagicMock()
        config.get_working_directory.side_effect = lambda rep: tmp_path / f"run_{rep}"
        config.output.effective_scratch_directory = tmp_path

        loader = TrajectoryLoader(config)
        topo = loader.find_topology(run_dir)
        assert topo == run_dir / "solvated_system.pdb"


# ---------------------------------------------------------------------------
# TrajectoryInfo validation
# ---------------------------------------------------------------------------


class TestTrajectoryInfoValidation:
    """TrajectoryInfo.validate() checks file existence."""

    def test_validate_passes_with_real_files(self, tmp_path):
        run_dir = tmp_path / "run_1"
        _create_openmm_daisy_chain(run_dir, n_segments=1)
        config = _make_openmm_config(tmp_path)
        loader = TrajectoryLoader(config)

        info = loader.get_trajectory_info(replicate=1)
        info.validate()  # Should not raise

    def test_validate_fails_missing_topology(self, tmp_path):
        from polyzymd.analyses.shared.loader import TrajectoryInfo

        info = TrajectoryInfo(
            topology_file=tmp_path / "nonexistent.pdb",
            trajectory_files=[],
        )
        with pytest.raises(FileNotFoundError, match="Topology not found"):
            info.validate()
