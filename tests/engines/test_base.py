"""Tests for engine base models."""

from pathlib import Path

from polyzymd.engines.base import EngineSubmitRequest, TrajectoryLayout


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
