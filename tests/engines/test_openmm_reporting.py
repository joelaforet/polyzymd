"""Tests for OpenMM production reporting cadence."""

from pathlib import Path
from unittest.mock import MagicMock, patch

from polyzymd.engines.openmm.engine import OpenMMEngine


def test_run_local_derives_report_interval_from_samples() -> None:
    """Requested trajectory samples determine the reporter interval."""
    config = MagicMock()
    production = config.simulation_phases.production
    production.duration = 1000.0
    production.samples = 2500
    production.time_step = 2.0
    production.checkpoint_interval = 60.0

    with patch("polyzymd.cli.main._run_initial_segment") as run_initial_segment:
        OpenMMEngine(config=config).run_local(
            replicate=1,
            working_dir=Path("run"),
        )

    run_initial_segment.assert_called_once_with(
        sim_config=config,
        working_dir=Path("run"),
        replicate=1,
        skip_build=False,
        duration_ns=1000.0,
        num_samples=2500,
        timestep_fs=2.0,
        report_interval=200000,
        checkpoint_interval_s=60.0,
    )
