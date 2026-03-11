"""Tests for the recover and check-progress CLI commands.

Covers:
- _find_topology_pdb helper
- recover CLI with new self-resubmitting interface
- check-progress CLI (exit code 0 = complete, 1 = work remains)
- Self-resubmitting model (no afterany dependencies)
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from polyzymd.cli.main import _find_topology_pdb, cli

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_completed_segment(seg_dir: Path, seg_idx: int) -> None:
    """Write minimal files to simulate a completed segment."""
    seg_dir.mkdir(parents=True, exist_ok=True)
    (seg_dir / f"production_{seg_idx}_state.xml").write_text("<mock/>")
    (seg_dir / f"production_{seg_idx}_system.xml").write_text("<mock/>")
    (seg_dir / f"production_{seg_idx}_parameters.json").write_text("{}")


# ---------------------------------------------------------------------------
# _find_topology_pdb
# ---------------------------------------------------------------------------


class TestFindTopologyPdb:
    """_find_topology_pdb locates a suitable PDB file."""

    def test_finds_solvated_pdb(self, tmp_path):
        (tmp_path / "system_solvated.pdb").write_text("ATOM ...")
        result = _find_topology_pdb(tmp_path)
        assert result.name == "system_solvated.pdb"

    def test_finds_equilibration_topology(self, tmp_path):
        eq_dir = tmp_path / "equilibration"
        eq_dir.mkdir()
        (eq_dir / "equil_topology.pdb").write_text("ATOM ...")
        result = _find_topology_pdb(tmp_path)
        assert "topology.pdb" in result.name

    def test_finds_production_0_topology(self, tmp_path):
        prod_dir = tmp_path / "production_0"
        prod_dir.mkdir()
        (prod_dir / "production_0_topology.pdb").write_text("ATOM ...")
        result = _find_topology_pdb(tmp_path)
        assert result.name == "production_0_topology.pdb"

    def test_raises_if_no_pdb(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Could not find topology PDB"):
            _find_topology_pdb(tmp_path)


# ---------------------------------------------------------------------------
# recover CLI (new self-resubmitting interface)
# ---------------------------------------------------------------------------


def _mock_sim_config(working_dir: Path):
    """Create a mock SimulationConfig that returns working_dir."""
    mock = MagicMock()
    mock.get_working_directory.return_value = working_dir
    mock.simulation_phases.production.time_step = 2.0
    mock.simulation_phases.production.duration = 20.0  # ns
    mock.simulation_phases.production.samples = 250
    return mock


def _mock_progress(
    total_steps=10000000,
    completed_steps=5000000,
    is_complete=False,
    n_segments=1,
):
    """Create a mock SimulationProgress."""
    from polyzymd.simulation.progress import (
        SegmentRecord,
        SegmentStatus,
        SimulationProgress,
        SimulationStatus,
    )

    segments = []
    for i in range(n_segments):
        seg = SegmentRecord(
            index=i,
            steps_completed=completed_steps // max(n_segments, 1),
            steps_requested=completed_steps // max(n_segments, 1),
            samples_written=125,
            status=SegmentStatus.COMPLETED,
            duration_ns=10.0 / max(n_segments, 1),
        )
        segments.append(seg)

    return SimulationProgress(
        config_path="/tmp/config.yaml",
        total_steps_requested=total_steps,
        total_samples_requested=250,
        timestep_fs=2.0,
        segments=segments,
        status=SimulationStatus.COMPLETED if is_complete else SimulationStatus.INTERRUPTED,
        replicate=1,
    )


class TestRecoverStatusReport:
    """recover CLI shows progress status."""

    @patch("polyzymd.simulation.progress.save_progress")
    @patch("polyzymd.simulation.progress.load_or_scan_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_shows_progress_report(self, mock_from_yaml, mock_load, mock_save, tmp_path):
        """recover should report steps completed and remaining."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()

        sim_config = _mock_sim_config(working_dir)
        mock_from_yaml.return_value = sim_config
        mock_load.return_value = _mock_progress(
            total_steps=10000000, completed_steps=5000000, n_segments=1
        )

        runner = CliRunner()
        result = runner.invoke(cli, ["recover", "-c", str(config_file), "-r", "1"])
        assert result.exit_code == 0
        assert "5000000/10000000" in result.output
        assert "50.0%" in result.output

    @patch("polyzymd.simulation.progress.save_progress")
    @patch("polyzymd.simulation.progress.load_or_scan_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_complete_simulation_shows_nothing_to_recover(
        self, mock_from_yaml, mock_load, mock_save, tmp_path
    ):
        """recover should show 'nothing to recover' when simulation is done."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()

        sim_config = _mock_sim_config(working_dir)
        mock_from_yaml.return_value = sim_config
        mock_load.return_value = _mock_progress(
            total_steps=10000000, completed_steps=10000000, is_complete=True, n_segments=2
        )

        runner = CliRunner()
        result = runner.invoke(cli, ["recover", "-c", str(config_file), "-r", "1"])
        assert result.exit_code == 0
        assert "nothing to recover" in result.output.lower()

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_missing_working_dir_exits_with_error(self, mock_from_yaml, tmp_path):
        """recover should fail when working directory does not exist."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")

        mock_from_yaml.return_value = _mock_sim_config(tmp_path / "nonexistent")

        runner = CliRunner()
        result = runner.invoke(cli, ["recover", "-c", str(config_file), "-r", "1"])
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# check-progress CLI
# ---------------------------------------------------------------------------


class TestCheckProgress:
    """check-progress reports status and exits with correct code."""

    @patch("polyzymd.simulation.progress.load_or_scan_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_complete_exits_zero(self, mock_from_yaml, mock_load, tmp_path):
        """check-progress should exit 0 when simulation is complete."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()

        sim_config = _mock_sim_config(working_dir)
        mock_from_yaml.return_value = sim_config
        mock_load.return_value = _mock_progress(
            total_steps=10000000, completed_steps=10000000, is_complete=True, n_segments=2
        )

        runner = CliRunner()
        result = runner.invoke(cli, ["check-progress", "-c", str(config_file), "-r", "1"])
        assert result.exit_code == 0
        assert "COMPLETE" in result.output

    @patch("polyzymd.simulation.progress.load_or_scan_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_incomplete_exits_one(self, mock_from_yaml, mock_load, tmp_path):
        """check-progress should exit 1 when work remains."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()

        sim_config = _mock_sim_config(working_dir)
        mock_from_yaml.return_value = sim_config
        mock_load.return_value = _mock_progress(
            total_steps=10000000, completed_steps=5000000, n_segments=1
        )

        runner = CliRunner()
        result = runner.invoke(cli, ["check-progress", "-c", str(config_file), "-r", "1"])
        assert result.exit_code == 1
        assert "remaining" in result.output.lower()

    @patch("polyzymd.simulation.progress.load_or_scan_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_shows_segment_count(self, mock_from_yaml, mock_load, tmp_path):
        """check-progress should report the number of segments."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()

        sim_config = _mock_sim_config(working_dir)
        mock_from_yaml.return_value = sim_config
        mock_load.return_value = _mock_progress(
            total_steps=10000000, completed_steps=5000000, n_segments=3
        )

        runner = CliRunner()
        result = runner.invoke(cli, ["check-progress", "-c", str(config_file), "-r", "1"])
        assert "3 segment(s)" in result.output


# ---------------------------------------------------------------------------
# Self-resubmitting model (no dependency chains)
# ---------------------------------------------------------------------------


class TestSelfResubmittingModel:
    """DaisyChainSubmitter uses self-resubmitting jobs, not afterany chains."""

    def test_no_afterany_in_submit_job(self):
        """Verify _submit_job does NOT use afterany dependencies."""
        import inspect

        from polyzymd.workflow.daisy_chain import DaisyChainSubmitter

        source = inspect.getsource(DaisyChainSubmitter._submit_job)
        assert "afterany" not in source
        assert "afterok" not in source

    def test_generate_job_script_produces_self_resubmitting(self):
        """The unified job template should resubmit via sbatch $SLURM_JOB_SCRIPT."""
        from polyzymd.workflow.slurm import SlurmConfig, SlurmScriptGenerator

        gen = SlurmScriptGenerator(SlurmConfig.from_preset("aa100"), pixi_env="cuda-12-4")
        script = gen.generate_job_script(
            config_path="/tmp/config.yaml",
            replicate=1,
            working_dir="/tmp/work",
        )
        assert "SLURM_JOB_SCRIPT" in script
        assert "sbatch" in script
        assert "check-progress" in script
        assert "run-segment" in script
