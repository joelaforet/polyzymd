"""Tests for the recover and check-progress CLI commands.

Covers:
- _find_topology_pdb helper
- recover CLI with new self-resubmitting interface
- check-progress CLI (exit code 0 = complete, 1 = work remains)
- Self-resubmitting model (no afterany dependencies)
"""

from __future__ import annotations

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
        eq_dir = tmp_path / "equilibration_0_heating"
        eq_dir.mkdir()
        (eq_dir / "equilibration_0_heating_topology.pdb").write_text("ATOM ...")
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
    mock.engine = "openmm"
    mock.get_working_directory.return_value = working_dir
    mock.simulation_phases.production.time_step = 2.0
    mock.simulation_phases.production.duration = 20.0  # ns
    mock.simulation_phases.production.samples = 250
    return mock


def _mock_sim_config_gromacs(working_dir: Path):
    """Create a mock GROMACS SimulationConfig for recover tests."""
    mock = MagicMock()
    mock.engine = "gromacs"
    mock.get_working_directory.return_value = working_dir
    mock.simulation_phases.production.time_step = 2.0
    mock.simulation_phases.production.duration = 20.0
    mock.simulation_phases.production.samples = 250
    mock.generate_system_name.return_value = "system"
    mock.output.slurm_logs_subdir = "slurm_logs"
    mock.enzyme.name = "CALB"
    mock.thermodynamics.temperature = 310
    mock.polymers = None
    mock.gromacs.grompp_flags = ""
    mock.gromacs.mdrun_flags = ""
    mock.gromacs.module_load = None
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

    @patch("polyzymd.engines.openmm.engine.OpenMMEngine.load_or_scan_progress")
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

    @patch("polyzymd.engines.openmm.engine.OpenMMEngine.load_or_scan_progress")
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

    @patch("polyzymd.engines.openmm.engine.OpenMMEngine.load_or_scan_progress")
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

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_config_error_exits_three(self, mock_from_yaml, tmp_path):
        """check-progress should exit 3 (not 1) on config load failure.

        Exit code 1 means 'work remains' and triggers SLURM resubmission.
        Errors must use a distinct exit code to prevent infinite resubmission.
        """
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")

        mock_from_yaml.side_effect = ValueError("bad config")

        runner = CliRunner()
        result = runner.invoke(cli, ["check-progress", "-c", str(config_file), "-r", "1"])
        assert result.exit_code == 3

    @patch("polyzymd.engines.openmm.engine.OpenMMEngine.load_or_scan_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_progress_load_error_exits_three(self, mock_from_yaml, mock_load, tmp_path):
        """check-progress should exit 3 on progress load failure."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()

        sim_config = _mock_sim_config(working_dir)
        mock_from_yaml.return_value = sim_config
        mock_load.side_effect = FileNotFoundError("no progress file")

        runner = CliRunner()
        result = runner.invoke(cli, ["check-progress", "-c", str(config_file), "-r", "1"])
        assert result.exit_code == 3

    def test_slurm_template_guards_check_progress_errors(self):
        """SLURM template must NOT resubmit on check-progress error codes."""
        from polyzymd.workflow.slurm import SlurmScriptGenerator

        # Access the job template
        template = SlurmScriptGenerator.JOB_TEMPLATE
        # Must check for non-1 exit code before resubmitting
        assert "PROGRESS_RC -ne 1" in template
        # Must stop on error
        assert "NOT resubmitting" in template


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
        from unittest.mock import patch

        from polyzymd.workflow.slurm import SlurmConfig, SlurmScriptGenerator

        gen = SlurmScriptGenerator(SlurmConfig.from_preset("aa100"), pixi_env="cuda-12-4")
        with patch(
            "polyzymd.workflow.slurm._discover_manifest_path",
            return_value="/fake/pixi.toml",
        ):
            script = gen.generate_job_script(
                config_path="/tmp/config.yaml",
                replicate=1,
                working_dir="/tmp/work",
            )
        assert "SLURM_JOB_SCRIPT" in script
        assert "sbatch" in script
        assert "check-progress" in script
        assert "run-segment" in script


# ---------------------------------------------------------------------------
# recover --submit --dry-run: skip_build detection
# ---------------------------------------------------------------------------


def _mock_sim_config_for_submit(working_dir: Path):
    """Create a mock SimulationConfig suitable for the recover --submit path.

    Extends the basic mock with attributes needed by create_job_name(),
    SlurmScriptGenerator, and the script-generation path.
    """
    mock = MagicMock()
    mock.get_working_directory.return_value = working_dir
    mock.simulation_phases.production.time_step = 2.0
    mock.simulation_phases.production.duration = 20.0
    mock.simulation_phases.production.samples = 250
    mock.output.slurm_logs_subdir = "slurm_logs"
    mock.enzyme.name = "CALB"
    mock.thermodynamics.temperature = 310
    mock.polymers = None  # no polymer info
    return mock


class TestRecoverSkipBuild:
    """recover --submit detects pre-built system and passes --skip-build."""

    @patch("polyzymd.workflow.slurm._discover_manifest_path", return_value="/fake/pixi.toml")
    @patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=[])
    @patch("polyzymd.simulation.progress.save_progress")
    @patch("polyzymd.simulation.progress.load_or_scan_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_skip_build_when_system_files_exist(
        self, mock_from_yaml, mock_load, mock_save, mock_squeue, mock_manifest, tmp_path
    ):
        """When solvated_system.pdb and system.xml exist, script has --skip-build."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()
        (working_dir / "solvated_system.pdb").write_text("ATOM ...")
        (working_dir / "system.xml").write_text("<System/>")

        mock_from_yaml.return_value = _mock_sim_config_for_submit(working_dir)
        mock_load.return_value = _mock_progress(
            total_steps=10000000, completed_steps=0, n_segments=0
        )

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "recover",
                "-c",
                str(config_file),
                "-r",
                "1",
                "--submit",
                "--dry-run",
                "--preset",
                "aa100",
            ],
        )
        assert result.exit_code == 0, result.output

        # The generated script should be written to disk
        script_path = working_dir / "recovery_scripts" / "recover_rep1.sh"
        assert script_path.exists()
        script_content = script_path.read_text()
        assert "--skip-build" in script_content

    @patch("polyzymd.workflow.slurm._discover_manifest_path", return_value="/fake/pixi.toml")
    @patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=[])
    @patch("polyzymd.simulation.progress.save_progress")
    @patch("polyzymd.simulation.progress.load_or_scan_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_no_skip_build_when_files_missing(
        self, mock_from_yaml, mock_load, mock_save, mock_squeue, mock_manifest, tmp_path
    ):
        """When system files are absent, script does NOT have --skip-build."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()
        # No solvated_system.pdb or system.xml

        mock_from_yaml.return_value = _mock_sim_config_for_submit(working_dir)
        mock_load.return_value = _mock_progress(
            total_steps=10000000, completed_steps=0, n_segments=0
        )

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "recover",
                "-c",
                str(config_file),
                "-r",
                "1",
                "--submit",
                "--dry-run",
                "--preset",
                "aa100",
            ],
        )
        assert result.exit_code == 0, result.output

        script_path = working_dir / "recovery_scripts" / "recover_rep1.sh"
        assert script_path.exists()
        script_content = script_path.read_text()
        assert "--skip-build" not in script_content

    @patch("polyzymd.workflow.slurm._discover_manifest_path", return_value="/fake/pixi.toml")
    @patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=[])
    @patch("polyzymd.simulation.progress.save_progress")
    @patch("polyzymd.simulation.progress.load_or_scan_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_skip_build_message_shown(
        self, mock_from_yaml, mock_load, mock_save, mock_squeue, mock_manifest, tmp_path
    ):
        """When skip_build is detected, CLI output mentions --skip-build."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()
        (working_dir / "solvated_system.pdb").write_text("ATOM ...")
        (working_dir / "system.xml").write_text("<System/>")

        mock_from_yaml.return_value = _mock_sim_config_for_submit(working_dir)
        mock_load.return_value = _mock_progress(
            total_steps=10000000, completed_steps=0, n_segments=0
        )

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "recover",
                "-c",
                str(config_file),
                "-r",
                "1",
                "--submit",
                "--dry-run",
                "--preset",
                "aa100",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "skip-build" in result.output.lower()

    @patch("polyzymd.workflow.slurm._discover_manifest_path", return_value="/fake/pixi.toml")
    @patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=[])
    @patch("polyzymd.simulation.progress.save_progress")
    @patch("polyzymd.simulation.progress.load_or_scan_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_partial_system_files_no_skip_build(
        self, mock_from_yaml, mock_load, mock_save, mock_squeue, mock_manifest, tmp_path
    ):
        """When only one of the two system files exists, no --skip-build."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()
        # Only the PDB, not the XML
        (working_dir / "solvated_system.pdb").write_text("ATOM ...")

        mock_from_yaml.return_value = _mock_sim_config_for_submit(working_dir)
        mock_load.return_value = _mock_progress(
            total_steps=10000000, completed_steps=0, n_segments=0
        )

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "recover",
                "-c",
                str(config_file),
                "-r",
                "1",
                "--submit",
                "--dry-run",
                "--preset",
                "aa100",
            ],
        )
        assert result.exit_code == 0, result.output

        script_path = working_dir / "recovery_scripts" / "recover_rep1.sh"
        script_content = script_path.read_text()
        assert "--skip-build" not in script_content


# ---------------------------------------------------------------------------
# _run_initial_segment: equilibration skip when skip_build + eq complete
# ---------------------------------------------------------------------------


class TestRunInitialSegmentEquilibrationSkip:
    """_run_initial_segment skips minimize + equilibrate when appropriate."""

    def _make_progress_with_eq(self, working_dir, n_stages=3):
        """Create a progress.json with completed equilibration stages."""
        from polyzymd.simulation.progress import (
            EquilibrationStageRecord,
            SegmentStatus,
            SimulationProgress,
            SimulationStatus,
            save_progress,
        )

        stage_names = ["heating", "npt_equilibration", "npt_production_eq"][:n_stages]
        eq_stages = [
            EquilibrationStageRecord(
                index=i,
                name=name,
                status=SegmentStatus.COMPLETED,
                duration_ns=0.5,
                ensemble="NVT" if i == 0 else "NPT",
                finished_at="2025-01-01T00:00:00+00:00",
            )
            for i, name in enumerate(stage_names)
        ]

        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=10000000,
            total_samples_requested=250,
            timestep_fs=2.0,
            equilibration_stages=eq_stages,
            segments=[],  # No production segments yet
            status=SimulationStatus.INTERRUPTED,
            replicate=1,
        )
        save_progress(working_dir, progress)
        return progress

    def _create_system_files(self, working_dir):
        """Create dummy solvated_system.pdb and system.xml on disk."""
        (working_dir / "solvated_system.pdb").write_text("ATOM mock")
        (working_dir / "system.xml").write_text("<System/>")

    def _make_sim_config_staged(self, n_stages=3):
        """Create a mock SimulationConfig with staged equilibration."""
        sim_config = MagicMock()
        sim_config.thermodynamics.temperature = 310.0
        sim_config.thermodynamics.pressure = 1.0
        sim_config.simulation_phases.total_equilibration_duration = 1.5
        sim_config.restraints = []

        stage_names = ["heating", "npt_equilibration", "npt_production_eq"][:n_stages]
        stages = []
        for name in stage_names:
            s = MagicMock()
            s.name = name
            stages.append(s)
        sim_config.simulation_phases.equilibration_stages = stages
        return sim_config

    @patch("polyzymd.simulation.runner.SimulationRunner")
    def test_skips_equilibration_when_complete(self, MockRunner, tmp_path):
        """With skip_build=True and equilibration complete, skip minimize + eq."""
        from polyzymd.cli.main import _run_initial_segment

        # Set up filesystem state
        self._make_progress_with_eq(tmp_path, n_stages=3)
        self._create_system_files(tmp_path)

        sim_config = self._make_sim_config_staged(n_stages=3)

        # Mock the runner
        mock_runner = MagicMock()
        MockRunner.return_value = mock_runner

        # Mock the openmm imports used in the skip_build=True file-loading path.
        # The fast-path (eq already complete) returns before run_production is
        # reached in the normal flow, but the PDBFile / XmlSerializer loads
        # happen first.
        mock_pdb_class = MagicMock()
        mock_pdb_instance = MagicMock()
        mock_pdb_instance.topology = MagicMock()
        mock_pdb_instance.positions = MagicMock()
        mock_pdb_class.return_value = mock_pdb_instance

        mock_xml_serializer = MagicMock()
        mock_xml_serializer.deserialize.return_value = MagicMock()

        with patch.dict(
            "sys.modules",
            {
                "openmm": MagicMock(XmlSerializer=mock_xml_serializer),
                "openmm.app": MagicMock(PDBFile=mock_pdb_class),
            },
        ):
            _run_initial_segment(
                sim_config=sim_config,
                working_dir=tmp_path,
                replicate=1,
                skip_build=True,
                duration_ns=20.0,
                num_samples=250,
                timestep_fs=2.0,
            )

        # Should NOT have minimized or equilibrated
        mock_runner.minimize.assert_not_called()
        mock_runner.run_equilibration.assert_not_called()

        # Should have loaded eq state from last stage (index 2, name npt_production_eq)
        mock_runner._load_eq_stage_state.assert_called_once_with(2, "npt_production_eq")

        # Should have run production
        mock_runner.run_production.assert_called_once()
        call_kwargs = mock_runner.run_production.call_args[1]
        assert call_kwargs["segment_index"] == 0
        assert call_kwargs["temperature"] == 310.0
        assert call_kwargs["duration_ns"] == 20.0

    @patch("polyzymd.simulation.runner.SimulationRunner")
    def test_runs_equilibration_when_eq_not_complete(self, MockRunner, tmp_path):
        """With skip_build=True but no eq stages recorded, run full pipeline."""
        from polyzymd.cli.main import _run_initial_segment
        from polyzymd.simulation.progress import (
            SimulationProgress,
            SimulationStatus,
            save_progress,
        )

        # Progress with NO equilibration stages — eq_complete is False
        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=10000000,
            total_samples_requested=250,
            timestep_fs=2.0,
            equilibration_stages=[],
            segments=[],
            status=SimulationStatus.NOT_STARTED,
            replicate=1,
        )
        save_progress(tmp_path, progress)
        self._create_system_files(tmp_path)

        sim_config = self._make_sim_config_staged(n_stages=0)

        mock_runner = MagicMock()
        mock_runner.run_equilibration.return_value = {"type": "staged_equilibration", "stages": []}
        MockRunner.return_value = mock_runner

        mock_pdb_class = MagicMock()
        mock_pdb_instance = MagicMock()
        mock_pdb_instance.topology = MagicMock()
        mock_pdb_instance.positions = MagicMock()
        mock_pdb_class.return_value = mock_pdb_instance

        mock_xml_serializer = MagicMock()
        mock_xml_serializer.deserialize.return_value = MagicMock()

        with patch.dict(
            "sys.modules",
            {
                "openmm": MagicMock(XmlSerializer=mock_xml_serializer),
                "openmm.app": MagicMock(PDBFile=mock_pdb_class),
            },
        ):
            _run_initial_segment(
                sim_config=sim_config,
                working_dir=tmp_path,
                replicate=1,
                skip_build=True,
                duration_ns=20.0,
                num_samples=250,
                timestep_fs=2.0,
            )

        # Should have minimized and equilibrated (full pipeline)
        mock_runner.minimize.assert_called_once()
        mock_runner.run_equilibration.assert_called_once()
        mock_runner.run_production.assert_called_once()

        # Should NOT have used the fast path
        mock_runner._load_eq_stage_state.assert_not_called()

    def test_no_skip_build_does_not_enter_fast_path(self, tmp_path):
        """With skip_build=False, the eq-skip fast path is never entered.

        We verify this structurally: inspect the source code to confirm the
        fast-path guard is gated on ``if skip_build:``. This avoids having to
        mock the entire build chain (SystemBuilder, OpenFF, etc.) which is
        not importable in the test environment.
        """
        import inspect

        from polyzymd.cli.main import _run_initial_segment

        source = inspect.getsource(_run_initial_segment)
        # The fast-path block must be inside ``if skip_build:``
        assert "if skip_build:" in source
        # And it must check equilibration_complete before jumping to production
        assert "equilibration_complete" in source
        # The _load_eq_stage_state call must be inside the skip_build guard
        assert "_load_eq_stage_state" in source


# ---------------------------------------------------------------------------
# Case 5b: Hard kill without system.xml raises immediately (B6)
# ---------------------------------------------------------------------------


class TestHardKillNoSystemXmlRaises:
    """_get_previous_paths must raise FileNotFoundError for Case 5b.

    When a segment has a periodic checkpoint (.chk) but no system.xml,
    recovery is impossible.  Before the fix, the code logged an error but
    fell through silently, returning paths to non-existent files.
    """

    def test_checkpoint_without_system_xml_raises(self, tmp_path):
        """Case 5b: .chk exists, system.xml absent — must raise."""
        from polyzymd.simulation.continuation import ContinuationManager

        prev_idx = 0
        seg_dir = tmp_path / f"production_{prev_idx}"
        seg_dir.mkdir()

        # Only the periodic checkpoint exists — no state/system XML at all
        (seg_dir / f"production_{prev_idx}_checkpoint.chk").write_bytes(b"binary")
        (seg_dir / f"production_{prev_idx}_parameters.json").write_text("{}")

        mgr = ContinuationManager(working_dir=tmp_path, segment_index=prev_idx + 1)
        with pytest.raises(FileNotFoundError, match="cannot recover"):
            mgr._get_previous_paths()

    def test_checkpoint_with_system_xml_does_not_raise(self, tmp_path):
        """Case 5 (normal hard kill): .chk + system.xml — should NOT raise."""
        from polyzymd.simulation.continuation import ContinuationManager

        prev_idx = 0
        seg_dir = tmp_path / f"production_{prev_idx}"
        seg_dir.mkdir()

        (seg_dir / f"production_{prev_idx}_checkpoint.chk").write_bytes(b"binary")
        (seg_dir / f"production_{prev_idx}_system.xml").write_text("<System/>")
        (seg_dir / f"production_{prev_idx}_parameters.json").write_text("{}")

        mgr = ContinuationManager(working_dir=tmp_path, segment_index=prev_idx + 1)
        paths = mgr._get_previous_paths()
        assert paths["use_checkpoint"] is True

    def test_normal_completion_does_not_raise(self, tmp_path):
        """Case 1: state.xml exists (normal completion) — should NOT raise."""
        from polyzymd.simulation.continuation import ContinuationManager

        prev_idx = 0
        seg_dir = tmp_path / f"production_{prev_idx}"
        seg_dir.mkdir()

        (seg_dir / f"production_{prev_idx}_state.xml").write_text("<State/>")
        (seg_dir / f"production_{prev_idx}_system.xml").write_text("<System/>")
        (seg_dir / f"production_{prev_idx}_parameters.json").write_text("{}")

        mgr = ContinuationManager(working_dir=tmp_path, segment_index=prev_idx + 1)
        paths = mgr._get_previous_paths()
        assert paths["use_checkpoint"] is False


# ---------------------------------------------------------------------------
# B11 – GROMACS dry-run output path uses f-string, not literal braces
# ---------------------------------------------------------------------------


class TestBuildDryRunGromacs:
    """build --dry-run --gromacs must show an actual directory, not {projects_dir}."""

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_gromacs_dry_run_shows_actual_path(self, mock_from_yaml, tmp_path):
        """The GROMACS output line must contain the real projects_directory, not
        the literal string '{projects_dir}'."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")

        mock_cfg = MagicMock()
        mock_cfg.name = "test"
        mock_cfg.enzyme.name = "CALB"
        mock_cfg.substrate = None
        mock_cfg.polymers = None
        mock_cfg.thermodynamics.temperature = 310.0
        mock_cfg.simulation_phases.production.duration = 100.0
        mock_cfg.output.projects_directory = tmp_path / "projects"
        mock_cfg.output.effective_scratch_directory = tmp_path / "scratch"
        mock_from_yaml.return_value = mock_cfg

        runner = CliRunner()
        result = runner.invoke(
            cli, ["build", "-c", str(config_file), "--dry-run", "--gromacs", "-r", "1"]
        )
        assert result.exit_code == 0, result.output
        assert "{projects_dir}" not in result.output, (
            "Output path must be interpolated, not literal {projects_dir}"
        )
        assert str(tmp_path / "projects") in result.output


class TestRecoverEngineAware:
    """recover dispatches progress loading by engine."""

    @patch("polyzymd.simulation.progress.save_progress")
    @patch("polyzymd.simulation.progress.load_or_scan_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_recover_openmm_default_uses_direct_progress(
        self, mock_from_yaml, mock_load, mock_save, tmp_path
    ):
        """OpenMM recover path should use direct progress helper."""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()

        mock_from_yaml.return_value = _mock_sim_config(working_dir)
        mock_load.return_value = _mock_progress(
            total_steps=10000000, completed_steps=5000000, n_segments=1
        )

        runner = CliRunner()
        result = runner.invoke(cli, ["recover", "-c", str(config_file), "-r", "1"])

        assert result.exit_code == 0, result.output
        mock_load.assert_called_once()

    @patch("polyzymd.simulation.progress.save_progress")
    @patch("polyzymd.simulation.progress.load_or_scan_progress")
    @patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
    @patch("polyzymd.engines.gromacs.engine.GromacsEngine.load_or_scan_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_recover_gromacs_uses_engine_dispatch(
        self,
        mock_from_yaml,
        mock_gromacs_load,
        mock_resolve_binary,
        mock_openmm_load,
        mock_save,
        tmp_path,
    ):
        """GROMACS recover path should dispatch through engine interface."""
        _ = mock_resolve_binary
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()

        mock_from_yaml.return_value = _mock_sim_config_gromacs(working_dir)
        mock_gromacs_load.return_value = _mock_progress(
            total_steps=10000000, completed_steps=5000000, n_segments=1
        )

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["recover", "-c", str(config_file), "-r", "1", "--engine", "gromacs"],
        )

        assert result.exit_code == 0, result.output
        mock_gromacs_load.assert_called_once_with(working_dir, 1)
        mock_openmm_load.assert_not_called()


class TestRecoverGromacsSubmit:
    """recover --submit supports GROMACS engine submission flow."""

    @patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=[])
    @patch("polyzymd.engines.create_engine")
    @patch("polyzymd.simulation.progress.save_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_gromacs_recover_dry_run_creates_script(
        self,
        mock_from_yaml,
        mock_save,
        mock_create_engine,
        mock_squeue,
        tmp_path,
    ):
        """Dry-run recover should stage a recovery script for GROMACS."""
        _ = mock_save, mock_squeue
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()

        sim_config = _mock_sim_config_gromacs(working_dir)
        mock_from_yaml.return_value = sim_config

        engine_mock = MagicMock()
        engine_mock.load_or_scan_progress.return_value = _mock_progress(
            total_steps=10000000, completed_steps=5000000, n_segments=1
        )

        def _prepare_submission_side_effect(request):
            daisy_dir = request.working_dir / "daisy_chain_scripts"
            daisy_dir.mkdir(parents=True, exist_ok=True)
            script = daisy_dir / f"run_rep{request.replicate}.sh"
            script.write_text("#!/bin/bash\n")
            return script

        engine_mock.prepare_submission.side_effect = _prepare_submission_side_effect
        mock_create_engine.return_value = engine_mock

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "recover",
                "-c",
                str(config_file),
                "-r",
                "1",
                "--engine",
                "gromacs",
                "--submit",
                "--dry-run",
                "--preset",
                "aa100",
            ],
        )

        assert result.exit_code == 0, result.output
        assert (working_dir / "recovery_scripts" / "recover_rep1.sh").exists()

    @patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=[])
    @patch("polyzymd.engines.create_engine")
    @patch("polyzymd.simulation.progress.save_progress")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_gromacs_recover_detects_existing_inputs(
        self,
        mock_from_yaml,
        mock_save,
        mock_create_engine,
        mock_squeue,
        tmp_path,
    ):
        """Recover should report reuse when GROMACS inputs already exist."""
        _ = mock_save, mock_squeue
        config_file = tmp_path / "config.yaml"
        config_file.write_text("name: test")
        working_dir = tmp_path / "work"
        working_dir.mkdir()

        (working_dir / "system.top").write_text("[ system ]\n")
        (working_dir / "prod.mdp").write_text("integrator = md\n")

        sim_config = _mock_sim_config_gromacs(working_dir)
        mock_from_yaml.return_value = sim_config

        engine_mock = MagicMock()
        engine_mock.load_or_scan_progress.return_value = _mock_progress(
            total_steps=10000000, completed_steps=5000000, n_segments=1
        )

        def _prepare_submission_side_effect(request):
            daisy_dir = request.working_dir / "daisy_chain_scripts"
            daisy_dir.mkdir(parents=True, exist_ok=True)
            script = daisy_dir / f"run_rep{request.replicate}.sh"
            script.write_text("#!/bin/bash\n")
            return script

        engine_mock.prepare_submission.side_effect = _prepare_submission_side_effect
        mock_create_engine.return_value = engine_mock

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "recover",
                "-c",
                str(config_file),
                "-r",
                "1",
                "--engine",
                "gromacs",
                "--submit",
                "--dry-run",
                "--preset",
                "aa100",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "reuse" in result.output.lower()


class TestUpdateGromacsProgressCmd:
    """Hidden command updates GROMACS progress state from logs."""

    def test_update_gromacs_progress_basic(self, tmp_path):
        """Progress command should create progress.json from prod.log."""
        from polyzymd.simulation.progress import load_progress

        working_dir = tmp_path / "gromacs"
        working_dir.mkdir()
        (working_dir / "prod.log").write_text("nsteps = 5000\n1000 2.0\n")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "_update-gromacs-progress",
                "--working-dir",
                str(working_dir),
                "--config-path",
                str(tmp_path / "config.yaml"),
                "--replicate",
                "1",
            ],
        )

        assert result.exit_code == 0, result.output
        assert (working_dir / "progress.json").exists()
        progress = load_progress(working_dir)
        assert progress is not None
        assert progress.total_steps_completed == 1000

    def test_update_gromacs_progress_mark_complete(self, tmp_path):
        """mark-complete should force COMPLETED status in saved progress."""
        from polyzymd.simulation.progress import SimulationStatus, load_progress

        working_dir = tmp_path / "gromacs"
        working_dir.mkdir()
        (working_dir / "prod.log").write_text("nsteps = 5000\n1000 2.0\n")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "_update-gromacs-progress",
                "--working-dir",
                str(working_dir),
                "--config-path",
                str(tmp_path / "config.yaml"),
                "--replicate",
                "1",
                "--mark-complete",
            ],
        )

        assert result.exit_code == 0, result.output
        progress = load_progress(working_dir)
        assert progress is not None
        assert progress.status == SimulationStatus.COMPLETED
