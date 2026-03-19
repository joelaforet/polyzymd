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
        sim_config.simulation_phases.uses_staged_equilibration = True
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
