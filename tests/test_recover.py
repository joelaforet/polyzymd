"""Tests for the recover and recover-chain CLI commands.

Covers:
- INTERRUPTED marker parsing
- recover --dry-run output
- recover-chain directory scanning and status reporting
- _find_topology_pdb helper
- SLURM template recovery preamble content
- afterany dependency in daisy chain
"""

from pathlib import Path

import pytest
from click.testing import CliRunner

from polyzymd.cli.main import _find_topology_pdb, cli

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_interrupted_marker(seg_dir: Path, steps_completed: int, total_steps: int) -> None:
    """Write an INTERRUPTED marker file to a segment directory."""
    seg_dir.mkdir(parents=True, exist_ok=True)
    remaining = total_steps - steps_completed
    (seg_dir / "INTERRUPTED").write_text(
        f"segment_index={seg_dir.name.split('_')[1]}\n"
        f"steps_completed={steps_completed}\n"
        f"total_steps={total_steps}\n"
        f"remaining_steps={remaining}\n"
    )


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
# recover --dry-run
# ---------------------------------------------------------------------------


class TestRecoverDryRun:
    """recover --dry-run reads the INTERRUPTED marker and exits without running."""

    def test_dry_run_shows_steps(self, tmp_path):
        seg_dir = tmp_path / "production_3"
        _write_interrupted_marker(seg_dir, steps_completed=50000, total_steps=100000)
        # Write parameters file (needed for dry run)
        import json

        params = {
            "__values__": {
                "integ_params": {
                    "__values__": {
                        "num_samples": 250,
                        "time_step": {
                            "__class__": "Quantity",
                            "__values__": {"value": 2.0, "unit": "femtosecond"},
                        },
                        "total_time": {
                            "__class__": "Quantity",
                            "__values__": {"value": 20.0, "unit": "nanosecond"},
                        },
                    }
                },
                "thermo_params": {"__values__": {}},
            }
        }
        (seg_dir / "production_3_parameters.json").write_text(json.dumps(params))

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["recover", "-w", str(tmp_path), "-s", "3", "--dry-run"],
        )
        assert result.exit_code == 0
        assert "50000/100000" in result.output
        assert "DRY RUN" in result.output

    def test_no_marker_exits_with_error(self, tmp_path):
        seg_dir = tmp_path / "production_3"
        seg_dir.mkdir(parents=True)
        # No INTERRUPTED marker

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["recover", "-w", str(tmp_path), "-s", "3"],
        )
        assert result.exit_code != 0
        assert "No INTERRUPTED marker" in result.output


# ---------------------------------------------------------------------------
# recover-chain
# ---------------------------------------------------------------------------


class TestRecoverChain:
    """recover-chain scans production directories and reports status."""

    def test_all_completed(self, tmp_path):
        for i in range(3):
            _write_completed_segment(tmp_path / f"production_{i}", i)

        runner = CliRunner()
        result = runner.invoke(cli, ["recover-chain", "-w", str(tmp_path), "--dry-run"])
        assert result.exit_code == 0
        assert "3 completed" in result.output
        assert "0 interrupted" in result.output
        assert "All segments healthy" in result.output

    def test_one_interrupted(self, tmp_path):
        _write_completed_segment(tmp_path / "production_0", 0)
        _write_completed_segment(tmp_path / "production_1", 1)
        seg2 = tmp_path / "production_2"
        _write_interrupted_marker(seg2, 30000, 100000)

        runner = CliRunner()
        result = runner.invoke(cli, ["recover-chain", "-w", str(tmp_path), "--dry-run"])
        assert result.exit_code == 0
        assert "1 interrupted" in result.output
        assert "INTERRUPTED" in result.output
        assert "DRY RUN" in result.output

    def test_missing_segment(self, tmp_path):
        _write_completed_segment(tmp_path / "production_0", 0)
        # production_1 exists but has no state.xml and no INTERRUPTED marker
        (tmp_path / "production_1").mkdir()

        runner = CliRunner()
        result = runner.invoke(cli, ["recover-chain", "-w", str(tmp_path), "--dry-run"])
        assert result.exit_code == 0
        assert "1 missing" in result.output
        assert "MISSING" in result.output

    def test_no_segments(self, tmp_path):
        runner = CliRunner()
        result = runner.invoke(cli, ["recover-chain", "-w", str(tmp_path), "--dry-run"])
        assert result.exit_code == 0
        assert "No production segments found" in result.output

    def test_resume_dir_noted(self, tmp_path):
        _write_completed_segment(tmp_path / "production_0", 0)
        # Create a resume directory
        (tmp_path / "production_0_resume").mkdir()

        runner = CliRunner()
        result = runner.invoke(cli, ["recover-chain", "-w", str(tmp_path), "--dry-run"])
        assert result.exit_code == 0
        assert "has resume/" in result.output


# ---------------------------------------------------------------------------
# SLURM template content
# ---------------------------------------------------------------------------


class TestSlurmTemplateRecovery:
    """SLURM templates include signal and recovery directives."""

    def test_signal_directive_in_initial_template(self):
        from polyzymd.workflow.slurm import JobContext, SlurmConfig, SlurmScriptGenerator

        gen = SlurmScriptGenerator(SlurmConfig.from_preset("aa100"), conda_env="test-env")
        ctx = JobContext(
            job_name="test",
            output_file="test.out",
            scratch_dir="/scratch/test",
            projects_dir="/projects/test",
        )
        script = gen.generate_initial_job(
            context=ctx,
            config_path="config.yaml",
            replicate=1,
            segment_time=20.0,
            segment_frames=250,
        )
        assert "#SBATCH --signal=B:USR1@300" in script
        assert "#SBATCH --no-requeue" in script

    def test_signal_directive_in_continuation_template(self):
        from polyzymd.workflow.slurm import JobContext, SlurmConfig, SlurmScriptGenerator

        gen = SlurmScriptGenerator(SlurmConfig.from_preset("aa100"), conda_env="test-env")
        ctx = JobContext(
            job_name="test",
            output_file="test.out",
            scratch_dir="/scratch/test",
            projects_dir="/projects/test",
            segment_index=3,
        )
        script = gen.generate_continuation_job(
            context=ctx,
            segment_time=20.0,
            num_samples=250,
        )
        assert "#SBATCH --signal=B:USR1@300" in script
        assert "#SBATCH --no-requeue" in script

    def test_recovery_preamble_in_continuation_template(self):
        from polyzymd.workflow.slurm import JobContext, SlurmConfig, SlurmScriptGenerator

        gen = SlurmScriptGenerator(SlurmConfig.from_preset("aa100"), conda_env="test-env")
        ctx = JobContext(
            job_name="test",
            output_file="test.out",
            scratch_dir="/scratch/test",
            projects_dir="/projects/test",
            segment_index=3,
        )
        script = gen.generate_continuation_job(
            context=ctx,
            segment_time=20.0,
            num_samples=250,
        )
        assert "Recovery preamble" in script
        assert "INTERRUPTED" in script
        assert "polyzymd recover" in script


# ---------------------------------------------------------------------------
# afterany dependency
# ---------------------------------------------------------------------------


class TestAfteranyDependency:
    """DaisyChainSubmitter uses afterany instead of afterok."""

    def test_afterany_in_source(self):
        """Verify the dependency string uses afterany."""
        import inspect

        from polyzymd.workflow.daisy_chain import DaisyChainSubmitter

        source = inspect.getsource(DaisyChainSubmitter._submit_job)
        assert "afterany" in source
        assert "afterok" not in source
