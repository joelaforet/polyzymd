"""Tests for the ``polyzymd cancel`` command and the STOP marker.

A self-resubmitting chain cannot be stopped with ``scancel`` alone: SIGTERM
makes ``run-segment`` exit 99, which the job wrapper reads as "interrupted,
work remains", so it queues a successor within seconds.  ``polyzymd cancel``
writes a STOP marker that the wrapper honours, then cancels the jobs.

Covers:
- STOP marker contents and location
- job lookup and cancellation by run directory or job name
- --resume, --stop-only and --dry-run
- cancel_slurm_jobs best-effort behaviour outside SLURM
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from polyzymd.cli.main import cli, stop_file_path, write_stop_file

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _config_and_mock(tmp_path: Path, working_dir_for):
    """Return (config path, mock SimulationConfig) wired to *working_dir_for*."""
    config_path = tmp_path / "config.yaml"
    config_path.write_text("engine: openmm\n")
    sim_config = MagicMock()
    sim_config.engine = "openmm"
    sim_config.get_working_directory.side_effect = working_dir_for
    return config_path, sim_config


def _invoke(config_path, sim_config, args, job_ids=("4242",)):
    runner = CliRunner()
    with (
        patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=sim_config),
        patch("polyzymd.workflow.daisy_chain.create_job_name", return_value="pzmd_r1"),
        patch(
            "polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=list(job_ids)
        ) as lookup,
        patch(
            "polyzymd.workflow.daisy_chain.cancel_slurm_jobs", side_effect=lambda ids: list(ids)
        ) as cancel_jobs,
    ):
        result = runner.invoke(cli, ["cancel", "-c", str(config_path), *args])
    return result, lookup, cancel_jobs


# ---------------------------------------------------------------------------
# STOP marker
# ---------------------------------------------------------------------------


class TestStopMarker:
    """The STOP marker is human-obvious and says who wrote it and when."""

    def test_written_at_the_root_of_the_working_directory(self, tmp_path):
        marker = write_stop_file(tmp_path / "run_1", "/projects/run/config.yaml", 1)
        assert marker == stop_file_path(tmp_path / "run_1")
        assert marker.name == "STOP"

    def test_contents_identify_author_time_and_run(self, tmp_path):
        marker = write_stop_file(tmp_path / "run_2", "/projects/run/config.yaml", 2)
        text = marker.read_text()
        assert "PolyzyMD STOP marker" in text
        assert "written_by:" in text
        assert "written_at:" in text
        assert "/projects/run/config.yaml" in text
        assert "replicate:    2" in text
        assert "polyzymd cancel --resume" in text

    def test_creates_a_missing_working_directory(self, tmp_path):
        """A chain can be stopped before its scratch directory exists."""
        marker = write_stop_file(tmp_path / "not_yet" / "run_1", "config.yaml", 1)
        assert marker.exists()


# ---------------------------------------------------------------------------
# polyzymd cancel
# ---------------------------------------------------------------------------


class TestCancelCommand:
    """``polyzymd cancel`` stops a chain and can hand it back."""

    def test_writes_marker_and_cancels_matching_jobs(self, tmp_path):
        working_dir = tmp_path / "run_1"
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: working_dir)

        result, lookup, cancel_jobs = _invoke(config_path, sim_config, ["-r", "1"])

        assert result.exit_code == 0, result.output
        assert stop_file_path(working_dir).exists()
        lookup.assert_called_once_with(working_dir, "pzmd_r1")
        cancel_jobs.assert_called_once_with(["4242"])
        assert "4242" in result.output

    def test_cancels_a_chain_that_works_in_the_folder_it_was_submitted_from(self, tmp_path):
        """Chains from older versions work where they were submitted; their name matches."""
        working_dir = tmp_path / "run_1"
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: working_dir)
        squeue = MagicMock(returncode=0, stdout=f"4242|pzmd_r1|{tmp_path / 'projects'}\n")

        with (
            patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=sim_config),
            patch("polyzymd.workflow.daisy_chain.create_job_name", return_value="pzmd_r1"),
            patch("polyzymd.workflow.daisy_chain.subprocess.run", return_value=squeue),
            patch(
                "polyzymd.workflow.daisy_chain.cancel_slurm_jobs",
                side_effect=lambda ids: list(ids),
            ) as cancel_jobs,
        ):
            result = CliRunner().invoke(cli, ["cancel", "-c", str(config_path), "-r", "1"])

        assert result.exit_code == 0, result.output
        cancel_jobs.assert_called_once_with(["4242"])

    def test_handles_several_replicates(self, tmp_path):
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: tmp_path / f"run_{r}")

        result, _, cancel_jobs = _invoke(config_path, sim_config, ["-r", "1-3"])

        assert result.exit_code == 0, result.output
        for replicate in (1, 2, 3):
            assert stop_file_path(tmp_path / f"run_{replicate}").exists()
        assert cancel_jobs.call_count == 3

    def test_scratch_dir_override_is_honoured(self, tmp_path):
        override = tmp_path / "elsewhere"
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: tmp_path / "unused")

        result, _, _ = _invoke(config_path, sim_config, ["-r", "1", "--scratch-dir", str(override)])

        assert result.exit_code == 0, result.output
        assert stop_file_path(override).exists()
        assert not stop_file_path(tmp_path / "unused").exists()

    def test_marker_is_written_even_when_no_job_is_queued(self, tmp_path):
        """A chain with nothing in the queue can still be stopped in advance."""
        working_dir = tmp_path / "run_1"
        working_dir.mkdir()
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: working_dir)

        result, _, cancel_jobs = _invoke(config_path, sim_config, ["-r", "1"], job_ids=())

        assert result.exit_code == 0, result.output
        assert stop_file_path(working_dir).exists()
        assert "no queued or running job" in result.output
        cancel_jobs.assert_called_once_with([])

    def test_stop_only_leaves_the_running_job_alone(self, tmp_path):
        working_dir = tmp_path / "run_1"
        working_dir.mkdir()
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: working_dir)

        result, lookup, cancel_jobs = _invoke(config_path, sim_config, ["-r", "1", "--stop-only"])

        assert result.exit_code == 0, result.output
        assert stop_file_path(working_dir).exists()
        lookup.assert_not_called()
        cancel_jobs.assert_not_called()

    def test_dry_run_changes_nothing(self, tmp_path):
        working_dir = tmp_path / "run_1"
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: working_dir)

        result, _, cancel_jobs = _invoke(config_path, sim_config, ["-r", "1", "--dry-run"])

        assert result.exit_code == 0, result.output
        assert not stop_file_path(working_dir).exists()
        cancel_jobs.assert_not_called()
        assert "[dry-run]" in result.output

    def test_resume_removes_the_marker(self, tmp_path):
        working_dir = tmp_path / "run_1"
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: working_dir)
        write_stop_file(working_dir, str(config_path), 1)

        result, _, cancel_jobs = _invoke(config_path, sim_config, ["-r", "1", "--resume"])

        assert result.exit_code == 0, result.output
        assert not stop_file_path(working_dir).exists()
        cancel_jobs.assert_not_called()
        assert "polyzymd submit" in result.output

    def test_resume_without_a_marker_warns(self, tmp_path):
        working_dir = tmp_path / "run_1"
        working_dir.mkdir()
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: working_dir)

        result, _, _ = _invoke(config_path, sim_config, ["-r", "1", "--resume"])

        assert result.exit_code == 0, result.output
        assert "nothing to resume" in result.output

    def test_cancel_then_resume_round_trips(self, tmp_path):
        working_dir = tmp_path / "run_1"
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: working_dir)

        _invoke(config_path, sim_config, ["-r", "1"])
        assert stop_file_path(working_dir).exists()
        _invoke(config_path, sim_config, ["-r", "1", "--resume"])
        assert not stop_file_path(working_dir).exists()


# ---------------------------------------------------------------------------
# GROMACS chains
# ---------------------------------------------------------------------------


def _run_gromacs_job_script(tmp_path: Path, replicate_dir: Path, monkeypatch) -> list[str]:
    """Run a GROMACS job script for *replicate_dir* with stub tools; return the tools it called."""
    import os
    import subprocess

    from polyzymd.engines.gromacs.slurm import GromacsSlurmScriptGenerator
    from polyzymd.workflow.slurm import SlurmConfig

    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path", lambda: "/tmp/pixi.toml"
    )
    script = tmp_path / "job.sh"
    script.write_text(
        GromacsSlurmScriptGenerator(slurm_config=SlurmConfig()).generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            # `polyzymd submit` runs a GROMACS job in <replicate>/gromacs.
            working_dir=str(replicate_dir / "gromacs"),
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
    )
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    calls = tmp_path / "calls"
    calls.unlink(missing_ok=True)
    for tool in ("sbatch", "module", "gmx", "pixi", "nvidia-smi"):
        stub = bin_dir / tool
        stub.write_text(f'#!/bin/sh\necho {tool} >> "{calls}"\nexit 1\n')
        stub.chmod(0o755)
    subprocess.run(
        ["bash", "--noprofile", "--norc", str(script)],
        # A cluster's Lmod `module` (an exported function, or one BASH_ENV
        # defines) would shadow the stubs.
        env={
            **{
                k: v
                for k, v in os.environ.items()
                if not k.startswith("BASH_FUNC_") and k != "BASH_ENV"
            },
            "PATH": f"{bin_dir}:{os.environ['PATH']}",
        },
        capture_output=True,
        text=True,
        timeout=60,
    )
    return calls.read_text().split() if calls.exists() else []


class TestCancelGromacsChain:
    """`polyzymd cancel` stops and restarts a GROMACS chain, whose job works in <replicate>/gromacs."""

    def _gromacs_config(self, tmp_path: Path, replicate_dir: Path):
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: replicate_dir)
        sim_config.engine = "gromacs"
        return config_path, sim_config

    def test_cancel_stops_a_gromacs_job(self, tmp_path, monkeypatch):
        replicate_dir = tmp_path / "run_1"
        config_path, sim_config = self._gromacs_config(tmp_path, replicate_dir)

        result, _, _ = _invoke(config_path, sim_config, ["-r", "1", "--stop-only"])
        assert result.exit_code == 0, result.output

        calls = _run_gromacs_job_script(tmp_path, replicate_dir, monkeypatch)
        assert "gmx" not in calls and "sbatch" not in calls, calls

    def test_cancel_writes_stop_where_a_running_gromacs_job_looks(self, tmp_path):
        """A chain submitted earlier checks only <replicate>/gromacs/STOP."""
        replicate_dir = tmp_path / "run_1"
        (replicate_dir / "gromacs").mkdir(parents=True)
        config_path, sim_config = self._gromacs_config(tmp_path, replicate_dir)

        result, _, _ = _invoke(config_path, sim_config, ["-r", "1", "--stop-only"])
        assert result.exit_code == 0, result.output
        assert stop_file_path(replicate_dir / "gromacs").exists()

    def test_resume_lets_a_gromacs_job_run_again(self, tmp_path, monkeypatch):
        replicate_dir = tmp_path / "run_1"
        (replicate_dir / "gromacs").mkdir(parents=True)
        stop_file_path(replicate_dir).write_text("PolyzyMD STOP marker\n")
        stop_file_path(replicate_dir / "gromacs").write_text("PolyzyMD STOP marker\n")
        config_path, sim_config = self._gromacs_config(tmp_path, replicate_dir)

        result, _, _ = _invoke(config_path, sim_config, ["-r", "1", "--resume"])
        assert result.exit_code == 0, result.output
        assert not stop_file_path(replicate_dir).exists()
        assert not stop_file_path(replicate_dir / "gromacs").exists()

        calls = _run_gromacs_job_script(tmp_path, replicate_dir, monkeypatch)
        assert "pixi" in calls or "gmx" in calls, calls


# ---------------------------------------------------------------------------
# cancel_slurm_jobs
# ---------------------------------------------------------------------------


class TestCancelSlurmJobs:
    """scancel is best effort, exactly like the squeue duplicate check."""

    def test_empty_input_does_not_shell_out(self):
        from polyzymd.workflow import daisy_chain

        with patch.object(daisy_chain.subprocess, "run") as run:
            assert daisy_chain.cancel_slurm_jobs([]) == []
        run.assert_not_called()

    def test_successful_cancel_reports_the_ids(self):
        from polyzymd.workflow import daisy_chain

        completed = MagicMock(returncode=0, stdout="", stderr="")
        with patch.object(daisy_chain.subprocess, "run", return_value=completed) as run:
            assert daisy_chain.cancel_slurm_jobs(["1", "2"]) == ["1", "2"]
        assert run.call_args[0][0] == ["scancel", "1", "2"]

    def test_missing_scancel_is_not_fatal(self):
        from polyzymd.workflow import daisy_chain

        with patch.object(daisy_chain.subprocess, "run", side_effect=FileNotFoundError):
            assert daisy_chain.cancel_slurm_jobs(["1"]) == []

    def test_nonzero_exit_reports_nothing_cancelled(self):
        from polyzymd.workflow import daisy_chain

        completed = MagicMock(returncode=1, stdout="", stderr="Invalid job id")
        with patch.object(daisy_chain.subprocess, "run", return_value=completed):
            assert daisy_chain.cancel_slurm_jobs(["1"]) == []


class TestCancelNeverRun:
    """A run that never started and has no job leaves nothing to cancel."""

    def test_writes_nothing_when_no_folder_and_no_job(self, tmp_path):
        working_dir = tmp_path / "scratch" / "run_1"
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: working_dir)

        result, _, cancel_jobs = _invoke(config_path, sim_config, ["-r", "1"], job_ids=())

        assert result.exit_code == 0, result.output
        assert not working_dir.exists()
        assert not (tmp_path / "scratch").exists()
        assert "nothing to cancel" in result.output
        cancel_jobs.assert_not_called()
