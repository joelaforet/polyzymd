"""Tests for the ``polyzymd cancel`` command and the STOP marker.

A self-resubmitting chain cannot be stopped with ``scancel`` alone: SIGTERM
makes ``run-segment`` exit 99, which the job wrapper reads as "interrupted,
work remains", so it queues a successor within seconds.  ``polyzymd cancel``
writes a STOP marker that the wrapper honours, then cancels the jobs.

Covers:
- STOP marker contents and location
- job lookup and cancellation by job name
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
        lookup.assert_called_once_with("pzmd_r1")
        cancel_jobs.assert_called_once_with(["4242"])
        assert "4242" in result.output

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
        config_path, sim_config = _config_and_mock(tmp_path, lambda r: working_dir)

        result, _, cancel_jobs = _invoke(config_path, sim_config, ["-r", "1"], job_ids=())

        assert result.exit_code == 0, result.output
        assert stop_file_path(working_dir).exists()
        assert "no queued or running job" in result.output
        cancel_jobs.assert_called_once_with([])

    def test_stop_only_leaves_the_running_job_alone(self, tmp_path):
        working_dir = tmp_path / "run_1"
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
