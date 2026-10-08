"""Tests for slurm_submit helper."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from polyzymd.workflow.slurm_submit import run_sbatch


@pytest.fixture(autouse=True)
def _sbatch_on_path(monkeypatch):
    monkeypatch.setattr(
        "polyzymd.workflow.slurm_submit.shutil.which", lambda name: "/usr/bin/sbatch"
    )


class TestRunSbatch:
    """Tests for run_sbatch()."""

    @patch("polyzymd.workflow.slurm_submit.subprocess.run")
    def test_missing_sbatch_is_an_error(self, mock_run, tmp_path, monkeypatch):
        """Without sbatch on PATH, submission fails and names the module to load."""
        monkeypatch.setattr("polyzymd.workflow.slurm_submit.shutil.which", lambda name: None)
        script = tmp_path / "job.sh"
        script.touch()

        with pytest.raises(RuntimeError, match=r"sbatch not found.*ml slurm/blanca.*job\.sh"):
            run_sbatch(script)
        mock_run.assert_not_called()

    @patch("polyzymd.workflow.slurm_submit.subprocess.run")
    def test_job_starts_with_a_clean_environment(self, mock_run, tmp_path):
        """sbatch runs directly with --export=NONE; no login-node shell or module load."""
        script = tmp_path / "my job.sh"
        script.touch()
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="Submitted batch job 123\n",
            stderr="",
        )

        result = run_sbatch(script)

        mock_run.assert_called_once()
        assert mock_run.call_args[0][0] == ["sbatch", "--export=NONE", str(script)]
        assert result.returncode == 0

    @patch("polyzymd.workflow.slurm_submit.subprocess.run")
    def test_the_log_folder_is_made_before_sbatch(self, mock_run, tmp_path):
        """sbatch cannot write a log into a missing folder, so it is made first."""
        logs = tmp_path / "slurm_logs"
        script = tmp_path / "job.sh"
        script.write_text(f"#!/bin/bash\n#SBATCH --output={logs}/job.%j.out\n")
        mock_run.side_effect = lambda *a, **kw: MagicMock(returncode=0, stdout="", stderr="")

        run_sbatch(script)

        assert logs.is_dir()
