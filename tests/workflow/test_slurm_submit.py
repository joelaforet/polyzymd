"""Tests for slurm_submit helper."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

from polyzymd.workflow.slurm_submit import run_sbatch


class TestRunSbatch:
    """Tests for run_sbatch()."""

    @patch("polyzymd.workflow.slurm_submit.subprocess.run")
    def test_no_module_load_calls_sbatch_directly(self, mock_run, tmp_path):
        """Without module_load, sbatch is called as a plain subprocess."""
        script = tmp_path / "job.sh"
        script.touch()
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="Submitted batch job 123\n",
            stderr="",
        )

        result = run_sbatch(script)

        mock_run.assert_called_once()
        args = mock_run.call_args
        assert args[0][0] == ["sbatch", str(script)]
        assert result.returncode == 0

    @patch("polyzymd.workflow.slurm_submit.subprocess.run")
    def test_with_module_load_wraps_in_bash(self, mock_run, tmp_path):
        """With module_load, sbatch is wrapped in bash -lc."""
        script = tmp_path / "job.sh"
        script.touch()
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="Submitted batch job 456\n",
            stderr="",
        )

        result = run_sbatch(script, module_load="module load slurm/blanca gromacs/2024.2")

        mock_run.assert_called_once()
        args = mock_run.call_args
        cmd = args[0][0]
        assert cmd[0] == "bash"
        assert cmd[1] == "-lc"
        assert "module load slurm/blanca gromacs/2024.2" in cmd[2]
        assert "sbatch" in cmd[2]
        assert str(script) in cmd[2]
        assert result.returncode == 0

    @patch("polyzymd.workflow.slurm_submit.subprocess.run")
    def test_extra_args_passed_to_sbatch(self, mock_run, tmp_path):
        """Extra arguments are placed between sbatch and script path."""
        script = tmp_path / "job.sh"
        script.touch()
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="Submitted batch job 789\n",
            stderr="",
        )

        run_sbatch(script, extra_args=("--export=NONE",))

        args = mock_run.call_args
        assert args[0][0] == ["sbatch", "--export=NONE", str(script)]

    @patch("polyzymd.workflow.slurm_submit.subprocess.run")
    def test_module_load_with_extra_args(self, mock_run, tmp_path):
        """module_load and extra_args both appear in the bash command."""
        script = tmp_path / "job.sh"
        script.touch()
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="Submitted batch job 101\n",
            stderr="",
        )

        run_sbatch(script, module_load="module load slurm/blanca", extra_args=("--export=NONE",))

        args = mock_run.call_args
        cmd = args[0][0]
        assert cmd[0] == "bash"
        shell_cmd = cmd[2]
        assert "module load slurm/blanca" in shell_cmd
        assert "--export=NONE" in shell_cmd

    @patch("polyzymd.workflow.slurm_submit.subprocess.run")
    def test_script_path_with_spaces_is_quoted(self, mock_run, tmp_path):
        """Script paths with spaces are properly shell-quoted."""
        script = tmp_path / "my job.sh"
        script.touch()
        mock_run.return_value = MagicMock(
            returncode=0,
            stdout="Submitted batch job 202\n",
            stderr="",
        )

        run_sbatch(script, module_load="module load slurm/blanca")

        args = mock_run.call_args
        shell_cmd = args[0][0][2]
        assert "my job.sh" in shell_cmd
