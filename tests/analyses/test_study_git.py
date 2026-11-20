"""study_git: which files of a study count as outputs rather than inputs."""

from __future__ import annotations

import subprocess

from polyzymd.analyses.study_git import is_output, unignored_paths


def test_compiled_python_and_job_files_are_never_inputs() -> None:
    """Freeze refuses uncommitted inputs, so bytecode and job logs must not count as inputs."""
    for path in (
        "lipa/analyses/__pycache__/helper.cpython-311.pyc",
        "x.pyc",
        "results/a/slurm/j.sh",
    ):
        assert is_output(path), path
    for path in ("lipa/analyses/helper.py", "study.yaml", "stats/plan.py"):
        assert not is_output(path), path


def test_describe_says_no_commit_yet_in_full() -> None:
    """A folder without a commit is described as 'no commit yet', not cut short."""
    from polyzymd.analyses.study_git import describe

    text = describe({"commit": None, "inputs_uncommitted": [], "outputs_uncommitted": []})
    assert "no commit yet" in text


def test_runs_are_never_inputs() -> None:
    """A trajectory under runs/, of a study or of a project, is not an uncommitted input."""
    for path in ("runs/water/w_run1/traj.dcd", "runs/lipa/water/w_run1/traj.dcd"):
        assert is_output(path), path


def test_sandbox_files_are_never_inputs() -> None:
    assert is_output("sandbox/check_build.py") and is_output("lipa/sandbox/plot.ipynb")
    assert not is_output("analyses/lid.py")


def test_unignored_paths_names_the_rules_an_old_gitignore_lacks(tmp_path) -> None:
    assert unignored_paths(tmp_path) == []
    subprocess.run(["git", "-C", str(tmp_path), "init", "-q"], check=True)
    (tmp_path / ".gitignore").write_text("runs/\n")
    assert unignored_paths(tmp_path) == ["results/RUN/polyzymd_results/", "slurm_logs/"]
    (tmp_path / ".gitignore").write_text("results/*/polyzymd_results/\nslurm_logs/\n")
    assert unignored_paths(tmp_path) == []
