"""polyzymd project: the options of its subcommands."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner


def test_project_check_takes_production() -> None:
    """project check takes --production, as study check does."""
    from polyzymd.cli.main import cli

    assert "--production" in CliRunner().invoke(cli, ["project", "check", "--help"]).output


def test_project_check_prints_the_git_line_once(tmp_path: Path) -> None:
    """With two studies and 20 uncommitted files, the git line appears once and is short."""
    import subprocess

    from polyzymd.cli.main import cli
    from tests._support.analysis_testkit import write_simulation_config

    for study in ("one", "two"):
        write_simulation_config(tmp_path / study / "control", scratch=tmp_path / "data" / study)
        (tmp_path / study / "study.yaml").write_text(
            "equilibration: 0ns\nconditions: {Control: control}\n"
        )
    (tmp_path / "project.yaml").write_text("studies: {one: one, two: two}\n")
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    for index in range(20):
        (tmp_path / f"note_{index}.txt").write_text("x\n")
    output = CliRunner().invoke(cli, ["project", "check", str(tmp_path)]).output
    git_lines = [line for line in output.splitlines() if line.startswith("git:")]
    assert len(git_lines) == 1, output
    assert "more" in git_lines[0] and "note_19.txt" not in git_lines[0]
