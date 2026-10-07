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


def test_add_study_writes_the_study_and_lists_it(tmp_path: Path) -> None:
    """add-study makes the study folder project init would, and lists it after the others."""
    from polyzymd.analyses.project_file import load_project_file
    from polyzymd.cli.main import cli

    project = tmp_path / "paper"
    init = ["project", "init", str(project), "--study", "lipa363", "--no-git"]
    assert CliRunner().invoke(cli, init).exit_code == 0
    result = CliRunner().invoke(cli, ["project", "add-study", "calb343", "--project", str(project)])

    assert result.exit_code == 0, result.output
    assert (project / "calb343" / "study.yaml").read_text() == (
        project / "lipa363" / "study.yaml"
    ).read_text()
    assert not (project / "calb343" / "LICENSE-code").exists()
    assert list(load_project_file(project).studies) == ["lipa363", "calb343"]
    assert "# label: folder holding its study.yaml" in (project / "project.yaml").read_text()


def test_add_study_refuses_a_taken_label_and_one_that_is_no_folder_name(tmp_path: Path) -> None:
    """add-study exits 2, and writes nothing, for a listed label or one with capitals or spaces."""
    from polyzymd.cli.main import cli

    project = tmp_path / "paper"
    init = ["project", "init", str(project), "--study", "lipa363", "--no-git"]
    assert CliRunner().invoke(cli, init).exit_code == 0
    before = (project / "project.yaml").read_text()
    for label, message in (("lipa363", "already has a study"), ("CalB 343", "not a folder name")):
        result = CliRunner().invoke(cli, ["project", "add-study", label, "--project", str(project)])
        assert result.exit_code == 2 and message in result.output, result.output
    assert (project / "project.yaml").read_text() == before
    assert not (project / "CalB 343").exists()
