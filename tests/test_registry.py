"""The per-user list of projects: init records, register, unregister and list."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from polyzymd.cli.main import cli


@pytest.fixture()
def config(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """An empty settings folder of this test's own."""
    folder = tmp_path / "config"
    monkeypatch.setenv("POLYZYMD_CONFIG_DIR", str(folder))
    return folder


def _run(*args: str) -> str:
    result = CliRunner().invoke(cli, list(args))
    assert result.exit_code == 0, result.output
    return result.output


def test_the_settings_folder_follows_xdg(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    from polyzymd.registry import registry_path

    monkeypatch.delenv("POLYZYMD_CONFIG_DIR")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    assert registry_path() == tmp_path / "polyzymd" / "projects.yaml"
    monkeypatch.setenv("XDG_CONFIG_HOME", "relative")
    assert registry_path() == Path.home() / ".config" / "polyzymd" / "projects.yaml"


def test_init_records_projects_and_lone_studies(config: Path, tmp_path: Path) -> None:
    _run("project", "init", str(tmp_path / "Paper_1"), "--study", "a", "--study", "b", "--no-git")
    _run("study", "init", str(tmp_path / "lone"), "--new-condition", "Water", "--no-git")
    _run("study", "init", str(tmp_path / "Paper_1" / "c"), "--new-condition", "W", "--no-git")
    listed = json.loads(_run("project", "list", "--format", "json"))
    assert [(e["name"], e["kind"], e["studies"], e["exists"]) for e in listed] == [
        ("Paper_1", "project", ["a", "b"], True),
        ("lone", "study", [], True),
    ]
    assert listed[0]["path"] == str((tmp_path / "Paper_1").resolve())
    assert listed[0]["created"]


def test_a_broken_list_does_not_fail_init(config: Path, tmp_path: Path) -> None:
    config.mkdir()
    (config / "projects.yaml").write_text("projects: [unclosed\n")
    output = _run("project", "init", str(tmp_path / "P"), "--study", "a", "--no-git")
    assert "warning: could not add" in output
    assert (tmp_path / "P" / "project.yaml").is_file()


def test_register_names_collisions_refreshes_and_unregisters(config: Path, tmp_path: Path) -> None:
    for parent in ("one", "two"):
        study = tmp_path / parent / "Paper"
        study.mkdir(parents=True)
        (study / "study.yaml").write_text("conditions: {}\n")
    (tmp_path / "one" / "Paper" / "sub").mkdir()
    assert "registered study Paper:" in _run("project", "register", str(tmp_path / "one/Paper/sub"))
    assert "registered study Paper-2:" in _run("project", "register", str(tmp_path / "two/Paper"))
    assert "refreshed study Paper:" in _run("project", "register", str(tmp_path / "one/Paper"))
    assert "unregistered Paper-2" in _run("project", "unregister", "Paper-2")
    result = CliRunner().invoke(cli, ["project", "unregister", "Paper-2"])
    assert result.exit_code == 2 and "polyzymd project list" in result.output
    assert (tmp_path / "two" / "Paper").is_dir()


def test_list_shows_a_missing_folder_and_keeps_it(config: Path, tmp_path: Path) -> None:
    from polyzymd.registry import register

    (tmp_path / "gone").mkdir()
    (tmp_path / "gone" / "project.yaml").write_text("studies: {x: x}\n")
    register(tmp_path / "gone")
    (tmp_path / "gone" / "project.yaml").unlink()
    (tmp_path / "gone").rmdir()
    line = _run("project", "list", "--format", "agent").strip()
    assert line == f"gone kind=project exists=no studies=x path={(tmp_path / 'gone').resolve()}"
    table = _run("project", "list").splitlines()
    assert table[0].split() == ["NAME", "KIND", "EXISTS", "STUDIES", "PATH"]
    assert table[1].split()[:4] == ["gone", "project", "no", "x"]


def test_an_empty_list_says_how_to_add_one(config: Path) -> None:
    assert "polyzymd project register PATH" in _run("project", "list", "--format", "agent")
    assert json.loads(_run("project", "list", "--format", "json")) == []
