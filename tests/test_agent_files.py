"""The agent guidance PolyzyMD ships and writes into projects (polyzymd/agent_files.py)."""

from __future__ import annotations

import re
import shlex
import subprocess
import tomllib
from pathlib import Path

import click
import pytest
from click.testing import CliRunner

import polyzymd
from polyzymd import agent_files
from polyzymd.agent_files import (
    DEVELOPER_SKILLS,
    PACKAGED_SKILLS,
    read_marker,
    sync,
    user_skills,
    with_marker,
)
from polyzymd.analyses.exceptions import ProtocolError

PACKAGE = Path(polyzymd.__file__).parent
REPO = PACKAGE.parents[1]
TEMPLATE = PACKAGE / "templates" / "agents" / "AGENTS.md"
in_checkout = pytest.mark.skipif(
    not (REPO / "pyproject.toml").is_file(), reason="needs the PolyzyMD repository"
)


@in_checkout
def test_the_wheel_copies_the_repository_skills_and_nothing_else_holds_them() -> None:
    """One source: .claude/skills/ goes into the wheel, and no second copy is committed."""
    pyproject = tomllib.loads((REPO / "pyproject.toml").read_text())
    wheel = pyproject["tool"]["hatch"]["build"]["targets"]["wheel"]
    assert wheel["force-include"] == {".claude/skills": f"polyzymd/{PACKAGED_SKILLS}"}
    assert "/.claude/skills" in pyproject["tool"]["hatch"]["build"]["targets"]["sdist"]["include"]
    assert not (PACKAGE / PACKAGED_SKILLS).exists()
    assert agent_files.skills_source() == REPO / ".claude" / "skills"


@in_checkout
def test_every_skill_but_the_developer_ones_is_a_user_skill() -> None:
    """A skill added to .claude/skills/ ships without being listed anywhere."""
    folders = {path.parent.name for path in (REPO / ".claude" / "skills").glob("*/SKILL.md")}
    assert DEVELOPER_SKILLS <= folders
    assert set(user_skills()) == folders - DEVELOPER_SKILLS
    assert {"polyzymd-analyze", "polyzymd-simulate"} <= set(user_skills())


def test_the_user_guide_is_self_contained() -> None:
    """The project AGENTS.md names no skill file, MCP server or contributor rule."""
    text = TEMPLATE.read_text()
    for word in ("SKILL.md", ".claude", "MCP", "contributor", "worktree", "Joe"):
        assert word not in text, word
    assert (TEMPLATE.parent / "CLAUDE.md").read_text() == "@AGENTS.md\n"


def _commands() -> list[str]:
    """Every polyzymd command the user AGENTS.md names, in code lines or inline code."""
    text = TEMPLATE.read_text()
    found = re.findall(r"^(polyzymd [^#\n]*?)\s*(?:#.*)?$", text, flags=re.MULTILINE)
    found += re.findall(r"`(polyzymd [^`]+)`", text)
    return found


@pytest.mark.parametrize("command", _commands())
def test_every_command_of_the_user_guide_exists(command: str) -> None:
    """Each named subcommand answers --help, and takes each option (and choice) the guide gives."""
    from polyzymd.cli.main import cli

    words = shlex.split(command)[1:]
    path, current = [], cli
    while words and isinstance(current, click.Group) and words[0] in current.commands:
        path.append(words[0])
        current = current.commands[words.pop(0)]
    assert path and not isinstance(current, click.Group), command
    result = CliRunner().invoke(cli, [*path, "--help"])
    assert result.exit_code == 0, (command, result.output)
    options = {name: param for param in current.params for name in param.opts}
    for index, word in enumerate(words):
        if word.startswith("-"):
            assert word in options, (command, word)
            choice = options[word].type
            given = words[index + 1] if index + 1 < len(words) else "<"
            if isinstance(choice, click.Choice) and not given.startswith("<"):
                assert given in choice.choices, (command, word)


def test_the_marker_follows_the_front_matter_and_hashes_the_rest() -> None:
    """Claude Code needs the front matter first, so the marker comes right after it."""
    body = "---\nname: x\ndescription: y\n---\n\n# x\n"
    text = with_marker(body, "9.9")
    assert text.startswith("---\nname: x\ndescription: y\n---\n<!-- polyzymd-managed: version=9.9 ")
    version, sha, rest = read_marker(text)
    assert (version, rest) == ("9.9", body)
    assert read_marker(with_marker("@AGENTS.md\n", "9.9"))[2] == "@AGENTS.md\n"
    assert read_marker("no marker\n") is None


def test_sync_creates_keeps_skips_and_forces(tmp_path: Path) -> None:
    """Missing files are created, edited ones skipped unless forced, and nothing else changes."""
    (tmp_path / "study.yaml").write_text("conditions: {}\n")
    assert {action.action for action in sync(tmp_path, dry_run=True)} == {"created"}
    assert not (tmp_path / "AGENTS.md").exists()
    created = sync(tmp_path)
    paths = {action.path for action in created}
    assert {"AGENTS.md", "CLAUDE.md"} <= paths
    for name in user_skills():
        assert {f".claude/skills/{name}/SKILL.md", f".agents/skills/{name}/SKILL.md"} <= paths
    assert {action.action for action in sync(tmp_path)} == {"unchanged"}

    agents = tmp_path / "AGENTS.md"
    agents.write_text(agents.read_text() + "\nOur own rule.\n")
    (tmp_path / "CLAUDE.md").write_text("@AGENTS.md\n")  # no marker: the user's own file
    lines = {action.path: action.line() for action in sync(tmp_path)}
    assert lines["AGENTS.md"] == "skipped: AGENTS.md (edited; --force to overwrite)"
    assert lines["CLAUDE.md"] == "skipped: CLAUDE.md (edited; --force to overwrite)"
    assert "Our own rule." in agents.read_text()
    forced = {action.path: action.action for action in sync(tmp_path, force=True)}
    assert forced["AGENTS.md"] == forced["CLAUDE.md"] == "updated"
    assert "Our own rule." not in agents.read_text()


def test_sync_updates_an_unedited_file_to_a_new_version(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A file still as an older PolyzyMD wrote it takes the installed version's text."""
    sync(tmp_path)
    shipped = agent_files.managed_files()
    shipped["AGENTS.md"] += "\nA new rule.\n"
    monkeypatch.setattr(agent_files, "managed_files", lambda: shipped)
    actions = {action.path: action.action for action in sync(tmp_path)}
    assert actions["AGENTS.md"] == "updated" and actions["CLAUDE.md"] == "unchanged"
    assert "A new rule." in (tmp_path / "AGENTS.md").read_text()


@pytest.mark.usefixtures("git_identity")
def test_project_init_commits_the_agent_files_once_at_its_root(tmp_path: Path) -> None:
    """project init writes the guidance at the project root, in the first commit."""
    from polyzymd.analyses.project_scaffold import create_project

    created = create_project(tmp_path / "paper", ["lipa"])
    root = created.root
    assert read_marker((root / "AGENTS.md").read_text())[0] == polyzymd.__version__
    assert not (root / "lipa" / "AGENTS.md").exists()
    if created.commit is not None:
        tracked = subprocess.run(
            ["git", "-C", str(root), "ls-files"], capture_output=True, text=True
        ).stdout.split()
        assert {"AGENTS.md", "CLAUDE.md"} <= set(tracked)
        assert f".agents/skills/{next(iter(user_skills()))}/SKILL.md" in tracked


def test_agent_files_are_no_analysis_input() -> None:
    """A sync leaves uncommitted guides, which must not make analyze warn or freeze refuse."""
    from polyzymd.analyses.study_git import is_output

    for path in ("AGENTS.md", "lipa/CLAUDE.md", ".claude/skills/a/SKILL.md", ".agents/x"):
        assert is_output(path), path
    assert not is_output("analyses/AGENTS.py")
