"""The retired ``polyzymd init`` command names its replacement and exits 2."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from polyzymd.cli.main import cli


def test_retired_commands_are_hidden_from_help() -> None:
    """polyzymd --help does not list the retired init command."""
    result = CliRunner().invoke(cli, ["--help"])

    assert result.exit_code == 0
    commands = [line.split()[0] for line in result.output.splitlines() if line.startswith("  ")]
    assert "init" not in commands


@pytest.mark.parametrize("arguments", [["init"], ["init", "-n", "my_simulation"]])
def test_init_points_to_project_init_and_add_condition(arguments: list[str]) -> None:
    """polyzymd init exits 2, writing nothing, with one line naming its replacements."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli, arguments)
        assert not list(Path().iterdir())

    assert result.exit_code == 2
    assert len(result.stderr.splitlines()) == 1
    assert "polyzymd project init" in result.stderr
    assert "polyzymd study add-condition LABEL --new" in result.stderr
