"""The hidden retired commands exit 2 and name their replacements."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from polyzymd.cli.main import cli

SKILL = ".claude/skills/polyzymd-analyze/SKILL.md"


@pytest.mark.parametrize(
    "arguments",
    [
        ["compare"],
        ["compare", "run", "rmsf", "-f", "comparison.yaml", "--recompute"],
        ["compare", "submit-all", "--preset", "blanca-shirts", "--dry-run"],
        ["compare", "--help"],
    ],
)
def test_compare_names_the_analyze_commands_docs_and_skill(arguments: list[str]) -> None:
    """Any compare invocation exits 2 with the run and SLURM commands, the docs page and the skill."""
    result = CliRunner().invoke(cli, arguments)

    assert result.exit_code == 2
    assert result.stdout == ""
    assert "polyzymd compare is retired" in result.stderr
    assert "polyzymd analyze NAME -c config.yaml" in result.stderr
    assert "--submit --preset <cluster>" in result.stderr
    assert (
        "https://polyzymd.readthedocs.io/en/latest/how_to/analysis_agent_protocol.html"
        in result.stderr
    )
    assert SKILL in result.stderr


@pytest.mark.parametrize(
    "arguments",
    [["new-analysis"], ["new-analysis", "solvent_shell", "--advanced", "--force"]],
)
def test_new_analysis_names_the_study_api_docs_and_skill(arguments: list[str]) -> None:
    """Any new-analysis invocation exits 2 pointing at the study API page and the skill."""
    result = CliRunner().invoke(cli, arguments)

    assert result.exit_code == 2
    assert "polyzymd new-analysis is retired" in result.stderr
    assert "Study.per_replicate" in result.stderr and "Study.timeseries" in result.stderr
    assert (
        "https://polyzymd.readthedocs.io/en/latest/how_to/study_api.html" in result.stderr
    )
    assert "docs/source/how_to/study_api.md" in result.stderr
    assert SKILL in result.stderr


def test_retired_commands_are_hidden_from_help() -> None:
    """polyzymd --help lists no retired command."""
    result = CliRunner().invoke(cli, ["--help"])

    assert result.exit_code == 0
    commands = [line.split()[0] for line in result.output.splitlines() if line.startswith("  ")]
    assert "compare" not in commands
    assert "new-analysis" not in commands
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
