"""The hidden retired commands exit 2 and name their replacement, docs page and agent skill."""

from __future__ import annotations

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


def test_retired_commands_are_hidden_from_help() -> None:
    """polyzymd --help does not list compare."""
    result = CliRunner().invoke(cli, ["--help"])

    assert result.exit_code == 0
    commands = [line.split()[0] for line in result.output.splitlines() if line.startswith("  ")]
    assert "compare" not in commands
