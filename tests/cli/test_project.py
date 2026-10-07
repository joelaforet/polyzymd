"""polyzymd project: the options of its subcommands."""

from __future__ import annotations

from click.testing import CliRunner


def test_project_check_takes_production() -> None:
    """project check takes --production, as study check does."""
    from polyzymd.cli.main import cli

    assert "--production" in CliRunner().invoke(cli, ["project", "check", "--help"]).output
