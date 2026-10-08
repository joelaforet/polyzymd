"""polyzymd hash-trajectories: errors in what it is given."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from polyzymd.cli.main import cli


def test_a_study_that_does_not_exist_is_an_error_with_a_fix(tmp_path: Path) -> None:
    result = CliRunner().invoke(
        cli, ["hash-trajectories", "--study", str(tmp_path / "nope"), "--dry-run"]
    )
    assert result.exit_code == 2, result.output
    assert isinstance(result.exception, SystemExit)
    assert "error:" in result.output and "fix:" in result.output
