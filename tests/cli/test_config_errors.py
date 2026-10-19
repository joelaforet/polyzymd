"""Tests for reporting a config that cannot be loaded."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from polyzymd.cli.config_errors import config_load_lines
from polyzymd.cli.main import cli
from tests._support.analysis_testkit import write_simulation_config

UNKNOWN_PDB = "error: pdb: unknown key; did you mean enzyme.pdb_path?"


def _misplaced_pdb_config(tmp_path: Path) -> Path:
    """Write a config whose PDB path is a top-level ``pdb:`` instead of ``enzyme.pdb_path``."""
    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    data = yaml.safe_load(path.read_text())
    data["pdb"] = data["enzyme"].pop("pdb_path")
    path.write_text(yaml.safe_dump(data))
    return path


@pytest.mark.parametrize(
    ("args", "exit_code"),
    [
        (["validate"], 1),
        (["build"], 1),
        (["run"], 1),
        (["submit", "--dry-run"], 1),
        (["run-segment", "-r", "1"], 1),
        (["check-progress", "-r", "1"], 3),
        (["cancel"], 1),
        (["status"], 1),
        (["recover", "-r", "1"], 1),
        (["hash-trajectories"], 2),
    ],
    ids=lambda value: value[0] if isinstance(value, list) else str(value),
)
def test_every_command_reports_a_bad_config_the_same_way(
    tmp_path: Path, args: list[str], exit_code: int
) -> None:
    config = _misplaced_pdb_config(tmp_path)

    result = CliRunner().invoke(cli, [args[0], "-c", str(config), *args[1:]])

    assert result.exit_code == exit_code, result.output
    assert isinstance(result.exception, SystemExit), result.exception
    assert UNKNOWN_PDB in result.stderr
    assert f"fix: correct these keys in {config} and run polyzymd validate -c {config} again." in (
        result.stderr
    )
    assert "Failed to load config" not in result.output


def test_study_check_reports_a_bad_condition_config_with_its_label(tmp_path: Path) -> None:
    from polyzymd.analyses.study_scaffold import create_study

    created = create_study(tmp_path / "s", new_conditions=["A", "B"], git=False)
    config = _misplaced_pdb_config(tmp_path)
    created.conditions["B"].write_text(config.read_text())

    result = CliRunner().invoke(cli, ["study", "check", str(created.root)])

    assert result.exit_code == 2, result.output
    path = created.conditions["B"]
    assert (
        f"error: condition B: cannot read {path}: pdb: unknown key; "
        "did you mean enzyme.pdb_path?"
    ) in result.output
    assert f"fix: correct these keys in {path} and run polyzymd validate -c {path} again." in (
        result.output
    )


def test_a_file_error_keeps_its_fix_and_says_to_validate_again(tmp_path: Path) -> None:
    from polyzymd.config.loader import ConfigFileError

    error = ConfigFileError(tmp_path / "c.yaml", "not valid YAML", "Correct the YAML.")

    assert config_load_lines(error, "c.yaml") == [
        ("error", f"{tmp_path / 'c.yaml'}: not valid YAML"),
        ("fix", "Correct the YAML."),
        ("fix", "correct c.yaml and run polyzymd validate -c c.yaml again."),
    ]


def test_a_value_error_with_its_own_fix_keeps_it() -> None:
    error = ValueError("production: ensemble NVE is not run.\nfix: use ensemble NVT or NPT.")

    assert config_load_lines(error, "c.yaml") == [
        ("error", "production: ensemble NVE is not run."),
        ("fix", "use ensemble NVT or NPT."),
        ("fix", "correct c.yaml and run polyzymd validate -c c.yaml again."),
    ]
