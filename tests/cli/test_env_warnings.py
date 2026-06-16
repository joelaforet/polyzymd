"""Tests for CLI pixi environment advisory warnings."""

from __future__ import annotations

from click.testing import CliRunner

from polyzymd.cli.compare import compare
from polyzymd.cli.env_warnings import (
    env_warnings_disabled,
    get_active_pixi_environment,
    warn_if_wrong_pixi_env,
)
from polyzymd.cli.main import cli


def test_get_active_pixi_environment_returns_none_when_unset(monkeypatch) -> None:
    """Missing pixi environment metadata should return ``None``."""
    monkeypatch.delenv("PIXI_ENVIRONMENT_NAME", raising=False)

    assert get_active_pixi_environment() is None


def test_get_active_pixi_environment_strips_empty_value(monkeypatch) -> None:
    """Blank pixi environment metadata should behave like an inactive shell."""
    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "   ")

    assert get_active_pixi_environment() is None


def test_env_warnings_disabled_accepts_truthy_values(monkeypatch) -> None:
    """The suppression environment variable should accept simple truthy values."""
    for value in ("1", "true", "TRUE", "yes", "on"):
        monkeypatch.setenv("POLYZYMD_DISABLE_ENV_WARNINGS", value)
        assert env_warnings_disabled()


def test_warn_if_wrong_pixi_env_writes_stderr_only(monkeypatch, capsys) -> None:
    """Wrong active environments should emit an advisory warning on stderr."""

    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "analysis")
    monkeypatch.delenv("POLYZYMD_DISABLE_ENV_WARNINGS", raising=False)

    warn_if_wrong_pixi_env("demo", "build")
    captured = capsys.readouterr()

    assert captured.out == ""
    assert "Warning: polyzymd demo" in captured.err
    assert "recommended environment is 'build'" in captured.err


def test_warn_if_wrong_pixi_env_skips_accepted_environment(monkeypatch, capsys) -> None:
    """Accepted active environments should not warn."""

    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "analysis")

    warn_if_wrong_pixi_env("demo", "build", accepted=("build", "analysis"))
    captured = capsys.readouterr()

    assert captured.out == ""
    assert captured.err == ""


def test_warn_if_wrong_pixi_env_respects_suppressor(monkeypatch, capsys) -> None:
    """The suppressor should silence advisory warnings."""

    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "analysis")
    monkeypatch.setenv("POLYZYMD_DISABLE_ENV_WARNINGS", "yes")

    warn_if_wrong_pixi_env("demo", "build")
    captured = capsys.readouterr()

    assert captured.out == ""
    assert captured.err == ""


def test_main_validate_warns_for_analysis_environment(monkeypatch, tmp_path) -> None:
    """System validation should recommend the build environment."""
    config_path = tmp_path / "config.yaml"
    config_path.write_text("name: demo\n", encoding="utf-8")
    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "analysis")
    monkeypatch.delenv("POLYZYMD_DISABLE_ENV_WARNINGS", raising=False)
    runner = CliRunner()

    result = runner.invoke(cli, ["validate", "-c", str(config_path)])

    assert result.exit_code != 0
    assert "Warning: polyzymd validate" in result.output
    assert "recommended environment is 'build'" in result.output


def test_run_openmm_warns_for_build_environment(monkeypatch, tmp_path) -> None:
    """Local OpenMM execution should recommend CUDA simulation environments."""
    config_path = tmp_path / "config.yaml"
    config_path.write_text("name: demo\n", encoding="utf-8")
    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "build")
    monkeypatch.delenv("POLYZYMD_DISABLE_ENV_WARNINGS", raising=False)
    runner = CliRunner()

    result = runner.invoke(cli, ["run", "-c", str(config_path), "--engine", "openmm"])

    assert result.exit_code != 0
    assert "Warning: polyzymd run --engine openmm" in result.output
    assert "sim-cuda-12-4" in result.output
    assert "sim-cuda-12-6" in result.output


def test_compare_init_warns_for_analysis_environment(monkeypatch, tmp_path) -> None:
    """Comparison scaffolding should recommend the build environment."""
    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "analysis")
    monkeypatch.delenv("POLYZYMD_DISABLE_ENV_WARNINGS", raising=False)
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path):
        result = runner.invoke(compare, ["init", "-n", "study"])

    assert result.exit_code == 0
    assert "Created comparison project" in result.output
    assert "Warning: polyzymd compare init" in result.output
    assert "recommended environment is 'build'" in result.output


def test_compare_status_json_skips_environment_warning(monkeypatch, tmp_path) -> None:
    """Machine-readable analysis status output should not get advisory warnings."""
    config_path = tmp_path / "comparison.yaml"
    config_path.write_text("name: demo\nconditions: []\nplugins: {}\n", encoding="utf-8")
    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "build")
    runner = CliRunner()

    result = runner.invoke(compare, ["status", "rmsf", "-f", str(config_path), "--json"])

    assert "Warning: polyzymd compare status" not in result.output
