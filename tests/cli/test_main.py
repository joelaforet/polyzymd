"""Tests for CLI replicate flag handling and run/submit UX.

Tests the _resolve_replicates_option() helper and the updated
--replicates / --replicate flags on `build` and `run`.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import click
import pytest
from click.testing import CliRunner

from polyzymd.cli.main import _resolve_replicates_option, cli


def _make_dry_run_config() -> SimpleNamespace:
    """Create a minimal config-like object for build --dry-run tests."""

    # Use SimpleNamespace to avoid importing full SimulationConfig internals in unit tests
    # This keeps the test lightweight and avoids heavy dependency requirements
    return SimpleNamespace(
        name="test_sim",
        description="test description",
        enzyme=SimpleNamespace(name="TestEnzyme", pdb_path=Path("structures/enzyme.pdb")),
        substrate=None,
        polymers=None,
        solvent=SimpleNamespace(
            primary=SimpleNamespace(model="tip3p"),
            box=SimpleNamespace(padding=1.2),
            ions=SimpleNamespace(nacl_concentration=0.15),
        ),
        force_field=SimpleNamespace(protein="amber14", small_molecule="openff-2.2.0"),
        thermodynamics=SimpleNamespace(temperature=300.0, pressure=1.0),
        simulation_phases=SimpleNamespace(
            equilibration_stages=[
                SimpleNamespace(name="heating", duration=0.2, ensemble="NVT"),
                SimpleNamespace(name="free_equilibration", duration=0.8, ensemble="NPT"),
            ],
            production=SimpleNamespace(duration=10.0, samples=250),
        ),
        output=SimpleNamespace(
            projects_directory=Path("/tmp/projects"),
            effective_scratch_directory=Path("/tmp/scratch"),
            scratch_directory=Path("/tmp/scratch"),
        ),
        restraints=[],
        get_working_directory=lambda rep: Path(f"/tmp/scratch/run_{rep}"),
    )


class TestResolveReplicatesOption:
    """Unit tests for the _resolve_replicates_option() helper."""

    def test_replicates_single(self) -> None:
        """--replicates '1' resolves to [1]."""
        assert _resolve_replicates_option("1", None, "test") == [1]

    def test_replicates_range(self) -> None:
        """--replicates '1-3' resolves to [1, 2, 3]."""
        assert _resolve_replicates_option("1-3", None, "test") == [1, 2, 3]

    def test_replicates_comma(self) -> None:
        """--replicates '1,3,5' resolves to [1, 3, 5]."""
        assert _resolve_replicates_option("1,3,5", None, "test") == [1, 3, 5]

    def test_replicates_range_with_step(self) -> None:
        """--replicates '1-10:2' resolves to [1, 3, 5, 7, 9]."""
        assert _resolve_replicates_option("1-10:2", None, "test") == [1, 3, 5, 7, 9]

    def test_deprecated_replicate_single(self) -> None:
        """--replicate 1 resolves to [1] with stderr warning."""
        from click.testing import CliRunner

        runner = CliRunner()
        with runner.isolation() as (_stdin, _stdout, stderr):
            result = _resolve_replicates_option(None, 1, "test")
            assert result == [1]
            assert "--replicate is deprecated" in stderr.getvalue().decode("utf-8")

    def test_deprecated_replicate_preserves_value(self) -> None:
        """--replicate 5 resolves to [5]."""
        assert _resolve_replicates_option(None, 5, "test") == [5]

    def test_both_flags_raises(self) -> None:
        """Using both --replicates and --replicate raises UsageError."""
        with pytest.raises(click.UsageError, match="Cannot use both"):
            _resolve_replicates_option("1-3", 1, "test")

    def test_neither_flag_defaults_to_one(self) -> None:
        """Omitting both flags defaults to [1]."""
        assert _resolve_replicates_option(None, None, "test") == [1]

    def test_error_message_includes_command_name(self) -> None:
        """Error message mentions the command name."""
        with pytest.raises(click.UsageError, match="build"):
            _resolve_replicates_option("1-3", 1, "build")

    def test_invalid_range_raises(self) -> None:
        """Invalid range string is converted to Click BadParameter."""
        with pytest.raises(click.BadParameter):
            _resolve_replicates_option("abc", None, "test")

    def test_empty_range_raises(self) -> None:
        """Empty string is converted to Click BadParameter."""
        with pytest.raises(click.BadParameter):
            _resolve_replicates_option("", None, "test")


class TestBuildCommandReplicateFlags:
    """Test that the build command accepts the new flags via Click invocation."""

    def test_build_help_shows_replicates(self) -> None:
        """'polyzymd build --help' should show --replicates option."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["build", "--help"])
        assert result.exit_code == 0
        assert "--replicates" in result.output

    def test_build_help_hides_replicate(self) -> None:
        """'polyzymd build --help' should NOT show deprecated --replicate."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["build", "--help"])
        assert result.exit_code == 0
        # --replicate is hidden, so it should not appear in help output
        # (but --replicates will appear, so we check there's no standalone --replicate)
        lines = result.output.split("\n")
        for line in lines:
            # If a line contains --replicate but NOT --replicates, that's the deprecated flag showing
            if "--replicate" in line and "--replicates" not in line:
                pytest.fail(f"Deprecated --replicate visible in help: {line}")


class TestRunCommandReplicateFlags:
    """Test that the run command accepts replicate flags."""

    def test_run_help_shows_replicates(self) -> None:
        """'polyzymd run --help' should show --replicates option."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--help"])
        assert result.exit_code == 0
        assert "--replicates" in result.output

    def test_run_help_hides_replicate(self) -> None:
        """'polyzymd run --help' should NOT show deprecated --replicate."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--help"])
        assert result.exit_code == 0
        lines = result.output.split("\n")
        for line in lines:
            if "--replicate" in line and "--replicates" not in line:
                pytest.fail(f"Deprecated --replicate visible in help: {line}")

    def test_run_help_shows_engine_flag(self) -> None:
        """'polyzymd run --help' should include required --engine option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--help"])
        assert result.exit_code == 0
        assert "--engine" in result.output

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_run_openmm_rejects_gmx_path(self, mock_from_yaml, tmp_path: Path) -> None:
        """--gmx-path with --engine openmm should raise UsageError."""
        mock_from_yaml.return_value = _make_dry_run_config()
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            [
                "run",
                "-c",
                str(config_path),
                "--engine",
                "openmm",
                "--gmx-path",
                "gmx",
                "--dry-run",
            ],
        )

        assert result.exit_code != 0
        assert "--gmx-path can only be used with --engine gromacs" in result.output


class TestSubmitCommandUnchanged:
    """Verify submit command options."""

    def test_submit_still_has_replicates(self) -> None:
        """'polyzymd submit --help' should still show --replicates."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["submit", "--help"])
        assert result.exit_code == 0
        assert "--replicates" in result.output

    def test_submit_help_shows_generate_only(self) -> None:
        """'polyzymd submit --help' should show --generate-only."""
        runner = CliRunner()
        result = runner.invoke(cli, ["submit", "--help"])
        assert result.exit_code == 0
        assert "--generate-only" in result.output

    def test_submit_dry_run_and_generate_only_conflict(self, tmp_path: Path) -> None:
        """submit should reject --dry-run with --generate-only."""
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            [
                "submit",
                "-c",
                str(config_path),
                "--dry-run",
                "--generate-only",
            ],
        )

        assert result.exit_code != 0
        assert "Cannot use both --dry-run and --generate-only" in result.output


class TestInternalCommandsUnchanged:
    """Verify internal commands still use --replicate (singular, int)."""

    @pytest.mark.parametrize("cmd", ["run-segment", "check-progress", "recover"])
    def test_internal_command_has_replicate(self, cmd: str) -> None:
        """Internal commands should still show --replicate."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, [cmd, "--help"])
        assert result.exit_code == 0
        assert "--replicate" in result.output


class TestDryRunOutput:
    """Tests for the enhanced --dry-run validation report."""

    def test_build_help_shows_dry_run_flag(self) -> None:
        """'polyzymd build --help' should include the --dry-run flag."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["build", "--help"])
        assert "--dry-run" in result.output

    def test_run_help_shows_dry_run_flag(self) -> None:
        """'polyzymd run --help' should include the --dry-run flag."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--help"])
        assert "--dry-run" in result.output


class TestDeprecatedReplicateValidation:
    """Unit tests for deprecated --replicate validation."""

    def test_replicate_zero_raises(self) -> None:
        """Replicate zero should raise BadParameter."""
        with pytest.raises(click.BadParameter):
            _resolve_replicates_option(None, 0, "build")

    def test_replicate_negative_raises(self) -> None:
        """Negative replicate should raise BadParameter."""
        with pytest.raises(click.BadParameter):
            _resolve_replicates_option(None, -1, "build")

    def test_replicate_positive_still_works(self) -> None:
        """Positive replicate should still resolve correctly."""
        assert _resolve_replicates_option(None, 3, "build") == [3]


class TestBuildDryRunEndToEnd:
    """End-to-end CliRunner tests for build dry-run behavior."""

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_build_dry_run_succeeds(self, mock_from_yaml, tmp_path: Path) -> None:
        """Build dry-run should exit successfully and print dry-run header."""
        mock_from_yaml.return_value = _make_dry_run_config()
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(cli, ["build", "-c", str(config_path), "--dry-run"])

        assert result.exit_code == 0
        assert "DRY RUN" in result.output

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_build_dry_run_with_deprecated_replicate_warns(
        self,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """Deprecated --replicate should still work and emit warning."""
        mock_from_yaml.return_value = _make_dry_run_config()
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            ["build", "-c", str(config_path), "--replicate", "1", "--dry-run"],
        )

        assert result.exit_code == 0
        stderr = getattr(result, "stderr", "")
        assert "deprecated" in f"{result.output}\n{stderr}".lower()

    def test_build_with_replicate_zero_fails(self, tmp_path: Path) -> None:
        """Replicate zero should fail with positive-integer guidance."""
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(cli, ["build", "-c", str(config_path), "--replicate", "0"])

        assert result.exit_code != 0
        stderr = getattr(result, "stderr", "")
        message = f"{result.output}\n{stderr}".lower()
        assert "positive" in message or "must be" in message


class TestCliExceptionHandlingNarrowing:
    """Regression tests for narrowed run/submit exception handling."""

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.cli.main._run_openmm_impl")
    def test_run_catches_runtime_error(
        self,
        mock_run_openmm,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """run should catch RuntimeError and exit cleanly."""
        mock_from_yaml.return_value = _make_dry_run_config()
        mock_run_openmm.side_effect = RuntimeError("boom")
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["run", "-c", str(config_path), "--engine", "openmm"],
        )

        assert result.exit_code != 0
        assert "Unexpected error: boom" in result.output

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.cli.main._run_openmm_impl")
    def test_run_does_not_catch_type_error(
        self,
        mock_run_openmm,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """run should not catch TypeError programmer errors."""
        mock_from_yaml.return_value = _make_dry_run_config()
        mock_run_openmm.side_effect = TypeError("programmer bug")
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")

        runner = CliRunner()
        with pytest.raises(TypeError, match="programmer bug"):
            runner.invoke(
                cli,
                ["run", "-c", str(config_path), "--engine", "openmm"],
                catch_exceptions=False,
            )

    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_submit_catches_runtime_error(self, mock_submit_daisy_chain, tmp_path: Path) -> None:
        """submit should catch RuntimeError and exit cleanly."""
        mock_submit_daisy_chain.side_effect = RuntimeError("submission blew up")
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["submit", "-c", str(config_path)],
        )

        assert result.exit_code != 0
        assert "Submission failed: submission blew up" in result.output

    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_submit_does_not_catch_type_error(
        self, mock_submit_daisy_chain, tmp_path: Path
    ) -> None:
        """submit should not catch TypeError programmer errors."""
        mock_submit_daisy_chain.side_effect = TypeError("submit programmer bug")
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")

        runner = CliRunner()
        with pytest.raises(TypeError, match="submit programmer bug"):
            runner.invoke(
                cli,
                ["submit", "-c", str(config_path)],
                catch_exceptions=False,
            )

    def test_build_with_replicate_negative_fails(self, tmp_path: Path) -> None:
        """Negative replicate should fail with positive-integer guidance."""
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(cli, ["build", "-c", str(config_path), "--replicate", "-1"])

        assert result.exit_code != 0
        stderr = getattr(result, "stderr", "")
        message = f"{result.output}\n{stderr}".lower()
        assert "positive" in message or "must be" in message


class TestRunEngineGromacs:
    """Regression tests for ``run --engine gromacs``."""

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_run_gromacs_dry_run_succeeds(self, mock_from_yaml, tmp_path: Path) -> None:
        """run --engine gromacs --dry-run should exit successfully."""
        mock_from_yaml.return_value = _make_dry_run_config()
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            ["run", "-c", str(config_path), "--engine", "gromacs", "--dry-run"],
        )

        assert result.exit_code == 0
        assert "DRY RUN" in result.output
        assert "gromacs" in result.output.lower()

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.cli.main._run_gromacs_impl")
    def test_run_gromacs_delegates_to_impl(
        self,
        mock_run_gromacs,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """run --engine gromacs should delegate to _run_gromacs_impl."""
        mock_from_yaml.return_value = _make_dry_run_config()
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            ["run", "-c", str(config_path), "--engine", "gromacs"],
        )

        assert result.exit_code == 0
        mock_run_gromacs.assert_called_once()

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.cli.main._run_gromacs_impl")
    def test_run_gromacs_passes_gmx_path(
        self,
        mock_run_gromacs,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """--gmx-path should be forwarded to _run_gromacs_impl."""
        mock_from_yaml.return_value = _make_dry_run_config()
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            [
                "run",
                "-c",
                str(config_path),
                "--engine",
                "gromacs",
                "--gmx-path",
                "/usr/local/bin/gmx_mpi",
            ],
        )

        assert result.exit_code == 0
        assert mock_run_gromacs.call_args is not None
        assert mock_run_gromacs.call_args.kwargs["gmx_path"] == "/usr/local/bin/gmx_mpi"

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.cli.main._run_gromacs_impl")
    def test_run_gromacs_catches_runtime_error(
        self,
        mock_run_gromacs,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """run --engine gromacs should catch RuntimeError and exit cleanly."""
        mock_from_yaml.return_value = _make_dry_run_config()
        mock_run_gromacs.side_effect = RuntimeError("GROMACS not found")
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            ["run", "-c", str(config_path), "--engine", "gromacs"],
        )

        assert result.exit_code != 0
        assert "GROMACS not found" in result.output


class TestSubmitDryRunVsGenerateOnly:
    """Side-effect tests for submit --dry-run vs --generate-only."""

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_dry_run_writes_no_files(self, mock_from_yaml, tmp_path: Path) -> None:
        """submit --dry-run should not create any files."""
        mock_from_yaml.return_value = _make_dry_run_config()
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        output_dir = tmp_path / "output"

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "submit",
                "-c",
                str(config_path),
                "--dry-run",
                "--output-dir",
                str(output_dir),
            ],
        )

        assert result.exit_code == 0
        assert "DRY RUN" in result.output
        assert not output_dir.exists()

    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_generate_only_calls_submit_with_flag(self, mock_submit, tmp_path: Path) -> None:
        """submit --generate-only should pass generate_only=True to the backend."""
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["submit", "-c", str(config_path), "--generate-only"],
        )

        assert result.exit_code == 0
        mock_submit.assert_called_once()
        call_kwargs = mock_submit.call_args.kwargs
        assert call_kwargs["generate_only"] is True
        assert call_kwargs["dry_run"] is False
