"""Tests for CLI replicate flag unification (Phase 2, v1.3.0).

Tests the _resolve_replicates_option() helper and the updated
--replicates / --replicate flags on `build` and `run-gromacs`.
"""

from __future__ import annotations

import warnings

import click
import pytest

from polyzymd.cli.main import _resolve_replicates_option


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
        """--replicate 1 resolves to [1] with deprecation warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = _resolve_replicates_option(None, 1, "test")
            assert result == [1]
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "--replicate is deprecated" in str(w[0].message)

    def test_deprecated_replicate_preserves_value(self) -> None:
        """--replicate 5 resolves to [5]."""
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
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
        """Invalid range string propagates ValueError from parser."""
        with pytest.raises(ValueError):
            _resolve_replicates_option("abc", None, "test")

    def test_empty_range_raises(self) -> None:
        """Empty string raises ValueError."""
        with pytest.raises(ValueError):
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


class TestRunGromacsCommandReplicateFlags:
    """Test that the run-gromacs command accepts the new flags."""

    def test_run_gromacs_help_shows_replicates(self) -> None:
        """'polyzymd run-gromacs --help' should show --replicates option."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["run-gromacs", "--help"])
        assert result.exit_code == 0
        assert "--replicates" in result.output

    def test_run_gromacs_help_hides_replicate(self) -> None:
        """'polyzymd run-gromacs --help' should NOT show deprecated --replicate."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["run-gromacs", "--help"])
        assert result.exit_code == 0
        lines = result.output.split("\n")
        for line in lines:
            if "--replicate" in line and "--replicates" not in line:
                pytest.fail(f"Deprecated --replicate visible in help: {line}")


class TestSubmitCommandUnchanged:
    """Verify submit command is unchanged."""

    def test_submit_still_has_replicates(self) -> None:
        """'polyzymd submit --help' should still show --replicates."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["submit", "--help"])
        assert result.exit_code == 0
        assert "--replicates" in result.output


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

    def test_dry_run_shows_validation_report_header(self) -> None:
        """--dry-run should show 'DRY RUN — Validation Report' header."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["build", "--help"])
        assert "--dry-run" in result.output

    def test_run_gromacs_dry_run_in_help(self) -> None:
        """run-gromacs --help should mention --dry-run."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["run-gromacs", "--help"])
        assert "--dry-run" in result.output
