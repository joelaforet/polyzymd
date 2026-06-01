"""Tests for CLI replicate flag handling and run/submit UX.

Tests the _resolve_replicates_option() helper and the --replicates flags on
`build` and `run`.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import click
import pytest
import yaml
from click.testing import CliRunner
from jinja2 import UndefinedError

from polyzymd.cli.main import _resolve_replicates_option, cli
from polyzymd.utils.templates import render_package_template

BRANDING_LINE = "PolyzyMD: Created by Joseph R. Laforet Jr."


def _assert_no_jinja_markers(content: str) -> None:
    """Assert rendered content contains no unresolved Jinja syntax."""
    assert "{{" not in content
    assert "{%" not in content
    assert "%}" not in content


def _make_dry_run_config() -> SimpleNamespace:
    """Create a minimal config-like object for build --dry-run tests."""

    # Use SimpleNamespace to avoid importing full SimulationConfig internals in unit tests
    # This keeps the test lightweight and avoids heavy dependency requirements
    return SimpleNamespace(
        name="test_sim",
        description="test description",
        engine="openmm",
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
            get_job_scripts_directory=lambda: Path("/tmp/projects/job_scripts"),
        ),
        restraints=[],
        get_working_directory=lambda rep: Path(f"/tmp/scratch/run_{rep}"),
    )


class TestResolveReplicatesOption:
    """Unit tests for the _resolve_replicates_option() helper."""

    def test_replicates_single(self) -> None:
        """--replicates '1' resolves to [1]."""
        assert _resolve_replicates_option("1") == [1]

    def test_replicates_range(self) -> None:
        """--replicates '1-3' resolves to [1, 2, 3]."""
        assert _resolve_replicates_option("1-3") == [1, 2, 3]

    def test_replicates_comma(self) -> None:
        """--replicates '1,3,5' resolves to [1, 3, 5]."""
        assert _resolve_replicates_option("1,3,5") == [1, 3, 5]

    def test_replicates_range_with_step(self) -> None:
        """--replicates '1-10:2' resolves to [1, 3, 5, 7, 9]."""
        assert _resolve_replicates_option("1-10:2") == [1, 3, 5, 7, 9]

    def test_neither_flag_defaults_to_one(self) -> None:
        """Omitting both flags defaults to [1]."""
        assert _resolve_replicates_option(None) == [1]

    def test_invalid_range_raises(self) -> None:
        """Invalid range string is converted to Click BadParameter."""
        with pytest.raises(click.BadParameter):
            _resolve_replicates_option("abc")

    def test_empty_range_raises(self) -> None:
        """Empty string is converted to Click BadParameter."""
        with pytest.raises(click.BadParameter):
            _resolve_replicates_option("")


class TestResolveEngineName:
    """Unit tests for _resolve_engine_name()."""

    def test_override_takes_priority(self) -> None:
        """CLI --engine override should take priority over config."""
        from types import SimpleNamespace

        from polyzymd.cli.main import _resolve_engine_name

        config = SimpleNamespace(engine="openmm")
        assert _resolve_engine_name(config, override="gromacs") == "gromacs"

    def test_override_works_without_config_engine(self) -> None:
        """CLI --engine override should work when config has no engine."""
        from types import SimpleNamespace

        from polyzymd.cli.main import _resolve_engine_name

        config = SimpleNamespace()
        assert _resolve_engine_name(config, override="gromacs") == "gromacs"

    def test_reads_config_engine(self) -> None:
        """Should read engine from config when no override."""
        from types import SimpleNamespace

        from polyzymd.cli.main import _resolve_engine_name

        config = SimpleNamespace(engine="gromacs")
        assert _resolve_engine_name(config) == "gromacs"

    def test_missing_engine_raises_usage_error(self) -> None:
        """Missing config engine should raise a usage error."""
        from types import SimpleNamespace

        from polyzymd.cli.main import _resolve_engine_name

        config = SimpleNamespace()
        with pytest.raises(click.UsageError, match="Configure 'engine'"):
            _resolve_engine_name(config)

    def test_none_engine_raises_usage_error(self) -> None:
        """None config engine should raise a usage error."""
        from types import SimpleNamespace

        from polyzymd.cli.main import _resolve_engine_name

        config = SimpleNamespace(engine=None)
        with pytest.raises(click.UsageError, match="Configure 'engine'"):
            _resolve_engine_name(config)

    def test_empty_engine_raises_usage_error(self) -> None:
        """Empty config engine should raise a usage error."""
        from types import SimpleNamespace

        from polyzymd.cli.main import _resolve_engine_name

        config = SimpleNamespace(engine="")
        with pytest.raises(click.UsageError, match="Configure 'engine'"):
            _resolve_engine_name(config)

    def test_case_insensitive(self) -> None:
        """Override should be case-insensitive."""
        from types import SimpleNamespace

        from polyzymd.cli.main import _resolve_engine_name

        config = SimpleNamespace(engine="openmm")
        assert _resolve_engine_name(config, override="GROMACS") == "gromacs"


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
        """'polyzymd build --help' should NOT show removed --replicate."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["build", "--help"])
        assert result.exit_code == 0
        # --replicate is hidden, so it should not appear in help output
        # (but --replicates will appear, so we check there's no standalone --replicate)
        lines = result.output.split("\n")
        for line in lines:
            # Detect removed singular option without matching the plural option
            if "--replicate" in line and "--replicates" not in line:
                pytest.fail(f"Removed singular --replicate visible in help: {line}")

    def test_build_gromacs_alias_is_rejected(self, tmp_path: Path) -> None:
        """Build should reject the removed --gromacs alias as unknown."""
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(cli, ["build", "-c", str(config_path), "--gromacs"])

        assert result.exit_code != 0
        assert "No such option: --gromacs" in result.output

    def test_build_help_clarifies_gromacs_handoff(self) -> None:
        """Build help should describe GROMACS as build-only handoff."""
        runner = CliRunner()

        result = runner.invoke(cli, ["build", "--help"])

        assert result.exit_code == 0
        assert "build-only GROMACS handoff" in result.output

    def test_run_help_keeps_gromacs_engine(self) -> None:
        """Run help should retain GROMACS as a full workflow engine."""
        runner = CliRunner()

        result = runner.invoke(cli, ["run", "--help"])

        assert result.exit_code == 0
        assert "gromacs" in result.output
        assert "openmm" in result.output

    @pytest.mark.parametrize("option", ["--output-dir", "-o"])
    def test_build_output_dir_alias_is_rejected(self, option: str, tmp_path: Path) -> None:
        """Build should reject removed output-dir aliases for scratch output."""
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(cli, ["build", "-c", str(config_path), option, str(tmp_path)])

        assert result.exit_code != 0
        assert f"No such option: {option}" in result.output


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
        """'polyzymd run --help' should NOT show removed --replicate."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--help"])
        assert result.exit_code == 0
        lines = result.output.split("\n")
        for line in lines:
            if "--replicate" in line and "--replicates" not in line:
                pytest.fail(f"Removed singular --replicate visible in help: {line}")

    def test_run_replicate_alias_is_rejected(self, tmp_path: Path) -> None:
        """Run should reject the removed singular --replicate option."""
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(cli, ["run", "-c", str(config_path), "--replicate", "1"])

        assert result.exit_code != 0
        assert "No such option: --replicate" in result.output

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

    def test_submit_help_shows_engine_flag(self) -> None:
        """submit --help should show --engine option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["submit", "--help"])
        assert result.exit_code == 0
        assert "--engine" in result.output

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


class TestInitCommand:
    """Tests for the project initialization scaffold."""

    def test_init_creates_project_scaffold(self, tmp_path: Path) -> None:
        """init should render project files with the requested project name."""
        runner = CliRunner()

        with runner.isolated_filesystem(temp_dir=tmp_path):
            result = runner.invoke(cli, ["init", "-n", "foo"])

            assert result.exit_code == 0
            project_dir = Path("foo")
            assert project_dir.is_dir()
            assert (project_dir / "config.yaml").is_file()
            assert (project_dir / "structures").is_dir()
            assert (project_dir / "job_scripts").is_dir()
            assert (project_dir / "slurm_logs").is_dir()
            assert (project_dir / "structures" / "place_protein_here.placeholder.txt").is_file()
            assert (project_dir / "structures" / "place_ligand_here.placeholder.txt").is_file()

            config_content = (project_dir / "config.yaml").read_text(encoding="utf-8")
            config_data = yaml.safe_load(config_content)
            assert config_data["name"] == "foo"
            assert config_content.count(BRANDING_LINE) == 1
            _assert_no_jinja_markers(config_content)

            for path in project_dir.rglob("*"):
                if path.is_file() and path.suffix in {".yaml", ".txt"}:
                    content = path.read_text(encoding="utf-8")
                    assert content.count(BRANDING_LINE) == 1
                    _assert_no_jinja_markers(content)

    def test_init_existing_directory_still_errors(self, tmp_path: Path) -> None:
        """init should preserve the existing-directory error behavior."""
        runner = CliRunner()

        with runner.isolated_filesystem(temp_dir=tmp_path):
            Path("foo").mkdir()
            result = runner.invoke(cli, ["init", "-n", "foo"])

            assert result.exit_code != 0
            assert "already exists" in result.output

    def test_shared_renderer_uses_strict_undefined(self) -> None:
        """Missing template context values should fail fast."""
        with pytest.raises(UndefinedError):
            render_package_template("polyzymd.templates", "config_template.yaml", {})


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
    def test_build_dry_run_rejects_singular_replicate(
        self,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """Build should reject the removed singular --replicate option."""
        mock_from_yaml.return_value = _make_dry_run_config()
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            ["build", "-c", str(config_path), "--replicate", "1", "--dry-run"],
        )

        assert result.exit_code != 0
        assert "No such option: --replicate" in result.output

    def test_build_with_replicate_zero_fails(self, tmp_path: Path) -> None:
        """Removed singular --replicate should fail as an unknown option."""
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(cli, ["build", "-c", str(config_path), "--replicate", "0"])

        assert result.exit_code != 0
        assert "No such option: --replicate" in result.output


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

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_submit_catches_runtime_error(
        self,
        mock_submit_daisy_chain,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """submit should catch RuntimeError and exit cleanly."""
        mock_from_yaml.return_value = _make_dry_run_config()
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

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_submit_does_not_catch_type_error(
        self,
        mock_submit_daisy_chain,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """submit should not catch TypeError programmer errors."""
        mock_from_yaml.return_value = _make_dry_run_config()
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
        """Removed singular --replicate should fail as an unknown option."""
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(cli, ["build", "-c", str(config_path), "--replicate", "-1"])

        assert result.exit_code != 0
        assert "No such option: --replicate" in result.output


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

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_generate_only_calls_submit_with_flag(
        self,
        mock_submit,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """submit --generate-only should pass generate_only=True to the backend."""
        mock_from_yaml.return_value = _make_dry_run_config()
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


class TestSubmitEngineAware:
    """Tests for engine-aware submit command."""

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.engines.create_engine")
    def test_submit_dry_run_gromacs(
        self,
        mock_create_engine,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """submit --engine gromacs --dry-run should preview without files."""
        mock_engine = SimpleNamespace()
        mock_engine._resolve_slurm_config = lambda base: SimpleNamespace(
            partition="aa100",
            time_limit="04:00:00",
            memory="8G",
            account=None,
            email="",
            nodes=1,
            ntasks=1,
            cpus_per_task=4,
            gpus=1,
            constraint=None,
            qos=None,
            nodelist=None,
        )
        mock_engine._resolve_mdrun_flags = lambda _effective: "-pin on"
        mock_create_engine.return_value = mock_engine

        mock_config = _make_dry_run_config()
        mock_config.engine = "gromacs"
        mock_config.gromacs = SimpleNamespace(
            gmx_binary=None,
            gpu=True,
            ntmpi=1,
            slurm_ntasks=None,
            ntomp=4,
            module_load=None,
        )
        mock_from_yaml.return_value = mock_config
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            ["submit", "-c", str(config_path), "--engine", "gromacs", "--dry-run"],
        )

        assert result.exit_code == 0
        assert "DRY RUN" in result.output
        mock_create_engine.assert_called_once_with(
            mock_config,
            override="gromacs",
            defer_binary=True,
        )

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.engines.create_engine")
    def test_submit_dry_run_gromacs_shows_effective_slurm_details(
        self,
        mock_create_engine,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """GROMACS dry-run should include effective SLURM and engine details."""
        mock_engine = SimpleNamespace()
        mock_engine._resolve_slurm_config = lambda base: SimpleNamespace(
            partition="bridges2",
            time_limit="12:00:00",
            memory="16G",
            account="mcb200029p",
            email="",
            nodes=1,
            ntasks=2,
            cpus_per_task=8,
            gpus=1,
            constraint="A40",
            qos="normal",
            nodelist=None,
        )
        mock_engine._resolve_mdrun_flags = lambda _effective: "-update gpu -bonded gpu"
        mock_create_engine.return_value = mock_engine

        mock_config = _make_dry_run_config()
        mock_config.engine = "gromacs"
        mock_config.gromacs = SimpleNamespace(
            gmx_binary="gmx_mpi",
            gpu=True,
            ntmpi=2,
            slurm_ntasks=None,
            ntomp=8,
            module_load="gromacs/2024.1",
        )
        mock_from_yaml.return_value = mock_config

        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            [
                "submit",
                "-c",
                str(config_path),
                "--engine",
                "gromacs",
                "--dry-run",
                "--preset",
                "bridges2",
            ],
        )

        assert result.exit_code == 0
        assert "Engine:" in result.output
        assert "gromacs" in result.output
        assert "Partition:" in result.output
        assert "Time limit:" in result.output
        assert "Memory:" in result.output
        assert "GPUs:" in result.output
        assert "GPU mode:" in result.output
        assert "ntmpi:" in result.output
        assert "ntomp:" in result.output
        assert "mdrun flags:" in result.output

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.engines.create_engine")
    def test_submit_dry_run_gromacs_shows_partition_qos_overrides(
        self,
        mock_create_engine,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """Dry run for GROMACS should show partition and QoS overrides."""
        mock_engine = SimpleNamespace()
        mock_engine._resolve_slurm_config = lambda base: SimpleNamespace(
            partition=base.partition,
            time_limit=base.time_limit,
            memory=base.memory,
            account=base.account,
            email=base.email,
            nodes=base.nodes,
            ntasks=base.ntasks,
            cpus_per_task=base.cpus_per_task,
            gpus=base.gpus,
            constraint=base.constraint,
            qos=base.qos,
            nodelist=base.nodelist,
        )
        mock_engine._resolve_mdrun_flags = lambda _effective: "-pin on"
        mock_create_engine.return_value = mock_engine

        mock_config = _make_dry_run_config()
        mock_config.engine = "gromacs"
        mock_config.gromacs = SimpleNamespace(
            gmx_binary=None,
            gpu=False,
            ntmpi=1,
            slurm_ntasks=None,
            ntomp=4,
            module_load=None,
        )
        mock_from_yaml.return_value = mock_config

        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            [
                "submit",
                "-c",
                str(config_path),
                "--engine",
                "gromacs",
                "--dry-run",
                "--partition",
                "blanca",
                "--qos",
                "preemptable",
                "--email",
                "user@example.com",
            ],
        )

        assert result.exit_code == 0
        assert "Partition:" in result.output
        assert "blanca" in result.output
        assert "QoS:" in result.output
        assert "preemptable" in result.output
        assert "Email:" in result.output
        assert "user@example.com" in result.output

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_submit_openmm_path_unchanged(
        self, mock_submit, mock_from_yaml, tmp_path: Path
    ) -> None:
        """submit without --engine should still use OpenMM daisy-chain."""
        mock_config = _make_dry_run_config()
        mock_config.engine = "openmm"
        mock_from_yaml.return_value = mock_config
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(cli, ["submit", "-c", str(config_path)])

        assert result.exit_code == 0
        mock_submit.assert_called_once()

    @patch("polyzymd.engines.gromacs.engine.GromacsEngine.submit")
    @patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_submit_gromacs_calls_engine_submit(
        self,
        mock_from_yaml,
        mock_resolve_gmx,
        mock_engine_submit,
        tmp_path: Path,
    ) -> None:
        """submit --engine gromacs should call GromacsEngine.submit()."""
        _ = mock_resolve_gmx
        mock_config = _make_dry_run_config()
        mock_config.engine = "gromacs"
        mock_config.gromacs = SimpleNamespace(
            grompp_flags="",
            mdrun_flags="",
            module_load=None,
            gmx_binary=None,
            ntmpi=1,
            slurm_ntasks=None,
            ntomp=4,
            gpu=False,
            gpus=1,
            memory="16G",
        )
        mock_config.generate_system_name = lambda: "test_system"
        mock_from_yaml.return_value = mock_config
        mock_engine_submit.return_value = {
            "submitted": False,
            "script_path": "/tmp/script.sh",
            "reason": "sbatch_not_available",
        }

        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        with patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=[]):
            with patch("polyzymd.workflow.daisy_chain.create_job_name", return_value="test_job"):
                result = runner.invoke(
                    cli,
                    ["submit", "-c", str(config_path), "--engine", "gromacs"],
                )

        assert result.exit_code == 0
        mock_engine_submit.assert_called_once()

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_submit_openff_logs_warns_for_gromacs(self, mock_from_yaml, tmp_path: Path) -> None:
        """--openff-logs with --engine gromacs should warn."""
        mock_config = _make_dry_run_config()
        mock_config.engine = "gromacs"
        mock_config.gromacs = SimpleNamespace(
            grompp_flags="",
            mdrun_flags="",
            module_load=None,
            gmx_binary=None,
            ntmpi=1,
            slurm_ntasks=None,
            ntomp=4,
            gpu=False,
            gpus=1,
            memory="16G",
        )
        mock_from_yaml.return_value = mock_config
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            [
                "submit",
                "-c",
                str(config_path),
                "--engine",
                "gromacs",
                "--openff-logs",
                "--dry-run",
            ],
        )

        assert result.exit_code == 0
        assert "no effect" in result.output.lower()


class TestSubmitConstraintOption:
    """Tests for --constraint CLI option on submit command."""

    def test_submit_help_shows_nodelist(self) -> None:
        """'polyzymd submit --help' should show --nodelist option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["submit", "--help"])
        assert result.exit_code == 0
        assert "--nodelist" in result.output

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_submit_gpu_type_a100_accepted_by_click(
        self,
        mock_submit,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """submit should accept arbitrary sanitized --gpu-type values like a100."""
        mock_config = _make_dry_run_config()
        mock_config.engine = "openmm"
        mock_from_yaml.return_value = mock_config
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["submit", "-c", str(config_path), "--gpu-type", "a100"],
        )

        assert result.exit_code == 0
        mock_submit.assert_called_once()
        assert mock_submit.call_args.kwargs["gpu_type"] == "a100"

    def test_submit_help_shows_constraint(self) -> None:
        """'polyzymd submit --help' should show --constraint option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["submit", "--help"])
        assert result.exit_code == 0
        assert "--constraint" in result.output

    def test_submit_help_shows_partition_qos_email(self) -> None:
        """'polyzymd submit --help' should show partition/qos/email options."""
        runner = CliRunner()
        result = runner.invoke(cli, ["submit", "--help"])
        assert result.exit_code == 0
        assert "--partition" in result.output
        assert "--qos" in result.output
        assert "--email" in result.output

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_constraint_passed_to_submit_daisy_chain(
        self,
        mock_submit,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """--constraint should be passed to submit_daisy_chain."""
        mock_config = _make_dry_run_config()
        mock_config.engine = "openmm"
        mock_from_yaml.return_value = mock_config
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["submit", "-c", str(config_path), "--constraint", "A40|A100"],
        )

        assert result.exit_code == 0
        mock_submit.assert_called_once()
        assert mock_submit.call_args.kwargs["constraint"] == "A40|A100"

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_partition_qos_email_passed_to_submit_daisy_chain(
        self,
        mock_submit,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """submit should pass partition/qos/email overrides to OpenMM path."""
        mock_config = _make_dry_run_config()
        mock_config.engine = "openmm"
        mock_from_yaml.return_value = mock_config
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "submit",
                "-c",
                str(config_path),
                "--partition",
                "debug",
                "--qos",
                "normal",
                "--email",
                "user@example.com",
            ],
        )

        assert result.exit_code == 0
        mock_submit.assert_called_once()
        kwargs = mock_submit.call_args.kwargs
        assert kwargs["partition"] == "debug"
        assert kwargs["qos"] == "normal"
        assert kwargs["email"] == "user@example.com"

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_nodelist_passed_to_submit_daisy_chain(
        self,
        mock_submit,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """submit should pass nodelist override to OpenMM submission path."""
        mock_config = _make_dry_run_config()
        mock_config.engine = "openmm"
        mock_from_yaml.return_value = mock_config
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "submit",
                "-c",
                str(config_path),
                "--nodelist",
                "bgpu-shirts3",
            ],
        )

        assert result.exit_code == 0
        mock_submit.assert_called_once()
        kwargs = mock_submit.call_args.kwargs
        assert kwargs["nodelist"] == "bgpu-shirts3"

    @patch("polyzymd.engines.gromacs.engine.GromacsEngine.submit")
    @patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_constraint_set_on_gromacs_slurm_config(
        self,
        mock_from_yaml,
        mock_resolve_gmx,
        mock_engine_submit,
        tmp_path: Path,
    ) -> None:
        """--constraint with --engine gromacs should set constraint on SlurmConfig."""
        _ = mock_resolve_gmx
        mock_config = _make_dry_run_config()
        mock_config.engine = "gromacs"
        mock_config.gromacs = SimpleNamespace(
            grompp_flags="", mdrun_flags="", module_load=None, gmx_binary=None
        )
        mock_config.gromacs.slurm_ntasks = None
        mock_config.generate_system_name = lambda: "test_system"
        mock_from_yaml.return_value = mock_config
        mock_engine_submit.return_value = {
            "submitted": False,
            "script_path": "/tmp/script.sh",
            "reason": "sbatch_not_available",
        }

        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        with patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=[]):
            with patch("polyzymd.workflow.daisy_chain.create_job_name", return_value="test_job"):
                result = runner.invoke(
                    cli,
                    [
                        "submit",
                        "-c",
                        str(config_path),
                        "--engine",
                        "gromacs",
                        "--constraint",
                        "A40",
                    ],
                )

        assert result.exit_code == 0
        mock_engine_submit.assert_called_once()
        request = mock_engine_submit.call_args[0][0]
        assert request.slurm_config.constraint == "A40"

    @patch("polyzymd.engines.gromacs.engine.GromacsEngine.submit")
    @patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_partition_qos_email_set_on_gromacs_slurm_config(
        self,
        mock_from_yaml,
        mock_resolve_gmx,
        mock_engine_submit,
        tmp_path: Path,
    ) -> None:
        """submit --engine gromacs should pass partition/qos/email to SlurmConfig."""
        _ = mock_resolve_gmx
        mock_config = _make_dry_run_config()
        mock_config.engine = "gromacs"
        mock_config.gromacs = SimpleNamespace(
            grompp_flags="",
            mdrun_flags="",
            mdrun_flags_equilibration=None,
            mdrun_flags_production=None,
            command_prefix=None,
            mpi_launcher_flags="",
            module_load=None,
            gmx_binary=None,
            ntmpi=1,
            slurm_ntasks=None,
            ntomp=4,
            gpu=False,
            gpus=1,
            memory="16G",
            env_exports={},
            setup_commands=[],
        )
        mock_config.generate_system_name = lambda: "test_system"
        mock_from_yaml.return_value = mock_config
        mock_engine_submit.return_value = {
            "submitted": False,
            "script_path": "/tmp/script.sh",
            "reason": "sbatch_not_available",
        }

        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        with patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=[]):
            with patch("polyzymd.workflow.daisy_chain.create_job_name", return_value="test_job"):
                result = runner.invoke(
                    cli,
                    [
                        "submit",
                        "-c",
                        str(config_path),
                        "--engine",
                        "gromacs",
                        "--partition",
                        "debug",
                        "--qos",
                        "normal",
                        "--email",
                        "user@example.com",
                    ],
                )

        assert result.exit_code == 0
        mock_engine_submit.assert_called_once()
        request = mock_engine_submit.call_args[0][0]
        assert request.slurm_config.partition == "debug"
        assert request.slurm_config.qos == "normal"
        assert request.slurm_config.email == "user@example.com"

    @patch("polyzymd.engines.gromacs.engine.GromacsEngine.submit")
    @patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_nodelist_set_on_gromacs_slurm_config(
        self,
        mock_from_yaml,
        mock_resolve_gmx,
        mock_engine_submit,
        tmp_path: Path,
    ) -> None:
        """submit --engine gromacs should pass nodelist override to SlurmConfig."""
        _ = mock_resolve_gmx
        mock_config = _make_dry_run_config()
        mock_config.engine = "gromacs"
        mock_config.gromacs = SimpleNamespace(
            grompp_flags="",
            mdrun_flags="",
            mdrun_flags_equilibration=None,
            mdrun_flags_production=None,
            command_prefix=None,
            mpi_launcher_flags="",
            module_load=None,
            gmx_binary=None,
            ntmpi=1,
            slurm_ntasks=None,
            ntomp=4,
            gpu=False,
            gpus=1,
            memory="16G",
            env_exports={},
            setup_commands=[],
        )
        mock_config.generate_system_name = lambda: "test_system"
        mock_from_yaml.return_value = mock_config
        mock_engine_submit.return_value = {
            "submitted": False,
            "script_path": "/tmp/script.sh",
            "reason": "sbatch_not_available",
        }

        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        with patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=[]):
            with patch("polyzymd.workflow.daisy_chain.create_job_name", return_value="test_job"):
                result = runner.invoke(
                    cli,
                    [
                        "submit",
                        "-c",
                        str(config_path),
                        "--engine",
                        "gromacs",
                        "--nodelist",
                        "bgpu-shirts3",
                    ],
                )

        assert result.exit_code == 0
        mock_engine_submit.assert_called_once()
        request = mock_engine_submit.call_args[0][0]
        assert request.slurm_config.nodelist == "bgpu-shirts3"

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.engines.create_engine")
    def test_submit_dry_run_gromacs_displays_nodelist_override(
        self,
        mock_create_engine,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """GROMACS dry-run output should display nodelist when provided."""
        mock_engine = SimpleNamespace()
        mock_engine._resolve_slurm_config = lambda base: SimpleNamespace(
            partition=base.partition,
            time_limit=base.time_limit,
            memory=base.memory,
            account=base.account,
            email=base.email,
            nodes=base.nodes,
            ntasks=base.ntasks,
            cpus_per_task=base.cpus_per_task,
            gpus=base.gpus,
            constraint=base.constraint,
            qos=base.qos,
            nodelist=base.nodelist,
        )
        mock_engine._resolve_mdrun_flags = lambda _effective: "-pin on"
        mock_create_engine.return_value = mock_engine

        mock_config = _make_dry_run_config()
        mock_config.engine = "gromacs"
        mock_config.gromacs = SimpleNamespace(
            gmx_binary=None,
            gpu=False,
            ntmpi=1,
            slurm_ntasks=None,
            ntomp=4,
            module_load=None,
        )
        mock_from_yaml.return_value = mock_config

        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            [
                "submit",
                "-c",
                str(config_path),
                "--engine",
                "gromacs",
                "--dry-run",
                "--nodelist",
                "bgpu-shirts3",
            ],
        )

        assert result.exit_code == 0
        assert "Nodelist:" in result.output
        assert "bgpu-shirts3" in result.output

    def test_recover_help_shows_constraint(self) -> None:
        """'polyzymd recover --help' should show --constraint option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["recover", "--help"])
        assert result.exit_code == 0
        assert "--constraint" in result.output

    def test_recover_help_shows_email(self) -> None:
        """'polyzymd recover --help' should show --email option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["recover", "--help"])
        assert result.exit_code == 0
        assert "--email" in result.output

    def test_recover_help_shows_engine(self) -> None:
        """'polyzymd recover --help' should show --engine option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["recover", "--help"])
        assert result.exit_code == 0
        assert "--engine" in result.output


class TestSubmitGromacsDuplicateGuard:
    """Tests for duplicate-job detection in GROMACS submit path."""

    @patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=["12345"])
    @patch("polyzymd.workflow.daisy_chain.create_job_name", return_value="test_job")
    @patch("polyzymd.engines.gromacs.engine.GromacsEngine.submit")
    @patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_submit_gromacs_blocked_by_existing_job(
        self,
        mock_from_yaml,
        mock_resolve,
        mock_submit,
        mock_create_name,
        mock_check,
        tmp_path,
    ):
        """submit --engine gromacs should be blocked when a job already exists."""
        _ = mock_resolve, mock_create_name, mock_check
        mock_config = _make_dry_run_config()
        mock_config.engine = "gromacs"
        mock_config.gromacs = SimpleNamespace(
            grompp_flags="",
            mdrun_flags="",
            module_load=None,
            gmx_binary=None,
        )
        mock_config.gromacs.slurm_ntasks = None
        mock_config.generate_system_name = lambda: "test_system"
        mock_from_yaml.return_value = mock_config

        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            ["submit", "-c", str(config_path), "--engine", "gromacs"],
        )

        assert result.exit_code == 0
        assert "already has RUNNING/PENDING" in result.output
        mock_submit.assert_not_called()

    @patch("polyzymd.engines.gromacs.engine.GromacsEngine.submit")
    @patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=["12345"])
    @patch("polyzymd.workflow.daisy_chain.create_job_name", return_value="test_job")
    @patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_submit_gromacs_force_bypasses_guard(
        self,
        mock_from_yaml,
        mock_resolve,
        mock_create_name,
        mock_check,
        mock_submit,
        tmp_path,
    ):
        """submit --engine gromacs --force should bypass duplicate guard."""
        _ = mock_resolve, mock_create_name, mock_check
        mock_config = _make_dry_run_config()
        mock_config.engine = "gromacs"
        mock_config.gromacs = SimpleNamespace(
            grompp_flags="",
            mdrun_flags="",
            module_load=None,
            gmx_binary=None,
        )
        mock_config.gromacs.slurm_ntasks = None
        mock_config.generate_system_name = lambda: "test_system"
        mock_from_yaml.return_value = mock_config
        mock_submit.return_value = {
            "submitted": False,
            "script_path": "/tmp/script.sh",
            "reason": "sbatch_not_available",
        }

        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n")
        runner = CliRunner()

        result = runner.invoke(
            cli,
            ["submit", "-c", str(config_path), "--engine", "gromacs", "--force"],
        )
        assert result.exit_code == 0
        mock_submit.assert_called_once()
