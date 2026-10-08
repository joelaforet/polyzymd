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

from polyzymd.cli.main import _resolve_replicates_option, _run_openmm_impl, cli
from polyzymd.config.schema import Ensemble, WaterModel
from polyzymd.utils.templates import render_package_template


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
            primary=SimpleNamespace(model=WaterModel.TIP3P),
            box=SimpleNamespace(padding=1.2),
            ions=SimpleNamespace(nacl_concentration=0.15, neutralize=True),
            co_solvents=[],
        ),
        force_field=SimpleNamespace(protein="amber14", small_molecule="openff-2.2.0"),
        thermodynamics=SimpleNamespace(temperature=300.0, pressure=1.0),
        simulation_phases=SimpleNamespace(
            equilibration_stages=[
                SimpleNamespace(
                    name="heating",
                    resolved_duration=0.2,
                    ensemble=Ensemble.NVT,
                    is_temperature_ramping=False,
                ),
                SimpleNamespace(
                    name="free_equilibration",
                    resolved_duration=0.8,
                    ensemble=Ensemble.NPT,
                    is_temperature_ramping=False,
                ),
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


def _minimal_cli_config_data(pdb_path: str | Path) -> dict[str, object]:
    """Create minimal YAML-serializable config data for CLI tests."""

    return {
        "name": "cli_reference_validation",
        "engine": "openmm",
        "enzyme": {"name": "Enz", "pdb_path": str(pdb_path)},
        "thermodynamics": {"temperature": 300.0},
        "simulation_phases": {
            "equilibration_stages": [
                {
                    "name": "eq",
                    "duration": 0.1,
                    "temperature": 300.0,
                    "ensemble": "NVT",
                }
            ],
            "production": {
                "ensemble": "NPT",
                "duration": 1.0,
                "samples": 10,
                "checkpoint_interval": 60.0,
            },
        },
    }


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


class TestResolveSubmissionPixiEnv:
    """Unit tests for submit/recover pixi environment resolution."""

    def test_explicit_env_takes_priority(self) -> None:
        """An explicit --pixi-env should be preserved for any engine."""
        from polyzymd.cli.main import _resolve_submission_pixi_env

        assert _resolve_submission_pixi_env("aa100", "gromacs", "sim-cuda-12-6") == "sim-cuda-12-6"

    def test_gromacs_defaults_to_build(self) -> None:
        """GROMACS should use build when no explicit pixi env is provided."""
        from polyzymd.cli.main import _resolve_submission_pixi_env

        assert _resolve_submission_pixi_env("bridges2", "gromacs") == "build"

    def test_gromacs_auto_maps_to_build(self) -> None:
        """Shared auto selection must not create a nonexistent GROMACS environment."""
        from polyzymd.cli.main import _resolve_submission_pixi_env

        assert _resolve_submission_pixi_env("bridges2", "gromacs", "auto") == "build"

    def test_openmm_defaults_from_preset(self) -> None:
        """OpenMM should keep preset-specific CUDA environment defaults."""
        from polyzymd.cli.main import _resolve_submission_pixi_env

        assert _resolve_submission_pixi_env("blanca-shirts", "openmm") == "sim-cuda-12-4"
        assert _resolve_submission_pixi_env("blanca-chbe-rdi", "openmm") == "sim-cuda-12-4"
        assert _resolve_submission_pixi_env("bridges2", "openmm") == "sim-cuda-12-6"

    def test_openmm_auto_resolves_from_site_policy(self) -> None:
        """Explicit auto should use a known site's fixed environment."""
        from polyzymd.cli.main import _resolve_submission_pixi_env

        assert _resolve_submission_pixi_env("blanca-shirts", "openmm", "auto") == "sim-cuda-12-4"
        assert _resolve_submission_pixi_env("bridges2", "openmm", "auto") == "sim-cuda-12-6"

    def test_unpinned_site_retains_capability_auto(self) -> None:
        """A site without validated policy should retain runtime detection."""
        from polyzymd.cli.main import _resolve_submission_pixi_env

        assert _resolve_submission_pixi_env("aa100", "openmm", "auto") == "auto"

    def test_explicit_openmm_override_is_preserved(self) -> None:
        """An expert override should not be replaced by site policy."""
        from polyzymd.cli.main import _resolve_submission_pixi_env

        assert (
            _resolve_submission_pixi_env("blanca-shirts", "openmm", "sim-cuda-12-6")
            == "sim-cuda-12-6"
        )

    def test_site_override_prints_reproducibility_warning(self, monkeypatch) -> None:
        """An expert override should identify the preset environment policy."""
        from polyzymd.cli.main import _warn_for_site_pixi_override

        messages: list[str] = []
        monkeypatch.setattr(click, "secho", lambda message, **_kwargs: messages.append(message))

        _warn_for_site_pixi_override("blanca-shirts", "openmm", "sim-cuda-12-6", "sim-cuda-12-6")

        assert messages
        assert "sim-cuda-12-4" in messages[0]
        assert "all segments of a replica" in messages[0]

    def test_auto_does_not_print_override_warning(self, monkeypatch) -> None:
        """Site auto resolution is normal and should not warn as an override."""
        from polyzymd.cli.main import _warn_for_site_pixi_override

        messages: list[str] = []
        monkeypatch.setattr(click, "secho", lambda message, **_kwargs: messages.append(message))

        _warn_for_site_pixi_override("blanca-shirts", "openmm", "sim-cuda-12-4", "auto")

        assert messages == []


class TestValidateCommandReferenceWarnings:
    """Tests for validate command runtime reference warnings."""

    def test_validate_exits_zero_and_warns_for_missing_referenced_files(
        self,
        tmp_path: Path,
    ) -> None:
        """Validate should warn about missing PDB files without failing schema validation."""

        config_path = tmp_path / "config.yaml"
        config_path.write_text(
            yaml.safe_dump(_minimal_cli_config_data("missing.pdb")),
            encoding="utf-8",
        )
        runner = CliRunner()

        result = runner.invoke(cli, ["validate", "-c", str(config_path)])

        assert result.exit_code == 0
        assert "Configuration is valid!" in result.output
        assert "Referenced file warnings" in result.output
        assert "Missing enzyme PDB" in result.output

    def test_validate_reports_derived_temperature_ramp_duration(self, tmp_path: Path) -> None:
        """Validation clearly reports rate-based heating duration."""
        data = _minimal_cli_config_data("missing.pdb")
        data["simulation_phases"]["equilibration_stages"] = [
            {
                "name": "heating",
                "ensemble": "NVT",
                "temperature_start": 60.0,
                "temperature_end": 300.0,
                "temperature_increment": 1.0,
                "temperature_interval_steps": 600,
                "time_step": 2.0,
            }
        ]
        config_path = tmp_path / "config.yaml"
        config_path.write_text(yaml.safe_dump(data), encoding="utf-8")

        result = CliRunner().invoke(cli, ["validate", "-c", str(config_path)])

        assert result.exit_code == 0
        assert "60 -> 300 K, +1 K every 600 steps" in result.output
        assert "derived duration 0.288000 ns" in result.output

    @pytest.mark.parametrize("where", ["production", "equilibration"])
    def test_validate_refuses_anisotropic_barostat_on_openmm(self, tmp_path: Path, where) -> None:
        """OpenMM runs no anisotropic barostat, so MCA with engine openmm is an error."""
        data = _minimal_cli_config_data("missing.pdb")
        phases = data["simulation_phases"]
        phase = phases["production"] if where == "production" else phases["equilibration_stages"][0]
        phase.update(ensemble="NPT", barostat="MCA")
        config_path = tmp_path / "config.yaml"
        config_path.write_text(yaml.safe_dump(data), encoding="utf-8")

        result = CliRunner().invoke(cli, ["validate", "-c", str(config_path)])

        assert result.exit_code != 0
        assert "MCA" in result.output and "anisotropic" in result.output

    def test_missing_equilibration_stages_names_the_key_without_a_pydantic_url(
        self, tmp_path: Path
    ) -> None:
        data = _minimal_cli_config_data("missing.pdb")
        del data["simulation_phases"]["equilibration_stages"]
        config_path = tmp_path / "config.yaml"
        config_path.write_text(yaml.safe_dump(data), encoding="utf-8")

        result = CliRunner().invoke(cli, ["validate", "-c", str(config_path)])

        assert result.exit_code == 1
        assert "equilibration_stages is missing: list at least one stage" in result.output
        assert "legacy" not in result.output
        assert "errors.pydantic.dev" not in result.output

    def test_validate_prints_engine_and_cosolvents(self, tmp_path: Path) -> None:
        data = _minimal_cli_config_data("missing.pdb")
        data["solvent"] = {
            "primary": {"type": "water", "model": "tip3p"},
            "co_solvents": [{"name": "dmso", "mole_fraction": 0.1}],
        }
        config_path = tmp_path / "config.yaml"
        config_path.write_text(yaml.safe_dump(data), encoding="utf-8")

        result = CliRunner().invoke(cli, ["validate", "-c", str(config_path)])

        assert result.exit_code == 0, result.output
        assert "Engine: openmm" in result.output
        assert "Co-solvents: dmso" in result.output


@pytest.mark.parametrize("export_format", ["lammps", "amber"])
def test_build_refuses_unimplemented_export_formats(tmp_path: Path, export_format: str) -> None:
    """Only implemented export formats are offered by build --format."""
    config = tmp_path / "config.yaml"
    config.write_text("{}\n")

    result = CliRunner().invoke(cli, ["build", "-c", str(config), "--format", export_format])

    assert result.exit_code == 2
    assert "Invalid value for '--format'" in result.output


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

    def test_build_dry_run_warns_for_missing_referenced_files(self, tmp_path: Path) -> None:
        """Build dry-run should warn when schema-valid referenced files are absent."""

        config_path = tmp_path / "config.yaml"
        config_path.write_text(
            yaml.safe_dump(_minimal_cli_config_data("missing.pdb")),
            encoding="utf-8",
        )
        runner = CliRunner()

        result = runner.invoke(cli, ["build", "-c", str(config_path), "--dry-run"])

        assert result.exit_code == 0
        assert "Referenced file warnings" in result.output
        assert "Missing enzyme PDB" in result.output

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


class TestTemplates:
    """Tests for the package template renderer."""

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


class TestOpenMMRunImplementation:
    """Regression tests for local OpenMM production configuration."""

    @patch("polyzymd.cli.main._run_initial_segment")
    def test_derives_report_interval_from_samples(
        self, run_initial_segment, tmp_path: Path
    ) -> None:
        """Local OpenMM runs use samples as the trajectory frame control."""
        production = SimpleNamespace(
            duration=1000.0,
            samples=2500,
            time_step=2.0,
            checkpoint_interval=60.0,
        )
        config = SimpleNamespace(
            simulation_phases=SimpleNamespace(production=production),
            require_engine_barostats=lambda engine: None,
            get_working_directory=lambda replicate: tmp_path / f"run_{replicate}",
        )

        _run_openmm_impl(config, replicate=1)

        run_initial_segment.assert_called_once_with(
            sim_config=config,
            working_dir=tmp_path / "run_1",
            replicate=1,
            skip_build=False,
            duration_ns=1000.0,
            num_samples=2500,
            timestep_fs=2.0,
            report_interval=200000,
            checkpoint_interval_s=60.0,
        )


class TestRunBuildManifest:
    """The build of `polyzymd run` records its provenance like `polyzymd build`."""

    @patch("polyzymd.simulation.runner.SimulationRunner")
    @patch("polyzymd.builders.system_builder.SystemBuilder.from_config")
    def test_run_build_manifest_records_packmol_seed(
        self, from_config, _runner, tmp_path: Path
    ) -> None:
        import json
        from unittest.mock import MagicMock

        from openmm import System, Vec3, unit
        from openmm.app import Element, Topology

        from polyzymd.cli.main import _run_initial_segment
        from polyzymd.config.loader import load_config

        topology = Topology()
        residue = topology.addResidue("HOH", topology.addChain("A"))
        system = System()
        for index in range(2):
            topology.addAtom(f"H{index}", Element.getByAtomicNumber(1), residue)
            system.addParticle(1.0)
        positions = [Vec3(float(index), 0, 0) for index in range(2)] * unit.nanometer
        builder = MagicMock(build_provenance={"polymer_packmol_seed": 1})
        builder.get_openmm_components.return_value = (topology, system, positions)
        from_config.return_value = builder

        config = load_config(Path(__file__).parents[2] / "examples/quickstart/config.yaml")
        _run_initial_segment(
            sim_config=config,
            working_dir=tmp_path,
            replicate=1,
            skip_build=False,
            duration_ns=0.004,
            num_samples=4,
            timestep_fs=2.0,
            report_interval=1,
            checkpoint_interval_s=60.0,
        )

        manifest = json.loads((tmp_path / "build_manifest.json").read_text())
        assert manifest["provenance"] == {"polymer_packmol_seed": 1}


class TestRunReusesBuild:
    """`polyzymd run` reuses the build of an earlier `polyzymd build`."""

    @staticmethod
    def _config(tmp_path: Path) -> SimpleNamespace:
        production = SimpleNamespace(
            duration=0.004, samples=4, time_step=2.0, checkpoint_interval=60.0
        )
        return SimpleNamespace(
            simulation_phases=SimpleNamespace(production=production),
            require_engine_barostats=lambda engine: None,
            get_working_directory=lambda replicate: tmp_path / f"run_{replicate}",
        )

    @patch("polyzymd.cli.main._run_initial_segment")
    @patch("polyzymd.simulation.artifact_integrity.validate_build_bundle")
    def test_openmm_run_reuses_a_valid_build(
        self, validate, run_initial_segment, tmp_path: Path, capsys
    ) -> None:
        config = self._config(tmp_path)
        _run_openmm_impl(config, replicate=1)

        validate.assert_called_once_with(tmp_path / "run_1", config)
        assert run_initial_segment.call_args.kwargs["skip_build"] is True
        assert f"Reusing the build in {tmp_path / 'run_1'}" in capsys.readouterr().out

    @patch("polyzymd.cli.main._run_initial_segment")
    def test_openmm_run_builds_without_a_build(
        self, run_initial_segment, tmp_path: Path, capsys
    ) -> None:
        _run_openmm_impl(self._config(tmp_path), replicate=1)

        assert run_initial_segment.call_args.kwargs["skip_build"] is False
        assert "Building the system in" in capsys.readouterr().out

    @patch("polyzymd.exporters.gromacs.GromacsRunner")
    @patch("polyzymd.builders.system_builder.SystemBuilder.from_config")
    @patch("polyzymd.analyses.shared.gromacs.system_prefix", return_value="sys")
    def test_gromacs_run_reuses_exported_files(
        self, _prefix, from_config, gromacs_runner, tmp_path: Path, capsys
    ) -> None:
        from polyzymd.cli.main import _run_gromacs_impl

        gromacs_dir = tmp_path / "run_1" / "gromacs"
        gromacs_dir.mkdir(parents=True)
        for name in ("sys.top", "sys.gro", "em.mdp", "eq_00_nvt.mdp", "prod.mdp"):
            (gromacs_dir / name).write_text("")

        _run_gromacs_impl(self._config(tmp_path), replicate=1, gmx_path="gmx")

        from_config.assert_not_called()
        gromacs_runner.assert_called_once_with(
            working_dir=gromacs_dir,
            prefix="sys",
            equilibration_mdps=["eq_00_nvt.mdp"],
            gmx_command="gmx",
        )
        assert f"Reusing the GROMACS files in {gromacs_dir}" in capsys.readouterr().out

    @patch("polyzymd.exporters.gromacs.GromacsRunner")
    @patch("polyzymd.builders.system_builder.SystemBuilder.from_config")
    @patch("polyzymd.analyses.shared.gromacs.system_prefix", return_value="sys")
    def test_gromacs_run_does_not_rebuild_a_started_run(
        self, _prefix, from_config, gromacs_runner, tmp_path: Path
    ) -> None:
        from polyzymd.cli.main import _run_gromacs_impl
        from polyzymd.simulation.artifact_integrity import ArtifactIntegrityError

        gromacs_dir = tmp_path / "run_1" / "gromacs"
        gromacs_dir.mkdir(parents=True)
        (gromacs_dir / "em.log").write_text("")

        with pytest.raises(ArtifactIntegrityError, match="Refusing to rebuild"):
            _run_gromacs_impl(self._config(tmp_path), replicate=1, gmx_path="gmx")

        from_config.assert_not_called()
        gromacs_runner.assert_not_called()

    @patch("polyzymd.exporters.gromacs.GromacsRunner")
    @patch("polyzymd.analyses.shared.gromacs.system_prefix", return_value="sys")
    def test_gromacs_run_records_the_polyzymd_version(
        self, _prefix, gromacs_runner, tmp_path: Path
    ) -> None:
        from polyzymd import __version__
        from polyzymd.cli.main import _run_gromacs_impl
        from polyzymd.simulation.progress import load_progress

        gromacs_dir = tmp_path / "run_1" / "gromacs"
        gromacs_dir.mkdir(parents=True)
        for name in ("sys.top", "sys.gro", "em.mdp", "eq_01_nvt.mdp", "eq_01.gro", "prod.mdp"):
            (gromacs_dir / name).write_text("")
        (gromacs_dir / "eq_02.gro").write_text("")
        # eq_02 ran before this run; this run skipped it.
        (gromacs_dir / "eq_02.log").write_text("Started mdrun on rank 0 Wed Oct  7 00:00:00 2020\n")

        _run_gromacs_impl(self._config(tmp_path), replicate=1, gmx_path="gmx")

        first, earlier = load_progress(gromacs_dir).equilibration_stages
        assert first.polyzymd_version == __version__
        assert earlier.polyzymd_version is None


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

    @pytest.fixture(autouse=True)
    def _gromacs_builds_exist(self, monkeypatch):
        """These tests check submission options, not the build check."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.engine.GromacsEngine.check_build", lambda self, request: None
        )

    def test_submit_help_includes_build_pixi_env_choice(self) -> None:
        """submit --help should list build as an allowed pixi environment."""
        runner = CliRunner()

        result = runner.invoke(cli, ["submit", "--help"])

        assert result.exit_code == 0
        assert "build" in result.output
        assert "sim-cuda-12-4" in result.output
        assert "sim-cuda-12-6" in result.output

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
            exclude=None,
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
            command_prefix=None,
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
        assert "Pixi env:    build" in result.output
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
            exclude=None,
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
            exclude=base.exclude,
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
            command_prefix=None,
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

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_submit_refuses_a_replicate_stopped_by_cancel(
        self, mock_submit, mock_from_yaml, tmp_path: Path
    ) -> None:
        """A job submitted while STOP exists would exit at once, so submit says so instead."""
        run_dir = tmp_path / "run_2"
        run_dir.mkdir()
        (run_dir / "STOP").write_text("stopped\n")
        mock_config = _make_dry_run_config()
        mock_config.engine = "openmm"
        mock_config.get_working_directory = lambda rep: tmp_path / f"run_{rep}"
        mock_from_yaml.return_value = mock_config
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")

        result = CliRunner().invoke(cli, ["submit", "-c", str(config_path), "-r", "1-2"])

        assert result.exit_code == 1
        assert f"polyzymd cancel -c {config_path} -r 2 --resume" in result.output
        mock_submit.assert_not_called()

    @pytest.mark.parametrize("limit", ["0:05:00", "5", "0:04:30"])
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_submit_refuses_time_limit_within_the_stop_signal_margin(
        self, mock_submit, mock_from_yaml, limit, tmp_path: Path
    ) -> None:
        """SLURM signals OpenMM jobs 5 minutes before the limit; a shorter limit never runs a step."""
        mock_config = _make_dry_run_config()
        mock_config.engine = "openmm"
        mock_from_yaml.return_value = mock_config
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")

        result = CliRunner().invoke(
            cli, ["submit", "-c", str(config_path), "--time-limit", limit]
        )

        assert result.exit_code == 2
        assert "--time-limit" in result.output
        mock_submit.assert_not_called()

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
            "submitted": True,
            "script_path": "/tmp/script.sh",
            "stdout": "Submitted batch job 1",
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
        request = mock_engine_submit.call_args.args[0]
        assert request.extra["pixi_env"] == "build"

    @patch("polyzymd.engines.gromacs.engine.GromacsEngine.submit")
    @patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_submit_gromacs_allows_explicit_cuda_pixi_env(
        self,
        mock_from_yaml,
        mock_resolve_gmx,
        mock_engine_submit,
        tmp_path: Path,
    ) -> None:
        """GROMACS submit should preserve explicit CUDA pixi env overrides."""
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
            gpu=True,
            gpus=1,
            memory="16G",
        )
        mock_from_yaml.return_value = mock_config
        mock_engine_submit.return_value = {
            "submitted": True,
            "script_path": "/tmp/script.sh",
            "stdout": "Submitted batch job 1",
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
                        "--pixi-env",
                        "sim-cuda-12-6",
                    ],
                )

        assert result.exit_code == 0
        request = mock_engine_submit.call_args.args[0]
        assert request.extra["pixi_env"] == "sim-cuda-12-6"

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_submit_openff_logs_warns_for_gromacs(self, mock_from_yaml, tmp_path: Path) -> None:
        """--openff-logs with --engine gromacs should warn."""
        mock_config = _make_dry_run_config()
        mock_config.engine = "gromacs"
        mock_config.gromacs = SimpleNamespace(
            grompp_flags="",
            mdrun_flags="",
            module_load=None,
            command_prefix=None,
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

    @pytest.fixture(autouse=True)
    def _gromacs_builds_exist(self, monkeypatch):
        """These tests check submission options, not the build check."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.engine.GromacsEngine.check_build", lambda self, request: None
        )

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
            "submitted": True,
            "script_path": "/tmp/script.sh",
            "stdout": "Submitted batch job 1",
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
            "submitted": True,
            "script_path": "/tmp/script.sh",
            "stdout": "Submitted batch job 1",
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

    def test_submit_help_shows_exclude(self) -> None:
        """'polyzymd submit --help' should show the --exclude option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["submit", "--help"])
        assert result.exit_code == 0
        assert "--exclude" in result.output

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_exclude_passed_to_submit_daisy_chain(
        self,
        mock_submit,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """submit --exclude should reach the OpenMM submission path."""
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
                "--exclude",
                "bgpu-g4-u20,bgpu-g4-u24",
            ],
        )

        assert result.exit_code == 0
        mock_submit.assert_called_once()
        assert mock_submit.call_args.kwargs["exclude"] == "bgpu-g4-u20,bgpu-g4-u24"

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.workflow.daisy_chain.submit_daisy_chain")
    def test_exclude_defaults_to_none_so_preset_wins(
        self,
        mock_submit,
        mock_from_yaml,
        tmp_path: Path,
    ) -> None:
        """Without --exclude the preset's excluded-node list is left alone."""
        mock_config = _make_dry_run_config()
        mock_config.engine = "openmm"
        mock_from_yaml.return_value = mock_config
        config_path = tmp_path / "fake.yaml"
        config_path.write_text("name: test\n", encoding="utf-8")

        runner = CliRunner()
        result = runner.invoke(cli, ["submit", "-c", str(config_path)])

        assert result.exit_code == 0
        assert mock_submit.call_args.kwargs["exclude"] is None

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
            "submitted": True,
            "script_path": "/tmp/script.sh",
            "stdout": "Submitted batch job 1",
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
            exclude=base.exclude,
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
            command_prefix=None,
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

    @pytest.fixture(autouse=True)
    def _gromacs_builds_exist(self, monkeypatch):
        """These tests check submission options, not the build check."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.engine.GromacsEngine.check_build", lambda self, request: None
        )

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
            "submitted": True,
            "script_path": "/tmp/script.sh",
            "stdout": "Submitted batch job 1",
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


def test_build_follows_the_config_engine(tmp_path: Path) -> None:
    """A GROMACS config builds GROMACS inputs, and run takes its engine from the config."""
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    (tmp_path / "c" / "test.pdb").write_text("END\n")
    data = yaml.safe_load(path.read_text())
    data["engine"] = "gromacs"
    data["solvent"] = {
        "co_solvents": [{"name": "sds", "smiles": "CCCCCCCCCCCCOS(=O)(=O)[O-]", "count": 8}]
    }
    path.write_text(yaml.safe_dump(data))
    result = CliRunner().invoke(cli, ["build", "-c", str(path), "--dry-run"])
    assert "Files to Generate (GROMACS)" in result.output, result.output
    assert "Co-solvent sds (SDS): 8 molecules" in result.output
    assert "Polymer seeds" not in result.output
    dry = CliRunner().invoke(cli, ["run", "-c", str(path), "--dry-run"])
    assert "Missing option '--engine'" not in dry.output


def test_clean_pdb_runs_on_the_cpu(tmp_path: Path, monkeypatch) -> None:
    """clean-pdb places hydrogens on the CPU platform, so a GPU driver mismatch cannot stop it."""
    import os

    pytest.importorskip("pdbfixer")

    monkeypatch.delenv("OPENMM_DEFAULT_PLATFORM", raising=False)
    source = Path(__file__).resolve().parents[2] / "examples" / "quickstart" / "trpcage.pdb"
    result = CliRunner().invoke(
        cli, ["clean-pdb", "-i", str(source), "-o", str(tmp_path / "clean.pdb")]
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "clean.pdb").is_file()
    assert os.environ["OPENMM_DEFAULT_PLATFORM"] == "CPU"


@pytest.mark.parametrize("text", ["garbage\n", ""], ids=["not a PDB", "empty"])
def test_clean_pdb_refuses_a_file_without_atoms(tmp_path: Path, text: str) -> None:
    """A PDB file with no atoms gives error: and fix:, not a traceback."""
    pytest.importorskip("pdbfixer")
    source = tmp_path / "in.pdb"
    source.write_text(text)
    result = CliRunner().invoke(cli, ["clean-pdb", "-i", str(source)])
    assert result.exit_code == 1, result.output
    assert isinstance(result.exception, SystemExit)
    assert f"error: {source} has no atoms" in result.output and "fix:" in result.output
    assert not (tmp_path / "in_clean.pdb").exists()


class TestSubmitDryRunHardwareWarnings:
    """The submit dry run warns when the job's hardware does not fit the config."""

    QUICKSTART = Path(__file__).resolve().parents[2] / "examples" / "quickstart"

    def _dry_run(self, config_name: str) -> str:
        result = CliRunner().invoke(
            cli, ["submit", "-c", str(self.QUICKSTART / config_name), "--dry-run"]
        )
        assert result.exit_code == 0, result.output
        return result.output

    def test_cpu_openmm_config_is_told_job_scripts_need_cuda(self) -> None:
        output = self._dry_run("config.yaml")
        assert "generated OpenMM job scripts need CUDA" in output
        assert "run on a CPU with `polyzymd run`" in output
        assert "choose a CPU partition" not in output

    def test_cpu_gromacs_config_is_warned_about_gpu_partition_and_gmx(self) -> None:
        output = self._dry_run("config_gromacs.yaml")
        assert "the job asks for no GPU, but preset aa100 is set up for GPU jobs" in output
        assert "pass --partition with a partition that has CPU nodes" in output
        assert "choose a CPU partition" not in output
        assert "neither gromacs.module_load nor gromacs.command_prefix is set" in output


class TestNoPolymerWording:
    """Runs without polymers do not mention polymers; console text shows plain values."""

    QUICKSTART = Path(__file__).resolve().parents[2] / "examples" / "quickstart"

    def test_build_dry_run_shows_plain_values_and_no_polymer(self) -> None:
        result = CliRunner().invoke(
            cli, ["build", "-c", str(self.QUICKSTART / "config.yaml"), "--dry-run"]
        )
        assert result.exit_code == 0, result.output
        assert "TIP3P" in result.output
        assert "WaterModel." not in result.output
        assert "Ensemble." not in result.output
        assert "Polymer" not in result.output

    def test_info_reports_packmol_and_gmx(self) -> None:
        result = CliRunner().invoke(cli, ["info"])
        assert result.exit_code == 0, result.output
        assert "packmol:" in result.output
        assert "gmx:" in result.output
        assert "Enzyme-Polymer" not in result.output


@pytest.mark.parametrize("command", ["compare", "new-analysis"])
def test_removed_commands_are_unknown(command: str) -> None:
    """polyzymd compare and polyzymd new-analysis are not commands."""
    result = CliRunner().invoke(cli, [command, "run"])

    assert result.exit_code == 2
    assert f"No such command '{command}'" in result.output


def _write_elongated_pdb(path: Path) -> None:
    """Six atoms spanning 37.5 x 41.8 x 50 Angstrom, long along z."""
    coords = [
        (-18.75, 0.0, 0.0),
        (18.75, 0.0, 0.0),
        (0.0, -20.9, 0.0),
        (0.0, 20.9, 0.0),
        (0.0, 0.0, -25.0),
        (0.0, 0.0, 25.0),
    ]
    lines = [
        f"ATOM  {i + 1:5d}  CA  ALA A{i + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        for i, (x, y, z) in enumerate(coords)
    ]
    path.write_text("\n".join(lines) + "\nEND\n")


@pytest.mark.parametrize("command", ["build", "run"])
def test_dry_run_reports_the_box_and_its_clearances(tmp_path: Path, command: str) -> None:
    """Dry runs print the box the build will make, so a box that is too small shows before Packmol."""
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    _write_elongated_pdb(tmp_path / "c" / "test.pdb")

    result = CliRunner().invoke(cli, [command, "-c", str(path), "--dry-run"])

    assert result.exit_code == 0, result.output
    # diameter 5.00 nm + 2 x 1.2 nm, grown so the 5 nm z extent fits the 0.707-edge brick
    assert "Box (rhombic_dodecahedron): edge 7.64 nm" in result.output, result.output
    assert "clearance to the brick faces 1.94 / 1.73 / 0.20 nm" in result.output


def test_dry_run_box_report_never_blocks(tmp_path: Path) -> None:
    """A box that cannot be estimated gives a note, not a failed dry run."""
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    _write_elongated_pdb(tmp_path / "c" / "test.pdb")

    with patch(
        "polyzymd.builders.solvent.SolventBuilder._get_box_shape_matrix",
        side_effect=AttributeError("no such shape"),
    ):
        result = CliRunner().invoke(cli, ["build", "-c", str(path), "--dry-run"])

    assert result.exit_code == 0, result.output
    assert "Box: not estimated (no such shape)" in result.output


def test_dry_run_box_report_says_a_substrate_is_left_out(tmp_path: Path) -> None:
    """The dry-run box comes from the enzyme PDB only; the line says so when a substrate is set."""
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    _write_elongated_pdb(tmp_path / "c" / "test.pdb")
    (tmp_path / "c" / "lig.sdf").write_text("\n")
    data = yaml.safe_load(path.read_text())
    data["substrate"] = {"name": "lig", "sdf_path": "lig.sdf"}
    path.write_text(yaml.safe_dump(data))

    result = CliRunner().invoke(cli, ["build", "-c", str(path), "--dry-run"])

    assert "from the enzyme PDB only; the substrate is not included" in result.output, result.output


def test_failed_build_keeps_a_run_folder_with_files(tmp_path: Path) -> None:
    """A failed build keeps the folder it made when files such as packmol_error.log are in it."""
    from polyzymd.config.loader import load_config
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    _write_elongated_pdb(tmp_path / "c" / "test.pdb")
    working_dir = load_config(path).get_working_directory(1)

    def write_log_and_fail(self, config, working_dir, polymer_seed, publish_topology=True):
        (Path(working_dir) / "packmol_error.log").write_text("Packmol failed\n")
        raise ValueError("Packmol failed; see packmol_error.log")

    with patch(
        "polyzymd.builders.system_builder.SystemBuilder.build_from_config", write_log_and_fail
    ):
        result = CliRunner().invoke(cli, ["build", "-c", str(path), "-r", "1"])

    assert result.exit_code == 1, result.output
    assert (working_dir / "packmol_error.log").is_file()


def test_build_that_cannot_take_the_lock_removes_nothing(tmp_path: Path) -> None:
    """A second build of a replicate that is being built leaves the first one's folder alone."""
    from polyzymd.config.loader import load_config
    from polyzymd.simulation.artifact_integrity import ArtifactIntegrityError
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    _write_elongated_pdb(tmp_path / "c" / "test.pdb")
    working_dir = load_config(path).get_working_directory(1)

    def first_build_makes_the_folder(working_dir):
        working_dir.mkdir(parents=True)
        raise ArtifactIntegrityError("Another PolyzyMD build or run holds the replicate lock")

    with patch(
        "polyzymd.simulation.artifact_integrity.replicate_lock", first_build_makes_the_folder
    ):
        result = CliRunner().invoke(cli, ["build", "-c", str(path), "-r", "1"])

    assert result.exit_code == 1, result.output
    assert working_dir.is_dir()


def test_failed_build_leaves_no_run_folder(tmp_path: Path) -> None:
    """A build that fails removes the run folder it made."""
    from polyzymd.config.loader import load_config
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    _write_elongated_pdb(tmp_path / "c" / "test.pdb")
    working_dir = load_config(path).get_working_directory(1)

    with patch(
        "polyzymd.builders.system_builder.SystemBuilder.build_from_config",
        side_effect=ValueError("atoms lie within 1.00 A of a periodic image"),
    ):
        result = CliRunner().invoke(cli, ["build", "-c", str(path), "-r", "1"])

    assert result.exit_code == 1, result.output
    assert "periodic image" in result.output
    assert not working_dir.exists()


def test_failed_gromacs_run_build_leaves_no_run_folder(tmp_path: Path) -> None:
    """`run --engine gromacs` removes the run folder its failed build made, and keeps an old one."""
    from polyzymd.cli.main import _run_gromacs_impl
    from polyzymd.config.loader import load_config
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    _write_elongated_pdb(tmp_path / "c" / "test.pdb")
    config = load_config(path)

    def make_folder_and_fail(self, config, working_dir, polymer_seed):
        Path(working_dir).mkdir(parents=True)
        raise ValueError("atoms lie within 1.00 A of a periodic image")

    with patch(
        "polyzymd.builders.system_builder.SystemBuilder.build_from_config", make_folder_and_fail
    ):
        with pytest.raises(ValueError):
            _run_gromacs_impl(config, replicate=1, gmx_path="gmx")
    assert not config.get_working_directory(1).exists()

    kept = config.get_working_directory(2)
    kept.mkdir(parents=True)
    with patch(
        "polyzymd.builders.system_builder.SystemBuilder.build_from_config",
        side_effect=ValueError("failed"),
    ):
        with pytest.raises(ValueError):
            _run_gromacs_impl(config, replicate=2, gmx_path="gmx")
    assert kept.is_dir()

    def write_log_and_fail(self, config, working_dir, polymer_seed):
        Path(working_dir).mkdir(parents=True)
        (Path(working_dir) / "packmol_error.log").write_text("Packmol failed\n")
        raise ValueError("Packmol failed; see packmol_error.log")

    with patch(
        "polyzymd.builders.system_builder.SystemBuilder.build_from_config", write_log_and_fail
    ):
        with pytest.raises(ValueError):
            _run_gromacs_impl(config, replicate=3, gmx_path="gmx")
    assert (config.get_working_directory(3) / "packmol_error.log").is_file()


@patch("polyzymd.engines.gromacs.engine.GromacsEngine.submit")
@patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
@patch("polyzymd.config.schema.SimulationConfig.from_yaml")
def test_gromacs_submit_checks_every_build_before_submitting(
    mock_from_yaml, _resolve, mock_engine_submit, tmp_path, monkeypatch
):
    """With replicate 2 unbuilt, replicate 1 is not submitted either."""
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

    def check_build(self, request):
        if request.replicate == 2:
            raise FileNotFoundError("No GROMACS build for replicate 2")

    monkeypatch.setattr("polyzymd.engines.gromacs.engine.GromacsEngine.check_build", check_build)
    config_path = tmp_path / "fake.yaml"
    config_path.write_text("name: test\n", encoding="utf-8")

    with patch("polyzymd.workflow.daisy_chain.check_existing_slurm_jobs", return_value=[]):
        with patch("polyzymd.workflow.daisy_chain.create_job_name", return_value="test_job"):
            result = CliRunner().invoke(
                cli, ["submit", "-c", str(config_path), "--engine", "gromacs", "-r", "1-2"]
            )

    assert result.exit_code == 1
    assert "No GROMACS build for replicate 2" in result.output
    mock_engine_submit.assert_not_called()
