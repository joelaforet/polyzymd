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
            water=SimpleNamespace(model=WaterModel.TIP3P),
            box=SimpleNamespace(padding=1.2),
            ions=SimpleNamespace(nacl_concentration=0.15, neutralize=True, ionic_strength=None),
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
        require_buildable=lambda engine=None, inputs=True: None,
    )


_ONE_ATOM_PDB = (
    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
)


def _config_with_pdb(tmp_path: Path) -> dict[str, object]:
    """Minimal config data whose enzyme PDB exists, with one atom."""
    (tmp_path / "enz.pdb").write_text(_ONE_ATOM_PDB)
    return _minimal_cli_config_data(tmp_path / "enz.pdb")


_SMALL_MOLECULE_FF = ["openff-2.3.0.offxml"]

#: NAGL charges with the model the config must name.
_NAGL = {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"}


def _solvent(**extra: object) -> dict[str, object]:
    """Solvent section with the small-molecule force field, TIP3P water and any ``extra`` keys."""
    return {"force_field": list(_SMALL_MOLECULE_FF), "water": {"model": "tip3p"}, **extra}


def _solute(type_: str, **keys: object) -> dict[str, object]:
    """A ``solutes`` entry of ``type_``; non-protein entries use the solvent force field.

    A protein takes its library (ff14SB) charges unless ``keys`` name others.
    """
    if type_ == "protein":
        keys = {"charges": {"method": "library"}, **keys}
    ff = ["ff14sb_off_impropers_0.0.4.offxml"] if type_ == "protein" else _SMALL_MOLECULE_FF
    return {"type": type_, **keys, "force_field": list(ff)}


def _minimal_cli_config_data(pdb_path: str | Path) -> dict[str, object]:
    """Create minimal YAML-serializable config data for CLI tests."""

    return {
        "name": "cli_reference_validation",
        "engine": "openmm",
        "solutes": [_solute("protein", name="Enz", pdb_path=str(pdb_path))],
        "solvent": _solvent(),
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

    @pytest.mark.parametrize("command", [["validate"], ["build", "--dry-run"]])
    def test_a_missing_or_empty_enzyme_pdb_is_an_error(self, tmp_path: Path, command) -> None:
        """validate and build --dry-run exit non-zero when the enzyme PDB is missing or empty."""
        config_path = tmp_path / "config.yaml"
        config_path.write_text(yaml.safe_dump(_minimal_cli_config_data("missing.pdb")))

        result = CliRunner().invoke(cli, [command[0], "-c", str(config_path), *command[1:]])

        assert result.exit_code == 1
        assert "missing.pdb does not exist" in result.output
        (tmp_path / "empty.pdb").write_text("END\n")
        config_path.write_text(yaml.safe_dump(_minimal_cli_config_data(tmp_path / "empty.pdb")))
        result = CliRunner().invoke(cli, [command[0], "-c", str(config_path), *command[1:]])
        assert result.exit_code == 1
        assert "has no ATOM or HETATM records" in result.output

    def test_validate_reports_derived_temperature_ramp_duration(self, tmp_path: Path) -> None:
        """Validation clearly reports rate-based heating duration."""
        data = _config_with_pdb(tmp_path)
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
        assert "Protocol warnings" not in result.output

    def test_validate_warns_when_equilibration_ends_below_production(self, tmp_path: Path) -> None:
        """A condition moved to 330 K whose stages still end at 300 K is valid but warned."""
        data = _config_with_pdb(tmp_path)
        data["thermodynamics"]["temperature"] = 330.0
        data["simulation_phases"]["equilibration_stages"] = [
            {
                "name": "heating",
                "ensemble": "NVT",
                "temperature_start": 60.0,
                "temperature_end": 300.0,
                "temperature_increment": 24.0,
                "temperature_interval_steps": 100,
            },
            {"name": "hold", "ensemble": "NVT", "temperature": 300.0, "duration": 0.002},
        ]
        config_path = tmp_path / "config.yaml"
        config_path.write_text(yaml.safe_dump(data), encoding="utf-8")

        result = CliRunner().invoke(cli, ["validate", "-c", str(config_path)])

        assert result.exit_code == 0
        assert "Configuration is valid!" in result.output
        assert "Protocol warnings:" in result.output
        assert (
            "equilibration stage 'hold' runs at 300 K, but production runs at "
            "thermodynamics.temperature = 330 K" in result.output
        )
        assert "simulation_phases.equilibration_stages.1.temperature" in result.output
        # A heating ramp that is not the last stage is not compared.
        assert "stage 'heating'" not in result.output
        assert "- hold: 300 K, 0.002 ns (NVT); WARNING: production runs at 330 K" in (result.output)

        # A copy moved to 330 K with the ramp still ending at 300 K jumps 30 K into 'hold'.
        data["simulation_phases"]["equilibration_stages"][1]["temperature"] = 330.0
        config_path.write_text(yaml.safe_dump(data), encoding="utf-8")
        result = CliRunner().invoke(cli, ["validate", "-c", str(config_path)])

        assert result.exit_code == 0
        assert (
            "equilibration stage 'heating' ends at 300 K, but the next stage 'hold' starts "
            "at 330 K" in result.output
        )
        assert "simulation_phases.equilibration_stages.0.temperature_end" in result.output
        assert "WARNING: the next stage starts at another temperature" in result.output

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
        data = _config_with_pdb(tmp_path)
        data["solvent"] = _solvent(
            water={"model": "tip3p"},
            co_solvents=[{"name": "dmso", "mole_fraction": 0.1, "charges": {"method": "preset"}}],
        )
        config_path = tmp_path / "config.yaml"
        config_path.write_text(yaml.safe_dump(data), encoding="utf-8")

        result = CliRunner().invoke(cli, ["validate", "-c", str(config_path)])

        assert result.exit_code == 0, result.output
        assert "Engine: openmm" in result.output
        assert "Co-solvents: dmso" in result.output


_DISTANCE_RESTRAINT = {
    "type": "flat_bottom",
    "name": "substrate_active_site",
    "atom1": {"component": "Enz", "selection": "resid 76 and name OG"},
    "atom2": {"component": "LIG", "selection": "name C1"},
    "distance": 3.5,
}


def _add_distance_restraint(data: dict) -> None:
    """Add the LIG substrate and a restraint from the enzyme to it."""
    data["solutes"].append(
        _solute(
            "small_molecule",
            name="LIG",
            sdf_path="missing.sdf",
            residue_name="LIG",
            charges=_NAGL,
        )
    )
    data["restraints"] = [_DISTANCE_RESTRAINT]


@pytest.mark.parametrize("command", [["validate"], ["build", "--dry-run"], ["run", "--dry-run"]])
def test_new_builds_refuse_a_water_model_the_build_cannot_make(
    tmp_path: Path, command: list[str]
) -> None:
    """validate, build and run refuse spce, which would build SPC/E charges on TIP3P water."""
    data = _minimal_cli_config_data("missing.pdb")
    data["solvent"] = _solvent(water={"model": "spce"})
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(data), encoding="utf-8")

    result = CliRunner().invoke(cli, [command[0], "-c", str(config_path), *command[1:]])

    assert result.exit_code != 0
    assert "model: tip3p for new systems" in result.output


def test_validate_refuses_distance_restraints_with_engine_gromacs(tmp_path: Path) -> None:
    """GROMACS runs apply no distance restraints, so validate refuses them."""
    data = _minimal_cli_config_data("missing.pdb")
    data["engine"] = "gromacs"
    _add_distance_restraint(data)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(data), encoding="utf-8")

    result = CliRunner().invoke(cli, ["validate", "-c", str(config_path)])

    assert result.exit_code != 0
    assert "set engine: openmm" in result.output


@pytest.mark.parametrize(
    "args",
    [
        ["build", "--format", "gromacs", "--dry-run"],
        ["run", "--engine", "gromacs", "--dry-run"],
        ["submit", "--engine", "gromacs", "--dry-run"],
    ],
)
def test_gromacs_overrides_refuse_distance_restraints(tmp_path: Path, args: list[str]) -> None:
    """A GROMACS build, run or submit of an OpenMM config never lists a restraint it drops."""
    data = _minimal_cli_config_data("missing.pdb")
    _add_distance_restraint(data)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(data), encoding="utf-8")

    result = CliRunner().invoke(cli, [args[0], "-c", str(config_path), *args[1:]])

    assert result.exit_code != 0
    assert "set engine: openmm" in result.output
    assert "ENABLED" not in result.output


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
        """Build dry-run warns about missing polymer files, which a build checks later."""
        data = _config_with_pdb(tmp_path)
        data["solutes"].append(
            _solute(
                "polymer",
                type_prefix="P",
                length=5,
                count=1,
                sdf_directory=str(tmp_path / "missing_sdfs"),
                monomers=[{"label": "A", "probability": 1.0}],
                charges={"method": "preset"},
            )
        )
        config_path = tmp_path / "config.yaml"
        config_path.write_text(yaml.safe_dump(data), encoding="utf-8")

        result = CliRunner().invoke(cli, ["build", "-c", str(config_path), "--dry-run"])

        assert result.exit_code == 0, result.output
        assert "Referenced file warnings" in result.output
        assert "Missing polymer SDF directory" in result.output

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
        self, from_config, _runner, tmp_path: Path, passing_audit
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

        manifest = json.loads((tmp_path / "build_provenance.json").read_text())
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

        validate.assert_called_once_with(tmp_path / "run_1", config, config_path="", replicate=1)
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
        self, _prefix, from_config, gromacs_runner, tmp_path: Path, capsys, monkeypatch
    ) -> None:
        from polyzymd.cli.main import _run_gromacs_impl

        monkeypatch.setattr(
            "polyzymd.simulation.artifact_integrity.require_audited_build", lambda *a, **k: {}
        )

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
            mdrun_flags={},  # the stub config has no gromacs section
        )
        assert f"Reusing the GROMACS files in {gromacs_dir}" in capsys.readouterr().out

    @patch("polyzymd.exporters.gromacs.GromacsRunner")
    @patch("polyzymd.builders.system_builder.SystemBuilder.from_config")
    @patch("polyzymd.analyses.shared.gromacs.system_prefix", return_value="sys")
    def test_gromacs_run_refuses_exported_files_no_audit_vouches_for(
        self, _prefix, from_config, gromacs_runner, tmp_path: Path
    ) -> None:
        """Exported files without a passing build record are not run."""
        from polyzymd.cli.main import _run_gromacs_impl
        from polyzymd.simulation.artifact_integrity import AuditGateError

        gromacs_dir = tmp_path / "run_1" / "gromacs"
        gromacs_dir.mkdir(parents=True)
        for name in ("sys.top", "sys.gro", "em.mdp", "eq_00_nvt.mdp", "prod.mdp"):
            (gromacs_dir / name).write_text("")

        with pytest.raises(AuditGateError, match="has no build record"):
            _run_gromacs_impl(self._config(tmp_path), replicate=1, gmx_path="gmx")
        gromacs_runner.assert_not_called()

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
        self, _prefix, gromacs_runner, tmp_path: Path, monkeypatch
    ) -> None:
        monkeypatch.setattr(
            "polyzymd.simulation.artifact_integrity.require_audited_build", lambda *a, **k: {}
        )
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

        result = CliRunner().invoke(cli, ["submit", "-c", str(config_path), "--time-limit", limit])

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
    (tmp_path / "c" / "test.pdb").write_text(_ONE_ATOM_PDB)
    data = yaml.safe_load(path.read_text())
    data["engine"] = "gromacs"
    data["solvent"] = _solvent(
        co_solvents=[
            {
                "name": "sds",
                "smiles": "CCCCCCCCCCCCOS(=O)(=O)[O-]",
                "count": 8,
                "charges": _NAGL,
            }
        ]
    )
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
    from rdkit import Chem

    (tmp_path / "c" / "lig.sdf").write_text(
        Chem.MolToMolBlock(Chem.MolFromSmiles("CCO")) + "$$$$\n"
    )
    data = yaml.safe_load(path.read_text())
    data["solutes"].append(_solute("small_molecule", name="lig", sdf_path="lig.sdf", charges=_NAGL))
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


@pytest.mark.parametrize(
    "command",
    [
        ["submit", "--dry-run", "--preset", "testing"],
        ["run", "--dry-run"],
        ["hash-trajectories", "--dry-run"],
        ["status"],
    ],
)
@pytest.mark.parametrize("text", ["- a\n- b\n", "name: x\nenzyme:\n\tpdb_path: a.pdb\n"])
def test_a_config_that_is_not_a_yaml_mapping_is_a_handled_error(tmp_path: Path, command, text):
    """A list or broken YAML gives an error line, not a traceback, in every command."""
    path = tmp_path / "config.yaml"
    path.write_text(text)
    result = CliRunner().invoke(cli, [command[0], "-c", str(path), *command[1:]])
    assert result.exit_code != 0
    assert isinstance(result.exception, SystemExit), result.exception
    assert "mapping" in result.output or "not valid YAML" in result.output


def test_validate_prints_error_and_fix_lines(tmp_path: Path) -> None:
    data = _config_with_pdb(tmp_path)
    data["simulation_phases"]["production"]["ensemble"] = "NVE"
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(data))
    result = CliRunner().invoke(cli, ["validate", "-c", str(path)])
    assert result.exit_code == 1
    assert "error: production: ensemble NVE is not run" in result.output
    assert "fix: use ensemble NVT or NPT." in result.output


@pytest.mark.parametrize("command", [["run", "--dry-run"], ["submit", "--dry-run"]])
def test_dry_runs_on_openmm_refuse_the_anisotropic_barostat(tmp_path: Path, command) -> None:
    """The dry runs make the same engine check as the run, so MCA on OpenMM fails early."""
    data = _config_with_pdb(tmp_path)
    data["engine"] = "gromacs"
    data["simulation_phases"]["production"]["barostat"] = "MCA"
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(data))
    result = CliRunner().invoke(
        cli, [command[0], "-c", str(path), *command[1:], "--engine", "openmm"]
    )
    assert result.exit_code != 0
    assert "MCA" in result.output


def test_submit_refuses_a_time_limit_slurm_cannot_read(tmp_path: Path) -> None:
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(_config_with_pdb(tmp_path)))
    result = CliRunner().invoke(
        cli, ["submit", "-c", str(path), "--dry-run", "--time-limit", "banana"]
    )
    assert result.exit_code == 2
    assert "is not a SLURM time" in result.output


def test_submit_does_not_check_the_inputs_of_a_build_it_does_not_make(tmp_path: Path) -> None:
    """submit runs existing builds, so a moved enzyme PDB does not stop it."""
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(_minimal_cli_config_data(tmp_path / "moved.pdb")))
    result = CliRunner().invoke(cli, ["submit", "-c", str(path), "--dry-run"])
    assert "moved.pdb does not exist" not in result.output, result.output


def _remove_thermodynamics(data: dict) -> None:
    del data["thermodynamics"]


@pytest.mark.parametrize("command", ["validate", "build", "run"])
@pytest.mark.parametrize(
    ("edit", "expected"),
    [
        (
            lambda data: data.update(pdb=data["solutes"][0].pop("pdb_path")),
            ["error: pdb: unknown key; did you mean solutes[].pdb_path?"],
        ),
        (
            _remove_thermodynamics,
            [
                "error: thermodynamics: required section missing; expected a mapping with "
                "temperature and pressure",
                "fix: start from the complete commented config that polyzymd study "
                "add-condition LABEL --new writes",
            ],
        ),
        (
            lambda data: data.update(simulation_phases="10 ns production"),
            [
                "error: simulation_phases: expected a mapping with minimization, "
                "equilibration_stages and production, not text"
            ],
        ),
    ],
    ids=["unknown key", "missing section", "wrong type"],
)
def test_config_errors_say_how_to_fix_them(tmp_path: Path, command, edit, expected) -> None:
    """validate, build and run print the same error and fix lines for a bad config."""
    data = _config_with_pdb(tmp_path)
    edit(data)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(data))

    result = CliRunner().invoke(cli, [command, "-c", str(config_path)])

    assert result.exit_code == 1, result.output
    for line in expected:
        assert line in result.stderr
    assert f"fix: correct these keys in {config_path}" in result.stderr
    assert "Configuration error:" not in result.output


class TestRunContinuesStartedProduction:
    """`polyzymd run` continues a started production instead of starting it again."""

    @staticmethod
    def _config(tmp_path: Path) -> SimpleNamespace:
        return TestRunReusesBuild._config(tmp_path)

    @staticmethod
    def _progress(steps_done: int, total: int = 2000):
        from polyzymd.simulation.progress import SegmentRecord, SegmentStatus, SimulationProgress

        status = SegmentStatus.COMPLETED if steps_done >= total else SegmentStatus.INTERRUPTED
        return SimulationProgress(
            total_steps_requested=total,
            segments=[SegmentRecord(index=0, status=status, steps_completed=steps_done)],
        )

    @patch("polyzymd.simulation.artifact_integrity.validate_build_bundle")
    @patch("polyzymd.cli.main._run_initial_segment")
    @patch("polyzymd.cli.main._run_segment_locked")
    def test_interrupted_production_continues_from_its_saved_state(
        self, run_segment, run_initial_segment, _validate, tmp_path: Path, capsys, monkeypatch
    ) -> None:
        from polyzymd.simulation import progress as progress_module

        monkeypatch.setattr(
            progress_module, "load_or_scan_progress", lambda **_: self._progress(540)
        )
        monkeypatch.setattr(progress_module, "load_progress", lambda _: self._progress(2000))

        _run_openmm_impl(self._config(tmp_path), replicate=1, config_path="/c/config.yaml")

        run_initial_segment.assert_not_called()
        run_segment.assert_called_once()
        assert run_segment.call_args.kwargs["skip_build"] is True
        assert "stopped at 540/2000 steps; continuing it" in capsys.readouterr().out

    @patch("polyzymd.cli.main._run_initial_segment")
    @patch("polyzymd.cli.main._run_segment_locked")
    def test_complete_production_is_left_alone(
        self, run_segment, run_initial_segment, tmp_path: Path, capsys, monkeypatch
    ) -> None:
        from polyzymd.simulation import progress as progress_module

        monkeypatch.setattr(
            progress_module, "load_or_scan_progress", lambda **_: self._progress(2000)
        )
        _run_openmm_impl(self._config(tmp_path), replicate=1)

        run_initial_segment.assert_not_called()
        run_segment.assert_not_called()
        assert "Replicate 1 is already complete" in capsys.readouterr().out

    @patch("polyzymd.simulation.artifact_integrity.validate_build_bundle")
    @patch("polyzymd.cli.main._run_segment_locked")
    def test_segment_without_progress_is_an_error(
        self, run_segment, _validate, tmp_path: Path, monkeypatch
    ) -> None:
        from polyzymd.simulation import progress as progress_module

        monkeypatch.setattr(
            progress_module, "load_or_scan_progress", lambda **_: self._progress(540)
        )
        monkeypatch.setattr(progress_module, "load_progress", lambda _: self._progress(540))
        with pytest.raises(RuntimeError, match="made no progress"):
            _run_openmm_impl(self._config(tmp_path), replicate=1)

    @patch("polyzymd.cli.main._run_openmm_impl")
    def test_interrupted_run_exits_99_and_names_the_command(self, run_impl, tmp_path) -> None:
        from polyzymd.simulation.signals import EXIT_CODE_INTERRUPTED, GracefulExit

        run_impl.side_effect = GracefulExit(15, 13500)
        config = Path(__file__).parents[2] / "examples/quickstart/config.yaml"
        result = CliRunner().invoke(cli, ["run", "-c", str(config), "-r", "2"])

        assert result.exit_code == EXIT_CODE_INTERRUPTED
        assert "Continue it with: polyzymd run -c" in result.output
        assert "-r 2" in result.output

    @patch("polyzymd.cli.main._run_segment_locked")
    def test_changed_config_is_refused_with_the_reason(
        self, run_segment, tmp_path: Path, monkeypatch
    ) -> None:
        from polyzymd.simulation import progress as progress_module
        from polyzymd.simulation.artifact_integrity import ArtifactIntegrityError

        def mismatch(*args, **kwargs):
            raise ArtifactIntegrityError("Configuration does not match build_provenance.json")

        monkeypatch.setattr(
            "polyzymd.simulation.artifact_integrity.validate_build_bundle", mismatch
        )
        monkeypatch.setattr("polyzymd.cli.main._build_matches_recorded_length", lambda *a: False)
        monkeypatch.setattr(
            progress_module, "load_or_scan_progress", lambda **_: self._progress(540)
        )
        (tmp_path / "run_1" / "production_0").mkdir(parents=True)

        with pytest.raises(ArtifactIntegrityError, match="Cannot continue replicate 1") as raised:
            _run_openmm_impl(self._config(tmp_path), replicate=1)
        assert "Configuration does not match" in str(raised.value)
        assert "study add-condition" in str(raised.value)
        run_segment.assert_not_called()

    @patch("polyzymd.cli.main._run_segment_locked")
    def test_changed_production_length_only_continues(
        self, run_segment, tmp_path: Path, monkeypatch, capsys
    ) -> None:
        from polyzymd.simulation import progress as progress_module
        from polyzymd.simulation.artifact_integrity import ArtifactIntegrityError

        def mismatch(*args, **kwargs):
            raise ArtifactIntegrityError("Configuration does not match")

        monkeypatch.setattr(
            "polyzymd.simulation.artifact_integrity.validate_build_bundle", mismatch
        )
        monkeypatch.setattr("polyzymd.cli.main._build_matches_recorded_length", lambda *a: True)
        monkeypatch.setattr(
            progress_module, "load_or_scan_progress", lambda **_: self._progress(540)
        )
        monkeypatch.setattr(progress_module, "load_progress", lambda _: self._progress(2000))
        (tmp_path / "run_1" / "production_0").mkdir(parents=True)

        _run_openmm_impl(self._config(tmp_path), replicate=1)
        run_segment.assert_called_once()
        assert "changes only the production length" in capsys.readouterr().out


class TestBuildMatchesRecordedLength:
    """A started chain may change its production length and nothing else."""

    @staticmethod
    def _built(tmp_path: Path):
        from polyzymd.config.loader import load_config
        from polyzymd.simulation.artifact_integrity import publish_build_bundle
        from polyzymd.simulation.progress import SimulationProgress, save_progress
        from tests.simulation.test_artifact_integrity import _tiny_openmm_bundle

        config = load_config(Path(__file__).parents[2] / "examples/quickstart/config.yaml")
        publish_build_bundle(tmp_path, *_tiny_openmm_bundle(), config)
        production = config.simulation_phases.production
        save_progress(
            tmp_path,
            SimulationProgress(
                total_steps_requested=int(production.duration * 1e6 / production.time_step),
                total_samples_requested=production.samples,
                timestep_fs=production.time_step,
            ),
        )
        return config

    def test_a_longer_production_matches(self, tmp_path: Path, passing_audit) -> None:
        from polyzymd.cli.main import _build_matches_recorded_length

        config = self._built(tmp_path)
        production = config.simulation_phases.production.model_copy(
            update={"duration": config.simulation_phases.production.duration * 2, "samples": 7}
        )
        longer = config.model_copy(
            update={
                "simulation_phases": config.simulation_phases.model_copy(
                    update={"production": production}
                )
            }
        )
        assert _build_matches_recorded_length(tmp_path, longer, "", 1)

    def test_the_first_segments_parameters_win_over_a_rewritten_progress(
        self, tmp_path: Path, passing_audit
    ) -> None:
        """A failed attempt with other samples rewrites progress.json; segment 0 still knows."""
        import json

        from polyzymd.cli.main import _build_matches_recorded_length
        from polyzymd.simulation.progress import SimulationProgress, save_progress

        config = self._built(tmp_path)
        production = config.simulation_phases.production
        params = tmp_path / "production_0" / "production_0_parameters.json"
        params.parent.mkdir()
        params.write_text(
            json.dumps(
                {
                    "__values__": {
                        "integ_params": {
                            "__values__": {
                                "total_time": {"__values__": {"value": production.duration}},
                                "num_samples": production.samples,
                            }
                        }
                    }
                }
            )
        )
        save_progress(
            tmp_path,
            SimulationProgress(
                total_steps_requested=123, total_samples_requested=1, timestep_fs=2.0
            ),
        )
        longer = config.model_copy(
            update={
                "simulation_phases": config.simulation_phases.model_copy(
                    update={
                        "production": production.model_copy(
                            update={"duration": production.duration * 2}
                        )
                    }
                )
            }
        )
        assert _build_matches_recorded_length(tmp_path, longer, "", 1)

    def test_an_extension_interrupted_without_segment_parameters_still_matches(
        self, tmp_path: Path, passing_audit
    ) -> None:
        """progress.json keeps the length production started with, whatever the totals say now."""
        from polyzymd.cli.main import _build_matches_recorded_length
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            load_or_scan_progress,
            load_progress,
            save_progress,
        )

        config = self._built(tmp_path)
        production = config.simulation_phases.production
        steps = int(production.duration * 1e6 / production.time_step)
        # The first run records the length production starts with.
        progress = load_or_scan_progress(
            tmp_path,
            total_steps=steps,
            total_samples=production.samples,
            timestep_fs=production.time_step,
        )
        progress.segments = [
            SegmentRecord(index=0, status=SegmentStatus.COMPLETED, steps_completed=steps),
            SegmentRecord(index=1, status=SegmentStatus.INTERRUPTED, steps_completed=10),
        ]
        save_progress(tmp_path, progress)
        # The extended run syncs the new totals, then is interrupted.
        extended = load_or_scan_progress(
            tmp_path,
            total_steps=steps * 2,
            total_samples=production.samples * 2,
            timestep_fs=production.time_step,
        )
        save_progress(tmp_path, extended)
        assert load_progress(tmp_path).started_total_steps == steps
        assert not (tmp_path / "production_0" / "production_0_parameters.json").exists()

        longer = config.model_copy(
            update={
                "simulation_phases": config.simulation_phases.model_copy(
                    update={
                        "production": production.model_copy(
                            update={
                                "duration": production.duration * 2,
                                "samples": production.samples * 2,
                            }
                        )
                    }
                )
            }
        )
        assert _build_matches_recorded_length(tmp_path, longer, "", 1)

    def test_the_old_length_written_into_the_folder_template_matches(
        self, tmp_path: Path, passing_audit
    ) -> None:
        """With {duration} in the folder name, the replicate keeps its folder this way."""
        from polyzymd.cli.main import _build_matches_recorded_length
        from polyzymd.config.loader import load_config
        from polyzymd.simulation.artifact_integrity import publish_build_bundle
        from polyzymd.simulation.progress import SimulationProgress, save_progress
        from tests.simulation.test_artifact_integrity import _tiny_openmm_bundle

        base = load_config(Path(__file__).parents[2] / "examples/quickstart/config.yaml")
        config = base.model_copy(
            update={
                "output": base.output.model_copy(
                    update={
                        "naming_template": "{enzyme}_{duration}ns_run{replicate}",
                        "scratch_directory": tmp_path,
                    }
                )
            }
        )
        working_dir = config.get_working_directory(1)
        publish_build_bundle(working_dir, *_tiny_openmm_bundle(), config)
        production = config.simulation_phases.production
        save_progress(
            working_dir,
            SimulationProgress(
                total_steps_requested=int(production.duration * 1e6 / production.time_step),
                total_samples_requested=production.samples,
                timestep_fs=production.time_step,
            ),
        )
        old = f"{production.duration:g}ns"
        longer = config.model_copy(
            update={
                "simulation_phases": config.simulation_phases.model_copy(
                    update={
                        "production": production.model_copy(
                            update={"duration": production.duration * 2}
                        )
                    }
                ),
                "output": config.output.model_copy(
                    update={"naming_template": "{enzyme}_" + old + "_run{replicate}"}
                ),
            }
        )
        assert longer.get_working_directory(1) == working_dir
        assert _build_matches_recorded_length(working_dir, longer, "", 1)

    def test_another_temperature_does_not(self, tmp_path: Path, passing_audit) -> None:
        from polyzymd.cli.main import _build_matches_recorded_length

        config = self._built(tmp_path)
        thermo = config.thermodynamics.model_copy(
            update={"temperature": config.thermodynamics.temperature + 10}
        )
        warmer = config.model_copy(update={"thermodynamics": thermo})
        assert not _build_matches_recorded_length(tmp_path, warmer, "", 1)


class TestGromacsRunStartedReplicate:
    """Local GROMACS `run` never repeats a started or finished production."""

    @staticmethod
    def _gromacs_dir(tmp_path: Path, *outputs: str) -> Path:
        gromacs_dir = tmp_path / "run_1" / "gromacs"
        gromacs_dir.mkdir(parents=True)
        for name in ("sys.top", "sys.gro", "em.mdp", "prod.mdp", *outputs):
            (gromacs_dir / name).write_text("")
        return gromacs_dir

    @staticmethod
    def _progress(monkeypatch, steps_done: int, total: int = 2000) -> list:
        from polyzymd.simulation.progress import SegmentRecord, SegmentStatus, SimulationProgress

        status = SegmentStatus.COMPLETED if steps_done >= total else SegmentStatus.INTERRUPTED
        progress = SimulationProgress(
            total_steps_requested=total,
            segments=[SegmentRecord(index=0, status=status, steps_completed=steps_done)],
        )
        gates: list = []
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.progress.load_or_scan_gromacs_progress",
            lambda **_: progress,
        )
        monkeypatch.setattr(
            "polyzymd.simulation.artifact_integrity.require_audited_build",
            lambda *a, **k: gates.append(a) or {},
        )
        return gates

    @patch("polyzymd.exporters.gromacs.GromacsRunner")
    @patch("polyzymd.analyses.shared.gromacs.system_prefix", return_value="sys")
    def test_complete_replicate_is_not_run_again(
        self, _prefix, gromacs_runner, tmp_path: Path, monkeypatch, capsys
    ) -> None:
        from polyzymd.cli.main import _run_gromacs_impl

        self._gromacs_dir(
            tmp_path, "em.gro", "prod.xtc", "prod.cpt", "prod_nojump.xtc", "prod_centered.xtc"
        )
        self._progress(monkeypatch, 2000)

        _run_gromacs_impl(TestRunReusesBuild._config(tmp_path), replicate=1, gmx_path="gmx")

        gromacs_runner.assert_not_called()
        assert "Replicate 1 is already complete" in capsys.readouterr().out

    @patch("polyzymd.exporters.gromacs.GromacsRunner")
    @patch("polyzymd.analyses.shared.gromacs.system_prefix", return_value="sys")
    def test_a_complete_replicate_missing_its_postprocessing_gets_it_again(
        self, _prefix, gromacs_runner, tmp_path: Path, monkeypatch, capsys
    ) -> None:
        """A stop during trjconv leaves no post-processed file; the next run writes them."""
        from polyzymd.cli.main import _run_gromacs_impl

        self._gromacs_dir(tmp_path, "em.gro", "prod.xtc", "prod.cpt")
        gates = self._progress(monkeypatch, 2000)
        monkeypatch.setattr("polyzymd.cli.main._save_gromacs_progress", lambda *a: None)

        _run_gromacs_impl(TestRunReusesBuild._config(tmp_path), replicate=1, gmx_path="gmx")

        gromacs_runner.return_value.run_postprocessing.assert_called_once()
        gromacs_runner.return_value.run_full_workflow.assert_not_called()
        assert gates == []
        assert "writing them again" in capsys.readouterr().out

    @patch("polyzymd.exporters.gromacs.GromacsRunner")
    @patch("polyzymd.analyses.shared.gromacs.system_prefix", return_value="sys")
    def test_started_production_continues_from_its_checkpoint_after_the_gates(
        self, _prefix, gromacs_runner, tmp_path: Path, monkeypatch
    ) -> None:
        from polyzymd.cli.main import _run_gromacs_impl

        self._gromacs_dir(tmp_path, "em.gro", "prod.tpr", "prod.cpt", "prod.xtc")
        gates = self._progress(monkeypatch, 500)
        monkeypatch.setattr("polyzymd.cli.main._save_gromacs_progress", lambda *a: None)

        _run_gromacs_impl(TestRunReusesBuild._config(tmp_path), replicate=1, gmx_path="gmx")

        assert len(gates) == 1
        gromacs_runner.return_value.continue_production.assert_called_once()
        gromacs_runner.return_value.run_full_workflow.assert_not_called()

    @patch("polyzymd.exporters.gromacs.GromacsRunner")
    @patch("polyzymd.analyses.shared.gromacs.system_prefix", return_value="sys")
    def test_started_production_without_a_checkpoint_says_how_to_start_it_again(
        self, _prefix, gromacs_runner, tmp_path: Path, monkeypatch
    ) -> None:
        from polyzymd.cli.main import _run_gromacs_impl
        from polyzymd.simulation.artifact_integrity import ArtifactIntegrityError

        gromacs_dir = self._gromacs_dir(tmp_path, "em.gro", "prod.tpr", "prod.xtc")
        (gromacs_dir / "prod.log").write_text(
            "   nsteps = 15000\n           Step           Time\n           3500        7.00000\n"
        )
        self._progress(monkeypatch, 0)

        with pytest.raises(ArtifactIntegrityError, match="wrote no checkpoint") as raised:
            _run_gromacs_impl(TestRunReusesBuild._config(tmp_path), replicate=1, gmx_path="gmx")
        assert "ran about 3500 of" in str(raised.value)
        assert "-cpi" not in str(raised.value)
        gromacs_runner.assert_not_called()

    @patch("polyzymd.exporters.gromacs.GromacsRunner")
    @patch("polyzymd.analyses.shared.gromacs.system_prefix", return_value="sys")
    def test_run_stopped_before_production_without_an_audit_record_is_refused(
        self, _prefix, gromacs_runner, tmp_path: Path, monkeypatch
    ) -> None:
        """A build from before 1.3.0rc5 has no record; it used to run again from minimization."""
        from polyzymd.cli.main import _run_gromacs_impl
        from polyzymd.simulation.artifact_integrity import AuditGateError
        from polyzymd.simulation.progress import SimulationProgress

        self._gromacs_dir(tmp_path, "em.log", "em.tpr")
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.progress.load_or_scan_gromacs_progress",
            lambda **_: SimulationProgress(total_steps_requested=2000),
        )
        with pytest.raises(AuditGateError, match="has no build record"):
            _run_gromacs_impl(TestRunReusesBuild._config(tmp_path), replicate=1, gmx_path="gmx")
        gromacs_runner.assert_not_called()

    @patch("polyzymd.exporters.gromacs.GromacsRunner")
    @patch("polyzymd.analyses.shared.gromacs.system_prefix", return_value="sys")
    def test_run_stopped_before_production_passes_the_audit_gate(
        self, _prefix, gromacs_runner, tmp_path: Path, monkeypatch
    ) -> None:
        """A run that stopped in minimization still checks the audited files and the config."""
        from polyzymd.cli.main import _run_gromacs_impl

        self._gromacs_dir(tmp_path, "em.log", "em.tpr")
        gates = self._progress(monkeypatch, 0)

        _run_gromacs_impl(TestRunReusesBuild._config(tmp_path), replicate=1, gmx_path="gmx")

        assert len(gates) == 1
        gromacs_runner.assert_called_once()


class TestRunWithoutBuildEnvironment:
    """In an environment that cannot build, `run` of an unbuilt replicate says where to build it."""

    @patch("polyzymd.cli.main._run_initial_segment")
    def test_unbuilt_replicate_names_the_build_command_and_writes_nothing(
        self, run_initial_segment, tmp_path: Path, monkeypatch
    ) -> None:
        from polyzymd.simulation.artifact_integrity import ArtifactIntegrityError

        monkeypatch.setattr("polyzymd.cli.main._missing_build_modules", lambda: ["rdkit"])
        config = TestRunReusesBuild._config(tmp_path)

        with pytest.raises(ArtifactIntegrityError) as raised:
            _run_openmm_impl(config, replicate=3, config_path="/p/config.yaml")

        message = str(raised.value)
        assert "Replicate 3 has not been built" in message
        assert "it has no rdkit" in message
        assert "pixi run -e build polyzymd build -c /p/config.yaml -r 3" in message
        assert not (tmp_path / "run_3").exists()
        run_initial_segment.assert_not_called()

    @patch("polyzymd.cli.main._run_initial_segment")
    def test_changed_build_says_restore_or_rebuild(
        self, run_initial_segment, tmp_path: Path, monkeypatch
    ) -> None:
        from polyzymd.simulation.artifact_integrity import ArtifactIntegrityError

        monkeypatch.setattr("polyzymd.cli.main._missing_build_modules", lambda: ["rdkit"])
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        (run_dir / "system.xml").write_text("<System/>")

        with pytest.raises(ArtifactIntegrityError, match="Restore the build's files"):
            _run_openmm_impl(TestRunReusesBuild._config(tmp_path), replicate=1)
        assert sorted(p.name for p in run_dir.iterdir()) == ["system.xml"]


class TestRunReportsOpenMMErrors:
    """An OpenMM platform error gets a hint and is recorded for `status`."""

    def test_ptx_error_names_the_sim_environment_and_is_recorded(self, tmp_path, monkeypatch):
        from polyzymd.cli.status_report import LOCAL_RUN_ERROR

        class OpenMMException(Exception):
            pass

        run_dir = tmp_path / "run_2"

        def fail(sim_config, replicate, config_path=""):
            run_dir.mkdir()
            raise OpenMMException(
                "Error loading CUDA module: CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)"
            )

        monkeypatch.setattr("polyzymd.cli.main._run_openmm_impl", fail)
        monkeypatch.setattr(
            "polyzymd.config.schema.SimulationConfig.get_working_directory",
            lambda self, replicate=1: tmp_path / f"run_{replicate}",
        )
        config = Path(__file__).parents[2] / "examples/quickstart/config.yaml"
        result = CliRunner().invoke(cli, ["run", "-c", str(config), "-r", "2"])

        assert result.exit_code == 1
        assert "OpenMM error: Error loading CUDA module" in result.output
        assert "sim-cuda" in result.output
        assert "UNSUPPORTED_PTX_VERSION" in (run_dir / LOCAL_RUN_ERROR).read_text()


class TestRunDryRunPlan:
    """`run --dry-run` says what the run would do, and fails when it could not."""

    def test_an_unbuilt_replicate_in_a_sim_environment_fails_the_dry_run(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.setattr("polyzymd.cli.main._missing_build_modules", lambda: ["openff.toolkit"])
        monkeypatch.setattr(
            "polyzymd.config.schema.SimulationConfig.get_working_directory",
            lambda self, replicate=1: tmp_path / f"run_{replicate}",
        )
        config = Path(__file__).parents[2] / "examples/quickstart/config.yaml"
        result = CliRunner().invoke(cli, ["run", "-c", str(config), "-r", "1", "--dry-run"])

        assert result.exit_code == 1
        assert "needs a build, which this environment cannot make" in result.output
        assert "run polyzymd build in the build environment first" in result.output
        assert not (tmp_path / "run_1").exists()

    def test_a_complete_replicate_is_reported_as_such(self, tmp_path, monkeypatch):
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            SimulationProgress,
            save_progress,
        )

        monkeypatch.setattr(
            "polyzymd.config.schema.SimulationConfig.get_working_directory",
            lambda self, replicate=1: tmp_path / f"run_{replicate}",
        )
        save_progress(
            tmp_path / "run_1",
            SimulationProgress(
                total_steps_requested=2000,
                segments=[
                    SegmentRecord(index=0, status=SegmentStatus.COMPLETED, steps_completed=2000)
                ],
            ),
        )
        config = Path(__file__).parents[2] / "examples/quickstart/config.yaml"
        result = CliRunner().invoke(cli, ["run", "-c", str(config), "-r", "1", "--dry-run"])
        assert result.exit_code == 0, result.output
        assert "already complete; nothing to run" in result.output

    def test_a_raised_production_length_is_previewed_as_a_continuation(self, tmp_path, monkeypatch):
        """After duration is raised, the preview says from which step to which."""
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            SimulationProgress,
            save_progress,
        )

        monkeypatch.setattr(
            "polyzymd.config.schema.SimulationConfig.get_working_directory",
            lambda self, replicate=1: tmp_path / f"run_{replicate}",
        )
        monkeypatch.setattr(
            "polyzymd.simulation.artifact_integrity.validate_build_bundle", lambda *a, **k: None
        )
        save_progress(
            tmp_path / "run_1",
            SimulationProgress(
                total_steps_requested=1000,
                segments=[
                    SegmentRecord(
                        index=0,
                        status=SegmentStatus.COMPLETED,
                        steps_completed=1000,
                        report_interval=500,
                    )
                ],
            ),
        )
        (tmp_path / "run_1" / "production_0").mkdir()
        (tmp_path / "run_1" / "production_0" / "production_0_state.xml").write_text("<State/>")
        config = Path(__file__).parents[2] / "examples/quickstart/config.yaml"
        result = CliRunner().invoke(cli, ["run", "-c", str(config), "-r", "1", "--dry-run"])
        assert result.exit_code == 0, result.output
        assert "continue production from its saved state, from step 1000 to 2000" in result.output


class TestRefusedRunWritesNothing:
    """A refused `polyzymd run` leaves progress.json as it was."""

    @patch("polyzymd.simulation.artifact_integrity.validate_build_bundle")
    @patch("polyzymd.cli.main._run_segment_locked")
    def test_a_frame_spacing_change_is_refused_before_any_write(
        self, run_segment, _validate, tmp_path: Path
    ) -> None:
        from polyzymd.simulation.artifact_integrity import ArtifactIntegrityError
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            SimulationProgress,
            save_progress,
        )

        run_dir = tmp_path / "run_1"
        save_progress(
            run_dir,
            SimulationProgress(
                total_steps_requested=2000,
                total_samples_requested=4,
                timestep_fs=2.0,
                segments=[
                    SegmentRecord(
                        index=0,
                        status=SegmentStatus.COMPLETED,
                        steps_completed=2000,
                        report_interval=500,
                    )
                ],
            ),
        )
        (run_dir / "production_0").mkdir()
        (run_dir / "production_0" / "production_0_state.xml").write_text("<State/>")
        before = (run_dir / "progress.json").read_text()
        production = SimpleNamespace(
            duration=0.008, samples=4, time_step=2.0, checkpoint_interval=60.0
        )
        config = SimpleNamespace(
            simulation_phases=SimpleNamespace(production=production),
            require_engine_barostats=lambda engine: None,
            get_working_directory=lambda replicate: tmp_path / f"run_{replicate}",
        )

        with pytest.raises(ArtifactIntegrityError, match="Cannot continue replicate 1") as raised:
            _run_openmm_impl(config, replicate=1)

        assert "samples: 8" in str(raised.value)
        assert "run-segment" not in str(raised.value)
        assert (run_dir / "progress.json").read_text() == before
        run_segment.assert_not_called()


class TestGromacsStopMessages:
    def test_each_stage_says_what_a_stop_left(self, tmp_path):
        from polyzymd.cli.main import _gromacs_stop_outcome

        assert (
            "no checkpoint"
            in _gromacs_stop_outcome("Running production MD (mdrun)", tmp_path).lower()
        )
        (tmp_path / "prod.cpt").write_text("")
        assert "latest checkpoint" in _gromacs_stop_outcome(
            "Running production MD (mdrun)", tmp_path
        )
        assert "post-processed" in _gromacs_stop_outcome("Removing PBC jumps (trjconv)", tmp_path)
        assert "starts them again" in _gromacs_stop_outcome(
            "Running energy minimization (mdrun)", tmp_path
        )
