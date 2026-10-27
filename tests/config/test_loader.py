"""Tests for config loader path handling and YAML dumping."""

from pathlib import Path

import pytest


class TestSaveConfigYamlDumper:
    """save_config uses a local Dumper subclass, not the global yaml.Dumper."""

    def test_global_dumper_unmodified_after_save(self, tmp_path):
        """Calling save_config must not add representers to yaml.Dumper."""
        import yaml

        from polyzymd.config.loader import save_config
        from polyzymd.config.schema import SimulationConfig

        original_representer = yaml.Dumper.yaml_representers.get(str)

        config = SimulationConfig.model_construct(
            name="test",
            _fields_set=set(),
        )
        out = tmp_path / "test_config.yaml"
        save_config(config, out)

        after_representer = yaml.Dumper.yaml_representers.get(str)
        message = "save_config must not mutate yaml.Dumper.yaml_representers"
        assert after_representer is original_representer, message

    def test_multiline_strings_use_block_style(self, tmp_path):
        """Multiline strings in saved YAML use literal block (|) style."""
        from polyzymd.config.loader import save_config
        from polyzymd.config.schema import SimulationConfig

        config = SimulationConfig.model_construct(
            name="multi\nline\ntest",
            _fields_set=set(),
        )
        out = tmp_path / "test_config.yaml"
        save_config(config, out)

        content = out.read_text()
        assert "|-" in content or "|\n" in content


class TestPolymerPathResolution:
    """_expand_paths and _convert_paths_to_relative handle the polymer directories."""

    @pytest.mark.parametrize("key", ["sdf_directory", "cache_directory"])
    def test_polymer_directories_are_expanded_and_relativized(self, key):
        """Relative polymer directories are resolved against the config and saved relative."""
        from polyzymd.config.loader import _convert_paths_to_relative, _expand_paths

        expanded = _expand_paths({"solutes": [{key: "polymers"}]}, Path("/configs"))
        assert expanded["solutes"][0][key] == "/configs/polymers"
        saved = _convert_paths_to_relative(expanded, Path("/configs"))
        assert saved["solutes"][0][key] == "polymers"

    def test_removed_reaction_keys_are_not_paths(self):
        """The old reactions block is left as written for the schema to refuse."""
        from polyzymd.config.loader import PATH_KEYS, _expand_paths

        assert not {"initiation", "polymerization", "termination"} & PATH_KEYS
        data = {"reactions": {"initiation": "templates/init.rxn"}}
        assert _expand_paths(data, Path("/configs")) == data


def test_output_folders_are_relative_to_the_config(tmp_path: Path) -> None:
    """projects_directory '.' is the config's folder, from any shell folder."""
    import os

    import yaml

    from polyzymd.config.schema import SimulationConfig
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "water", scratch=Path("."))
    data = yaml.safe_load(path.read_text())
    data["output"]["projects_directory"] = "."
    path.write_text(yaml.safe_dump(data))
    old = Path.cwd()
    os.chdir(tmp_path)
    try:
        config = SimulationConfig.from_yaml(path)
    finally:
        os.chdir(old)
    assert Path(config.output.projects_directory) == (tmp_path / "water").resolve()


def test_an_output_folder_with_dots_is_kept_as_written(tmp_path: Path) -> None:
    """projects_directory '../sims' keeps its '..', so the config hash of existing builds holds."""
    import yaml

    from polyzymd.config.schema import SimulationConfig
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "water", scratch=Path("."))
    data = yaml.safe_load(path.read_text())
    data["output"]["projects_directory"] = "../sims"
    data["output"]["scratch_directory"] = "../sims"
    path.write_text(yaml.safe_dump(data))

    config = SimulationConfig.from_yaml(path)

    assert str(config.output.projects_directory) == str(tmp_path / "water" / ".." / "sims")
    assert str(config.output.scratch_directory) == str(tmp_path / "water" / ".." / "sims")


def test_configs_of_existing_runs_load_with_values_new_builds_refuse(tmp_path: Path) -> None:
    """Old spce, KCl and GROMACS-restraint configs load, so status and analysis still work."""
    import yaml

    from polyzymd.config.loader import load_config
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "water", scratch=tmp_path / "s")
    data = yaml.safe_load(path.read_text())
    data["engine"] = "gromacs"
    data["solvent"]["water"] = {"model": "spce"}
    data["solvent"]["ions"] = {"kcl_concentration": 0.15}
    data["restraints"] = [
        {
            "type": "flat_bottom",
            "name": "site",
            "atom1": {"component": "TestEnzyme", "selection": "resid 1 and name CA"},
            "atom2": {"component": "TestEnzyme", "selection": "resid 2 and name CA"},
            "distance": 3.5,
        }
    ]
    path.write_text(yaml.safe_dump(data))

    config = load_config(path)

    assert config.solvent.water.model.value == "spce"
    assert config.solvent.ions.kcl_concentration == 0.15


def _config_file(tmp_path: Path, edit=None) -> Path:
    """Write the minimal test config, changed by ``edit(data)``, and return its path."""
    import yaml

    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    if edit is not None:
        data = yaml.safe_load(path.read_text())
        edit(data)
        path.write_text(yaml.safe_dump(data))
    return path


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("- a\n- b\n", "must be a mapping of config sections, not a list"),
        ("42\n", "must be a mapping of config sections, not a int"),
        ("name: x\nenzyme:\n\tpdb_path: a.pdb\n", "not valid YAML"),
        ("name: [x\n", "not valid YAML"),
        ("name: x\nengine: openmm\nengine: gromacs\n", "duplicate key 'engine'"),
    ],
)
def test_a_file_that_is_not_a_config_mapping_is_refused(tmp_path: Path, text, message) -> None:
    """Not-a-mapping, broken YAML and keys given twice are refused with a fix line."""
    from polyzymd.config.loader import ConfigFileError, load_config

    path = tmp_path / "config.yaml"
    path.write_text(text)
    with pytest.raises(ConfigFileError, match=message) as error:
        load_config(path)
    assert "fix:" in str(error.value)
    if "YAML" in message:
        assert "line" in str(error.value)


def test_home_and_variables_expand_in_paths(tmp_path: Path, monkeypatch) -> None:
    """'~' is the home folder; an unset $VAR is kept as written and a new build names it."""
    from polyzymd.config.loader import load_config

    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("POLYZYMD_TEST_UNSET", raising=False)

    def edit(data):
        data["output"]["scratch_directory"] = "~/sims"
        data["output"]["projects_directory"] = "$POLYZYMD_TEST_UNSET/projects"

    config = load_config(_config_file(tmp_path, edit))
    assert config.output.scratch_directory == tmp_path / "home" / "sims"
    assert str(config.output.projects_directory) == "$POLYZYMD_TEST_UNSET/projects"
    with pytest.raises(ValueError, match=r"\$POLYZYMD_TEST_UNSET, which is not set"):
        config.require_buildable()


def test_yaml_merge_keys_load_and_may_be_overridden(tmp_path: Path) -> None:
    """An anchor merged with '<<:' loads, and a key the mapping sets again overrides it."""
    from polyzymd.config.loader import read_yaml_mapping

    path = tmp_path / "config.yaml"
    path.write_text(
        "base: &base {ensemble: NVT, duration: 0.1}\n" "stage: {<<: *base, duration: 0.2}\n"
    )
    assert read_yaml_mapping(path)["stage"] == {"ensemble": "NVT", "duration": 0.2}


def test_runs_written_under_a_literal_tilde_folder_are_still_found(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    """Earlier versions wrote '~/sims' runs to <config dir>/~/sims; that folder is kept."""
    from polyzymd.config.loader import load_config

    monkeypatch.setenv("HOME", str(tmp_path / "home"))

    def edit(data):
        data["output"]["scratch_directory"] = "~/sims"

    path = _config_file(tmp_path, edit)
    (tmp_path / "c" / "~" / "sims").mkdir(parents=True)
    assert load_config(path).output.scratch_directory == tmp_path / "c" / "~" / "sims"
    assert "earlier versions put the runs" in caplog.text
    (tmp_path / "home" / "sims").mkdir(parents=True)
    assert load_config(path).output.scratch_directory == tmp_path / "home" / "sims"


def _error_lines(tmp_path: Path, edit) -> list[tuple[str, str]]:
    """Return the error and fix lines of the minimal config changed by ``edit``."""
    from pydantic import ValidationError

    from polyzymd.config.loader import config_error_lines, load_config

    path = _config_file(tmp_path, edit)
    with pytest.raises(ValidationError) as error:
        load_config(path)
    return config_error_lines(error.value, path)


def test_an_unknown_key_names_the_closest_field_in_any_section(tmp_path: Path) -> None:
    def edit(data):
        data["pdb"] = data["solutes"][0].pop("pdb_path")
        data["system"] = {"pdb": "x.pdb", "water": "tip3p"}
        data["thermodynamics"]["temprature"] = 300.0

    lines = _error_lines(tmp_path, edit)

    assert ("error", "pdb: unknown key; did you mean solutes[].pdb_path?") in lines
    assert ("error", "system.pdb: unknown key; did you mean solutes[].pdb_path?") in lines
    assert ("error", "system.water: unknown key; did you mean solvent.water?") in lines
    assert (
        "error",
        "thermodynamics.temprature: unknown key; did you mean thermodynamics.temperature?",
    ) in lines


def test_an_unknown_key_without_a_close_field_lists_the_keys_of_its_section(
    tmp_path: Path,
) -> None:
    lines = _error_lines(tmp_path, lambda data: data["thermodynamics"].update(zzz=1))

    assert ("error", "thermodynamics.zzz: unknown key; valid keys here: temperature, pressure") in (
        lines
    )


def test_a_missing_section_points_to_the_commented_template(tmp_path: Path) -> None:
    from polyzymd.config.loader import TEMPLATE_FIX

    def edit(data):
        del data["thermodynamics"]
        del data["simulation_phases"]

    lines = _error_lines(tmp_path, edit)

    assert (
        "error",
        "thermodynamics: required section missing; expected a mapping with temperature and "
        "pressure",
    ) in lines
    assert [line for line in lines if line == ("fix", TEMPLATE_FIX)] == [("fix", TEMPLATE_FIX)]
    assert "polyzymd study add-condition LABEL --new" in TEMPLATE_FIX


def test_a_value_of_the_wrong_type_names_the_expected_shape(tmp_path: Path) -> None:
    lines = _error_lines(tmp_path, lambda data: data.update(simulation_phases="10 ns production"))

    assert (
        "error",
        "simulation_phases: expected a mapping with minimization, equilibration_stages and "
        "production, not text",
    ) in lines
    assert lines[-1][0] == "fix"


def test_a_misspelled_section_names_the_section(tmp_path: Path) -> None:
    from polyzymd.config.loader import TEMPLATE_FIX

    lines = _error_lines(tmp_path, lambda data: data.update(solute=data.pop("solutes")))

    assert ("error", "solute: unknown key; did you mean solutes?") in lines
    assert ("fix", TEMPLATE_FIX) not in lines


def test_a_mapping_given_for_a_list_names_the_list(tmp_path: Path) -> None:
    def edit(data):
        phases = data["simulation_phases"]
        phases["equilibration_stages"] = phases["equilibration_stages"][0]

    lines = _error_lines(tmp_path, edit)

    assert any(
        text.startswith("simulation_phases.equilibration_stages: expected a list of mappings with")
        and text.endswith("not a mapping")
        for _, text in lines
    ), lines


def test_save_config_leaves_out_keys_no_code_reads(tmp_path: Path) -> None:
    """Unread keys at their default are not saved; any other value is kept, so the hash holds."""
    import yaml

    from polyzymd.config.loader import load_config, save_config
    from polyzymd.simulation.artifact_integrity import config_hash

    path = _config_file(
        tmp_path,
        lambda data: data.setdefault("output", {}).update(
            trajectory_format="xtc", save_checkpoint=True
        ),
    )
    out = tmp_path / "saved.yaml"
    loaded = load_config(path)
    save_config(loaded, out)

    output = yaml.safe_load(out.read_text())["output"]
    assert not {"save_checkpoint", "save_state_data"} & set(output)
    assert output["trajectory_format"] == "xtc"
    assert config_hash(load_config(out)) == config_hash(loaded)
