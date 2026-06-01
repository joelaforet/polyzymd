"""Tests for the legacy conversion utility."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest
import yaml

from polyzymd.config.comparison import ComparisonConfig
from polyzymd.config.schema import SimulationConfig

SCRIPT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "convert_legacy.py"


@pytest.fixture(scope="module")
def convert_legacy() -> ModuleType:
    """Load ``scripts/convert_legacy.py`` as an importable module.

    Returns
    -------
    ModuleType
        Imported conversion script module.
    """
    spec = importlib.util.spec_from_file_location("convert_legacy", SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load script module from {SCRIPT_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_parse_folder_name_with_polymer(convert_legacy: ModuleType) -> None:
    """Folder parsing extracts replicate, polymer, and phase metadata."""
    folder = (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_SBMA-EGMA-50%_38x_"
        "363.0K_0.5ns-NVT_1000.0ns-NPT_run2"
    )

    metadata = convert_legacy.parse_folder_name(folder)

    assert metadata.enzyme_name == "LipA"
    assert metadata.replicate == 2
    assert metadata.temperature_K == 363.0
    assert metadata.production.duration_ns == 1000.0
    assert metadata.polymer is not None
    assert metadata.polymer.type_prefix == "SBMA-EGMA"
    assert metadata.polymer.chain_count == 38


def test_discover_files_uses_exact_segment_layout(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Daisy-chain discovery accepts exact segment topology and trajectories."""
    sim_dir = tmp_path / "legacy"
    prod0 = sim_dir / "production_0"
    prod1 = sim_dir / "production_1"
    prod0.mkdir(parents=True)
    prod1.mkdir()
    (prod0 / "production_0_topology.pdb").write_text("ATOM\n", encoding="utf-8")
    (prod0 / "production_0_trajectory.dcd").write_bytes(b"DCD")
    (prod1 / "production_1_trajectory.dcd").write_bytes(b"DCD")

    metadata = convert_legacy.SimMetadata(
        folder_name="legacy",
        restraint_distance="10A",
        enzyme_name="LipA",
        substrate_name="Resorufin-Butyrate",
        temperature_K=363.0,
        replicate=1,
    )

    convert_legacy.discover_files(sim_dir, metadata)

    assert metadata.topology_path == prod0 / "production_0_topology.pdb"
    assert metadata.trajectory_paths == [
        prod0 / "production_0_trajectory.dcd",
        prod1 / "production_1_trajectory.dcd",
    ]


def test_discover_files_rejects_incomplete_daisy_chain_segment(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Daisy-chain discovery fails when an expected segment trajectory is missing."""
    sim_dir = tmp_path / "legacy"
    prod0 = sim_dir / "production_0"
    prod1 = sim_dir / "production_1"
    prod0.mkdir(parents=True)
    prod1.mkdir()
    (prod0 / "production_0_topology.pdb").write_text("ATOM\n", encoding="utf-8")
    (prod0 / "production_0_trajectory.dcd").write_bytes(b"DCD")

    metadata = convert_legacy.SimMetadata(
        folder_name="legacy",
        restraint_distance="10A",
        enzyme_name="LipA",
        substrate_name="Resorufin-Butyrate",
        temperature_K=363.0,
        replicate=1,
    )

    with pytest.raises(FileNotFoundError, match="production_1_trajectory.dcd"):
        convert_legacy.discover_files(sim_dir, metadata)

    assert metadata.trajectory_paths == []


def test_discover_files_rejects_arbitrary_names(convert_legacy: ModuleType, tmp_path: Path) -> None:
    """Discovery rejects non-exact topology and trajectory names."""
    sim_dir = tmp_path / "legacy"
    prod0 = sim_dir / "production_0"
    prod0.mkdir(parents=True)
    (prod0 / "custom_topology.pdb").write_text("ATOM\n", encoding="utf-8")
    (prod0 / "production_custom_trajectory.dcd").write_bytes(b"DCD")

    metadata = convert_legacy.SimMetadata(
        folder_name="legacy",
        restraint_distance="10A",
        enzyme_name="LipA",
        substrate_name="Resorufin-Butyrate",
        temperature_K=363.0,
        replicate=1,
    )

    with pytest.raises(FileNotFoundError, match="No topology PDB found"):
        convert_legacy.discover_files(sim_dir, metadata)


def test_generate_config_yaml_includes_engine_and_loads(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Generated simulation config should include canonical OpenMM engine."""
    metadata = convert_legacy.SimMetadata(
        folder_name="10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
        restraint_distance="10A",
        enzyme_name="LipA",
        substrate_name="Resorufin-Butyrate",
        temperature_K=363.0,
        replicate=1,
        equilibration=convert_legacy.PhaseParams(
            ensemble="NVT",
            duration_ns=0.5,
            samples=10,
            time_step_fs=2.0,
            temperature_K=363.0,
            pressure_atm=1.0,
        ),
        production=convert_legacy.PhaseParams(
            ensemble="NPT",
            duration_ns=1000.0,
            samples=2500,
            time_step_fs=2.0,
            temperature_K=363.0,
            pressure_atm=1.0,
        ),
    )
    config_path = tmp_path / "config.yaml"

    convert_legacy.generate_config_yaml(metadata, tmp_path, config_path)

    data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    assert data["engine"] == "openmm"
    config = SimulationConfig.from_yaml(config_path)
    assert config.engine == "openmm"


def test_generate_comparison_yaml_uses_canonical_plot_key(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Generated comparison config should use catalytic_triad plot settings."""
    converted = tmp_path / "converted"
    sim_dir = converted / (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1"
    )
    sim_dir.mkdir(parents=True)
    (sim_dir / "config.yaml").write_text("name: legacy\n", encoding="utf-8")

    comparison_yaml = convert_legacy.generate_comparison_yaml(
        output_dir=converted,
        comparison_dir=tmp_path / "comparison",
        project_name="legacy_compare",
        control_label="No Polymer (Control)",
    )

    data = yaml.safe_load(comparison_yaml.read_text(encoding="utf-8"))
    assert "catalytic_triad" in data["plot_settings"]
    assert "triad" not in data["plot_settings"]
    config = ComparisonConfig.from_yaml(comparison_yaml)
    assert config.validate_config() == []


def _write_converted_config(root: Path, folder_name: str) -> Path:
    """Create a synthetic converted simulation directory with config.yaml.

    Parameters
    ----------
    root : Path
        Converted-output root directory.
    folder_name : str
        Legacy-format simulation folder name.

    Returns
    -------
    Path
        Path to the generated config file.
    """
    sim_dir = root / folder_name
    sim_dir.mkdir(parents=True)
    config_path = sim_dir / "config.yaml"
    config_path.write_text("name: legacy\nengine: openmm\n", encoding="utf-8")
    return config_path


def test_group_conditions_separates_duplicate_control_labels(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Control-like conditions with different restraint or conformation stay distinct."""
    converted = tmp_path / "converted"
    config_a = _write_converted_config(
        converted,
        "10A_RESTRAINT_CALB_Resorufin-Butyrate_conf1_363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )
    config_b = _write_converted_config(
        converted,
        "12A_RESTRAINT_CALB_Resorufin-Butyrate_conf2_363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )
    config_b_rep2 = _write_converted_config(
        converted,
        "12A_RESTRAINT_CALB_Resorufin-Butyrate_conf2_363.0K_0.5ns-NVT_1000.0ns-NPT_run2",
    )

    conditions = convert_legacy.group_conditions(converted)

    assert len(conditions) == 2
    assert all(label.startswith("No Polymer (Control) (") for label in conditions)
    assert any("10A" in label and "conf1" in label for label in conditions)
    assert any("12A" in label and "conf2" in label for label in conditions)
    assert sorted(cond["config"] for cond in conditions.values()) == sorted([config_a, config_b])
    assert sorted(cond["replicates"] for cond in conditions.values()) == [[1], [1, 2]]
    assert config_b_rep2 not in [cond["config"] for cond in conditions.values()]


def test_group_conditions_separates_duplicate_polymer_chain_counts(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Polymer conditions with identical friendly labels but chain counts stay distinct."""
    converted = tmp_path / "converted"
    config_38 = _write_converted_config(
        converted,
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_SBMA-EGMA-50%_38x_"
        "363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )
    config_77 = _write_converted_config(
        converted,
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_SBMA-EGMA-50%_77x_"
        "363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )
    _write_converted_config(
        converted,
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_SBMA-EGMA-50%_77x_"
        "363.0K_0.5ns-NVT_1000.0ns-NPT_run2",
    )

    conditions = convert_legacy.group_conditions(converted)

    assert len(conditions) == 2
    assert all(label.startswith("SBMA-EGMA 50% (") for label in conditions)
    assert any("38 chains" in label for label in conditions)
    assert any("77 chains" in label for label in conditions)
    assert sorted(cond["config"] for cond in conditions.values()) == sorted([config_38, config_77])
    assert sorted(cond["replicates"] for cond in conditions.values()) == [[1], [1, 2]]


def test_generate_comparison_yaml_preserves_disambiguated_entries(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Comparison YAML should include all disambiguated duplicate-label conditions."""
    converted = tmp_path / "converted"
    control_config = _write_converted_config(
        converted,
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )
    polymer_38 = _write_converted_config(
        converted,
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_SBMA-EGMA-50%_38x_"
        "363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )
    polymer_77 = _write_converted_config(
        converted,
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_SBMA-EGMA-50%_77x_"
        "363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )

    comparison_yaml = convert_legacy.generate_comparison_yaml(
        output_dir=converted,
        comparison_dir=tmp_path / "comparison",
        project_name="legacy_compare",
        control_label="No Polymer (Control)",
    )

    data = yaml.safe_load(comparison_yaml.read_text(encoding="utf-8"))
    entries = {entry["label"]: entry for entry in data["conditions"]}
    assert list(entries)[0] == "No Polymer (Control)"
    assert len(entries) == 3
    assert entries["No Polymer (Control)"]["replicates"] == [1]
    assert any(label.startswith("SBMA-EGMA 50% (") for label in entries)
    assert any("38 chains" in label for label in entries)
    assert any("77 chains" in label for label in entries)
    config_values = {entry["config"] for entry in entries.values()}
    assert f"../{control_config.relative_to(tmp_path)}" in config_values
    assert f"../{polymer_38.relative_to(tmp_path)}" in config_values
    assert f"../{polymer_77.relative_to(tmp_path)}" in config_values
    config = ComparisonConfig.from_yaml(comparison_yaml)
    assert config.validate_config() == []


def test_generate_comparison_yaml_rejects_ambiguous_control_base_label(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Ambiguous base control labels should ask for an exact disambiguated label."""
    converted = tmp_path / "converted"
    _write_converted_config(
        converted,
        "10A_RESTRAINT_CALB_Resorufin-Butyrate_conf1_363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )
    _write_converted_config(
        converted,
        "12A_RESTRAINT_CALB_Resorufin-Butyrate_conf2_363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )

    with pytest.raises(ValueError, match="--control") as exc_info:
        convert_legacy.generate_comparison_yaml(
            output_dir=converted,
            comparison_dir=tmp_path / "comparison",
            project_name="legacy_compare",
            control_label="No Polymer (Control)",
        )

    message = str(exc_info.value)
    assert "No Polymer (Control) (10A, CALB" in message
    assert "No Polymer (Control) (12A, CALB" in message
    assert 'Pass --control "<exact label>"' in message


def test_generate_comparison_yaml_accepts_explicit_disambiguated_control(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """An exact disambiguated control label should be accepted and written."""
    converted = tmp_path / "converted"
    _write_converted_config(
        converted,
        "10A_RESTRAINT_CALB_Resorufin-Butyrate_conf1_363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )
    _write_converted_config(
        converted,
        "12A_RESTRAINT_CALB_Resorufin-Butyrate_conf2_363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )
    control_label = (
        "No Polymer (Control) " "(12A, CALB, Resorufin-Butyrate, conf2, 363K, no polymer)"
    )

    comparison_yaml = convert_legacy.generate_comparison_yaml(
        output_dir=converted,
        comparison_dir=tmp_path / "comparison",
        project_name="legacy_compare",
        control_label=control_label,
    )

    data = yaml.safe_load(comparison_yaml.read_text(encoding="utf-8"))
    assert data["control"] == control_label
    assert data["conditions"][0]["label"] == control_label


def test_validate_output_config_validation_failure_returns_false(
    convert_legacy: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Invalid generated simulation config should fail output validation."""

    class FakeAtomGroup:
        def __init__(self, size: int) -> None:
            self._size = size

        def __len__(self) -> int:
            return self._size

    class FakeUniverse:
        def __init__(self, _path: str) -> None:
            self.atoms = type("FakeAtoms", (), {"chainIDs": ["A"]})()

        def select_atoms(self, selection: str) -> FakeAtomGroup:
            if selection in {"protein and chainID A", "protein and name CA"}:
                return FakeAtomGroup(1)
            return FakeAtomGroup(0)

    fake_mda = type("FakeMDAnalysis", (), {"Universe": FakeUniverse})()
    monkeypatch.setitem(sys.modules, "MDAnalysis", fake_mda)

    output_sim_dir = tmp_path / "converted"
    output_sim_dir.mkdir()
    (output_sim_dir / "solvated_system.pdb").write_text("ATOM\n", encoding="utf-8")
    (output_sim_dir / "config.yaml").write_text("name: invalid\n", encoding="utf-8")
    metadata = convert_legacy.SimMetadata(
        folder_name="legacy",
        restraint_distance="10A",
        enzyme_name="LipA",
        substrate_name="Resorufin-Butyrate",
        temperature_K=363.0,
        replicate=1,
    )

    assert convert_legacy.validate_output(output_sim_dir, metadata) is False


def test_validate_output_missing_expected_segment_trajectory_returns_false(
    convert_legacy: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Validation should fail when a daisy-chain segment trajectory is missing."""

    class FakeAtomGroup:
        def __init__(self, size: int) -> None:
            self._size = size

        def __len__(self) -> int:
            return self._size

    class FakeUniverse:
        def __init__(self, _path: str) -> None:
            self.atoms = type("FakeAtoms", (), {"chainIDs": ["A"]})()

        def select_atoms(self, selection: str) -> FakeAtomGroup:
            if selection in {"protein and chainID A", "protein and name CA"}:
                return FakeAtomGroup(1)
            return FakeAtomGroup(0)

    class FakeConfig:
        def get_working_directory(self, _replicate: int) -> Path:
            return output_sim_dir

    class FakeSimulationConfig:
        @staticmethod
        def from_yaml(_path: Path) -> FakeConfig:
            return FakeConfig()

    fake_mda = type("FakeMDAnalysis", (), {"Universe": FakeUniverse})()
    monkeypatch.setitem(sys.modules, "MDAnalysis", fake_mda)
    monkeypatch.setattr("polyzymd.config.schema.SimulationConfig", FakeSimulationConfig)

    output_sim_dir = tmp_path / "converted"
    output_sim_dir.mkdir()
    (output_sim_dir / "solvated_system.pdb").write_text("ATOM\n", encoding="utf-8")
    (output_sim_dir / "config.yaml").write_text("name: legacy\n", encoding="utf-8")
    prod0 = output_sim_dir / "production_0"
    prod1 = output_sim_dir / "production_1"
    prod0.mkdir()
    prod1.mkdir()
    (prod0 / "production_0_trajectory.dcd").write_bytes(b"DCD")

    metadata = convert_legacy.SimMetadata(
        folder_name="legacy",
        restraint_distance="10A",
        enzyme_name="LipA",
        substrate_name="Resorufin-Butyrate",
        temperature_K=363.0,
        replicate=1,
        n_segments=2,
        production_dirs=[Path("production_0"), Path("production_1")],
    )

    assert convert_legacy.validate_output(output_sim_dir, metadata) is False


def test_convert_simulation_returns_false_when_validation_fails(
    convert_legacy: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Conversion should fail when converted output validation fails."""

    def fake_discover_files(_sim_dir: Path, metadata: object) -> None:
        metadata.production_dirs = []
        metadata.topology_path = tmp_path / "topology.pdb"

    monkeypatch.setattr(convert_legacy, "discover_files", fake_discover_files)
    monkeypatch.setattr(convert_legacy, "read_parameters_json", lambda *_args: None)
    monkeypatch.setattr(convert_legacy, "rewrite_topology", lambda *_args: None)
    monkeypatch.setattr(convert_legacy, "generate_config_yaml", lambda *_args: None)
    monkeypatch.setattr(convert_legacy, "create_symlinks", lambda *_args: None)
    monkeypatch.setattr(convert_legacy, "validate_output", lambda *_args: False)
    sim_dir = tmp_path / (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1"
    )
    sim_dir.mkdir()

    success = convert_legacy.convert_simulation(
        sim_dir=sim_dir,
        output_dir=tmp_path / "converted",
        reference_pdb=tmp_path / "reference.pdb",
    )

    assert success is False


def test_convert_simulation_returns_false_when_daisy_chain_segment_missing(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Conversion should fail before writing output for incomplete daisy chains."""
    sim_dir = tmp_path / (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1"
    )
    prod0 = sim_dir / "production_0"
    prod1 = sim_dir / "production_1"
    prod0.mkdir(parents=True)
    prod1.mkdir()
    (prod0 / "production_0_topology.pdb").write_text("ATOM\n", encoding="utf-8")
    (prod0 / "production_0_trajectory.dcd").write_bytes(b"DCD")

    success = convert_legacy.convert_simulation(
        sim_dir=sim_dir,
        output_dir=tmp_path / "converted",
        reference_pdb=tmp_path / "reference.pdb",
    )

    assert success is False


def test_conversion_guidance_uses_canonical_compare_commands(
    convert_legacy: ModuleType, caplog: pytest.LogCaptureFixture, tmp_path: Path
) -> None:
    """Generated guidance uses current comparison commands only."""
    caplog.set_level("INFO", logger=convert_legacy.logger.name)
    converted = tmp_path / "converted"
    sim_dir = converted / (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1"
    )
    sim_dir.mkdir(parents=True)
    (sim_dir / "config.yaml").write_text("name: legacy\n", encoding="utf-8")

    convert_legacy.generate_comparison_yaml(
        output_dir=converted,
        comparison_dir=tmp_path / "comparison",
        project_name="legacy_compare",
        control_label="No Polymer (Control)",
    )

    guidance = caplog.text
    assert "polyzymd compare validate -f" in guidance
    assert "polyzymd compare run rmsf -f" in guidance
    assert "polyzymd compare plot-all -f" in guidance
    assert "polyzymd analyze" not in guidance
    assert "analysis.yaml" not in guidance
