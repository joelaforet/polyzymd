"""Tests for the legacy conversion utility."""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path
from types import ModuleType

import pytest
import yaml

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
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_SBMA-EGMA-50%_38x_363.0K_0.5ns-NVT_1000.0ns-NPT_run2"
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


def test_discover_files_rejects_empty_daisy_chain_segment(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Daisy-chain discovery fails before assigning empty segment trajectories."""
    sim_dir = tmp_path / "legacy"
    prod0 = sim_dir / "production_0"
    prod1 = sim_dir / "production_1"
    prod0.mkdir(parents=True)
    prod1.mkdir()
    (prod0 / "production_0_topology.pdb").write_text("ATOM\n", encoding="utf-8")
    (prod0 / "production_0_trajectory.dcd").write_bytes(b"DCD")
    (prod1 / "production_1_trajectory.dcd").write_bytes(b"")

    metadata = convert_legacy.SimMetadata(
        folder_name="legacy",
        restraint_distance="10A",
        enzyme_name="LipA",
        substrate_name="Resorufin-Butyrate",
        temperature_K=363.0,
        replicate=1,
    )

    with pytest.raises(ValueError, match="production_1_trajectory.dcd"):
        convert_legacy.discover_files(sim_dir, metadata)

    assert metadata.trajectory_paths == []


def test_discover_files_rejects_non_contiguous_daisy_chain_segments(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Discovery should reject daisy-chain segment gaps before assigning metadata."""
    sim_dir = tmp_path / "legacy"
    prod0 = sim_dir / "production_0"
    prod2 = sim_dir / "production_2"
    prod0.mkdir(parents=True)
    prod2.mkdir()
    (prod0 / "production_0_topology.pdb").write_text("ATOM\n", encoding="utf-8")
    (prod0 / "production_0_trajectory.dcd").write_bytes(b"DCD")
    (prod2 / "production_2_trajectory.dcd").write_bytes(b"DCD")

    metadata = convert_legacy.SimMetadata(
        folder_name="legacy",
        restraint_distance="10A",
        enzyme_name="LipA",
        substrate_name="Resorufin-Butyrate",
        temperature_K=363.0,
        replicate=1,
    )

    with pytest.raises(FileNotFoundError) as exc_info:
        convert_legacy.discover_files(sim_dir, metadata)

    assert "production_1/production_1_trajectory.dcd" in str(exc_info.value)
    assert metadata.production_dirs == []
    assert metadata.trajectory_paths == []


def test_discover_files_rejects_daisy_chain_missing_production_zero(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Discovery should require daisy chains to start at production_0."""
    sim_dir = tmp_path / "legacy"
    prod1 = sim_dir / "production_1"
    prod2 = sim_dir / "production_2"
    prod1.mkdir(parents=True)
    prod2.mkdir()
    (prod1 / "production_1_trajectory.dcd").write_bytes(b"DCD")
    (prod2 / "production_2_trajectory.dcd").write_bytes(b"DCD")

    metadata = convert_legacy.SimMetadata(
        folder_name="legacy",
        restraint_distance="10A",
        enzyme_name="LipA",
        substrate_name="Resorufin-Butyrate",
        temperature_K=363.0,
        replicate=1,
    )

    with pytest.raises(FileNotFoundError) as exc_info:
        convert_legacy.discover_files(sim_dir, metadata)

    message = str(exc_info.value)
    assert "production_0/production_0_topology.pdb" in message
    assert "production_0/production_0_trajectory.dcd" in message


def test_discover_files_rejects_missing_first_segment_topology(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Discovery should require the exact production_0 topology path."""
    sim_dir = tmp_path / "legacy"
    prod0 = sim_dir / "production_0"
    prod1 = sim_dir / "production_1"
    prod0.mkdir(parents=True)
    prod1.mkdir()
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

    with pytest.raises(FileNotFoundError, match="production_0_topology.pdb"):
        convert_legacy.discover_files(sim_dir, metadata)


def test_discover_files_rejects_empty_single_production_trajectory(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Single-production discovery should reject empty trajectory files."""
    sim_dir = tmp_path / "legacy"
    prod = sim_dir / "production"
    prod.mkdir(parents=True)
    (prod / "production_topology.pdb").write_text("ATOM\n", encoding="utf-8")
    (prod / "production_trajectory.dcd").write_bytes(b"")

    metadata = convert_legacy.SimMetadata(
        folder_name="legacy",
        restraint_distance="10A",
        enzyme_name="LipA",
        substrate_name="Resorufin-Butyrate",
        temperature_K=363.0,
        replicate=1,
    )

    with pytest.raises(ValueError, match="Empty production trajectory"):
        convert_legacy.discover_files(sim_dir, metadata)


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

    with pytest.raises(FileNotFoundError, match="production_0_topology.pdb"):
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
    assert "report_interval" not in data["simulation_phases"]["production"]
    config = SimulationConfig.from_yaml(config_path)
    assert config.engine == "openmm"


def test_analyze_commands_write_the_triad_pairs_and_no_comparison_yaml(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """The commands read the converted configs; only the triad pairs file is written."""
    converted = tmp_path / "converted"
    sim_dir = converted / (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1"
    )
    sim_dir.mkdir(parents=True)
    (sim_dir / "config.yaml").write_text("name: legacy\n", encoding="utf-8")
    pairs_path = tmp_path / "analysis" / "triad_pairs.yaml"

    commands = convert_legacy.analyze_commands(
        output_dir=converted,
        control_label="No Polymer (Control)",
        pairs_path=pairs_path,
    )

    pairs = yaml.safe_load(pairs_path.read_text())
    assert [pair["label"] for pair in pairs] == ["Ser77-His156", "His156-Asp133"]
    assert sorted(path.name for path in tmp_path.rglob("*.yaml")) == [
        "config.yaml",
        "triad_pairs.yaml",
    ]
    config = (sim_dir / "config.yaml").resolve()
    assert commands[1] == (
        f"polyzymd analyze distances -c {config} --label 'No Polymer (Control)' "
        f"--replicates 1 --eq 10ns --set pairs={pairs_path}"
    )


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


def test_analyze_commands_keep_every_disambiguated_condition(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Every disambiguated duplicate-label condition gets a -c pair, control first."""
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

    rmsf, distances, contacts = convert_legacy.analyze_commands(
        output_dir=converted,
        control_label="No Polymer (Control)",
    )

    assert contacts.startswith(
        f"polyzymd analyze contacts -c {control_config.resolve()} --label 'No Polymer (Control)' "
    )
    assert f"-c {polymer_38.resolve()} --label 'SBMA-EGMA 50% (" in contacts
    assert f"-c {polymer_77.resolve()} --label 'SBMA-EGMA 50% (" in contacts
    assert "38 chains" in contacts and "77 chains" in contacts
    assert contacts.endswith("--replicates 1 --eq 10ns")
    assert rmsf.endswith("--set reference_mode=average --set 'highlight_residues=[77,133,156]'")
    assert distances.endswith("--set pairs=<pairs.yaml>")


def test_analyze_commands_reject_ambiguous_control_base_label(
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
        convert_legacy.analyze_commands(output_dir=converted, control_label="No Polymer (Control)")

    message = str(exc_info.value)
    assert "No Polymer (Control) (10A, CALB" in message
    assert "No Polymer (Control) (12A, CALB" in message
    assert 'Pass --control "<exact label>"' in message


def test_analyze_commands_accept_explicit_disambiguated_control(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """An exact disambiguated control label is accepted and comes first."""
    converted = tmp_path / "converted"
    _write_converted_config(
        converted,
        "10A_RESTRAINT_CALB_Resorufin-Butyrate_conf1_363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )
    control = _write_converted_config(
        converted,
        "12A_RESTRAINT_CALB_Resorufin-Butyrate_conf2_363.0K_0.5ns-NVT_1000.0ns-NPT_run1",
    )
    control_label = "No Polymer (Control) (12A, CALB, Resorufin-Butyrate, conf2, 363K, no polymer)"

    commands = convert_legacy.analyze_commands(output_dir=converted, control_label=control_label)

    assert commands[2].startswith(
        f"polyzymd analyze contacts -c {control.resolve()} --label '{control_label}' "
    )


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
    convert_legacy: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
    tmp_path: Path,
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
    caplog.set_level("ERROR", logger=convert_legacy.logger.name)

    output_sim_dir = tmp_path / "converted"
    output_sim_dir.mkdir()
    (output_sim_dir / "solvated_system.pdb").write_text("ATOM\n", encoding="utf-8")
    (output_sim_dir / "config.yaml").write_text("name: legacy\n", encoding="utf-8")
    prod0 = output_sim_dir / "production_0"
    prod2 = output_sim_dir / "production_2"
    prod0.mkdir()
    prod2.mkdir()
    (prod0 / "production_0_trajectory.dcd").write_bytes(b"DCD")
    (prod2 / "production_2_trajectory.dcd").write_bytes(b"DCD")

    metadata = convert_legacy.SimMetadata(
        folder_name="legacy",
        restraint_distance="10A",
        enzyme_name="LipA",
        substrate_name="Resorufin-Butyrate",
        temperature_K=363.0,
        replicate=1,
        n_segments=3,
        production_dirs=[Path("production_0"), Path("production_2")],
    )

    assert convert_legacy.validate_output(output_sim_dir, metadata) is False
    assert "production_1/production_1_trajectory.dcd" in caplog.text


def test_validate_output_empty_expected_segment_trajectory_returns_false(
    convert_legacy: ModuleType,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
    tmp_path: Path,
) -> None:
    """Validation should fail when an expected segment trajectory is empty."""

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
    caplog.set_level("ERROR", logger=convert_legacy.logger.name)

    output_sim_dir = tmp_path / "converted"
    output_sim_dir.mkdir()
    (output_sim_dir / "solvated_system.pdb").write_text("ATOM\n", encoding="utf-8")
    (output_sim_dir / "config.yaml").write_text("name: legacy\n", encoding="utf-8")
    prod0 = output_sim_dir / "production_0"
    prod1 = output_sim_dir / "production_1"
    prod0.mkdir()
    prod1.mkdir()
    (prod0 / "production_0_trajectory.dcd").write_bytes(b"DCD")
    (prod1 / "production_1_trajectory.dcd").write_bytes(b"")

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
    assert "Trajectory is empty" in caplog.text
    assert "production_1/production_1_trajectory.dcd" in caplog.text


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


def test_convert_simulation_propagates_topology_validation_error(
    convert_legacy: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Topology validation errors should propagate after discovery succeeds."""

    def fake_discover_files(_sim_dir: Path, metadata: object) -> None:
        metadata.production_dirs = []
        metadata.topology_path = tmp_path / "topology.pdb"

    def fake_rewrite_topology(*_args: object) -> None:
        raise ValueError("invalid chain assignment")

    monkeypatch.setattr(convert_legacy, "discover_files", fake_discover_files)
    monkeypatch.setattr(convert_legacy, "read_parameters_json", lambda *_args: None)
    monkeypatch.setattr(convert_legacy, "rewrite_topology", fake_rewrite_topology)
    sim_dir = tmp_path / (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1"
    )
    sim_dir.mkdir()

    with pytest.raises(ValueError, match="invalid chain assignment"):
        convert_legacy.convert_simulation(
            sim_dir=sim_dir,
            output_dir=tmp_path / "converted",
            reference_pdb=tmp_path / "reference.pdb",
        )


def test_convert_simulation_propagates_config_validation_error(
    convert_legacy: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Config validation errors should propagate after discovery succeeds."""

    def fake_discover_files(_sim_dir: Path, metadata: object) -> None:
        metadata.production_dirs = []
        metadata.topology_path = tmp_path / "topology.pdb"

    def fake_generate_config_yaml(*_args: object) -> None:
        raise ValueError("invalid config")

    monkeypatch.setattr(convert_legacy, "discover_files", fake_discover_files)
    monkeypatch.setattr(convert_legacy, "read_parameters_json", lambda *_args: None)
    monkeypatch.setattr(convert_legacy, "rewrite_topology", lambda *_args: None)
    monkeypatch.setattr(convert_legacy, "generate_config_yaml", fake_generate_config_yaml)
    sim_dir = tmp_path / (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1"
    )
    sim_dir.mkdir()

    with pytest.raises(ValueError, match="invalid config"):
        convert_legacy.convert_simulation(
            sim_dir=sim_dir,
            output_dir=tmp_path / "converted",
            reference_pdb=tmp_path / "reference.pdb",
        )


def test_convert_simulation_returns_false_for_discovery_validation_error(
    convert_legacy: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Discovery validation errors should fail before creating output."""

    def fake_discover_files(*_args: object) -> None:
        raise ValueError("empty trajectory")

    monkeypatch.setattr(convert_legacy, "discover_files", fake_discover_files)
    sim_dir = tmp_path / (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1"
    )
    sim_dir.mkdir()
    converted = tmp_path / "converted"

    success = convert_legacy.convert_simulation(
        sim_dir=sim_dir,
        output_dir=converted,
        reference_pdb=tmp_path / "reference.pdb",
    )

    assert success is False
    assert not (converted / sim_dir.name).exists()


def test_convert_simulation_returns_false_for_gapped_daisy_chain_before_output(
    convert_legacy: ModuleType, tmp_path: Path
) -> None:
    """Conversion should reject non-contiguous daisy chains before output writes."""
    sim_dir = tmp_path / (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1"
    )
    prod0 = sim_dir / "production_0"
    prod2 = sim_dir / "production_2"
    prod0.mkdir(parents=True)
    prod2.mkdir()
    (prod0 / "production_0_topology.pdb").write_text("ATOM\n", encoding="utf-8")
    (prod0 / "production_0_trajectory.dcd").write_bytes(b"DCD")
    (prod2 / "production_2_trajectory.dcd").write_bytes(b"DCD")

    success = convert_legacy.convert_simulation(
        sim_dir=sim_dir,
        output_dir=tmp_path / "converted",
        reference_pdb=tmp_path / "reference.pdb",
    )

    assert success is False
    assert not (tmp_path / "converted" / sim_dir.name).exists()


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


def test_convert_simulation_returns_false_for_empty_trajectory_daisy_chain_segment(
    convert_legacy: ModuleType,
    caplog: pytest.LogCaptureFixture,
    tmp_path: Path,
) -> None:
    """Conversion should skip empty daisy-chain trajectories before output writes."""
    caplog.set_level("ERROR", logger=convert_legacy.logger.name)
    sim_dir = tmp_path / (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1"
    )
    prod0 = sim_dir / "production_0"
    prod1 = sim_dir / "production_1"
    prod0.mkdir(parents=True)
    prod1.mkdir()
    (prod0 / "production_0_topology.pdb").write_text("ATOM\n", encoding="utf-8")
    (prod0 / "production_0_trajectory.dcd").write_bytes(b"DCD")
    (prod1 / "production_1_trajectory.dcd").write_bytes(b"")
    converted = tmp_path / "converted"

    success = convert_legacy.convert_simulation(
        sim_dir=sim_dir,
        output_dir=converted,
        reference_pdb=tmp_path / "reference.pdb",
    )

    assert success is False
    assert not (converted / sim_dir.name).exists()
    assert "SKIP" in caplog.text
    assert "Empty daisy-chain trajectory" in caplog.text
    assert "production_1/production_1_trajectory.dcd" in caplog.text


def test_convert_simulation_returns_false_for_empty_trajectory_single_production(
    convert_legacy: ModuleType,
    caplog: pytest.LogCaptureFixture,
    tmp_path: Path,
) -> None:
    """Conversion should skip empty single-production trajectories before output writes."""
    caplog.set_level("ERROR", logger=convert_legacy.logger.name)
    sim_dir = tmp_path / (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1"
    )
    prod = sim_dir / "production"
    prod.mkdir(parents=True)
    (prod / "production_topology.pdb").write_text("ATOM\n", encoding="utf-8")
    (prod / "production_trajectory.dcd").write_bytes(b"")
    converted = tmp_path / "converted"

    success = convert_legacy.convert_simulation(
        sim_dir=sim_dir,
        output_dir=converted,
        reference_pdb=tmp_path / "reference.pdb",
    )

    assert success is False
    assert not (converted / sim_dir.name).exists()
    assert "SKIP" in caplog.text
    assert "Empty production trajectory" in caplog.text


def test_print_commands_logs_the_analyze_commands_only(
    convert_legacy: ModuleType,
    caplog: pytest.LogCaptureFixture,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """--print-commands logs polyzymd analyze commands and the triad routine, nothing else."""
    caplog.set_level("INFO", logger=convert_legacy.logger.name)
    converted = tmp_path / "converted"
    sim_dir = converted / (
        "10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1"
    )
    sim_dir.mkdir(parents=True)
    (sim_dir / "config.yaml").write_text("name: legacy\n", encoding="utf-8")
    pairs = tmp_path / "triad_pairs.yaml"
    monkeypatch.setattr(
        "sys.argv",
        [
            "convert_legacy.py",
            "--output-dir",
            str(converted),
            "--print-commands",
            "--triad-pairs",
            str(pairs),
        ],
    )

    with pytest.raises(SystemExit) as exit_info:
        convert_legacy.main()

    assert exit_info.value.code == 0
    guidance = caplog.text
    for name in ("rmsf", "distances", "contacts"):
        assert f"polyzymd analyze {name} -c " in guidance
    assert f"--set pairs={pairs}" in guidance
    assert "how_to/analysis_triad_quickstart.html" in guidance
    assert "compare" not in guidance
    assert "comparison.yaml" not in guidance
    assert not (tmp_path / "comparison").exists()
