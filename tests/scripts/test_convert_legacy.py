"""Tests for the legacy conversion utility."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

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
