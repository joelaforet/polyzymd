"""EnzymeBuilder: what it passes to OpenFF when it loads the enzyme."""

from __future__ import annotations

from pathlib import Path

_CYSTINE_TEMPLATES = (
    Path(__file__).resolve().parents[2]
    / "examples"
    / "pdb_preparation"
    / "4cha"
    / "nterminal_cystine_substructure.json"
)


def _config_with_templates(tmp_path: Path, templates):
    """Load a config whose enzyme.custom_substructures_path is templates.json beside it."""
    import shutil

    import yaml

    from polyzymd.config.schema import SimulationConfig
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    (tmp_path / "c" / "test.pdb").write_text("END\n")
    if isinstance(templates, Path):
        shutil.copy(templates, tmp_path / "c" / "templates.json")
    else:
        (tmp_path / "c" / "templates.json").write_text(templates)
    data = yaml.safe_load(path.read_text())
    data["enzyme"]["custom_substructures_path"] = "templates.json"
    path.write_text(yaml.safe_dump(data))
    return SimulationConfig.from_yaml(path)


def test_the_templates_reach_openff(tmp_path: Path, monkeypatch) -> None:
    """The custom substructures file reaches Topology.from_pdb as _custom_substructures."""
    from openff.toolkit import Topology

    from polyzymd.builders.enzyme import EnzymeBuilder

    seen = {}

    def from_pdb(path, **kwargs):
        seen.update(kwargs)
        return Topology()

    monkeypatch.setattr(Topology, "from_pdb", staticmethod(from_pdb))
    config = _config_with_templates(tmp_path, _CYSTINE_TEMPLATES)
    (tmp_path / "c" / "test.pdb").write_text("END\n")
    EnzymeBuilder().build_from_config(config.enzyme)
    assert list(seen["_custom_substructures"]) == ["NCYX"]
