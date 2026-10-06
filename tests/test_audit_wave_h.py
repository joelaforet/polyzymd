"""Regressions for the 1.3 audit round, wave H (findings of the docs rewrite).

Each test names its finding in
``PAPERS/polyzymd_v1.3_refactor_handoff/audit_2026-10-06/AUDIT_LOG.md``.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest
import yaml

from tests._support.analysis_testkit import write_simulation_config

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "pdb_preparation" / "4cha"


def _config_with_templates(tmp_path: Path, templates):
    from polyzymd.config.schema import SimulationConfig

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


class TestCustomSubstructures:
    def test_the_config_takes_the_file_relative_to_itself(self, tmp_path: Path) -> None:
        """G-1: the documented key exists and resolves like pdb_path."""
        config = _config_with_templates(tmp_path, EXAMPLE / "nterminal_cystine_substructure.json")
        assert config.enzyme.custom_substructures_path == (tmp_path / "c" / "templates.json").resolve()

    def test_a_file_of_another_shape_is_refused(self, tmp_path: Path) -> None:
        from pydantic import ValidationError

        with pytest.raises(ValidationError, match="must map each residue name"):
            _config_with_templates(tmp_path, json.dumps({"NCYX": ["N", "CA"]}))

    def test_the_templates_reach_openff(self, tmp_path: Path, monkeypatch) -> None:
        from openff.toolkit import Topology

        from polyzymd.builders.enzyme import EnzymeBuilder

        seen = {}

        def from_pdb(path, **kwargs):
            seen.update(kwargs)
            return Topology()

        monkeypatch.setattr(Topology, "from_pdb", staticmethod(from_pdb))
        config = _config_with_templates(tmp_path, EXAMPLE / "nterminal_cystine_substructure.json")
        (tmp_path / "c" / "test.pdb").write_text("END\n")
        EnzymeBuilder().build_from_config(config.enzyme)
        assert list(seen["_custom_substructures"]) == ["NCYX"]

    def test_the_templates_are_part_of_the_config_hash(self, tmp_path: Path) -> None:
        from polyzymd.analyses.identity import compute_config_hash

        one = _config_with_templates(tmp_path / "a", json.dumps({"X": {"[#6:1]": ["C1"]}}))
        two = _config_with_templates(tmp_path / "b", json.dumps({"X": {"[#6:1]": ["C2"]}}))
        assert compute_config_hash(one) != compute_config_hash(two)
