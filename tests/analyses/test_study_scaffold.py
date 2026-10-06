"""study init: where create_study records the runs of each condition."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from polyzymd.analyses.study_scaffold import create_study
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

pytest.importorskip("MDAnalysis")
pytestmark = [pytest.mark.filterwarnings("ignore"), pytest.mark.usefixtures("git_identity")]


def test_a_relative_scratch_is_relative_to_its_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """study init records where the runs are, not '.'."""
    folder = tmp_path / "water"
    config = write_simulation_config(folder, scratch=Path("."))
    (folder / "test.pdb").write_text("END\n")
    monkeypatch.chdir(folder)  # where a user runs it, so the runs land beside the config
    write_openmm_replicate(config, 1, [1.0, 1.1, 1.2])
    monkeypatch.chdir(tmp_path)
    root = tmp_path / "st"
    create_study(root, conditions={"Water": config}, equilibration="0ns")
    recorded = yaml.safe_load((root / "data.local.yaml").read_text())["Water"]
    assert Path(recorded) == folder.resolve()
