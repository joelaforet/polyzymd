"""Study, Condition and Replicate: how a replicate reports a universe it cannot load."""

from __future__ import annotations

from pathlib import Path

import pytest


def test_a_topology_that_cannot_load_is_named_as_such(tmp_path: Path) -> None:
    """A broken topology is not reported as an equilibration problem."""
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.study import Condition, Replicate

    class Broken:
        def load_universe(self, index):
            raise ValueError("Length of charges does not match number of atoms")

    replicate = Replicate.__new__(Replicate)
    replicate._universe = None
    replicate.index = 1
    replicate.condition = type("C", (), {"_provider": Broken(), "label": "A"})()
    with pytest.raises(ProtocolError, match="Cannot load the topology") as info:
        Replicate.universe(replicate)
    assert "analysis-topology" in info.value.hint
    assert Condition  # the class used above is the study's own


def test_labels_that_share_a_folder_name_are_refused(tmp_path: Path) -> None:
    """Study.from_configs refuses two labels that give one results folder name."""
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.study import Study
    from tests._support.analysis_testkit import write_simulation_config

    configs = {
        label: write_simulation_config(tmp_path / name, scratch=tmp_path / "scratch" / name)
        for label, name in (("SBMA 50", "a"), ("SBMA 50%", "b"))
    }
    with pytest.raises(ProtocolError, match="SBMA 50.*SBMA 50%"):
        Study.from_configs(configs, equilibration="0ns")


def test_labels_differing_in_case_are_refused_only_in_a_study_file() -> None:
    """WT and wt keep separate results folders, but share one conditions/ folder name."""
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.study_file import check_folder_names

    check_folder_names(["WT", "wt"], "Study.from_configs", study_folders=False)
    with pytest.raises(ProtocolError, match="'WT' and 'wt'"):
        check_folder_names(["WT", "wt"], "study.yaml")
