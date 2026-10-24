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


def test_a_replicate_folder_that_was_never_run_is_left_out(tmp_path: Path, caplog) -> None:
    """A folder a failed `polyzymd run` left behind does not drop the whole condition."""
    import logging

    pytest.importorskip("MDAnalysis")
    from polyzymd.analyses.study import Study
    from polyzymd.config.schema import SimulationConfig
    from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

    config = write_simulation_config(tmp_path / "cond", scratch=tmp_path / "scratch")
    for replicate in (1, 2):
        write_openmm_replicate(config, replicate, [1.0 + 0.01 * k for k in range(10)])
    empty = SimulationConfig.from_yaml(config).get_working_directory(3)
    empty.mkdir(parents=True)
    (empty / "progress.json").write_text("{}")

    with caplog.at_level(logging.WARNING):
        condition = Study.from_configs({"A": config}, equilibration="0ns")["A"]
    assert [replicate.index for replicate in condition.replicates] == [1, 2]
    assert condition.not_run == [3]
    assert "replicate 3 has no production output yet" in caplog.text


def test_a_named_replicate_that_was_never_run_is_still_used(tmp_path: Path) -> None:
    """Replicates named on purpose are kept, so a missing run is an error, not a gap."""
    pytest.importorskip("MDAnalysis")
    from polyzymd.analyses.study import Study
    from polyzymd.config.schema import SimulationConfig
    from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

    config = write_simulation_config(tmp_path / "cond", scratch=tmp_path / "scratch")
    write_openmm_replicate(config, 1, [1.0 + 0.01 * k for k in range(10)])
    SimulationConfig.from_yaml(config).get_working_directory(2).mkdir(parents=True)

    study = Study.from_configs({"A": config}, equilibration="0ns", replicates=[1, 2])
    assert [replicate.index for replicate in study["A"].replicates] == [1, 2]


def test_segments_that_do_not_join_are_not_blamed_on_the_topology() -> None:
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.shared.loader import TrajectoryLineageError
    from polyzymd.analyses.study import Replicate

    class Branched:
        def load_universe(self, index):
            raise TrajectoryLineageError("production_1 overlaps production_0")

    replicate = Replicate.__new__(Replicate)
    replicate._universe = None
    replicate.index = 1
    replicate.condition = type("C", (), {"_provider": Branched(), "label": "A"})()
    with pytest.raises(ProtocolError, match="Cannot join the production segments") as info:
        Replicate.universe(replicate)
    assert "analysis-topology" not in info.value.hint
