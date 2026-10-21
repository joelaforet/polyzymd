"""A repeated or a missing frame at a segment boundary is repaired; other defects still raise.

Each replicate is an OpenMM run directory of two production segments of four
atoms, each holding five frames 100 ps apart whose first frame is step
``istart``. Segment 1 starting at step 5 joins cleanly. Starting at step 4, it
records the time of segment 0's last frame again, as a restart from a state
saved shortly before that frame does. Starting at step 6, it leaves one frame
out, as a report write interrupted by a file-system error does.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest

from polyzymd.analyses.shared.loader import TrajectoryLineageError, TrajectoryLoader
from polyzymd.config.schema import SimulationConfig
from tests._support.analysis_testkit import CROSS, write_simulation_config

mda = pytest.importorskip("MDAnalysis")
pytestmark = [pytest.mark.filterwarnings("ignore::UserWarning")]

DT_PS = 100.0


def _write_segment(run_dir: Path, index: int, istart: int, n_frames: int = 5) -> None:
    """Write production_<index> with ``n_frames`` frames starting at step ``istart``."""
    segment = run_dir / f"production_{index}"
    segment.mkdir(parents=True, exist_ok=True)
    universe = mda.Universe.empty(4, n_residues=1, atom_resindex=[0] * 4, trajectory=True)
    universe.add_TopologyAttr("names", ["C1", "C2", "C3", "C4"])
    universe.add_TopologyAttr("resnames", ["MOL"])
    universe.atoms.positions = np.asarray(CROSS, dtype=np.float32)
    universe.atoms.write(str(run_dir / "solvated_system.pdb"))
    path = segment / f"production_{index}_trajectory.dcd"
    with mda.Writer(str(path), n_atoms=4, dt=DT_PS, istart=istart, nsavc=1) as writer:
        for k in range(n_frames):
            universe.atoms.positions = np.asarray(CROSS, dtype=np.float32) * (1.0 + 0.01 * k)
            writer.write(universe.atoms)


def _config(tmp_path: Path, starts: list[int]) -> Path:
    config = write_simulation_config(tmp_path / "cond", scratch=tmp_path / "scratch")
    run_dir = SimulationConfig.from_yaml(config).get_working_directory(1)
    for index, istart in enumerate(starts):
        _write_segment(run_dir, index, istart)
    return config


def _replicate(config: Path, equilibration: str = "0ns"):
    from polyzymd.analyses.study import Study

    return Study.from_configs({"toy": config}, equilibration=equilibration)["toy"].replicates[0]


def test_clean_chain_has_no_repairs(tmp_path) -> None:
    replicate = _replicate(_config(tmp_path, [0, 5]))
    assert replicate.frames.tolist() == list(range(10))
    assert replicate.times * 1000 == pytest.approx(DT_PS * np.arange(10), abs=1e-3)
    assert replicate.condition._provider.provenance_for(1).segment_join is None


def test_repeated_boundary_step_drops_the_earlier_segments_last_frame(tmp_path, caplog) -> None:
    with caplog.at_level(logging.WARNING):
        replicate = _replicate(_config(tmp_path, [0, 4]))
        frames = replicate.frames
    assert len(replicate.universe().trajectory) == 10
    assert frames.tolist() == [0, 1, 2, 3, 5, 6, 7, 8, 9]
    assert replicate.times * 1000 == pytest.approx(DT_PS * np.arange(9), abs=1e-3)
    assert any("leaving that frame out" in record.message for record in caplog.records)


def test_one_missing_frame_is_accepted_and_later_times_stay_recorded(tmp_path) -> None:
    replicate = _replicate(_config(tmp_path, [0, 6]))
    assert replicate.frames.tolist() == list(range(10))
    expected = DT_PS * np.array([0, 1, 2, 3, 4, 6, 7, 8, 9, 10])
    assert replicate.times * 1000 == pytest.approx(expected, abs=1e-3)


def test_equilibration_cut_uses_recorded_times_after_a_missing_frame(tmp_path) -> None:
    """At 600 ps, index 5 is the first production frame, not index 6 as 600/100 gives."""
    replicate = _replicate(_config(tmp_path, [0, 6]), equilibration="600ps")
    assert replicate.frames.tolist() == [5, 6, 7, 8, 9]
    assert replicate.times[0] * 1000 == pytest.approx(600.0, abs=1e-3)


def test_provenance_records_the_repair(tmp_path) -> None:
    replicate = _replicate(_config(tmp_path, [0, 4]))
    assert len(replicate.frames) == 9
    record = replicate.condition._provider.provenance_for(1, refresh=True).as_dict()
    assert record["segment_join"]["dropped_frames"] == [4]
    assert record["segment_join"]["dropped_frame_times_ps"] == pytest.approx([400.0], abs=1e-3)
    assert record["segment_join"]["missing_before"] == []
    assert any("Left out frame 4" in text for text in record["warnings"])


@pytest.mark.parametrize("second_start", [3, 7])
def test_larger_overlaps_and_gaps_still_raise(tmp_path, second_start) -> None:
    loader = TrajectoryLoader(SimulationConfig.from_yaml(_config(tmp_path, [0, second_start])))
    with pytest.raises(TrajectoryLineageError, match="repaired"):
        loader.load_universe(1, cache=False)
