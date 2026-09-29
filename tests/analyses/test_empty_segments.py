"""A production segment whose DCD holds no frame is skipped, and the rest must still join.

Each replicate is an OpenMM run directory of three production segments of
four atoms. Segment ``k`` holds five frames 100 ps apart whose first frame is
step ``istart``, so the segments join without a gap when each one starts one
frame after the previous one ends. Segment 1 is either a 0-byte file, as a
segment interrupted at start-up leaves it, or a DCD header with no frame.
"""

from __future__ import annotations

import importlib
import logging
from pathlib import Path

import numpy as np
import pytest

from polyzymd.analyses.shared.loader import TrajectoryLineageError, TrajectoryLoader
from polyzymd.config.schema import SimulationConfig
from tests._support.analysis_testkit import CROSS, write_simulation_config

mda = pytest.importorskip("MDAnalysis")
pytestmark = [pytest.mark.filterwarnings("ignore::UserWarning")]


def _write_segment(run_dir: Path, index: int, istart: int, n_frames: int, dt_ps: float) -> Path:
    """Write production_<index> with ``n_frames`` frames starting at step ``istart``."""
    segment = run_dir / f"production_{index}"
    segment.mkdir(parents=True, exist_ok=True)
    universe = mda.Universe.empty(4, n_residues=1, atom_resindex=[0] * 4, trajectory=True)
    universe.add_TopologyAttr("names", ["C1", "C2", "C3", "C4"])
    universe.add_TopologyAttr("resnames", ["MOL"])
    universe.atoms.positions = np.asarray(CROSS, dtype=np.float32)
    universe.atoms.write(str(run_dir / "solvated_system.pdb"))
    path = segment / f"production_{index}_trajectory.dcd"
    with mda.Writer(str(path), n_atoms=4, dt=dt_ps, istart=istart, nsavc=1) as writer:
        for k in range(n_frames):
            universe.atoms.positions = np.asarray(CROSS, dtype=np.float32) * (1.0 + 0.01 * k)
            writer.write(universe.atoms)
    return path


def _empty(path: Path, kind: str) -> None:
    """Leave ``path`` as a 0-byte file or as its DCD header with no frame."""
    if kind == "zero_bytes":
        path.write_bytes(b"")
        return
    data = path.read_bytes()
    title = int.from_bytes(data[92:96], "little")
    path.write_bytes(data[: 92 + 4 + title + 4 + 12])


def _loader(tmp_path: Path, starts: list[int], empty: set[int], kind: str = "zero_bytes"):
    config = write_simulation_config(tmp_path / "cond", scratch=tmp_path / "scratch")
    run_dir = SimulationConfig.from_yaml(config).get_working_directory(1)
    for index, istart in enumerate(starts):
        path = _write_segment(run_dir, index, istart, 5, 100.0)
        if index in empty:
            _empty(path, kind)
    return TrajectoryLoader(SimulationConfig.from_yaml(config))


@pytest.mark.parametrize("kind", ["zero_bytes", "header_only"])
def test_empty_middle_segment_is_skipped_with_a_warning(tmp_path, caplog, kind) -> None:
    """Segments 0 and 2 join, so the empty segment 1 is skipped and ten frames remain."""
    loader = _loader(tmp_path, [0, 999, 5], {1}, kind)
    with caplog.at_level(logging.WARNING):
        info = loader.get_trajectory_info(replicate=1)
    assert [path.parent.name for path in info.trajectory_files] == ["production_0", "production_2"]
    assert info.empty_segments == [1] and info.excluded_segments == []
    assert any("production_1_trajectory.dcd holds no frames" in r.message for r in caplog.records)
    assert any("Skipped production segment(s) [1]" in text for text in info.warnings)
    universe = loader.load_universe(1, cache=False)
    assert len(universe.trajectory) == 10
    times = [ts.time for ts in universe.trajectory]
    assert times == pytest.approx([100.0 * k for k in range(10)], abs=1e-3)


def test_empty_segment_between_segments_that_do_not_join_still_raises(tmp_path) -> None:
    """Segment 2 starts two frames late, so skipping segment 1 leaves a gap."""
    loader = _loader(tmp_path, [0, 999, 7], {1})
    with pytest.raises(TrajectoryLineageError, match=r"segment\(s\) \[1\] were skipped"):
        loader.load_universe(1, cache=False)


def test_all_segments_empty_still_raises(tmp_path) -> None:
    loader = _loader(tmp_path, [0, 5, 10], {0, 1, 2})
    with pytest.raises(FileNotFoundError, match="No production trajectory files found"):
        loader.get_trajectory_info(replicate=1)


def test_mixed_frame_intervals_still_raise_across_a_skipped_segment(tmp_path) -> None:
    config = write_simulation_config(tmp_path / "cond", scratch=tmp_path / "scratch")
    run_dir = SimulationConfig.from_yaml(config).get_working_directory(1)
    _write_segment(run_dir, 0, 0, 5, 100.0)
    _empty(_write_segment(run_dir, 1, 5, 5, 100.0), "zero_bytes")
    _write_segment(run_dir, 2, 3, 5, 200.0)
    loader = TrajectoryLoader(SimulationConfig.from_yaml(config))
    with pytest.raises(TrajectoryLineageError, match="frame interval"):
        loader.load_universe(1, cache=False)


def test_universe_provenance_records_the_empty_segment(tmp_path) -> None:
    from polyzymd.analyses.mda.universe import UniverseProvider

    config = write_simulation_config(tmp_path / "cond", scratch=tmp_path / "scratch")
    run_dir = SimulationConfig.from_yaml(config).get_working_directory(1)
    for index, istart in enumerate([0, 999, 5]):
        path = _write_segment(run_dir, index, istart, 5, 100.0)
        if index == 1:
            _empty(path, "zero_bytes")
    provider = UniverseProvider.from_config(SimulationConfig.from_yaml(config))
    record = provider.provenance_for(1).as_dict()
    assert record["empty_segments"] == [1] and record["n_segments"] == 2
    assert importlib.import_module("polyzymd.engines.openmm.engine")._dcd_has_no_frames(
        run_dir / "production_1" / "production_1_trajectory.dcd"
    )
