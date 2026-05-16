"""Tests for the shared pair-distance MDAnalysis primitive."""

from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace
from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.mda import PairDistanceSpec, build_pair_distance_analysis
from polyzymd.analyses.shared.selections import SelectionMode


class _FakeTrajectory:
    """Trajectory fake that exposes frame-specific timestep metadata."""

    def __init__(self, dimensions: list[Any], timestep_ps: float = 2.0) -> None:
        """Store fake timestep dimensions for each frame."""

        self.dimensions = dimensions
        self.timestep_ps = timestep_ps
        self.current_frame = 0

    def __len__(self) -> int:
        """Return the number of fake frames."""

        return len(self.dimensions)

    def __getitem__(self, frame: int) -> SimpleNamespace:
        """Return a fake timestep and update the current frame."""

        self.current_frame = int(frame)
        return SimpleNamespace(
            frame=int(frame),
            time=float(frame) * self.timestep_ps,
            dimensions=self.dimensions[int(frame)],
        )


class _FakeAnalysisBase:
    """Minimal ``AnalysisBase`` fake that drives the analysis frame loop."""

    def __init__(self, trajectory: _FakeTrajectory) -> None:
        """Store the trajectory and allocate a results namespace."""

        self._trajectory = trajectory
        self.results = SimpleNamespace()

    def run(self, start: int = 0, stop: int | None = None, step: int = 1, **kwargs: Any) -> Any:
        """Run the fake AnalysisBase lifecycle for selected frames."""

        if kwargs:
            raise ValueError(f"Unexpected backend kwargs: {sorted(kwargs)}")
        if stop is None:
            stop = len(self._trajectory)
        self.frames = np.asarray(list(range(start, stop, step)), dtype=np.int64)
        self.times = self.frames.astype(np.float64) * float(self._trajectory.timestep_ps)
        self._prepare()
        for frame in self.frames:
            self._ts = self._trajectory[int(frame)]
            self._single_frame()
        self._conclude()
        return self


class _FakeAtomGroup:
    """One-atom group with frame-dependent positions."""

    def __init__(self, trajectory: _FakeTrajectory, positions: list[list[float]]) -> None:
        """Store positions indexed by the trajectory current frame."""

        self._trajectory = trajectory
        self._positions = np.asarray(positions, dtype=np.float64)

    def __len__(self) -> int:
        """Return one atom for single-position reductions."""

        return 1

    @property
    def positions(self) -> np.ndarray:
        """Return the current frame position as an atom-position array."""

        return np.asarray([self._positions[self._trajectory.current_frame]], dtype=np.float64)


@pytest.fixture
def fake_mdanalysis_pair_modules(monkeypatch: pytest.MonkeyPatch) -> list[Any]:
    """Install fake MDAnalysis modules and record ``calc_bonds`` box arguments."""

    calc_bonds_boxes: list[Any] = []

    def calc_bonds(positions_a: np.ndarray, positions_b: np.ndarray, box: Any = None) -> np.ndarray:
        """Record the forwarded PBC box and compute Euclidean pair distances."""

        calc_bonds_boxes.append(None if box is None else np.asarray(box).copy())
        return np.linalg.norm(np.asarray(positions_b) - np.asarray(positions_a), axis=1)

    mda_module = ModuleType("MDAnalysis")
    analysis_module = ModuleType("MDAnalysis.analysis")
    base_module = ModuleType("MDAnalysis.analysis.base")
    lib_module = ModuleType("MDAnalysis.lib")
    distances_module = ModuleType("MDAnalysis.lib.distances")
    base_module.AnalysisBase = _FakeAnalysisBase
    distances_module.calc_bonds = calc_bonds

    monkeypatch.setitem(sys.modules, "MDAnalysis", mda_module)
    monkeypatch.setitem(sys.modules, "MDAnalysis.analysis", analysis_module)
    monkeypatch.setitem(sys.modules, "MDAnalysis.analysis.base", base_module)
    monkeypatch.setitem(sys.modules, "MDAnalysis.lib", lib_module)
    monkeypatch.setitem(sys.modules, "MDAnalysis.lib.distances", distances_module)
    return calc_bonds_boxes


def test_pair_distance_analysis_collects_pair_frame_matrix_and_forwards_boxes(
    fake_mdanalysis_pair_modules: list[Any],
) -> None:
    """The custom AnalysisBase should collect pair x frame arrays via calc_bonds."""

    boxes = [np.asarray([10.0, 11.0, 12.0, 90.0, 90.0, 90.0]) + idx for idx in range(3)]
    trajectory = _FakeTrajectory(boxes)
    universe = SimpleNamespace(trajectory=trajectory)
    pairs = [
        PairDistanceSpec(
            label="Configured A",
            selection_a="sel a1",
            selection_b="sel b1",
            atoms_a=_FakeAtomGroup(trajectory, [[0.0, 0.0, 0.0]] * 3),
            atoms_b=_FakeAtomGroup(trajectory, [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]),
            mode_a=SelectionMode.SINGLE,
            mode_b=SelectionMode.SINGLE,
            threshold=2.5,
        ),
        PairDistanceSpec(
            label="Configured B",
            selection_a="sel a2",
            selection_b="sel b2",
            atoms_a=_FakeAtomGroup(trajectory, [[0.0, 1.0, 0.0]] * 3),
            atoms_b=_FakeAtomGroup(trajectory, [[0.0, 3.0, 0.0], [0.0, 4.0, 0.0], [0.0, 5.0, 0.0]]),
            mode_a=SelectionMode.SINGLE,
            mode_b=SelectionMode.SINGLE,
            threshold=None,
        ),
    ]

    analysis = build_pair_distance_analysis(universe=universe, pairs=pairs, use_pbc=True).run(
        start=0, stop=3, step=1
    )

    np.testing.assert_allclose(
        analysis.results.distance_matrix,
        np.asarray([[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]], dtype=np.float64),
    )
    assert analysis.results.distance_matrix.shape == (2, 3)
    np.testing.assert_array_equal(analysis.results.frames, np.asarray([0, 1, 2], dtype=np.int64))
    np.testing.assert_allclose(analysis.results.times_ps, np.asarray([0.0, 2.0, 4.0]))
    assert analysis.results.warnings == []
    assert len(fake_mdanalysis_pair_modules) == 3
    for observed_box, expected_box in zip(fake_mdanalysis_pair_modules, boxes, strict=True):
        np.testing.assert_allclose(observed_box, expected_box)


@pytest.mark.parametrize(
    ("dimensions", "message"),
    [
        ([None, None], "no box dimensions"),
        ([np.asarray([0.0, 10.0, 10.0, 90.0, 90.0, 90.0])] * 2, "box is invalid"),
    ],
)
def test_pair_distance_analysis_warns_once_and_disables_invalid_pbc_boxes(
    fake_mdanalysis_pair_modules: list[Any],
    dimensions: list[Any],
    message: str,
) -> None:
    """Missing or invalid boxes should disable PBC for affected frames with one warning."""

    trajectory = _FakeTrajectory(dimensions)
    universe = SimpleNamespace(trajectory=trajectory)
    pairs = [
        PairDistanceSpec(
            label="Configured label",
            selection_a="sel a",
            selection_b="sel b",
            atoms_a=_FakeAtomGroup(trajectory, [[0.0, 0.0, 0.0]] * 2),
            atoms_b=_FakeAtomGroup(trajectory, [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            mode_a=SelectionMode.SINGLE,
            mode_b=SelectionMode.SINGLE,
        )
    ]

    analysis = build_pair_distance_analysis(universe=universe, pairs=pairs, use_pbc=True).run(
        start=0, stop=2, step=1
    )

    assert fake_mdanalysis_pair_modules == [None, None]
    assert len(analysis.results.warnings) == 1
    assert message in analysis.results.warnings[0]
