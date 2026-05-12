"""Tests for shared trajectory alignment helpers."""

from __future__ import annotations

import sys
import types
from unittest.mock import MagicMock

from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory


class _FakeTrajectory:
    """Minimal trajectory object with a finite frame count."""

    def __len__(self) -> int:
        """Return the fake trajectory length."""

        return 100


class _FakeUniverse:
    """Minimal universe object for alignment tests."""

    trajectory = _FakeTrajectory()


def test_average_reference_alignment_uses_selected_production_window(monkeypatch) -> None:
    """Average-reference alignment should pass post-equilibration frame bounds."""

    average_runs: list[dict[str, int | None]] = []
    align_runs: list[dict[str, int | None]] = []

    class FakeAverageStructure:
        """Capture AverageStructure construction and run keyword arguments."""

        def __init__(self, *args, **kwargs) -> None:
            """Store construction arguments for assertions."""

            self.args = args
            self.kwargs = kwargs
            self.results = MagicMock(universe=_FakeUniverse())

        def run(self, *, start: int, stop: int | None, step: int) -> "FakeAverageStructure":
            """Record selected-window run arguments."""

            average_runs.append({"start": start, "stop": stop, "step": step})
            return self

    class FakeAlignTraj:
        """Capture AlignTraj construction and run keyword arguments."""

        def __init__(self, *args, **kwargs) -> None:
            """Store construction arguments for assertions."""

            self.args = args
            self.kwargs = kwargs

        def run(self, *, start: int, stop: int | None, step: int) -> "FakeAlignTraj":
            """Record selected-window run arguments."""

            align_runs.append({"start": start, "stop": stop, "step": step})
            return self

    mdanalysis_module = types.ModuleType("MDAnalysis")
    analysis_module = types.ModuleType("MDAnalysis.analysis")
    align_module = types.ModuleType("MDAnalysis.analysis.align")
    align_module.AverageStructure = FakeAverageStructure
    align_module.AlignTraj = FakeAlignTraj
    analysis_module.align = align_module
    mdanalysis_module.analysis = analysis_module
    monkeypatch.setitem(sys.modules, "MDAnalysis", mdanalysis_module)
    monkeypatch.setitem(sys.modules, "MDAnalysis.analysis", analysis_module)
    monkeypatch.setitem(sys.modules, "MDAnalysis.analysis.align", align_module)

    ref_frame = align_trajectory(
        _FakeUniverse(),
        AlignmentConfig(reference_mode="average"),
        start_frame=10,
        stop_frame=60,
        step_frame=5,
    )

    assert ref_frame is None
    assert average_runs == [{"start": 10, "stop": 60, "step": 5}]
    assert align_runs == [{"start": 10, "stop": 60, "step": 5}]
