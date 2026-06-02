"""Tests for shared trajectory alignment helpers."""

from __future__ import annotations

import sys
import types
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory


class _FakeTrajectory:
    """Minimal trajectory object with a finite frame count."""

    def __len__(self) -> int:
        """Return the fake trajectory length."""

        return 100


class _FakeUniverse:
    """Minimal universe object for alignment tests."""

    trajectory = _FakeTrajectory()

    def select_atoms(self, selection: str) -> "_FakeAtomGroup":
        """Return fake atoms for external-reference validation."""

        del selection
        return _FakeAtomGroup()


class _FakeResidue:
    """Minimal residue with a residue name."""

    resname = "ALA"


class _FakeAtomGroup:
    """Minimal atom group with one matching residue."""

    residues = [_FakeResidue()]

    def __len__(self) -> int:
        """Return the fake atom count."""

        return 1


@pytest.fixture
def fake_mdanalysis(monkeypatch: pytest.MonkeyPatch):
    """Install fake MDAnalysis alignment classes and capture run windows."""

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

    def make_universe(*args, **kwargs) -> _FakeUniverse:
        """Return a fake Universe for external references."""

        del args, kwargs
        return _FakeUniverse()

    mdanalysis_module = types.ModuleType("MDAnalysis")
    mdanalysis_module.Universe = make_universe
    analysis_module = types.ModuleType("MDAnalysis.analysis")
    align_module = types.ModuleType("MDAnalysis.analysis.align")
    align_module.AverageStructure = FakeAverageStructure
    align_module.AlignTraj = FakeAlignTraj
    analysis_module.align = align_module
    mdanalysis_module.analysis = analysis_module
    monkeypatch.setitem(sys.modules, "MDAnalysis", mdanalysis_module)
    monkeypatch.setitem(sys.modules, "MDAnalysis.analysis", analysis_module)
    monkeypatch.setitem(sys.modules, "MDAnalysis.analysis.align", align_module)

    return average_runs, align_runs


def test_average_reference_alignment_uses_selected_production_window(fake_mdanalysis) -> None:
    """Average-reference alignment should pass post-equilibration frame bounds."""

    average_runs, align_runs = fake_mdanalysis

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


def test_centroid_reference_alignment_uses_selected_production_window(
    fake_mdanalysis,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Centroid-reference alignment should pass post-equilibration frame bounds."""

    average_runs, align_runs = fake_mdanalysis

    from polyzymd.analyses.shared import centroid

    centroid_calls: list[dict[str, int | None]] = []

    def fake_find_centroid_frame(
        universe,
        *,
        selection: str,
        start_frame: int,
        stop_frame: int | None,
        verbose: bool,
    ) -> int:
        """Capture centroid selection arguments and return a frame."""

        del universe, selection, verbose
        centroid_calls.append({"start": start_frame, "stop": stop_frame})
        return 12

    monkeypatch.setattr(centroid, "find_centroid_frame", fake_find_centroid_frame)

    ref_frame = align_trajectory(
        _FakeUniverse(),
        AlignmentConfig(reference_mode="centroid"),
        start_frame=10,
        stop_frame=60,
        step_frame=5,
    )

    assert ref_frame == 12
    assert average_runs == []
    assert centroid_calls == [{"start": 10, "stop": 60}]
    assert align_runs == [{"start": 10, "stop": 60, "step": 5}]


def test_frame_reference_alignment_uses_selected_production_window(fake_mdanalysis) -> None:
    """Frame-reference alignment should pass post-equilibration frame bounds."""

    average_runs, align_runs = fake_mdanalysis

    ref_frame = align_trajectory(
        _FakeUniverse(),
        AlignmentConfig(reference_mode="frame", reference_frame=20),
        start_frame=10,
        stop_frame=60,
        step_frame=5,
    )

    assert ref_frame == 19
    assert average_runs == []
    assert align_runs == [{"start": 10, "stop": 60, "step": 5}]


def test_external_reference_alignment_uses_selected_production_window(
    fake_mdanalysis,
    tmp_path: Path,
) -> None:
    """External-reference alignment should pass post-equilibration frame bounds."""

    average_runs, align_runs = fake_mdanalysis
    reference_file = tmp_path / "reference.pdb"
    reference_file.write_text("MODEL\nENDMDL\n", encoding="utf-8")

    ref_frame = align_trajectory(
        _FakeUniverse(),
        AlignmentConfig(reference_mode="external", reference_file=reference_file),
        start_frame=10,
        stop_frame=60,
        step_frame=5,
    )

    assert ref_frame is None
    assert average_runs == []
    assert align_runs == [{"start": 10, "stop": 60, "step": 5}]
