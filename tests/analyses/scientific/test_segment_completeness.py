"""Tests for excluding production segments that are still being written.

A segment marked ``running`` or ``failed`` in ``progress.json`` has a
trajectory file that ends wherever the last flush landed. Reading it gives an
analysis window that is short for a reason the result does not record.
"""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace
from typing import Sequence

import pytest

from polyzymd.engines.openmm.engine import OpenMMEngine
from polyzymd.simulation.progress import (
    SegmentRecord,
    SegmentStatus,
    SimulationProgress,
    SimulationStatus,
    save_progress,
)


def _make_engine() -> OpenMMEngine:
    """Create an OpenMM engine with the production settings the layout needs.

    Returns
    -------
    OpenMMEngine
        Engine bound to a minimal stand-in configuration.
    """

    config = SimpleNamespace(
        simulation_phases=SimpleNamespace(
            production=SimpleNamespace(duration=100.0, time_step=2.0, samples=250)
        )
    )
    return OpenMMEngine(config=config)


def _write_segments(working_dir: Path, statuses: Sequence[SegmentStatus]) -> None:
    """Write segment directories, trajectories, and a matching ``progress.json``.

    Parameters
    ----------
    working_dir : Path
        Replicate working directory to populate.
    statuses : sequence of SegmentStatus
        Status recorded for segment ``0`` through ``len(statuses) - 1``.
    """

    segments = []
    for index, status in enumerate(statuses):
        segment_dir = working_dir / f"production_{index}"
        segment_dir.mkdir(parents=True, exist_ok=True)
        (segment_dir / f"production_{index}_trajectory.dcd").write_bytes(b"DCD" * (index + 1))
        segments.append(
            SegmentRecord(
                index=index,
                steps_completed=500_000,
                steps_requested=500_000,
                samples_written=250,
                status=status,
                duration_ns=1.0,
            )
        )
    save_progress(
        working_dir,
        SimulationProgress(
            total_steps_requested=1_000_000,
            total_samples_requested=500,
            timestep_fs=2.0,
            segments=segments,
            status=SimulationStatus.RUNNING,
        ),
    )


class TestSegmentStatusFiltering:
    """The trajectory resolver must not hand analyses a segment being written."""

    def test_running_segment_is_excluded_and_recorded(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A RUNNING segment is dropped, logged, and listed on the layout."""

        _write_segments(tmp_path, [SegmentStatus.COMPLETED, SegmentStatus.RUNNING])

        with caplog.at_level(logging.WARNING, logger="polyzymd.engines.openmm.engine"):
            layout = _make_engine().resolve_trajectory_layout(tmp_path, replicate=1)

        assert [path.name for path in layout.trajectory_paths] == ["production_0_trajectory.dcd"]
        assert layout.excluded_segments == [1]
        assert layout.segment_status == {0: "completed", 1: "running"}
        assert "production_1" in caplog.text

    def test_require_complete_false_includes_and_flags_segment(self, tmp_path: Path) -> None:
        """With ``require_complete=False`` the segment loads but stays flagged."""

        _write_segments(tmp_path, [SegmentStatus.COMPLETED, SegmentStatus.RUNNING])

        layout = _make_engine().resolve_trajectory_layout(
            tmp_path, replicate=1, require_complete=False
        )

        assert [path.name for path in layout.trajectory_paths] == [
            "production_0_trajectory.dcd",
            "production_1_trajectory.dcd",
        ]
        assert layout.excluded_segments == []
        assert layout.incomplete_segments == [1]
        assert layout.segment_status == {0: "completed", 1: "running"}

    def test_completed_segments_are_all_included(self, tmp_path: Path) -> None:
        """Completed chains are unchanged by the status check."""

        _write_segments(tmp_path, [SegmentStatus.COMPLETED, SegmentStatus.COMPLETED])

        layout = _make_engine().resolve_trajectory_layout(tmp_path, replicate=1)

        assert len(layout.trajectory_paths) == 2
        assert layout.excluded_segments == []

    def test_interrupted_segment_is_kept(self, tmp_path: Path) -> None:
        """Interrupted segments stay in the chain the continuation resumes from."""

        _write_segments(tmp_path, [SegmentStatus.INTERRUPTED, SegmentStatus.COMPLETED])

        layout = _make_engine().resolve_trajectory_layout(tmp_path, replicate=1)

        assert len(layout.trajectory_paths) == 2
        assert layout.excluded_segments == []
        assert layout.segment_status == {0: "interrupted", 1: "completed"}

    def test_missing_progress_file_keeps_every_segment(self, tmp_path: Path) -> None:
        """Runs with no ``progress.json`` keep the previous inclusive behaviour."""

        for index in (0, 1):
            segment_dir = tmp_path / f"production_{index}"
            segment_dir.mkdir()
            (segment_dir / f"production_{index}_trajectory.dcd").write_bytes(b"DCD")

        layout = _make_engine().resolve_trajectory_layout(tmp_path, replicate=1)

        assert len(layout.trajectory_paths) == 2
        assert layout.excluded_segments == []
        assert layout.segment_status == {}


class TestExclusionWarningText:
    """The warning must say what the exclusion did to the time line."""

    def test_middle_exclusion_warning_describes_the_hole(self, tmp_path: Path) -> None:
        """Excluding a segment other than the last one leaves a gap, and says so."""

        from polyzymd.analyses.shared.loader import _segment_completeness_warning
        from polyzymd.engines.base import TrajectoryLayout

        layout = TrajectoryLayout(
            trajectory_format="dcd",
            topology_format="pdb",
            trajectory_paths=[Path("production_0.dcd"), Path("production_2.dcd")],
            segment_status={0: "completed", 1: "failed", 2: "completed"},
            excluded_segments=[1],
        )

        warning = _segment_completeness_warning(layout)

        assert warning is not None
        assert "gap" in warning
        assert "[1]" in warning

    def test_trailing_exclusion_warning_describes_a_short_window(self, tmp_path: Path) -> None:
        """Excluding the last segment shortens the window instead."""

        from polyzymd.analyses.shared.loader import _segment_completeness_warning
        from polyzymd.engines.base import TrajectoryLayout

        layout = TrajectoryLayout(
            trajectory_format="dcd",
            topology_format="pdb",
            trajectory_paths=[Path("production_0.dcd")],
            segment_status={0: "completed", 1: "running"},
            excluded_segments=[1],
        )

        warning = _segment_completeness_warning(layout)

        assert warning is not None
        assert "gap" not in warning
        assert "ends before" in warning


class TestLoaderRecordsExcludedSegments:
    """The exclusion must reach provenance, not only the engine layout."""

    def test_trajectory_info_and_provenance_record_the_exclusion(self, tmp_path: Path) -> None:
        """Excluded segments and segment status reach ``UniverseProvenance``."""

        from unittest.mock import MagicMock

        from polyzymd.analyses.mda.universe import UniverseProvider
        from polyzymd.analyses.shared.loader import TrajectoryLoader

        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        (run_dir / "solvated_system.pdb").write_text("ATOM")
        _write_segments(run_dir, [SegmentStatus.COMPLETED, SegmentStatus.RUNNING])
        config = MagicMock()
        config.engine = "openmm"
        config.get_working_directory.side_effect = lambda rep: tmp_path / f"run_{rep}"
        config.output.effective_scratch_directory = tmp_path

        info = TrajectoryLoader(config).get_trajectory_info(1)
        assert info.excluded_segments == [1]
        assert info.n_segments == 1
        assert any("Excluded production segment" in warning for warning in info.warnings)

        provenance = UniverseProvider.from_config(config).provenance_for(1)
        payload = provenance.as_dict()
        assert payload["excluded_segments"] == [1]
        assert payload["segment_status"] == {"0": "completed", "1": "running"}
