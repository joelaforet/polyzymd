"""Tests for simulation progress tracking.

Covers:
- SegmentRecord and SimulationProgress Pydantic models
- Computed properties (total_steps_completed, is_complete, etc.)
- save_progress / load_progress (atomic write + JSON round-trip)
- scan_filesystem (reconstruct progress from production_N/ dirs)
- get_next_segment_info (determine what work remains)
- validate_progress (reconcile file vs filesystem)
- load_or_scan_progress (primary entry point)
- Unit conversion helpers (_convert_to_ns, _convert_to_fs)
"""

import json
from pathlib import Path

import pytest

from polyzymd.simulation.progress import (
    SegmentRecord,
    SegmentStatus,
    SimulationProgress,
    SimulationStatus,
    get_next_segment_info,
    load_or_scan_progress,
    load_progress,
    save_progress,
    scan_filesystem,
    validate_progress,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_segment(
    index: int,
    steps: int = 5000000,
    status: SegmentStatus = SegmentStatus.COMPLETED,
    duration_ns: float = 10.0,
    samples: int = 125,
) -> SegmentRecord:
    """Create a SegmentRecord for testing."""
    return SegmentRecord(
        index=index,
        steps_completed=steps,
        steps_requested=steps,
        samples_written=samples,
        status=status,
        duration_ns=duration_ns,
    )


def _make_progress(
    segments: list[SegmentRecord] | None = None,
    total_steps: int = 10000000,
    total_samples: int = 250,
) -> SimulationProgress:
    """Create a SimulationProgress for testing."""
    return SimulationProgress(
        config_path="/tmp/config.yaml",
        total_steps_requested=total_steps,
        total_samples_requested=total_samples,
        timestep_fs=2.0,
        segments=segments or [],
        replicate=1,
    )


def _write_completed_segment_on_disk(
    working_dir: Path,
    seg_idx: int,
    duration_ns: float = 10.0,
    num_samples: int = 125,
) -> None:
    """Write minimal files simulating a completed production segment."""
    seg_dir = working_dir / f"production_{seg_idx}"
    seg_dir.mkdir(parents=True, exist_ok=True)
    (seg_dir / f"production_{seg_idx}_state.xml").write_text("<mock-state/>")
    (seg_dir / f"production_{seg_idx}_system.xml").write_text("<mock-system/>")

    params = {
        "__values__": {
            "integ_params": {
                "__values__": {
                    "num_samples": num_samples,
                    "time_step": {
                        "__class__": "Quantity",
                        "__values__": {"value": 2.0, "unit": "femtosecond"},
                    },
                    "total_time": {
                        "__class__": "Quantity",
                        "__values__": {"value": duration_ns, "unit": "nanosecond"},
                    },
                }
            },
            "thermo_params": {"__values__": {}},
        }
    }
    (seg_dir / f"production_{seg_idx}_parameters.json").write_text(json.dumps(params))


def _write_interrupted_segment_on_disk(
    working_dir: Path,
    seg_idx: int,
    steps_completed: int = 3000000,
    total_steps: int = 5000000,
) -> None:
    """Write an INTERRUPTED marker in a production segment directory."""
    seg_dir = working_dir / f"production_{seg_idx}"
    seg_dir.mkdir(parents=True, exist_ok=True)
    remaining = total_steps - steps_completed
    (seg_dir / "INTERRUPTED").write_text(
        f"segment_index={seg_idx}\n"
        f"steps_completed={steps_completed}\n"
        f"total_steps={total_steps}\n"
        f"remaining_steps={remaining}\n"
    )


# ---------------------------------------------------------------------------
# SegmentRecord model
# ---------------------------------------------------------------------------


class TestSegmentRecord:
    """SegmentRecord Pydantic model."""

    def test_defaults(self):
        rec = SegmentRecord(index=0)
        assert rec.steps_completed == 0
        assert rec.status == SegmentStatus.RUNNING
        assert rec.duration_ns == 0.0

    def test_completed_segment(self):
        rec = _make_segment(0, steps=5000000, status=SegmentStatus.COMPLETED)
        assert rec.status == SegmentStatus.COMPLETED
        assert rec.steps_completed == 5000000

    def test_interrupted_segment(self):
        rec = _make_segment(1, steps=3000000, status=SegmentStatus.INTERRUPTED)
        assert rec.status == SegmentStatus.INTERRUPTED


# ---------------------------------------------------------------------------
# SimulationProgress model — computed properties
# ---------------------------------------------------------------------------


class TestSimulationProgress:
    """SimulationProgress computed properties."""

    def test_empty_progress(self):
        p = _make_progress()
        assert p.total_steps_completed == 0
        assert p.total_samples_completed == 0
        assert p.next_segment_index == 0
        assert p.is_complete is False
        assert p.fraction_complete() == 0.0
        assert p.last_completed_segment is None
        assert p.last_segment is None

    def test_one_completed_segment(self):
        seg = _make_segment(0, steps=5000000)
        p = _make_progress(segments=[seg], total_steps=10000000)
        assert p.total_steps_completed == 5000000
        assert p.next_segment_index == 1
        assert p.is_complete is False
        assert p.fraction_complete() == pytest.approx(0.5)
        assert p.steps_remaining == 5000000

    def test_two_segments_complete(self):
        segs = [_make_segment(0, steps=5000000), _make_segment(1, steps=5000000)]
        p = _make_progress(segments=segs, total_steps=10000000)
        assert p.total_steps_completed == 10000000
        assert p.is_complete is True
        assert p.fraction_complete() == pytest.approx(1.0)
        assert p.steps_remaining == 0

    def test_interrupted_segment_counted_toward_total(self):
        """Interrupted segments SHOULD count toward total_steps_completed."""
        segs = [
            _make_segment(0, steps=5000000, status=SegmentStatus.COMPLETED),
            _make_segment(1, steps=3000000, status=SegmentStatus.INTERRUPTED),
        ]
        p = _make_progress(segments=segs, total_steps=10000000)
        assert p.total_steps_completed == 8000000  # Both segments counted
        assert p.next_segment_index == 2
        assert p.is_complete is False

    def test_last_completed_segment(self):
        segs = [
            _make_segment(0, steps=5000000),
            _make_segment(1, steps=3000000),
        ]
        p = _make_progress(segments=segs)
        last = p.last_completed_segment
        assert last is not None
        assert last.index == 1

    def test_last_segment_includes_interrupted(self):
        segs = [
            _make_segment(0, steps=5000000),
            _make_segment(1, steps=3000000, status=SegmentStatus.INTERRUPTED),
        ]
        p = _make_progress(segments=segs)
        assert p.last_segment is not None
        assert p.last_segment.index == 1

    def test_time_completed_ns(self):
        segs = [
            _make_segment(0, duration_ns=10.0),
            _make_segment(1, duration_ns=5.0),
        ]
        p = _make_progress(segments=segs)
        assert p.time_completed_ns() == pytest.approx(15.0)

    def test_fraction_complete_clamped_to_one(self):
        """fraction_complete should not exceed 1.0."""
        seg = _make_segment(0, steps=12000000)
        p = _make_progress(segments=[seg], total_steps=10000000)
        assert p.fraction_complete() == 1.0


# ---------------------------------------------------------------------------
# save_progress / load_progress round-trip
# ---------------------------------------------------------------------------


class TestProgressIO:
    """Atomic save and load of progress files."""

    def test_save_creates_file(self, tmp_path):
        p = _make_progress()
        result = save_progress(tmp_path, p)
        assert result.exists()
        assert result.name == "progress.json"

    def test_round_trip(self, tmp_path):
        seg = _make_segment(0, steps=5000000)
        p = _make_progress(segments=[seg], total_steps=10000000)
        p.config_path = "/test/config.yaml"

        save_progress(tmp_path, p)
        loaded = load_progress(tmp_path)

        assert loaded is not None
        assert loaded.config_path == "/test/config.yaml"
        assert loaded.total_steps_requested == 10000000
        assert len(loaded.segments) == 1
        assert loaded.segments[0].steps_completed == 5000000
        assert loaded.segments[0].status == SegmentStatus.COMPLETED

    def test_load_returns_none_when_missing(self, tmp_path):
        assert load_progress(tmp_path) is None

    def test_load_returns_none_on_corrupt_json(self, tmp_path):
        (tmp_path / "progress.json").write_text("not valid json{{{")
        assert load_progress(tmp_path) is None

    def test_save_creates_parent_dirs(self, tmp_path):
        nested = tmp_path / "a" / "b" / "c"
        p = _make_progress()
        result = save_progress(nested, p)
        assert result.exists()

    def test_atomic_write_no_tmp_left(self, tmp_path):
        """After save, no .tmp file should remain."""
        save_progress(tmp_path, _make_progress())
        tmp_files = list(tmp_path.glob("*.tmp"))
        assert len(tmp_files) == 0


# ---------------------------------------------------------------------------
# scan_filesystem
# ---------------------------------------------------------------------------


class TestScanFilesystem:
    """Reconstruct progress from production_N/ directories."""

    def test_empty_directory(self, tmp_path):
        p = scan_filesystem(tmp_path)
        assert len(p.segments) == 0
        assert p.status == SimulationStatus.NOT_STARTED

    def test_one_completed_segment(self, tmp_path):
        _write_completed_segment_on_disk(tmp_path, 0, duration_ns=10.0, num_samples=125)
        p = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert len(p.segments) == 1
        assert p.segments[0].status == SegmentStatus.COMPLETED
        assert p.segments[0].steps_completed == 5000000  # 10ns at 2fs timestep

    def test_multiple_completed_segments(self, tmp_path):
        _write_completed_segment_on_disk(tmp_path, 0, duration_ns=10.0)
        _write_completed_segment_on_disk(tmp_path, 1, duration_ns=10.0)
        p = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert len(p.segments) == 2
        assert p.total_steps_completed == 10000000

    def test_interrupted_segment(self, tmp_path):
        _write_interrupted_segment_on_disk(
            tmp_path, 0, steps_completed=3000000, total_steps=5000000
        )
        p = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert len(p.segments) == 1
        assert p.segments[0].status == SegmentStatus.INTERRUPTED
        assert p.segments[0].steps_completed == 3000000
        assert p.status == SimulationStatus.INTERRUPTED

    def test_mixed_completed_and_interrupted(self, tmp_path):
        _write_completed_segment_on_disk(tmp_path, 0, duration_ns=10.0)
        _write_interrupted_segment_on_disk(tmp_path, 1, steps_completed=3000000)
        p = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert len(p.segments) == 2
        assert p.segments[0].status == SegmentStatus.COMPLETED
        assert p.segments[1].status == SegmentStatus.INTERRUPTED

    def test_empty_segment_dir_treated_as_failed(self, tmp_path):
        """A production_N/ dir with no state.xml or INTERRUPTED → FAILED."""
        (tmp_path / "production_0").mkdir()
        p = scan_filesystem(tmp_path)
        assert len(p.segments) == 1
        assert p.segments[0].status == SegmentStatus.FAILED

    def test_ignores_non_production_dirs(self, tmp_path):
        """Directories not matching production_N/ should be ignored."""
        (tmp_path / "equilibration").mkdir()
        (tmp_path / "production_0_resume").mkdir()
        (tmp_path / "random_dir").mkdir()
        p = scan_filesystem(tmp_path)
        assert len(p.segments) == 0

    def test_segments_sorted_by_index(self, tmp_path):
        """Segments should be sorted by index regardless of filesystem order."""
        _write_completed_segment_on_disk(tmp_path, 2, duration_ns=5.0)
        _write_completed_segment_on_disk(tmp_path, 0, duration_ns=5.0)
        _write_completed_segment_on_disk(tmp_path, 1, duration_ns=5.0)
        p = scan_filesystem(tmp_path, timestep_fs=2.0)
        indices = [s.index for s in p.segments]
        assert indices == [0, 1, 2]


# ---------------------------------------------------------------------------
# get_next_segment_info
# ---------------------------------------------------------------------------


class TestGetNextSegmentInfo:
    """Determine what work the next segment should run."""

    def test_no_segments_returns_full_run(self):
        p = _make_progress(total_steps=10000000, total_samples=250)
        info = get_next_segment_info(p, total_steps=10000000, total_samples=250)
        assert info is not None
        assert info["segment_index"] == 0
        assert info["steps_to_run"] == 10000000
        assert info["samples_to_write"] == 250

    def test_half_done_returns_remainder(self):
        seg = _make_segment(0, steps=5000000)
        p = _make_progress(segments=[seg], total_steps=10000000, total_samples=250)
        info = get_next_segment_info(p, total_steps=10000000, total_samples=250)
        assert info is not None
        assert info["segment_index"] == 1
        assert info["steps_to_run"] == 5000000
        assert info["samples_to_write"] == 125  # Half the samples

    def test_complete_returns_none(self):
        segs = [_make_segment(0, steps=5000000), _make_segment(1, steps=5000000)]
        p = _make_progress(segments=segs, total_steps=10000000)
        info = get_next_segment_info(p, total_steps=10000000, total_samples=250)
        assert info is None

    def test_overshoot_returns_none(self):
        seg = _make_segment(0, steps=12000000)
        p = _make_progress(segments=[seg], total_steps=10000000)
        info = get_next_segment_info(p, total_steps=10000000, total_samples=250)
        assert info is None


# ---------------------------------------------------------------------------
# validate_progress
# ---------------------------------------------------------------------------


class TestValidateProgress:
    """Reconcile progress file against filesystem state."""

    def test_consistent_state_unchanged(self, tmp_path):
        """When progress file and filesystem agree, nothing changes."""
        _write_completed_segment_on_disk(tmp_path, 0, duration_ns=10.0)
        seg = _make_segment(0, steps=5000000)
        p = _make_progress(segments=[seg], total_steps=10000000)
        validated = validate_progress(tmp_path, p, timestep_fs=2.0)
        assert len(validated.segments) == 1
        assert validated.segments[0].status == SegmentStatus.COMPLETED

    def test_filesystem_adds_missing_segment(self, tmp_path):
        """Segment on disk but not in progress file should be added."""
        _write_completed_segment_on_disk(tmp_path, 0, duration_ns=10.0)
        _write_completed_segment_on_disk(tmp_path, 1, duration_ns=10.0)
        # Progress file only knows about segment 0
        seg = _make_segment(0, steps=5000000)
        p = _make_progress(segments=[seg], total_steps=10000000)
        validated = validate_progress(tmp_path, p, timestep_fs=2.0)
        assert len(validated.segments) == 2

    def test_filesystem_status_is_authoritative(self, tmp_path):
        """If filesystem says interrupted but file says completed, use filesystem."""
        _write_interrupted_segment_on_disk(tmp_path, 0, steps_completed=3000000)
        # Progress file says completed
        seg = _make_segment(0, steps=5000000, status=SegmentStatus.COMPLETED)
        p = _make_progress(segments=[seg], total_steps=10000000)
        validated = validate_progress(tmp_path, p, timestep_fs=2.0)
        assert validated.segments[0].status == SegmentStatus.INTERRUPTED


# ---------------------------------------------------------------------------
# load_or_scan_progress (primary entry point)
# ---------------------------------------------------------------------------


class TestLoadOrScanProgress:
    """Primary entry point for progress loading."""

    def test_creates_progress_from_empty_dir(self, tmp_path):
        p = load_or_scan_progress(
            working_dir=tmp_path,
            config_path="/test/config.yaml",
            total_steps=10000000,
            total_samples=250,
            timestep_fs=2.0,
            replicate=1,
        )
        assert p.config_path == "/test/config.yaml"
        assert p.total_steps_requested == 10000000
        assert p.status == SimulationStatus.NOT_STARTED
        assert len(p.segments) == 0

    def test_loads_existing_progress_file(self, tmp_path):
        """When progress.json exists, it should be loaded and validated."""
        seg = _make_segment(0, steps=5000000)
        original = _make_progress(segments=[seg], total_steps=10000000)
        original.config_path = "/test/config.yaml"
        save_progress(tmp_path, original)

        # Also create the matching directory on disk
        _write_completed_segment_on_disk(tmp_path, 0, duration_ns=10.0)

        p = load_or_scan_progress(
            working_dir=tmp_path,
            config_path="/test/config.yaml",
            total_steps=10000000,
            total_samples=250,
            timestep_fs=2.0,
        )
        assert p.total_steps_requested == 10000000
        assert len(p.segments) == 1
        assert p.segments[0].status == SegmentStatus.COMPLETED

    def test_scans_filesystem_when_no_progress_file(self, tmp_path):
        """Without progress.json, should reconstruct from production dirs."""
        _write_completed_segment_on_disk(tmp_path, 0, duration_ns=10.0)
        _write_completed_segment_on_disk(tmp_path, 1, duration_ns=10.0)

        p = load_or_scan_progress(
            working_dir=tmp_path,
            config_path="/test/config.yaml",
            total_steps=10000000,
            total_samples=250,
            timestep_fs=2.0,
        )
        assert len(p.segments) == 2
        assert p.total_steps_completed == 10000000

    def test_complete_simulation_marked_as_such(self, tmp_path):
        _write_completed_segment_on_disk(tmp_path, 0, duration_ns=20.0, num_samples=250)

        p = load_or_scan_progress(
            working_dir=tmp_path,
            config_path="/test/config.yaml",
            total_steps=10000000,
            total_samples=250,
            timestep_fs=2.0,
        )
        assert p.is_complete is True
        assert p.status == SimulationStatus.COMPLETED


# ---------------------------------------------------------------------------
# Unit conversion helpers
# ---------------------------------------------------------------------------


class TestUnitConversion:
    """Internal unit conversion functions."""

    def test_convert_to_ns_nanoseconds(self):
        from polyzymd.simulation.progress import _convert_to_ns

        assert _convert_to_ns(10.0, "nanosecond") == pytest.approx(10.0)
        assert _convert_to_ns(10.0, "nanoseconds") == pytest.approx(10.0)

    def test_convert_to_ns_picoseconds(self):
        from polyzymd.simulation.progress import _convert_to_ns

        assert _convert_to_ns(1000.0, "picosecond") == pytest.approx(1.0)

    def test_convert_to_ns_femtoseconds(self):
        from polyzymd.simulation.progress import _convert_to_ns

        assert _convert_to_ns(1e6, "femtosecond") == pytest.approx(1.0)

    def test_convert_to_fs_femtoseconds(self):
        from polyzymd.simulation.progress import _convert_to_fs

        assert _convert_to_fs(2.0, "femtosecond") == pytest.approx(2.0)

    def test_convert_to_fs_picoseconds(self):
        from polyzymd.simulation.progress import _convert_to_fs

        assert _convert_to_fs(0.002, "picosecond") == pytest.approx(2.0)
