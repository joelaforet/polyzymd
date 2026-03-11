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
    EquilibrationStageRecord,
    SegmentRecord,
    SegmentStatus,
    SimulationProgress,
    SimulationStatus,
    _estimate_steps_from_csv,
    get_next_segment_info,
    load_or_scan_progress,
    load_progress,
    save_progress,
    scan_equilibration_stages,
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
        assert info["report_interval"] == 40000  # 10_000_000 // 250

    def test_half_done_returns_remainder(self):
        seg = _make_segment(0, steps=5000000)
        p = _make_progress(segments=[seg], total_steps=10000000, total_samples=250)
        info = get_next_segment_info(p, total_steps=10000000, total_samples=250)
        assert info is not None
        assert info["segment_index"] == 1
        assert info["steps_to_run"] == 5000000
        assert info["samples_to_write"] == 125  # Half the samples
        assert info["report_interval"] == 40000  # Same interval as full run

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

    def test_report_interval_is_uniform_across_segments(self):
        """Report interval should be identical regardless of how many steps remain."""
        total_steps = 750000
        total_samples = 75
        expected_interval = 10000  # 750_000 // 75

        # First segment: full run
        p0 = _make_progress(total_steps=total_steps, total_samples=total_samples)
        info0 = get_next_segment_info(p0, total_steps, total_samples)
        assert info0 is not None
        assert info0["report_interval"] == expected_interval

        # After segment 0 ran 200_000 steps
        seg0 = _make_segment(0, steps=200000)
        p1 = _make_progress(segments=[seg0], total_steps=total_steps, total_samples=total_samples)
        info1 = get_next_segment_info(p1, total_steps, total_samples)
        assert info1 is not None
        assert info1["report_interval"] == expected_interval
        assert info1["samples_to_write"] == 55  # 550_000 // 10_000

        # After two segments: 200_000 + 300_000 = 500_000
        seg1 = _make_segment(1, steps=300000)
        p2 = _make_progress(
            segments=[seg0, seg1], total_steps=total_steps, total_samples=total_samples
        )
        info2 = get_next_segment_info(p2, total_steps, total_samples)
        assert info2 is not None
        assert info2["report_interval"] == expected_interval
        assert info2["samples_to_write"] == 25  # 250_000 // 10_000


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


# ---------------------------------------------------------------------------
# Helpers — equilibration stages
# ---------------------------------------------------------------------------


def _write_completed_eq_stage_on_disk(
    working_dir: Path,
    stage_idx: int,
    stage_name: str,
) -> None:
    """Write minimal files simulating a completed equilibration stage."""
    dir_name = f"equilibration_{stage_idx}_{stage_name}"
    stage_dir = working_dir / dir_name
    stage_dir.mkdir(parents=True, exist_ok=True)
    # Checkpoint existence is the completion indicator
    (stage_dir / f"{dir_name}_checkpoint.chk").write_bytes(b"mock-chk")
    (stage_dir / f"{dir_name}_trajectory.dcd").write_bytes(b"mock-dcd")
    (stage_dir / f"{dir_name}_topology.pdb").write_text("<mock-pdb/>")
    (stage_dir / f"{dir_name}_state_data.csv").write_text("step,time\n1,0.001\n")


def _write_partial_eq_stage_on_disk(
    working_dir: Path,
    stage_idx: int,
    stage_name: str,
) -> None:
    """Write an equilibration stage dir WITHOUT a checkpoint (partial/failed)."""
    dir_name = f"equilibration_{stage_idx}_{stage_name}"
    stage_dir = working_dir / dir_name
    stage_dir.mkdir(parents=True, exist_ok=True)
    # Some files exist but NO checkpoint → stage didn't finish
    (stage_dir / f"{dir_name}_trajectory.dcd").write_bytes(b"mock-dcd")


# ---------------------------------------------------------------------------
# EquilibrationStageRecord model
# ---------------------------------------------------------------------------


class TestEquilibrationStageRecord:
    """EquilibrationStageRecord Pydantic model."""

    def test_defaults(self):
        rec = EquilibrationStageRecord(index=0)
        assert rec.name == ""
        assert rec.status == SegmentStatus.COMPLETED
        assert rec.duration_ns == 0.0
        assert rec.ensemble == "NVT"
        assert rec.finished_at is None

    def test_custom_values(self):
        rec = EquilibrationStageRecord(
            index=2,
            name="density_equilibration",
            status=SegmentStatus.COMPLETED,
            duration_ns=1.0,
            ensemble="NPT",
        )
        assert rec.index == 2
        assert rec.name == "density_equilibration"
        assert rec.ensemble == "NPT"
        assert rec.duration_ns == 1.0


# ---------------------------------------------------------------------------
# SimulationProgress — equilibration properties
# ---------------------------------------------------------------------------


class TestEquilibrationProgressProperties:
    """Equilibration-related properties on SimulationProgress."""

    def test_no_eq_stages_not_complete(self):
        p = _make_progress()
        assert p.equilibration_complete is False
        assert p.num_eq_stages_completed == 0

    def test_all_eq_stages_completed(self):
        stages = [
            EquilibrationStageRecord(index=0, name="heating", status=SegmentStatus.COMPLETED),
            EquilibrationStageRecord(index=1, name="npt", status=SegmentStatus.COMPLETED),
        ]
        p = _make_progress()
        p.equilibration_stages = stages
        assert p.equilibration_complete is True
        assert p.num_eq_stages_completed == 2

    def test_partial_eq_stages_not_complete(self):
        stages = [
            EquilibrationStageRecord(index=0, name="heating", status=SegmentStatus.COMPLETED),
            EquilibrationStageRecord(index=1, name="npt", status=SegmentStatus.FAILED),
        ]
        p = _make_progress()
        p.equilibration_stages = stages
        assert p.equilibration_complete is False
        assert p.num_eq_stages_completed == 1


# ---------------------------------------------------------------------------
# scan_equilibration_stages
# ---------------------------------------------------------------------------


class TestScanEquilibrationStages:
    """Filesystem scanning for equilibration stage directories."""

    def test_empty_directory(self, tmp_path):
        records = scan_equilibration_stages(tmp_path)
        assert records == []

    def test_nonexistent_directory(self, tmp_path):
        records = scan_equilibration_stages(tmp_path / "nonexistent")
        assert records == []

    def test_completed_stages(self, tmp_path):
        _write_completed_eq_stage_on_disk(tmp_path, 0, "heating")
        _write_completed_eq_stage_on_disk(tmp_path, 1, "npt_equilibration")
        records = scan_equilibration_stages(tmp_path)
        assert len(records) == 2
        assert records[0].index == 0
        assert records[0].name == "heating"
        assert records[0].status == SegmentStatus.COMPLETED
        assert records[1].index == 1
        assert records[1].name == "npt_equilibration"
        assert records[1].status == SegmentStatus.COMPLETED

    def test_partial_stage_marked_failed(self, tmp_path):
        _write_partial_eq_stage_on_disk(tmp_path, 0, "heating")
        records = scan_equilibration_stages(tmp_path)
        assert len(records) == 1
        assert records[0].status == SegmentStatus.FAILED

    def test_mixed_completed_and_partial(self, tmp_path):
        _write_completed_eq_stage_on_disk(tmp_path, 0, "heating")
        _write_partial_eq_stage_on_disk(tmp_path, 1, "npt_equilibration")
        records = scan_equilibration_stages(tmp_path)
        assert len(records) == 2
        assert records[0].status == SegmentStatus.COMPLETED
        assert records[1].status == SegmentStatus.FAILED

    def test_stages_sorted_by_index(self, tmp_path):
        _write_completed_eq_stage_on_disk(tmp_path, 2, "production_prep")
        _write_completed_eq_stage_on_disk(tmp_path, 0, "heating")
        _write_completed_eq_stage_on_disk(tmp_path, 1, "npt_equilibration")
        records = scan_equilibration_stages(tmp_path)
        indices = [r.index for r in records]
        assert indices == [0, 1, 2]

    def test_ignores_non_eq_directories(self, tmp_path):
        """Directories not matching equilibration_N_name/ should be ignored."""
        (tmp_path / "production_0").mkdir()
        (tmp_path / "equilibration").mkdir()  # No index/name
        (tmp_path / "random_dir").mkdir()
        records = scan_equilibration_stages(tmp_path)
        assert records == []


# ---------------------------------------------------------------------------
# scan_filesystem — equilibration stages integration
# ---------------------------------------------------------------------------


class TestScanFilesystemEquilibration:
    """scan_filesystem picks up equilibration stages."""

    def test_eq_stages_populated_in_scan(self, tmp_path):
        _write_completed_eq_stage_on_disk(tmp_path, 0, "heating")
        _write_completed_eq_stage_on_disk(tmp_path, 1, "npt_equilibration")
        p = scan_filesystem(tmp_path)
        assert len(p.equilibration_stages) == 2
        assert p.num_eq_stages_completed == 2

    def test_eq_stages_and_production_segments(self, tmp_path):
        _write_completed_eq_stage_on_disk(tmp_path, 0, "heating")
        _write_completed_segment_on_disk(tmp_path, 0, duration_ns=10.0)
        p = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert len(p.equilibration_stages) == 1
        assert len(p.segments) == 1


# ---------------------------------------------------------------------------
# Round-trip serialization with equilibration stages
# ---------------------------------------------------------------------------


class TestEquilibrationRoundTrip:
    """Save/load progress with equilibration stage records."""

    def test_round_trip_with_eq_stages(self, tmp_path):
        stages = [
            EquilibrationStageRecord(index=0, name="heating", duration_ns=0.5, ensemble="NVT"),
            EquilibrationStageRecord(index=1, name="npt_eq", duration_ns=1.0, ensemble="NPT"),
        ]
        p = _make_progress()
        p.equilibration_stages = stages

        save_progress(tmp_path, p)
        loaded = load_progress(tmp_path)

        assert loaded is not None
        assert len(loaded.equilibration_stages) == 2
        assert loaded.equilibration_stages[0].name == "heating"
        assert loaded.equilibration_stages[0].ensemble == "NVT"
        assert loaded.equilibration_stages[1].name == "npt_eq"
        assert loaded.equilibration_stages[1].ensemble == "NPT"
        assert loaded.equilibration_complete is True


# ---------------------------------------------------------------------------
# _estimate_steps_from_csv
# ---------------------------------------------------------------------------


class TestEstimateStepsFromCsv:
    """_estimate_steps_from_csv reads the last data line of a CSV."""

    def test_typical_csv(self, tmp_path):
        csv = tmp_path / "state_data.csv"
        csv.write_text(
            '#"Step","Time (ps)","Potential Energy (kJ/mole)"\n'
            "40000,80.0,-123456.0\n"
            "80000,160.0,-123500.0\n"
            "120000,240.0,-123550.0\n"
        )
        assert _estimate_steps_from_csv(csv) == 120000

    def test_quoted_columns(self, tmp_path):
        csv = tmp_path / "state_data.csv"
        csv.write_text(
            '#"Step","Time (ps)"\n'
            '"40000","80.0"\n'
            '"80000","160.0"\n'
        )
        assert _estimate_steps_from_csv(csv) == 80000

    def test_header_only_returns_zero(self, tmp_path):
        csv = tmp_path / "state_data.csv"
        csv.write_text('#"Step","Time (ps)"\n')
        assert _estimate_steps_from_csv(csv) == 0

    def test_empty_file_returns_zero(self, tmp_path):
        csv = tmp_path / "state_data.csv"
        csv.write_text("")
        assert _estimate_steps_from_csv(csv) == 0

    def test_missing_file_returns_zero(self, tmp_path):
        csv = tmp_path / "nonexistent.csv"
        assert _estimate_steps_from_csv(csv) == 0

    def test_trailing_empty_lines_skipped(self, tmp_path):
        csv = tmp_path / "state_data.csv"
        csv.write_text(
            '#"Step","Time (ps)"\n'
            "40000,80.0\n"
            "80000,160.0\n"
            "\n"
            "\n"
        )
        assert _estimate_steps_from_csv(csv) == 80000

    def test_single_data_line(self, tmp_path):
        csv = tmp_path / "state_data.csv"
        csv.write_text(
            '#"Step","Time (ps)"\n'
            "50000,100.0\n"
        )
        assert _estimate_steps_from_csv(csv) == 50000


# ---------------------------------------------------------------------------
# Checkpoint-only segment classification (hard-kill recovery)
# ---------------------------------------------------------------------------


def _write_checkpoint_only_segment(
    working_dir: Path,
    seg_idx: int,
    with_csv: bool = False,
    csv_steps: int = 80000,
) -> None:
    """Write files simulating a hard-killed segment (checkpoint but no state.xml)."""
    seg_dir = working_dir / f"production_{seg_idx}"
    seg_dir.mkdir(parents=True, exist_ok=True)
    # Checkpoint exists (from periodic CheckpointReporter)
    (seg_dir / f"production_{seg_idx}_checkpoint.chk").write_bytes(b"\x00" * 32)
    # system.xml exists (saved at segment start)
    (seg_dir / f"production_{seg_idx}_system.xml").write_text("<mock-system/>")
    if with_csv:
        (seg_dir / f"production_{seg_idx}_state_data.csv").write_text(
            '#"Step","Time (ps)","PE (kJ/mole)"\n'
            f"{csv_steps},{csv_steps * 0.002},-100000.0\n"
        )


class TestCheckpointOnlySegment:
    """Segments with checkpoint but no state.xml should be INTERRUPTED (not FAILED)."""

    def test_checkpoint_only_is_interrupted(self, tmp_path):
        _write_checkpoint_only_segment(tmp_path, 0)
        p = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert len(p.segments) == 1
        assert p.segments[0].status == SegmentStatus.INTERRUPTED
        assert p.segments[0].steps_completed == 0  # No CSV to estimate from

    def test_checkpoint_with_csv_estimates_steps(self, tmp_path):
        _write_checkpoint_only_segment(tmp_path, 0, with_csv=True, csv_steps=120000)
        p = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert len(p.segments) == 1
        assert p.segments[0].status == SegmentStatus.INTERRUPTED
        assert p.segments[0].steps_completed == 120000

    def test_checkpoint_only_counted_in_total_steps(self, tmp_path):
        """Checkpoint-only segments should count toward total_steps_completed."""
        _write_completed_segment_on_disk(tmp_path, 0, duration_ns=10.0)
        _write_checkpoint_only_segment(tmp_path, 1, with_csv=True, csv_steps=200000)
        p = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert p.segments[0].status == SegmentStatus.COMPLETED
        assert p.segments[1].status == SegmentStatus.INTERRUPTED
        # Segment 0: 10ns at 2fs = 5_000_000 steps. Segment 1: ~200_000 steps
        assert p.total_steps_completed == 5200000

    def test_no_checkpoint_no_state_is_failed(self, tmp_path):
        """Segment with directory but no checkpoint and no state.xml → FAILED."""
        seg_dir = tmp_path / "production_0"
        seg_dir.mkdir()
        # Only system.xml, no checkpoint
        (seg_dir / "production_0_system.xml").write_text("<mock/>")
        p = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert len(p.segments) == 1
        assert p.segments[0].status == SegmentStatus.FAILED

    def test_checkpoint_only_next_segment_increments(self, tmp_path):
        """After a checkpoint-only (INTERRUPTED) segment, next_segment_index increments."""
        _write_checkpoint_only_segment(tmp_path, 0, with_csv=True, csv_steps=50000)
        p = scan_filesystem(tmp_path, timestep_fs=2.0)
        progress = _make_progress(segments=p.segments, total_steps=10000000)
        assert progress.next_segment_index == 1


# ---------------------------------------------------------------------------
# FAILED segment cleanup in progress
# ---------------------------------------------------------------------------


class TestFailedSegmentHandling:
    """FAILED segments contribute 0 steps and get cleaned up."""

    def test_failed_segment_zero_steps(self):
        """FAILED segments should contribute 0 to total_steps_completed."""
        segs = [
            _make_segment(0, steps=0, status=SegmentStatus.FAILED),
        ]
        p = _make_progress(segments=segs, total_steps=10000000)
        assert p.total_steps_completed == 0

    def test_failed_after_completed(self):
        """A FAILED segment after a COMPLETED one doesn't affect completed count."""
        segs = [
            _make_segment(0, steps=5000000, status=SegmentStatus.COMPLETED),
            _make_segment(1, steps=0, status=SegmentStatus.FAILED),
        ]
        p = _make_progress(segments=segs, total_steps=10000000)
        assert p.total_steps_completed == 5000000
        assert p.next_segment_index == 2

    def test_removing_failed_segments_resets_index(self):
        """After removing FAILED segments, next_segment_index recalculates correctly."""
        segs = [
            _make_segment(0, steps=0, status=SegmentStatus.FAILED),
        ]
        p = _make_progress(segments=segs, total_steps=10000000)
        # Remove FAILED segments (simulating what main.py does)
        p.segments = [s for s in p.segments if s.status != SegmentStatus.FAILED]
        assert p.next_segment_index == 0  # Can retry from scratch
