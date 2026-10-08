"""Tests that production segments resume from the state of the last written frame.

Reproduces two restart defects seen on Blanca in September 2026:

- B, re-simulated boundary step: after a hard kill the next segment resumed
  from a wall-clock ``restart_state.xml`` older than the last trajectory
  frame and wrote that frame's step a second time.
- C, lost frame: a trajectory write raised ``OSError`` (stale file handle),
  the crash handler saved ``interrupted_state.xml`` at that step, and the
  next segment started one frame later, so the step never had a frame.

Also covers the per-segment frame interval check (a chain must keep one
interval) and the overlap and gap bookkeeping in ``progress.json``.
"""

import json
import logging
import re
import struct
from pathlib import Path

import pytest

from polyzymd.simulation.progress import (
    ReportIntervalChangeError,
    SegmentRecord,
    SegmentStatus,
    SimulationProgress,
    check_report_interval_unchanged,
    load_or_scan_progress,
    load_progress,
    save_progress,
)
from polyzymd.simulation.report_state import (
    first_report_step,
    last_reported_frame,
    read_state_step,
    segment_frames,
    write_text_atomic,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_dcd_header(path: Path, frames: int, interval: int) -> None:
    """Write the 100-byte header OpenMM's DCDReporter starts every trajectory with.

    DCDReporter always passes ``firstStep=reportInterval``, whatever step the
    simulation starts at, so the header steps are relative to the segment start.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    first_step = interval
    last_step = first_step + (frames - 1) * interval
    header = struct.pack(
        "<i4c9if",
        84,
        b"C",
        b"O",
        b"R",
        b"D",
        frames,
        first_step,
        interval,
        last_step,
        0,
        0,
        0,
        0,
        0,
        0.0409,
    )
    header += struct.pack("<13i", 1, 0, 0, 0, 0, 0, 0, 0, 0, 24, 84, 164, 2)
    path.write_bytes(header)


def _write_state(path: Path, step: int | None) -> None:
    """Write a minimal complete State XML carrying *step* as ``stepCount``."""
    path.parent.mkdir(parents=True, exist_ok=True)
    attr = f' stepCount="{step}"' if step is not None else ""
    path.write_text(f'<?xml version="1.0" ?>\n<State{attr} time="0" type="State" version="1"/>\n')


def _prev_segment(tmp_path: Path, index: int, *, frames: int, first: int, interval: int) -> Path:
    """Create ``production_<index>`` with a system XML, a DCD header and a state-data CSV.

    *first* is the step of the first frame; the CSV holds one row per frame.
    """
    seg_dir = tmp_path / f"production_{index}"
    seg_dir.mkdir(parents=True, exist_ok=True)
    (seg_dir / f"production_{index}_system.xml").write_text("<System/>")
    _write_dcd_header(seg_dir / f"production_{index}_trajectory.dcd", frames, interval)
    rows = "".join(f"{first + i * interval},0.0\n" for i in range(frames))
    (seg_dir / f"production_{index}_state_data.csv").write_text('#"Step","Time (ps)"\n' + rows)
    return seg_dir


def _manager(working_dir: Path, prev_segment: int):
    from polyzymd.simulation.continuation import ContinuationManager

    mgr = ContinuationManager.__new__(ContinuationManager)
    mgr._working_dir = Path(working_dir)
    mgr._prev_segment = prev_segment
    return mgr


def _csv_steps(path: Path) -> list[int]:
    rows = [line for line in path.read_text().splitlines() if line and not line.startswith("#")]
    return [int(float(row.split(",")[0])) for row in rows]


# ---------------------------------------------------------------------------
# File-header readers
# ---------------------------------------------------------------------------


class TestHeaderReaders:
    def test_read_state_step(self, tmp_path):
        _write_state(tmp_path / "s.xml", 250598907)
        assert read_state_step(tmp_path / "s.xml") == 250598907
        _write_state(tmp_path / "old.xml", None)
        assert read_state_step(tmp_path / "old.xml") is None
        assert read_state_step(tmp_path / "missing.xml") is None

    def test_segment_frames_uses_absolute_steps_of_a_later_segment(self, tmp_path):
        """The DCD header counts steps from the segment start; the CSV holds the real steps."""
        _prev_segment(tmp_path, 9, frames=154, first=220000000, interval=200000)
        frame = segment_frames(tmp_path, 9)
        assert frame.last_step == 250600000
        assert frame.report_interval == 200000
        assert frame.frames == 154

    def test_segment_frames_counts_a_frame_missing_from_the_csv(self, tmp_path):
        """A kill between the DCD write and the CSV row of one step leaves one more frame."""
        seg_dir = _prev_segment(tmp_path, 1, frames=3, first=1200, interval=200)
        _write_dcd_header(seg_dir / "production_1_trajectory.dcd", 4, 200)
        frame = segment_frames(tmp_path, 1)
        assert (frame.last_step, frame.frames) == (1800, 4)

    def test_segment_frames_falls_back_to_csv(self, tmp_path):
        seg_dir = tmp_path / "production_0"
        seg_dir.mkdir()
        (seg_dir / "production_0_state_data.csv").write_text(
            '#"Step","Time (ps)"\n200,0.4\n400,0.8\n600,1.2\n'
        )
        frame = segment_frames(tmp_path, 0)
        assert (frame.last_step, frame.report_interval, frame.frames) == (600, 200, 3)

    def test_last_reported_frame_passes_over_empty_segments(self, tmp_path):
        _prev_segment(tmp_path, 0, frames=5, first=200, interval=200)
        seg1 = tmp_path / "production_1"
        seg1.mkdir()
        (seg1 / "production_1_trajectory.dcd").write_bytes(b"")
        frame = last_reported_frame(tmp_path, 2)
        assert frame.segment_index == 0
        assert frame.last_step == 1000

    def test_first_report_step(self):
        assert first_report_step(250598907, 200000) == 250600000
        assert first_report_step(210400000, 200000) == 210600000

    def test_write_text_atomic_leaves_no_temporary(self, tmp_path):
        target = tmp_path / "restart_state.xml"
        write_text_atomic(target, "<State/>")
        write_text_atomic(target, "<State stepCount='2'/>")
        assert target.read_text() == "<State stepCount='2'/>"
        assert [p.name for p in tmp_path.iterdir()] == ["restart_state.xml"]


# ---------------------------------------------------------------------------
# Choosing the state to resume from
# ---------------------------------------------------------------------------


class TestResumeStateSelection:
    """``_find_portable_state`` picks the state that continues the frames."""

    def test_lost_frame_state_is_not_chosen(self, tmp_path):
        """Pattern C: interrupted_state.xml at the step whose frame failed."""
        seg_dir = _prev_segment(tmp_path, 2, frames=122, first=186000000, interval=200000)
        _write_state(seg_dir / "interrupted_state.xml", 210400000)
        (seg_dir / "interrupted_system.xml").write_text("<System/>")
        _write_state(seg_dir / "restart_state.xml", 210362569)

        state, _ = _manager(tmp_path, 2)._find_portable_state()
        assert state.name == "restart_state.xml"

    def test_latest_state_before_the_next_frame_wins(self, tmp_path):
        seg_dir = _prev_segment(tmp_path, 0, frames=10, first=200, interval=200)
        _write_state(seg_dir / "interrupted_state.xml", 2050)
        (seg_dir / "interrupted_system.xml").write_text("<System/>")
        _write_state(seg_dir / "restart_state.xml", 2150)

        state, _ = _manager(tmp_path, 0)._find_portable_state()
        assert state.name == "restart_state.xml"

    def test_completed_segment_keeps_its_final_state(self, tmp_path):
        seg_dir = _prev_segment(tmp_path, 0, frames=10, first=200, interval=200)
        _write_state(seg_dir / "production_0_state.xml", 2000)
        _write_state(seg_dir / "restart_state.xml", 2000)

        state, _ = _manager(tmp_path, 0)._find_portable_state()
        assert state.name == "production_0_state.xml"

    def test_all_states_past_a_lost_frame_choose_the_earliest(self, tmp_path):
        seg_dir = _prev_segment(tmp_path, 0, frames=10, first=200, interval=200)
        _write_state(seg_dir / "interrupted_state.xml", 2400)
        (seg_dir / "interrupted_system.xml").write_text("<System/>")
        _write_state(seg_dir / "restart_state.xml", 2200)

        state, _ = _manager(tmp_path, 0)._find_portable_state()
        assert state.name == "restart_state.xml"

    def test_states_without_step_counts_keep_the_file_order(self, tmp_path):
        seg_dir = _prev_segment(tmp_path, 0, frames=10, first=200, interval=200)
        _write_state(seg_dir / "interrupted_state.xml", None)
        (seg_dir / "interrupted_system.xml").write_text("<System/>")
        _write_state(seg_dir / "restart_state.xml", 2150)

        state, _ = _manager(tmp_path, 0)._find_portable_state()
        assert state.name == "interrupted_state.xml"

    def test_truncated_restart_state_is_still_skipped(self, tmp_path):
        seg_dir = _prev_segment(tmp_path, 0, frames=10, first=200, interval=200)
        _write_state(seg_dir / "interrupted_state.xml", 1900)
        (seg_dir / "interrupted_system.xml").write_text("<System/>")
        (seg_dir / "restart_state.xml").write_text('<State stepCount="2100">')

        state, _ = _manager(tmp_path, 0)._find_portable_state()
        assert state.name == "interrupted_state.xml"


# ---------------------------------------------------------------------------
# Frame interval kept constant along a chain
# ---------------------------------------------------------------------------


class TestReportIntervalCheck:
    def _progress(self, **fields) -> SimulationProgress:
        return SimulationProgress(
            segments=[SegmentRecord(index=0, status=SegmentStatus.INTERRUPTED, **fields)]
        )

    def test_recorded_interval_change_is_refused(self, tmp_path):
        progress = self._progress(report_interval=20000)
        with pytest.raises(ReportIntervalChangeError, match="every 20000 steps"):
            check_report_interval_unchanged(progress, tmp_path, 1, 200000)

    def test_override_allows_the_change(self, tmp_path, caplog):
        progress = self._progress(report_interval=20000)
        # `polyzymd status` raises this logger to ERROR for the rest of the process.
        with caplog.at_level(logging.WARNING, logger="polyzymd.simulation.progress"):
            check_report_interval_unchanged(progress, tmp_path, 1, 200000, allow_change=True)
        assert "explicitly allowed" in caplog.text

    def test_unrecorded_interval_falls_back_to_the_dcd_header(self, tmp_path):
        _prev_segment(tmp_path, 0, frames=10, first=20000, interval=20000)
        with pytest.raises(ReportIntervalChangeError):
            check_report_interval_unchanged(self._progress(), tmp_path, 1, 200000)
        check_report_interval_unchanged(self._progress(), tmp_path, 1, 20000)

    def test_no_earlier_frames_means_nothing_to_compare(self, tmp_path):
        check_report_interval_unchanged(SimulationProgress(), tmp_path, 0, 200000)


# ---------------------------------------------------------------------------
# End-to-end on the OpenMM Reference platform
# ---------------------------------------------------------------------------

REPORT_INTERVAL = 5


def _toy_system():
    """One uncharged particle in a periodic box: the smallest NPT-capable system."""
    from openmm import NonbondedForce, System, Vec3, unit
    from openmm.app import Element, Topology

    topology = Topology()
    chain = topology.addChain("A")
    residue = topology.addResidue("MOL", chain)
    topology.addAtom("C", Element.getBySymbol("C"), residue)
    system = System()
    system.addParticle(12.0 * unit.dalton)
    system.setDefaultPeriodicBoxVectors(
        Vec3(2.0, 0.0, 0.0), Vec3(0.0, 2.0, 0.0), Vec3(0.0, 0.0, 2.0)
    )
    nonbonded = NonbondedForce()
    nonbonded.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
    nonbonded.setCutoffDistance(0.9 * unit.nanometer)
    nonbonded.addParticle(0.0, 0.3, 0.1)
    system.addForce(nonbonded)
    positions = [Vec3(1.0, 1.0, 1.0)] * unit.nanometer
    return topology, system, positions


def _run_segment_zero(tmp_path: Path, total_steps: int = 20) -> None:
    from polyzymd.simulation.runner import SimulationRunner

    progress = load_or_scan_progress(
        tmp_path, total_steps=100, total_samples=100 // REPORT_INTERVAL, timestep_fs=1.0
    )
    save_progress(tmp_path, progress)
    topology, system, positions = _toy_system()
    runner = SimulationRunner(
        topology=topology,
        system=system,
        positions=positions,
        working_dir=tmp_path,
        platform="Reference",
    )
    runner.run_production(
        temperature=300.0,
        duration_ns=total_steps * 1e-6,
        num_samples=total_steps // REPORT_INTERVAL,
        timestep_fs=1.0,
        segment_index=0,
        report_interval=REPORT_INTERVAL,
        checkpoint_interval_s=1e-9,  # wall-clock restart save after every chunk
    )


def _run_segment_one(tmp_path: Path, total_steps: int = 10, *, index: int = 1):
    from polyzymd.simulation.continuation import ContinuationManager

    manager = ContinuationManager(working_dir=tmp_path, segment_index=index, platform="Reference")
    manager.load_previous_state()
    manager.run_segment(
        duration_ns=total_steps * 1e-6,
        num_samples=total_steps // REPORT_INTERVAL,
        timestep_fs=1.0,
        report_interval=REPORT_INTERVAL,
        checkpoint_interval_s=3600.0,
    )
    return manager


def _hard_kill(seg_dir: Path) -> None:
    """Remove what only a clean finish or graceful stop leaves behind."""
    for name in ("production_0_state.xml", "INTERRUPTED", "interrupted_state.xml"):
        (seg_dir / name).unlink(missing_ok=True)


@pytest.fixture
def _no_signal_state():
    from polyzymd.simulation import signals

    signals.reset()
    yield
    signals.reset()


@pytest.mark.usefixtures("_no_signal_state")
class TestRestartEndToEnd:
    def test_restart_state_tracks_every_frame(self, tmp_path):
        _run_segment_zero(tmp_path)
        seg_dir = tmp_path / "production_0"
        assert _csv_steps(seg_dir / "production_0_state_data.csv") == [5, 10, 15, 20]
        assert read_state_step(seg_dir / "restart_state.xml") == 20

        seg = load_progress(tmp_path).segments[0]
        assert seg.report_interval == REPORT_INTERVAL
        assert seg.start_step == 0
        assert seg.last_reported_step == 20
        assert seg.samples_written == 4

    def test_hard_kill_resumes_after_the_last_frame(self, tmp_path):
        """Pattern B: a hard-killed segment is continued without a repeated step."""
        _run_segment_zero(tmp_path)
        _hard_kill(tmp_path / "production_0")

        _run_segment_one(tmp_path)
        steps = _csv_steps(tmp_path / "production_1" / "production_1_state_data.csv")
        assert steps == [25, 30]
        seg1 = next(s for s in load_progress(tmp_path).segments if s.index == 1)
        assert seg1.resumed_from == "restart_state.xml"
        assert seg1.start_step == 20
        assert seg1.overlap_frames == 0
        assert seg1.gap_frames == 0

    def test_resume_before_the_last_frame_writes_no_frame_twice(self, tmp_path, caplog):
        """A state older than the last frame does not repeat that frame's steps."""
        _run_segment_zero(tmp_path)
        seg_dir = tmp_path / "production_0"
        _hard_kill(seg_dir)
        # Stand-in for the old wall-clock restart state: relabel it step 12.
        restart = seg_dir / "restart_state.xml"
        restart.write_text(re.sub(r'stepCount="\d+"', 'stepCount="12"', restart.read_text()))

        _run_segment_one(tmp_path, total_steps=10)
        assert _csv_steps(tmp_path / "production_1" / "production_1_state_data.csv") == [25, 30]
        seg1 = next(s for s in load_progress(tmp_path).segments if s.index == 1)
        assert (seg1.overlap_frames, seg1.samples_written) == (0, 2)
        assert "not writing its first 2 frame(s)" in caplog.text

    def test_hard_killed_segment_resumed_from_an_older_state_completes_the_chain(self, tmp_path):
        """The chain ends at the requested step with one frame per report step.

        The killed segment wrote its frame at step 20 but its restart state
        is still at step 15, and the scan estimated its steps as 20.
        """
        _run_segment_zero(tmp_path)
        seg_dir = tmp_path / "production_0"
        _hard_kill(seg_dir)
        restart = seg_dir / "restart_state.xml"
        restart.write_text(re.sub(r'stepCount="\d+"', 'stepCount="15"', restart.read_text()))
        progress = load_progress(tmp_path)
        seg0 = progress.segments[0]
        seg0.status = SegmentStatus.INTERRUPTED
        seg0.steps_completed = 20
        save_progress(tmp_path, progress)

        _run_segment_one(tmp_path, total_steps=progress.total_steps_requested - 20)

        frames = _csv_steps(seg_dir / "production_0_state_data.csv") + _csv_steps(
            tmp_path / "production_1" / "production_1_state_data.csv"
        )
        assert frames == list(range(5, 101, 5))
        reloaded = load_or_scan_progress(tmp_path, total_steps=100, timestep_fs=1.0)
        assert [s.steps_completed for s in reloaded.segments] == [15, 85]
        assert reloaded.total_steps_completed == 100
        params = json.loads((tmp_path / "production_1" / "production_1_parameters.json").read_text())
        assert params["__values__"]["integ_params"]["__values__"]["num_samples"] == 16

    def test_failed_frame_write_is_not_skipped(self, tmp_path, monkeypatch):
        """Pattern C: the frame whose write raised is written by the next segment."""
        from openmm.app import DCDReporter

        original = DCDReporter.report
        calls = {"n": 0}

        def failing_report(self, simulation, state):
            calls["n"] += 1
            if calls["n"] == 3:
                raise OSError(116, "Stale file handle")
            return original(self, simulation, state)

        monkeypatch.setattr(DCDReporter, "report", failing_report)
        with pytest.raises(OSError, match="Stale file handle"):
            _run_segment_zero(tmp_path)
        monkeypatch.setattr(DCDReporter, "report", original)

        seg_dir = tmp_path / "production_0"
        assert not (seg_dir / "interrupted_state.xml").exists()
        assert read_state_step(seg_dir / "restart_state.xml") == 10
        marker = dict(
            line.split("=") for line in (seg_dir / "INTERRUPTED").read_text().splitlines()
        )
        assert marker["steps_completed"] == "10"
        seg0 = load_progress(tmp_path).segments[0]
        assert seg0.status == SegmentStatus.INTERRUPTED
        assert (seg0.steps_completed, seg0.last_reported_step, seg0.samples_written) == (10, 10, 2)

        _run_segment_one(tmp_path)
        steps = _csv_steps(tmp_path / "production_1" / "production_1_state_data.csv")
        assert steps[0] == 15
        seg1 = next(s for s in load_progress(tmp_path).segments if s.index == 1)
        assert (seg1.overlap_frames, seg1.gap_frames) == (0, 0)
        params = json.loads(
            (tmp_path / "production_1" / "production_1_parameters.json").read_text()
        )
        reporter = params["__values__"]["reporter_params"]["__values__"]
        assert reporter["report_interval"] == REPORT_INTERVAL

    def test_segment_after_a_resumed_segment_continues_its_frames(self, tmp_path, caplog):
        """Frames of a segment that began at a non-zero step are read at their real steps."""
        _run_segment_zero(tmp_path)
        _run_segment_one(tmp_path)
        assert segment_frames(tmp_path, 1).last_step == 30

        _run_segment_one(tmp_path, index=2)
        assert _csv_steps(tmp_path / "production_2" / "production_2_state_data.csv") == [35, 40]
        seg2 = next(s for s in load_progress(tmp_path).segments if s.index == 2)
        assert (seg2.start_step, seg2.overlap_frames, seg2.gap_frames) == (30, 0, 0)
        assert "never written" not in caplog.text
        assert "have no frame" not in caplog.text

    def test_hard_killed_segment_is_followed_by_one_that_ends_at_the_total(self, tmp_path):
        """The next segment stops at the requested total, whatever the killed one was estimated at."""
        _run_segment_zero(tmp_path)
        _hard_kill(tmp_path / "production_0")
        # What the filesystem scan records for a hard-killed segment: an
        # estimated step count and no frame bookkeeping.
        progress = load_progress(tmp_path)
        seg0 = progress.segments[0]
        seg0.status = SegmentStatus.INTERRUPTED
        seg0.steps_completed = 15
        seg0.last_reported_step = None
        seg0.samples_written = 0
        save_progress(tmp_path, progress)

        _run_segment_one(tmp_path, total_steps=progress.total_steps_requested - 15)
        steps = _csv_steps(tmp_path / "production_1" / "production_1_state_data.csv")
        assert steps[-1] == 100
        assert load_progress(tmp_path).total_steps_completed == 100

    def test_hard_killed_segment_frames_are_read_from_its_files(self, tmp_path):
        """The next segment fills in the frame fields a hard-killed segment never wrote."""
        _run_segment_zero(tmp_path)
        _hard_kill(tmp_path / "production_0")
        progress = load_progress(tmp_path)
        seg0 = progress.segments[0]
        seg0.status = SegmentStatus.INTERRUPTED
        seg0.last_reported_step = None
        seg0.samples_written = 0
        save_progress(tmp_path, progress)

        _run_segment_one(tmp_path)
        seg0 = load_progress(tmp_path).segments[0]
        assert (seg0.last_reported_step, seg0.samples_written) == (20, 4)
