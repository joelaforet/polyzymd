"""Tests for the automatic restart and hard-kill recovery infrastructure.

Covers:
- Checkpoint fallback recovery paths in ContinuationManager._get_previous_paths()
- File validation for checkpoint vs state-based recovery
- _estimate_steps_from_csv edge cases (covered more extensively in test_progress.py)
- EQ_INTERRUPTED marker parsing in _find_interrupted_eq_stage()
- _find_completed_eq_stages with interrupted stages
- FAILED segment cleanup logic
- End-to-end restart flow for checkpoint-only segments
"""

import json
import shutil
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_normal_segment(working_dir: Path, seg_idx: int) -> None:
    """Write minimal files simulating a normally completed segment."""
    seg_dir = working_dir / f"production_{seg_idx}"
    seg_dir.mkdir(parents=True, exist_ok=True)
    (seg_dir / f"production_{seg_idx}_state.xml").write_text("<State/>")
    (seg_dir / f"production_{seg_idx}_system.xml").write_text("<System/>")
    (seg_dir / f"production_{seg_idx}_checkpoint.chk").write_bytes(b"\x00" * 16)
    params = {
        "__values__": {
            "integ_params": {
                "__values__": {
                    "num_samples": 100,
                    "time_step": {
                        "__class__": "Quantity",
                        "__values__": {"value": 2.0, "unit": "femtosecond"},
                    },
                    "total_time": {
                        "__class__": "Quantity",
                        "__values__": {"value": 10.0, "unit": "nanosecond"},
                    },
                }
            },
            "thermo_params": {"__values__": {}},
        }
    }
    (seg_dir / f"production_{seg_idx}_parameters.json").write_text(json.dumps(params))


def _write_interrupted_segment(
    working_dir: Path, seg_idx: int, steps_completed: int = 3000000
) -> None:
    """Write files for a gracefully interrupted segment (SIGUSR1/SIGTERM)."""
    seg_dir = working_dir / f"production_{seg_idx}"
    seg_dir.mkdir(parents=True, exist_ok=True)
    (seg_dir / "interrupted_state.xml").write_text("<State/>")
    (seg_dir / "interrupted_system.xml").write_text("<System/>")
    (seg_dir / "interrupted_checkpoint.chk").write_bytes(b"\x00" * 16)
    (seg_dir / "INTERRUPTED").write_text(
        f"segment_index={seg_idx}\n"
        f"steps_completed={steps_completed}\n"
        f"total_steps=5000000\n"
        f"remaining_steps={5000000 - steps_completed}\n"
    )
    params = {
        "__values__": {
            "integ_params": {
                "__values__": {
                    "num_samples": 100,
                    "time_step": {
                        "__class__": "Quantity",
                        "__values__": {"value": 2.0, "unit": "femtosecond"},
                    },
                    "total_time": {
                        "__class__": "Quantity",
                        "__values__": {"value": 10.0, "unit": "nanosecond"},
                    },
                }
            },
            "thermo_params": {"__values__": {}},
        }
    }
    (seg_dir / f"production_{seg_idx}_parameters.json").write_text(json.dumps(params))


def _write_hard_killed_segment(
    working_dir: Path,
    seg_idx: int,
    with_system_xml: bool = True,
    with_csv: bool = False,
    csv_steps: int = 80000,
) -> None:
    """Write files for a hard-killed segment (SIGKILL/OOM): checkpoint + system.xml only."""
    import os
    import time as _time

    seg_dir = working_dir / f"production_{seg_idx}"
    seg_dir.mkdir(parents=True, exist_ok=True)
    # Periodic checkpoint from CheckpointReporter
    chk_path = seg_dir / f"production_{seg_idx}_checkpoint.chk"
    chk_path.write_bytes(b"\x00" * 16)
    if with_system_xml:
        (seg_dir / f"production_{seg_idx}_system.xml").write_text("<System/>")
    if with_csv:
        # Write two data rows so _estimate_steps_from_csv computes the
        # correct per-segment delta (last_step - first_step = csv_steps).
        (seg_dir / f"production_{seg_idx}_state_data.csv").write_text(
            f'#"Step","Time (ps)","PE (kJ/mole)"\n'
            f"0,0.0,-100000.0\n"
            f"{csv_steps},{csv_steps * 0.002},-100000.0\n"
        )
    # Parameters file (continuation manager needs this)
    params = {
        "__values__": {
            "integ_params": {
                "__values__": {
                    "num_samples": 100,
                    "time_step": {
                        "__class__": "Quantity",
                        "__values__": {"value": 2.0, "unit": "femtosecond"},
                    },
                    "total_time": {
                        "__class__": "Quantity",
                        "__values__": {"value": 10.0, "unit": "nanosecond"},
                    },
                }
            },
            "thermo_params": {"__values__": {}},
        }
    }
    (seg_dir / f"production_{seg_idx}_parameters.json").write_text(json.dumps(params))

    # Backdate checkpoint mtime so the recency heuristic classifies this
    # as INTERRUPTED (hard-killed) rather than RUNNING.
    old_time = _time.time() - 1200  # 20 minutes ago
    os.utime(chk_path, (old_time, old_time))


def _write_eq_stage(
    working_dir: Path,
    stage_index: int,
    stage_name: str,
    completed: bool = True,
    interrupted: bool = False,
    steps_completed: int = 50000,
    total_steps: int = 100000,
    current_temperature: float = 300.0,
    is_temperature_ramping: bool = False,
) -> None:
    """Write equilibration stage files on disk."""
    dir_name = f"equilibration_{stage_index}_{stage_name}"
    stage_dir = working_dir / dir_name
    stage_dir.mkdir(parents=True, exist_ok=True)

    if completed or interrupted:
        # Checkpoint exists in both completed and interrupted cases
        (stage_dir / f"{dir_name}_checkpoint.chk").write_bytes(b"\x00" * 16)

    if interrupted:
        (stage_dir / "EQ_INTERRUPTED").write_text(
            f"stage_index={stage_index}\n"
            f"stage_name={stage_name}\n"
            f"steps_completed={steps_completed}\n"
            f"total_steps={total_steps}\n"
            f"current_temperature={current_temperature}\n"
            f"is_temperature_ramping={is_temperature_ramping}\n"
        )


# ---------------------------------------------------------------------------
# ContinuationManager._get_previous_paths — recovery path selection
# ---------------------------------------------------------------------------


class TestGetPreviousPaths:
    """ContinuationManager._get_previous_paths selects the right recovery path."""

    def _make_manager(self, working_dir, prev_segment):
        """Create a minimal ContinuationManager-like object for testing."""
        # We import here to test in isolation; mock out openmm if unavailable
        try:
            from polyzymd.simulation.continuation import ContinuationManager

            mgr = ContinuationManager.__new__(ContinuationManager)
            mgr._working_dir = Path(working_dir)
            mgr._prev_segment = prev_segment
            return mgr
        except ImportError:
            pytest.skip("polyzymd.simulation.continuation not importable")

    def test_normal_completion_path(self, tmp_path):
        """Normal segment: state.xml exists → no checkpoint recovery."""
        _write_normal_segment(tmp_path, 0)
        mgr = self._make_manager(tmp_path, 0)
        paths = mgr._get_previous_paths()
        assert paths["use_checkpoint"] is False
        assert paths["state"].exists()
        assert "state.xml" in str(paths["state"])

    def test_graceful_interrupt_prefers_state_xml(self, tmp_path):
        """Gracefully interrupted with interrupted_state.xml → portable state recovery."""
        _write_interrupted_segment(tmp_path, 0)
        mgr = self._make_manager(tmp_path, 0)
        paths = mgr._get_previous_paths()
        # New behaviour: interrupted_state.xml exists → use loadState (portable)
        assert paths["use_checkpoint"] is False
        assert "interrupted_state.xml" in str(paths["state"])
        assert "interrupted_system.xml" in str(paths["system"])

    def test_graceful_interrupt_legacy_chk_only(self, tmp_path):
        """Gracefully interrupted with only .chk (no state.xml) → checkpoint recovery."""
        seg_dir = tmp_path / "production_0"
        seg_dir.mkdir(parents=True, exist_ok=True)
        # Legacy format: only .chk and system.xml, no interrupted_state.xml
        (seg_dir / "interrupted_checkpoint.chk").write_bytes(b"\x00" * 16)
        (seg_dir / "interrupted_system.xml").write_text("<System/>")
        (seg_dir / "INTERRUPTED").write_text(
            "segment_index=0\nsteps_completed=3000000\n"
            "total_steps=5000000\nremaining_steps=2000000\n"
        )
        params = {
            "__values__": {
                "integ_params": {
                    "__values__": {
                        "num_samples": 100,
                        "time_step": {
                            "__class__": "Quantity",
                            "__values__": {"value": 2.0, "unit": "femtosecond"},
                        },
                        "total_time": {
                            "__class__": "Quantity",
                            "__values__": {"value": 10.0, "unit": "nanosecond"},
                        },
                    }
                },
                "thermo_params": {"__values__": {}},
            }
        }
        (seg_dir / "production_0_parameters.json").write_text(json.dumps(params))

        mgr = self._make_manager(tmp_path, 0)
        paths = mgr._get_previous_paths()
        assert paths["use_checkpoint"] is True
        assert "interrupted_checkpoint.chk" in str(paths["checkpoint"])
        assert "interrupted_system.xml" in str(paths["system"])

    def test_restart_checkpoint_recovery(self, tmp_path):
        """Wall-time restart checkpoint: restart_state.xml → portable state recovery."""
        seg_dir = tmp_path / "production_0"
        seg_dir.mkdir(parents=True, exist_ok=True)
        # No normal state.xml, no interrupted files, but restart checkpoint exists
        (seg_dir / "restart_state.xml").write_text("<State/>")
        (seg_dir / "restart_system.xml").write_text("<System/>")
        params = {
            "__values__": {
                "integ_params": {
                    "__values__": {
                        "num_samples": 100,
                        "time_step": {
                            "__class__": "Quantity",
                            "__values__": {"value": 2.0, "unit": "femtosecond"},
                        },
                        "total_time": {
                            "__class__": "Quantity",
                            "__values__": {"value": 10.0, "unit": "nanosecond"},
                        },
                    }
                },
                "thermo_params": {"__values__": {}},
            }
        }
        (seg_dir / "production_0_parameters.json").write_text(json.dumps(params))

        mgr = self._make_manager(tmp_path, 0)
        paths = mgr._get_previous_paths()
        assert paths["use_checkpoint"] is False
        assert "restart_state.xml" in str(paths["state"])
        assert "restart_system.xml" in str(paths["system"])

    def test_interrupted_state_preferred_over_restart(self, tmp_path):
        """When both interrupted_state.xml and restart_state.xml exist,
        interrupted_state.xml is preferred (more recent)."""
        seg_dir = tmp_path / "production_0"
        seg_dir.mkdir(parents=True, exist_ok=True)
        (seg_dir / "interrupted_state.xml").write_text("<State-interrupted/>")
        (seg_dir / "interrupted_system.xml").write_text("<System-interrupted/>")
        (seg_dir / "restart_state.xml").write_text("<State-restart/>")
        (seg_dir / "restart_system.xml").write_text("<System-restart/>")
        params = {
            "__values__": {
                "integ_params": {
                    "__values__": {
                        "num_samples": 100,
                        "time_step": {
                            "__class__": "Quantity",
                            "__values__": {"value": 2.0, "unit": "femtosecond"},
                        },
                        "total_time": {
                            "__class__": "Quantity",
                            "__values__": {"value": 10.0, "unit": "nanosecond"},
                        },
                    }
                },
                "thermo_params": {"__values__": {}},
            }
        }
        (seg_dir / "production_0_parameters.json").write_text(json.dumps(params))

        mgr = self._make_manager(tmp_path, 0)
        paths = mgr._get_previous_paths()
        assert paths["use_checkpoint"] is False
        assert "interrupted_state.xml" in str(paths["state"])
        assert "interrupted_system.xml" in str(paths["system"])

    def test_hard_kill_checkpoint_recovery(self, tmp_path):
        """Hard-killed: only checkpoint + system.xml → checkpoint recovery."""
        _write_hard_killed_segment(tmp_path, 0)
        mgr = self._make_manager(tmp_path, 0)
        paths = mgr._get_previous_paths()
        assert paths["use_checkpoint"] is True
        assert "production_0_system.xml" in str(paths["system"])
        assert "production_0_checkpoint.chk" in str(paths["checkpoint"])

    def test_hard_kill_no_system_xml_not_recoverable(self, tmp_path):
        """Hard-killed with no system.xml: checkpoint exists but not recoverable."""
        _write_hard_killed_segment(tmp_path, 0, with_system_xml=False)
        mgr = self._make_manager(tmp_path, 0)
        paths = mgr._get_previous_paths()
        # use_checkpoint should still be False because there's no system.xml
        # (the error branch logs an error but doesn't set use_checkpoint)
        assert paths["use_checkpoint"] is False

    def test_hard_kill_recovery_after_completed_segment(self, tmp_path):
        """Normal segment 0, then hard-killed segment 1 → recovery from seg 1."""
        _write_normal_segment(tmp_path, 0)
        _write_hard_killed_segment(tmp_path, 1)
        mgr = self._make_manager(tmp_path, 1)
        paths = mgr._get_previous_paths()
        assert paths["use_checkpoint"] is True


# ---------------------------------------------------------------------------
# File validation in load_previous_state
# ---------------------------------------------------------------------------


class TestLoadPreviousStateValidation:
    """load_previous_state validates files based on recovery mode."""

    def _make_manager(self, working_dir, prev_segment):
        try:
            from polyzymd.simulation.continuation import ContinuationManager

            mgr = ContinuationManager.__new__(ContinuationManager)
            mgr._working_dir = Path(working_dir)
            mgr._prev_segment = prev_segment
            return mgr
        except ImportError:
            pytest.skip("polyzymd.simulation.continuation not importable")

    def test_normal_segment_validates_state_xml(self, tmp_path):
        """For normal recovery, state.xml must exist."""
        seg_dir = tmp_path / "production_0"
        seg_dir.mkdir()
        # No state.xml, no checkpoint, no interrupted files
        (seg_dir / "production_0_system.xml").write_text("<System/>")
        (seg_dir / "production_0_parameters.json").write_text("{}")
        mgr = self._make_manager(tmp_path, 0)
        with pytest.raises(FileNotFoundError, match="state"):
            mgr.load_previous_state()

    def test_checkpoint_recovery_validates_checkpoint(self, tmp_path):
        """Legacy format: only interrupted_system.xml + interrupted_checkpoint.chk
        (no interrupted_state.xml) → checkpoint recovery validates .chk exists."""
        seg_dir = tmp_path / "production_0"
        seg_dir.mkdir()
        # Has interrupted_system.xml and interrupted_checkpoint.chk but no state XML
        (seg_dir / "interrupted_system.xml").write_text("<System/>")
        # No interrupted_checkpoint.chk — should fail validation
        (seg_dir / "production_0_parameters.json").write_text("{}")
        mgr = self._make_manager(tmp_path, 0)
        # Without interrupted_state.xml or interrupted_checkpoint.chk,
        # the recovery path doesn't match case 4.  Falls through to
        # no-match → state_path points to non-existent production_0_state.xml.
        with pytest.raises(FileNotFoundError, match="state"):
            mgr.load_previous_state()

    def test_legacy_chk_recovery_validates_checkpoint_exists(self, tmp_path):
        """Legacy format: interrupted_system.xml + interrupted_checkpoint.chk
        present → checkpoint recovery mode."""
        seg_dir = tmp_path / "production_0"
        seg_dir.mkdir()
        (seg_dir / "interrupted_system.xml").write_text("<System/>")
        (seg_dir / "interrupted_checkpoint.chk").write_bytes(b"\x00")
        (seg_dir / "production_0_parameters.json").write_text("{}")
        mgr = self._make_manager(tmp_path, 0)
        paths = mgr._get_previous_paths()
        # Verify case 4: legacy .chk recovery
        assert paths["use_checkpoint"] is True
        assert "interrupted_checkpoint.chk" in str(paths["checkpoint"])
        assert "interrupted_system.xml" in str(paths["system"])


# ---------------------------------------------------------------------------
# EQ_INTERRUPTED marker parsing
# ---------------------------------------------------------------------------


class TestEqInterruptedMarker:
    """Tests for EQ_INTERRUPTED marker write/read cycle."""

    def test_marker_format(self, tmp_path):
        """Verify the marker file format written by run_equilibration_stage."""
        stage_dir = tmp_path / "equilibration_2_npt_eq"
        stage_dir.mkdir()
        marker = stage_dir / "EQ_INTERRUPTED"
        marker.write_text(
            "stage_index=2\n"
            "stage_name=npt_eq\n"
            "steps_completed=75000\n"
            "total_steps=100000\n"
            "current_temperature=350.5\n"
            "is_temperature_ramping=True\n"
        )

        # Parse manually (same logic as _find_interrupted_eq_stage)
        info = {}
        for line in marker.read_text().strip().splitlines():
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip()
            if key == "steps_completed":
                info["steps_completed"] = int(value)
            elif key == "total_steps":
                info["total_steps"] = int(value)
            elif key == "current_temperature":
                info["current_temperature"] = float(value)
            elif key == "is_temperature_ramping":
                info["is_temperature_ramping"] = value.lower() == "true"

        assert info["steps_completed"] == 75000
        assert info["total_steps"] == 100000
        assert info["current_temperature"] == 350.5
        assert info["is_temperature_ramping"] is True


class TestFindCompletedEqStages:
    """_find_completed_eq_stages handles completed/interrupted stages."""

    def _make_runner(self, working_dir):
        """Create a minimal mock SimulationRunner for eq stage detection."""
        try:
            from polyzymd.simulation.runner import SimulationRunner

            runner = SimulationRunner.__new__(SimulationRunner)
            runner._working_dir = Path(working_dir)
            return runner
        except ImportError:
            pytest.skip("polyzymd.simulation.runner not importable")

    def _make_stage(self, name):
        """Create a minimal mock EquilibrationStageConfig."""
        stage = MagicMock()
        stage.name = name
        return stage

    def test_all_completed(self, tmp_path):
        stages = [self._make_stage("heating"), self._make_stage("npt_eq")]
        _write_eq_stage(tmp_path, 0, "heating", completed=True)
        _write_eq_stage(tmp_path, 1, "npt_eq", completed=True)
        runner = self._make_runner(tmp_path)
        completed = runner._find_completed_eq_stages(stages)
        assert completed == [0, 1]

    def test_stops_at_interrupted(self, tmp_path):
        """Interrupted stage (with EQ_INTERRUPTED marker) is NOT completed."""
        stages = [self._make_stage("heating"), self._make_stage("npt_eq")]
        _write_eq_stage(tmp_path, 0, "heating", completed=True)
        _write_eq_stage(tmp_path, 1, "npt_eq", interrupted=True)
        runner = self._make_runner(tmp_path)
        completed = runner._find_completed_eq_stages(stages)
        assert completed == [0]  # Only stage 0

    def test_stops_at_first_gap(self, tmp_path):
        stages = [
            self._make_stage("heating"),
            self._make_stage("npt_eq"),
            self._make_stage("final"),
        ]
        _write_eq_stage(tmp_path, 0, "heating", completed=True)
        # Stage 1 missing
        _write_eq_stage(tmp_path, 2, "final", completed=True)
        runner = self._make_runner(tmp_path)
        completed = runner._find_completed_eq_stages(stages)
        assert completed == [0]  # Stops at gap

    def test_empty(self, tmp_path):
        stages = [self._make_stage("heating")]
        runner = self._make_runner(tmp_path)
        completed = runner._find_completed_eq_stages(stages)
        assert completed == []


class TestFindInterruptedEqStage:
    """_find_interrupted_eq_stage detects mid-stage interrupts."""

    def _make_runner(self, working_dir):
        try:
            from polyzymd.simulation.runner import SimulationRunner

            runner = SimulationRunner.__new__(SimulationRunner)
            runner._working_dir = Path(working_dir)
            return runner
        except ImportError:
            pytest.skip("polyzymd.simulation.runner not importable")

    def _make_stage(self, name):
        stage = MagicMock()
        stage.name = name
        return stage

    def test_finds_interrupted_stage(self, tmp_path):
        stages = [self._make_stage("heating"), self._make_stage("npt_eq")]
        _write_eq_stage(tmp_path, 0, "heating", completed=True)
        _write_eq_stage(
            tmp_path,
            1,
            "npt_eq",
            interrupted=True,
            steps_completed=75000,
            total_steps=100000,
            current_temperature=300.0,
        )
        runner = self._make_runner(tmp_path)
        info = runner._find_interrupted_eq_stage(stages, completed_indices=[0])
        assert info is not None
        assert info["stage_index"] == 1
        assert info["steps_completed"] == 75000
        assert info["total_steps"] == 100000
        assert info["current_temperature"] == 300.0

    def test_no_interrupted_stage(self, tmp_path):
        stages = [self._make_stage("heating"), self._make_stage("npt_eq")]
        _write_eq_stage(tmp_path, 0, "heating", completed=True)
        runner = self._make_runner(tmp_path)
        info = runner._find_interrupted_eq_stage(stages, completed_indices=[0])
        assert info is None

    def test_all_completed_returns_none(self, tmp_path):
        stages = [self._make_stage("heating")]
        _write_eq_stage(tmp_path, 0, "heating", completed=True)
        runner = self._make_runner(tmp_path)
        info = runner._find_interrupted_eq_stage(stages, completed_indices=[0])
        assert info is None

    def test_interrupted_without_checkpoint_returns_none(self, tmp_path):
        """EQ_INTERRUPTED marker without checkpoint → can't resume → None."""
        stages = [self._make_stage("heating")]
        dir_name = "equilibration_0_heating"
        stage_dir = tmp_path / dir_name
        stage_dir.mkdir()
        (stage_dir / "EQ_INTERRUPTED").write_text(
            "stage_index=0\nstage_name=heating\n"
            "steps_completed=5000\ntotal_steps=10000\n"
            "current_temperature=300.0\nis_temperature_ramping=False\n"
        )
        # No checkpoint file
        runner = self._make_runner(tmp_path)
        info = runner._find_interrupted_eq_stage(stages, completed_indices=[])
        assert info is None

    def test_temperature_ramping_metadata(self, tmp_path):
        stages = [self._make_stage("ramp")]
        _write_eq_stage(
            tmp_path,
            0,
            "ramp",
            interrupted=True,
            steps_completed=25000,
            total_steps=50000,
            current_temperature=200.0,
            is_temperature_ramping=True,
        )
        runner = self._make_runner(tmp_path)
        info = runner._find_interrupted_eq_stage(stages, completed_indices=[])
        assert info is not None
        assert info["is_temperature_ramping"] is True
        assert info["current_temperature"] == 200.0


# ---------------------------------------------------------------------------
# FAILED segment cleanup in main.py
# ---------------------------------------------------------------------------


class TestFailedSegmentCleanup:
    """FAILED segment directories get removed before retry."""

    def test_failed_dir_removed(self, tmp_path):
        """Simulates the cleanup logic from main.py run_segment()."""
        from polyzymd.simulation.progress import SegmentRecord, SegmentStatus

        # Create a failed segment on disk
        seg_dir = tmp_path / "production_0"
        seg_dir.mkdir()
        (seg_dir / "some_partial_file.txt").write_text("partial")

        # Create a progress with a FAILED record
        from polyzymd.simulation.progress import SimulationProgress

        failed_seg = SegmentRecord(index=0, steps_completed=0, status=SegmentStatus.FAILED)
        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=10000000,
            total_samples_requested=250,
            timestep_fs=2.0,
            segments=[failed_seg],
            replicate=1,
        )

        # Simulate the cleanup logic from main.py
        failed_segments = [s for s in progress.segments if s.status == SegmentStatus.FAILED]
        for failed in failed_segments:
            failed_dir = tmp_path / f"production_{failed.index}"
            if failed_dir.exists():
                shutil.rmtree(failed_dir)
            progress.segments = [s for s in progress.segments if s.index != failed.index]

        # Verify cleanup
        assert not seg_dir.exists()
        assert len(progress.segments) == 0
        assert progress.next_segment_index == 0

    def test_failed_after_completed_cleanup(self, tmp_path):
        """Only the FAILED segment is removed; completed segments survive."""
        from polyzymd.simulation.progress import SegmentRecord, SegmentStatus, SimulationProgress

        # Completed segment 0
        seg0_dir = tmp_path / "production_0"
        seg0_dir.mkdir()
        (seg0_dir / "production_0_state.xml").write_text("<State/>")

        # Failed segment 1
        seg1_dir = tmp_path / "production_1"
        seg1_dir.mkdir()
        (seg1_dir / "partial.txt").write_text("x")

        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=10000000,
            total_samples_requested=250,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=5000000,
                    status=SegmentStatus.COMPLETED,
                    samples_written=125,
                ),
                SegmentRecord(index=1, steps_completed=0, status=SegmentStatus.FAILED),
            ],
            replicate=1,
        )

        # Cleanup
        failed_segments = [s for s in progress.segments if s.status == SegmentStatus.FAILED]
        for failed in failed_segments:
            failed_dir = tmp_path / f"production_{failed.index}"
            if failed_dir.exists():
                shutil.rmtree(failed_dir)
            progress.segments = [s for s in progress.segments if s.index != failed.index]

        # Verify
        assert seg0_dir.exists()
        assert not seg1_dir.exists()
        assert len(progress.segments) == 1
        assert progress.segments[0].index == 0
        assert progress.next_segment_index == 1


# ---------------------------------------------------------------------------
# End-to-end: scan → recovery path selection
# ---------------------------------------------------------------------------


class TestEndToEndRecoveryFlow:
    """Integration test: filesystem scan → correct recovery classification."""

    def test_hard_killed_seg0_becomes_interrupted(self, tmp_path):
        """Segment 0 hard-killed → classified as INTERRUPTED → continuation resumes."""
        from polyzymd.simulation.progress import (
            SegmentStatus,
            get_next_segment_info,
            scan_filesystem,
        )

        _write_hard_killed_segment(tmp_path, 0, with_csv=True, csv_steps=200000)

        progress = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert len(progress.segments) == 1
        assert progress.segments[0].status == SegmentStatus.INTERRUPTED
        assert progress.segments[0].steps_completed == 200000

        # Configure progress for next segment calculation
        progress.total_steps_requested = 10000000
        progress.total_samples_requested = 250

        info = get_next_segment_info(progress, total_steps=10000000, total_samples=250)
        assert info is not None
        assert info["segment_index"] == 1  # Continuation from seg 0
        assert info["steps_to_run"] == 10000000 - 200000

    def test_truly_empty_seg_is_failed(self, tmp_path):
        """Empty directory → FAILED → next run should clean up and retry."""
        from polyzymd.simulation.progress import SegmentStatus, scan_filesystem

        (tmp_path / "production_0").mkdir()
        progress = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert progress.segments[0].status == SegmentStatus.FAILED

    def test_completed_then_hard_killed(self, tmp_path):
        """Segment 0 completed, segment 1 hard-killed → correct recovery."""
        from polyzymd.simulation.progress import (
            SegmentStatus,
            get_next_segment_info,
            scan_filesystem,
        )

        _write_normal_segment(tmp_path, 0)
        _write_hard_killed_segment(tmp_path, 1, with_csv=True, csv_steps=300000)

        progress = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert progress.segments[0].status == SegmentStatus.COMPLETED
        assert progress.segments[1].status == SegmentStatus.INTERRUPTED

        progress.total_steps_requested = 10000000
        progress.total_samples_requested = 250
        # Segment 0: 10ns at 2fs = 5_000_000. Segment 1: 300_000 estimated
        total_done = progress.total_steps_completed
        assert total_done == 5300000

        info = get_next_segment_info(progress, total_steps=10000000, total_samples=250)
        assert info is not None
        assert info["segment_index"] == 2
        assert info["steps_to_run"] == 10000000 - 5300000


# ---------------------------------------------------------------------------
# Concurrency guard: RUNNING segment blocks new segment start
# ---------------------------------------------------------------------------


class TestRunningSegmentConcurrencyGuard:
    """Prevent starting a new segment while one is still RUNNING."""

    def test_running_segment_blocks_new_segment(self, tmp_path):
        """If any segment is RUNNING, the guard prevents advancing."""
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            SimulationProgress,
        )

        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=10000000,
            total_samples_requested=250,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=5000000,
                    status=SegmentStatus.COMPLETED,
                    samples_written=125,
                ),
                SegmentRecord(
                    index=1,
                    steps_completed=1000000,
                    status=SegmentStatus.RUNNING,
                ),
            ],
            replicate=1,
        )

        # Simulate the guard logic from main.py run_segment()
        running_segments = [s for s in progress.segments if s.status == SegmentStatus.RUNNING]
        assert len(running_segments) == 1
        assert running_segments[0].index == 1

    def test_no_running_segments_allows_advance(self, tmp_path):
        """With no RUNNING segments, the guard is a no-op."""
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            SimulationProgress,
            get_next_segment_info,
        )

        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=10000000,
            total_samples_requested=250,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=5000000,
                    status=SegmentStatus.COMPLETED,
                    samples_written=125,
                ),
                SegmentRecord(
                    index=1,
                    steps_completed=1000000,
                    status=SegmentStatus.INTERRUPTED,
                ),
            ],
            replicate=1,
        )

        running_segments = [s for s in progress.segments if s.status == SegmentStatus.RUNNING]
        assert len(running_segments) == 0

        # Advance is allowed — get_next_segment_info returns the next segment
        info = get_next_segment_info(progress, total_steps=10000000, total_samples=250)
        assert info is not None
        assert info["segment_index"] == 2

    def test_scan_filesystem_classifies_recent_checkpoint_as_running(self, tmp_path):
        """A segment with a recent checkpoint (no INTERRUPTED marker) is RUNNING."""
        import os
        import time as _time

        from polyzymd.simulation.progress import SegmentStatus, scan_filesystem

        # Create a segment with a very recent checkpoint
        seg_dir = tmp_path / "production_0"
        seg_dir.mkdir()
        chk = seg_dir / "production_0_checkpoint.chk"
        chk.write_bytes(b"\x00" * 16)
        (seg_dir / "production_0_system.xml").write_text("<System/>")
        # Touch checkpoint to be recent (now)
        now = _time.time()
        os.utime(chk, (now, now))

        progress = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert len(progress.segments) == 1
        assert progress.segments[0].status == SegmentStatus.RUNNING

    def test_scan_filesystem_classifies_stale_checkpoint_as_interrupted(self, tmp_path):
        """A segment with a stale checkpoint (no INTERRUPTED marker) is INTERRUPTED."""
        import os
        import time as _time

        from polyzymd.simulation.progress import SegmentStatus, scan_filesystem

        seg_dir = tmp_path / "production_0"
        seg_dir.mkdir()
        chk = seg_dir / "production_0_checkpoint.chk"
        chk.write_bytes(b"\x00" * 16)
        (seg_dir / "production_0_system.xml").write_text("<System/>")
        # Backdate checkpoint to be stale (20 min ago)
        old = _time.time() - 1200
        os.utime(chk, (old, old))

        progress = scan_filesystem(tmp_path, timestep_fs=2.0)
        assert len(progress.segments) == 1
        assert progress.segments[0].status == SegmentStatus.INTERRUPTED


# ---------------------------------------------------------------------------
# Hard-kill retry-in-place: clean up and retry hard-killed segments
# ---------------------------------------------------------------------------


class TestHardKillRetryInPlace:
    """Hard-killed segments (no INTERRUPTED marker) get cleaned up and retried."""

    def test_hard_killed_segment_cleaned_up(self, tmp_path):
        """Hard-killed last segment is removed so the index is retried."""
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            SimulationProgress,
            get_next_segment_info,
            save_progress,
        )

        # Segment 0 completed normally
        _write_normal_segment(tmp_path, 0)
        # Segment 1 hard-killed (no INTERRUPTED marker)
        _write_hard_killed_segment(tmp_path, 1, with_csv=True, csv_steps=300000)

        # Build progress with segment 1 as INTERRUPTED (from filesystem scan)
        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=10000000,
            total_samples_requested=250,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=5000000,
                    status=SegmentStatus.COMPLETED,
                    samples_written=125,
                ),
                SegmentRecord(
                    index=1,
                    steps_completed=300000,
                    status=SegmentStatus.INTERRUPTED,
                ),
            ],
            replicate=1,
        )
        save_progress(tmp_path, progress)

        # --- Simulate the hard-kill retry logic from main.py ---
        last_seg = max(progress.segments, key=lambda s: s.index)
        assert last_seg.status == SegmentStatus.INTERRUPTED
        last_seg_dir = tmp_path / f"production_{last_seg.index}"
        interrupted_marker = last_seg_dir / "INTERRUPTED"
        assert last_seg_dir.exists()
        assert not interrupted_marker.exists()  # hard-killed: no marker

        # Clean up
        shutil.rmtree(last_seg_dir)
        progress.segments = [s for s in progress.segments if s.index != last_seg.index]
        save_progress(tmp_path, progress)

        # Verify: directory removed, progress updated
        assert not last_seg_dir.exists()
        assert len(progress.segments) == 1
        assert progress.segments[0].index == 0

        # The same segment index (1) should be assigned for retry
        info = get_next_segment_info(progress, total_steps=10000000, total_samples=250)
        assert info is not None
        assert info["segment_index"] == 1  # Retry, not advance to 2!

    def test_graceful_interrupted_segment_not_cleaned_up(self, tmp_path):
        """Gracefully interrupted segment (with INTERRUPTED marker) is preserved."""
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            SimulationProgress,
            get_next_segment_info,
        )

        # Segment 0 completed, segment 1 gracefully interrupted
        _write_normal_segment(tmp_path, 0)
        _write_interrupted_segment(tmp_path, 1, steps_completed=3000000)

        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=10000000,
            total_samples_requested=250,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=5000000,
                    status=SegmentStatus.COMPLETED,
                    samples_written=125,
                ),
                SegmentRecord(
                    index=1,
                    steps_completed=3000000,
                    status=SegmentStatus.INTERRUPTED,
                ),
            ],
            replicate=1,
        )

        # Check that the INTERRUPTED marker exists
        last_seg = max(progress.segments, key=lambda s: s.index)
        last_seg_dir = tmp_path / f"production_{last_seg.index}"
        interrupted_marker = last_seg_dir / "INTERRUPTED"
        assert interrupted_marker.exists()

        # The hard-kill guard should NOT trigger
        should_cleanup = (
            last_seg.status == SegmentStatus.INTERRUPTED
            and last_seg_dir.exists()
            and not interrupted_marker.exists()
        )
        assert not should_cleanup

        # Segment advances normally to 2
        info = get_next_segment_info(progress, total_steps=10000000, total_samples=250)
        assert info is not None
        assert info["segment_index"] == 2

    def test_hard_killed_seg0_retried_from_scratch(self, tmp_path):
        """Hard-killed segment 0 is cleaned up so the run restarts from scratch."""
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            SimulationProgress,
            get_next_segment_info,
            save_progress,
        )

        _write_hard_killed_segment(tmp_path, 0, with_csv=True, csv_steps=100000)

        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=10000000,
            total_samples_requested=250,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=100000,
                    status=SegmentStatus.INTERRUPTED,
                ),
            ],
            replicate=1,
        )
        save_progress(tmp_path, progress)

        # Simulate hard-kill retry logic
        last_seg = max(progress.segments, key=lambda s: s.index)
        last_seg_dir = tmp_path / f"production_{last_seg.index}"
        interrupted_marker = last_seg_dir / "INTERRUPTED"
        assert not interrupted_marker.exists()

        shutil.rmtree(last_seg_dir)
        progress.segments = [s for s in progress.segments if s.index != last_seg.index]
        save_progress(tmp_path, progress)

        assert not last_seg_dir.exists()
        assert len(progress.segments) == 0

        # Segment 0 should be retried
        info = get_next_segment_info(progress, total_steps=10000000, total_samples=250)
        assert info is not None
        assert info["segment_index"] == 0

    def test_completed_segment_not_affected_by_hard_kill_guard(self, tmp_path):
        """A COMPLETED last segment is not affected by the hard-kill logic."""
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            SimulationProgress,
            get_next_segment_info,
        )

        _write_normal_segment(tmp_path, 0)

        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=10000000,
            total_samples_requested=250,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=5000000,
                    status=SegmentStatus.COMPLETED,
                    samples_written=125,
                ),
            ],
            replicate=1,
        )

        last_seg = max(progress.segments, key=lambda s: s.index)
        # Guard only applies to INTERRUPTED status
        assert last_seg.status != SegmentStatus.INTERRUPTED

        # Normal advance to segment 1
        info = get_next_segment_info(progress, total_steps=10000000, total_samples=250)
        assert info is not None
        assert info["segment_index"] == 1
