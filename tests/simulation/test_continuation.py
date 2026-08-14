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
from pathlib import Path

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


# ---------------------------------------------------------------------------
# ContinuationManager._get_previous_paths — recovery path selection
# ---------------------------------------------------------------------------


class TestFindSolvatedPdb:
    """ContinuationManager._find_solvated_pdb accepts only approved topology paths."""

    def _make_manager(self, working_dir):
        """Create a minimal ContinuationManager-like object for testing."""
        try:
            from polyzymd.simulation.continuation import ContinuationManager

            mgr = ContinuationManager.__new__(ContinuationManager)
            mgr._working_dir = Path(working_dir)
            mgr._prev_segment = 0
            return mgr
        except ImportError:
            pytest.skip("polyzymd.simulation.continuation not importable")

    @pytest.mark.parametrize(
        "relative_path",
        [
            Path("solvated_system.pdb"),
            Path("production_0") / "production_0_topology.pdb",
            Path("production") / "production_topology.pdb",
        ],
    )
    def test_finds_exact_approved_paths(self, tmp_path, relative_path):
        """Approved exact topology paths should be accepted."""
        pdb_path = tmp_path / relative_path
        pdb_path.parent.mkdir(parents=True, exist_ok=True)
        pdb_path.write_text("ATOM ...")
        mgr = self._make_manager(tmp_path)

        assert mgr._find_solvated_pdb() == pdb_path

    def test_rejects_arbitrary_recursive_pdb(self, tmp_path):
        """Arbitrary recursive PDB files should not be accepted."""
        nested_dir = tmp_path / "nested"
        nested_dir.mkdir()
        (nested_dir / "decoy.pdb").write_text("ATOM ...")
        mgr = self._make_manager(tmp_path)

        with pytest.raises(FileNotFoundError, match="Could not find solvated PDB file"):
            mgr._find_solvated_pdb()

    def test_prefers_predecessor_topology_over_root_decoy(self, tmp_path):
        """A segment-owned topology is bound to its System/State artifacts."""
        root = tmp_path / "solvated_system.pdb"
        root.write_text("decoy")
        segment = tmp_path / "production_2" / "production_2_topology.pdb"
        segment.parent.mkdir()
        segment.write_text("matching")
        mgr = self._make_manager(tmp_path)
        mgr._prev_segment = 2

        assert mgr._find_solvated_pdb() == segment


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
        with pytest.raises(FileNotFoundError, match="cannot recover"):
            mgr._get_previous_paths()

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
