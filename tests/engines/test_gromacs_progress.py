"""Tests for GROMACS progress tracking helpers."""

from pathlib import Path

from polyzymd.engines.gromacs.progress import (
    _parse_gromacs_log,
    _scan_equilibration_gromacs,
    scan_gromacs_progress,
    update_gromacs_progress,
)
from polyzymd.simulation.progress import load_progress


def test_parse_gromacs_log_extracts_steps_and_completion(tmp_path: Path) -> None:
    """Log parser should extract nsteps, latest step, and finished state."""
    log_path = tmp_path / "prod.log"
    log_path.write_text("""
nsteps                   = 500000

           Step           Time
              0          0.000
         250000        500.000
         500000       1000.000

Finished mdrun
""")

    parsed = _parse_gromacs_log(log_path)
    assert parsed["nsteps_requested"] == 500000
    assert parsed["steps_completed"] == 500000
    assert parsed["time_completed_ps"] == 1000.0
    assert parsed["is_finished"] is True


def test_scan_gromacs_progress_reads_eq_and_prod(tmp_path: Path) -> None:
    """Scanner should detect completed equilibration and production state."""
    (tmp_path / "eq_01.gro").write_text("eq1")
    (tmp_path / "eq_02.gro").write_text("eq2")
    (tmp_path / "prod.log").write_text("""
nsteps = 100000
 Step Time
 50000 100.0
""")

    progress = scan_gromacs_progress(
        working_dir=tmp_path,
        config_path="/path/config.yaml",
        replicate=3,
        total_steps=100000,
        total_samples=200,
        timestep_fs=2.0,
    )

    assert progress.replicate == 3
    assert progress.config_path == "/path/config.yaml"
    assert progress.num_eq_stages_completed == 2
    assert progress.total_steps_completed == 50000


def test_update_gromacs_progress_creates_progress_json(tmp_path: Path) -> None:
    """Updater should create and then update progress.json records."""
    (tmp_path / "prod.log").write_text("""
nsteps = 100000
 Step Time
 25000 50.0
""")

    first = update_gromacs_progress(
        working_dir=tmp_path,
        config_path="/path/config.yaml",
        replicate=1,
        mark_complete=False,
    )
    assert first.total_steps_completed == 25000

    (tmp_path / "prod.log").write_text("""
nsteps = 100000
 Step Time
 100000 200.0
Finished mdrun
""")

    second = update_gromacs_progress(
        working_dir=tmp_path,
        config_path="/path/config.yaml",
        replicate=1,
        mark_complete=True,
    )
    assert second.is_complete
    loaded = load_progress(tmp_path)
    assert loaded is not None
    assert loaded.is_complete


def test_scan_equilibration_gromacs_finds_ordered_eq_files(tmp_path: Path) -> None:
    """Equilibration scan should find eq_NN.gro files in stage order."""
    (tmp_path / "eq_02.gro").write_text("two")
    (tmp_path / "eq_01.gro").write_text("one")
    (tmp_path / "eq_10.gro").write_text("ten")

    records = _scan_equilibration_gromacs(tmp_path)

    assert [r.name for r in records] == ["eq_01", "eq_02", "eq_10"]
    assert all(r.status.value == "completed" for r in records)
