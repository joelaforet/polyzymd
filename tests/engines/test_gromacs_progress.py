"""Tests for GROMACS progress tracking helpers."""

from pathlib import Path

import pytest

from polyzymd.engines.gromacs.progress import (
    _parse_gromacs_log,
    _scan_equilibration_gromacs,
    scan_gromacs_progress,
    update_gromacs_progress,
)
from polyzymd.simulation.progress import SegmentStatus, SimulationStatus, load_progress

#: The lines mdrun writes to its log when it writes a checkpoint.
_CHECKPOINT = "Writing checkpoint, step {step} at Thu Oct  8 00:23:40 2026\n"


def test_parse_gromacs_log_extracts_steps_and_completion(tmp_path: Path) -> None:
    """Log parser should extract nsteps, latest step, and the last checkpoint step."""
    log_path = tmp_path / "prod.log"
    log_path.write_text("""
nsteps                   = 500000

           Step           Time
              0          0.000
         250000        500.000
         500000       1000.000

Writing checkpoint, step 500000 at Thu Oct  8 00:54:26 2026

Finished mdrun
""")

    parsed = _parse_gromacs_log(log_path)
    assert parsed["nsteps_requested"] == 500000
    assert parsed["steps_completed"] == 500000
    assert parsed["time_completed_ps"] == 1000.0
    assert parsed["checkpoint_step"] == 500000


def test_scan_gromacs_progress_reads_eq_and_prod(tmp_path: Path) -> None:
    """Scanner should detect completed equilibration and production state."""
    (tmp_path / "eq_01.gro").write_text("eq1")
    (tmp_path / "eq_02.gro").write_text("eq2")
    (tmp_path / "prod.log").write_text("""
nsteps = 100000
 Step Time
 50000 100.0
Writing checkpoint, step 50000 at Thu Oct  8 00:54:26 2026
""")
    (tmp_path / "state.cpt").write_text("cpt")

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
Writing checkpoint, step 25000 at Thu Oct  8 00:54:26 2026
""")
    (tmp_path / "state.cpt").write_text("cpt")

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
Writing checkpoint, step 100000 at Thu Oct  8 00:54:26 2026
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


def _logged(stamp: str) -> str:
    """Return the UTC ISO time of a time stamp mdrun writes in local time."""
    from datetime import datetime, timezone

    when = datetime.strptime(" ".join(stamp.split()), "%a %b %d %H:%M:%S %Y")
    return when.astimezone(timezone.utc).isoformat()


def test_gromacs_records_take_times_and_seeds_from_mdrun(tmp_path: Path) -> None:
    """Stage and segment records hold when mdrun ran and the seeds it used."""
    (tmp_path / "eq_01.gro").write_text("eq1")
    (tmp_path / "eq_01_heating.mdp").write_text(
        "ld_seed = -1\ngen_vel = yes ; draw velocities\ngen_seed = 2122911565\n"
    )
    (tmp_path / "eq_01.log").write_text(
        "   ld-seed                        = 504874113\n"
        "Started mdrun on rank 0 Wed Oct  7 00:16:48 2026\n"
        "Finished mdrun on rank 0 Wed Oct  7 00:17:17 2026\n"
    )
    (tmp_path / "prod.mdp").write_text("ld_seed = 1482133866\ngen_vel = no\ngen_seed = 5\n")
    (tmp_path / "prod.log").write_text(
        "nsteps = 100\n"
        "   ld-seed                        = 1482133866\n"
        "Started mdrun on rank 0 Wed Oct  7 00:18:20 2026\n"
        " Step Time\n 100 0.2\n"
        "Writing checkpoint, step 100 at Wed Oct  7 00:18:43 2026\n"
        "Finished mdrun on rank 0 Wed Oct  7 00:18:43 2026\n"
    )
    (tmp_path / "state.cpt").write_text("cpt")

    progress = scan_gromacs_progress(tmp_path, total_steps=100)

    (stage,) = progress.equilibration_stages
    assert stage.started_at == _logged("Wed Oct  7 00:16:48 2026")
    assert stage.finished_at == _logged("Wed Oct  7 00:17:17 2026")
    assert stage.seeds == {"ld_seed": 504874113, "gen_seed": 2122911565}
    (segment,) = progress.segments
    assert segment.started_at == _logged("Wed Oct  7 00:18:20 2026")
    assert segment.finished_at == _logged("Wed Oct  7 00:18:43 2026")
    assert segment.seeds == {"ld_seed": 1482133866}


def test_killed_gromacs_run_ends_at_its_last_log_write(tmp_path: Path) -> None:
    """A run with no Finished line after its last start ends at the log's mtime, not now."""
    import os
    from datetime import datetime, timezone

    log = tmp_path / "prod.log"
    log.write_text(
        "nsteps = 100000\n"
        "Started mdrun on rank 0 Wed Oct  7 00:00:00 2026\n"
        "Finished mdrun on rank 0 Wed Oct  7 01:00:00 2026\n"
        "Started mdrun on rank 0 Wed Oct  7 05:00:00 2026\n"
        " Step Time\n 50000 100.0\n"
        "Writing checkpoint, step 50000 at Wed Oct  7 05:30:00 2026\n"
    )
    (tmp_path / "state.cpt").write_text("cpt")
    os.utime(log, (1_790_000_000, 1_790_000_000))

    (segment,) = scan_gromacs_progress(tmp_path, total_steps=100000).segments
    assert segment.started_at == _logged("Wed Oct  7 00:00:00 2026")
    assert segment.finished_at == datetime.fromtimestamp(1_790_000_000, timezone.utc).isoformat()


def test_unreadable_mdp_and_unlogged_random_seeds(tmp_path: Path) -> None:
    """An unreadable MDP is skipped; a seed of -1 that GROMACS chose is recorded as None."""
    from polyzymd.engines.gromacs.progress import _seeds

    unreadable = tmp_path / "eq_01_heating.mdp"
    unreadable.write_text("gen_vel = yes\ngen_seed = 5\n")
    unreadable.chmod(0)
    assert _seeds(unreadable, tmp_path / "eq_01.log") is None
    mdp = tmp_path / "prod.mdp"
    mdp.write_text("ld_seed = -1\ngen_vel = yes\ngen_seed = -1\n")
    assert _seeds(mdp, tmp_path / "prod.log") == {"ld_seed": None, "gen_seed": None}


def test_stage_records_hold_the_ensemble_and_length_that_ran(tmp_path: Path) -> None:
    """A stage with pressure coupling is NPT, and its length is the logged steps times dt."""
    for idx, pcoupl in ((1, "no"), (2, "c-rescale")):
        (tmp_path / f"eq_{idx:02d}.gro").write_text("eq")
        (tmp_path / f"eq_{idx:02d}_stage.mdp").write_text(f"dt = 0.002\npcoupl = {pcoupl}\n")
        (tmp_path / f"eq_{idx:02d}.log").write_text(" Step Time\n 0 0.0\n 5000 10.0\n")

    nvt, npt = _scan_equilibration_gromacs(tmp_path)

    assert (nvt.ensemble, npt.ensemble) == ("NVT", "NPT")
    assert nvt.duration_ns == npt.duration_ns == 0.01


def test_production_records_count_the_frames_written(tmp_path: Path) -> None:
    """samples_written counts the prod.xtc frames, step 0 included, and a restart adds only new ones."""
    (tmp_path / "prod.mdp").write_text("nstxout-compressed = 100\n")
    (tmp_path / "state.cpt").write_text("cpt")
    (tmp_path / "prod.log").write_text("nsteps = 5000\n" + _CHECKPOINT.format(step=2000))
    first = update_gromacs_progress(tmp_path)
    assert first.segments[-1].samples_written == 21

    (tmp_path / "prod.log").write_text("nsteps = 5000\n" + _CHECKPOINT.format(step=5000))
    second = update_gromacs_progress(tmp_path)
    assert second.segments[-1].samples_written == 30
    assert sum(segment.samples_written for segment in second.segments) == 51


def test_records_keep_the_polyzymd_version_that_ran_mdrun(tmp_path: Path, monkeypatch) -> None:
    """The process that ran mdrun records its version; a later scan by another version keeps it."""
    from polyzymd.engines.gromacs.progress import load_or_scan_gromacs_progress

    (tmp_path / "eq_01.gro").write_text("eq")
    (tmp_path / "state.cpt").write_text("cpt")
    (tmp_path / "prod.log").write_text("nsteps = 100\n" + _CHECKPOINT.format(step=100))
    monkeypatch.setattr("polyzymd.utils.version.get_polyzymd_version", lambda: "1.3.0")
    progress = update_gromacs_progress(tmp_path)
    assert [r.polyzymd_version for r in progress.equilibration_stages] == ["1.3.0"]
    assert [r.polyzymd_version for r in progress.segments] == ["1.3.0"]
    assert progress.segments[0].openmm_version is None

    monkeypatch.setattr("polyzymd.utils.version.get_polyzymd_version", lambda: "9.9.9")
    rescanned = load_or_scan_gromacs_progress(tmp_path)
    assert rescanned.equilibration_stages[0].polyzymd_version == "1.3.0"
    assert scan_gromacs_progress(tmp_path).segments[0].polyzymd_version is None


def test_update_records_its_version_only_on_what_it_ran(tmp_path: Path, monkeypatch) -> None:
    """Records saved earlier without a version stay without one; new and re-run ones get it."""
    from polyzymd.engines.gromacs.progress import load_or_scan_gromacs_progress
    from polyzymd.simulation.progress import save_progress

    (tmp_path / "eq_01.gro").write_text("eq")
    (tmp_path / "eq_01.log").write_text("Started mdrun on rank 0 Wed Oct  7 00:00:00 2026\n")
    (tmp_path / "state.cpt").write_text("cpt")
    (tmp_path / "prod.log").write_text("nsteps = 200\n" + _CHECKPOINT.format(step=100))
    save_progress(tmp_path, load_or_scan_gromacs_progress(tmp_path))
    monkeypatch.setattr("polyzymd.utils.version.get_polyzymd_version", lambda: "1.3.0")

    (tmp_path / "prod.log").write_text("nsteps = 200\n" + _CHECKPOINT.format(step=200))
    progress = update_gromacs_progress(tmp_path)
    assert [r.polyzymd_version for r in progress.segments] == [None, "1.3.0"]
    assert progress.equilibration_stages[0].polyzymd_version is None

    progress.equilibration_stages[0].polyzymd_version = "1.2.0"
    save_progress(tmp_path, progress)
    (tmp_path / "eq_01.log").write_text(
        "Started mdrun on rank 0 Wed Oct  7 02:00:00 2026\n"
        "Finished mdrun on rank 0 Wed Oct  7 02:01:00 2026\n"
    )
    rerun = update_gromacs_progress(tmp_path)
    assert rerun.equilibration_stages[0].polyzymd_version == "1.3.0"


def test_unreadable_stage_log_does_not_stop_a_scan(tmp_path: Path) -> None:
    """A stage log that cannot be read gives a stage record of length 0."""
    (tmp_path / "eq_01.gro").write_text("eq")
    (tmp_path / "eq_01.log").mkdir()

    (stage,) = _scan_equilibration_gromacs(tmp_path)
    assert stage.duration_ns == 0.0


def test_update_records_its_version_on_what_status_saved_during_the_job(
    tmp_path: Path, monkeypatch
) -> None:
    """Records that status saved during the job get its version; records of earlier jobs do not."""
    from polyzymd.engines.gromacs.progress import load_or_scan_gromacs_progress
    from polyzymd.simulation.progress import save_progress

    for idx, hour in ((1, "01"), (2, "03")):
        (tmp_path / f"eq_{idx:02d}.gro").write_text("eq")
        (tmp_path / f"eq_{idx:02d}.log").write_text(
            f"Started mdrun on rank 0 Wed Oct  7 {hour}:00:00 2026\n"
            f"Finished mdrun on rank 0 Wed Oct  7 {hour}:30:00 2026\n"
        )
    (tmp_path / "prod.log").write_text(
        "nsteps = 200\nStarted mdrun on rank 0 Wed Oct  7 04:00:00 2026\n"
        + _CHECKPOINT.format(step=200)
    )
    (tmp_path / "state.cpt").write_text("cpt")
    save_progress(tmp_path, load_or_scan_gromacs_progress(tmp_path))
    monkeypatch.setattr("polyzymd.utils.version.get_polyzymd_version", lambda: "1.3.0")

    progress = update_gromacs_progress(tmp_path, since=_logged("Wed Oct  7 02:00:00 2026"))

    assert [r.polyzymd_version for r in progress.equilibration_stages] == [None, "1.3.0"]
    assert [r.polyzymd_version for r in progress.segments] == ["1.3.0"]


#: prod.log of a production run that ``-maxh`` stopped (GROMACS 2024.2 on Blanca).
_MAXH_STOPPED_LOG = """\
   nsteps                         = 10000000
Started mdrun on rank 0 Thu Oct  8 00:20:42 2026
           Step           Time
        1700000     3400.00000

Step 1703100: Run time exceeded 0.050 hours, will terminate the run within 200 steps

Writing checkpoint, step 1703150 at Thu Oct  8 00:23:40 2026

Finished mdrun on rank 0 Thu Oct  8 00:23:40 2026
"""

#: prod.log of a production run that a TERM signal stopped.
_TERM_STOPPED_LOG = """\
   nsteps                         = 10000000
Started mdrun on rank 0 Thu Oct  8 00:20:42 2026
           Step           Time
        1700000     3400.00000

Received the TERM signal, stopping within 200 steps

Writing checkpoint, step 1703150 at Thu Oct  8 00:23:40 2026

Finished mdrun on rank 0 Thu Oct  8 00:23:40 2026
"""

#: prod.log of a production run that reached nsteps.
_FINISHED_LOG = """\
   nsteps                         = 10000000
Started mdrun on rank 0 Thu Oct  8 00:32:30 2026
           Step           Time
       10000000    20000.00000

Writing checkpoint, step 10000000 at Thu Oct  8 00:54:26 2026

Finished mdrun on rank 0 Thu Oct  8 00:54:26 2026
"""


@pytest.mark.parametrize("log", [_MAXH_STOPPED_LOG, _TERM_STOPPED_LOG])
def test_a_run_that_mdrun_stopped_early_is_interrupted(tmp_path: Path, log: str) -> None:
    """mdrun writes "Finished mdrun" also when -maxh or TERM stops it; the run is not complete."""
    (tmp_path / "state.cpt").write_text("cpt")
    (tmp_path / "prod.log").write_text(log)

    progress = update_gromacs_progress(tmp_path)

    assert progress.status == SimulationStatus.INTERRUPTED
    assert progress.total_steps_completed == 1703150
    assert [segment.status for segment in progress.segments] == [SegmentStatus.INTERRUPTED]


def test_a_run_at_nsteps_is_complete(tmp_path: Path) -> None:
    """A run whose last checkpoint is at nsteps is complete."""
    (tmp_path / "state.cpt").write_text("cpt")
    (tmp_path / "prod.log").write_text(_FINISHED_LOG)

    progress = update_gromacs_progress(tmp_path)

    assert progress.status == SimulationStatus.COMPLETED
    assert progress.total_steps_completed == 10000000


def test_steps_after_the_last_checkpoint_are_not_counted(tmp_path: Path) -> None:
    """A killed run counts up to its last checkpoint, which is where the restart begins."""
    (tmp_path / "state.cpt").write_text("cpt")
    (tmp_path / "prod.log").write_text(
        "nsteps = 2500000\n" + _CHECKPOINT.format(step=300000) + " Step Time\n 375000 750.0\n"
    )

    assert scan_gromacs_progress(tmp_path).total_steps_completed == 300000
    (tmp_path / "state.cpt").unlink()
    assert scan_gromacs_progress(tmp_path).total_steps_completed == 0


def test_a_production_restart_from_step_0_replaces_the_lost_segments(tmp_path: Path) -> None:
    """Without state.cpt production starts again from step 0; the lost steps leave the records."""
    from polyzymd.engines.gromacs.progress import load_or_scan_gromacs_progress
    from polyzymd.simulation.progress import save_progress

    first = "nsteps = 20000\nStarted mdrun on rank 0 Thu Oct  8 06:46:46 2026\n"
    (tmp_path / "state.cpt").write_text("cpt")
    (tmp_path / "prod.log").write_text(first + _CHECKPOINT.format(step=5000))
    save_progress(tmp_path, update_gromacs_progress(tmp_path))

    # Those 5000 steps were then lost: the next job found no state.cpt and
    # started again from step 0, and mdrun backed up the old prod.log.
    (tmp_path / "state.cpt").unlink()
    (tmp_path / "prod.log").write_text(first + " Step Time\n 3000 6.0\n")
    assert load_or_scan_gromacs_progress(tmp_path).total_steps_completed == 0

    (tmp_path / "state.cpt").write_text("cpt")
    (tmp_path / "prod.log").write_text(
        "nsteps = 20000\nStarted mdrun on rank 0 Thu Oct  8 06:48:29 2026\n"
        + _CHECKPOINT.format(step=20000)
    )
    progress = update_gromacs_progress(tmp_path)

    assert [(s.steps_completed, s.status) for s in progress.segments] == [
        (20000, SegmentStatus.COMPLETED)
    ]
    assert progress.status == SimulationStatus.COMPLETED
