"""Tests for atomic lifecycle phase records."""

from polyzymd.simulation.phase_state import load_phase_record, phase_completed, write_phase_record


def test_only_completed_record_skips_phase(tmp_path):
    path = tmp_path / "phase.json"
    write_phase_record(path, phase="equilibration_0", status="started", total_steps=100)
    assert not phase_completed(path)
    write_phase_record(
        path,
        phase="equilibration_0",
        status="recovery",
        step=40,
        total_steps=100,
        temperature=280.0,
    )
    assert not phase_completed(path)
    write_phase_record(path, phase="equilibration_0", status="completed", step=100)
    assert phase_completed(path)


def test_invalid_record_is_not_completion(tmp_path):
    path = tmp_path / "phase.json"
    path.write_text('{"status":"completed"')
    assert load_phase_record(path) is None
    assert not phase_completed(path)


def test_minimization_record_carries_the_frozen_and_hydrogen_metrics(tmp_path):
    """The minimization record reports what was frozen and how far hydrogens moved."""
    import json

    path = tmp_path / "phase.json"
    write_phase_record(
        path,
        phase="minimization",
        status="completed",
        frozen_atoms=2497,
        frozen_rmsd_angstrom=0.0,
        hydrogen_max_displacement_angstrom=0.234,
    )
    record = load_phase_record(path)
    assert record is not None
    assert record.frozen_atoms == 2497
    assert record.frozen_rmsd_angstrom == 0.0
    assert record.hydrogen_max_displacement_angstrom == 0.234
    assert json.loads(path.read_text())["hydrogen_max_displacement_angstrom"] == 0.234


def test_hydrogen_displacement_defaults_to_none(tmp_path):
    """Records written before the field existed, and unfrozen runs, stay valid."""
    path = tmp_path / "phase.json"
    write_phase_record(path, phase="minimization", status="completed")
    record = load_phase_record(path)
    assert record is not None
    assert record.hydrogen_max_displacement_angstrom is None
