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
