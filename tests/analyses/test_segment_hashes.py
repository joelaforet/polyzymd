"""Trajectories identified by content: segment hashes in progress.json, the hash cache, and records.

Conditions are synthetic OpenMM runs (four unit-mass atoms on a cross).
"""

from __future__ import annotations

import hashlib
import io
import os
import shutil
from pathlib import Path

import pytest

import polyzymd as pz
from polyzymd.analyses.shared.file_hashes import file_sha256, recorded_segment_hashes
from polyzymd.simulation.progress import (
    SegmentRecord,
    SegmentStatus,
    SimulationProgress,
    flush_reporters,
    load_progress,
    save_progress,
    trajectory_digest,
)
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]


def _progress(path: Path, **segment: object) -> None:
    save_progress(
        path,
        SimulationProgress(
            config_path="config.yaml",
            total_steps_requested=1000,
            total_samples_requested=10,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=1000,
                    steps_requested=1000,
                    samples_written=10,
                    status=SegmentStatus.COMPLETED,
                    **segment,
                )
            ],
        ),
    )


class TestRunnerBookkeeping:
    def test_trajectory_digest(self, tmp_path: Path) -> None:
        trajectory = tmp_path / "t.dcd"
        trajectory.write_bytes(b"frames")
        digest = trajectory_digest(trajectory)
        assert digest == {
            "trajectory_sha256": hashlib.sha256(b"frames").hexdigest(),
            "trajectory_bytes": 6,
        }
        assert trajectory_digest(tmp_path / "missing.dcd") == {
            "trajectory_sha256": None,
            "trajectory_bytes": None,
        }

    def test_flush_reporters_writes_buffered_frames(self, tmp_path: Path) -> None:
        path = tmp_path / "t.dcd"
        handle = path.open("wb", buffering=1 << 20)
        handle.write(b"x" * 100)
        reporter = type("Reporter", (), {"_out": handle})()
        assert path.stat().st_size == 0
        flush_reporters(type("Simulation", (), {"reporters": [reporter]})())
        assert path.stat().st_size == 100
        handle.close()

    def test_continuation_records_hash_and_versions(self, tmp_path: Path) -> None:
        from polyzymd.simulation.continuation import ContinuationManager

        working = tmp_path / "run"
        (working / "production_1").mkdir(parents=True)
        (working / "production_1" / "production_1_trajectory.dcd").write_bytes(b"segment one")
        _progress(working)
        manager = object.__new__(ContinuationManager)
        manager._working_dir = working
        manager._segment_index = 1
        manager._update_progress_completed(
            total_steps=1000, num_samples=10, duration_ns=1.0, timestep_fs=2.0
        )
        segment = next(s for s in load_progress(working).segments if s.index == 1)
        assert segment.trajectory_sha256 == hashlib.sha256(b"segment one").hexdigest()
        assert segment.trajectory_bytes == len(b"segment one")
        assert segment.polyzymd_version == pz.__version__

    def test_first_segment_records_hash(self, tmp_path: Path) -> None:
        from polyzymd.simulation.runner import SimulationRunner

        working = tmp_path / "run"
        (working / "production_0").mkdir(parents=True)
        (working / "production_0" / "production_0_trajectory.dcd").write_bytes(b"segment zero")
        save_progress(
            working,
            SimulationProgress(
                config_path="config.yaml",
                total_steps_requested=1000,
                total_samples_requested=10,
                timestep_fs=2.0,
            ),
        )
        runner = object.__new__(SimulationRunner)
        runner._working_dir = working
        runner._update_progress_completed(
            segment_index=0, total_steps=1000, num_samples=10, duration_ns=1.0, timestep_fs=2.0
        )
        segment = load_progress(working).segments[0]
        assert segment.trajectory_sha256 == hashlib.sha256(b"segment zero").hexdigest()


class TestHashCache:
    def test_cache_is_reused_until_the_file_changes(self, tmp_path: Path, monkeypatch) -> None:
        file = tmp_path / "t.dcd"
        file.write_bytes(b"one")
        first = file_sha256(file)
        assert first == hashlib.sha256(b"one").hexdigest()
        reads = []
        real_open = Path.open

        def counting_open(self, *args, **kwargs):
            if self == file.resolve():
                reads.append(1)
            return real_open(self, *args, **kwargs)

        monkeypatch.setattr(Path, "open", counting_open)
        assert file_sha256(file) == first and reads == []
        file.write_bytes(b"two")
        os.utime(file, ns=(1, 1))
        assert file_sha256(file) == hashlib.sha256(b"two").hexdigest()

    def test_recorded_hash_is_trusted_when_the_size_matches(self, tmp_path: Path) -> None:
        file = tmp_path / "t.dcd"
        file.write_bytes(b"abc")
        assert file_sha256(file, ("recorded", 3)) == "recorded"
        assert file_sha256(file, ("recorded", 99)) == hashlib.sha256(b"abc").hexdigest()

    def test_recorded_segment_hashes(self, tmp_path: Path) -> None:
        _progress(tmp_path, trajectory_sha256="f" * 64, trajectory_bytes=12)
        assert recorded_segment_hashes(tmp_path) == {"production_0_trajectory.dcd": ("f" * 64, 12)}
        assert recorded_segment_hashes(tmp_path / "nothing") == {}


CALLS: list[int] = []


def counted_rg(atoms):
    CALLS.append(1)
    return atoms.radius_of_gyration()


class TestRecords:
    def _config(self, tmp_path: Path) -> Path:
        config = write_simulation_config(tmp_path / "A", scratch=tmp_path / "data")
        write_openmm_replicate(config, 1, [1.0 + 0.01 * k for k in range(10)])
        return config

    def test_record_names_files_by_content(self, tmp_path: Path) -> None:
        study = pz.Study.from_configs({"A": self._config(tmp_path)}, equilibration="0ns")
        identity = study["A"].replicates[0].identity
        for item in (identity["topology"], *identity["trajectories"]):
            assert len(item["sha256"]) == 64 and "mtime_ns" not in item
            assert not Path(item["path"]).is_absolute()

    def test_copied_data_reuses_stored_results(self, tmp_path: Path) -> None:
        config = self._config(tmp_path)
        out = tmp_path / "results"
        CALLS.clear()
        pz.Study.from_configs({"A": config}, equilibration="0ns").timeseries(
            counted_rg, pz.select("all"), unit="A", output_dir=out
        )
        measured = len(CALLS)
        assert measured > 0
        shutil.copytree(tmp_path / "data", tmp_path / "copy")  # new modification times
        for path in (tmp_path / "copy").rglob("*"):
            if path.is_file():
                os.utime(path, None)
        CALLS.clear()
        pz.Study.from_configs(
            {"A": config}, equilibration="0ns", data={"A": tmp_path / "copy"}
        ).timeseries(counted_rg, pz.select("all"), unit="A", output_dir=out)
        assert CALLS == []

    def test_changed_trajectory_recomputes(self, tmp_path: Path) -> None:
        config = self._config(tmp_path)
        out = tmp_path / "results"
        pz.Study.from_configs({"A": config}, equilibration="0ns").timeseries(
            counted_rg, pz.select("all"), unit="A", output_dir=out
        )
        write_openmm_replicate(config, 1, [5.0 + 0.01 * k for k in range(10)])
        CALLS.clear()
        series = pz.Study.from_configs({"A": config}, equilibration="0ns").timeseries(
            counted_rg, pz.select("all"), unit="A", output_dir=out
        )
        assert CALLS and series.reduce().summary().conditions[0].mean > 4.0


def test_buffered_reporter_type_is_tolerated() -> None:
    flush_reporters(
        type("Simulation", (), {"reporters": [object(), type("R", (), {"_out": io.BytesIO()})()]})()
    )


class TestHashExistingRuns:
    """record_trajectory_hashes: hashing a finished run is idempotent and never overwrites."""

    def _run(
        self, tmp_path: Path, status=SegmentStatus.COMPLETED, overall=None, **recorded
    ) -> Path:
        from polyzymd.simulation.progress import SimulationStatus

        working = tmp_path / "run"
        for index in (0, 1):
            (working / f"production_{index}").mkdir(parents=True, exist_ok=True)
            (working / f"production_{index}" / f"production_{index}_trajectory.dcd").write_bytes(
                f"segment {index}".encode()
            )
        save_progress(
            working,
            SimulationProgress(
                config_path="config.yaml",
                total_steps_requested=2000,
                total_samples_requested=20,
                timestep_fs=2.0,
                status=overall or SimulationStatus.COMPLETED,
                segments=[
                    SegmentRecord(index=0, steps_completed=1000, steps_requested=1000,
                                  samples_written=10, status=status, **recorded),
                    SegmentRecord(index=1, steps_completed=1000, steps_requested=1000,
                                  samples_written=10, status=SegmentStatus.INTERRUPTED),
                    SegmentRecord(index=2, steps_completed=1000, steps_requested=1000,
                                  samples_written=10, status=SegmentStatus.COMPLETED),
                ],
            ),
        )  # fmt: skip
        return working

    def test_records_then_changes_nothing(self, tmp_path: Path) -> None:
        from polyzymd.simulation.progress import record_trajectory_hashes

        working = self._run(tmp_path)
        first = record_trajectory_hashes(working)
        assert first["hashed"] == [0, 1] and first["missing"] == [2]
        segments = load_progress(working).segments
        assert segments[0].trajectory_sha256 == hashlib.sha256(b"segment 0").hexdigest()
        assert segments[1].trajectory_bytes == len(b"segment 1")
        before = (
            (working / "progress.json").read_bytes(),
            (working / "progress.json").stat().st_mtime_ns,
        )
        second = record_trajectory_hashes(working)
        assert second["hashed"] == [] and second["recorded"] == [0, 1]
        after = (
            (working / "progress.json").read_bytes(),
            (working / "progress.json").stat().st_mtime_ns,
        )
        assert after == before
        assert record_trajectory_hashes(working, verify=True)["verified"] == [0, 1]

    def test_never_overwrites_a_recorded_hash(self, tmp_path: Path) -> None:
        from polyzymd.simulation.progress import record_trajectory_hashes

        working = self._run(tmp_path, trajectory_sha256="0" * 64, trajectory_bytes=99)
        report = record_trajectory_hashes(working)
        assert any("records 99 bytes" in c for c in report["conflicts"])
        assert load_progress(working).segments[0].trajectory_sha256 == "0" * 64
        working = self._run(
            tmp_path / "b", trajectory_sha256="0" * 64, trajectory_bytes=len(b"segment 0")
        )
        report = record_trajectory_hashes(working, verify=True)
        assert any("differs from the one recorded" in c for c in report["conflicts"])
        assert load_progress(working).segments[0].trajectory_sha256 == "0" * 64

    def test_running_runs_are_left_alone(self, tmp_path: Path) -> None:
        from polyzymd.simulation.progress import SimulationStatus, record_trajectory_hashes

        working = self._run(tmp_path, overall=SimulationStatus.RUNNING)
        assert "running" in record_trajectory_hashes(working)["skipped"]
        assert load_progress(working).segments[0].trajectory_sha256 is None
        assert record_trajectory_hashes(working, force=True)["hashed"] == [0, 1]

    def test_dry_run_writes_nothing(self, tmp_path: Path) -> None:
        from polyzymd.simulation.progress import record_trajectory_hashes

        working = self._run(tmp_path)
        before = (working / "progress.json").read_bytes()
        assert record_trajectory_hashes(working, dry_run=True)["hashed"] == [0, 1]
        assert (working / "progress.json").read_bytes() == before

    def test_cli_reports_and_exits_2_on_conflict(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        config = write_simulation_config(tmp_path / "A", scratch=tmp_path / "data")
        write_openmm_replicate(config, 1, [1.0 + 0.01 * k for k in range(5)])
        working = pz.Study.from_configs({"A": config}, equilibration="0ns")[
            "A"
        ].config.get_working_directory(1)
        _progress(working)
        result = CliRunner().invoke(cli, ["hash-trajectories", "-c", str(config)])
        assert result.exit_code == 0, result.output
        assert "A replicate 1: hashed 1, already recorded 0" in result.output
        again = CliRunner().invoke(cli, ["hash-trajectories", "-c", str(config)])
        assert "hashed 0, already recorded 1" in again.output
        progress = load_progress(working)
        progress.segments[0].trajectory_bytes = 1
        save_progress(working, progress)
        conflict = CliRunner().invoke(cli, ["hash-trajectories", "-c", str(config)])
        assert conflict.exit_code == 2 and "conflict:" in conflict.output
