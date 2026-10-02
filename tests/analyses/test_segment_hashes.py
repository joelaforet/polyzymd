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
from polyzymd.analyses.shared.file_hashes import file_sha256
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

    def test_openmm_engine_reads_segment_hashes(self, tmp_path: Path) -> None:
        engine = _engine(tmp_path, "openmm")
        _progress(tmp_path, trajectory_sha256="f" * 64, trajectory_bytes=12)
        dcd = (tmp_path / "production_0" / "production_0_trajectory.dcd").resolve()
        assert engine.recorded_trajectory_hashes(tmp_path) == {dcd: ("f" * 64, 12)}
        assert engine.recorded_trajectory_hashes(tmp_path / "nothing") == {}


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


def _engine(tmp_path: Path, name: str):
    """Return the ``name`` engine of a minimal config, without its binary."""
    import yaml

    from polyzymd.config.schema import SimulationConfig
    from polyzymd.engines import create_engine

    config = write_simulation_config(tmp_path / f"config_{name}", scratch=tmp_path / "data")
    data = yaml.safe_load(config.read_text())
    data["engine"] = name
    config.write_text(yaml.safe_dump(data))
    return create_engine(SimulationConfig.from_yaml(config), defer_binary=True)


DCD = "production_{0}/production_{0}_trajectory.dcd"


class TestHashExistingRuns:
    """SimulationEngine.record_trajectory_hashes: idempotent, never overwrites, for every engine."""

    def _run(
        self, tmp_path: Path, status=SegmentStatus.COMPLETED, overall=None, **recorded
    ) -> Path:
        from polyzymd.simulation.progress import SimulationStatus

        working = tmp_path / "run"
        for index in (0, 1):
            (working / f"production_{index}").mkdir(parents=True, exist_ok=True)
            (working / DCD.format(index)).write_bytes(f"segment {index}".encode())
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
        engine = _engine(tmp_path, "openmm")
        working = self._run(tmp_path)
        first = engine.record_trajectory_hashes(working, 1)
        assert first["hashed"] == [DCD.format(0), DCD.format(1)]
        progress = load_progress(working)
        assert progress.segments[0].trajectory_sha256 == hashlib.sha256(b"segment 0").hexdigest()
        assert progress.segments[1].trajectory_bytes == len(b"segment 1")
        assert progress.trajectory_hashes[DCD.format(1)].bytes == len(b"segment 1")
        # Only hashes are added: segment 2 has no file and no state.xml, and stays completed.
        assert [x.status for x in progress.segments] == [
            SegmentStatus.COMPLETED,
            SegmentStatus.INTERRUPTED,
            SegmentStatus.COMPLETED,
        ]
        before = (
            (working / "progress.json").read_bytes(),
            (working / "progress.json").stat().st_mtime_ns,
        )
        second = engine.record_trajectory_hashes(working, 1)
        assert second["hashed"] == [] and second["recorded"] == [DCD.format(0), DCD.format(1)]
        after = (
            (working / "progress.json").read_bytes(),
            (working / "progress.json").stat().st_mtime_ns,
        )
        assert after == before
        verified = engine.record_trajectory_hashes(working, 1, verify=True)["verified"]
        assert verified == [DCD.format(0), DCD.format(1)]

    def test_never_overwrites_a_recorded_hash(self, tmp_path: Path) -> None:
        engine = _engine(tmp_path, "openmm")
        working = self._run(tmp_path, trajectory_sha256="0" * 64, trajectory_bytes=99)
        report = engine.record_trajectory_hashes(working, 1)
        assert any("records 99 bytes" in c for c in report["conflicts"])
        assert load_progress(working).segments[0].trajectory_sha256 == "0" * 64
        working = self._run(
            tmp_path / "b", trajectory_sha256="0" * 64, trajectory_bytes=len(b"segment 0")
        )
        report = engine.record_trajectory_hashes(working, 1, verify=True)
        assert any("differs from the one recorded" in c for c in report["conflicts"])
        assert load_progress(working).segments[0].trajectory_sha256 == "0" * 64

    def test_running_runs_are_left_alone(self, tmp_path: Path) -> None:
        from polyzymd.simulation.progress import SimulationStatus

        engine = _engine(tmp_path, "openmm")
        working = self._run(tmp_path, overall=SimulationStatus.RUNNING)
        assert "running" in engine.record_trajectory_hashes(working, 1)["skipped"]
        assert load_progress(working).segments[0].trajectory_sha256 is None
        forced = engine.record_trajectory_hashes(working, 1, force=True)
        assert forced["hashed"] == [DCD.format(0), DCD.format(1)]

    def test_dry_run_writes_nothing(self, tmp_path: Path) -> None:
        engine = _engine(tmp_path, "openmm")
        working = self._run(tmp_path)
        before = (working / "progress.json").read_bytes()
        assert engine.record_trajectory_hashes(working, 1, dry_run=True)["hashed"] == [
            DCD.format(0),
            DCD.format(1),
        ]
        assert (working / "progress.json").read_bytes() == before

    def test_gromacs_hashes_its_trajectories(self, tmp_path: Path) -> None:
        engine = _engine(tmp_path, "gromacs")
        working = tmp_path / "run" / "gromacs"
        working.mkdir(parents=True)
        for name in ("prod.xtc", "prod_nojump.xtc"):
            (working / name).write_bytes(name.encode())
        (working / "prod.edr").write_bytes(b"not a trajectory")
        first = engine.record_trajectory_hashes(working, 1)
        assert first["hashed"] == ["prod.xtc", "prod_nojump.xtc"] and first["created"]
        recorded = load_progress(working).trajectory_hashes
        assert recorded["prod.xtc"].sha256 == hashlib.sha256(b"prod.xtc").hexdigest()
        assert engine.recorded_trajectory_hashes(working)[(working / "prod.xtc").resolve()] == (
            hashlib.sha256(b"prod.xtc").hexdigest(),
            len(b"prod.xtc"),
        )
        before = (working / "progress.json").read_bytes()
        second = engine.record_trajectory_hashes(working, 1)
        assert second["hashed"] == [] and not second["created"]
        assert (working / "progress.json").read_bytes() == before
        (working / "prod.xtc").write_bytes(b"appended by a later mdrun")
        conflict = engine.record_trajectory_hashes(working, 1)["conflicts"]
        assert any(c.startswith("prod.xtc:") for c in conflict)

    def test_missing_working_directory_is_skipped(self, tmp_path: Path) -> None:
        for name in ("openmm", "gromacs"):
            report = _engine(tmp_path, name).record_trajectory_hashes(tmp_path / "absent", 1)
            assert report["skipped"] and report["hashed"] == []

    def test_cli_reports_and_exits_2_on_conflict(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        config = write_simulation_config(tmp_path / "A", scratch=tmp_path / "data")
        write_openmm_replicate(config, 1, [1.0 + 0.01 * k for k in range(5)])
        condition = pz.Study.from_configs({"A": config}, equilibration="0ns")["A"]
        working = condition.config.get_working_directory(1)
        _progress(working)
        result = CliRunner().invoke(cli, ["hash-trajectories", "-c", str(config)])
        assert result.exit_code == 0, result.output
        assert "A replicate 1 (openmm): hashed 1, already recorded 0" in result.output
        again = CliRunner().invoke(cli, ["hash-trajectories", "-c", str(config)])
        assert "hashed 0, already recorded 1" in again.output
        dcd = (working / DCD.format(0)).resolve()
        assert dcd in condition._provider.recorded_trajectory_hashes(1)
        progress = load_progress(working)
        progress.segments[0].trajectory_bytes = 1
        progress.trajectory_hashes.clear()
        save_progress(working, progress)
        conflict = CliRunner().invoke(cli, ["hash-trajectories", "-c", str(config)])
        assert conflict.exit_code == 2 and "conflict:" in conflict.output
        none = CliRunner().invoke(
            cli, ["hash-trajectories", "-c", str(config), "--replicates", "7"]
        )
        assert none.exit_code == 0 and "A: no runs found in" in none.output


class TestFreezeNudge:
    def test_freeze_names_the_command_until_hashes_are_recorded(
        self, tmp_path: Path, monkeypatch
    ) -> None:
        import subprocess

        from polyzymd.analyses.study_freeze import freeze
        from polyzymd.analyses.study_scaffold import create_study

        for key in ("GIT_AUTHOR_NAME", "GIT_COMMITTER_NAME"):
            monkeypatch.setenv(key, "Test")
        for key in ("GIT_AUTHOR_EMAIL", "GIT_COMMITTER_EMAIL"):
            monkeypatch.setenv(key, "test@example.com")
        config = write_simulation_config(tmp_path / "A", scratch=tmp_path / "data")
        (config.parent / "test.pdb").write_text("REMARK\nEND\n")
        write_openmm_replicate(config, 1, [1.0 + 0.01 * k for k in range(5)])
        root = tmp_path / "study"
        create_study(root, conditions={"A": config}, equilibration="0ns")
        working = pz.Study(root)["A"].config.get_working_directory(1)
        assert (
            any("hash-trajectories" in w for w in freeze(root).warnings) is False
        )  # no progress.json
        _progress(working)
        subprocess.run(["git", "-C", str(root), "add", "-A"], check=False)
        nudged = freeze(root)
        assert any("polyzymd hash-trajectories --study ." in w for w in nudged.warnings)
        assert str(tmp_path) not in (root / "manifest.json").read_text()
        _engine(tmp_path, "openmm").record_trajectory_hashes(working, 1)
        assert not any("hash-trajectories" in w for w in freeze(root).warnings)
