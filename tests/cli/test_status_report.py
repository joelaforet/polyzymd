"""Tests for the agent-oriented ``polyzymd status --format agent|json`` path."""

from __future__ import annotations

import json
from datetime import datetime, timedelta, timezone
from pathlib import Path
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from polyzymd.cli.colors import TerminalColorSupport, set_color_support
from polyzymd.cli.main import cli
from polyzymd.cli.status_report import (
    ReplicateReport,
    SlurmJob,
    SystemReport,
    classify,
    discover_config_files,
    estimate_eta_days,
    estimate_rate_ns_per_day,
    find_latest_log,
    jobs_by_name,
    last_error_line,
    parse_squeue_output,
    render_agent,
)
from polyzymd.simulation.progress import (
    SegmentRecord,
    SegmentStatus,
    SimulationProgress,
    SimulationStatus,
)

NOW = datetime(2026, 9, 11, 16, 0, tzinfo=timezone.utc)


# ---------------------------------------------------------------------------
# squeue parsing
# ---------------------------------------------------------------------------


def test_parse_squeue_output_fields_and_grouping():
    text = (
        "28248388|CALB_none_run5|R|3:27|bgpu-shirts2|None\n"
        "28248422|CALB_none_run3|PD|0:00||(Priority)\n"
        "\n"
    )
    jobs = parse_squeue_output(text)
    assert [j.job_id for j in jobs] == ["28248388", "28248422"]
    assert jobs[0].node == "bgpu-shirts2"
    assert jobs[1].reason == "(Priority)"
    grouped = jobs_by_name(jobs)
    assert set(grouped) == {"CALB_none_run5", "CALB_none_run3"}


# ---------------------------------------------------------------------------
# classification
# ---------------------------------------------------------------------------


def test_classify_verdicts():
    running = [SlurmJob("1", "n", "R")]
    pending = [SlurmJob("2", "n", "PD")]
    assert classify("completed", [], True) == "completed"
    assert classify("interrupted", running, True) == "running"
    assert classify("interrupted", pending, True) == "queued"
    assert classify("interrupted", [], True) == "dead"
    assert classify("running", [], True) == "dead"  # stale "running" with no job
    assert classify("not_started", [], True) == "not_started"
    assert classify("interrupted", running, False) == "not_found"


# ---------------------------------------------------------------------------
# throughput / ETA
# ---------------------------------------------------------------------------


def _progress(segments, timestep_fs=2.0):
    return SimulationProgress(
        config_path="c.yaml",
        total_steps_requested=500_000_000,
        total_samples_requested=2500,
        timestep_fs=timestep_fs,
        segments=segments,
        status=SimulationStatus.INTERRUPTED,
    )


def test_rate_prefers_finished_segment():
    start = NOW - timedelta(days=2)
    seg0 = SegmentRecord(
        index=0,
        steps_completed=100_000_000,  # 200 ns
        started_at=start.isoformat(),
        finished_at=(start + timedelta(days=1)).isoformat(),
        status=SegmentStatus.COMPLETED,
    )
    seg1 = SegmentRecord(
        index=1,
        steps_completed=10_000,
        started_at=(NOW - timedelta(minutes=1)).isoformat(),
        status=SegmentStatus.RUNNING,
    )
    rate = estimate_rate_ns_per_day(_progress([seg0, seg1]), NOW, live=True)
    assert rate is not None and abs(rate - 200.0) < 1e-6


def test_rate_falls_back_to_live_segment_only_when_live():
    seg = SegmentRecord(
        index=0,
        steps_completed=50_000_000,  # 100 ns
        started_at=(NOW - timedelta(hours=12)).isoformat(),
        status=SegmentStatus.RUNNING,
    )
    assert estimate_rate_ns_per_day(_progress([seg]), NOW, live=False) is None
    rate = estimate_rate_ns_per_day(_progress([seg]), NOW, live=True)
    assert rate is not None and abs(rate - 200.0) < 1e-6


def test_rate_ignores_tiny_windows():
    seg = SegmentRecord(
        index=0,
        steps_completed=5_000_000,
        started_at=(NOW - timedelta(minutes=2)).isoformat(),
        status=SegmentStatus.RUNNING,
    )
    assert estimate_rate_ns_per_day(_progress([seg]), NOW, live=True) is None


def test_eta():
    assert estimate_eta_days(600.0, 200.0) == 3.0
    assert estimate_eta_days(600.0, None) is None
    assert estimate_eta_days(0.0, 200.0) == 0.0


# ---------------------------------------------------------------------------
# log inspection
# ---------------------------------------------------------------------------


def test_find_latest_log_and_error_line(tmp_path: Path):
    logs = tmp_path / "slurm_logs"
    logs.mkdir()
    old = logs / "sys_run1.100.out"
    old.write_text("run-segment exited with code 0\n")
    new = logs / "sys_run1.200.out"
    new.write_text(
        "ROUTING: sim-cuda-12-4 is incompatible with driver 525.147\n"
        "FATAL: CUDA routing failed after 3 retries\n"
    )
    import os
    import time

    t = time.time()
    os.utime(old, (t - 100, t - 100))
    os.utime(new, (t, t))

    assert find_latest_log(logs, "sys_run1") == new
    assert last_error_line(new) == "FATAL: CUDA routing failed after 3 retries"
    assert last_error_line(old) is None
    assert find_latest_log(logs, "other") is None
    assert find_latest_log(None, "sys_run1") is None


def test_error_line_concurrent_and_segment_failed(tmp_path: Path):
    p = tmp_path / "a.out"
    p.write_text(
        "Segment 0 failed: Particle coordinate is NaN.\n"
        "run-segment exited with code 1\n"
        "FATAL: run-segment failed (exit code 1) — NOT resubmitting\n"
    )
    assert last_error_line(p).startswith("FATAL: run-segment failed")
    p.write_text("run-segment exited with code 2\nCONCURRENT: Another job is already running\n")
    assert last_error_line(p).startswith("CONCURRENT:")


# ---------------------------------------------------------------------------
# rendering
# ---------------------------------------------------------------------------


def _rep(n, verdict, ns=100.0, **kw):
    return ReplicateReport(
        replicate=n,
        directory=f"/scratch/sys_run{n}",
        progress_status="interrupted" if verdict != "completed" else "completed",
        completed_ns=ns,
        total_ns=1000.0,
        fraction=ns / 1000.0,
        verdict=verdict,
        **kw,
    )


def test_render_agent_lines_and_resubmit_hint():
    report = SystemReport(
        name="SYS",
        config_path="a/config.yaml",
        scratch_directory="/scratch",
        replicates=[
            _rep(1, "completed", ns=1000.0),
            _rep(
                2,
                "running",
                ns=250.0,
                jobs=[SlurmJob("11", "SYS_run2", "R", "3:00:00", "node1")],
                rate_ns_per_day=200.0,
                eta_days=3.75,
            ),
            _rep(
                3, "queued", ns=50.0, jobs=[SlurmJob("12", "SYS_run3", "PD", "", "", "(Priority)")]
            ),
            _rep(
                4,
                "dead",
                ns=112.3,
                last_error="FATAL: CUDA routing failed after 3 retries",
                last_log="SYS_run4.99.out",
            ),
        ],
    )
    text = render_agent([report], now=NOW, preset_hint="blanca-shirts")
    assert "1 system(s)  4 replicate(s): 1 completed, 1 running, 1 queued, 1 dead" in text
    assert (
        "run2   250.0/1000ns   25%  RUNNING      job 11 R 3:00:00 node1  200ns/d  eta 3.8d" in text
    )
    assert "run3    50.0/1000ns    5%  QUEUED       job 12 PD ((Priority))  eta ?" in text
    assert (
        "DEAD         no job  last: FATAL: CUDA routing failed after 3 retries [SYS_run4.99.out]"
        in text
    )
    assert "polyzymd submit -c a/config.yaml -r 4 --preset blanca-shirts" in text
    assert "\x1b[" not in text  # no colour codes


def test_render_agent_warns_when_slurm_unavailable():
    report = SystemReport("S", "c.yaml", "/s", [_rep(1, "dead", note="slurm unavailable")])
    text = render_agent([report], now=NOW, slurm_available=False)
    assert "WARNING: squeue unavailable" in text
    assert "slurm unavailable" in text


# ---------------------------------------------------------------------------
# config discovery
# ---------------------------------------------------------------------------


def test_discover_config_files_depth_limited(tmp_path: Path):
    (tmp_path / "a").mkdir()
    (tmp_path / "a" / "config.yaml").write_text("x")
    deep = tmp_path / "b" / "c" / "d" / "e"
    deep.mkdir(parents=True)
    (deep / "config.yaml").write_text("x")
    (tmp_path / ".hidden").mkdir()
    (tmp_path / ".hidden" / "config.yaml").write_text("x")
    found = discover_config_files([tmp_path])
    assert found == [tmp_path / "a" / "config.yaml"]
    assert discover_config_files([tmp_path / "a" / "config.yaml"]) == [
        tmp_path / "a" / "config.yaml"
    ]


# ---------------------------------------------------------------------------
# CLI end-to-end (squeue patched)
# ---------------------------------------------------------------------------


def _mock_cfg(scratch: Path, logs: Path):
    cfg = MagicMock()
    cfg.simulation_phases.production.duration = 1000.0
    cfg.output.effective_scratch_directory = scratch
    cfg.output.get_slurm_logs_directory.return_value = logs
    cfg._format_run_directory_name.side_effect = lambda r: f"SYS_run{r}"
    cfg.format_run_directory_name.side_effect = lambda replicate=1: f"SYS_run{replicate}"
    cfg.engine = "openmm"
    cfg.simulation_engine = None
    return cfg


def _write_progress(rep_dir: Path, steps: int, status: SimulationStatus, seg_status):
    rep_dir.mkdir(parents=True, exist_ok=True)
    start = NOW - timedelta(days=1)
    prog = SimulationProgress(
        config_path="c.yaml",
        total_steps_requested=500_000_000,
        total_samples_requested=2500,
        timestep_fs=2.0,
        segments=[
            SegmentRecord(
                index=0,
                steps_completed=steps,
                steps_requested=500_000_000,
                started_at=start.isoformat(),
                finished_at=(start + timedelta(hours=12)).isoformat(),
                status=seg_status,
            )
        ],
        status=status,
    )
    (rep_dir / "progress.json").write_text(prog.model_dump_json())
    return prog


class TestStatusAgentCli:
    def setup_method(self):
        set_color_support(TerminalColorSupport.NONE)

    def test_agent_format_joins_squeue_and_logs(self, tmp_path: Path):
        scratch = tmp_path / "scratch"
        logs = tmp_path / "slurm_logs"
        logs.mkdir()
        rep1 = scratch / "SYS_run1"
        rep2 = scratch / "SYS_run2"
        p1 = _write_progress(
            rep1, 125_000_000, SimulationStatus.INTERRUPTED, SegmentStatus.INTERRUPTED
        )
        p2 = _write_progress(
            rep2, 60_000_000, SimulationStatus.INTERRUPTED, SegmentStatus.INTERRUPTED
        )
        (logs / "SYS_run2.77.out").write_text("FATAL: CUDA routing failed after 3 retries\n")

        cfg = _mock_cfg(scratch, logs)
        cfg.discover_replicate_dirs.return_value = [(1, rep1), (2, rep2)]
        config_path = tmp_path / "config.yaml"
        config_path.write_text("name: x\n")

        engine = MagicMock()
        engine.resolve_engine_working_directory.side_effect = lambda p: p
        engine.load_or_scan_progress.side_effect = lambda d, r: {rep1: p1, rep2: p2}[d]

        squeue_text = "555|SYS_run1|R|2:00:00|nodeA|None\n"
        runner = CliRunner()
        with (
            patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=cfg),
            patch("polyzymd.engines.create_engine", return_value=engine),
            patch("polyzymd.cli.main._resolve_engine_name", return_value="openmm"),
            patch("polyzymd.cli.main.warn_if_wrong_pixi_env"),
            patch("polyzymd.simulation.progress.save_progress"),
            patch(
                "polyzymd.workflow.daisy_chain.create_job_name",
                side_effect=lambda c, r: f"SYS_run{r}",
            ),
            patch("polyzymd.cli.status_report.subprocess.run") as run,
        ):
            run.return_value = MagicMock(returncode=0, stdout=squeue_text, stderr="")
            result = runner.invoke(
                cli, ["status", "--format", "agent", "-c", str(config_path), "--preset", "p"]
            )

        assert result.exit_code == 0, result.output
        out = result.output
        assert "1 running, 1 dead" in out
        assert (
            "run1   250.0/1000ns   25%  RUNNING      job 555 R 2:00:00 nodeA  500ns/d  eta 1.5d"
            in out
        )
        assert (
            "run2   120.0/1000ns   12%  DEAD         no job  last: FATAL: CUDA routing failed after 3 retries [SYS_run2.77.out]"
            in out
        )
        assert f"polyzymd submit -c {config_path} -r 2 --preset p" in out
        # exactly one squeue call regardless of replicate count
        assert run.call_count == 1

    def test_json_format(self, tmp_path: Path):
        scratch = tmp_path / "scratch"
        rep1 = scratch / "SYS_run1"
        p1 = _write_progress(rep1, 500_000_000, SimulationStatus.COMPLETED, SegmentStatus.COMPLETED)
        cfg = _mock_cfg(scratch, tmp_path / "logs")
        cfg.discover_replicate_dirs.return_value = [(1, rep1)]
        config_path = tmp_path / "config.yaml"
        config_path.write_text("name: x\n")
        engine = MagicMock()
        engine.resolve_engine_working_directory.side_effect = lambda p: p
        engine.load_or_scan_progress.return_value = p1

        with (
            patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=cfg),
            patch("polyzymd.engines.create_engine", return_value=engine),
            patch("polyzymd.cli.main._resolve_engine_name", return_value="openmm"),
            patch("polyzymd.cli.main.warn_if_wrong_pixi_env"),
            patch("polyzymd.simulation.progress.save_progress"),
            patch(
                "polyzymd.workflow.daisy_chain.create_job_name",
                side_effect=lambda c, r: f"SYS_run{r}",
            ),
        ):
            result = CliRunner().invoke(
                cli, ["status", "--format", "json", "--no-slurm", "-c", str(config_path)]
            )
        assert result.exit_code == 0, result.output
        data = json.loads(result.output)
        assert data["slurm_available"] is False
        assert data["systems"][0]["counts"] == {"completed": 1}
        assert data["systems"][0]["replicates"][0]["verdict"] == "completed"

    def test_table_rejects_multiple_configs(self, tmp_path: Path):
        a = tmp_path / "a.yaml"
        b = tmp_path / "b.yaml"
        a.write_text("x")
        b.write_text("x")
        result = CliRunner().invoke(cli, ["status", "-c", str(a), "-c", str(b)])
        assert result.exit_code != 0
        assert "one config at a time" in result.output

    def test_requires_config_or_all(self):
        result = CliRunner().invoke(cli, ["status"])
        assert result.exit_code != 0
        assert "at least one" in result.output
