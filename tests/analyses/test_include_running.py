"""``include_running``: read production segments still being written, or leave the replicate out.

A replicate whose only production segment ``progress.json`` records as
running has no finished frames. Without ``include_running`` it is left out of
its condition with a warning and listed in the report; with it, the frames
written so far are read and the report says the values are provisional.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.simulation.progress import (
    SegmentRecord,
    SegmentStatus,
    SimulationProgress,
    SimulationStatus,
    save_progress,
)
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config


def _mark(run_dir: Path, status: SegmentStatus, frames: int = 10) -> None:
    """Record segment 0 of ``run_dir`` with ``status`` in its ``progress.json``."""
    done = status == SegmentStatus.COMPLETED
    save_progress(
        run_dir,
        SimulationProgress(
            total_steps_requested=500_000,
            total_samples_requested=10,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=frames * 50_000,
                    steps_requested=500_000,
                    samples_written=frames,
                    status=status,
                    duration_ns=frames * 0.1,
                )
            ],
            status=SimulationStatus.COMPLETED if done else SimulationStatus.RUNNING,
        ),
    )


def _config(tmp_path: Path, statuses: list[SegmentStatus]) -> Path:
    """Write one condition whose replicate ``i`` has one segment with ``statuses[i - 1]``."""
    pytest.importorskip("MDAnalysis")
    config = write_simulation_config(tmp_path / "noPoly", scratch=tmp_path / "scratch")
    (tmp_path / "noPoly" / "test.pdb").write_text("REMARK input\nEND\n")
    for replicate, status in enumerate(statuses, start=1):
        scales = [1.0 + 0.1 * replicate + 0.01 * k for k in range(10)]
        _mark(write_openmm_replicate(config, replicate, scales), status)
    return config


def _study(config: Path, **options):
    from polyzymd.analyses.study import Study

    return Study.from_configs({"No polymer": config}, equilibration="0ns", **options)


def test_a_replicate_with_only_a_running_segment_is_left_out(tmp_path: Path, caplog) -> None:
    """The other replicates are analysed; the running one is named, with the flag that reads it."""
    config = _config(tmp_path, [SegmentStatus.COMPLETED, SegmentStatus.RUNNING])
    with caplog.at_level(logging.WARNING, logger="polyzymd.analyses.study"):
        condition = _study(config)["No polymer"]
    assert [replicate.index for replicate in condition.replicates] == [1]
    assert condition.only_running == [2]
    assert "--include-running" in caplog.text


def test_every_replicate_running_says_how_to_read_them(tmp_path: Path) -> None:
    """With no finished frames anywhere, the error gives --include-running, not a missing file."""
    config = _config(tmp_path, [SegmentStatus.RUNNING, SegmentStatus.RUNNING])
    with pytest.raises(ProtocolError, match="still writing the first production segment") as error:
        _study(config)
    assert "--include-running" in error.value.hint


def test_include_running_reads_the_running_segment_and_marks_the_record(tmp_path: Path) -> None:
    """Both replicates are read; the stored-result identity says the frames were in flight."""
    config = _config(tmp_path, [SegmentStatus.COMPLETED, SegmentStatus.RUNNING])
    condition = _study(config, include_running=True)["No polymer"]
    assert [replicate.index for replicate in condition.replicates] == [1, 2]
    assert condition.only_running == []
    assert len(condition.replicates[1].frames) == 10
    assert condition.replicates[1].identity["include_running"] is True
    assert "include_running" not in _study(config)["No polymer"].replicates[0].identity


def test_the_report_lists_the_left_out_replicate_or_says_it_is_provisional(tmp_path: Path) -> None:
    """Without the flag the report is partial and names replicate 2; with it, it is provisional."""
    from polyzymd.analyses.protocols import ONLY_RUNNING_REASON, PROVISIONAL_WARNING, analyze

    config = _config(tmp_path, [SegmentStatus.COMPLETED, SegmentStatus.RUNNING])
    common = {
        "equilibration": "0ns",
        "settings": {"selection": "all"},
        "labels": ["No polymer"],
        "plots": False,
    }
    left_out = analyze("rg", [config], output_dir=tmp_path / "final", **common)
    assert left_out.status == "partial"
    assert left_out.excluded_replicates == {"No polymer": {"2": ONLY_RUNNING_REASON}}
    assert PROVISIONAL_WARNING not in left_out.warnings

    provisional = analyze(
        "rg", [config], output_dir=tmp_path / "live", include_running=True, **common
    )
    assert provisional.excluded_replicates == {}
    assert PROVISIONAL_WARNING in provisional.warnings


@pytest.mark.usefixtures("git_identity")
def test_study_yaml_include_running(tmp_path: Path) -> None:
    """include_running: true is read from study.yaml; anything but true or false is refused."""
    pytest.importorskip("MDAnalysis")
    from polyzymd.analyses.study_file import load_study_file
    from tests._support.analysis_testkit import write_committed_study

    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    assert load_study_file(root).include_running is False
    text = (root / "study.yaml").read_text()
    (root / "study.yaml").write_text(text + "include_running: true\n")
    assert load_study_file(root).include_running is True
    (root / "study.yaml").write_text(text + "include_running: yes please\n")
    with pytest.raises(ProtocolError, match="include_running must be true or false"):
        load_study_file(root)


@pytest.mark.filterwarnings("ignore")
def test_the_command_line_flag_reads_running_segments_and_reaches_submit_tasks(
    tmp_path: Path,
) -> None:
    """analyze --include-running reports provisional values; --submit passes the flag on."""
    from click.testing import CliRunner

    from polyzymd.cli.analyze import EXIT_PARTIAL
    from polyzymd.cli.main import cli

    config = _config(tmp_path, [SegmentStatus.COMPLETED, SegmentStatus.RUNNING])
    base = ["analyze", "rg", "-c", str(config), "--eq", "0ns", "--set", "selection=all"]
    result = CliRunner().invoke(
        cli, [*base, "--include-running", "--no-plots", "--output-dir", str(tmp_path / "live")]
    )
    # Production is not finished, so the report is partial: exit code 3.
    assert result.exit_code == EXIT_PARTIAL, result.output
    assert "provisional" in result.output

    dry = CliRunner().invoke(
        cli, [*base, "--include-running", "--dry-run", "--output-dir", str(tmp_path / "jobs")]
    )
    assert dry.exit_code == 0, dry.output
    scripts = [path.read_text() for path in (tmp_path / "jobs").rglob("*.sbatch")]
    assert scripts and any("--include-running" in text for text in scripts)


def test_include_running_leaves_out_a_replicate_that_has_no_frame_yet(tmp_path: Path) -> None:
    """A run that just started (a segment folder, no frame) is left out, not fatal for the rest."""
    from polyzymd.analyses.protocols import NO_FRAMES_YET_REASON, analyze
    from polyzymd.config.schema import SimulationConfig

    config = _config(tmp_path, [SegmentStatus.COMPLETED, SegmentStatus.RUNNING])
    just_started = SimulationConfig.from_yaml(config).get_working_directory(3)
    (just_started / "production_0").mkdir(parents=True)
    _mark(just_started, SegmentStatus.RUNNING, frames=0)
    condition = _study(config, include_running=True)["No polymer"]
    assert [replicate.index for replicate in condition.replicates] == [1, 2]
    assert condition.no_frames_yet == [3]
    report = analyze(
        "rg",
        [config],
        equilibration="0ns",
        settings={"selection": "all"},
        labels=["No polymer"],
        plots=False,
        output_dir=tmp_path / "live",
        include_running=True,
    )
    assert report.excluded_replicates == {"No polymer": {"3": NO_FRAMES_YET_REASON}}
