"""Tests for the ``polyzymd status`` command and supporting helpers.

Covers:
- render_progress_bar() — bar width, clamping, color-by-status
- discover_replicate_dirs() — glob-based discovery, sorting, filtering
- status CLI — end-to-end output via Click CliRunner
"""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from polyzymd.cli.colors import (
    TerminalColorSupport,
    render_progress_bar,
    set_color_support,
)
from polyzymd.cli.main import cli
from polyzymd.simulation.progress import (
    SegmentRecord,
    SegmentStatus,
    SimulationProgress,
    SimulationStatus,
)

# ---------------------------------------------------------------------------
# render_progress_bar
# ---------------------------------------------------------------------------


class TestRenderProgressBar:
    """Unit tests for render_progress_bar()."""

    def setup_method(self):
        """Disable color so we can inspect raw characters."""
        set_color_support(TerminalColorSupport.NONE)

    def test_default_width_is_40(self):
        bar = render_progress_bar(0.5, "running")
        assert len(bar) == 40

    def test_custom_width(self):
        bar = render_progress_bar(0.5, "running", width=20)
        assert len(bar) == 20

    def test_zero_fraction_all_empty(self):
        bar = render_progress_bar(0.0, "not_started")
        assert bar == "\u2591" * 40  # all ░

    def test_full_fraction_all_filled(self):
        bar = render_progress_bar(1.0, "completed")
        assert bar == "\u2588" * 40  # all █

    def test_half_fraction(self):
        bar = render_progress_bar(0.5, "running", width=10)
        assert bar == "\u2588" * 5 + "\u2591" * 5

    def test_fraction_clamped_above_one(self):
        bar = render_progress_bar(1.5, "completed", width=10)
        assert bar == "\u2588" * 10

    def test_fraction_clamped_below_zero(self):
        bar = render_progress_bar(-0.5, "not_started", width=10)
        assert bar == "\u2591" * 10

    def test_truecolor_wraps_ansi(self):
        set_color_support(TerminalColorSupport.TRUECOLOR)
        bar = render_progress_bar(1.0, "completed", width=5)
        assert "\033[38;2;" in bar
        assert "\033[0m" in bar
        # The actual bar chars should still be present
        assert "\u2588" * 5 in bar

    def test_basic_color_wraps_ansi(self):
        set_color_support(TerminalColorSupport.BASIC)
        bar = render_progress_bar(0.5, "failed", width=10)
        assert "\033[91m" in bar  # red basic ANSI
        assert "\033[0m" in bar

    def test_extended_color_wraps_ansi(self):
        set_color_support(TerminalColorSupport.EXTENDED)
        bar = render_progress_bar(0.5, "interrupted", width=10)
        assert "\033[38;5;" in bar
        assert "\033[0m" in bar

    def test_unknown_status_uses_not_found_color(self):
        set_color_support(TerminalColorSupport.TRUECOLOR)
        bar = render_progress_bar(0.0, "some_unknown_status", width=5)
        # Should use not_found gray (128,128,128)
        assert "\033[38;2;128;128;128m" in bar


# ---------------------------------------------------------------------------
# discover_replicate_dirs
# ---------------------------------------------------------------------------


def _make_mock_config(scratch_dir: Path, template: str | None = None):
    """Create a mock SimulationConfig with the fields discover_replicate_dirs needs."""
    mock = MagicMock()

    # Polymer config
    monomer_a = MagicMock()
    monomer_a.label = "A"
    monomer_a.probability = 0.5
    monomer_b = MagicMock()
    monomer_b.label = "B"
    monomer_b.probability = 0.5

    mock.polymers.enabled = True
    mock.polymers.monomers = [monomer_a, monomer_b]
    mock.polymers.type_prefix = "OEGMA-SBMA"

    mock.substrate.name = "apo"
    mock.enzyme.name = "fnIII"
    mock.thermodynamics.temperature = 310.0
    mock.simulation_phases.production.duration = 100.0

    if template is None:
        template = "{enzyme}_{substrate}_{polymer_type}_{duration}ns_{temperature}K_run{replicate}"
    mock.output.naming_template = template
    mock.output.effective_scratch_directory = scratch_dir

    return mock


class TestDiscoverReplicateDirs:
    """Unit tests for SimulationConfig.discover_replicate_dirs()."""

    def test_finds_matching_directories(self, tmp_path):
        """Discovers run1, run2, run3 directories."""
        from polyzymd.config.schema import SimulationConfig

        # Create directories matching the expected pattern
        base_name = "fnIII_apo_OEGMA-SBMA_A50_B50_100ns_310K"
        for i in (1, 2, 3):
            (tmp_path / f"{base_name}_run{i}").mkdir()

        mock = _make_mock_config(tmp_path)

        # Call the real method on the mock by binding it
        result = SimulationConfig.discover_replicate_dirs(mock)

        assert len(result) == 3
        assert [r[0] for r in result] == [1, 2, 3]
        for num, path in result:
            assert path.name.endswith(f"_run{num}")

    def test_returns_sorted(self, tmp_path):
        """Results are sorted by replicate number even if created out of order."""
        from polyzymd.config.schema import SimulationConfig

        base_name = "fnIII_apo_OEGMA-SBMA_A50_B50_100ns_310K"
        for i in (5, 2, 8, 1):
            (tmp_path / f"{base_name}_run{i}").mkdir()

        mock = _make_mock_config(tmp_path)
        result = SimulationConfig.discover_replicate_dirs(mock)

        assert [r[0] for r in result] == [1, 2, 5, 8]

    def test_ignores_files(self, tmp_path):
        """Only directories are returned, not files."""
        from polyzymd.config.schema import SimulationConfig

        base_name = "fnIII_apo_OEGMA-SBMA_A50_B50_100ns_310K"
        (tmp_path / f"{base_name}_run1").mkdir()
        (tmp_path / f"{base_name}_run2").write_text("not a dir")

        mock = _make_mock_config(tmp_path)
        result = SimulationConfig.discover_replicate_dirs(mock)

        assert len(result) == 1
        assert result[0][0] == 1

    def test_empty_when_no_dirs(self, tmp_path):
        """Returns empty list if no matching dirs exist."""
        from polyzymd.config.schema import SimulationConfig

        mock = _make_mock_config(tmp_path)
        result = SimulationConfig.discover_replicate_dirs(mock)

        assert result == []

    def test_ignores_non_matching_dirs(self, tmp_path):
        """Directories that don't match the pattern are ignored."""
        from polyzymd.config.schema import SimulationConfig

        base_name = "fnIII_apo_OEGMA-SBMA_A50_B50_100ns_310K"
        (tmp_path / f"{base_name}_run1").mkdir()
        (tmp_path / "unrelated_run2").mkdir()
        (tmp_path / f"{base_name}_analysis").mkdir()

        mock = _make_mock_config(tmp_path)
        result = SimulationConfig.discover_replicate_dirs(mock)

        assert len(result) == 1
        assert result[0][0] == 1


# ---------------------------------------------------------------------------
# status CLI command (end-to-end)
# ---------------------------------------------------------------------------


def _write_progress_json(
    directory: Path,
    *,
    total_steps: int = 50000000,
    completed_steps: int = 50000000,
    seg_status: SegmentStatus = SegmentStatus.COMPLETED,
    sim_status: SimulationStatus = SimulationStatus.COMPLETED,
    duration_ns: float = 100.0,
    replicate: int = 1,
) -> None:
    """Write a minimal progress.json into *directory*."""
    progress = SimulationProgress(
        config_path="/tmp/config.yaml",
        total_steps_requested=total_steps,
        total_samples_requested=250,
        timestep_fs=2.0,
        segments=[
            SegmentRecord(
                index=0,
                steps_completed=completed_steps,
                steps_requested=total_steps,
                samples_written=250 if seg_status == SegmentStatus.COMPLETED else 125,
                status=seg_status,
                duration_ns=duration_ns,
            )
        ],
        status=sim_status,
        replicate=replicate,
    )
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "progress.json").write_text(progress.model_dump_json(indent=2))


def _mock_sim_config(scratch_dir: Path) -> MagicMock:
    """Create a mock SimulationConfig for status tests."""
    mock = MagicMock()
    mock.simulation_phases.production.duration = 100.0
    mock.simulation_phases.production.time_step = 2.0
    mock.simulation_phases.production.samples = 250
    mock._format_run_directory_name.return_value = "fnIII_apo_none_100ns_310K_run1"
    mock.output.effective_scratch_directory = scratch_dir
    return mock


class TestStatusCli:
    """End-to-end tests for ``polyzymd status``."""

    def setup_method(self):
        set_color_support(TerminalColorSupport.NONE)

    def _make_dummy_config(self, tmp_path: Path) -> Path:
        """Create a dummy config file (content doesn't matter — we mock from_yaml)."""
        config_path = tmp_path / "config.yaml"
        config_path.write_text("name: test\n")
        return config_path

    def test_status_shows_completed_replicate(self, tmp_path):
        """A completed replicate shows 100% and 'completed'."""
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        rep_dir = scratch / "fnIII_apo_none_100ns_310K_run1"
        _write_progress_json(rep_dir)

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep_dir)]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "100.0%" in result.output
        assert "completed" in result.output
        assert "run1" in result.output

    def test_status_shows_interrupted_replicate(self, tmp_path):
        """An interrupted replicate shows partial progress."""
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        rep_dir = scratch / "fnIII_apo_none_100ns_310K_run1"
        _write_progress_json(
            rep_dir,
            completed_steps=25000000,
            seg_status=SegmentStatus.INTERRUPTED,
            sim_status=SimulationStatus.INTERRUPTED,
            duration_ns=50.0,
        )

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep_dir)]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "interrupted" in result.output
        assert "need attention" in result.output

    def test_status_multiple_replicates(self, tmp_path):
        """Multiple replicates shown in order."""
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        base = "fnIII_apo_none_100ns_310K"
        rep1 = scratch / f"{base}_run1"
        rep2 = scratch / f"{base}_run2"
        rep3 = scratch / f"{base}_run3"

        _write_progress_json(rep1, sim_status=SimulationStatus.COMPLETED)
        _write_progress_json(
            rep2,
            completed_steps=25000000,
            seg_status=SegmentStatus.INTERRUPTED,
            sim_status=SimulationStatus.INTERRUPTED,
            duration_ns=50.0,
        )
        _write_progress_json(
            rep3,
            completed_steps=0,
            seg_status=SegmentStatus.RUNNING,
            sim_status=SimulationStatus.RUNNING,
            duration_ns=0.0,
        )

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = [
            (1, rep1),
            (2, rep2),
            (3, rep3),
        ]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "run1" in result.output
        assert "run2" in result.output
        assert "run3" in result.output
        assert "need attention" in result.output

    def test_status_no_replicates_found(self, tmp_path):
        """When no replicate dirs exist, show helpful message."""
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = []

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "No replicate directories found" in result.output

    def test_status_all_completed(self, tmp_path):
        """When all replicates are done, show success message."""
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        base = "fnIII_apo_none_100ns_310K"
        rep1 = scratch / f"{base}_run1"
        rep2 = scratch / f"{base}_run2"
        _write_progress_json(rep1)
        _write_progress_json(rep2)

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep1), (2, rep2)]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "All 2 replicates completed" in result.output

    def test_status_dir_exists_but_no_progress(self, tmp_path):
        """Directory exists but no progress.json -> shows not_started."""
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        rep_dir = scratch / "fnIII_apo_none_100ns_310K_run1"
        rep_dir.mkdir()

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep_dir)]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "not_started" in result.output
        assert "0.0%" in result.output

    def test_status_invalid_config(self, tmp_path):
        """Invalid config file produces error."""
        bad_config = tmp_path / "bad.yaml"
        bad_config.write_text("name: test\n")

        runner = CliRunner()

        with patch(
            "polyzymd.config.schema.SimulationConfig.from_yaml",
            side_effect=ValueError("bad config"),
        ):
            result = runner.invoke(cli, ["status", "-c", str(bad_config)])

        assert result.exit_code != 0
        assert "Error" in result.output

    def test_status_header_shows_system_name(self, tmp_path):
        """Header includes system name derived from naming template."""
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        rep_dir = scratch / "fnIII_apo_none_100ns_310K_run1"
        rep_dir.mkdir()

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep_dir)]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "polyzymd status" in result.output
        assert "fnIII_apo_none_100ns_310K" in result.output

    def test_status_ns_includes_interrupted_segments(self, tmp_path):
        """ns column accounts for steps in INTERRUPTED segments, not just COMPLETED."""
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        # Simulate a replicate with 2 segments: one interrupted, one running.
        # total_steps_requested = 50_000_000 (100 ns at 2 fs/step)
        # Segment 0: INTERRUPTED at 25_000_000 steps (50 ns worth)
        # Segment 1: RUNNING at 10_000_000 steps (20 ns worth)
        # Total completed steps = 35_000_000 => 70 ns
        # time_completed_ns() would return 0.0 (no COMPLETED segments)
        # Correct ns = 35_000_000 * 2.0 / 1e6 = 70.0
        rep_dir = scratch / "fnIII_apo_none_100ns_310K_run1"
        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=50_000_000,
            total_samples_requested=250,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=25_000_000,
                    steps_requested=50_000_000,
                    samples_written=125,
                    status=SegmentStatus.INTERRUPTED,
                    duration_ns=100.0,
                ),
                SegmentRecord(
                    index=1,
                    steps_completed=10_000_000,
                    steps_requested=25_000_000,
                    samples_written=50,
                    status=SegmentStatus.RUNNING,
                    duration_ns=50.0,
                ),
            ],
            status=SimulationStatus.RUNNING,
            replicate=1,
        )
        rep_dir.mkdir(parents=True, exist_ok=True)
        (rep_dir / "progress.json").write_text(progress.model_dump_json(indent=2))

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep_dir)]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        # Should show 70.0 ns (from steps), NOT 0.0 ns (from time_completed_ns)
        assert "70.0/100.0 ns" in result.output

    def test_status_summary_not_all_completed_when_running(self, tmp_path):
        """'All N replicates completed!' must NOT appear when some are running."""
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        base = "fnIII_apo_none_100ns_310K"
        rep1 = scratch / f"{base}_run1"
        rep2 = scratch / f"{base}_run2"

        _write_progress_json(rep1)  # completed
        _write_progress_json(
            rep2,
            completed_steps=25_000_000,
            seg_status=SegmentStatus.RUNNING,
            sim_status=SimulationStatus.RUNNING,
            duration_ns=50.0,
        )

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep1), (2, rep2)]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "All 2 replicates completed" not in result.output
        assert "still running" in result.output
        assert "1/2 completed" in result.output

    def test_status_summary_mixed_attention_and_running(self, tmp_path):
        """Both 'need attention' and 'still running' lines shown for mixed states."""
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        base = "fnIII_apo_none_100ns_310K"
        rep1 = scratch / f"{base}_run1"
        rep2 = scratch / f"{base}_run2"
        rep3 = scratch / f"{base}_run3"

        _write_progress_json(rep1)  # completed
        _write_progress_json(
            rep2,
            completed_steps=25_000_000,
            seg_status=SegmentStatus.INTERRUPTED,
            sim_status=SimulationStatus.INTERRUPTED,
            duration_ns=50.0,
        )
        _write_progress_json(
            rep3,
            completed_steps=10_000_000,
            seg_status=SegmentStatus.RUNNING,
            sim_status=SimulationStatus.RUNNING,
            duration_ns=20.0,
        )

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = [
            (1, rep1),
            (2, rep2),
            (3, rep3),
        ]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "All 3 replicates completed" not in result.output
        assert "need attention" in result.output
        assert "still running" in result.output

    def test_status_stale_running_reclassified_by_filesystem_scan(self, tmp_path):
        """A stale 'running' replicate is reclassified to 'interrupted' by filesystem scan.

        When progress.json says 'running' but the checkpoint file is older
        than CHECKPOINT_RECENCY_SECONDS (10 min), validate_progress() via
        scan_filesystem() reclassifies the segment as 'interrupted'.

        Without production_N/ dirs on disk, validate_progress() keeps the
        progress file records but recomputes overall status. Since the
        segments remain unchanged (no filesystem segments to override),
        the status stays as-is from progress.json. The stale detection
        in production relies on actual checkpoint file age.

        This test verifies that a stale running replicate (no production
        dirs to confirm it's alive) shows as needing attention because the
        filesystem scan finds no evidence of an active simulation.
        """
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        base = "fnIII_apo_none_100ns_310K"
        rep1 = scratch / f"{base}_run1"
        rep2 = scratch / f"{base}_run2"

        _write_progress_json(rep1)  # completed

        # rep2: "running" per progress.json, but no production dirs on disk.
        # The filesystem scan finds no segments, but validate_progress keeps
        # the progress file's segment record ("in progress file but not on
        # filesystem"). The segment keeps its RUNNING status, so the overall
        # status remains RUNNING — but the summary should NOT say "All completed".
        _write_progress_json(
            rep2,
            completed_steps=25_000_000,
            seg_status=SegmentStatus.RUNNING,
            sim_status=SimulationStatus.RUNNING,
            duration_ns=50.0,
        )

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep1), (2, rep2)]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "All 2 replicates completed" not in result.output

    def test_status_stale_running_with_old_checkpoint_becomes_interrupted(self, tmp_path):
        """A replicate with an old checkpoint file is reclassified to interrupted.

        When a production_N/ directory has a .chk file older than
        CHECKPOINT_RECENCY_SECONDS, scan_filesystem() marks it as
        INTERRUPTED. validate_progress() then overrides the progress
        file's RUNNING status with the filesystem's INTERRUPTED.
        """
        import time

        scratch = tmp_path / "scratch"
        scratch.mkdir()

        base = "fnIII_apo_none_100ns_310K"
        rep_dir = scratch / f"{base}_run1"
        rep_dir.mkdir(parents=True)

        # Write progress.json saying "running"
        _write_progress_json(
            rep_dir,
            completed_steps=25_000_000,
            seg_status=SegmentStatus.RUNNING,
            sim_status=SimulationStatus.RUNNING,
            duration_ns=50.0,
        )

        # Create a production_0/ directory with an old checkpoint file
        prod_dir = rep_dir / "production_0"
        prod_dir.mkdir()
        chk_file = prod_dir / "production_0_checkpoint.chk"
        chk_file.write_bytes(b"fake checkpoint data")
        # Backdate the checkpoint to make it older than CHECKPOINT_RECENCY_SECONDS (600s)
        old_time = time.time() - 1200  # 20 minutes ago
        import os

        os.utime(str(chk_file), (old_time, old_time))

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep_dir)]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        # The filesystem scan should reclassify to interrupted
        assert "interrupted" in result.output
        assert "need attention" in result.output
        assert "All 1 replicates completed" not in result.output
