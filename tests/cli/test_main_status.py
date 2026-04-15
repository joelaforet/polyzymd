"""Tests for the ``polyzymd status`` command and supporting helpers."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from polyzymd.cli.colors import TerminalColorSupport, set_color_support
from polyzymd.cli.main import cli
from polyzymd.simulation.progress import (
    SegmentRecord,
    SegmentStatus,
    SimulationProgress,
    SimulationStatus,
)


def _make_mock_config(scratch_dir: Path, template: str | None = None):
    """Create a mock SimulationConfig with status helper fields."""
    mock = MagicMock()

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
        from polyzymd.config.schema import SimulationConfig

        base_name = "fnIII_apo_OEGMA-SBMA_A50_B50_100ns_310K"
        for i in (1, 2, 3):
            (tmp_path / f"{base_name}_run{i}").mkdir()

        mock = _make_mock_config(tmp_path)
        result = SimulationConfig.discover_replicate_dirs(mock)

        assert len(result) == 3
        assert [replicate[0] for replicate in result] == [1, 2, 3]
        for number, path in result:
            assert path.name.endswith(f"_run{number}")

    def test_returns_sorted(self, tmp_path):
        from polyzymd.config.schema import SimulationConfig

        base_name = "fnIII_apo_OEGMA-SBMA_A50_B50_100ns_310K"
        for i in (5, 2, 8, 1):
            (tmp_path / f"{base_name}_run{i}").mkdir()

        mock = _make_mock_config(tmp_path)
        result = SimulationConfig.discover_replicate_dirs(mock)

        assert [replicate[0] for replicate in result] == [1, 2, 5, 8]

    def test_ignores_files(self, tmp_path):
        from polyzymd.config.schema import SimulationConfig

        base_name = "fnIII_apo_OEGMA-SBMA_A50_B50_100ns_310K"
        (tmp_path / f"{base_name}_run1").mkdir()
        (tmp_path / f"{base_name}_run2").write_text("not a dir")

        mock = _make_mock_config(tmp_path)
        result = SimulationConfig.discover_replicate_dirs(mock)

        assert len(result) == 1
        assert result[0][0] == 1

    def test_empty_when_no_dirs(self, tmp_path):
        from polyzymd.config.schema import SimulationConfig

        mock = _make_mock_config(tmp_path)
        result = SimulationConfig.discover_replicate_dirs(mock)

        assert result == []

    def test_ignores_non_matching_dirs(self, tmp_path):
        from polyzymd.config.schema import SimulationConfig

        base_name = "fnIII_apo_OEGMA-SBMA_A50_B50_100ns_310K"
        (tmp_path / f"{base_name}_run1").mkdir()
        (tmp_path / "unrelated_run2").mkdir()
        (tmp_path / f"{base_name}_analysis").mkdir()

        mock = _make_mock_config(tmp_path)
        result = SimulationConfig.discover_replicate_dirs(mock)

        assert len(result) == 1
        assert result[0][0] == 1


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
    mock.engine = "openmm"
    mock.simulation_phases.production.duration = 100.0
    mock.simulation_phases.production.time_step = 2.0
    mock.simulation_phases.production.samples = 250
    mock._format_run_directory_name.return_value = "fnIII_apo_none_100ns_310K_run1"
    mock.output.effective_scratch_directory = scratch_dir
    return mock


class TestCheckProgressEngineAware:
    """Tests for engine-aware check-progress command."""

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.engines.openmm.engine.OpenMMEngine.load_or_scan_progress")
    def test_check_progress_uses_engine_dispatch(
        self, mock_load_progress, mock_from_yaml, tmp_path
    ):
        """check-progress should use engine dispatch for progress loading."""
        mock_cfg = _mock_sim_config(tmp_path / "scratch")
        mock_cfg.engine = "openmm"
        mock_from_yaml.return_value = mock_cfg

        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=50_000_000,
            total_samples_requested=250,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=50_000_000,
                    steps_requested=50_000_000,
                    samples_written=250,
                    status=SegmentStatus.COMPLETED,
                    duration_ns=100.0,
                )
            ],
            status=SimulationStatus.COMPLETED,
            replicate=1,
        )
        mock_load_progress.return_value = progress

        config_path = tmp_path / "config.yaml"
        config_path.write_text("name: test\n")
        runner = CliRunner()

        result = runner.invoke(cli, ["check-progress", "-c", str(config_path)])

        assert result.exit_code == 0
        assert "COMPLETE" in result.output
        mock_load_progress.assert_called_once()

    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    @patch("polyzymd.engines.openmm.engine.OpenMMEngine.load_or_scan_progress")
    def test_check_progress_incomplete_exits_1(self, mock_load_progress, mock_from_yaml, tmp_path):
        """check-progress should exit 1 when work remains."""
        mock_cfg = _mock_sim_config(tmp_path / "scratch")
        mock_cfg.engine = "openmm"
        mock_from_yaml.return_value = mock_cfg

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
                    duration_ns=50.0,
                )
            ],
            status=SimulationStatus.INTERRUPTED,
            replicate=1,
        )
        mock_load_progress.return_value = progress

        config_path = tmp_path / "config.yaml"
        config_path.write_text("name: test\n")
        runner = CliRunner()

        result = runner.invoke(cli, ["check-progress", "-c", str(config_path)])

        assert result.exit_code == 1
        assert "remaining" in result.output.lower()


class TestCheckProgressGromacs:
    """Tests for check-progress with GROMACS engine."""

    @patch("polyzymd.engines.gromacs.engine.GromacsEngine.load_or_scan_progress")
    @patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_check_progress_gromacs_complete(
        self,
        mock_from_yaml,
        mock_resolve,
        mock_load,
        tmp_path,
    ):
        """check-progress --engine gromacs should exit 0 when complete."""
        _ = mock_resolve
        mock_cfg = _mock_sim_config(tmp_path / "scratch")
        mock_cfg.engine = "gromacs"
        mock_cfg.gromacs = SimpleNamespace(
            gmx_binary=None,
            grompp_flags="",
            mdrun_flags="",
            module_load=None,
        )
        mock_from_yaml.return_value = mock_cfg

        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=50_000_000,
            total_samples_requested=250,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=50_000_000,
                    steps_requested=50_000_000,
                    samples_written=0,
                    status=SegmentStatus.COMPLETED,
                    duration_ns=100.0,
                )
            ],
            status=SimulationStatus.COMPLETED,
            replicate=1,
        )
        mock_load.return_value = progress

        config_path = tmp_path / "config.yaml"
        config_path.write_text("name: test\n")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["check-progress", "-c", str(config_path), "--engine", "gromacs"],
        )

        assert result.exit_code == 0
        assert "COMPLETE" in result.output
        mock_load.assert_called_once()

    @patch("polyzymd.engines.gromacs.engine.GromacsEngine.load_or_scan_progress")
    @patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_check_progress_gromacs_incomplete_exits_1(
        self,
        mock_from_yaml,
        mock_resolve,
        mock_load,
        tmp_path,
    ):
        """check-progress --engine gromacs should exit 1 when work remains."""
        _ = mock_resolve
        mock_cfg = _mock_sim_config(tmp_path / "scratch")
        mock_cfg.engine = "gromacs"
        mock_cfg.gromacs = SimpleNamespace(
            gmx_binary=None,
            grompp_flags="",
            mdrun_flags="",
            module_load=None,
        )
        mock_from_yaml.return_value = mock_cfg

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
                    samples_written=0,
                    status=SegmentStatus.INTERRUPTED,
                    duration_ns=50.0,
                )
            ],
            status=SimulationStatus.INTERRUPTED,
            replicate=1,
        )
        mock_load.return_value = progress

        config_path = tmp_path / "config.yaml"
        config_path.write_text("name: test\n")

        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["check-progress", "-c", str(config_path), "--engine", "gromacs"],
        )

        assert result.exit_code == 1
        assert "remaining" in result.output.lower()


class TestStatusGromacs:
    """Tests for status command with GROMACS engine."""

    def setup_method(self):
        set_color_support(TerminalColorSupport.NONE)

    @patch("polyzymd.engines.gromacs.engine.GromacsEngine.load_or_scan_progress")
    @patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_status_gromacs_shows_progress(self, mock_from_yaml, mock_resolve, mock_load, tmp_path):
        """status should work with engine=gromacs in config."""
        _ = mock_resolve
        scratch = tmp_path / "scratch"
        rep_dir = scratch / "system_run1"
        rep_dir.mkdir(parents=True)

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.engine = "gromacs"
        mock_cfg.gromacs = SimpleNamespace(
            gmx_binary=None,
            grompp_flags="",
            mdrun_flags="",
            module_load=None,
        )
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep_dir)]
        mock_from_yaml.return_value = mock_cfg

        progress = SimulationProgress(
            config_path="/tmp/config.yaml",
            total_steps_requested=50_000_000,
            total_samples_requested=250,
            timestep_fs=2.0,
            segments=[
                SegmentRecord(
                    index=0,
                    steps_completed=50_000_000,
                    steps_requested=50_000_000,
                    samples_written=0,
                    status=SegmentStatus.COMPLETED,
                    duration_ns=100.0,
                )
            ],
            status=SimulationStatus.COMPLETED,
            replicate=1,
        )
        mock_load.return_value = progress

        config_path = tmp_path / "config.yaml"
        config_path.write_text("name: test\n")

        runner = CliRunner()
        result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0


class TestStatusGromacsReplicateDiscovery:
    """Tests for status command GROMACS replicate directory discovery."""

    def setup_method(self):
        set_color_support(TerminalColorSupport.NONE)

    @patch("polyzymd.engines.gromacs.binary.resolve_gromacs_binary", return_value="gmx")
    @patch("polyzymd.config.schema.SimulationConfig.from_yaml")
    def test_no_fallback_when_scratch_has_replicates(self, mock_from_yaml, mock_resolve, tmp_path):
        """Legacy fallback should NOT be used when scratch has replicates."""
        _ = mock_resolve
        scratch = tmp_path / "scratch"
        rep_dir = scratch / "system_run1"
        rep_dir.mkdir(parents=True)
        gromacs_dir = rep_dir / "gromacs"
        gromacs_dir.mkdir()

        projects = tmp_path / "projects"
        (projects / "replicate_1" / "gromacs").mkdir(parents=True)

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.engine = "gromacs"
        mock_cfg.gromacs = SimpleNamespace(
            gmx_binary=None,
            grompp_flags="",
            mdrun_flags="",
            module_load=None,
        )
        mock_cfg.output.projects_directory = projects
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep_dir)]
        mock_from_yaml.return_value = mock_cfg

        _write_progress_json(gromacs_dir)

        config_path = tmp_path / "config.yaml"
        config_path.write_text("name: test\n")

        runner = CliRunner()
        result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "run1" in result.output


class TestStatusCli:
    """End-to-end tests for ``polyzymd status``."""

    def setup_method(self):
        set_color_support(TerminalColorSupport.NONE)

    def _make_dummy_config(self, tmp_path: Path) -> Path:
        config_path = tmp_path / "config.yaml"
        config_path.write_text("name: test\n")
        return config_path

    def test_status_shows_completed_replicate(self, tmp_path):
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
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep1), (2, rep2), (3, rep3)]

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
        scratch = tmp_path / "scratch"
        scratch.mkdir()

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
        assert "70.0/100.0 ns" in result.output

    def test_status_summary_not_all_completed_when_running(self, tmp_path):
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        base = "fnIII_apo_none_100ns_310K"
        rep1 = scratch / f"{base}_run1"
        rep2 = scratch / f"{base}_run2"

        _write_progress_json(rep1)
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
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        base = "fnIII_apo_none_100ns_310K"
        rep1 = scratch / f"{base}_run1"
        rep2 = scratch / f"{base}_run2"
        rep3 = scratch / f"{base}_run3"

        _write_progress_json(rep1)
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
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep1), (2, rep2), (3, rep3)]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "All 3 replicates completed" not in result.output
        assert "need attention" in result.output
        assert "still running" in result.output

    def test_status_stale_running_reclassified_by_filesystem_scan(self, tmp_path):
        scratch = tmp_path / "scratch"
        scratch.mkdir()

        base = "fnIII_apo_none_100ns_310K"
        rep1 = scratch / f"{base}_run1"
        rep2 = scratch / f"{base}_run2"

        _write_progress_json(rep1)
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
        import os
        import time

        scratch = tmp_path / "scratch"
        scratch.mkdir()

        base = "fnIII_apo_none_100ns_310K"
        rep_dir = scratch / f"{base}_run1"
        rep_dir.mkdir(parents=True)

        _write_progress_json(
            rep_dir,
            completed_steps=25_000_000,
            seg_status=SegmentStatus.RUNNING,
            sim_status=SimulationStatus.RUNNING,
            duration_ns=50.0,
        )

        prod_dir = rep_dir / "production_0"
        prod_dir.mkdir()
        chk_file = prod_dir / "production_0_checkpoint.chk"
        chk_file.write_bytes(b"fake checkpoint data")
        old_time = time.time() - 1200
        os.utime(str(chk_file), (old_time, old_time))

        mock_cfg = _mock_sim_config(scratch)
        mock_cfg.discover_replicate_dirs.return_value = [(1, rep_dir)]

        config_path = self._make_dummy_config(tmp_path)
        runner = CliRunner()

        with patch("polyzymd.config.schema.SimulationConfig.from_yaml", return_value=mock_cfg):
            result = runner.invoke(cli, ["status", "-c", str(config_path)])

        assert result.exit_code == 0, f"Output: {result.output}"
        assert "interrupted" in result.output
        assert "need attention" in result.output
        assert "All 1 replicates completed" not in result.output
