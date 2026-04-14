"""Tests for GROMACS engine override resolution."""

from __future__ import annotations

from unittest.mock import MagicMock

from polyzymd.engines.gromacs.engine import GromacsEngine
from polyzymd.workflow.slurm import SlurmConfig


def _make_config(
    *,
    ntmpi: int = 1,
    ntomp: int = 8,
    gpu: bool = False,
    memory: str = "16G",
    mdrun_flags: str = "",
    grompp_flags: str = "-maxwarn 1",
    gmx_binary: str | None = None,
    module_load: str | None = None,
) -> MagicMock:
    """Build a minimal mock SimulationConfig with GROMACS settings."""
    config = MagicMock()
    config.gromacs.ntmpi = ntmpi
    config.gromacs.ntomp = ntomp
    config.gromacs.gpu = gpu
    config.gromacs.memory = memory
    config.gromacs.mdrun_flags = mdrun_flags
    config.gromacs.grompp_flags = grompp_flags
    config.gromacs.gmx_binary = gmx_binary
    config.gromacs.module_load = module_load
    config.enzyme.name = "CALB"
    config.polymers.enabled = True
    config.polymers.type_prefix = "PEG"
    return config


class TestResolveSlurmConfig:
    """Test GROMACS hardware overrides on SlurmConfig."""

    def test_cpu_only_forces_zero_gpus(self):
        config = _make_config(gpu=False)
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        base = SlurmConfig(gpus=1)
        resolved = engine._resolve_slurm_config(base)
        assert resolved.gpus == 0

    def test_gpu_enabled_preserves_gpus(self):
        config = _make_config(gpu=True)
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        base = SlurmConfig(gpus=2)
        resolved = engine._resolve_slurm_config(base)
        assert resolved.gpus == 2

    def test_gpu_enabled_sets_minimum_one(self):
        config = _make_config(gpu=True)
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        base = SlurmConfig(gpus=0)
        resolved = engine._resolve_slurm_config(base)
        assert resolved.gpus == 1

    def test_ntasks_from_ntmpi(self):
        config = _make_config(ntmpi=4)
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        base = SlurmConfig(ntasks=1)
        resolved = engine._resolve_slurm_config(base)
        assert resolved.ntasks == 4

    def test_cpus_per_task_from_ntomp(self):
        config = _make_config(ntomp=16)
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        base = SlurmConfig()
        resolved = engine._resolve_slurm_config(base)
        assert resolved.cpus_per_task == 16

    def test_memory_override(self):
        config = _make_config(memory="32G")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        base = SlurmConfig(memory="3G")
        resolved = engine._resolve_slurm_config(base)
        assert resolved.memory == "32G"


class TestResolveMdrunFlags:
    """Test auto-composition of -ntmpi/-ntomp flags."""

    def test_auto_adds_ntmpi_ntomp(self):
        config = _make_config(ntmpi=2, ntomp=16, mdrun_flags="")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=2, cpus_per_task=16)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-ntmpi 2" in flags
        assert "-ntomp 16" in flags

    def test_preserves_explicit_ntmpi(self):
        config = _make_config(mdrun_flags="-ntmpi 4")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=8)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-ntmpi 4" in flags
        assert flags.count("-ntmpi") == 1
        assert "-ntomp 8" in flags

    def test_preserves_explicit_ntomp(self):
        config = _make_config(mdrun_flags="-ntomp 32")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=8)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-ntomp 32" in flags
        assert flags.count("-ntomp") == 1
        assert "-ntmpi 1" in flags

    def test_preserves_extra_flags(self):
        config = _make_config(mdrun_flags="-nb gpu -pme gpu")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=8)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-nb gpu" in flags
        assert "-pme gpu" in flags
        assert "-ntmpi 1" in flags
        assert "-ntomp 8" in flags

    def test_empty_flags_gets_defaults(self):
        config = _make_config(mdrun_flags="", ntmpi=1, ntomp=8)
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=8)
        flags = engine._resolve_mdrun_flags(slurm)
        assert flags == "-ntmpi 1 -ntomp 8"


class TestResolveMdrunFlagsMPIBinary:
    """Test -ntmpi suppression for real-MPI GROMACS builds (gmx_mpi)."""

    def test_mpi_binary_omits_ntmpi(self):
        """gmx_mpi does not support -ntmpi; only -ntomp should be added."""
        config = _make_config(ntmpi=1, ntomp=16, mdrun_flags="")
        engine = GromacsEngine(config=config, gmx_binary="gmx_mpi")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=16)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-ntmpi" not in flags
        assert "-ntomp 16" in flags

    def test_mpi_binary_full_path_omits_ntmpi(self):
        """Even a full path like /usr/bin/gmx_mpi should suppress -ntmpi."""
        config = _make_config(ntmpi=2, ntomp=8, mdrun_flags="")
        engine = GromacsEngine(config=config, gmx_binary="/usr/local/bin/gmx_mpi")
        slurm = SlurmConfig(ntasks=2, cpus_per_task=8)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-ntmpi" not in flags
        assert "-ntomp 8" in flags

    def test_mpi_double_binary_omits_ntmpi(self):
        """gmx_mpi_d (double precision MPI) should also suppress -ntmpi."""
        config = _make_config(ntmpi=1, ntomp=4, mdrun_flags="")
        engine = GromacsEngine(config=config, gmx_binary="gmx_mpi_d")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=4)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-ntmpi" not in flags
        assert "-ntomp 4" in flags

    def test_thread_mpi_binary_includes_ntmpi(self):
        """Plain gmx (thread-MPI) should still get -ntmpi."""
        config = _make_config(ntmpi=2, ntomp=8, mdrun_flags="")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=2, cpus_per_task=8)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-ntmpi 2" in flags
        assert "-ntomp 8" in flags

    def test_mpi_binary_preserves_explicit_ntomp(self):
        """User-specified -ntomp should be preserved with gmx_mpi."""
        config = _make_config(mdrun_flags="-ntomp 32")
        engine = GromacsEngine(config=config, gmx_binary="gmx_mpi")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=8)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-ntomp 32" in flags
        assert flags.count("-ntomp") == 1
        assert "-ntmpi" not in flags

    def test_mpi_binary_with_extra_flags(self):
        """Extra flags like -nb gpu should pass through with gmx_mpi."""
        config = _make_config(mdrun_flags="-nb gpu -pme gpu")
        engine = GromacsEngine(config=config, gmx_binary="gmx_mpi")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=16)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-nb gpu" in flags
        assert "-pme gpu" in flags
        assert "-ntomp 16" in flags
        assert "-ntmpi" not in flags


class TestDiscoverLegacyReplicates:
    """Tests for GromacsEngine.discover_legacy_replicates()."""

    def test_finds_legacy_gromacs_dirs(self, tmp_path):
        """Should find replicate_N/gromacs/ directories in projects_directory."""
        projects = tmp_path / "projects"
        for i in (1, 2, 3):
            (projects / f"replicate_{i}" / "gromacs").mkdir(parents=True)

        config = _make_config()
        config.output.projects_directory = projects
        engine = GromacsEngine(config=config, gmx_binary="gmx")

        result = engine.discover_legacy_replicates(config)
        assert len(result) == 3
        assert [r[0] for r in result] == [1, 2, 3]
        for num, path in result:
            assert path == projects / f"replicate_{num}" / "gromacs"

    def test_returns_sorted(self, tmp_path):
        """Should return replicates in numeric order."""
        projects = tmp_path / "projects"
        for i in (5, 2, 8):
            (projects / f"replicate_{i}" / "gromacs").mkdir(parents=True)

        config = _make_config()
        config.output.projects_directory = projects
        engine = GromacsEngine(config=config, gmx_binary="gmx")

        result = engine.discover_legacy_replicates(config)
        assert [r[0] for r in result] == [2, 5, 8]

    def test_ignores_replicate_without_gromacs_subdir(self, tmp_path):
        """Should skip replicate dirs that lack a gromacs/ subdirectory."""
        projects = tmp_path / "projects"
        (projects / "replicate_1" / "gromacs").mkdir(parents=True)
        (projects / "replicate_2").mkdir(parents=True)

        config = _make_config()
        config.output.projects_directory = projects
        engine = GromacsEngine(config=config, gmx_binary="gmx")

        result = engine.discover_legacy_replicates(config)
        assert len(result) == 1
        assert result[0][0] == 1

    def test_ignores_non_matching_dirs(self, tmp_path):
        """Should ignore directories that don't match replicate_N pattern."""
        projects = tmp_path / "projects"
        (projects / "replicate_1" / "gromacs").mkdir(parents=True)
        (projects / "other_stuff").mkdir(parents=True)
        (projects / "replicate_bad" / "gromacs").mkdir(parents=True)

        config = _make_config()
        config.output.projects_directory = projects
        engine = GromacsEngine(config=config, gmx_binary="gmx")

        result = engine.discover_legacy_replicates(config)
        assert len(result) == 1
        assert result[0][0] == 1

    def test_empty_when_projects_dir_missing(self, tmp_path):
        """Should return empty list when projects_directory doesn't exist."""
        config = _make_config()
        config.output.projects_directory = tmp_path / "nonexistent"
        engine = GromacsEngine(config=config, gmx_binary="gmx")

        result = engine.discover_legacy_replicates(config)
        assert result == []

    def test_empty_when_no_replicates(self, tmp_path):
        """Should return empty when projects_directory exists but has no replicate dirs."""
        projects = tmp_path / "projects"
        projects.mkdir()

        config = _make_config()
        config.output.projects_directory = projects
        engine = GromacsEngine(config=config, gmx_binary="gmx")

        result = engine.discover_legacy_replicates(config)
        assert result == []
