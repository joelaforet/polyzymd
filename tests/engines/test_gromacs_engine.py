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
