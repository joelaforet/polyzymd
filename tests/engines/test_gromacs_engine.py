"""Tests for GROMACS engine override resolution."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

from polyzymd.engines.base import EngineSubmitRequest
from polyzymd.engines.gromacs.engine import GromacsEngine
from polyzymd.workflow.slurm import SlurmConfig


def _make_config(
    *,
    ntmpi: int = 1,
    ntomp: int = 8,
    gpu: bool = False,
    gpus: int = 1,
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
    config.gromacs.gpus = gpus
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

    def test_gpu_enabled_uses_config_gpus(self):
        config = _make_config(gpu=True, gpus=2)
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        base = SlurmConfig(gpus=0)
        resolved = engine._resolve_slurm_config(base)
        assert resolved.gpus == 2

    def test_gpu_enabled_single_gpu_default(self):
        config = _make_config(gpu=True)
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        base = SlurmConfig(gpus=0)
        resolved = engine._resolve_slurm_config(base)
        assert resolved.gpus == 1

    def test_gpu_multi_gpu_from_config(self):
        """gpu=True with gpus=4 should request 4 GPUs."""
        config = _make_config(gpu=True, gpus=4)
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        base = SlurmConfig(gpus=0)
        resolved = engine._resolve_slurm_config(base)
        assert resolved.gpus == 4

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


class TestResolveMdrunFlagsGPU:
    """Test GPU offload flag auto-composition in _resolve_mdrun_flags."""

    def test_gpu_flags_auto_added(self):
        """gpu:true should auto-add all four GPU offload flags."""
        config = _make_config(gpu=True, ntmpi=1, ntomp=12, mdrun_flags="")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=12)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-nb gpu" in flags
        assert "-pme gpu" in flags
        assert "-bonded gpu" in flags
        assert "-update gpu" in flags

    def test_gpu_flags_not_duplicated(self):
        """User-provided GPU flags should not be duplicated."""
        config = _make_config(
            gpu=True,
            ntmpi=1,
            ntomp=12,
            mdrun_flags="-nb gpu -pme gpu -bonded gpu -update gpu",
        )
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=12)
        flags = engine._resolve_mdrun_flags(slurm)
        assert flags.count("-nb") == 1
        assert flags.count("-pme") == 1
        assert flags.count("-bonded") == 1
        assert flags.count("-update") == 1

    def test_user_override_nb_cpu_respected(self):
        """User specifying -nb cpu should prevent auto-add of -nb gpu."""
        config = _make_config(gpu=True, ntmpi=1, ntomp=12, mdrun_flags="-nb cpu")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=12)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-nb cpu" in flags
        assert "-nb gpu" not in flags
        # Other GPU flags should still be auto-added
        assert "-pme gpu" in flags
        assert "-bonded gpu" in flags
        assert "-update gpu" in flags

    def test_user_override_update_cpu_respected(self):
        """-update cpu should prevent -update gpu but leave others."""
        config = _make_config(gpu=True, ntmpi=1, ntomp=12, mdrun_flags="-update cpu")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=12)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-update cpu" in flags
        assert "-update gpu" not in flags
        assert "-nb gpu" in flags
        assert "-pme gpu" in flags
        assert "-bonded gpu" in flags

    def test_cpu_mode_no_gpu_flags(self):
        """gpu:false should not add any GPU offload flags."""
        config = _make_config(gpu=False, ntmpi=1, ntomp=8, mdrun_flags="")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=8)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-nb gpu" not in flags
        assert "-pme gpu" not in flags
        assert "-bonded gpu" not in flags
        assert "-update gpu" not in flags
        assert "-ntmpi 1" in flags
        assert "-ntomp 8" in flags

    def test_gpu_flags_with_thread_counts(self):
        """GPU flags should coexist with -ntmpi and -ntomp."""
        config = _make_config(gpu=True, ntmpi=1, ntomp=12, mdrun_flags="")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=12)
        flags = engine._resolve_mdrun_flags(slurm)
        assert "-ntmpi 1" in flags
        assert "-ntomp 12" in flags
        assert "-nb gpu" in flags
        assert "-pme gpu" in flags
        assert "-bonded gpu" in flags
        assert "-update gpu" in flags

    def test_gpu_partial_user_override(self):
        """User specifies some GPU flags; others auto-added."""
        config = _make_config(gpu=True, ntmpi=1, ntomp=12, mdrun_flags="-nb gpu -bonded cpu")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        slurm = SlurmConfig(ntasks=1, cpus_per_task=12)
        flags = engine._resolve_mdrun_flags(slurm)
        # User-specified flags preserved
        assert "-nb gpu" in flags
        assert "-bonded cpu" in flags
        # Non-specified GPU flags auto-added
        assert "-pme gpu" in flags
        assert "-update gpu" in flags
        # No duplicates
        assert flags.count("-nb") == 1
        assert flags.count("-bonded") == 1


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


class TestEngineSubmitModuleLoad:
    """Tests that submit() passes module_load to run_sbatch."""

    @patch("polyzymd.workflow.slurm_submit.run_sbatch")
    def test_submit_passes_module_load(self, mock_run_sbatch):
        """submit() should pass gromacs.module_load to run_sbatch."""
        config = _make_config(module_load="module load slurm/blanca gromacs/2024.2")
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        script_path = Path("/tmp/run_rep1.sh")
        engine.prepare_submission = MagicMock(return_value=script_path)

        mock_run_sbatch.return_value = MagicMock(
            returncode=0,
            stdout="Submitted batch job 123",
            stderr="",
        )

        request = EngineSubmitRequest(
            replicate=1,
            config_path=Path("/tmp/config.yaml"),
            working_dir=Path("/tmp/work"),
            slurm_config=SlurmConfig(),
        )

        result = engine.submit(request)

        mock_run_sbatch.assert_called_once_with(script_path, module_load=config.gromacs.module_load)
        assert result["submitted"] is True

    @patch("polyzymd.workflow.slurm_submit.run_sbatch")
    @patch("polyzymd.engines.gromacs.engine.shutil.which")
    def test_submit_no_module_load(self, mock_which, mock_run_sbatch):
        """submit() should pass module_load=None when not configured."""
        config = _make_config(module_load=None)
        engine = GromacsEngine(config=config, gmx_binary="gmx")
        script_path = Path("/tmp/run_rep2.sh")
        engine.prepare_submission = MagicMock(return_value=script_path)
        mock_which.return_value = "/usr/bin/sbatch"
        mock_run_sbatch.return_value = MagicMock(
            returncode=0,
            stdout="Submitted batch job 456",
            stderr="",
        )

        request = EngineSubmitRequest(
            replicate=2,
            config_path=Path("/tmp/config.yaml"),
            working_dir=Path("/tmp/work"),
            slurm_config=SlurmConfig(),
        )

        result = engine.submit(request)

        mock_run_sbatch.assert_called_once_with(script_path, module_load=None)
        assert result["submitted"] is True


class TestPrepareSubmissionPassThrough:
    """Tests for passing config fields to script generator."""

    @patch("polyzymd.engines.gromacs.engine.GromacsSlurmScriptGenerator")
    def test_prepare_submission_passes_env_exports_and_setup_commands(
        self, mock_generator_cls, tmp_path
    ):
        """prepare_submission should pass env_exports and setup_commands to generator."""
        config = _make_config()
        config.gromacs.env_exports = {
            "GMX_GPU_DD_COMMS": "true",
            "GMX_FORCE_UPDATE_DEFAULT_GPU": "true",
        }
        config.gromacs.setup_commands = [
            "source /projects/software/gromacs/bin/GMXRC",
            "export PATH=$PATH:/projects/software/plumed/bin",
        ]
        config.gromacs.mdrun_flags_equilibration = None
        config.gromacs.mdrun_flags_production = None
        config.gromacs.command_prefix = None
        config.gromacs.mpi_launcher_flags = ""

        engine = GromacsEngine(config=config, gmx_binary="gmx")

        working_dir = tmp_path / "gromacs"
        working_dir.mkdir(parents=True)
        (working_dir / "CALB_PEG.top").write_text("[ system ]\n")
        (working_dir / "CALB_PEG.gro").write_text("test\n")
        (working_dir / "em.mdp").write_text("integrator = steep\n")
        (working_dir / "prod.mdp").write_text("integrator = md\n")

        mock_generator = MagicMock()
        mock_generator.generate_job_script.return_value = "#!/bin/bash\n"
        mock_generator.save_script.return_value = (
            working_dir / "daisy_chain_scripts" / "run_rep1.sh"
        )
        mock_generator_cls.return_value = mock_generator

        request = EngineSubmitRequest(
            replicate=1,
            config_path=tmp_path / "config.yaml",
            working_dir=working_dir,
            slurm_config=SlurmConfig(),
            extra={"skip_build": True},
        )

        engine.prepare_submission(request)

        kwargs = mock_generator_cls.call_args.kwargs
        assert kwargs["env_exports"] == config.gromacs.env_exports
        assert kwargs["setup_commands"] == config.gromacs.setup_commands

    @patch("polyzymd.engines.gromacs.engine.GromacsSlurmScriptGenerator")
    def test_prepare_submission_passes_stage_specific_mdrun_flags(
        self, mock_generator_cls, tmp_path
    ):
        """prepare_submission should pass resolved stage-specific mdrun flags."""
        config = _make_config(mdrun_flags="-ntomp 8")
        config.gromacs.mdrun_flags_equilibration = "-ntomp 4"
        config.gromacs.mdrun_flags_production = "-ntomp 8 -plumed plumed_setup.dat"
        config.gromacs.command_prefix = None
        config.gromacs.mpi_launcher_flags = ""
        config.gromacs.env_exports = {}
        config.gromacs.setup_commands = []

        engine = GromacsEngine(config=config, gmx_binary="gmx")

        working_dir = tmp_path / "gromacs"
        working_dir.mkdir(parents=True)
        (working_dir / "CALB_PEG.top").write_text("[ system ]\n")
        (working_dir / "CALB_PEG.gro").write_text("test\n")
        (working_dir / "em.mdp").write_text("integrator = steep\n")
        (working_dir / "prod.mdp").write_text("integrator = md\n")

        mock_generator = MagicMock()
        mock_generator.generate_job_script.return_value = "#!/bin/bash\n"
        mock_generator.save_script.return_value = (
            working_dir / "daisy_chain_scripts" / "run_rep1.sh"
        )
        mock_generator_cls.return_value = mock_generator

        request = EngineSubmitRequest(
            replicate=1,
            config_path=tmp_path / "config.yaml",
            working_dir=working_dir,
            slurm_config=SlurmConfig(ntasks=1, cpus_per_task=8),
            extra={"skip_build": True},
        )

        engine.prepare_submission(request)

        kwargs = mock_generator_cls.call_args.kwargs
        assert kwargs["mdrun_flags_eq"] == "-ntomp 4 -ntmpi 1"
        assert kwargs["mdrun_flags_prod"] == "-ntomp 8 -plumed plumed_setup.dat -ntmpi 1"

    @patch("polyzymd.engines.gromacs.engine.GromacsSlurmScriptGenerator")
    def test_prepare_submission_passes_wrapper_and_mpi_launcher(self, mock_generator_cls, tmp_path):
        """prepare_submission should pass command wrapper and MPI launcher flags."""
        config = _make_config()
        config.gromacs.command_prefix = (
            "singularity exec --rocm --bind /scratch /path/to/gromacs.sif"
        )
        config.gromacs.mpi_launcher_flags = "-genv I_MPI_FABRICS shm:tcp"
        config.gromacs.mdrun_flags_equilibration = None
        config.gromacs.mdrun_flags_production = None
        config.gromacs.env_exports = {}
        config.gromacs.setup_commands = []

        engine = GromacsEngine(config=config, gmx_binary="gmx_mpi")

        working_dir = tmp_path / "gromacs"
        working_dir.mkdir(parents=True)
        (working_dir / "CALB_PEG.top").write_text("[ system ]\n")
        (working_dir / "CALB_PEG.gro").write_text("test\n")
        (working_dir / "em.mdp").write_text("integrator = steep\n")
        (working_dir / "prod.mdp").write_text("integrator = md\n")

        mock_generator = MagicMock()
        mock_generator.generate_job_script.return_value = "#!/bin/bash\n"
        mock_generator.save_script.return_value = (
            working_dir / "daisy_chain_scripts" / "run_rep1.sh"
        )
        mock_generator_cls.return_value = mock_generator

        request = EngineSubmitRequest(
            replicate=1,
            config_path=tmp_path / "config.yaml",
            working_dir=working_dir,
            slurm_config=SlurmConfig(),
            extra={"skip_build": True},
        )

        engine.prepare_submission(request)

        kwargs = mock_generator_cls.call_args.kwargs
        assert kwargs["command_prefix"] == config.gromacs.command_prefix
        assert kwargs["mpi_launcher_flags"] == config.gromacs.mpi_launcher_flags


class TestFromConfigBinaryResolution:
    """Tests for deferred GROMACS binary resolution."""

    @patch("polyzymd.engines.gromacs.engine.resolve_gromacs_binary")
    def test_from_config_deferred_uses_configured_binary(self, mock_resolve):
        """defer_binary=True should skip PATH probing and use configured binary."""
        config = _make_config(gmx_binary="gmx_mpi")
        engine = GromacsEngine.from_config(config, defer_binary=True)

        assert engine._gmx_binary == "gmx_mpi"
        mock_resolve.assert_not_called()

    @patch("polyzymd.engines.gromacs.engine.resolve_gromacs_binary")
    def test_from_config_deferred_uses_default_when_missing(self, mock_resolve):
        """defer_binary=True should fallback to gmx when unset in config."""
        config = _make_config(gmx_binary=None)
        engine = GromacsEngine.from_config(config, defer_binary=True)

        assert engine._gmx_binary == "gmx"
        mock_resolve.assert_not_called()
