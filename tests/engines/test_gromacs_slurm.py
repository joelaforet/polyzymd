"""Tests for GROMACS SLURM script generation."""

import re

from polyzymd.engines.gromacs.slurm import GromacsSlurmScriptGenerator
from polyzymd.workflow.slurm import SlurmConfig


def _generator() -> GromacsSlurmScriptGenerator:
    return GromacsSlurmScriptGenerator(
        slurm_config=SlurmConfig(
            partition="aa100",
            qos="normal",
            account="ucb625_asc1",
            time_limit="23:59:59",
            email="",
            nodes=1,
            ntasks=1,
            memory="3G",
            gpus=1,
        ),
        pixi_env="cuda-12-4",
        gmx_binary="gmx",
        grompp_flags="-maxwarn 1",
        mdrun_flags="-nb gpu",
        module_load="module load gromacs/2024",
    )


def test_script_contains_expected_slurm_and_restart_logic(monkeypatch) -> None:
    """Generated script should include directives and restart features."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    script = _generator().generate_job_script(
        config_path="/path/config.yaml",
        replicate=2,
        working_dir="/scratch/run2/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp", "eq_02_npt.mdp"],
        job_name="gmx_r2",
        output_file="slurm_logs/gmx_r2.%j.out",
    )

    assert "#SBATCH --partition=aa100" in script
    assert "#SBATCH --job-name=gmx_r2" in script
    assert "#SBATCH --time=23:59:59" in script
    assert "-maxh $MAXH" in script
    assert "-cpi state.cpt" in script
    assert 'sbatch "$THIS_SCRIPT"' in script
    assert "if [ ! -f eq_01.gro ]; then" in script
    assert "if [ ! -f eq_02.gro ]; then" in script


def test_parse_wall_time_hours_hms() -> None:
    """Wall time parser should parse HH:MM:SS format."""
    hours = _generator()._parse_wall_time_hours("23:30:00")
    assert hours == 23.5


def test_parse_wall_time_hours_with_days() -> None:
    """Wall time parser should parse D-HH:MM:SS format."""
    hours = _generator()._parse_wall_time_hours("1-12:00:00")
    assert hours == 36.0


def test_em_block_single_stage() -> None:
    """Single-stage EM block should include health check and no soft EM artifacts."""
    block = _generator()._generate_energy_minimization_block()

    assert "em_soft" not in block
    assert "em.mdp" in block
    assert "em.gro" in block
    assert "force.*not finite" in block

    em_if_idx = block.index("if [ ! -f em.gro ]; then")
    health_idx = block.index("force.*not finite")
    em_else_idx = block.index('else\n    echo "Skipping energy minimization (em.gro exists)."')
    assert em_if_idx < health_idx < em_else_idx


def test_script_never_contains_soft_em(monkeypatch) -> None:
    """Generated scripts should never reference soft EM artifacts."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    script = _generator().generate_job_script(
        config_path="/path/config.yaml",
        replicate=1,
        working_dir="/scratch/run1/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp"],
    )

    assert "em_soft" not in script


def test_cpu_only_omits_gpu_directive(monkeypatch) -> None:
    """CPU-only scripts should omit GPU SBATCH directives."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    generator = GromacsSlurmScriptGenerator(
        slurm_config=SlurmConfig(
            partition="aa100",
            qos="normal",
            account="ucb625_asc1",
            time_limit="23:59:59",
            email="",
            nodes=1,
            ntasks=1,
            cpus_per_task=1,
            memory="3G",
            gpus=0,
        ),
        pixi_env="cuda-12-4",
        gmx_binary="gmx",
    )

    script = generator.generate_job_script(
        config_path="/path/config.yaml",
        replicate=1,
        working_dir="/scratch/run1/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp"],
    )

    assert "--gres" not in script
    assert "--gpus" not in script


def test_cpus_per_task_in_script(monkeypatch) -> None:
    """Script should include cpus-per-task when configured."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    generator = GromacsSlurmScriptGenerator(
        slurm_config=SlurmConfig(
            partition="aa100",
            qos="normal",
            account="ucb625_asc1",
            time_limit="23:59:59",
            email="",
            nodes=1,
            ntasks=1,
            cpus_per_task=16,
            memory="3G",
            gpus=0,
        ),
        pixi_env="cuda-12-4",
        gmx_binary="gmx",
    )

    script = generator.generate_job_script(
        config_path="/path/config.yaml",
        replicate=1,
        working_dir="/scratch/run1/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp"],
    )

    assert "#SBATCH --cpus-per-task=16" in script


def test_mdrun_flags_in_all_stages(monkeypatch) -> None:
    """EM should use filtered flags while equilibration and production use full flags."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    script = _generator().generate_job_script(
        config_path="/path/config.yaml",
        replicate=2,
        working_dir="/scratch/run2/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp", "eq_02_npt.mdp"],
        job_name="gmx_r2",
        output_file="slurm_logs/gmx_r2.%j.out",
    )

    assert "$MDRUN -deffnm em -cpo em.cpt $MDRUN_FLAGS_EM -v" in script
    assert "$MDRUN -deffnm eq_01 -cpo eq_01.cpt $MDRUN_FLAGS_EQ -v" in script
    assert "$MDRUN -deffnm eq_02 -cpo eq_02.cpt $MDRUN_FLAGS_EQ -v" in script
    assert "$MDRUN -deffnm prod -cpo state.cpt -maxh $MAXH $MDRUN_FLAGS_PROD -v" in script
    assert script.count("$MDRUN_FLAGS") >= 1


def test_mdrun_flags_variable_defined_before_em(monkeypatch) -> None:
    """MDRUN flag variables should be declared before EM block commands."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    script = _generator().generate_job_script(
        config_path="/path/config.yaml",
        replicate=1,
        working_dir="/scratch/run1/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp"],
    )

    assert script.index("MDRUN_FLAGS=") < script.index("em.mdp")
    assert script.index("MDRUN_FLAGS_EQ=") < script.index("em.mdp")
    assert script.index("MDRUN_FLAGS_PROD=") < script.index("em.mdp")
    assert script.index("MDRUN_FLAGS_EM=") < script.index("em.mdp")


def test_stage_specific_mdrun_flags_override(monkeypatch) -> None:
    """Equilibration and production should use stage-specific flags when set."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )
    generator = GromacsSlurmScriptGenerator(
        slurm_config=SlurmConfig(
            partition="aa100",
            qos="normal",
            account="ucb625_asc1",
            time_limit="23:59:59",
            email="",
            nodes=1,
            ntasks=1,
            memory="3G",
            gpus=1,
        ),
        pixi_env="cuda-12-4",
        gmx_binary="gmx",
        mdrun_flags="-ntomp 8",
        mdrun_flags_eq="-ntomp 4",
        mdrun_flags_prod="-ntomp 8 -plumed plumed_setup.dat",
    )
    script = generator.generate_job_script(
        config_path="/path/config.yaml",
        replicate=1,
        working_dir="/scratch/run1/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp"],
    )
    assert 'MDRUN_FLAGS_EQ="-ntomp 4"' in script
    assert 'MDRUN_FLAGS_PROD="-ntomp 8 -plumed plumed_setup.dat"' in script
    assert "$MDRUN -deffnm eq_01 -cpo eq_01.cpt $MDRUN_FLAGS_EQ -v" in script
    assert "$MDRUN -deffnm prod -cpo state.cpt -maxh $MAXH $MDRUN_FLAGS_PROD -v" in script


def test_mdrun_variable_uses_mpirun_for_mpi_binary(monkeypatch) -> None:
    """MPI binaries should set MDRUN to use mpirun."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    generator = GromacsSlurmScriptGenerator(
        slurm_config=SlurmConfig(
            partition="aa100",
            qos="normal",
            account="ucb625_asc1",
            time_limit="23:59:59",
            email="",
            nodes=1,
            ntasks=1,
            memory="3G",
            gpus=1,
        ),
        pixi_env="cuda-12-4",
        gmx_binary="gmx_mpi",
        grompp_flags="-maxwarn 1",
        mdrun_flags="-nb gpu",
        module_load="module load gromacs/2024",
    )

    script = generator.generate_job_script(
        config_path="/path/config.yaml",
        replicate=1,
        working_dir="/scratch/run1/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp"],
    )

    assert 'MDRUN="mpirun $GMX mdrun"' in script
    assert "$MDRUN -deffnm em -cpo em.cpt $MDRUN_FLAGS_EM -v" in script
    assert "$MDRUN -deffnm eq_01 -cpo eq_01.cpt $MDRUN_FLAGS_EQ -v" in script
    assert "$MDRUN -deffnm prod -cpo state.cpt -maxh $MAXH $MDRUN_FLAGS_PROD -v" in script


def test_mpirun_uses_launcher_flags_when_configured(monkeypatch) -> None:
    """MPI launcher flags should be inserted for real-MPI binaries."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    generator = GromacsSlurmScriptGenerator(
        slurm_config=SlurmConfig(
            partition="aa100",
            qos="normal",
            account="ucb625_asc1",
            time_limit="23:59:59",
            email="",
            nodes=1,
            ntasks=1,
            memory="3G",
            gpus=1,
        ),
        pixi_env="cuda-12-4",
        gmx_binary="gmx_mpi",
        mpi_launcher_flags="-genv I_MPI_FABRICS shm:tcp",
    )
    script = generator.generate_job_script(
        config_path="/path/config.yaml",
        replicate=1,
        working_dir="/scratch/run1/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp"],
    )

    assert 'MDRUN="mpirun -genv I_MPI_FABRICS shm:tcp $GMX mdrun"' in script


def test_mdrun_variable_omits_mpirun_for_non_mpi_binary(monkeypatch) -> None:
    """Non-MPI binaries should use direct GMX mdrun."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    script = _generator().generate_job_script(
        config_path="/path/config.yaml",
        replicate=1,
        working_dir="/scratch/run1/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp"],
    )

    assert 'MDRUN="$GMX mdrun"' in script
    assert 'MDRUN="mpirun $GMX mdrun"' not in script


def test_command_prefix_wraps_gmx_commands(monkeypatch) -> None:
    """command_prefix should prepend wrapper to GMX command variable."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    generator = GromacsSlurmScriptGenerator(
        slurm_config=SlurmConfig(
            partition="aa100",
            qos="normal",
            account="ucb625_asc1",
            time_limit="23:59:59",
            email="",
            nodes=1,
            ntasks=1,
            memory="3G",
            gpus=1,
        ),
        pixi_env="cuda-12-4",
        gmx_binary="gmx",
        command_prefix="singularity exec --rocm --bind /scratch /path/to/gromacs.sif",
    )
    script = generator.generate_job_script(
        config_path="/path/config.yaml",
        replicate=1,
        working_dir="/scratch/run1/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp"],
    )

    assert 'GMX="singularity exec --rocm --bind /scratch /path/to/gromacs.sif gmx"' in script


def test_nodelist_and_typed_gres_render(monkeypatch) -> None:
    """typed --gres and nodelist should be rendered for gres style."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    generator = GromacsSlurmScriptGenerator(
        slurm_config=SlurmConfig(
            partition="aa100",
            qos="normal",
            account="ucb625_asc1",
            time_limit="23:59:59",
            email="",
            nodes=1,
            ntasks=1,
            memory="3G",
            gpus=2,
            gpu_type="a100",
            gpu_directive_style="gres",
            nodelist="nodeA40-01",
        ),
        pixi_env="cuda-12-4",
        gmx_binary="gmx",
    )
    script = generator.generate_job_script(
        config_path="/path/config.yaml",
        replicate=1,
        working_dir="/scratch/run1/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp"],
    )

    assert "#SBATCH --gres=gpu:a100:2" in script
    assert "#SBATCH --nodelist=nodeA40-01" in script


def test_grompp_calls_never_use_mpirun(monkeypatch) -> None:
    """grompp commands should never be prefixed with mpirun."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    generator = GromacsSlurmScriptGenerator(
        slurm_config=SlurmConfig(
            partition="aa100",
            qos="normal",
            account="ucb625_asc1",
            time_limit="23:59:59",
            email="",
            nodes=1,
            ntasks=1,
            memory="3G",
            gpus=1,
        ),
        pixi_env="cuda-12-4",
        gmx_binary="gmx_mpi",
        grompp_flags="-maxwarn 1",
        mdrun_flags="-nb gpu",
        module_load="module load gromacs/2024",
    )

    script = generator.generate_job_script(
        config_path="/path/config.yaml",
        replicate=1,
        working_dir="/scratch/run1/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp"],
    )

    assert "mpirun $GMX grompp" not in script


def test_env_exports_and_setup_commands_rendered_in_order(monkeypatch) -> None:
    """env_exports and setup_commands should render after module load."""
    monkeypatch.setattr(
        "polyzymd.engines.gromacs.slurm._discover_manifest_path",
        lambda: "/tmp/pixi.toml",
    )

    generator = GromacsSlurmScriptGenerator(
        slurm_config=SlurmConfig(
            partition="aa100",
            qos="normal",
            account="ucb625_asc1",
            time_limit="23:59:59",
            email="",
            nodes=1,
            ntasks=1,
            memory="3G",
            gpus=1,
        ),
        pixi_env="cuda-12-4",
        gmx_binary="gmx",
        module_load="module load gromacs/2024",
        env_exports={
            "GMX_GPU_DD_COMMS": "true",
            "GMX_GPU_PME_PP_COMMS": "true",
        },
        setup_commands=[
            "source /projects/software/gromacs/bin/GMXRC",
            "export PATH=$PATH:/projects/software/plumed/bin",
        ],
    )

    script = generator.generate_job_script(
        config_path="/path/config.yaml",
        replicate=1,
        working_dir="/scratch/run1/gromacs",
        system_prefix="enzyme_polymer",
        equilibration_mdps=["eq_01_nvt.mdp"],
    )

    module_idx = script.index("module load gromacs/2024")
    env_idx = script.index('export GMX_GPU_DD_COMMS="true"')
    setup_idx = script.index("source /projects/software/gromacs/bin/GMXRC")
    resolve_idx = script.index("THIS_SCRIPT=")
    assert module_idx < env_idx < setup_idx < resolve_idx


class TestSignalInfrastructure:
    """Signal handling and preemption resilience tests."""

    def test_constraint_directive_rendered_when_configured(self, monkeypatch) -> None:
        """Constraint directive should appear when configured."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        generator = GromacsSlurmScriptGenerator(
            slurm_config=SlurmConfig(
                partition="blanca,blanca-shirts",
                qos="preemptable",
                account="blanca-shirts",
                time_limit="23:59:59",
                nodes=1,
                ntasks=1,
                cpus_per_task=12,
                memory="64G",
                gpus=1,
                constraint="A40|A100",
            ),
            pixi_env="cuda-12-4",
            gmx_binary="gmx",
            mdrun_flags="-nb gpu",
        )
        script = generator.generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert "#SBATCH --constraint=A40|A100" in script

    def test_constraint_directive_omitted_when_unset(self, monkeypatch) -> None:
        """No constraint line should appear when constraint is None."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert "#SBATCH --constraint=" not in script

    def test_script_contains_term_trap(self, monkeypatch) -> None:
        """Script should trap SIGTERM and define handle_term function."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert "trap 'handle_term' TERM" in script
        assert "handle_term()" in script

    def test_script_does_not_include_usr1_signal_directive(self, monkeypatch) -> None:
        """GROMACS scripts should not use USR1 signal (that's OpenMM-only)."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert "--signal=B:USR1" not in script
        assert "trap" in script
        assert "USR1" not in script

    def test_script_defines_resubmit_once_function(self, monkeypatch) -> None:
        """Script should define resubmit_once with double-submit guard."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert "resubmit_once()" in script
        assert "RESUBMITTED" in script

    def test_script_defines_run_mdrun_stage_function(self, monkeypatch) -> None:
        """run_mdrun_stage helper should be defined exactly once."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert script.count("run_mdrun_stage()") == 1
        assert script.count("run_mdrun_stage ") >= 3

    def test_mdrun_stages_are_backgrounded(self, monkeypatch) -> None:
        """run_mdrun_stage should background the child and capture PID."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert '"$@" &' in script
        assert "CHILD_PID=$!" in script

    def test_term_received_state_declared_and_checked(self, monkeypatch) -> None:
        """TERM_RECEIVED should be declared and checked in run_mdrun_stage."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert "TERM_RECEIVED=0" in script
        assert '"$TERM_RECEIVED" -eq 1' in script


class TestGPUSlurmScripts:
    """GPU-specific SLURM script regression tests."""

    def test_single_gpu_gres_directive(self, monkeypatch) -> None:
        """Single GPU should emit --gres=gpu:1."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )

        generator = GromacsSlurmScriptGenerator(
            slurm_config=SlurmConfig(
                partition="blanca-shirts",
                qos="blanca-shirts",
                account="blanca-shirts",
                time_limit="23:59:59",
                nodes=1,
                ntasks=1,
                cpus_per_task=12,
                memory="64G",
                gpus=1,
            ),
            pixi_env="cuda-12-4",
            gmx_binary="gmx",
            grompp_flags="-maxwarn 1",
            mdrun_flags="-ntmpi 1 -ntomp 12 -nb gpu -pme gpu -bonded gpu -update gpu",
            module_load="module load gromacs/2024.2",
        )

        script = generator.generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/rc_scratch/run1/gromacs",
            system_prefix="CALB_PEG",
            equilibration_mdps=["eq_01_nvt.mdp", "eq_02_npt.mdp"],
        )

        assert "#SBATCH --gres=gpu:1" in script

    def test_multi_gpu_gres_directive(self, monkeypatch) -> None:
        """Multiple GPUs should emit --gres=gpu:4."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )

        generator = GromacsSlurmScriptGenerator(
            slurm_config=SlurmConfig(
                partition="blanca-shirts",
                qos="blanca-shirts",
                account="blanca-shirts",
                time_limit="23:59:59",
                nodes=1,
                ntasks=4,
                cpus_per_task=12,
                memory="128G",
                gpus=4,
            ),
            pixi_env="cuda-12-4",
            gmx_binary="gmx",
            mdrun_flags="-ntmpi 4 -ntomp 12 -nb gpu -pme gpu -bonded gpu -update gpu",
        )

        script = generator.generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/rc_scratch/run1/gromacs",
            system_prefix="CALB_PEG",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )

        assert "#SBATCH --gres=gpu:4" in script

    def test_gpu_script_gmx_uses_direct_launch(self, monkeypatch) -> None:
        """GPU scripts with gmx (thread-MPI) should NOT use mpirun."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )

        generator = GromacsSlurmScriptGenerator(
            slurm_config=SlurmConfig(
                partition="blanca-shirts",
                qos="blanca-shirts",
                account="blanca-shirts",
                time_limit="23:59:59",
                nodes=1,
                ntasks=1,
                cpus_per_task=12,
                memory="64G",
                gpus=1,
            ),
            pixi_env="cuda-12-4",
            gmx_binary="gmx",
            mdrun_flags="-ntmpi 1 -ntomp 12 -nb gpu -pme gpu -bonded gpu -update gpu",
            module_load="module load gromacs/2024.2",
        )

        script = generator.generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/rc_scratch/run1/gromacs",
            system_prefix="CALB_PEG",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )

        assert 'MDRUN="$GMX mdrun"' in script
        assert "mpirun" not in script

    def test_gpu_script_gmx_mpi_uses_mpirun(self, monkeypatch) -> None:
        """GPU scripts with gmx_mpi should use mpirun launch."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )

        generator = GromacsSlurmScriptGenerator(
            slurm_config=SlurmConfig(
                partition="blanca-shirts",
                qos="blanca-shirts",
                account="blanca-shirts",
                time_limit="23:59:59",
                nodes=1,
                ntasks=4,
                cpus_per_task=12,
                memory="128G",
                gpus=4,
            ),
            pixi_env="cuda-12-4",
            gmx_binary="gmx_mpi",
            mdrun_flags="-ntomp 12 -nb gpu -pme gpu",
            module_load="module load gromacs/2024.2",
        )

        script = generator.generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/rc_scratch/run1/gromacs",
            system_prefix="CALB_PEG",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )

        assert 'MDRUN="mpirun $GMX mdrun"' in script

    def test_gpu_type_directive_style(self, monkeypatch) -> None:
        """gpu_directive_style='gpus' with gpu_type should emit --gpus=a100:2."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )

        generator = GromacsSlurmScriptGenerator(
            slurm_config=SlurmConfig(
                partition="gpu",
                qos="normal",
                account="myaccount",
                time_limit="23:59:59",
                nodes=1,
                ntasks=1,
                cpus_per_task=8,
                memory="64G",
                gpus=2,
                gpu_type="a100",
                gpu_directive_style="gpus",
            ),
            pixi_env="cuda-12-4",
            gmx_binary="gmx",
            mdrun_flags="-nb gpu",
        )

        script = generator.generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="CALB_PEG",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )

        assert "#SBATCH --gpus=a100:2" in script
        assert "--gres" not in script

    def test_gpu_script_contains_module_load(self, monkeypatch) -> None:
        """GPU scripts should include the module load command."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )

        generator = GromacsSlurmScriptGenerator(
            slurm_config=SlurmConfig(
                partition="blanca-shirts",
                qos="blanca-shirts",
                account="blanca-shirts",
                time_limit="23:59:59",
                nodes=1,
                ntasks=1,
                cpus_per_task=12,
                memory="64G",
                gpus=1,
            ),
            pixi_env="cuda-12-4",
            gmx_binary="gmx",
            mdrun_flags="-ntmpi 1 -ntomp 12 -nb gpu -pme gpu -bonded gpu -update gpu",
            module_load="module load gromacs/2024.2",
        )

        script = generator.generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/rc_scratch/run1/gromacs",
            system_prefix="CALB_PEG",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )

        assert "module load gromacs/2024.2" in script

    def test_gpu_mdrun_flags_in_script_variable(self, monkeypatch) -> None:
        """GPU offload flags should appear in the MDRUN_FLAGS variable."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )

        mdrun_flags = "-ntmpi 1 -ntomp 12 -nb gpu -pme gpu -bonded gpu -update gpu"
        generator = GromacsSlurmScriptGenerator(
            slurm_config=SlurmConfig(
                partition="blanca-shirts",
                qos="blanca-shirts",
                account="blanca-shirts",
                time_limit="23:59:59",
                nodes=1,
                ntasks=1,
                cpus_per_task=12,
                memory="64G",
                gpus=1,
            ),
            pixi_env="cuda-12-4",
            gmx_binary="gmx",
            mdrun_flags=mdrun_flags,
            module_load="module load gromacs/2024.2",
        )

        script = generator.generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/rc_scratch/run1/gromacs",
            system_prefix="CALB_PEG",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )

        assert f'MDRUN_FLAGS="{mdrun_flags}"' in script

    def test_em_uses_filtered_gpu_flags(self, monkeypatch) -> None:
        """EM stage should use MDRUN_FLAGS_EM with unsupported GPU flags removed."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )

        mdrun_flags = "-ntmpi 1 -ntomp 12 -nb gpu -pme gpu -bonded gpu -update gpu"
        generator = GromacsSlurmScriptGenerator(
            slurm_config=SlurmConfig(
                partition="blanca-shirts",
                qos="blanca-shirts",
                account="blanca-shirts",
                time_limit="23:59:59",
                nodes=1,
                ntasks=1,
                cpus_per_task=12,
                memory="64G",
                gpus=1,
            ),
            pixi_env="cuda-12-4",
            gmx_binary="gmx",
            mdrun_flags=mdrun_flags,
            module_load="module load gromacs/2024.2",
        )

        script = generator.generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/rc_scratch/run1/gromacs",
            system_prefix="CALB_PEG",
            equilibration_mdps=["eq_01_nvt.mdp", "eq_02_npt.mdp"],
        )

        assert "$MDRUN -deffnm em -cpo em.cpt $MDRUN_FLAGS_EM -v" in script
        assert "$MDRUN -deffnm eq_01 -cpo eq_01.cpt $MDRUN_FLAGS_EQ -v" in script
        assert "$MDRUN -deffnm eq_02 -cpo eq_02.cpt $MDRUN_FLAGS_EQ -v" in script
        assert "$MDRUN -deffnm prod -cpo state.cpt -maxh $MAXH $MDRUN_FLAGS_PROD -v" in script
        assert "MDRUN_FLAGS_EM=" in script
        assert "-pme" in script
        assert "-update" in script

    def test_cpu_only_em_flags_still_filtered(self, monkeypatch) -> None:
        """CPU-only scripts should still define MDRUN_FLAGS_EM for a no-op filter."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )

        generator = GromacsSlurmScriptGenerator(
            slurm_config=SlurmConfig(
                partition="aa100",
                qos="normal",
                account="ucb625_asc1",
                time_limit="23:59:59",
                nodes=1,
                ntasks=1,
                cpus_per_task=16,
                memory="3G",
                gpus=0,
            ),
            pixi_env="cuda-12-4",
            gmx_binary="gmx_mpi",
            mdrun_flags="-ntomp 16",
        )

        script = generator.generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )

        assert "$MDRUN -deffnm em -cpo em.cpt $MDRUN_FLAGS_EM -v" in script
        assert "MDRUN_FLAGS_EM=" in script


class TestCheckpointResume:
    """Tests for stage-local checkpoint resume logic."""

    def test_em_checkpoint_resume_logic(self, monkeypatch) -> None:
        """EM block should check for em.cpt and resume if found."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        block = _generator()._generate_energy_minimization_block()
        assert "-cpi em.cpt" in block
        assert "-cpo em.cpt" in block
        assert "Resuming energy minimization from checkpoint" in block

    def test_em_fresh_run_creates_checkpoint(self, monkeypatch) -> None:
        """Fresh EM run should use -cpo to create checkpoint."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        block = _generator()._generate_energy_minimization_block()
        assert "-cpo em.cpt" in block

    def test_em_skips_grompp_when_tpr_exists(self, monkeypatch) -> None:
        """EM should skip grompp if em.tpr already exists."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        block = _generator()._generate_energy_minimization_block()
        assert "if [ ! -f em.tpr ]; then" in block

    def test_equilibration_checkpoint_resume_logic(self, monkeypatch) -> None:
        """Equilibration stages should check for stage-local checkpoints."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        block = _generator()._generate_equilibration_block(["eq_01_nvt.mdp", "eq_02_npt.mdp"])
        assert "-cpi eq_01.cpt" in block
        assert "-cpo eq_01.cpt" in block
        assert "-cpi eq_02.cpt" in block
        assert "-cpo eq_02.cpt" in block
        assert "Resuming equilibration stage 1" in block
        assert "Resuming equilibration stage 2" in block

    def test_equilibration_skips_grompp_when_tpr_exists(self, monkeypatch) -> None:
        """Equilibration should skip grompp if stage.tpr exists."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        block = _generator()._generate_equilibration_block(["eq_01_nvt.mdp"])
        assert "if [ ! -f eq_01.tpr ]; then" in block

    def test_append_flag_only_on_resume(self, monkeypatch) -> None:
        """The -append flag should only appear in checkpoint resume branches."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        block = _generator()._generate_equilibration_block(["eq_01_nvt.mdp"])
        lines = block.split("\n")
        for line in lines:
            if "-append" in line:
                assert "-cpi" in line, "'-append' found without '-cpi' on same line"


class TestSetEOrdering:
    """Tests for set -e placement relative to environment setup."""

    def test_set_e_before_pixi_shell_hook(self, monkeypatch) -> None:
        """set -e must appear before pixi shell-hook."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert script.index("set -e") < script.index("pixi shell-hook")

    def test_set_e_before_module_load(self, monkeypatch) -> None:
        """set -e must appear before module load."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert script.index("set -e") < script.index("module load")


class TestGlobalTermHandling:
    """Tests for script-global SIGTERM handling."""

    def test_handle_term_function_defined(self, monkeypatch) -> None:
        """Script must define handle_term with phase-aware logic."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert "handle_term()" in script
        assert "CURRENT_PHASE=" in script
        assert "JOB_COMPLETE=" in script

    def test_trap_installed_before_simulation_work(self, monkeypatch) -> None:
        """SIGTERM trap must be installed before any grompp or mdrun."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        trap_idx = script.index("trap 'handle_term' TERM")
        assert trap_idx < script.index("run_foreground $GMX grompp")
        assert trap_idx < script.index('run_mdrun_stage "energy minimization"')

    def test_job_complete_flag_prevents_spurious_resubmit(self, monkeypatch) -> None:
        """JOB_COMPLETE=1 must be set before final exit 0."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        complete_idx = script.index("JOB_COMPLETE=1")
        done_idx = script.index("All done.")
        assert complete_idx < done_idx

    def test_foreground_term_resubmit_failure_exits_non_zero(self, monkeypatch) -> None:
        """foreground/idle TERM path should propagate resubmit_once failure code."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert "resubmit_once || true" not in script
        assert "resubmit_once\n            exit $?" in script

    def test_run_foreground_wrapper_defined(self, monkeypatch) -> None:
        """run_foreground wrapper must be defined for tracked foreground commands."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        assert "run_foreground()" in script
        assert 'CURRENT_PHASE="foreground"' in script

    def test_grompp_calls_use_run_foreground(self, monkeypatch) -> None:
        """All grompp calls should be wrapped with run_foreground."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp", "eq_02_npt.mdp"],
        )
        grompp_calls = [m.start() for m in re.finditer(r"\$GMX grompp", script)]
        assert len(grompp_calls) >= 3
        for pos in grompp_calls:
            line_start = script.rfind("\n", 0, pos) + 1
            line = script[line_start : script.find("\n", pos)]
            assert "run_foreground" in line, f"grompp call not wrapped: {line.strip()}"

    def test_trjconv_calls_use_foreground_phase(self, monkeypatch) -> None:
        """trjconv section should execute inside foreground phase tracking."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        post_idx = script.index("=== Post-processing ===")
        start_idx = script.index('CURRENT_PHASE="foreground"', post_idx)
        trjconv_idx = script.index("trjconv", post_idx)
        end_idx = script.index('CURRENT_PHASE="idle"', trjconv_idx)
        assert start_idx < trjconv_idx < end_idx

    def test_progress_update_uses_foreground_phase(self, monkeypatch) -> None:
        """Progress update calls should run inside foreground phase tracking."""
        monkeypatch.setattr(
            "polyzymd.engines.gromacs.slurm._discover_manifest_path",
            lambda: "/tmp/pixi.toml",
        )
        script = _generator().generate_job_script(
            config_path="/path/config.yaml",
            replicate=1,
            working_dir="/scratch/run1/gromacs",
            system_prefix="enzyme_polymer",
            equilibration_mdps=["eq_01_nvt.mdp"],
        )
        first_idx = script.index("_update-gromacs-progress")
        pre_idx = script.rfind('CURRENT_PHASE="foreground"', 0, first_idx)
        post_idx = script.find('CURRENT_PHASE="idle"', first_idx)
        assert pre_idx != -1
        assert post_idx != -1
        assert pre_idx < first_idx < post_idx
