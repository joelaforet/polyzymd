"""Tests for GROMACS SLURM script generation."""

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


def test_em_block_two_stage() -> None:
    """Two-stage EM block should include soft + standard EM and in-block health check."""
    block = _generator()._generate_energy_minimization_block(use_soft_em=True)

    assert "em_soft.mdp" in block
    assert "em_soft.gro" in block
    assert "em.mdp" in block
    assert "em.gro" in block
    assert "force.*not finite" in block

    em_if_idx = block.index("if [ ! -f em.gro ]; then")
    health_idx = block.index("force.*not finite")
    em_else_idx = block.index(
        'else\n    echo "Skipping standard energy minimization (em.gro exists)."'
    )
    assert em_if_idx < health_idx < em_else_idx


def test_em_block_legacy() -> None:
    """Legacy EM block should not include soft-stage artifacts."""
    block = _generator()._generate_energy_minimization_block(use_soft_em=False)

    assert "em_soft" not in block
    assert "em.mdp" in block
    assert "em.gro" in block
    assert "force.*not finite" in block


def test_use_soft_em_parameter(monkeypatch) -> None:
    """Job script should omit soft EM when use_soft_em is False."""
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
        use_soft_em=False,
    )

    assert "em_soft" not in script


def test_use_soft_em_default(monkeypatch) -> None:
    """Job script should include soft EM by default."""
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

    assert "em_soft" in script


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
    """MDRUN_FLAGS should be used in EM, equilibration, and production."""
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

    assert "$GMX mdrun -deffnm em_soft $MDRUN_FLAGS -v" in script
    assert "$GMX mdrun -deffnm em $MDRUN_FLAGS -v" in script
    assert "$GMX mdrun -deffnm eq_01 $MDRUN_FLAGS -v" in script
    assert "$GMX mdrun -deffnm eq_02 $MDRUN_FLAGS -v" in script
    assert "$GMX mdrun -deffnm prod -cpo state.cpt -maxh $MAXH $MDRUN_FLAGS -v" in script
    assert script.count("$MDRUN_FLAGS") >= 5


def test_mdrun_flags_variable_defined_before_em(monkeypatch) -> None:
    """MDRUN_FLAGS variable should be declared before EM block commands."""
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

    assert script.index("MDRUN_FLAGS=") < script.index("em_soft.mdp")
