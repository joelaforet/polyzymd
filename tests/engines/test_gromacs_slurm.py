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
