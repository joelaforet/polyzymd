"""cuPotential job scripts: Apptainer image, plugin channel, GPU constraint, exit-86 routing."""

import os
import subprocess

import pytest

import polyzymd.workflow.slurm as slurm_module
from polyzymd.workflow.cupotential import BLANCA_CUPOTENTIAL_FEATURES, resolve_cupotential_image
from polyzymd.workflow.slurm import SlurmConfig, SlurmScriptGenerator


@pytest.fixture
def checkout(tmp_path, monkeypatch):
    """A checkout with envs/openmm-cupotential, an image and a plugin channel."""
    env_dir = tmp_path / "repo" / "envs" / "openmm-cupotential"
    env_dir.mkdir(parents=True)
    (env_dir / "pixi.toml").write_text("")
    (tmp_path / "repo.sif").write_text("")
    channel = tmp_path / "channel"
    (channel / "linux-64").mkdir(parents=True)
    (channel / "linux-64" / "repodata.json").write_text("{}")
    monkeypatch.setattr(
        slurm_module, "_discover_manifest_path", lambda: str(tmp_path / "repo" / "pixi.toml")
    )
    monkeypatch.setenv("POLYZYMD_CUPOTENTIAL_IMAGE", str(tmp_path / "repo.sif"))
    monkeypatch.setenv("POLYZYMD_CUPOTENTIAL_CHANNEL", str(channel))
    return tmp_path


def _render(tmp_path, pixi_env, preset="blanca-shirts", **config):
    config = SlurmConfig(**{**vars(SlurmConfig.from_preset(preset)), **config})
    return SlurmScriptGenerator(config, pixi_env=pixi_env, preset=preset).generate_job_script(
        config_path=str(tmp_path / "config.yaml"), replicate=1, working_dir=str(tmp_path / "run")
    )


def test_cupotential_script_runs_in_the_image_on_capable_gpus(checkout):
    script = _render(checkout, "openmm-cupotential")

    assert "#SBATCH --constraint=" + "|".join(BLANCA_CUPOTENTIAL_FEATURES) in script
    assert "V100" not in script and "P100" not in script
    assert "REQUESTED_PIXI_ENV=openmm-cupotential" in script
    assert "sim-cuda" not in script and "pixi shell-hook -e" not in script
    assert '"${CUPO_EXEC[@]}" polyzymd run-segment' in script
    assert '"${CUPO_EXEC[@]}" polyzymd check-progress' in script
    assert f"CUPO_IMAGE={checkout / 'repo.sif'}" in script
    channel_link = checkout / "repo" / "envs" / "openmm-cupotential" / "nv-channel"
    assert channel_link.resolve() == (checkout / "channel").resolve()


def test_stock_script_is_untouched(checkout):
    script = _render(checkout, "sim-cuda-12-4")

    assert "apptainer" not in script and "CUPO" not in script
    assert "--constraint" not in script
    assert "polyzymd run-segment" in script and '"${CUPO_EXEC' not in script


def test_explicit_constraint_wins(checkout):
    assert "#SBATCH --constraint=A100\n" in _render(
        checkout, "openmm-cupotential", constraint="A100"
    )


def test_missing_image_is_refused_at_submit_time(checkout, monkeypatch):
    monkeypatch.delenv("POLYZYMD_CUPOTENTIAL_IMAGE")
    with pytest.raises(ValueError, match="POLYZYMD_CUPOTENTIAL_IMAGE"):
        _render(checkout, "openmm-cupotential", preset="aa100")
    monkeypatch.setenv("POLYZYMD_CUPOTENTIAL_IMAGE", str(checkout / "absent.sif"))
    with pytest.raises(ValueError, match="does not exist"):
        resolve_cupotential_image("blanca-shirts")


def test_missing_channel_is_refused(checkout, monkeypatch):
    monkeypatch.setenv("POLYZYMD_CUPOTENTIAL_CHANNEL", str(checkout / "nowhere"))
    with pytest.raises(ValueError, match="POLYZYMD_CUPOTENTIAL_CHANNEL"):
        _render(checkout, "openmm-cupotential")


@pytest.mark.parametrize("probe_rc", [86, 0])
def test_exit_86_reroutes_and_a_good_node_runs_in_the_container(checkout, probe_rc):
    """``cupo-run`` exit 86 excludes the node; otherwise every command goes through apptainer."""
    script = checkout / "job.sh"
    script.write_text(_render(checkout, "openmm-cupotential"))
    bin_dir = checkout / "bin"
    bin_dir.mkdir()
    fakes = {
        "nvidia-smi": "echo '550.54.15, 8.6'",
        "sbatch": 'printf "%s\\n" "$@" >> "$LOG"; echo "Submitted batch job 42"',
        # Drop the options and the image; cupo-run probes, then the command runs here.
        "apptainer": 'echo "apptainer $*" >> "$LOG"\n'
        'while [ "$1" != "$IMAGE" ]; do shift; done; shift\n'
        'if [ "$1" = cupo-run ]; then [ "$2" = true ] && exit "$PROBE_RC"; shift; fi\n'
        'exec "$@"',
        "pixi": "exit 0",
        "python": 'cat > /dev/null; [ -n "$3" ] && echo "{}" > "$3"; exit 0',
        "polyzymd": 'echo "polyzymd $1" >> "$LOG"',
    }
    for name, body in fakes.items():
        (bin_dir / name).write_text(f"#!/bin/bash\n{body}\n")
        (bin_dir / name).chmod(0o755)
    log = checkout / "log"
    env = {
        **os.environ,
        "PATH": f"{bin_dir}:{os.environ['PATH']}",
        "LOG": str(log),
        "IMAGE": str(checkout / "repo.sif"),
        "PROBE_RC": str(probe_rc),
        "SLURM_JOB_ID": "7",
        "SLURMD_NODENAME": "bgpu-r550",
        "SLURM_JOB_SCRIPT": str(script),
    }

    result = subprocess.run(["bash", str(script)], env=env, capture_output=True, text=True)

    assert result.returncode == 0, result.stderr
    text = log.read_text()
    if probe_rc == 86:
        assert "cupo-run found no CUDA >= 13.4 driver API" in result.stderr
        assert "--exclude=bgpu-bortz1,bgpu-g4-u20,bgpu-g4-u24,bgpu-r550" in text
        assert "POLYZYMD_ROUTING_RETRY_COUNT=1" in text
        assert "polyzymd run-segment" not in text
    else:
        assert "pixi install --manifest-path" in text
        assert "polyzymd run-segment" in text and "polyzymd check-progress" in text
        assert "--dependency" not in text


def test_cupotential_jobs_get_a_larger_routing_retry_budget(checkout):
    """Each R525/R550 bounce excludes one node; cuPotential jobs may step past all of them."""
    from polyzymd.workflow.cupotential import CUPOTENTIAL_ROUTING_RETRY_LIMIT

    assert CUPOTENTIAL_ROUTING_RETRY_LIMIT > 3
    assert f"ROUTING_RETRY_LIMIT={CUPOTENTIAL_ROUTING_RETRY_LIMIT}\n" in _render(
        checkout, "openmm-cupotential"
    )
    assert "ROUTING_RETRY_LIMIT=3\n" in _render(checkout, "sim-cuda-12-4")


def test_container_activation_names_the_environment(checkout):
    """run-segment inside the image sees PIXI_ENVIRONMENT_NAME=openmm-cupotential, not 'default'."""
    assert "export PIXI_ENVIRONMENT_NAME=openmm-cupotential" in _render(
        checkout, "openmm-cupotential"
    )
