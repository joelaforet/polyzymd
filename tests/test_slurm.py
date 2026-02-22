"""Tests for SLURM job script generation and preset configuration.

Covers:
- SlurmConfig.from_preset() for all named presets
- Conditional SBATCH directive generation (gpu_line, qos_line, mem_line)
- Bridges2-specific GPU type and directive style
- BRIDGES2_GPU_TYPES registry
- Script generation produces well-formed output for Alpine and Bridges2 presets
- account validation guard in submit_daisy_chain()
"""

import pytest

from polyzymd.workflow.slurm import (
    BRIDGES2_GPU_TYPES,
    JobContext,
    SlurmConfig,
    SlurmScriptGenerator,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_generator(config: SlurmConfig) -> SlurmScriptGenerator:
    return SlurmScriptGenerator(config, conda_env="test-env")


def _make_context(**kwargs) -> JobContext:
    defaults = dict(
        job_name="test_job",
        output_file="logs/test.out",
        scratch_dir="/scratch/test",
        projects_dir="/projects/test",
        segment_index=0,
    )
    defaults.update(kwargs)
    return JobContext(**defaults)


# ---------------------------------------------------------------------------
# SlurmConfig.from_preset — all presets load without error
# ---------------------------------------------------------------------------


class TestPresetLoading:
    """SlurmConfig.from_preset() creates valid configs for every named preset."""

    def test_aa100_preset(self):
        cfg = SlurmConfig.from_preset("aa100")
        assert cfg.partition == "aa100"
        assert cfg.qos == "normal"
        assert cfg.account == "ucb625_asc1"
        assert cfg.time_limit == "23:59:59"
        assert cfg.memory == "3G"
        assert cfg.gpu_directive_style == "gres"
        assert cfg.gpu_type is None

    def test_al40_preset(self):
        cfg = SlurmConfig.from_preset("al40")
        assert cfg.partition == "al40"
        assert cfg.qos == "normal"
        assert cfg.account == "ucb625_asc1"
        assert cfg.gpu_directive_style == "gres"

    def test_blanca_shirts_preset(self):
        cfg = SlurmConfig.from_preset("blanca-shirts")
        assert cfg.partition == "blanca,blanca-shirts"
        assert cfg.qos == "preemptable"
        assert cfg.account == "blanca-shirts"
        assert cfg.exclude == "bgpu-bortz1"

    def test_testing_preset(self):
        cfg = SlurmConfig.from_preset("testing")
        assert cfg.partition == "atesting_a100"
        assert cfg.qos == "testing"
        assert cfg.time_limit == "0:05:59"

    def test_bridges2_preset(self):
        cfg = SlurmConfig.from_preset("bridges2")
        assert cfg.partition == "GPU-shared"
        assert cfg.qos == ""  # No QoS on Bridges2
        assert cfg.account == ""  # No shared default — user must supply
        assert cfg.time_limit == "24:00:00"
        assert cfg.memory is None  # Per-GPU allocation; omit --mem
        assert cfg.gpu_type == "v100-32"
        assert cfg.gpu_directive_style == "gpus"

    def test_preset_accepts_email(self):
        cfg = SlurmConfig.from_preset("aa100", email="user@example.com")
        assert cfg.email == "user@example.com"

    def test_unknown_preset_falls_back_to_aa100(self):
        """Unrecognised preset names silently fall back to aa100 (existing behaviour)."""
        cfg = SlurmConfig.from_preset("nonexistent")  # type: ignore[arg-type]
        assert cfg.partition == "aa100"


# ---------------------------------------------------------------------------
# BRIDGES2_GPU_TYPES registry
# ---------------------------------------------------------------------------


class TestBridges2GpuTypes:
    """BRIDGES2_GPU_TYPES list is the single source of truth for valid GPU types."""

    def test_registry_is_non_empty_list(self):
        assert isinstance(BRIDGES2_GPU_TYPES, list)
        assert len(BRIDGES2_GPU_TYPES) > 0

    def test_registry_contains_expected_types(self):
        for gpu in ("v100-16", "v100-32", "l40s-48", "h100-80"):
            assert gpu in BRIDGES2_GPU_TYPES, f"{gpu} missing from BRIDGES2_GPU_TYPES"

    def test_registry_contains_only_strings(self):
        assert all(isinstance(t, str) for t in BRIDGES2_GPU_TYPES)


# ---------------------------------------------------------------------------
# Conditional directive helpers (_gpu_line, _qos_line, _mem_line)
# ---------------------------------------------------------------------------


class TestConditionalDirectives:
    """The three conditional SBATCH helper methods produce correct output."""

    # --- GPU line ---

    def test_gpu_line_gres_style(self):
        """Default (Alpine) style emits --gres=gpu:N."""
        cfg = SlurmConfig.from_preset("aa100")
        gen = _make_generator(cfg)
        assert gen._gpu_line() == "#SBATCH --gres=gpu:1"

    def test_gpu_line_gres_style_multi_gpu(self):
        cfg = SlurmConfig(gpu_directive_style="gres", gpus=4)
        gen = _make_generator(cfg)
        assert gen._gpu_line() == "#SBATCH --gres=gpu:4"

    def test_gpu_line_gpus_style_bridges2_default(self):
        """Bridges2 preset emits --gpus=v100-32:1."""
        cfg = SlurmConfig.from_preset("bridges2")
        gen = _make_generator(cfg)
        assert gen._gpu_line() == "#SBATCH --gpus=v100-32:1"

    def test_gpu_line_gpus_style_custom_type(self):
        cfg = SlurmConfig(gpu_directive_style="gpus", gpu_type="h100-80", gpus=2)
        gen = _make_generator(cfg)
        assert gen._gpu_line() == "#SBATCH --gpus=h100-80:2"

    def test_gpu_line_gpus_style_without_gpu_type_falls_back_to_gres(self):
        """If gpu_directive_style=='gpus' but gpu_type is None, fall back to gres."""
        cfg = SlurmConfig(gpu_directive_style="gpus", gpu_type=None, gpus=1)
        gen = _make_generator(cfg)
        assert gen._gpu_line() == "#SBATCH --gres=gpu:1"

    # --- QoS line ---

    def test_qos_line_present_when_non_empty(self):
        cfg = SlurmConfig(qos="normal")
        gen = _make_generator(cfg)
        assert gen._qos_line() == "#SBATCH --qos=normal"

    def test_qos_line_omitted_when_empty(self):
        """Bridges2 sets qos='' so the directive must be absent."""
        cfg = SlurmConfig(qos="")
        gen = _make_generator(cfg)
        assert gen._qos_line() == ""

    def test_qos_line_bridges2_preset(self):
        cfg = SlurmConfig.from_preset("bridges2")
        gen = _make_generator(cfg)
        assert gen._qos_line() == ""

    # --- Memory line ---

    def test_mem_line_present_when_set(self):
        cfg = SlurmConfig(memory="4G")
        gen = _make_generator(cfg)
        assert gen._mem_line() == "#SBATCH --mem=4G"

    def test_mem_line_omitted_when_none(self):
        """Bridges2 sets memory=None so the directive must be absent."""
        cfg = SlurmConfig(memory=None)
        gen = _make_generator(cfg)
        assert gen._mem_line() == ""

    def test_mem_line_bridges2_preset(self):
        cfg = SlurmConfig.from_preset("bridges2")
        gen = _make_generator(cfg)
        assert gen._mem_line() == ""

    # --- Exclude line ---

    def test_exclude_line_present(self):
        cfg = SlurmConfig.from_preset("blanca-shirts")
        gen = _make_generator(cfg)
        assert gen._exclude_line() == "#SBATCH --exclude=bgpu-bortz1"

    def test_exclude_line_absent(self):
        cfg = SlurmConfig.from_preset("aa100")
        gen = _make_generator(cfg)
        assert gen._exclude_line() == ""


# ---------------------------------------------------------------------------
# Script generation — Alpine (gres) presets produce unchanged output
# ---------------------------------------------------------------------------


class TestAlpineScriptGeneration:
    """Generated scripts for Alpine presets are structurally correct."""

    def _gen_initial(self, preset: str = "aa100") -> str:
        cfg = SlurmConfig.from_preset(preset, email="test@example.com")
        gen = _make_generator(cfg)
        ctx = _make_context()
        return gen.generate_initial_job(
            context=ctx,
            config_path="config.yaml",
            replicate=1,
            segment_time=10.0,
            segment_frames=250,
        )

    def test_initial_script_has_gres_directive(self):
        script = self._gen_initial("aa100")
        assert "#SBATCH --gres=gpu:1" in script
        assert "--gpus=" not in script

    def test_initial_script_has_qos(self):
        script = self._gen_initial("aa100")
        assert "#SBATCH --qos=normal" in script

    def test_initial_script_has_mem(self):
        script = self._gen_initial("aa100")
        assert "#SBATCH --mem=3G" in script

    def test_initial_script_has_account(self):
        script = self._gen_initial("aa100")
        assert "#SBATCH --account=ucb625_asc1" in script

    def test_initial_script_has_shebang(self):
        script = self._gen_initial("aa100")
        assert script.strip().startswith("#!/bin/bash")

    def test_continuation_script_has_gres_directive(self):
        cfg = SlurmConfig.from_preset("aa100")
        gen = _make_generator(cfg)
        ctx = _make_context(segment_index=1)
        script = gen.generate_continuation_job(context=ctx, segment_time=10.0, num_samples=250)
        assert "#SBATCH --gres=gpu:1" in script
        assert "--gpus=" not in script


# ---------------------------------------------------------------------------
# Script generation — Bridges2 preset
# ---------------------------------------------------------------------------


class TestBridges2ScriptGeneration:
    """Generated scripts for the bridges2 preset use correct Bridges2 directives."""

    def _gen_initial(self, account: str = "abc123_gpu", gpu_type: str | None = None) -> str:
        cfg = SlurmConfig.from_preset("bridges2", email="collab@pitt.edu")
        cfg.account = account
        if gpu_type:
            cfg.gpu_type = gpu_type
        gen = _make_generator(cfg)
        ctx = _make_context()
        return gen.generate_initial_job(
            context=ctx,
            config_path="config.yaml",
            replicate=1,
            segment_time=10.0,
            segment_frames=250,
        )

    def test_uses_gpus_directive_not_gres(self):
        script = self._gen_initial()
        assert "#SBATCH --gpus=v100-32:1" in script
        assert "--gres=" not in script

    def test_no_qos_directive(self):
        script = self._gen_initial()
        assert "--qos" not in script

    def test_no_mem_directive(self):
        script = self._gen_initial()
        assert "--mem" not in script

    def test_account_in_script(self):
        script = self._gen_initial(account="myalloc_gpu")
        assert "#SBATCH --account=myalloc_gpu" in script

    def test_gpu_type_override_v100_16(self):
        script = self._gen_initial(gpu_type="v100-16")
        assert "#SBATCH --gpus=v100-16:1" in script

    def test_gpu_type_override_h100_80(self):
        script = self._gen_initial(gpu_type="h100-80")
        assert "#SBATCH --gpus=h100-80:1" in script

    def test_partition_is_gpu_shared(self):
        script = self._gen_initial()
        assert "#SBATCH --partition=GPU-shared" in script

    def test_time_limit_24h(self):
        script = self._gen_initial()
        assert "#SBATCH --time=24:00:00" in script

    def test_continuation_script_uses_gpus_directive(self):
        cfg = SlurmConfig.from_preset("bridges2")
        cfg.account = "abc123_gpu"
        gen = _make_generator(cfg)
        ctx = _make_context(segment_index=1)
        script = gen.generate_continuation_job(context=ctx, segment_time=10.0, num_samples=250)
        assert "#SBATCH --gpus=v100-32:1" in script
        assert "--gres=" not in script
        assert "--qos" not in script
        assert "--mem" not in script


# ---------------------------------------------------------------------------
# Account validation in submit_daisy_chain()
# ---------------------------------------------------------------------------


class TestAccountValidation:
    """submit_daisy_chain() enforces non-empty account for actual submission."""

    def test_empty_account_raises_on_real_submission(self, tmp_path):
        """Attempting real submission with empty account raises ValueError."""
        from unittest.mock import MagicMock, patch

        from polyzymd.workflow.daisy_chain import submit_daisy_chain

        # Minimal fake SimulationConfig with required attributes
        fake_output = MagicMock()
        fake_output.get_job_scripts_directory.return_value = tmp_path
        fake_sim_config = MagicMock()
        fake_sim_config.output = fake_output

        with (
            patch(
                "polyzymd.workflow.daisy_chain.SimulationConfig.from_yaml",
                return_value=fake_sim_config,
            ),
            pytest.raises(ValueError, match="SLURM account is required"),
        ):
            submit_daisy_chain(
                config_path=tmp_path / "config.yaml",
                slurm_preset="bridges2",
                dry_run=False,  # Real submission — must raise
            )

    def test_empty_account_warns_on_dry_run(self, tmp_path, caplog):
        """Dry-run with empty account logs a warning but does not raise."""
        import logging
        from unittest.mock import MagicMock, patch

        from polyzymd.workflow.daisy_chain import submit_daisy_chain

        fake_output = MagicMock()
        fake_output.get_job_scripts_directory.return_value = tmp_path
        fake_sim_config = MagicMock()
        fake_sim_config.output = fake_output

        fake_submitter = MagicMock()
        fake_submitter.submit_all.return_value = {}

        with (
            patch(
                "polyzymd.workflow.daisy_chain.SimulationConfig.from_yaml",
                return_value=fake_sim_config,
            ),
            patch(
                "polyzymd.workflow.daisy_chain.DaisyChainConfig.from_simulation_config",
                return_value=MagicMock(),
            ),
            patch(
                "polyzymd.workflow.daisy_chain.DaisyChainSubmitter",
                return_value=fake_submitter,
            ),
            caplog.at_level(logging.WARNING, logger="polyzymd.workflow.daisy_chain"),
        ):
            # Should NOT raise
            submit_daisy_chain(
                config_path=tmp_path / "config.yaml",
                slurm_preset="bridges2",
                dry_run=True,
            )

        assert any("SLURM account is required" in r.message for r in caplog.records)

    def test_provided_account_overrides_preset_empty(self, tmp_path):
        """Supplying --account clears the validation error for bridges2."""
        from unittest.mock import MagicMock, patch

        from polyzymd.workflow.daisy_chain import submit_daisy_chain

        fake_output = MagicMock()
        fake_output.get_job_scripts_directory.return_value = tmp_path
        fake_sim_config = MagicMock()
        fake_sim_config.output = fake_output

        fake_submitter = MagicMock()
        fake_submitter.submit_all.return_value = {}

        with (
            patch(
                "polyzymd.workflow.daisy_chain.SimulationConfig.from_yaml",
                return_value=fake_sim_config,
            ),
            patch(
                "polyzymd.workflow.daisy_chain.DaisyChainConfig.from_simulation_config",
                return_value=MagicMock(),
            ),
            patch(
                "polyzymd.workflow.daisy_chain.DaisyChainSubmitter",
                return_value=fake_submitter,
            ),
        ):
            # Should not raise — account is provided
            submit_daisy_chain(
                config_path=tmp_path / "config.yaml",
                slurm_preset="bridges2",
                account="abc123_gpu",
                dry_run=True,
            )
