"""Tests for SLURM job script generation and preset configuration.

Covers:
- SlurmConfig.from_preset() for all named presets
- Conditional SBATCH directive generation (gpu_line, qos_line, mem_line,
  nodes_line, account_line, mail_line)
- Bridges2-specific GPU type, directive style
- BRIDGES2_GPU_TYPES registry
- Job name generation for self-resubmitting model
- Pixi environment activation in generated scripts
"""

import pytest

from polyzymd.workflow.slurm import (
    BRIDGES2_GPU_TYPES,
    PRESET_DEFAULT_PIXI_ENV,
    SlurmConfig,
    SlurmScriptGenerator,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_generator(config: SlurmConfig) -> SlurmScriptGenerator:
    return SlurmScriptGenerator(config, pixi_env="cuda-12-4")


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
        assert cfg.account == ""  # Inferred from login; omit directive
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
# Conditional directive helpers
# ---------------------------------------------------------------------------


class TestConditionalDirectives:
    """All conditional SBATCH helper methods produce correct output."""

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

    # --- Nodes line ---

    def test_nodes_line_alpine_emits_two_directives(self):
        """Alpine (gres) style emits --nodes and --ntasks on separate lines."""
        cfg = SlurmConfig.from_preset("aa100")
        gen = _make_generator(cfg)
        line = gen._nodes_line()
        assert "#SBATCH --nodes=1" in line
        assert "#SBATCH --ntasks=1" in line

    def test_nodes_line_bridges2_emits_single_N_flag(self):
        """Bridges2 (gpus) style emits only '#SBATCH -N N'."""
        cfg = SlurmConfig.from_preset("bridges2")
        gen = _make_generator(cfg)
        line = gen._nodes_line()
        assert line == "#SBATCH -N 1"
        assert "--nodes" not in line
        assert "--ntasks" not in line

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

    # --- Account line ---

    def test_account_line_present_when_non_empty(self):
        cfg = SlurmConfig(account="ucb625_asc1")
        gen = _make_generator(cfg)
        assert gen._account_line() == "#SBATCH --account=ucb625_asc1"

    def test_account_line_omitted_when_empty(self):
        """Bridges2 infers allocation from login; account line must be absent."""
        cfg = SlurmConfig(account="")
        gen = _make_generator(cfg)
        assert gen._account_line() == ""

    def test_account_line_bridges2_preset_omitted(self):
        cfg = SlurmConfig.from_preset("bridges2")
        gen = _make_generator(cfg)
        assert gen._account_line() == ""

    def test_account_line_bridges2_with_user_account(self):
        """User-supplied account via --account flag appears in directive."""
        cfg = SlurmConfig.from_preset("bridges2")
        cfg.account = "chm250017p"
        gen = _make_generator(cfg)
        assert gen._account_line() == "#SBATCH --account=chm250017p"

    # --- Mail line ---

    def test_mail_line_present_when_email_set(self):
        cfg = SlurmConfig(email="user@pitt.edu")
        gen = _make_generator(cfg)
        line = gen._mail_line()
        assert "#SBATCH --mail-type=FAIL" in line
        assert "#SBATCH --mail-user=user@pitt.edu" in line

    def test_mail_line_omitted_when_no_email(self):
        """Both --mail-type and --mail-user omitted when email is empty."""
        cfg = SlurmConfig(email="")
        gen = _make_generator(cfg)
        assert gen._mail_line() == ""

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
# Job name generation
# ---------------------------------------------------------------------------


class TestJobNameGeneration:
    """Tests for DaisyChainSubmitter._create_job_name().

    Ensures the polymer composition suffix in SLURM --job-name matches the
    directory naming convention (e.g. SBMA-OEGMA_A75_B25) rather than the
    old minority-only format (SBMA-OEGMA-25%).

    The new signature is _create_job_name(replicate) — segment index was
    removed because each replicate uses a single self-resubmitting job.
    """

    def _make_submitter(self, enzyme_name, temperature, monomers_by_label):
        """Return a DaisyChainSubmitter with a mocked SimulationConfig.

        Parameters
        ----------
        enzyme_name : str
            e.g. "Fibronectin_8_to_10"
        temperature : float
            Simulation temperature in K.
        monomers_by_label : dict[str, float] | None
            Mapping of monomer label → probability, e.g. {"A": 0.75, "B": 0.25}.
            Pass None to simulate no polymer.
        """
        from unittest.mock import MagicMock

        from polyzymd.workflow.daisy_chain import DaisyChainSubmitter

        sim_config = MagicMock()
        sim_config.enzyme.name = enzyme_name
        sim_config.thermodynamics.temperature = temperature

        if monomers_by_label is None:
            sim_config.polymers = None
        else:
            sim_config.polymers.enabled = True
            sim_config.polymers.type_prefix = "SBMA-OEGMA"
            monomers = []
            for label, prob in monomers_by_label.items():
                m = MagicMock()
                m.label = label
                m.probability = prob
                monomers.append(m)
            sim_config.polymers.monomers = monomers

        dc_config = MagicMock()
        return DaisyChainSubmitter(sim_config=sim_config, dc_config=dc_config)

    def test_copolymer_uses_label_composition_format(self):
        """75/25 copolymer should produce SBMA-OEGMA_A75_B25, not SBMA-OEGMA-25%."""
        submitter = self._make_submitter("Fibronectin_8_to_10", 310.0, {"A": 0.75, "B": 0.25})
        name = submitter._create_job_name(1)
        assert "SBMA-OEGMA_A75_B25" in name
        assert "25%" not in name
        assert "-25" not in name

    def test_copolymer_label_order_is_alphabetical(self):
        """Labels must be sorted alphabetically regardless of input order."""
        submitter = self._make_submitter("Fibronectin_8_to_10", 310.0, {"B": 0.25, "A": 0.75})
        name = submitter._create_job_name(1)
        # A must come before B
        assert name.index("_A75") < name.index("_B25")

    def test_full_name_format_two_monomers(self):
        """Full job name should be r{rep}_{T}K_{enzyme}_{prefix}_{comp}."""
        submitter = self._make_submitter("Fibronectin_8_to_10", 310.0, {"A": 0.75, "B": 0.25})
        assert submitter._create_job_name(1) == ("r1_310K_Fibronectin_8_to_10_SBMA-OEGMA_A75_B25")

    def test_replicate_encoded(self):
        """Replicate index must appear in the job name."""
        submitter = self._make_submitter("Fibronectin_8_to_10", 310.0, {"A": 0.75, "B": 0.25})
        assert submitter._create_job_name(2) == ("r2_310K_Fibronectin_8_to_10_SBMA-OEGMA_A75_B25")

    def test_no_polymer_omits_composition_suffix(self):
        """When no polymer is configured, name should have no polymer suffix."""
        submitter = self._make_submitter("Fibronectin_8_to_10", 310.0, None)
        assert submitter._create_job_name(1) == "r1_310K_Fibronectin_8_to_10"

    def test_three_monomer_composition(self):
        """Three-monomer system should list all three labels in sorted order."""
        submitter = self._make_submitter("LipA", 300.0, {"A": 0.50, "B": 0.30, "C": 0.20})
        name = submitter._create_job_name(1)
        assert "SBMA-OEGMA_A50_B30_C20" in name

    def test_integer_rounding_of_probabilities(self):
        """Floating-point probabilities must round to clean integer percentages."""
        submitter = self._make_submitter("LipA", 300.0, {"A": 0.333, "B": 0.667})
        name = submitter._create_job_name(1)
        # int(0.333 * 100) = 33, int(0.667 * 100) = 66
        assert "A33" in name
        assert "B66" in name
