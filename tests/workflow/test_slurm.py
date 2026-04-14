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

    def test_unknown_preset_raises_value_error(self):
        """Unrecognised preset names raise ValueError with available presets."""
        with pytest.raises(ValueError, match="Unknown SLURM preset 'nonexistent'"):
            SlurmConfig.from_preset("nonexistent")  # type: ignore[arg-type]

    def test_unknown_preset_error_lists_valid_presets(self):
        """The ValueError message includes all valid preset names."""
        with pytest.raises(ValueError, match="aa100") as exc_info:
            SlurmConfig.from_preset("bogus")  # type: ignore[arg-type]
        msg = str(exc_info.value)
        for name in ("aa100", "al40", "blanca-shirts", "bridges2", "testing"):
            assert name in msg, f"Missing preset {name!r} from error message"


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


class TestConstraintDirective:
    """Tests for the --constraint SBATCH directive."""

    def test_constraint_none_omits_line(self):
        """When constraint is None, no --constraint line appears."""
        cfg = SlurmConfig(constraint=None)
        gen = _make_generator(cfg)
        assert gen._constraint_line() == ""

    def test_constraint_single_value(self):
        """A single constraint value renders correctly."""
        cfg = SlurmConfig(constraint="A40")
        gen = _make_generator(cfg)
        assert gen._constraint_line() == "#SBATCH --constraint=A40"

    def test_constraint_or_expression(self):
        """SLURM OR expressions (pipe) are accepted."""
        cfg = SlurmConfig(constraint="A40|A100")
        gen = _make_generator(cfg)
        assert gen._constraint_line() == "#SBATCH --constraint=A40|A100"

    def test_constraint_and_expression(self):
        """SLURM AND expressions (ampersand) are accepted."""
        cfg = SlurmConfig(constraint="avx2&rh8")
        gen = _make_generator(cfg)
        assert gen._constraint_line() == "#SBATCH --constraint=avx2&rh8"

    def test_constraint_complex_expression(self):
        """Mixed constraint expressions are accepted."""
        cfg = SlurmConfig(constraint="A40|A100|H100")
        gen = _make_generator(cfg)
        assert gen._constraint_line() == "#SBATCH --constraint=A40|A100|H100"


# ---------------------------------------------------------------------------
# Exit code handling in JOB_TEMPLATE
# ---------------------------------------------------------------------------


class TestJobTemplateExitCodeHandling:
    """Verify that JOB_TEMPLATE handles exit codes correctly.

    The bash wrapper must:
    - Exit code 2 (concurrent): terminate cleanly without resubmitting
    - Exit code 0 or 99: proceed to check-progress / resubmit logic
    - Other non-zero: abort with error
    """

    def test_template_contains_exit_code_2_guard(self):
        """JOB_TEMPLATE must intercept exit code 2 before the generic fatal check."""
        template = SlurmScriptGenerator.JOB_TEMPLATE
        assert "if [ $RC -eq 2 ]" in template

    def test_exit_code_2_does_not_resubmit(self):
        """The exit-code-2 block must exit 0 (no resubmit) and log CONCURRENT."""
        template = SlurmScriptGenerator.JOB_TEMPLATE
        # Find the exit-code-2 block
        lines = template.splitlines()
        in_block = False
        block_lines = []
        for line in lines:
            if "if [ $RC -eq 2 ]" in line:
                in_block = True
            if in_block:
                block_lines.append(line)
                if line.strip() == "fi":
                    break

        block = "\n".join(block_lines)
        assert "exit 0" in block, "Exit code 2 block must exit 0 (clean termination)"
        assert "CONCURRENT" in block, "Exit code 2 block must log CONCURRENT message"

    def test_exit_code_2_guard_precedes_fatal_check(self):
        """Exit code 2 handling must appear BEFORE the generic non-zero check.

        If the order were reversed, exit code 2 would be caught by the
        'RC -ne 0 && RC -ne 99' guard and treated as a fatal error.
        """
        template = SlurmScriptGenerator.JOB_TEMPLATE
        pos_concurrent = template.index("if [ $RC -eq 2 ]")
        pos_fatal = template.index("if [ $RC -ne 0 ] && [ $RC -ne 99 ]")
        assert pos_concurrent < pos_fatal, (
            "Exit code 2 guard must appear before the generic fatal-error guard"
        )

    def test_template_fatal_check_excludes_code_2(self):
        """The fatal-error guard checks 'RC -ne 0 && RC -ne 99'.

        Since exit code 2 is intercepted earlier, it never reaches this guard.
        This test documents the expected pattern.
        """
        template = SlurmScriptGenerator.JOB_TEMPLATE
        assert "if [ $RC -ne 0 ] && [ $RC -ne 99 ]" in template


class TestScriptValueValidation:
    """C8-H4: SLURM script values should reject shell metacharacters."""

    def test_semicolon_injection_rejected(self):
        from polyzymd.workflow.slurm import _validate_script_value

        with pytest.raises(ValueError, match="unsafe characters"):
            _validate_script_value("gpu; rm -rf /", "partition")

    def test_pipe_injection_rejected(self):
        from polyzymd.workflow.slurm import _validate_script_value

        with pytest.raises(ValueError, match="unsafe characters"):
            _validate_script_value("gpu | cat /etc/passwd", "partition")

    def test_dollar_injection_rejected(self):
        from polyzymd.workflow.slurm import _validate_script_value

        with pytest.raises(ValueError, match="unsafe characters"):
            _validate_script_value("$HOME", "pixi_env")

    def test_backtick_injection_rejected(self):
        from polyzymd.workflow.slurm import _validate_script_value

        with pytest.raises(ValueError, match="unsafe characters"):
            _validate_script_value("`whoami`", "job_name")

    def test_newline_injection_rejected(self):
        from polyzymd.workflow.slurm import _validate_script_value

        with pytest.raises(ValueError, match="unsafe characters"):
            _validate_script_value("job\n#SBATCH --partition=evil", "job_name")

    def test_valid_partition_accepted(self):
        from polyzymd.workflow.slurm import _validate_script_value

        assert _validate_script_value("blanca-shirts", "partition") == "blanca-shirts"

    def test_comma_partition_accepted(self):
        """Comma-separated partitions like 'blanca,blanca-shirts' should work."""
        from polyzymd.workflow.slurm import _validate_script_value

        assert _validate_script_value("blanca,blanca-shirts", "partition") == "blanca,blanca-shirts"

    def test_valid_path_accepted(self):
        from polyzymd.workflow.slurm import _validate_script_value

        assert (
            _validate_script_value("/home/user/project/config.yaml", "config_path")
            == "/home/user/project/config.yaml"
        )

    def test_empty_string_accepted(self):
        """Empty strings should pass validation (e.g. empty mdrun_flags)."""
        from polyzymd.workflow.slurm import _validate_script_value

        assert _validate_script_value("", "mdrun_flags") == ""

    def test_constraint_pipe_accepted(self):
        """Pipe is valid in constraint expressions (OR)."""
        from polyzymd.workflow.slurm import _validate_constraint_value

        assert _validate_constraint_value("A40|A100", "constraint") == "A40|A100"

    def test_constraint_ampersand_accepted(self):
        """Ampersand is valid in constraint expressions (AND)."""
        from polyzymd.workflow.slurm import _validate_constraint_value

        assert _validate_constraint_value("avx2&rh8", "constraint") == "avx2&rh8"

    def test_constraint_semicolon_rejected(self):
        """Semicolons are still rejected in constraint values."""
        from polyzymd.workflow.slurm import _validate_constraint_value

        with pytest.raises(ValueError, match="unsafe characters"):
            _validate_constraint_value("A40; rm -rf /", "constraint")

    def test_constraint_dollar_rejected(self):
        """Dollar signs are rejected in constraint values."""
        from polyzymd.workflow.slurm import _validate_constraint_value

        with pytest.raises(ValueError, match="unsafe characters"):
            _validate_constraint_value("$HOME", "constraint")

    def test_constraint_backtick_rejected(self):
        """Backticks are rejected in constraint values."""
        from polyzymd.workflow.slurm import _validate_constraint_value

        with pytest.raises(ValueError, match="unsafe characters"):
            _validate_constraint_value("`whoami`", "constraint")

    def test_constraint_space_rejected(self):
        """Spaces are rejected in constraint values (not valid in SLURM constraints)."""
        from polyzymd.workflow.slurm import _validate_constraint_value

        with pytest.raises(ValueError, match="unsafe characters"):
            _validate_constraint_value("A40 A100", "constraint")
