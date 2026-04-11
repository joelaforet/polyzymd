"""Tests for the ``polyzymd new-analysis`` scaffold command."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from polyzymd.cli.scaffold import (
    VALID_STYLES,
    generate_scaffold,
    to_pascal_case,
    validate_class_name,
    validate_name,
)

# ---------------------------------------------------------------------------
# Unit tests — validate_name
# ---------------------------------------------------------------------------


class TestValidateName:
    """Validation of proposed plugin names."""

    def test_valid_simple(self):
        assert validate_name("solvent_shell", check_existing=False) is None

    def test_valid_single_word(self):
        assert validate_name("density", check_existing=False) is None

    def test_valid_with_digits(self):
        assert validate_name("rdf2d", check_existing=False) is None

    def test_reject_uppercase(self):
        err = validate_name("SolventShell", check_existing=False)
        assert err is not None
        assert "snake_case" in err

    def test_reject_leading_digit(self):
        err = validate_name("2drdf", check_existing=False)
        assert err is not None

    def test_reject_hyphen(self):
        err = validate_name("solvent-shell", check_existing=False)
        assert err is not None

    def test_reject_python_keyword(self):
        err = validate_name("class", check_existing=False)
        assert err is not None
        assert "keyword" in err

    def test_reject_reserved_name(self):
        for reserved in ("base", "discovery", "orchestrator", "shared", "stats"):
            err = validate_name(reserved, check_existing=False)
            assert err is not None
            assert "reserved" in err

    def test_reject_existing_plugin(self):
        """Collision check against registered plugins (requires discovery)."""
        err = validate_name("rmsf", check_existing=True)
        assert err is not None
        assert "already exists" in err

    def test_accept_novel_name_with_collision_check(self):
        assert validate_name("my_brand_new_analysis_xyz", check_existing=True) is None


# ---------------------------------------------------------------------------
# Unit tests — to_pascal_case
# ---------------------------------------------------------------------------


class TestToPascalCase:
    def test_single_word(self):
        assert to_pascal_case("density") == "Density"

    def test_two_words(self):
        assert to_pascal_case("solvent_shell") == "SolventShell"

    def test_three_words(self):
        assert to_pascal_case("binding_free_energy") == "BindingFreeEnergy"


# ---------------------------------------------------------------------------
# Integration tests — generate_scaffold (default dict style)
# ---------------------------------------------------------------------------


class TestGenerateScaffold:
    """Scaffold file generation using the default dict style."""

    def test_creates_two_files(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        created = generate_scaffold("solvent_shell", tmp_path)

        assert len(created) == 2
        for p in created:
            assert p.exists(), f"{p} was not created"

    def test_file_paths(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        created = generate_scaffold("solvent_shell", tmp_path)
        names = {p.name for p in created}

        assert "__init__.py" in names
        assert "test_solvent_shell.py" in names

    def test_plugin_init_content(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        init = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        text = init.read_text()

        # Dict-style: class exists, has name, Settings, lifecycle methods
        assert "class SolventShellAnalysis(Analysis):" in text
        assert 'name: ClassVar[str] = "solvent_shell"' in text
        assert "class Settings(BaseModel):" in text
        assert "def compute_replicate(" in text
        assert "def aggregate(" in text
        assert "def extract_metrics(" in text
        assert "def plot(self, ctx: PlotContext)" in text
        assert "_build_plot_data" in text

    def test_dict_style_uses_plain_dicts(self, tmp_path: Path):
        """Dict-style scaffold should NOT have AggregatedResultClass or typed result models."""
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        init = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        text = init.read_text()

        # No Pydantic result models in dict style
        assert "class SolventShellReplicateResult" not in text
        assert "class SolventShellAggregatedResult" not in text
        # No AggregatedResultClass ClassVar assignment (may appear in comments)
        assert "AggregatedResultClass: ClassVar" not in text
        # Return type annotations should be dict
        assert "-> dict[str, Any]:" in text

    def test_test_file_content(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        tf = tmp_path / "tests" / "analyses" / "plugins" / "test_solvent_shell.py"
        text = tf.read_text()

        assert "class TestDiscovery:" in text
        assert "class TestComputeReplicate:" in text
        assert "class TestAggregate:" in text
        assert "class TestExtractMetrics:" in text
        assert "class TestPlot:" in text

    def test_custom_class_name(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("density", tmp_path, class_name="MassDensity")
        init = tmp_path / "src" / "polyzymd" / "analyses" / "density" / "__init__.py"
        text = init.read_text()

        assert "class MassDensityAnalysis(Analysis):" in text
        assert "class Settings(BaseModel):" in text

    def test_refuses_overwrite_without_force(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        with pytest.raises(FileExistsError, match="already exists"):
            generate_scaffold("solvent_shell", tmp_path)

    def test_force_overwrites(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        # Should not raise
        created = generate_scaffold("solvent_shell", tmp_path, force=True)
        assert len(created) == 2

    def test_dry_run_creates_no_files(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        created = generate_scaffold("solvent_shell", tmp_path, dry_run=True)
        assert len(created) == 2
        for p in created:
            assert not p.exists(), f"{p} should not exist in dry-run mode"

    def test_default_style_is_dict(self, tmp_path: Path):
        """Verify that omitting --style gives dict-style output."""
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        init = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        text = init.read_text()

        assert "Uses plain dicts for results" in text
        assert "AggregatedResultClass: ClassVar" not in text


# ---------------------------------------------------------------------------
# Integration tests — generate_scaffold (pydantic style)
# ---------------------------------------------------------------------------


class TestGenerateScaffoldPydantic:
    """Scaffold file generation using the pydantic style."""

    def test_creates_two_files(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        created = generate_scaffold("solvent_shell", tmp_path, style="pydantic")

        assert len(created) == 2
        for p in created:
            assert p.exists(), f"{p} was not created"

    def test_plugin_init_content(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path, style="pydantic")
        init = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        text = init.read_text()

        assert "class SolventShellAnalysis(Analysis):" in text
        assert 'name: ClassVar[str] = "solvent_shell"' in text
        assert "AggregatedResultClass: ClassVar[type] = SolventShellAggregatedResult" in text
        assert "class Settings(BaseModel):" in text
        assert "def compute_replicate(" in text
        assert "def aggregate(" in text
        assert "def extract_metrics(" in text
        assert "def plot(self, ctx: PlotContext)" in text
        assert "_build_plot_data" in text

    def test_init_contains_result_models(self, tmp_path: Path):
        """Result models (ReplicateResult, AggregatedResult) are inline in __init__.py."""
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path, style="pydantic")
        init = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        text = init.read_text()

        assert "class SolventShellReplicateResult(BaseModel):" in text
        assert "class SolventShellAggregatedResult(BaseModel):" in text

    def test_test_file_content(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path, style="pydantic")
        tf = tmp_path / "tests" / "analyses" / "plugins" / "test_solvent_shell.py"
        text = tf.read_text()

        assert "class TestDiscovery:" in text
        assert "class TestComputeReplicate:" in text
        assert "class TestAggregate:" in text
        assert "class TestExtractMetrics:" in text
        assert "class TestPlot:" in text
        # Pydantic test imports the typed result models
        assert "SolventShellReplicateResult" in text
        assert "SolventShellAggregatedResult" in text

    def test_custom_class_name(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("density", tmp_path, class_name="MassDensity", style="pydantic")
        init = tmp_path / "src" / "polyzymd" / "analyses" / "density" / "__init__.py"
        text = init.read_text()

        assert "class MassDensityAnalysis(Analysis):" in text
        assert "class Settings(BaseModel):" in text
        assert "class MassDensityReplicateResult(BaseModel):" in text
        assert "class MassDensityAggregatedResult(BaseModel):" in text

    def test_force_overwrites(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path, style="pydantic")
        created = generate_scaffold("solvent_shell", tmp_path, style="pydantic", force=True)
        assert len(created) == 2


# ---------------------------------------------------------------------------
# Style validation
# ---------------------------------------------------------------------------


class TestStyleValidation:
    """Test style parameter validation."""

    def test_valid_styles_constant(self):
        assert "dict" in VALID_STYLES
        assert "pydantic" in VALID_STYLES

    def test_invalid_style_raises(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        with pytest.raises(ValueError, match="Invalid style"):
            generate_scaffold("solvent_shell", tmp_path, style="bad_style")


# ---------------------------------------------------------------------------
# Integration test — generated code is syntactically valid Python
# ---------------------------------------------------------------------------


class TestGeneratedCodeQuality:
    """Ensure generated code compiles and passes basic checks."""

    def test_dict_init_compiles(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        init = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        source = init.read_text()
        compile(source, str(init), "exec")  # raises SyntaxError if invalid

    def test_dict_test_compiles(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        tf = tmp_path / "tests" / "analyses" / "plugins" / "test_solvent_shell.py"
        source = tf.read_text()
        compile(source, str(tf), "exec")

    def test_pydantic_init_compiles(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path, style="pydantic")
        init = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        source = init.read_text()
        compile(source, str(init), "exec")

    def test_pydantic_init_contains_result_models(self, tmp_path: Path):
        """Pydantic result models are inline in __init__.py; verify it compiles."""
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path, style="pydantic")
        init = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        source = init.read_text()
        assert "ReplicateResult" in source
        assert "AggregatedResult" in source
        compile(source, str(init), "exec")

    def test_pydantic_test_compiles(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path, style="pydantic")
        tf = tmp_path / "tests" / "analyses" / "plugins" / "test_solvent_shell.py"
        source = tf.read_text()
        compile(source, str(tf), "exec")


# ---------------------------------------------------------------------------
# Unit tests — validate_class_name
# ---------------------------------------------------------------------------


class TestValidateClassName:
    """Validation of custom class name prefixes."""

    def test_valid_pascal_case(self):
        assert validate_class_name("SolventShell") is None

    def test_valid_single_word(self):
        assert validate_class_name("Density") is None

    def test_valid_with_digits(self):
        assert validate_class_name("Rdf2D") is None

    def test_reject_leading_digit(self):
        err = validate_class_name("2DRdf")
        assert err is not None
        assert "identifier" in err

    def test_reject_python_keyword(self):
        err = validate_class_name("class")
        assert err is not None
        assert "keyword" in err

    def test_reject_lowercase_start(self):
        err = validate_class_name("solventShell")
        assert err is not None
        assert "uppercase" in err

    def test_reject_hyphenated(self):
        err = validate_class_name("Solvent-Shell")
        assert err is not None
        assert "identifier" in err

    def test_reject_spaces(self):
        err = validate_class_name("Solvent Shell")
        assert err is not None
        assert "identifier" in err

    def test_reject_empty_string(self):
        err = validate_class_name("")
        assert err is not None


class TestGenerateScaffoldClassNameValidation:
    """Ensure generate_scaffold rejects invalid class names."""

    def test_invalid_class_name_raises(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        with pytest.raises(ValueError, match="identifier"):
            generate_scaffold("density", tmp_path, class_name="2Bad")

    def test_lowercase_class_name_raises(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        with pytest.raises(ValueError, match="uppercase"):
            generate_scaffold("density", tmp_path, class_name="massDensity")


# ---------------------------------------------------------------------------
# CLI integration tests — CliRunner
# ---------------------------------------------------------------------------


class TestNewAnalysisCLI:
    """Test the ``new-analysis`` Click command via CliRunner."""

    @pytest.fixture(autouse=True)
    def _setup_project_dir(self, tmp_path: Path):
        """Create a minimal project structure and store the root."""
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()
        self.root = tmp_path

    @pytest.fixture
    def runner(self):
        return CliRunner()

    @pytest.fixture
    def cli(self):
        from polyzymd.cli.main import new_analysis

        return new_analysis

    def test_success(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["solvent_shell", "--project-root", str(self.root)],
        )
        assert result.exit_code == 0, result.output
        assert "scaffolded successfully" in result.output
        init = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        assert init.exists()

    def test_default_style_is_dict(self, runner: CliRunner, cli):
        """Default invocation (no --style) produces dict-style output."""
        result = runner.invoke(
            cli,
            ["solvent_shell", "--project-root", str(self.root)],
        )
        assert result.exit_code == 0, result.output
        init = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        text = init.read_text()
        assert "Uses plain dicts for results" in text
        assert "AggregatedResultClass: ClassVar" not in text

    def test_style_pydantic(self, runner: CliRunner, cli):
        """--style pydantic produces pydantic-style output."""
        result = runner.invoke(
            cli,
            ["solvent_shell", "--style", "pydantic", "--project-root", str(self.root)],
        )
        assert result.exit_code == 0, result.output
        init = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        text = init.read_text()
        assert "AggregatedResultClass" in text
        assert "class SolventShellReplicateResult(BaseModel):" in text

    def test_style_dict_explicit(self, runner: CliRunner, cli):
        """--style dict explicitly produces dict-style output."""
        result = runner.invoke(
            cli,
            ["solvent_shell", "--style", "dict", "--project-root", str(self.root)],
        )
        assert result.exit_code == 0, result.output
        init = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        text = init.read_text()
        assert "Uses plain dicts for results" in text

    def test_invalid_style_rejected(self, runner: CliRunner, cli):
        """--style with an invalid value is rejected by Click."""
        result = runner.invoke(
            cli,
            ["solvent_shell", "--style", "invalid", "--project-root", str(self.root)],
        )
        assert result.exit_code != 0
        # Click Choice validation produces an error message
        assert "Invalid value" in result.output or "invalid" in result.output.lower()

    def test_dry_run(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["solvent_shell", "--project-root", str(self.root), "--dry-run"],
        )
        assert result.exit_code == 0, result.output
        assert "Would create" in result.output
        init = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        assert not init.exists()

    def test_invalid_name(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["BadName", "--project-root", str(self.root)],
        )
        assert result.exit_code != 0
        assert "snake_case" in result.output

    def test_existing_plugin_rejected(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["rmsf", "--project-root", str(self.root)],
        )
        assert result.exit_code != 0
        assert "already exists" in result.output

    def test_invalid_class_name(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["density", "--class-name", "2Bad", "--project-root", str(self.root)],
        )
        assert result.exit_code != 0
        assert "identifier" in result.output

    def test_lowercase_class_name_rejected(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["density", "--class-name", "massDensity", "--project-root", str(self.root)],
        )
        assert result.exit_code != 0
        assert "uppercase" in result.output

    def test_custom_class_name_success(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["density", "--class-name", "MassDensity", "--project-root", str(self.root)],
        )
        assert result.exit_code == 0, result.output
        init = self.root / "src" / "polyzymd" / "analyses" / "density" / "__init__.py"
        assert init.exists()
        assert "class MassDensityAnalysis(Analysis):" in init.read_text()

    def test_force_overwrites(self, runner: CliRunner, cli):
        # First run
        result = runner.invoke(
            cli,
            ["solvent_shell", "--project-root", str(self.root)],
        )
        assert result.exit_code == 0

        # Second run without force — should fail
        result = runner.invoke(
            cli,
            ["solvent_shell", "--project-root", str(self.root)],
        )
        assert result.exit_code != 0
        assert "already exists" in result.output

        # Third run with force — should succeed
        result = runner.invoke(
            cli,
            ["solvent_shell", "--project-root", str(self.root), "--force"],
        )
        assert result.exit_code == 0, result.output

    def test_help_shows_correct_test_path(self, runner: CliRunner, cli):
        """Help text should reference the nested test layout."""
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "tests/analyses/plugins/test_<NAME>.py" in result.output

    def test_success_message_shows_correct_test_path(self, runner: CliRunner, cli):
        """Success message should reference the correct pytest command."""
        result = runner.invoke(
            cli,
            ["solvent_shell", "--project-root", str(self.root)],
        )
        assert result.exit_code == 0, result.output
        assert "tests/analyses/plugins/test_solvent_shell.py" in result.output
