"""Tests for the ``polyzymd new-analysis`` scaffold command."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from polyzymd.cli.scaffold import generate_scaffold, to_pascal_case, validate_name

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
# Integration tests — generate_scaffold
# ---------------------------------------------------------------------------


class TestGenerateScaffold:
    """Scaffold file generation into a temporary directory."""

    def test_creates_three_files(self, tmp_path: Path):
        # Set up the expected directory structure
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        created = generate_scaffold("solvent_shell", tmp_path)

        assert len(created) == 3
        for p in created:
            assert p.exists(), f"{p} was not created"

    def test_file_paths(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        created = generate_scaffold("solvent_shell", tmp_path)
        names = {p.name for p in created}

        assert "__init__.py" in names
        assert "_comparison_results.py" in names
        assert "test_solvent_shell_plugin.py" in names

    def test_plugin_init_content(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        init = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        text = init.read_text()

        assert "class SolventShellAnalysis(Analysis):" in text
        assert 'name: ClassVar[str] = "solvent_shell"' in text
        assert "class SolventShellSettings(BaseModel):" in text
        assert "def compute_replicate(" in text
        assert "def aggregate(" in text
        assert "def extract_metrics(" in text

    def test_comparison_results_content(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        cr = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "_comparison_results.py"
        text = cr.read_text()

        assert "class SolventShellComparisonResult(ComparisonResult):" in text
        assert 'analysis_type: str = "solvent_shell"' in text

    def test_test_file_content(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        tf = tmp_path / "tests" / "test_solvent_shell_plugin.py"
        text = tf.read_text()

        assert "class TestDiscovery:" in text
        assert "class TestComputeReplicate:" in text
        assert "class TestAggregate:" in text
        assert "class TestExtractMetrics:" in text

    def test_custom_class_name(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("density", tmp_path, class_name="MassDensity")
        init = tmp_path / "src" / "polyzymd" / "analyses" / "density" / "__init__.py"
        text = init.read_text()

        assert "class MassDensityAnalysis(Analysis):" in text
        assert "class MassDensitySettings(BaseModel):" in text

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
        assert len(created) == 3

    def test_dry_run_creates_no_files(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        created = generate_scaffold("solvent_shell", tmp_path, dry_run=True)
        assert len(created) == 3
        for p in created:
            assert not p.exists(), f"{p} should not exist in dry-run mode"


# ---------------------------------------------------------------------------
# Integration test — generated plugin is syntactically valid Python
# ---------------------------------------------------------------------------


class TestGeneratedCodeQuality:
    """Ensure generated code compiles and passes basic checks."""

    def test_generated_init_compiles(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        init = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        source = init.read_text()
        compile(source, str(init), "exec")  # raises SyntaxError if invalid

    def test_generated_results_compiles(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        cr = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "_comparison_results.py"
        source = cr.read_text()
        compile(source, str(cr), "exec")

    def test_generated_test_compiles(self, tmp_path: Path):
        (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
        (tmp_path / "tests").mkdir()

        generate_scaffold("solvent_shell", tmp_path)
        tf = tmp_path / "tests" / "test_solvent_shell_plugin.py"
        source = tf.read_text()
        compile(source, str(tf), "exec")
