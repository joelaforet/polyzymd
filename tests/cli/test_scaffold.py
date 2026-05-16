"""Tests for the ``polyzymd new-analysis`` scaffold command."""

from __future__ import annotations

import subprocess
import sys
from importlib.resources import files
from pathlib import Path
from types import SimpleNamespace

import black
import pytest
from click.testing import CliRunner
from jinja2.exceptions import UndefinedError

from polyzymd.cli._scaffold.models import ScaffoldSpec
from polyzymd.cli._scaffold.renderer import create_environment, render_template
from polyzymd.cli.scaffold import (
    VALID_STYLES,
    generate_scaffold,
    to_pascal_case,
    validate_class_name,
    validate_name,
)

FakeSimulationConfig = SimpleNamespace
REPO_ROOT = Path(__file__).resolve().parents[2]


def _prepare_project(tmp_path: Path) -> None:
    """Create the minimal project directories required by the scaffold."""
    (tmp_path / "src" / "polyzymd" / "analyses").mkdir(parents=True)
    (tmp_path / "tests").mkdir()


def _generated_py_files(tmp_path: Path) -> list[Path]:
    """Return generated Python files in a deterministic order."""
    return sorted(tmp_path.rglob("*.py"))


def _assert_black_compliant(path: Path) -> None:
    """Assert that a Python file already satisfies Black formatting."""
    source = path.read_text(encoding="utf-8")
    mode = black.FileMode(line_length=100)
    try:
        black.format_file_contents(source, fast=True, mode=mode)
    except black.NothingChanged:
        return
    pytest.fail(f"Generated file is not Black formatted: {path}")


def _assert_ruff_compliant(paths: list[Path]) -> None:
    """Assert that generated Python files already satisfy the project Ruff config."""
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "ruff",
            "check",
            "--config",
            str(REPO_ROOT / "pyproject.toml"),
            "--config",
            'lint.isort.known-first-party=["polyzymd"]',
            *(str(path) for path in paths),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(
            "Generated files are not Ruff compliant:\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


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
        assert validate_name("2drdf", check_existing=False) is not None

    def test_reject_hyphen(self):
        assert validate_name("solvent-shell", check_existing=False) is not None

    def test_reject_python_keyword(self):
        err = validate_name("class", check_existing=False)
        assert err is not None
        assert "keyword" in err

    def test_reject_reserved_name(self):
        for reserved in ("base", "discovery", "mda", "orchestrator", "shared", "stats"):
            err = validate_name(reserved, check_existing=False)
            assert err is not None
            assert "reserved" in err

    def test_reject_existing_plugin(self):
        err = validate_name("rmsf", check_existing=True)
        assert err is not None
        assert "already exists" in err

    def test_accept_novel_name_with_collision_check(self):
        assert validate_name("my_brand_new_analysis_xyz", check_existing=True) is None


class TestToPascalCase:
    """Conversion from plugin names to class prefixes."""

    def test_single_word(self):
        assert to_pascal_case("density") == "Density"

    def test_two_words(self):
        assert to_pascal_case("solvent_shell") == "SolventShell"

    def test_three_words(self):
        assert to_pascal_case("solvent_shell_density") == "SolventShellDensity"


class TestGenerateScaffold:
    """Scaffold file generation using the default measurement style."""

    def test_creates_two_files(self, tmp_path: Path):
        _prepare_project(tmp_path)

        created = generate_scaffold("solvent_shell", tmp_path)

        assert len(created) == 2
        for path in created:
            assert path.exists(), f"{path} was not created"

    def test_file_paths(self, tmp_path: Path):
        _prepare_project(tmp_path)

        created = generate_scaffold("solvent_shell", tmp_path)
        names = {path.name for path in created}

        assert names == {"solvent_shell.py", "test_solvent_shell.py"}

    def test_measurement_plugin_content(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        plugin_path = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell.py"
        text = plugin_path.read_text(encoding="utf-8")

        assert "class SolventShellSettings(BaseModel):" in text
        assert "scale: float = Field(" in text
        assert "class SolventShellMeasurement(ScalarMeasurement):" in text
        assert "metric: ClassVar[MetricSpec] = MetricSpec(" in text
        assert "def measure(" in text
        assert "class SolventShellAnalysis(ScalarMeasurementAnalysis):" in text
        assert 'name: ClassVar[str] = "solvent_shell"' in text
        assert "Settings: ClassVar[type[BaseModel]] = SolventShellSettings" in text
        assert "measurement: ClassVar[type[ScalarMeasurement]] = SolventShellMeasurement" in text
        assert "def build_runner(" not in text
        assert "def aggregate(" not in text
        assert "def plot(" not in text
        assert "._runner" not in text

    def test_default_does_not_create_advanced_package_files(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        analyses_root = tmp_path / "src" / "polyzymd" / "analyses"

        assert (analyses_root / "solvent_shell.py").exists()
        assert not (analyses_root / "solvent_shell").exists()

    def test_test_file_content(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        test_path = tmp_path / "tests" / "analyses" / "plugins" / "test_solvent_shell.py"
        text = test_path.read_text(encoding="utf-8")

        assert "class TestDiscovery:" in text
        assert "class TestSettings:" in text
        assert "class TestMeasurementLogic:" in text
        assert "class TestMetricMetadata:" in text
        assert "test_measure_counts_windowed_frames_with_scale" in text
        assert "AggregateContext" not in text
        assert "ReplicateContext" not in text
        assert "PlotContext" not in text

    def test_custom_class_name(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("density", tmp_path, class_name="MassDensity")
        plugin_path = tmp_path / "src" / "polyzymd" / "analyses" / "density.py"

        assert "class MassDensitySettings(BaseModel):" in plugin_path.read_text(encoding="utf-8")
        assert "class MassDensityMeasurement(ScalarMeasurement):" in plugin_path.read_text(
            encoding="utf-8"
        )
        assert "class MassDensityAnalysis(ScalarMeasurementAnalysis):" in plugin_path.read_text(
            encoding="utf-8"
        )

    def test_refuses_overwrite_without_force(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        with pytest.raises(FileExistsError, match="already exists"):
            generate_scaffold("solvent_shell", tmp_path)

    def test_preflights_default_test_file_collision_before_writing_source(
        self,
        tmp_path: Path,
    ):
        _prepare_project(tmp_path)
        test_path = tmp_path / "tests" / "analyses" / "plugins" / "test_solvent_shell.py"
        test_path.parent.mkdir(parents=True)
        test_path.write_text("# existing test\n", encoding="utf-8")
        plugin_path = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell.py"

        with pytest.raises(FileExistsError, match="test_solvent_shell.py"):
            generate_scaffold("solvent_shell", tmp_path)

        assert not plugin_path.exists()
        assert test_path.read_text(encoding="utf-8") == "# existing test\n"

    @pytest.mark.parametrize(
        "collision_parts",
        [
            ("analyses",),
            ("analyses", "plugins"),
        ],
    )
    def test_preflights_default_test_parent_collision_before_writing_source(
        self,
        tmp_path: Path,
        collision_parts: tuple[str, ...],
    ):
        _prepare_project(tmp_path)
        collision_path = tmp_path / "tests" / Path(*collision_parts)
        collision_path.parent.mkdir(parents=True, exist_ok=True)
        collision_path.write_text("# parent collision\n", encoding="utf-8")
        plugin_path = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell.py"
        test_path = tmp_path / "tests" / "analyses" / "plugins" / "test_solvent_shell.py"

        with pytest.raises(FileExistsError, match="parent path is not a directory"):
            generate_scaffold("solvent_shell", tmp_path)

        assert not plugin_path.exists()
        assert not test_path.exists()
        assert collision_path.read_text(encoding="utf-8") == "# parent collision\n"

    def test_force_overwrites(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        plugin_path = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell.py"
        plugin_path.write_text(
            plugin_path.read_text(encoding="utf-8") + "\n# stale local edit\n",
            encoding="utf-8",
        )
        created = generate_scaffold("solvent_shell", tmp_path, force=True)
        assert len(created) == 2
        assert "stale local edit" not in plugin_path.read_text(encoding="utf-8")

    def test_dry_run_creates_no_files(self, tmp_path: Path):
        _prepare_project(tmp_path)

        created = generate_scaffold("solvent_shell", tmp_path, dry_run=True)
        assert len(created) == 2
        for path in created:
            assert not path.exists(), f"{path} should not exist in dry-run mode"

    def test_default_style_is_measurement(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        plugin_path = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell.py"
        text = plugin_path.read_text(encoding="utf-8")

        assert "ScalarMeasurementAnalysis" in text
        assert "SolventShellReplicateRunner" not in text


class TestGenerateScaffoldAdvancedDict:
    """Advanced dict scaffolds should not generate old runner templates."""

    def test_style_dict_is_temporarily_unavailable(self, tmp_path: Path):
        _prepare_project(tmp_path)

        with pytest.raises(
            ValueError, match="Advanced analysis scaffolds are temporarily unavailable"
        ):
            generate_scaffold("solvent_shell", tmp_path, style="dict")

    def test_advanced_flag_defaults_to_dict_style(self, tmp_path: Path):
        _prepare_project(tmp_path)

        with pytest.raises(
            ValueError, match="Advanced analysis scaffolds are temporarily unavailable"
        ):
            generate_scaffold("solvent_shell", tmp_path, advanced=True)


class TestGenerateScaffoldPydantic:
    """Pydantic scaffolds should wait for the MDA-native advanced template."""

    def test_style_pydantic_is_temporarily_unavailable(self, tmp_path: Path):
        _prepare_project(tmp_path)

        with pytest.raises(
            ValueError, match="Advanced analysis scaffolds are temporarily unavailable"
        ):
            generate_scaffold("solvent_shell", tmp_path, style="pydantic")


class TestStyleValidation:
    """Style parameter validation."""

    def test_valid_styles_constant(self):
        assert "measurement" in VALID_STYLES
        assert "dict" in VALID_STYLES
        assert "pydantic" in VALID_STYLES

    def test_invalid_style_raises(self, tmp_path: Path):
        _prepare_project(tmp_path)

        with pytest.raises(ValueError, match="Invalid style"):
            generate_scaffold("solvent_shell", tmp_path, style="bad_style")


class TestLayoutCollisionChecks:
    """Cross-layout collision behavior for module and package scaffolds."""

    def test_measurement_refuses_existing_package_even_with_force(self, tmp_path: Path):
        _prepare_project(tmp_path)
        package_dir = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell"
        package_dir.mkdir()

        with pytest.raises(FileExistsError, match="package layout"):
            generate_scaffold("solvent_shell", tmp_path, force=True)

    def test_advanced_refuses_existing_module_even_with_force(self, tmp_path: Path):
        _prepare_project(tmp_path)
        module_path = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell.py"
        module_path.write_text("", encoding="utf-8")

        with pytest.raises(
            ValueError, match="Advanced analysis scaffolds are temporarily unavailable"
        ):
            generate_scaffold("solvent_shell", tmp_path, style="dict", force=True)


class TestTemplateResources:
    """Jinja package-resource template behavior."""

    def test_template_resources_are_present(self):
        template_root = files("polyzymd.cli._scaffold") / "templates"
        names = {path.name for path in template_root.iterdir()}

        assert {
            "measurement_plugin.py.jinja",
            "test_measurement_plugin.py.jinja",
        }.issubset(names)

    def test_renderer_uses_strict_undefined(self):
        env = create_environment()

        with pytest.raises(UndefinedError):
            env.from_string("{{ missing_value }}").render()

    def test_templates_render_with_sample_context(self):
        measurement_spec = ScaffoldSpec(
            name="strict_sample",
            class_name="StrictSample",
            style="measurement",
        )
        for template_name in ("measurement_plugin.py.jinja", "test_measurement_plugin.py.jinja"):
            rendered = render_template(template_name, measurement_spec)
            assert "StrictSample" in rendered or "strict_sample" in rendered


class TestGeneratedCodeQuality:
    """Ensure generated code compiles and satisfies formatters."""

    @pytest.mark.parametrize(
        ("name", "style", "expected_count"),
        [
            ("solvent_shell", "measurement", 2),
            ("scaffold_pydantic_e2e", "measurement", 2),
        ],
    )
    def test_generated_files_compile_and_are_formatter_clean(
        self,
        tmp_path: Path,
        name: str,
        style: str,
        expected_count: int,
    ):
        _prepare_project(tmp_path)

        generate_scaffold(name, tmp_path, style=style)

        generated = _generated_py_files(tmp_path)
        assert len(generated) == expected_count
        for path in generated:
            compile(path.read_text(encoding="utf-8"), str(path), "exec")
            _assert_black_compliant(path)
        _assert_ruff_compliant(generated)


class TestGeneratedPluginEndToEnd:
    """End-to-end checks for scaffolded plugin discovery and lifecycle."""

    class FakeTrajectory:
        """Minimal trajectory with a length for measurement tests."""

        def __len__(self) -> int:
            return 5

    class FakeUniverse:
        """Minimal Universe exposing a trajectory."""

        def __init__(self) -> None:
            self.trajectory = TestGeneratedPluginEndToEnd.FakeTrajectory()

    class FakeWindow:
        """Minimal trajectory window for scalar measurements."""

        warning_message: str | None = None

        def run_kwargs(self) -> dict[str, int | None]:
            return {"start": 1, "stop": 5, "step": 2}

    class FakeLoader:
        """Minimal loader that returns a fake Universe."""

        def __init__(self, sim_config):
            self.sim_config = sim_config

        def load_universe(self, replicate: int):
            self.replicate = replicate
            return TestGeneratedPluginEndToEnd.FakeUniverse()

    def test_measurement_scaffold_is_discoverable_and_measurable(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        analyses_root = tmp_path / "src" / "polyzymd" / "analyses"
        _prepare_project(tmp_path)

        plugin_name = "scaffold_e2e"
        generate_scaffold(plugin_name, tmp_path)

        import polyzymd.analyses as analyses_pkg
        from polyzymd.analyses.base import AggregateContext, Condition, ReplicateContext
        from polyzymd.analyses.discovery import clear_cache, get_analysis, list_analyses

        original_path = list(analyses_pkg.__path__)
        monkeypatch.setattr(analyses_pkg, "__path__", [*original_path, str(analyses_root)])
        clear_cache()

        try:
            discovered = list_analyses()
            assert plugin_name in discovered

            analysis_cls = get_analysis(plugin_name)
            analysis = analysis_cls()
            settings = analysis_cls.Settings()
            assert "run_replicate" not in analysis_cls.__dict__
            assert "build_runner" not in analysis_cls.__dict__
            assert analysis_cls.measurement.metric.name == "scaffold_e2e_value"

            condition = Condition(
                label="Scaffold Condition",
                config_path=tmp_path / "config.yaml",
                replicates=(1, 2),
                sim_config=FakeSimulationConfig(name="scaffold_condition"),
            )

            rep_ctx = ReplicateContext(
                condition=condition,
                replicate=1,
                sim_config=condition.sim_config,
                output_dir=tmp_path / "run_1",
                equilibration="0ns",
                recompute=True,
                settings=settings,
            )
            monkeypatch.setattr(analysis, "_trajectory_loader_factory", lambda: self.FakeLoader)
            monkeypatch.setattr(analysis, "get_trajectory_window", lambda *args: self.FakeWindow())
            replicate_result = analysis.run_replicate(rep_ctx, replicate=1)
            assert isinstance(replicate_result, dict)
            assert replicate_result["analysis"] == plugin_name
            assert replicate_result["measurement"] == "scaffold_e2e_measurement"
            assert replicate_result["metric"] == "scaffold_e2e_value"
            assert replicate_result["value"] == pytest.approx(2.0)

            agg_ctx = AggregateContext(
                condition=condition,
                replicates=(1, 2),
                output_dir=tmp_path / "aggregated",
                equilibration="0ns",
                settings=settings,
            )
            aggregated = analysis.aggregate(
                agg_ctx,
                results=[replicate_result],
            )
            assert isinstance(aggregated, dict)
            assert aggregated["mean_value"] == pytest.approx(2.0)
            assert aggregated["n_replicates"] == 1

            metrics = analysis.extract_metrics(aggregated)
            assert "scaffold_e2e_value" in metrics
        finally:
            clear_cache()
            sys.modules.pop(f"polyzymd.analyses.{plugin_name}", None)
            plugin_path = analyses_root / f"{plugin_name}.py"
            if plugin_path.exists():
                plugin_path.unlink()


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
        assert validate_class_name("") is not None


class TestGenerateScaffoldClassNameValidation:
    """Ensure generate_scaffold rejects invalid class names."""

    def test_invalid_class_name_raises(self, tmp_path: Path):
        _prepare_project(tmp_path)

        with pytest.raises(ValueError, match="identifier"):
            generate_scaffold("density", tmp_path, class_name="2Bad")

    def test_lowercase_class_name_raises(self, tmp_path: Path):
        _prepare_project(tmp_path)

        with pytest.raises(ValueError, match="uppercase"):
            generate_scaffold("density", tmp_path, class_name="massDensity")


class TestNewAnalysisCLI:
    """Test the ``new-analysis`` Click command via CliRunner."""

    @pytest.fixture(autouse=True)
    def _setup_project_dir(self, tmp_path: Path):
        """Create a minimal project structure and store the root."""
        _prepare_project(tmp_path)
        self.root = tmp_path

    @pytest.fixture
    def runner(self):
        return CliRunner()

    @pytest.fixture
    def cli(self):
        from polyzymd.cli.main import new_analysis

        return new_analysis

    def test_success(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["solvent_shell", "--project-root", str(self.root)])
        assert result.exit_code == 0, result.output
        assert "scaffolded successfully" in result.output
        assert (self.root / "src" / "polyzymd" / "analyses" / "solvent_shell.py").exists()
        assert not (self.root / "src" / "polyzymd" / "analyses" / "solvent_shell").exists()

    def test_default_style_is_measurement(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["solvent_shell", "--project-root", str(self.root)])
        assert result.exit_code == 0, result.output
        plugin_path = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell.py"
        text = plugin_path.read_text(encoding="utf-8")
        assert "ScalarMeasurementAnalysis" in text
        assert "class SolventShellMeasurement" in text
        assert "SolventShellReplicateRunner" not in text

    def test_style_pydantic(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["solvent_shell", "--style", "pydantic", "--project-root", str(self.root)],
        )
        assert result.exit_code != 0
        assert "Advanced analysis scaffolds are temporarily unavailable" in result.output

    def test_style_dict_explicit(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["solvent_shell", "--style", "dict", "--project-root", str(self.root)],
        )
        assert result.exit_code != 0
        assert "Advanced analysis scaffolds are temporarily unavailable" in result.output

    def test_advanced_flag(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["solvent_shell", "--advanced", "--project-root", str(self.root)],
        )
        assert result.exit_code != 0
        assert "Advanced analysis scaffolds are temporarily unavailable" in result.output

    def test_invalid_style_rejected(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["solvent_shell", "--style", "invalid", "--project-root", str(self.root)],
        )
        assert result.exit_code != 0
        assert "Invalid value" in result.output or "invalid" in result.output.lower()

    def test_dry_run(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["solvent_shell", "--project-root", str(self.root), "--dry-run"],
        )
        assert result.exit_code == 0, result.output
        assert result.output.count("Would create") == 2
        assert "src/polyzymd/analyses/solvent_shell.py" in result.output
        plugin_path = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell.py"
        assert not plugin_path.exists()

    def test_invalid_name(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["BadName", "--project-root", str(self.root)])
        assert result.exit_code != 0
        assert "snake_case" in result.output

    def test_existing_plugin_rejected(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["rmsf", "--project-root", str(self.root)])
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
        plugin_path = self.root / "src" / "polyzymd" / "analyses" / "density.py"
        assert "class MassDensityAnalysis(ScalarMeasurementAnalysis):" in plugin_path.read_text(
            encoding="utf-8"
        )

    def test_force_overwrites(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["solvent_shell", "--project-root", str(self.root)])
        assert result.exit_code == 0

        result = runner.invoke(cli, ["solvent_shell", "--project-root", str(self.root)])
        assert result.exit_code != 0
        assert "already exists" in result.output

        result = runner.invoke(
            cli,
            ["solvent_shell", "--project-root", str(self.root), "--force"],
        )
        assert result.exit_code == 0, result.output

    @pytest.mark.parametrize(
        ("plugin_name", "style_args", "stale_path_parts"),
        [
            ("scaffold_force_measurement", [], ("scaffold_force_measurement.py",)),
        ],
    )
    def test_force_overwrites_discoverable_scaffold_after_fresh_discovery(
        self,
        runner: CliRunner,
        cli,
        monkeypatch: pytest.MonkeyPatch,
        plugin_name: str,
        style_args: list[str],
        stale_path_parts: tuple[str, ...],
    ):
        analyses_root = self.root / "src" / "polyzymd" / "analyses"
        result = runner.invoke(cli, [plugin_name, *style_args, "--project-root", str(self.root)])
        assert result.exit_code == 0, result.output

        plugin_path = analyses_root / Path(*stale_path_parts)
        plugin_path.write_text(
            plugin_path.read_text(encoding="utf-8") + "\n# stale local edit\n",
            encoding="utf-8",
        )

        import polyzymd.analyses as analyses_pkg
        from polyzymd.analyses.discovery import clear_cache, list_analyses

        original_path = list(analyses_pkg.__path__)
        monkeypatch.setattr(analyses_pkg, "__path__", [*original_path, str(analyses_root)])
        clear_cache()

        try:
            assert plugin_name in list_analyses()

            result = runner.invoke(
                cli,
                [plugin_name, *style_args, "--project-root", str(self.root), "--force"],
            )
            assert result.exit_code == 0, result.output
            assert "stale local edit" not in plugin_path.read_text(encoding="utf-8")
        finally:
            clear_cache()
            sys.modules.pop(f"polyzymd.analyses.{plugin_name}", None)

    def test_force_rejects_registered_builtin_name(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["rmsf", "--project-root", str(self.root), "--force"],
        )
        assert result.exit_code != 0
        assert "registered analysis plugin" in result.output

    def test_force_rejects_registered_external_name(
        self,
        runner: CliRunner,
        cli,
        monkeypatch: pytest.MonkeyPatch,
    ):
        plugin_name = "external_force"
        external_root = self.root / "external_plugins"
        external_root.mkdir()
        external_module = external_root / f"{plugin_name}.py"
        external_module.write_text(
            "from typing import ClassVar\n\n"
            "from pydantic import BaseModel\n\n"
            "from polyzymd.analyses.base import Analysis\n\n\n"
            "class ExternalForceSettings(BaseModel):\n"
            "    pass\n\n\n"
            "class ExternalForceAnalysis(Analysis):\n"
            f"    name: ClassVar[str] = {plugin_name!r}\n"
            "    Settings: ClassVar[type[BaseModel]] = ExternalForceSettings\n"
            "    has_compute_stage: ClassVar[bool] = False\n"
            "    has_aggregate_stage: ClassVar[bool] = False\n",
            encoding="utf-8",
        )

        import polyzymd.analyses as analyses_pkg
        from polyzymd.analyses.discovery import clear_cache, list_analyses

        original_path = list(analyses_pkg.__path__)
        monkeypatch.setattr(analyses_pkg, "__path__", [*original_path, str(external_root)])
        clear_cache()

        try:
            assert plugin_name in list_analyses()

            result = runner.invoke(
                cli,
                [plugin_name, "--project-root", str(self.root), "--force"],
            )

            assert result.exit_code != 0
            assert "registered analysis plugin" in result.output
            assert not (self.root / "src" / "polyzymd" / "analyses" / f"{plugin_name}.py").exists()
        finally:
            clear_cache()
            sys.modules.pop(f"polyzymd.analyses.{plugin_name}", None)

    def test_help_shows_correct_test_path(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "tests/analyses/plugins/test_<NAME>.py" in result.output

    def test_help_shows_measurement_and_advanced_layouts(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "src/polyzymd/analyses/<NAME>.py" in result.output
        assert "_runner.py" not in result.output
        assert "Advanced MDAnalysis-native package scaffolds" in result.output

    def test_success_message_shows_correct_test_path(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["solvent_shell", "--project-root", str(self.root)])
        assert result.exit_code == 0, result.output
        assert "tests/analyses/plugins/test_solvent_shell.py" in result.output
