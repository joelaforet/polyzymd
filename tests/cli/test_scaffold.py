"""Tests for the ``polyzymd new-analysis`` scaffold command."""

from __future__ import annotations

import shutil
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
        for reserved in ("base", "discovery", "orchestrator", "shared", "stats"):
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
    """Scaffold file generation using the default dict style."""

    def test_creates_three_files(self, tmp_path: Path):
        _prepare_project(tmp_path)

        created = generate_scaffold("solvent_shell", tmp_path)

        assert len(created) == 3
        for path in created:
            assert path.exists(), f"{path} was not created"

    def test_file_paths(self, tmp_path: Path):
        _prepare_project(tmp_path)

        created = generate_scaffold("solvent_shell", tmp_path)
        names = {path.name for path in created}

        assert names == {"__init__.py", "_runner.py", "test_solvent_shell.py"}

    def test_dict_style_facade_and_runner_content(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        plugin_dir = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell"
        init_text = (plugin_dir / "__init__.py").read_text(encoding="utf-8")
        runner_text = (plugin_dir / "_runner.py").read_text(encoding="utf-8")

        assert "class SolventShellAnalysis(Analysis):" in init_text
        assert "class SolventShellReplicateRunner:" not in init_text
        assert "from ._runner import SolventShellReplicateRunner" in init_text
        assert 'name: ClassVar[str] = "solvent_shell"' in init_text
        assert "class Settings(BaseModel):" in init_text
        assert "def build_runner(" in init_text
        assert "def summarize_replicate(" in init_text
        assert "def run_replicate(" not in init_text
        assert "def aggregate(" in init_text
        assert "def extract_metrics(" in init_text
        assert "def plot(self, ctx: PlotContext)" in init_text
        assert "_build_plot_data" in init_text
        assert "class SolventShellReplicateRunner:" in runner_text
        assert "def run(" in runner_text

    def test_dict_style_uses_plain_dicts(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        plugin_dir = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell"
        init_text = (plugin_dir / "__init__.py").read_text(encoding="utf-8")

        assert not (plugin_dir / "_results.py").exists()
        assert "class SolventShellReplicateResult" not in init_text
        assert "class SolventShellAggregatedResult" not in init_text
        assert "AggregatedResultClass: ClassVar" not in init_text
        assert "-> dict[str, Any]:" in init_text

    def test_test_file_content(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        test_path = tmp_path / "tests" / "analyses" / "plugins" / "test_solvent_shell.py"
        text = test_path.read_text(encoding="utf-8")

        assert "class TestDiscovery:" in text
        assert "class TestRunnerBackedReplicate:" in text
        assert "class TestAggregate:" in text
        assert "class TestExtractMetrics:" in text
        assert "class TestPlot:" in text
        assert "test_run_replicate_uses_base_runner_dispatch" in text
        assert 'runner.__class__.__name__ == "SolventShellReplicateRunner"' in text

    def test_custom_class_name(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("density", tmp_path, class_name="MassDensity")
        plugin_dir = tmp_path / "src" / "polyzymd" / "analyses" / "density"

        assert "class MassDensityAnalysis(Analysis):" in (plugin_dir / "__init__.py").read_text(
            encoding="utf-8"
        )
        assert "class MassDensityReplicateRunner:" in (plugin_dir / "_runner.py").read_text(
            encoding="utf-8"
        )

    def test_refuses_overwrite_without_force(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        with pytest.raises(FileExistsError, match="already exists"):
            generate_scaffold("solvent_shell", tmp_path)

    def test_force_overwrites(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        created = generate_scaffold("solvent_shell", tmp_path, force=True)
        assert len(created) == 3

    def test_dry_run_creates_no_files(self, tmp_path: Path):
        _prepare_project(tmp_path)

        created = generate_scaffold("solvent_shell", tmp_path, dry_run=True)
        assert len(created) == 3
        for path in created:
            assert not path.exists(), f"{path} should not exist in dry-run mode"

    def test_default_style_is_dict(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        init_path = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        text = init_path.read_text(encoding="utf-8")

        assert "Uses plain dicts for result containers" in text
        assert "AggregatedResultClass: ClassVar" not in text


class TestGenerateScaffoldPydantic:
    """Scaffold file generation using the pydantic style."""

    def test_creates_four_files(self, tmp_path: Path):
        _prepare_project(tmp_path)

        created = generate_scaffold("solvent_shell", tmp_path, style="pydantic")

        assert len(created) == 4
        for path in created:
            assert path.exists(), f"{path} was not created"

    def test_file_paths(self, tmp_path: Path):
        _prepare_project(tmp_path)

        created = generate_scaffold("solvent_shell", tmp_path, style="pydantic")
        names = {path.name for path in created}

        assert names == {"__init__.py", "_runner.py", "_results.py", "test_solvent_shell.py"}

    def test_plugin_facade_content(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path, style="pydantic")
        plugin_dir = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell"
        init_text = (plugin_dir / "__init__.py").read_text(encoding="utf-8")
        runner_text = (plugin_dir / "_runner.py").read_text(encoding="utf-8")

        assert "class SolventShellAnalysis(Analysis):" in init_text
        assert 'name: ClassVar[str] = "solvent_shell"' in init_text
        assert "ReplicateResultClass: ClassVar[type] = SolventShellReplicateResult" in init_text
        assert "AggregatedResultClass: ClassVar[type] = SolventShellAggregatedResult" in init_text
        assert "class SolventShellReplicateRunner:" not in init_text
        assert "from ._runner import SolventShellReplicateRunner" in init_text
        assert "from ._results import (" in init_text
        assert "SolventShellAggregatedResult," in init_text
        assert "def build_runner(" in init_text
        assert "def summarize_replicate(" in init_text
        assert "def run_replicate(" not in init_text
        assert "def aggregate(" in init_text
        assert "def extract_metrics(" in init_text
        assert "def plot(self, ctx: PlotContext)" in init_text
        assert "class SolventShellReplicateRunner:" in runner_text

    def test_result_models_live_in_results_module(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path, style="pydantic")
        plugin_dir = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell"
        init_text = (plugin_dir / "__init__.py").read_text(encoding="utf-8")
        results_text = (plugin_dir / "_results.py").read_text(encoding="utf-8")

        assert "class SolventShellReplicateResult(BaseModel):" not in init_text
        assert "class SolventShellAggregatedResult(BaseModel):" not in init_text
        assert "class SolventShellReplicateResult(BaseModel):" in results_text
        assert "class SolventShellAggregatedResult(BaseModel):" in results_text
        assert "replicates: list[int] = Field(default_factory=list)" in results_text
        assert "settings_fingerprint: str | None = None" in results_text

    def test_test_file_content(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path, style="pydantic")
        test_path = tmp_path / "tests" / "analyses" / "plugins" / "test_solvent_shell.py"
        text = test_path.read_text(encoding="utf-8")

        assert "class TestDiscovery:" in text
        assert "class TestRunnerBackedReplicate:" in text
        assert "class TestAggregate:" in text
        assert "class TestExtractMetrics:" in text
        assert "class TestPlot:" in text
        assert "test_run_replicate_uses_base_runner_dispatch" in text
        assert "SolventShellReplicateResult" in text
        assert "SolventShellAggregatedResult" in text
        assert 'runner.__class__.__name__ == "SolventShellReplicateRunner"' in text

    def test_custom_class_name(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("density", tmp_path, class_name="MassDensity", style="pydantic")
        plugin_dir = tmp_path / "src" / "polyzymd" / "analyses" / "density"

        assert "class MassDensityAnalysis(Analysis):" in (plugin_dir / "__init__.py").read_text(
            encoding="utf-8"
        )
        assert "class MassDensityReplicateRunner:" in (plugin_dir / "_runner.py").read_text(
            encoding="utf-8"
        )
        assert "class MassDensityReplicateResult(BaseModel):" in (
            plugin_dir / "_results.py"
        ).read_text(encoding="utf-8")
        assert "class MassDensityAggregatedResult(BaseModel):" in (
            plugin_dir / "_results.py"
        ).read_text(encoding="utf-8")

    def test_force_overwrites(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path, style="pydantic")
        created = generate_scaffold("solvent_shell", tmp_path, style="pydantic", force=True)
        assert len(created) == 4


class TestStyleValidation:
    """Style parameter validation."""

    def test_valid_styles_constant(self):
        assert "dict" in VALID_STYLES
        assert "pydantic" in VALID_STYLES

    def test_invalid_style_raises(self, tmp_path: Path):
        _prepare_project(tmp_path)

        with pytest.raises(ValueError, match="Invalid style"):
            generate_scaffold("solvent_shell", tmp_path, style="bad_style")


class TestTemplateResources:
    """Jinja package-resource template behavior."""

    def test_template_resources_are_present(self):
        template_root = files("polyzymd.cli._scaffold") / "templates"
        names = {path.name for path in template_root.iterdir()}

        assert {
            "plugin_init.py.jinja",
            "runner.py.jinja",
            "results.py.jinja",
            "test_plugin.py.jinja",
        }.issubset(names)

    def test_renderer_uses_strict_undefined(self):
        env = create_environment()

        with pytest.raises(UndefinedError):
            env.from_string("{{ missing_value }}").render()

    def test_templates_render_with_sample_context(self):
        dict_spec = ScaffoldSpec(name="strict_sample", class_name="StrictSample", style="dict")
        pydantic_spec = ScaffoldSpec(
            name="strict_sample",
            class_name="StrictSample",
            style="pydantic",
        )

        for template_name in ("plugin_init.py.jinja", "runner.py.jinja", "test_plugin.py.jinja"):
            rendered = render_template(template_name, dict_spec)
            assert "StrictSample" in rendered or "strict_sample" in rendered

        rendered_results = render_template("results.py.jinja", pydantic_spec)
        assert "StrictSampleReplicateResult" in rendered_results


class TestGeneratedCodeQuality:
    """Ensure generated code compiles and satisfies formatters."""

    @pytest.mark.parametrize(
        ("name", "style", "expected_count"),
        [
            ("solvent_shell", "dict", 3),
            ("solvent_shell", "pydantic", 4),
            ("scaffold_pydantic_e2e", "dict", 3),
            ("scaffold_pydantic_e2e", "pydantic", 4),
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
        """Minimal trajectory with a length for runner dispatch tests."""

        def __len__(self) -> int:
            return 5

    class FakeUniverse:
        """Minimal Universe exposing a trajectory."""

        def __init__(self) -> None:
            self.trajectory = TestGeneratedPluginEndToEnd.FakeTrajectory()

    class FakeWindow:
        """Minimal trajectory window for base runner dispatch."""

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

    def test_dict_scaffold_is_discoverable_and_comparable(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        analyses_root = tmp_path / "src" / "polyzymd" / "analyses"
        _prepare_project(tmp_path)

        plugin_name = "scaffold_e2e"
        generate_scaffold(plugin_name, tmp_path, style="dict")

        import polyzymd.analyses as analyses_pkg
        from polyzymd.analyses.base import (
            AggregateContext,
            AggregateValidationError,
            ComparisonContext,
            Condition,
            PlotContext,
            ReplicateContext,
        )
        from polyzymd.analyses.discovery import clear_cache, get_analysis, list_analyses
        from polyzymd.config.comparison import PlotSettings

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
            assert replicate_result["value"] == pytest.approx(2.0)
            assert replicate_result["n_frames"] == 2

            agg_ctx = AggregateContext(
                condition=condition,
                replicates=(1, 2),
                output_dir=tmp_path / "aggregated",
                equilibration="0ns",
                settings=settings,
            )
            aggregated = analysis.aggregate(
                agg_ctx,
                results=[
                    {"value": 1.0, "replicate": 1, "n_frames": 4},
                    {"value": 3.0, "replicate": 2, "n_frames": 4},
                ],
            )
            assert isinstance(aggregated, dict)
            assert aggregated["mean_value"] == pytest.approx(2.0)
            assert aggregated["n_replicates"] == 2

            metrics = analysis.extract_metrics(aggregated)
            assert "value" in metrics

            compare_ctx = ComparisonContext(
                name="scaffold-e2e",
                conditions=[condition],
                excluded_conditions=[],
                control_label=None,
                analysis_dirs={condition.label: tmp_path / "analysis_dir"},
                results_dir=tmp_path / "results",
                equilibration="0ns",
                settings=settings,
                recompute=True,
                aggregated_results={condition.label: aggregated},
            )
            comparison = analysis.compare(compare_ctx)
            assert comparison is not None
            assert comparison.analysis_type == plugin_name
            assert len(comparison.conditions) == 1
            assert comparison.conditions[0].label == condition.label

            stale_analysis_dir = tmp_path / "stale_analysis_dir"
            stale_agg_dir = stale_analysis_dir / "aggregated"
            stale_agg_dir.mkdir(parents=True)
            (stale_agg_dir / "result.json").write_text(
                '{"mean_value": 2.0, "sem_value": 0.1, "replicate_values": [1.0, 3.0], '
                '"n_replicates": 2, "replicates": [1, 2], '
                '"settings_fingerprint": "deadbeef"}',
                encoding="utf-8",
            )
            plot_ctx = PlotContext(
                conditions=[condition],
                analysis_dirs={condition.label: stale_analysis_dir},
                results_dir=tmp_path / "comparison" / plugin_name,
                output_dir=tmp_path / "figures" / plugin_name,
                settings=settings,
                plot_settings=PlotSettings(),
                equilibration="0ns",
            )
            with pytest.raises(AggregateValidationError, match="settings fingerprint mismatch"):
                analysis.plot(plot_ctx)
        finally:
            clear_cache()
            for module_name in [
                f"polyzymd.analyses.{plugin_name}",
                f"polyzymd.analyses.{plugin_name}._runner",
            ]:
                sys.modules.pop(module_name, None)
            plugin_dir = analyses_root / plugin_name
            if plugin_dir.exists():
                shutil.rmtree(plugin_dir)

    def test_pydantic_result_models_are_importable_from_facade(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        analyses_root = tmp_path / "src" / "polyzymd" / "analyses"
        _prepare_project(tmp_path)

        plugin_name = "scaffold_pydantic_e2e"
        generate_scaffold(plugin_name, tmp_path, style="pydantic")

        import polyzymd.analyses as analyses_pkg
        from polyzymd.analyses.discovery import clear_cache, get_analysis

        original_path = list(analyses_pkg.__path__)
        monkeypatch.setattr(analyses_pkg, "__path__", [*original_path, str(analyses_root)])
        clear_cache()

        try:
            analysis_cls = get_analysis(plugin_name)
            module = sys.modules[analysis_cls.__module__]

            assert hasattr(module, "ScaffoldPydanticE2eReplicateResult")
            assert hasattr(module, "ScaffoldPydanticE2eAggregatedResult")
            assert analysis_cls.ReplicateResultClass is module.ScaffoldPydanticE2eReplicateResult
            assert analysis_cls.AggregatedResultClass is module.ScaffoldPydanticE2eAggregatedResult
        finally:
            clear_cache()
            for module_name in [
                f"polyzymd.analyses.{plugin_name}",
                f"polyzymd.analyses.{plugin_name}._results",
                f"polyzymd.analyses.{plugin_name}._runner",
            ]:
                sys.modules.pop(module_name, None)
            plugin_dir = analyses_root / plugin_name
            if plugin_dir.exists():
                shutil.rmtree(plugin_dir)


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
        assert (
            self.root / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        ).exists()
        assert (
            self.root / "src" / "polyzymd" / "analyses" / "solvent_shell" / "_runner.py"
        ).exists()

    def test_default_style_is_dict(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["solvent_shell", "--project-root", str(self.root)])
        assert result.exit_code == 0, result.output
        plugin_dir = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell"
        text = (plugin_dir / "__init__.py").read_text(encoding="utf-8")
        assert "Uses plain dicts for result containers" in text
        assert "AggregatedResultClass: ClassVar" not in text
        assert not (plugin_dir / "_results.py").exists()

    def test_style_pydantic(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["solvent_shell", "--style", "pydantic", "--project-root", str(self.root)],
        )
        assert result.exit_code == 0, result.output
        plugin_dir = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell"
        init_text = (plugin_dir / "__init__.py").read_text(encoding="utf-8")
        results_text = (plugin_dir / "_results.py").read_text(encoding="utf-8")
        assert "AggregatedResultClass" in init_text
        assert "class SolventShellReplicateResult(BaseModel):" in results_text

    def test_style_dict_explicit(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["solvent_shell", "--style", "dict", "--project-root", str(self.root)],
        )
        assert result.exit_code == 0, result.output
        init_path = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        assert "Uses plain dicts for result containers" in init_path.read_text(encoding="utf-8")

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
        assert result.output.count("Would create") == 3
        assert "src/polyzymd/analyses/solvent_shell/_runner.py" in result.output
        init_path = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell" / "__init__.py"
        assert not init_path.exists()

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
        plugin_dir = self.root / "src" / "polyzymd" / "analyses" / "density"
        assert "class MassDensityAnalysis(Analysis):" in (plugin_dir / "__init__.py").read_text(
            encoding="utf-8"
        )
        assert "class MassDensityReplicateRunner:" in (plugin_dir / "_runner.py").read_text(
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

    def test_help_shows_correct_test_path(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "tests/analyses/plugins/test_<NAME>.py" in result.output

    def test_help_shows_facade_layout(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "src/polyzymd/analyses/<NAME>/_runner.py" in result.output
        assert "src/polyzymd/analyses/<NAME>/_results.py" in result.output

    def test_success_message_shows_correct_test_path(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["solvent_shell", "--project-root", str(self.root)])
        assert result.exit_code == 0, result.output
        assert "tests/analyses/plugins/test_solvent_shell.py" in result.output
