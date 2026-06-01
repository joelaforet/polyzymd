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
    """Scaffold file generation using the default MDAnalysis-native style."""

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

    def test_simple_mda_plugin_content(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        plugin_path = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell.py"
        text = plugin_path.read_text(encoding="utf-8")

        assert "class SolventShellSettings(BaseModel):" in text
        assert "selection: str = Field(" in text
        assert "scale: float = Field(" in text
        assert "def measure_solvent_shell(" in text
        assert "class SolventShellArtifactCollector:" in text
        assert "class SolventShellAnalysis(Analysis):" in text
        assert 'name: ClassVar[str] = "solvent_shell"' in text
        assert "Settings: ClassVar[type[BaseModel]] = SolventShellSettings" in text
        assert "def build_mda_jobs(" in text
        assert "def build_mda_collector(" in text
        assert "def extract_metrics(" in text
        assert "ScalarMeasurement" not in text
        assert "def build_runner(" not in text

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
        assert "class TestMDAJobs:" in text
        assert "class TestCollector:" in text
        assert "class TestDefaultAggregation:" in text
        assert "test_build_mda_jobs_returns_function_adapter_job" in text
        assert "AggregateContext" in text
        assert "ReplicateContext" in text
        assert "ReplicateArtifact" in text

    def test_custom_class_name(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("density", tmp_path, class_name="MassDensity")
        plugin_path = tmp_path / "src" / "polyzymd" / "analyses" / "density.py"

        assert "class MassDensitySettings(BaseModel):" in plugin_path.read_text(encoding="utf-8")
        assert "class MassDensityArtifactCollector:" in plugin_path.read_text(encoding="utf-8")
        assert "class MassDensityAnalysis(Analysis):" in plugin_path.read_text(encoding="utf-8")

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

    def test_default_uses_simple_mda_scaffold(self, tmp_path: Path):
        _prepare_project(tmp_path)

        generate_scaffold("solvent_shell", tmp_path)
        plugin_path = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell.py"
        text = plugin_path.read_text(encoding="utf-8")

        assert "MDAAnalysisJob.from_function" in text
        assert "ScalarMeasurementAnalysis" not in text
        assert "SolventShellReplicateRunner" not in text

    def test_simple_style_uses_simple_mda_scaffold(self, tmp_path: Path):
        _prepare_project(tmp_path)

        created = generate_scaffold("solvent_shell", tmp_path, style="simple")
        plugin_path = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell.py"
        analyses_root = tmp_path / "src" / "polyzymd" / "analyses"
        text = plugin_path.read_text(encoding="utf-8")

        assert len(created) == 2
        assert plugin_path.exists()
        assert not (analyses_root / "solvent_shell").exists()
        assert "MDAAnalysisJob.from_function" in text
        assert "class SolventShellAnalysis(Analysis):" in text
        assert "ScalarMeasurementAnalysis" not in text
        assert "SolventShellReplicateRunner" not in text


class TestGenerateScaffoldAdvancedDict:
    """Advanced dict scaffolds generate MDAnalysis-native packages."""

    def test_style_dict_creates_package_files(self, tmp_path: Path):
        _prepare_project(tmp_path)

        created = generate_scaffold("solvent_shell", tmp_path, style="dict")

        assert len(created) == 3
        package_dir = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell"
        assert (package_dir / "__init__.py").exists()
        assert (package_dir / "_mda.py").exists()
        assert not (package_dir / "_results.py").exists()

    def test_advanced_flag_defaults_to_dict_style(self, tmp_path: Path):
        _prepare_project(tmp_path)

        created = generate_scaffold("solvent_shell", tmp_path, advanced=True)

        assert len(created) == 3
        package_dir = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell"
        assert (package_dir / "__init__.py").exists()
        assert (package_dir / "_mda.py").exists()


class TestStyleValidation:
    """Style parameter validation."""

    def test_valid_styles_constant(self):
        assert "simple" in VALID_STYLES
        assert "measurement" not in VALID_STYLES
        assert "dict" in VALID_STYLES
        assert "pydantic" not in VALID_STYLES

    def test_measurement_style_raises(self, tmp_path: Path):
        _prepare_project(tmp_path)
        legacy_style = "measure" + "ment"

        with pytest.raises(ValueError, match="Invalid style"):
            generate_scaffold("solvent_shell", tmp_path, style=legacy_style)

    def test_pydantic_style_raises(self, tmp_path: Path):
        _prepare_project(tmp_path)

        with pytest.raises(ValueError, match="Invalid style"):
            generate_scaffold("solvent_shell", tmp_path, style="pydantic")

    def test_invalid_style_raises(self, tmp_path: Path):
        _prepare_project(tmp_path)

        with pytest.raises(ValueError, match="Invalid style"):
            generate_scaffold("solvent_shell", tmp_path, style="bad_style")


class TestLayoutCollisionChecks:
    """Cross-layout collision behavior for module and package scaffolds."""

    def test_compatibility_style_refuses_existing_package_even_with_force(self, tmp_path: Path):
        _prepare_project(tmp_path)
        package_dir = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell"
        package_dir.mkdir()

        with pytest.raises(FileExistsError, match="package layout"):
            generate_scaffold("solvent_shell", tmp_path, force=True)

    def test_advanced_refuses_existing_module_even_with_force(self, tmp_path: Path):
        _prepare_project(tmp_path)
        module_path = tmp_path / "src" / "polyzymd" / "analyses" / "solvent_shell.py"
        module_path.write_text("", encoding="utf-8")

        with pytest.raises(FileExistsError, match="single-file layout"):
            generate_scaffold("solvent_shell", tmp_path, style="dict", force=True)


class TestTemplateResources:
    """Jinja package-resource template behavior."""

    def test_template_resources_are_present(self):
        template_root = files("polyzymd.cli._scaffold") / "templates"
        names = {path.name for path in template_root.iterdir()}

        assert {
            "simple_mda_plugin.py.jinja",
            "test_simple_mda_plugin.py.jinja",
            "advanced_plugin_init.py.jinja",
            "advanced_mda.py.jinja",
            "test_advanced_plugin.py.jinja",
        }.issubset(names)

    def test_renderer_uses_strict_undefined(self):
        env = create_environment()

        with pytest.raises(UndefinedError):
            env.from_string("{{ missing_value }}").render()

    def test_templates_render_with_sample_context(self):
        simple_spec = ScaffoldSpec(
            name="strict_sample",
            class_name="StrictSample",
            style="simple",
        )
        advanced_spec = ScaffoldSpec(
            name="strict_sample",
            class_name="StrictSample",
            style="dict",
        )
        for template_name, spec in (
            ("simple_mda_plugin.py.jinja", simple_spec),
            ("test_simple_mda_plugin.py.jinja", simple_spec),
            ("advanced_plugin_init.py.jinja", advanced_spec),
            ("advanced_mda.py.jinja", advanced_spec),
            ("test_advanced_plugin.py.jinja", advanced_spec),
        ):
            rendered = render_template(template_name, spec)
            assert "StrictSample" in rendered or "strict_sample" in rendered


class TestGeneratedCodeQuality:
    """Ensure generated code compiles and satisfies formatters."""

    @pytest.mark.parametrize(
        ("name", "style", "expected_count"),
        [
            ("solvent_shell", "simple", 2),
            ("scaffold_advanced_dict_e2e", "dict", 3),
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

    class FakeAtomGroup:
        """Minimal atom group with a length."""

        def __len__(self) -> int:
            return 3

    class FakeTrajectory:
        """Minimal trajectory with a length for MDA job tests."""

        def __len__(self) -> int:
            return 10

    class FakeUniverse:
        """Minimal Universe exposing selection and trajectory behavior."""

        def __init__(self) -> None:
            self.trajectory = TestGeneratedPluginEndToEnd.FakeTrajectory()
            self.selections: list[str] = []

        def select_atoms(self, selection: str):
            self.selections.append(selection)
            return TestGeneratedPluginEndToEnd.FakeAtomGroup()

    def test_simple_mda_scaffold_is_discoverable_and_aggregates_artifacts(
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
        from polyzymd.analyses.mda import (
            ArtifactStore,
            FrameSelection,
            MDACollectorContext,
            MDAUniversePolicy,
            ReplicateArtifact,
        )

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

            condition = Condition(
                label="Scaffold Condition",
                config_path=tmp_path / "config.yaml",
                replicates=(1, 2),
                sim_config=FakeSimulationConfig(name="scaffold_condition"),
            )

            frame_selection = FrameSelection(start=1, stop=5, step=2, n_frames_total=10)
            universe_policy = MDAUniversePolicy(condition_label=condition.label, replicate=1)
            universe = self.FakeUniverse()
            job_ctx = SimpleNamespace(
                universe=universe,
                frame_selection=frame_selection,
                universe_policy=universe_policy,
                settings=settings,
            )
            completed = analysis.build_mda_jobs(job_ctx)[0].run()
            assert completed.results["metrics"]["scaffold_e2e_value"] == pytest.approx(6.0)

            explicit_frame_ctx = SimpleNamespace(
                universe=self.FakeUniverse(),
                frame_selection=FrameSelection(frames=(0, 2, 4), n_frames_total=10),
                universe_policy=universe_policy,
                settings=settings,
            )
            explicit_frame_result = analysis.build_mda_jobs(explicit_frame_ctx)[0].run()
            assert explicit_frame_result.results["n_frames"] == 3
            assert explicit_frame_result.results["metrics"]["scaffold_e2e_value"] == pytest.approx(
                9.0
            )

            artifacts: list[ReplicateArtifact] = []
            for replicate, value in ((1, 2.0), (2, 4.0)):
                output_dir = tmp_path / f"run_{replicate}"
                rep_ctx = ReplicateContext(
                    condition=condition,
                    replicate=replicate,
                    sim_config=condition.sim_config,
                    output_dir=output_dir,
                    equilibration="0ns",
                    recompute=True,
                    settings=settings,
                    result_path=output_dir / "result.json",
                )
                collector_ctx = MDACollectorContext(
                    analysis_name=plugin_name,
                    replicate_context=rep_ctx,
                    frame_selection=frame_selection,
                    universe_policy=MDAUniversePolicy(
                        condition_label=condition.label,
                        replicate=replicate,
                    ),
                    artifact_store=ArtifactStore(output_dir),
                )
                artifact = ReplicateArtifact(
                    analysis_name=plugin_name,
                    condition_label=condition.label,
                    replicate=replicate,
                    payload={"metrics": {"scaffold_e2e_value": value}},
                    provenance={"frame_selection": {"start": 0, "stop": None, "step": 1}},
                    metadata={
                        "settings_fingerprint": analysis.aggregate_settings_fingerprint(settings),
                    },
                )
                artifact_from_collector = analysis.build_mda_collector(collector_ctx)(
                    collector_ctx,
                    [completed],
                )
                assert artifact_from_collector.payload["metrics"]["scaffold_e2e_value"]
                ArtifactStore(output_dir).write_replicate_result(artifact)
                artifacts.append(artifact)

            agg_ctx = AggregateContext(
                condition=condition,
                replicates=(1, 2),
                output_dir=tmp_path / "aggregated",
                equilibration="0ns",
                settings=settings,
            )
            aggregated = analysis.aggregate(
                agg_ctx,
                results=artifacts,
            )
            assert aggregated.payload["metrics"]["scaffold_e2e_value"]["mean"] == pytest.approx(3.0)

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

    def test_default_uses_simple_mda_scaffold(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["solvent_shell", "--project-root", str(self.root)])
        assert result.exit_code == 0, result.output
        plugin_path = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell.py"
        text = plugin_path.read_text(encoding="utf-8")
        assert "MDAAnalysisJob.from_function" in text
        assert "class SolventShellAnalysis(Analysis):" in text
        assert "ScalarMeasurementAnalysis" not in text
        assert "SolventShellReplicateRunner" not in text

    def test_style_pydantic_rejected(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["solvent_shell", "--style", "pydantic", "--project-root", str(self.root)],
        )
        assert result.exit_code != 0
        assert "Invalid value" in result.output

    def test_style_measurement_rejected(self, runner: CliRunner, cli):
        legacy_style = "measure" + "ment"
        result = runner.invoke(
            cli,
            ["solvent_shell", "--style", legacy_style, "--project-root", str(self.root)],
        )
        assert result.exit_code != 0
        assert "Invalid value" in result.output

    def test_style_dict_explicit(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["solvent_shell", "--style", "dict", "--project-root", str(self.root)],
        )
        assert result.exit_code == 0, result.output
        package_dir = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell"
        assert (package_dir / "__init__.py").exists()
        assert (package_dir / "_mda.py").exists()

    def test_advanced_flag(self, runner: CliRunner, cli):
        result = runner.invoke(
            cli,
            ["solvent_shell", "--advanced", "--project-root", str(self.root)],
        )
        assert result.exit_code == 0, result.output
        package_dir = self.root / "src" / "polyzymd" / "analyses" / "solvent_shell"
        assert (package_dir / "__init__.py").exists()
        assert (package_dir / "_mda.py").exists()

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
        assert "class MassDensityAnalysis(Analysis):" in plugin_path.read_text(encoding="utf-8")

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
            ("scaffold_force_simple", [], ("scaffold_force_simple.py",)),
            (
                "scaffold_force_advanced",
                ["--style", "dict"],
                ("scaffold_force_advanced", "__init__.py"),
            ),
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

    def test_help_shows_default_and_advanced_layouts(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "src/polyzymd/analyses/<NAME>.py" in result.output
        assert "src/polyzymd/analyses/<NAME>/__init__.py" in result.output
        assert "only dict canonical" in result.output
        assert "artifacts are supported" in result.output
        assert f"--style {'measure' + 'ment'}" not in result.output
        assert "_runner.py" not in result.output
        assert "Advanced package scaffolds" in result.output

    def test_success_message_shows_correct_test_path(self, runner: CliRunner, cli):
        result = runner.invoke(cli, ["solvent_shell", "--project-root", str(self.root)])
        assert result.exit_code == 0, result.output
        assert "tests/analyses/plugins/test_solvent_shell.py" in result.output
