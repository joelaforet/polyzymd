"""The Analysis class and the context objects the runner passes it."""

from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar

import pytest
from pydantic import BaseModel, ValidationError

from polyzymd.analyses.base import ComparisonContext, Condition, PlotContext, SlurmResourceHint
from polyzymd.analyses.contract import Observable, contract_analysis


class ToySettings(BaseModel):
    threshold: float = 1.0


class _Toy:
    name: ClassVar[str] = "toy"
    Settings: ClassVar[type[BaseModel]] = ToySettings
    references: ClassVar[tuple[str, ...]] = ()

    def compute(self, universe: Any, frames: Any, settings: ToySettings) -> list[Observable]:
        return [Observable(name="value", kind="mean_of_timeseries", unit="A", values=[1.0])]


ToyAnalysis = contract_analysis(_Toy)


@pytest.fixture
def toy_analysis():
    return ToyAnalysis()


@pytest.fixture
def toy_condition(tmp_path):
    """Build a Condition with a stand-in sim_config."""
    return Condition(
        label="Test Condition",
        config_path=tmp_path / "config.yaml",
        replicates=(1, 2, 3),
        sim_config=object(),
    )


class TestSlurmResourceHint:
    def test_default_slurm_resource_hint_is_none(self, toy_analysis) -> None:
        """An analysis that states no SLURM needs has no hint."""
        assert toy_analysis.slurm_resource_hint is None

    def test_plugin_hint_is_copied_onto_the_analysis(self) -> None:
        """contract_analysis copies a plugin's SLURM hint, which the CLI reads."""

        class _Hungry(_Toy):
            name: ClassVar[str] = "hungry"
            slurm_resource_hint = SlurmResourceHint(mem="16G", time="04:00:00", cpus_per_task=4)

        hint = contract_analysis(_Hungry)().slurm_resource_hint
        assert (hint.mem, hint.time, hint.cpus_per_task) == ("16G", "04:00:00", 4)

    def test_slurm_resource_hint_model_validation(self) -> None:
        """SlurmResourceHint should validate declared field types."""
        with pytest.raises(ValidationError):
            SlurmResourceHint(cpus_per_task="four")


class TestContextObjects:
    """Test context dataclass construction and properties."""

    def test_condition_creation(self, toy_condition):
        assert toy_condition.label == "Test Condition"
        assert toy_condition.replicates == (1, 2, 3)

    def test_comparison_context_effective_control(self, toy_condition):
        cond2 = Condition(
            label="Control",
            config_path=Path("/tmp/ctrl.yaml"),
            replicates=(1, 2),
            sim_config=object(),
        )
        ctx = ComparisonContext(
            name="Test Project",
            conditions=[toy_condition, cond2],
            excluded_conditions=[],
            control_label="Control",
            analysis_dirs={},
            results_dir=Path("/tmp/results"),
            equilibration="10ns",
            settings=ToySettings(),
            recompute=False,
            result_path=Path("/tmp/results/result.json"),
        )
        assert ctx.effective_control == "Control"

    def test_comparison_context_excluded_control(self, toy_condition):
        """If control was excluded, effective_control returns None."""
        ctx = ComparisonContext(
            name="Test",
            conditions=[toy_condition],
            excluded_conditions=[],
            control_label="Missing Control",
            analysis_dirs={},
            results_dir=Path("/tmp"),
            equilibration="0ns",
            settings=ToySettings(),
            recompute=False,
            result_path=Path("/tmp/result.json"),
        )
        assert ctx.effective_control is None


# ============================================================================
# Tests: MetricValue and extract_metrics integration
# ============================================================================


# ============================================================================
# Tests: resolve_output_dir
# ============================================================================


class TestResolveOutputDir:
    """Test path resolution utilities."""

    def test_resolve_output_dir(self, toy_analysis, tmp_path):
        result = toy_analysis.resolve_output_dir(tmp_path / "analysis", "100% SBMA")
        assert "100pct_SBMA" in str(result)
        assert result.name == "toy"


def test_plot_context_default_plot_settings() -> None:
    """PlotContext without explicit plot_settings should get a real PlotSettings."""
    from polyzymd.config.comparison import PlotSettings

    ctx = PlotContext(
        conditions=[],
        analysis_dirs={},
        results_dir=Path("/fake"),
        output_dir=Path("/fake"),
        settings=ToySettings(),
    )
    assert isinstance(ctx.plot_settings, PlotSettings)


def test_plot_context_materializes_when_none_is_explicitly_passed() -> None:
    """PlotContext should materialize PlotSettings when None is passed explicitly."""
    from polyzymd.config.comparison import PlotSettings

    ctx = PlotContext(
        conditions=[],
        analysis_dirs={},
        results_dir=Path("/fake"),
        output_dir=Path("/fake"),
        settings=ToySettings(),
        plot_settings=None,
    )

    assert isinstance(ctx.plot_settings, PlotSettings)


def test_plot_context_keeps_explicit_plot_settings_instance() -> None:
    """PlotContext should keep a valid PlotSettings instance as-is."""
    from polyzymd.config.comparison import PlotSettings

    plot_settings = PlotSettings()
    ctx = PlotContext(
        conditions=[],
        analysis_dirs={},
        results_dir=Path("/fake"),
        output_dir=Path("/fake"),
        settings=ToySettings(),
        plot_settings=plot_settings,
    )

    assert ctx.plot_settings is plot_settings


@pytest.mark.parametrize("invalid_value", [False, object()])
def test_plot_context_rejects_invalid_plot_settings_type(invalid_value: object) -> None:
    """PlotContext should raise TypeError for non-PlotSettings values."""
    with pytest.raises(TypeError, match="plot_settings must be a PlotSettings instance"):
        PlotContext(
            conditions=[],
            analysis_dirs={},
            results_dir=Path("/fake"),
            output_dir=Path("/fake"),
            settings=ToySettings(),
            plot_settings=invalid_value,
        )
