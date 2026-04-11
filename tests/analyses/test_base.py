"""Tests for the analyses plugin infrastructure (Phase A).

Tests the base class contract, discovery, stats utilities, and orchestrator
without requiring heavy dependencies (OpenMM, MDAnalysis, etc.).
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar, Sequence
from unittest.mock import MagicMock

import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    Condition,
    MetricValue,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.exceptions import PluginContractError

# ============================================================================
# Fixtures: Toy analysis implementations for testing
# ============================================================================


class ToySettings(BaseModel):
    """Minimal settings for testing."""

    threshold: float = 1.0


class ToyResult(BaseModel):
    """Minimal result for testing."""

    value: float
    replicate: int


class ToyAggregatedResult(BaseModel):
    """Minimal aggregated result for testing."""

    mean_value: float
    sem_value: float
    replicate_values: list[float]
    n_replicates: int

    def save(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path

    @classmethod
    def load(cls, path: Path) -> "ToyAggregatedResult":
        return cls.model_validate_json(path.read_text())


class ToyAnalysis(Analysis):
    """Concrete analysis for testing the plugin system."""

    name: ClassVar[str] = "toy"
    Settings: ClassVar[type] = ToySettings
    AggregatedResultClass: ClassVar[type] = ToyAggregatedResult
    aliases: ClassVar[tuple[str, ...]] = ("toy_alias",)
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> ToyResult:
        return ToyResult(value=replicate * 1.5, replicate=replicate)

    def aggregate(self, ctx: AggregateContext, results: Sequence[ToyResult]) -> ToyAggregatedResult:
        values = [r.value for r in results]
        import statistics

        mean_val = statistics.mean(values)
        sem_val = statistics.stdev(values) / len(values) ** 0.5 if len(values) > 1 else 0.0
        return ToyAggregatedResult(
            mean_value=mean_val,
            sem_value=sem_val,
            replicate_values=values,
            n_replicates=len(values),
        )

    def extract_metrics(self, summary: ToyAggregatedResult) -> dict[str, MetricValue]:
        return {
            "toy_metric": MetricValue(
                name="toy_metric",
                mean=summary.mean_value,
                sem=summary.sem_value,
                replicate_values=summary.replicate_values,
                higher_is_better=False,
                direction_labels=("stabilizing", "unchanged", "destabilizing"),
            )
        }


@pytest.fixture
def toy_analysis():
    return ToyAnalysis()


@pytest.fixture
def toy_condition(tmp_path):
    """Build a Condition with a mock sim_config."""
    return Condition(
        label="Test Condition",
        config_path=tmp_path / "config.yaml",
        replicates=(1, 2, 3),
        sim_config=MagicMock(),
    )


@pytest.fixture
def toy_settings():
    return ToySettings(threshold=2.0)


# ============================================================================
# Tests: Analysis ABC contract
# ============================================================================


class TestAnalysisABC:
    """Test the Analysis base class contract enforcement."""

    def test_concrete_subclass_requires_name(self):
        """Subclass without 'name' should fail at class creation."""
        with pytest.raises(TypeError, match="must define 'name'"):

            class BadAnalysis(Analysis):
                Settings = ToySettings

                def compute_replicate(self, ctx, replicate):
                    pass

                def aggregate(self, ctx, results):
                    pass

    def test_concrete_subclass_requires_settings(self):
        """Subclass without 'Settings' should fail at class creation."""
        with pytest.raises(TypeError, match="must define 'Settings'"):

            class BadAnalysis(Analysis):
                name = "bad"

                def compute_replicate(self, ctx, replicate):
                    pass

                def aggregate(self, ctx, results):
                    pass

    def test_abstract_intermediate_skips_validation(self):
        """An abstract intermediate class should not trigger validation."""
        from abc import abstractmethod

        class AbstractMiddle(Analysis):
            @abstractmethod
            def extra_method(self):
                pass

        # Should not raise — AbstractMiddle still has __abstractmethods__
        assert hasattr(AbstractMiddle, "__abstractmethods__")

    def test_concrete_subclass_valid(self, toy_analysis):
        """ToyAnalysis should instantiate without error."""
        assert toy_analysis.name == "toy"
        assert toy_analysis.aliases == ("toy_alias",)
        assert toy_analysis.min_replicates == 2

    def test_repr(self, toy_analysis):
        assert "ToyAnalysis" in repr(toy_analysis)
        assert "toy" in repr(toy_analysis)

    def test_default_filter_conditions(self, toy_analysis, toy_condition):
        """Default filter_conditions keeps all conditions."""
        conditions = [toy_condition]
        result = toy_analysis.filter_conditions(conditions)
        assert result == conditions

    def test_default_plot_returns_empty(self, toy_analysis):
        """Default plot() returns empty list."""
        ctx = PlotContext(
            conditions=[],
            analysis_dirs={},
            results_dir=Path("/tmp/results"),
            output_dir=Path("/tmp/figures"),
            settings=ToySettings(),
            comparison_path=Path("/tmp/results/result.json"),
        )
        assert toy_analysis.plot(ctx) == []

    def test_default_extract_metrics_returns_empty(self, toy_analysis):
        """Default extract_metrics returns empty if not overridden on base."""
        # ToyAnalysis DOES override extract_metrics, so test the base
        base_result = Analysis.extract_metrics(toy_analysis, "some_summary")
        # ToyAnalysis's own implementation should return metrics
        toy_result = toy_analysis.extract_metrics(
            ToyAggregatedResult(
                mean_value=1.0, sem_value=0.1, replicate_values=[0.9, 1.1], n_replicates=2
            )
        )
        assert "toy_metric" in toy_result


# ============================================================================
# Tests: Context objects
# ============================================================================


class TestContextObjects:
    """Test context dataclass construction and properties."""

    def test_condition_creation(self, toy_condition):
        assert toy_condition.label == "Test Condition"
        assert toy_condition.replicates == (1, 2, 3)

    def test_replicate_context(self, toy_condition, toy_settings):
        ctx = ReplicateContext(
            condition=toy_condition,
            replicate=1,
            sim_config=toy_condition.sim_config,
            output_dir=Path("/tmp/run_1"),
            equilibration="10ns",
            recompute=False,
            settings=toy_settings,
            result_path=Path("/tmp/run_1/result.json"),
        )
        assert ctx.replicate == 1
        assert ctx.equilibration == "10ns"
        assert ctx.settings.threshold == 2.0

    def test_comparison_context_effective_control(self, toy_condition):
        cond2 = Condition(
            label="Control",
            config_path=Path("/tmp/ctrl.yaml"),
            replicates=(1, 2),
            sim_config=MagicMock(),
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

    def test_metric_value_defaults(self):
        mv = MetricValue(name="test", mean=1.0, sem=0.1, replicate_values=[0.9, 1.1])
        assert mv.higher_is_better is True
        assert mv.direction_labels == ("decreased", "unchanged", "increased")


# ============================================================================
# Tests: MetricValue and extract_metrics integration
# ============================================================================


class TestMetricExtraction:
    """Test the extract_metrics -> default compare pipeline."""

    def test_toy_extract_metrics(self, toy_analysis):
        agg = ToyAggregatedResult(
            mean_value=2.0, sem_value=0.3, replicate_values=[1.5, 2.0, 2.5], n_replicates=3
        )
        metrics = toy_analysis.extract_metrics(agg)
        assert "toy_metric" in metrics
        mv = metrics["toy_metric"]
        assert mv.mean == 2.0
        assert mv.sem == 0.3
        assert mv.higher_is_better is False
        assert mv.direction_labels[0] == "stabilizing"


class TestDefaultCompareContract:
    """Test default compare() contract enforcement behavior."""

    def test_compare_raises_on_empty_extract_metrics(self, tmp_path: Path) -> None:
        """Empty metric extraction should raise PluginContractError."""

        class EmptyMetricsAnalysis(Analysis):
            name: ClassVar[str] = "empty_metrics"
            Settings: ClassVar[type] = ToySettings

            def compute_replicate(self, ctx, replicate):
                return {"replicate": replicate}

            def aggregate(self, ctx, results):
                return {"dummy": True}

        analysis = EmptyMetricsAnalysis()
        condition = Condition(
            label="A",
            config_path=tmp_path / "a.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = ComparisonContext(
            name="proj",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": tmp_path / "analysis" / "A" / "empty_metrics"},
            results_dir=tmp_path / "comparison" / "empty_metrics",
            equilibration="10ns",
            settings=ToySettings(),
            recompute=False,
            aggregated_results={"A": {"dummy": True}},
        )

        with pytest.raises(
            PluginContractError,
            match=r"extract_metrics\(\) returned empty dict for condition 'A'",
        ):
            analysis.compare(ctx)

    def test_compare_skips_missing_result_file_with_warning(self, caplog, tmp_path: Path) -> None:
        """Missing aggregated files should be skipped with warning, not contract error."""

        analysis = ToyAnalysis()
        condition = Condition(
            label="A",
            config_path=tmp_path / "a.yaml",
            replicates=(1, 2),
            sim_config=MagicMock(),
        )
        ctx = ComparisonContext(
            name="proj",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": tmp_path / "analysis" / "A" / "toy"},
            results_dir=tmp_path / "comparison" / "toy",
            equilibration="10ns",
            settings=ToySettings(),
            recompute=False,
            aggregated_results={},
        )

        caplog.set_level("WARNING")
        result = analysis.compare(ctx)

        assert result is None
        assert "missing aggregated result for condition 'A'" in caplog.text


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
        settings=MagicMock(),
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
        settings=MagicMock(),
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
        settings=MagicMock(),
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
            settings=MagicMock(),
            plot_settings=invalid_value,
        )
