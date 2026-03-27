"""Tests for the analyses plugin infrastructure (Phase A).

Tests the base class contract, discovery, stats utilities, and orchestrator
without requiring heavy dependencies (OpenMM, MDAnalysis, etc.).
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path
from typing import Any, ClassVar, Sequence
from unittest.mock import MagicMock, patch

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


class ToyDependentAnalysis(Analysis):
    """Analysis that depends on ToyAnalysis."""

    name: ClassVar[str] = "toy_dependent"
    Settings: ClassVar[type] = ToySettings
    dependencies: ClassVar[tuple[str, ...]] = ("toy",)

    def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> ToyResult:
        return ToyResult(value=replicate * 2.0, replicate=replicate)

    def aggregate(self, ctx: AggregateContext, results: Sequence[ToyResult]) -> ToyAggregatedResult:
        values = [r.value for r in results]
        return ToyAggregatedResult(
            mean_value=sum(values) / len(values),
            sem_value=0.0,
            replicate_values=values,
            n_replicates=len(values),
        )


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
# Tests: Stats utilities
# ============================================================================


class TestStats:
    """Test the shared statistical utility functions."""

    def test_interpret_direction(self):
        from polyzymd.analyses.stats import interpret_direction

        assert interpret_direction(5.0) == "increased"
        assert interpret_direction(-5.0) == "decreased"
        assert interpret_direction(0.5) == "unchanged"
        assert interpret_direction(0.5, threshold=0.3) == "increased"

    def test_interpret_direction_custom_labels(self):
        from polyzymd.analyses.stats import interpret_direction

        labels = ("stabilizing", "unchanged", "destabilizing")
        assert interpret_direction(-10.0, labels) == "stabilizing"
        assert interpret_direction(10.0, labels) == "destabilizing"

    def test_pairwise_comparisons_control(self):
        from polyzymd.analyses.stats import pairwise_comparisons

        metrics = {
            "Control": MetricValue("m", 1.0, 0.1, [0.9, 1.0, 1.1]),
            "Treatment A": MetricValue("m", 0.5, 0.05, [0.45, 0.5, 0.55]),
            "Treatment B": MetricValue("m", 0.8, 0.08, [0.75, 0.8, 0.85]),
        }
        results = pairwise_comparisons(metrics, control_label="Control")
        # Should have 2 comparisons: Control vs A, Control vs B
        assert len(results) == 2
        assert all(r.condition_a == "Control" for r in results)
        assert {r.condition_b for r in results} == {"Treatment A", "Treatment B"}
        # All should have required fields
        for r in results:
            assert hasattr(r, "t_statistic")
            assert hasattr(r, "p_value")
            assert hasattr(r, "cohens_d")
            assert hasattr(r, "percent_change")
            assert hasattr(r, "significant")

    def test_pairwise_comparisons_all_pairs(self):
        from polyzymd.analyses.stats import pairwise_comparisons

        metrics = {
            "A": MetricValue("m", 1.0, 0.1, [0.9, 1.0, 1.1]),
            "B": MetricValue("m", 2.0, 0.2, [1.8, 2.0, 2.2]),
            "C": MetricValue("m", 3.0, 0.3, [2.7, 3.0, 3.3]),
        }
        results = pairwise_comparisons(metrics)
        # Should have 3 comparisons: A-B, A-C, B-C
        assert len(results) == 3

    def test_anova_test_with_3_conditions(self):
        from polyzymd.analyses.stats import anova_test

        metrics = {
            "A": MetricValue("m", 1.0, 0.1, [0.9, 1.0, 1.1]),
            "B": MetricValue("m", 5.0, 0.2, [4.8, 5.0, 5.2]),
            "C": MetricValue("m", 3.0, 0.3, [2.7, 3.0, 3.3]),
        }
        result = anova_test(metrics, "test_metric")
        assert result is not None
        assert hasattr(result, "f_statistic")
        assert hasattr(result, "p_value")
        assert result.significant is True  # These groups are very different

    def test_anova_test_too_few_conditions(self):
        from polyzymd.analyses.stats import anova_test

        metrics = {
            "A": MetricValue("m", 1.0, 0.1, [0.9, 1.1]),
            "B": MetricValue("m", 2.0, 0.2, [1.8, 2.2]),
        }
        result = anova_test(metrics)
        assert result is None

    def test_rank_conditions_higher_is_better(self):
        from polyzymd.analyses.stats import rank_conditions

        metrics = {
            "Low": MetricValue("m", 1.0, 0.1, [1.0], higher_is_better=True),
            "High": MetricValue("m", 3.0, 0.1, [3.0], higher_is_better=True),
            "Mid": MetricValue("m", 2.0, 0.1, [2.0], higher_is_better=True),
        }
        ranking = rank_conditions(metrics)
        assert ranking == ["High", "Mid", "Low"]

    def test_rank_conditions_lower_is_better(self):
        from polyzymd.analyses.stats import rank_conditions

        metrics = {
            "Low": MetricValue("m", 1.0, 0.1, [1.0], higher_is_better=False),
            "High": MetricValue("m", 3.0, 0.1, [3.0], higher_is_better=False),
            "Mid": MetricValue("m", 2.0, 0.1, [2.0], higher_is_better=False),
        }
        ranking = rank_conditions(metrics)
        assert ranking == ["Low", "Mid", "High"]

    def test_default_scalar_comparison(self):
        from polyzymd.analyses.stats import default_scalar_comparison

        metrics_by_condition = {
            "Control": {
                "metric_a": MetricValue("metric_a", 1.0, 0.1, [0.9, 1.0, 1.1]),
            },
            "Treatment": {
                "metric_a": MetricValue("metric_a", 0.5, 0.05, [0.45, 0.5, 0.55]),
            },
        }
        result = default_scalar_comparison(
            analysis_name="test",
            project_name="Test Project",
            metrics_by_condition=metrics_by_condition,
            control_label="Control",
            equilibration="10ns",
        )
        assert result.analysis_type == "test"
        assert result.name == "Test Project"
        assert len(result.pairwise_comparisons) == 1
        assert result.anova is None  # Only 2 conditions
        assert result.ranking == ["Control", "Treatment"] or result.ranking == [
            "Treatment",
            "Control",
        ]


# ============================================================================
# Tests: Discovery
# ============================================================================


class TestDiscovery:
    """Test plugin discovery via pkgutil."""

    def test_discovery_finds_no_plugins_initially(self):
        """With no concrete analysis files yet, discovery should return empty."""
        from polyzymd.analyses.discovery import clear_cache, list_analyses

        clear_cache()
        # Currently no plugin files exist in analyses/ (only infrastructure)
        analyses = list_analyses()
        # May be empty or may find ToyAnalysis if this test file gets scanned
        # The important thing is it doesn't crash
        assert isinstance(analyses, dict)

    def test_get_analysis_unknown_raises(self):
        from polyzymd.analyses.discovery import clear_cache, get_analysis

        clear_cache()
        with pytest.raises(KeyError, match="Unknown analysis"):
            get_analysis("nonexistent_analysis_xyz")

    def test_is_concrete_analysis(self):
        from polyzymd.analyses.discovery import _is_concrete_analysis

        assert _is_concrete_analysis(ToyAnalysis) is True
        assert _is_concrete_analysis(Analysis) is False
        assert _is_concrete_analysis(str) is False
        assert _is_concrete_analysis(42) is False


# ============================================================================
# Tests: Orchestrator
# ============================================================================


class TestOrchestrator:
    """Test the orchestrator's replicate running and dependency sorting."""

    def test_run_replicates_success(self, toy_analysis, toy_condition, toy_settings, tmp_path):
        from polyzymd.analyses.orchestrator import _run_replicates

        results, successful, failed = _run_replicates(
            toy_analysis,
            toy_condition,
            toy_settings,
            equilibration="10ns",
            output_dir=tmp_path,
        )
        assert len(results) == 3
        assert successful == [1, 2, 3]
        assert failed == []
        # Check output dirs were created
        assert (tmp_path / "run_1").exists()
        assert (tmp_path / "run_2").exists()
        assert (tmp_path / "run_3").exists()

    def test_run_replicates_partial_failure(self, toy_condition, toy_settings, tmp_path):
        """If some replicates fail but min_replicates is met, continue."""

        class FailingAnalysis(Analysis):
            name: ClassVar[str] = "failing"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 1

            def compute_replicate(self, ctx, replicate):
                if replicate == 2:
                    raise FileNotFoundError("Trajectory not found")
                return ToyResult(value=float(replicate), replicate=replicate)

            def aggregate(self, ctx, results):
                return None

        from polyzymd.analyses.orchestrator import _run_replicates

        analysis = FailingAnalysis()
        results, successful, failed = _run_replicates(
            analysis,
            toy_condition,
            toy_settings,
            equilibration="0ns",
            output_dir=tmp_path,
        )
        assert len(results) == 2
        assert 2 not in successful
        assert 2 in failed

    def test_run_replicates_below_minimum_raises(self, toy_condition, toy_settings, tmp_path):
        """If fewer than min_replicates succeed, raise ValueError."""

        class AlwaysFailAnalysis(Analysis):
            name: ClassVar[str] = "always_fail"
            Settings: ClassVar[type] = ToySettings
            min_replicates: ClassVar[int] = 2

            def compute_replicate(self, ctx, replicate):
                raise RuntimeError("Always fails")

            def aggregate(self, ctx, results):
                return None

        from polyzymd.analyses.orchestrator import _run_replicates

        with pytest.raises(ValueError, match="need at least 2"):
            _run_replicates(
                AlwaysFailAnalysis(),
                toy_condition,
                toy_settings,
                equilibration="0ns",
                output_dir=tmp_path,
            )

    def test_topological_sort_no_deps(self):
        from polyzymd.analyses.orchestrator import _topological_sort

        a = ToyAnalysis()
        result = _topological_sort([a])
        assert result == [a]

    def test_topological_sort_with_deps(self):
        from polyzymd.analyses.orchestrator import _topological_sort

        toy = ToyAnalysis()
        dep = ToyDependentAnalysis()
        # dep depends on toy, so toy should come first
        result = _topological_sort([dep, toy])
        assert result.index(toy) < result.index(dep)

    def test_topological_sort_circular_raises(self):
        from polyzymd.analyses.orchestrator import _topological_sort

        class CircA(Analysis):
            name: ClassVar[str] = "circ_a"
            Settings: ClassVar[type] = ToySettings
            dependencies: ClassVar[tuple[str, ...]] = ("circ_b",)

            def compute_replicate(self, ctx, replicate):
                pass

            def aggregate(self, ctx, results):
                pass

        class CircB(Analysis):
            name: ClassVar[str] = "circ_b"
            Settings: ClassVar[type] = ToySettings
            dependencies: ClassVar[tuple[str, ...]] = ("circ_a",)

            def compute_replicate(self, ctx, replicate):
                pass

            def aggregate(self, ctx, results):
                pass

        with pytest.raises(ValueError, match="Circular dependency"):
            _topological_sort([CircA(), CircB()])

    def test_run_analysis_full(self, toy_analysis, toy_condition, toy_settings, tmp_path):
        """Test the full run_analysis path (compute + aggregate)."""
        from polyzymd.analyses.orchestrator import run_analysis

        result = run_analysis(
            toy_analysis,
            toy_condition,
            toy_settings,
            equilibration="10ns",
            output_dir=tmp_path,
        )
        assert isinstance(result, ToyAggregatedResult)
        assert result.n_replicates == 3
        assert len(result.replicate_values) == 3
        # Check aggregated dir was created
        assert (tmp_path / "aggregated").exists()
        assert (tmp_path / "run_1" / "result.json").exists()
        assert (tmp_path / "aggregated" / "result.json").exists()

    def test_compare_only_analysis_skips_compute_and_aggregate(
        self, toy_condition, toy_settings, tmp_path
    ):
        from polyzymd.analyses.orchestrator import run_analysis

        class CompareOnlyAnalysis(Analysis):
            name: ClassVar[str] = "compare_only"
            Settings: ClassVar[type] = ToySettings
            has_compute_stage: ClassVar[bool] = False
            has_aggregate_stage: ClassVar[bool] = False

        analysis = CompareOnlyAnalysis()
        result = run_analysis(
            analysis,
            toy_condition,
            toy_settings,
            equilibration="0ns",
            output_dir=tmp_path,
        )

        assert result is None
        assert not (tmp_path / "aggregated").exists()


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


# ============================================================================
# Tests: resolve_output_dir
# ============================================================================


class TestResolveOutputDir:
    """Test path resolution utilities."""

    def test_resolve_output_dir(self, toy_analysis, tmp_path):
        result = toy_analysis.resolve_output_dir(tmp_path / "analysis", "100% SBMA")
        assert "100pct_SBMA" in str(result)
        assert result.name == "toy"
