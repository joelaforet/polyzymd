"""Tests for the analyses plugin infrastructure (Phase A).

Tests the base class contract, discovery, stats utilities, and orchestrator
without requiring heavy dependencies (OpenMM, MDAnalysis, etc.).
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar, Sequence

import pytest
from pydantic import BaseModel, ValidationError

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    Condition,
    MetricValue,
    PlotContext,
    ReplicateContext,
    SlurmResourceHint,
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
    replicates: list[int] | None = None
    settings_fingerprint: str | None = None

    def save(self, path: Path) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path

    @classmethod
    def load(cls, path: Path) -> "ToyAggregatedResult":
        return cls.model_validate_json(path.read_text())


class _MDAContractMixin:
    """Provide the required MDA lifecycle seam for direct compute fakes."""

    def build_mda_jobs(self, ctx):
        """Return no jobs for tests that override the internal dispatcher."""

        del ctx
        return []


class ToyAnalysis(_MDAContractMixin, Analysis):
    """Concrete analysis for testing the plugin system."""

    name: ClassVar[str] = "toy"
    Settings: ClassVar[type] = ToySettings
    AggregatedResultClass: ClassVar[type] = ToyAggregatedResult
    aliases: ClassVar[tuple[str, ...]] = ("toy_alias",)
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    def _run_compute_stage(self, ctx: ReplicateContext, replicate: int) -> ToyResult:
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
            replicates=list(ctx.replicates),
            settings_fingerprint=self.aggregate_settings_fingerprint(ctx.settings),
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
        sim_config=object(),
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

            class BadAnalysis(_MDAContractMixin, Analysis):
                Settings = ToySettings

                def aggregate(self, ctx, results):
                    pass

    def test_concrete_subclass_requires_settings(self):
        """Subclass without 'Settings' should fail at class creation."""
        with pytest.raises(TypeError, match="must define 'Settings'"):

            class BadAnalysis(_MDAContractMixin, Analysis):
                name = "bad"

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

    def test_concrete_subclass_requires_compute_contract(self) -> None:
        """Compute-stage plugins must provide MDA jobs."""

        with pytest.raises(
            TypeError,
            match=r"public plugins must implement build_mda_jobs\(\) when has_compute_stage=True",
        ):

            class BadComputeContractAnalysis(Analysis):
                name: ClassVar[str] = "bad_compute_contract"
                Settings: ClassVar[type] = ToySettings

                def aggregate(self, ctx, results):
                    return {"dummy": True}

    def test_noncanonical_compute_replicate_hook_is_not_supported(self) -> None:
        """Non-canonical compute_replicate-only plugins should not satisfy the contract."""

        assert not hasattr(Analysis, "compute_replicate")

        with pytest.raises(
            TypeError,
            match=r"public plugins must implement build_mda_jobs\(\)",
        ):

            class NoncanonicalComputeOnlyAnalysis(Analysis):
                name: ClassVar[str] = "noncanonical_compute_only"
                Settings: ClassVar[type] = ToySettings

                def compute_replicate(self, ctx, replicate):
                    return {"replicate": replicate}

                def aggregate(self, ctx, results):
                    return {"dummy": True}

    def test_run_replicate_override_is_rejected(self) -> None:
        """Concrete plugins cannot define the removed replicate hook."""

        with pytest.raises(TypeError, match=r"defines removed hook run_replicate\(\)"):

            class RunReplicateAnalysis(_MDAContractMixin, Analysis):
                name: ClassVar[str] = "run_replicate_removed"
                Settings: ClassVar[type] = ToySettings

                def run_replicate(self, ctx: ReplicateContext, replicate: int) -> dict[str, float]:
                    del ctx
                    return {"value": float(replicate)}

                def aggregate(self, ctx, results):
                    del ctx, results
                    return {"dummy": True}

    def test_run_replicate_inherited_from_mixin_is_rejected(self) -> None:
        """Concrete plugins cannot inherit the removed replicate hook from mixins."""

        class RemovedHookMixin:
            def run_replicate(self, ctx: ReplicateContext, replicate: int) -> dict[str, float]:
                """Non-canonical hook that should be rejected after MRO resolution."""

                del ctx
                return {"value": float(replicate)}

        with pytest.raises(TypeError, match=r"inherits removed hook run_replicate\(\)"):

            class InheritedRunReplicateAnalysis(_MDAContractMixin, RemovedHookMixin, Analysis):
                name: ClassVar[str] = "inherited_run_replicate_removed"
                Settings: ClassVar[type] = ToySettings

                def aggregate(self, ctx, results):
                    del ctx, results
                    return {"dummy": True}

    def test_run_replicate_inherited_from_abstract_intermediate_is_rejected(self) -> None:
        """Concrete plugins cannot inherit the removed hook from abstract bases."""

        from abc import abstractmethod

        class AbstractNoncanonicalRunReplicate(Analysis):
            @abstractmethod
            def extra_method(self) -> None:
                """Keep the intermediate class abstract during definition."""

            def run_replicate(self, ctx: ReplicateContext, replicate: int) -> dict[str, float]:
                """Non-canonical hook that should be rejected once concrete."""

                del ctx
                return {"value": float(replicate)}

        assert getattr(AbstractNoncanonicalRunReplicate, "__abstractmethods__", None)

        with pytest.raises(TypeError, match=r"inherits removed hook run_replicate\(\)"):

            class ConcreteNoncanonicalRunReplicate(_MDAContractMixin, AbstractNoncanonicalRunReplicate):
                name: ClassVar[str] = "concrete_inherited_run_replicate_removed"
                Settings: ClassVar[type] = ToySettings

                def extra_method(self) -> None:
                    """Implement the abstract method to trigger concrete validation."""

                def aggregate(self, ctx, results):
                    del ctx, results
                    return {"dummy": True}

    def test_internal_compute_dispatch_uses_mda_lifecycle(
        self, toy_analysis, toy_condition
    ) -> None:
        """The internal compute dispatcher should run the configured compute seam."""
        ctx = ReplicateContext(
            condition=toy_condition,
            replicate=2,
            sim_config=toy_condition.sim_config,
            output_dir=Path("/tmp/run_2"),
            equilibration="10ns",
            recompute=False,
            settings=ToySettings(),
        )

        result = toy_analysis._run_compute_stage(ctx, replicate=2)

        assert result == ToyResult(value=3.0, replicate=2)

    def test_runner_hooks_are_rejected(self) -> None:
        """Removed runner hooks should fail with migration guidance."""

        with pytest.raises(
            TypeError,
            match=r"defines removed runner hook\(s\): build_runner.*Implement build_mda_jobs",
        ):

            class RemovedRunnerHookAnalysis(Analysis):
                name: ClassVar[str] = "removed_runner_hook"
                Settings: ClassVar[type] = ToySettings

                def build_runner(self, ctx, replicate, universe, window):
                    """Non-canonical hook that should be rejected."""

                    del ctx, replicate, universe, window
                    return object()

                def aggregate(self, ctx, results):
                    """Return a simple aggregate result."""

                    del ctx, results
                    return {"dummy": True}

    def test_mda_job_subclass_satisfies_compute_contract(self) -> None:
        """MDA job-backed plugins should satisfy the compute contract."""

        class MDAJobOnlyAnalysis(Analysis):
            name: ClassVar[str] = "mda_job_only_contract"
            Settings: ClassVar[type] = ToySettings

            def build_mda_jobs(self, ctx):
                """Return no jobs for contract-only validation."""

                del ctx
                return []

            def aggregate(self, ctx, results):
                """Return a simple aggregate result."""

                del ctx, results
                return {"dummy": True}

        assert MDAJobOnlyAnalysis().name == "mda_job_only_contract"

    def test_compare_only_subclass_can_disable_compute_stage(self) -> None:
        """Compare-only plugins should remain valid with compute disabled."""

        class CompareOnlyAnalysis(Analysis):
            name: ClassVar[str] = "compare_only"
            Settings: ClassVar[type] = ToySettings
            has_compute_stage: ClassVar[bool] = False
            has_aggregate_stage: ClassVar[bool] = False

        plugin = CompareOnlyAnalysis()
        assert plugin.has_compute_stage is False
        assert plugin.has_aggregate_stage is False

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
        assert base_result == {}
        # ToyAnalysis's own implementation should return metrics
        toy_result = toy_analysis.extract_metrics(
            ToyAggregatedResult(
                mean_value=1.0, sem_value=0.1, replicate_values=[0.9, 1.1], n_replicates=2
            )
        )
        assert "toy_metric" in toy_result

    def test_default_slurm_resource_hint_is_none(self, toy_analysis) -> None:
        """Default slurm_resource_hint should be unset."""
        assert toy_analysis.slurm_resource_hint is None

    def test_subclass_can_set_slurm_resource_hint(self) -> None:
        """Subclass should be able to provide SLURM defaults."""

        class SlurmHintAnalysis(_MDAContractMixin, Analysis):
            name: ClassVar[str] = "slurm_hint"
            Settings: ClassVar[type] = ToySettings
            slurm_resource_hint: ClassVar[SlurmResourceHint | None] = SlurmResourceHint(
                mem="16G",
                time="04:00:00",
                cpus_per_task=4,
            )

            def _run_compute_stage(self, ctx: ReplicateContext, replicate: int) -> dict[str, float]:
                del ctx
                return {"value": float(replicate)}

            def aggregate(
                self,
                ctx: AggregateContext,
                results: Sequence[dict[str, float]],
            ) -> dict[str, float]:
                return {"mean": 1.0}

        plugin = SlurmHintAnalysis()
        assert plugin.slurm_resource_hint is not None
        assert plugin.slurm_resource_hint.mem == "16G"
        assert plugin.slurm_resource_hint.time == "04:00:00"
        assert plugin.slurm_resource_hint.cpus_per_task == 4

    def test_slurm_resource_hint_model_validation(self) -> None:
        """SlurmResourceHint should validate declared field types."""
        with pytest.raises(ValidationError):
            SlurmResourceHint(cpus_per_task="four")


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

    def test_aggregate_context_preserves_positional_api(self, toy_condition, toy_settings):
        """Old positional aggregate context construction should remain valid."""
        result_path = Path("/tmp/aggregated/result.json")

        ctx = AggregateContext(
            toy_condition,
            (1, 2),
            Path("/tmp/aggregated"),
            "10ns",
            toy_settings,
            result_path,
        )

        assert ctx.condition is toy_condition
        assert ctx.replicates == (1, 2)
        assert ctx.output_dir == Path("/tmp/aggregated")
        assert ctx.equilibration == "10ns"
        assert ctx.settings is toy_settings
        assert ctx.result_path == result_path
        assert ctx.recompute is False

        recompute_ctx = AggregateContext(
            toy_condition,
            (1, 2),
            Path("/tmp/aggregated"),
            "10ns",
            toy_settings,
            result_path,
            recompute=True,
        )
        assert recompute_ctx.result_path == result_path
        assert recompute_ctx.recompute is True

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

    def test_comparison_context_preserves_positional_api(self, toy_condition):
        """Old positional comparison context construction should remain valid."""
        result_path = Path("/tmp/results/result.json")
        failed_conditions = [toy_condition]
        aggregated_results = {"Test Condition": {"ok": True}}

        ctx = ComparisonContext(
            "Test Project",
            [toy_condition],
            [],
            None,
            {"Test Condition": Path("/tmp/analysis/toy")},
            Path("/tmp/results"),
            "10ns",
            ToySettings(),
            0.1,
            "welch",
            "tukey_hsd",
            result_path,
            failed_conditions,
            aggregated_results,
        )

        assert ctx.name == "Test Project"
        assert ctx.conditions == [toy_condition]
        assert ctx.excluded_conditions == []
        assert ctx.control_label is None
        assert ctx.analysis_dirs == {"Test Condition": Path("/tmp/analysis/toy")}
        assert ctx.results_dir == Path("/tmp/results")
        assert ctx.equilibration == "10ns"
        assert isinstance(ctx.settings, ToySettings)
        assert ctx.fdr_alpha == 0.1
        assert ctx.ttest_method == "welch"
        assert ctx.posthoc_method == "tukey_hsd"
        assert ctx.result_path == result_path
        assert ctx.failed_conditions == failed_conditions
        assert ctx.aggregated_results == aggregated_results
        assert ctx.recompute is False

        recompute_ctx = ComparisonContext(
            "Test Project",
            [toy_condition],
            [],
            None,
            {"Test Condition": Path("/tmp/analysis/toy")},
            Path("/tmp/results"),
            "10ns",
            ToySettings(),
            0.1,
            "welch",
            "tukey_hsd",
            result_path,
            failed_conditions,
            aggregated_results,
            recompute=True,
        )
        assert recompute_ctx.fdr_alpha == 0.1
        assert recompute_ctx.aggregated_results == aggregated_results
        assert recompute_ctx.recompute is True

    def test_plot_context_preserves_positional_api(self, toy_condition):
        """Old positional plot context construction should remain valid."""
        from polyzymd.config.comparison import PlotSettings

        plot_settings = PlotSettings()
        comparison_path = Path("/tmp/results/result.json")

        ctx = PlotContext(
            [toy_condition],
            {"Test Condition": Path("/tmp/analysis/toy")},
            Path("/tmp/results"),
            Path("/tmp/figures/toy"),
            ToySettings(),
            plot_settings,
            comparison_path,
            "Test Condition",
            "10ns",
        )

        assert ctx.conditions == [toy_condition]
        assert ctx.analysis_dirs == {"Test Condition": Path("/tmp/analysis/toy")}
        assert ctx.results_dir == Path("/tmp/results")
        assert ctx.output_dir == Path("/tmp/figures/toy")
        assert isinstance(ctx.settings, ToySettings)
        assert ctx.plot_settings is plot_settings
        assert ctx.comparison_path == comparison_path
        assert ctx.control_label == "Test Condition"
        assert ctx.equilibration == "10ns"
        assert ctx.recompute is False

        recompute_ctx = PlotContext(
            [toy_condition],
            {"Test Condition": Path("/tmp/analysis/toy")},
            Path("/tmp/results"),
            Path("/tmp/figures/toy"),
            ToySettings(),
            plot_settings,
            comparison_path,
            "Test Condition",
            "10ns",
            recompute=True,
        )
        assert recompute_ctx.comparison_path == comparison_path
        assert recompute_ctx.equilibration == "10ns"
        assert recompute_ctx.recompute is True

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

        class EmptyMetricsAnalysis(_MDAContractMixin, Analysis):
            name: ClassVar[str] = "empty_metrics"
            Settings: ClassVar[type] = ToySettings

            def _run_compute_stage(self, ctx, replicate):
                return {"replicate": replicate}

            def aggregate(self, ctx, results):
                return {"dummy": True}

        analysis = EmptyMetricsAnalysis()
        condition = Condition(
            label="A",
            config_path=tmp_path / "a.yaml",
            replicates=(1, 2),
            sim_config=object(),
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
            aggregated_results={
                "A": {
                    "dummy": True,
                    "n_replicates": 2,
                    "settings_fingerprint": analysis.aggregate_settings_fingerprint(ToySettings()),
                }
            },
        )

        with pytest.raises(
            PluginContractError,
            match=r"extract_metrics\(\) returned empty dict for condition 'A'",
        ):
            analysis.compare(ctx)

    def test_compare_raises_when_extract_metrics_returns_non_dict(self, tmp_path: Path) -> None:
        """extract_metrics must return a dict mapping metric names to MetricValue."""

        class BadTypeAnalysis(_MDAContractMixin, Analysis):
            name: ClassVar[str] = "bad_type"
            Settings: ClassVar[type] = ToySettings

            def _run_compute_stage(self, ctx, replicate):
                return {"replicate": replicate}

            def aggregate(self, ctx, results):
                return {"dummy": True}

            def extract_metrics(self, summary):
                del summary
                return ["not", "a", "dict"]

        analysis = BadTypeAnalysis()
        condition = Condition(
            label="A",
            config_path=tmp_path / "a.yaml",
            replicates=(1, 2),
            sim_config=object(),
        )
        ctx = ComparisonContext(
            name="proj",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": tmp_path / "analysis" / "A" / "bad_type"},
            results_dir=tmp_path / "comparison" / "bad_type",
            equilibration="10ns",
            settings=ToySettings(),
            recompute=False,
            aggregated_results={
                "A": {
                    "dummy": True,
                    "n_replicates": 2,
                    "settings_fingerprint": analysis.aggregate_settings_fingerprint(ToySettings()),
                }
            },
        )

        with pytest.raises(
            PluginContractError,
            match=r"extract_metrics\(\) must return dict\[str, MetricValue\]",
        ):
            analysis.compare(ctx)

    def test_compare_raises_when_extract_metrics_contains_non_metric_value(
        self, tmp_path: Path
    ) -> None:
        """extract_metrics values must be MetricValue instances."""

        class BadValueAnalysis(_MDAContractMixin, Analysis):
            name: ClassVar[str] = "bad_value"
            Settings: ClassVar[type] = ToySettings

            def _run_compute_stage(self, ctx, replicate):
                return {"replicate": replicate}

            def aggregate(self, ctx, results):
                return {"dummy": True}

            def extract_metrics(self, summary):
                del summary
                return {"bad_metric": 123}

        analysis = BadValueAnalysis()
        condition = Condition(
            label="A",
            config_path=tmp_path / "a.yaml",
            replicates=(1, 2),
            sim_config=object(),
        )
        ctx = ComparisonContext(
            name="proj",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={"A": tmp_path / "analysis" / "A" / "bad_value"},
            results_dir=tmp_path / "comparison" / "bad_value",
            equilibration="10ns",
            settings=ToySettings(),
            recompute=False,
            aggregated_results={
                "A": {
                    "dummy": True,
                    "n_replicates": 2,
                    "settings_fingerprint": analysis.aggregate_settings_fingerprint(ToySettings()),
                }
            },
        )

        with pytest.raises(
            PluginContractError,
            match=r"returned invalid value for key 'bad_metric'.*expected MetricValue",
        ):
            analysis.compare(ctx)

    def test_compare_skips_missing_result_file_with_warning(self, caplog, tmp_path: Path) -> None:
        """Missing aggregated files should be skipped with warning, not contract error."""

        analysis = ToyAnalysis()
        condition = Condition(
            label="A",
            config_path=tmp_path / "a.yaml",
            replicates=(1, 2),
            sim_config=object(),
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
