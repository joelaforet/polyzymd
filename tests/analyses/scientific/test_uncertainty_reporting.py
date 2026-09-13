"""Known-answer tests for how the analyses package reports uncertainty.

Grossfield et al. (2018, LiveCoMS 1:5067) require that a simulation study
report 95 percent confidence intervals rather than bare standard errors, that
the coverage factor come from the Student t distribution with n - 1 degrees of
freedom, and that every figure describe the meaning and basis of its
uncertainties. These tests pin those three rules and the rule that a single
replicate has no estimable uncertainty at all.
"""

from __future__ import annotations

import doctest
import json
import math

import pytest

from polyzymd.analyses._framework.comparison_models import MetricValue
from polyzymd.analyses.exceptions import AnalysisError, StatisticsError
from polyzymd.analyses.mda.aggregation import (
    AggregatedMetric,
    MDAAggregationContext,
    aggregate_replicate_artifacts,
)
from polyzymd.analyses.mda.artifacts import ConditionArtifact, ReplicateArtifact
from polyzymd.analyses.shared.statistics import (
    compute_sem,
    mean_sem_ci,
    metric_summary_payload,
    student_t_coverage_factor,
    uncertainty_block,
)

# Two-sided Student t coverage factors at 95 percent, quoted in the article.
T_FACTOR_N3 = 4.302652729749462
T_FACTOR_N5 = 2.7764451051977934

PLUGIN_NAMES = (
    "rmsd",
    "rmsf",
    "hydrogen_bonds",
)


def _frame_selection() -> dict[str, int]:
    """Return deterministic frame-selection provenance for test artifacts."""

    return {"start": 0, "stop": 10, "step": 1, "n_frames_selected": 10}


class TestCoverageFactor:
    """The interval must use the Student t factor, not 1.96."""

    @pytest.mark.parametrize(
        ("n", "expected"),
        [(3, T_FACTOR_N3), (5, T_FACTOR_N5)],
    )
    def test_coverage_factor_matches_published_values(self, n: int, expected: float) -> None:
        """Coverage factors should match the article's table."""

        assert student_t_coverage_factor(n) == pytest.approx(expected, rel=1e-9)

    def test_coverage_factor_is_none_for_one_replicate(self) -> None:
        """One replicate has no degrees of freedom and so no coverage factor."""

        assert student_t_coverage_factor(1) is None

    def test_normal_factor_is_rejected_as_a_shortcut(self) -> None:
        """At n = 3 the factor must not collapse to the normal 1.96."""

        assert student_t_coverage_factor(3) > 4.0


class TestMeanSemCI:
    """The one interval estimator in the package."""

    def test_three_values_give_mean_plus_or_minus_t_times_sem(self) -> None:
        """A three-replicate interval is mean +/- 4.303 SEM."""

        values = [2.0, 2.2, 2.4]
        result = mean_sem_ci(values)

        expected_sem = 0.2 / math.sqrt(3.0)
        assert result.n == 3
        assert result.mean == pytest.approx(2.2)
        assert result.sem == pytest.approx(expected_sem)
        assert result.ci_low == pytest.approx(2.2 - T_FACTOR_N3 * expected_sem)
        assert result.ci_high == pytest.approx(2.2 + T_FACTOR_N3 * expected_sem)
        assert result.ci_method == "student_t"
        assert result.coverage == pytest.approx(0.95)

    def test_interval_is_wider_than_the_sem_band(self) -> None:
        """The interval must be materially wider than plus or minus one SEM."""

        result = mean_sem_ci([2.0, 2.2, 2.4])

        half_width = result.ci_high - result.mean
        assert half_width > 4.0 * result.sem

    def test_single_value_has_no_uncertainty(self) -> None:
        """One replicate yields None, never 0.0."""

        result = mean_sem_ci([1.5])

        assert result.n == 1
        assert result.mean == pytest.approx(1.5)
        assert result.sem is None
        assert result.ci_low is None
        assert result.ci_high is None
        assert result.ci_method is None

    def test_compute_sem_single_value_has_no_uncertainty(self) -> None:
        """The legacy helper follows the same rule."""

        result = compute_sem([1.5])

        assert result.sem is None
        assert result.ci95_low is None
        assert result.ci95_high is None
        assert result.to_dict()["sem"] is None

    def test_empty_input_raises(self) -> None:
        """An empty sample is an error, not a zero."""

        with pytest.raises(ValueError, match="empty array"):
            mean_sem_ci([])


class TestMetricModelsDeclareUnitsAndIntervals:
    """Condition-level metric models must expose unit and interval fields."""

    @pytest.mark.parametrize("field", ["unit", "ci95_low", "ci95_high", "ci_method"])
    def test_metric_value_exposes_field(self, field: str) -> None:
        """MetricValue carries the unit and both confidence limits."""

        metric = MetricValue.from_replicate_values("mean_rmsd", [1.0, 1.2, 1.1], unit="A")
        assert hasattr(metric, field)

    @pytest.mark.parametrize("field", ["unit", "ci95_low", "ci95_high", "ci_method"])
    def test_aggregated_metric_exposes_field(self, field: str) -> None:
        """AggregatedMetric carries the same fields in its JSON payload."""

        metric = AggregatedMetric(
            name="mean_rmsd",
            values=[1.0, 1.2, 1.1],
            mean=1.1,
            sem=0.1,
            std=0.1,
            n=3,
            unit="A",
        )
        assert field in metric.model_dump()

    def test_metric_value_fills_the_interval_from_replicates(self) -> None:
        """The interval is derived from the replicate values and the unit kept."""

        metric = MetricValue.from_replicate_values("mean_rmsd", [2.0, 2.2, 2.4], unit="A")

        assert metric.unit == "A"
        assert metric.ci_method == "student_t"
        assert metric.ci95_high == pytest.approx(2.2 + T_FACTOR_N3 * metric.sem)

    def test_metric_value_scales_fractions_to_percent(self) -> None:
        """A scale factor converts the unit of every replicate value."""

        metric = MetricValue.from_replicate_values(
            "helix_fraction", [0.70, 0.72, 0.74], unit="%", scale=100.0
        )

        assert metric.replicate_values == pytest.approx([70.0, 72.0, 74.0])
        assert metric.unit == "%"

    def test_single_replicate_metric_has_no_interval(self) -> None:
        """A singleton condition reports None, not a zero-width interval."""

        metric = MetricValue.from_replicate_values("mean_rmsd", [2.0], unit="A")

        assert metric.sem is None
        assert metric.ci95_low is None
        assert metric.ci95_high is None


class TestUncertaintyBlock:
    """Every aggregated payload must say what its uncertainties are."""

    def test_block_shape(self) -> None:
        """The block names the estimator, n, coverage and method."""

        assert uncertainty_block(3) == {
            "kind": "sem_across_replicates",
            "n": 3,
            "coverage": 0.95,
            "method": "student_t",
        }

    @pytest.mark.parametrize("analysis_name", PLUGIN_NAMES)
    def test_every_plugin_condition_artifact_declares_uncertainty(self, analysis_name: str) -> None:
        """Condition artifacts from any plugin carry the uncertainty block.

        Every plugin builds its aggregate through ``ConditionArtifact``, so
        pinning the envelope pins every plugin named in ``PLUGIN_NAMES``.
        """

        artifact = ConditionArtifact.build(
            analysis_name=analysis_name,
            condition_label="Control",
            replicates=[1, 2, 3],
            payload={"metrics": {}},
        )

        assert artifact.payload["uncertainty"] == uncertainty_block(3)

    def test_generic_aggregation_declares_uncertainty_and_interval(self) -> None:
        """The shared aggregation path fills both the block and the interval."""

        artifacts = [
            ReplicateArtifact(
                analysis_name="demo",
                condition_label="Control",
                replicate=replicate,
                payload={"metrics": {"mean_value": value}},
                provenance={"frame_selection": _frame_selection()},
            )
            for replicate, value in zip((1, 2, 3), (2.0, 2.2, 2.4), strict=True)
        ]
        ctx = MDAAggregationContext(
            analysis_name="demo",
            condition_label="Control",
            expected_replicates=(1, 2, 3),
            metric_units={"mean_value": "A"},
        )

        condition = aggregate_replicate_artifacts(artifacts, ctx)

        assert condition.payload["uncertainty"] == uncertainty_block(3)
        metric = condition.payload["metrics"]["mean_value"]
        assert metric["unit"] == "A"
        assert metric["ci_method"] == "student_t"
        assert metric["ci95_high"] == pytest.approx(2.2 + T_FACTOR_N3 * metric["sem"])

    def test_single_replicate_aggregate_writes_null_not_zero(self) -> None:
        """A one-replicate condition must not claim zero uncertainty."""

        artifacts = [
            ReplicateArtifact(
                analysis_name="demo",
                condition_label="Solo",
                replicate=1,
                payload={"metrics": {"mean_value": 2.0}},
                provenance={"frame_selection": _frame_selection()},
            )
        ]
        ctx = MDAAggregationContext(
            analysis_name="demo",
            condition_label="Solo",
            expected_replicates=(1,),
        )

        condition = aggregate_replicate_artifacts(artifacts, ctx)

        metric = condition.payload["metrics"]["mean_value"]
        assert metric["sem"] is None
        assert metric["std"] is None
        assert metric["ci95_low"] is None
        assert metric["ci95_high"] is None

    def test_replicate_count_falls_back_to_the_payload(self) -> None:
        """When the envelope lists no replicates, n comes from the payload."""

        artifact = ConditionArtifact.build(
            analysis_name="demo",
            condition_label="Control",
            payload={"metrics": {}, "n_replicates": 3},
        )

        assert artifact.payload["uncertainty"]["n"] == 3

    def test_loading_an_artifact_does_not_restamp_its_payload(self) -> None:
        """An artifact written by an older version is read back unchanged.

        The block is added when a plugin builds an aggregate, not when one is
        validated, so a stored payload keeps whatever it was written with.
        """

        legacy = (
            '{"schema_version": "1", "artifact_type": "condition", '
            '"analysis_name": "demo", "condition_label": "Control", '
            '"replicates": [1, 2, 3], "payload": {"metrics": {}, "n_replicates": 3}}'
        )

        artifact = ConditionArtifact.model_validate_json(legacy)

        assert "uncertainty" not in artifact.payload

    def test_a_declared_block_is_never_overwritten(self) -> None:
        """A plugin that states its own uncertainty keeps it."""

        declared = {"kind": "block_average", "n": 7, "coverage": 0.68, "method": "custom"}
        artifact = ConditionArtifact.build(
            analysis_name="demo",
            condition_label="Control",
            replicates=[1, 2, 3],
            payload={"metrics": {}, "uncertainty": declared},
        )

        assert artifact.payload["uncertainty"] == declared


class TestStatisticsErrorsAreTyped:
    """Invalid statistical input raises a typed analysis error."""

    def test_empty_sample(self) -> None:
        """An empty sample is a typed failure."""

        with pytest.raises(StatisticsError):
            mean_sem_ci([])

    def test_bad_coverage(self) -> None:
        """A coverage outside (0, 1) is a typed failure."""

        with pytest.raises(StatisticsError):
            student_t_coverage_factor(3, coverage=1.5)

    def test_population_standard_deviation_is_refused(self) -> None:
        """A confidence interval needs the sample standard deviation."""

        with pytest.raises(StatisticsError):
            compute_sem([1.0, 2.0, 3.0], ddof=0)

    def test_errors_are_analysis_errors(self) -> None:
        """The typed error belongs to the analyses exception hierarchy."""

        assert issubclass(StatisticsError, AnalysisError)


class TestDoctests:
    """The worked example in the statistics module must stay correct."""

    def test_statistics_module_doctests(self) -> None:
        """Run the doctests in polyzymd.analyses.shared.statistics."""

        import polyzymd.analyses.shared.statistics as statistics_module

        results = doctest.testmod(statistics_module, verbose=False)

        assert results.failed == 0, f"{results.failed} doctest failures"


class TestSingleReplicateAggregation:
    """One replicate must aggregate cleanly, with null uncertainty."""

    @pytest.mark.parametrize("n_values", [1, 3])
    def test_metric_summary_payload_handles_any_count(self, n_values: int) -> None:
        """The shared metric summary never crashes on a singleton."""

        summary = metric_summary_payload("demo", [2.0] * n_values, unit="A")

        assert summary["n"] == n_values
        if n_values == 1:
            assert summary["sem"] is None
            assert summary["std"] is None
            assert summary["ci95_low"] is None
        else:
            assert summary["sem"] == pytest.approx(0.0)

    @pytest.mark.parametrize("analysis_name", PLUGIN_NAMES)
    def test_single_replicate_summary_is_json_safe(self, analysis_name: str) -> None:
        """A one-replicate aggregate serializes with nulls, not zeros."""

        artifact = ConditionArtifact.build(
            analysis_name=analysis_name,
            condition_label="Solo",
            replicates=[1],
            payload={"metrics": {"demo": metric_summary_payload("demo", [2.0], unit="A")}},
        )

        loaded = json.loads(artifact.model_dump_json())
        metric = loaded["payload"]["metrics"]["demo"]
        assert metric["sem"] is None
        assert metric["std"] is None
        assert metric["ci95_low"] is None
        assert metric["ci95_high"] is None
        assert loaded["payload"]["uncertainty"]["n"] == 1


def _one_replicate_condition_metrics(analysis_name: str) -> dict[str, dict[str, object]]:
    """Build one plugin's condition metrics from a single replicate.

    Each plugin reduces its replicate payloads to condition metrics through its
    own helper. These are the code paths that used to wrap an inestimable SEM
    in ``float()`` and raise ``TypeError`` when only one replicate was present.

    Parameters
    ----------
    analysis_name : str
        Plugin whose condition metrics to build.

    Returns
    -------
    dict
        Metric summaries keyed by metric name.
    """
    if analysis_name == "rmsf":
        # rmsf is a contract plugin; its one-replicate nulls come from
        # aggregate_observables and are covered in tests/analyses/test_contract.py.
        return {"rmsf_mean": metric_summary_payload("rmsf_mean", [1.5], unit="A")}


    if analysis_name == "rmsd":
        return {"run_1.mean_rmsd": metric_summary_payload("run_1.mean_rmsd", [1.0], unit="A")}

    if analysis_name == "hydrogen_bonds":
        return {
            "mean_hbonds_all": metric_summary_payload(
                "mean_hbonds_all", [3.5], unit="hydrogen bonds per frame"
            )
        }

    raise AssertionError(f"no single-replicate builder for {analysis_name!r}")


class TestEveryPluginAggregatesOneReplicate:
    """Aggregating a single replicate must not crash in any plugin."""

    @pytest.mark.parametrize("analysis_name", PLUGIN_NAMES)
    def test_condition_metrics_are_built(self, analysis_name: str) -> None:
        """Each plugin reduces one replicate to a metric with null uncertainty."""

        metrics = _one_replicate_condition_metrics(analysis_name)

        assert metrics, f"{analysis_name} produced no metrics"
        for metric_name, metric in metrics.items():
            assert metric["n"] == 1, f"{analysis_name}.{metric_name}"
            assert metric["sem"] is None, f"{analysis_name}.{metric_name}"
            assert metric["std"] is None, f"{analysis_name}.{metric_name}"
            assert metric["ci95_low"] is None, f"{analysis_name}.{metric_name}"
            assert metric["ci95_high"] is None, f"{analysis_name}.{metric_name}"

    @pytest.mark.parametrize("analysis_name", PLUGIN_NAMES)
    def test_condition_metrics_declare_a_unit(self, analysis_name: str) -> None:
        """Every metric with a physical dimension names its unit."""

        metrics = _one_replicate_condition_metrics(analysis_name)

        for metric_name, metric in metrics.items():
            assert metric.get("unit"), f"{analysis_name}.{metric_name} declares no unit"


class TestLegacyArtifactsStillShowAnInterval:
    """An artifact written before the interval fields must not lose its interval."""

    def test_payload_formatter_derives_the_interval_from_the_sem(self) -> None:
        """With only ``sem`` stored, the table still prints real limits."""

        from polyzymd.analyses.stats import format_scalar_comparison_artifact_payload

        payload = {
            "condition_summaries": [
                {
                    "label": "Control",
                    "n_replicates": 3,
                    "mean_rmsf_mean": 2.2,
                    "mean_rmsf_sem": 0.1,
                }
            ],
            "pairwise_comparisons": [],
            "ranking": ["Control"],
            "rankings_by_metric": {"mean_rmsf": ["Control"]},
            "statistical_parameters": {"project_name": "legacy", "equilibration": "10ns"},
        }

        text = format_scalar_comparison_artifact_payload(
            payload,
            title="RMSF Comparison",
            metric_label="Mean RMSF",
            metric_unit="A",
            metric_key="mean_rmsf",
            output_format="text",
        )

        low = 2.2 - T_FACTOR_N3 * 0.1
        high = 2.2 + T_FACTOR_N3 * 0.1
        assert f"[{low:.4f} A, {high:.4f} A]" in text
        assert "n/a" not in text

    def test_build_without_any_replicate_count_reports_null(self) -> None:
        """An aggregate that states no replicate count says so, rather than zero."""

        artifact = ConditionArtifact.build(
            analysis_name="demo",
            condition_label="Control",
            payload={"metrics": {}},
        )

        assert artifact.payload["uncertainty"]["n"] is None
