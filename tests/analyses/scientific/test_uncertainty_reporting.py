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

from polyzymd.analyses.contract import ObservableEstimate, aggregate_observables
from polyzymd.analyses.exceptions import AnalysisError, StatisticsError
from polyzymd.analyses.mda.artifacts import ConditionArtifact, ReplicateArtifact
from polyzymd.analyses.shared.statistics import (
    compute_sem,
    mean_sem_ci,
    student_t_coverage_factor,
    uncertainty_block,
)

# Two-sided Student t coverage factors at 95 percent, quoted in the article.
T_FACTOR_N3 = 4.302652729749462
T_FACTOR_N5 = 2.7764451051977934

PLUGIN_NAMES = (
    "rmsd",
    "rmsf",
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

    def test_contract_aggregation_declares_uncertainty_and_interval(self) -> None:
        """The shared aggregation path fills the interval and names its method."""

        replicates = [
            [
                ObservableEstimate(
                    name="mean_value", kind="mean_of_timeseries", unit="A", value=v, n_frames=10
                )
            ]
            for v in (2.0, 2.2, 2.4)
        ]

        aggregate = aggregate_observables(replicates)[0]

        assert aggregate.n_replicates == 3
        assert aggregate.unit == "A"
        assert aggregate.ci_method == "student_t"
        assert aggregate.mean == pytest.approx(2.2)
        assert aggregate.ci95_high == pytest.approx(2.2 + T_FACTOR_N3 * aggregate.sem)

    def test_single_replicate_aggregate_writes_null_not_zero(self) -> None:
        """A one-replicate condition must not claim zero uncertainty."""

        replicates = [
            [
                ObservableEstimate(
                    name="mean_value", kind="mean_of_timeseries", unit="A", value=2.0, n_frames=10
                )
            ]
        ]

        aggregate = aggregate_observables(replicates)[0]

        assert aggregate.sem is None
        assert aggregate.ci95_low is None
        assert aggregate.ci95_high is None

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


class TestAggregateWithoutReplicates:
    """An aggregate that names no replicates reports a null count."""

    def test_build_without_any_replicate_count_reports_null(self) -> None:
        """An aggregate that states no replicate count says so, rather than zero."""

        artifact = ConditionArtifact.build(
            analysis_name="demo",
            condition_label="Control",
            payload={"metrics": {}},
        )

        assert artifact.payload["uncertainty"]["n"] is None
