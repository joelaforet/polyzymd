"""Hypothesis testing must behave the same way in every comparison plugin.

Each plugin that still overrides ``compare`` runs its own pairwise tests, so
these tests pin the three properties that must not depend on which plugin was
asked. A plugin on the observable contract does not override ``compare``; the
same three properties are pinned once for all of them in
``tests/analyses/test_contract.py``. The properties are:

1. ``ttest_method`` from the comparison config selects the variance
   assumption. With equal group sizes and unequal variances Welch's test
   gives the same t statistic as Student's but fewer degrees of freedom,
   so the Welch p-value is strictly larger.
2. Every pairwise result carries ``p_value_adjusted``. With a single test
   in the Benjamini-Hochberg family the adjusted p-value equals the raw one.
3. Effect sizes carry Hedges' g and drop the Cohen adjectives when the
   combined sample is smaller than ten replicates.

The fixtures below build condition artifacts by hand. Group A has a
standard deviation of 0.1 and group B a standard deviation of 1.0, so the
variances differ by a factor of 100.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock

import numpy as np
import pytest

from polyzymd.analyses._framework.cache_identity import settings_fingerprint
from polyzymd.analyses.base import ComparisonContext, Condition
from polyzymd.analyses.mda import ConditionArtifact

LOW_VARIANCE = (10.0, 10.1, 9.9)
HIGH_VARIANCE = (12.0, 13.0, 11.0)


# ---------------------------------------------------------------------------
# Context helpers
# ---------------------------------------------------------------------------


def _condition(label: str) -> Condition:
    """Build a condition with three replicates and a mock simulation config."""
    return Condition(
        label=label,
        config_path=Path(f"/fake/{label}/config.yaml"),
        replicates=(1, 2, 3),
        sim_config=MagicMock(),
    )


def _context(
    tmp_path: Path,
    settings: Any,
    aggregated: dict[str, ConditionArtifact],
    ttest_method: str,
) -> ComparisonContext:
    """Build a comparison context that requests a specific t-test method."""
    labels = list(aggregated)
    analysis_dirs = {}
    for label in labels:
        analysis_dir = tmp_path / label
        (analysis_dir / "aggregated").mkdir(parents=True, exist_ok=True)
        analysis_dirs[label] = analysis_dir
    return ComparisonContext(
        name="hypothesis_testing_consistency",
        conditions=[_condition(label) for label in labels],
        excluded_conditions=[],
        control_label=labels[0],
        analysis_dirs=analysis_dirs,
        results_dir=tmp_path / "comparison",
        equilibration="10ns",
        settings=settings,
        ttest_method=ttest_method,
        aggregated_results=aggregated,
    )


def _base_metadata(settings: Any) -> dict[str, Any]:
    """Return artifact metadata accepted by the aggregate validators."""
    return {
        "settings_fingerprint": settings_fingerprint(settings),
        "config_hash": "hash",
        "polyzymd_version": "test",
        "equilibration_time": 10.0,
        "equilibration_unit": "ns",
    }


# ---------------------------------------------------------------------------
# Per-plugin fixtures
# ---------------------------------------------------------------------------


def _rmsd_case(tmp_path: Path, ttest_method: str) -> tuple[Any, ComparisonContext]:
    """The ported rmsd plugin, whose comparison runs through the contract."""
    from polyzymd.analyses.contract import Observable, aggregate_observables
    from polyzymd.analyses.rmsd import RMSDAnalysis, RMSDRunSettings, RMSDSettings

    settings = RMSDSettings(runs=[RMSDRunSettings(label="protein_backbone")])

    def artifact(label: str, values: tuple[float, ...]) -> ConditionArtifact:
        aggregates = aggregate_observables(
            [
                [
                    Observable(
                        name="rmsd_protein_backbone_ref_centroid",
                        kind="mean_of_timeseries",
                        unit="A",
                        values=[value],
                    )
                ]
                for value in values
            ]
        )
        return ConditionArtifact(
            analysis_name="rmsd",
            condition_label=label,
            replicates=[1, 2, 3],
            payload={
                "observables": [aggregate.model_dump(mode="json") for aggregate in aggregates]
            },
            metadata=_base_metadata(settings),
            provenance={"frame_selection": {"equilibration": "10ns"}},
        )

    aggregated = {
        "Control": artifact("Control", LOW_VARIANCE),
        "Treated": artifact("Treated", HIGH_VARIANCE),
    }
    return RMSDAnalysis(), _context(tmp_path, settings, aggregated, ttest_method)


CASES = {
    "rmsd": _rmsd_case,
}


def _pairwise_p_values(plugin: str, result: Any) -> list[tuple[float, float | None]]:
    """Return ``(p_value, p_value_adjusted)`` for every testable pairwise test."""
    pairs: list[tuple[float, float | None]] = []
    if plugin == "rmsd":
        for comparison in result.payload["comparisons"]:
            pairs.append((comparison["p_value"], comparison["p_adjusted"]))
    else:
        for comparison in result.pairwise_comparisons:
            pairs.append((comparison.p_value, comparison.p_value_adjusted))
    return pairs


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("plugin", sorted(CASES))
def test_welch_gives_larger_p_than_student(plugin: str, tmp_path: Path) -> None:
    """Welch's test must widen the p-value when the variances differ 100-fold."""
    analysis, ctx = CASES[plugin](tmp_path / "welch", "welch")
    welch = _pairwise_p_values(plugin, analysis.compare(ctx))

    analysis, ctx = CASES[plugin](tmp_path / "student", "student")
    student = _pairwise_p_values(plugin, analysis.compare(ctx))

    assert welch, f"{plugin} produced no pairwise comparisons"
    assert len(welch) == len(student)
    for (welch_p, _), (student_p, _) in zip(welch, student, strict=True):
        assert welch_p != pytest.approx(student_p)
        assert welch_p > student_p


@pytest.mark.parametrize("plugin", sorted(CASES))
def test_pairwise_results_carry_adjusted_p_values(plugin: str, tmp_path: Path) -> None:
    """Every plugin must report a Benjamini-Hochberg adjusted p-value."""
    analysis, ctx = CASES[plugin](tmp_path, "welch")
    pairs = _pairwise_p_values(plugin, analysis.compare(ctx))

    assert pairs, f"{plugin} produced no pairwise comparisons"
    for raw_p, adjusted_p in pairs:
        assert adjusted_p is not None
        assert adjusted_p >= raw_p


@pytest.mark.parametrize("plugin", ["rmsd"])
def test_single_comparison_leaves_p_value_unchanged(plugin: str, tmp_path: Path) -> None:
    """A family of one test must have an adjusted p-value equal to the raw one."""
    analysis, ctx = CASES[plugin](tmp_path, "welch")
    pairs = _pairwise_p_values(plugin, analysis.compare(ctx))

    assert len(pairs) == 1
    raw_p, adjusted_p = pairs[0]
    assert adjusted_p == pytest.approx(raw_p)


def test_hedges_g_corrects_cohens_d_downward() -> None:
    """Hedges' g must apply the 1981 bias correction J to Cohen's d."""
    from polyzymd.analyses.shared.inferential_statistics import cohens_d

    effect = cohens_d(LOW_VARIANCE, HIGH_VARIANCE)

    n1 = n2 = 3
    expected_j = 1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0)
    assert effect.hedges_g == pytest.approx(effect.cohens_d * expected_j)
    assert abs(effect.hedges_g) < abs(effect.cohens_d)


def test_effect_size_adjectives_dropped_for_small_samples() -> None:
    """Cohen's adjectives are noise below ten replicates, so they are withheld."""
    from polyzymd.analyses.shared.inferential_statistics import cohens_d

    small = cohens_d([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
    assert small.interpretation is None

    large_group_a = [float(value) for value in range(6)]
    large_group_b = [float(value) + 10.0 for value in range(6)]
    large = cohens_d(large_group_a, large_group_b)
    assert large.interpretation == "large"


def test_public_functions_have_resolvable_annotations() -> None:
    """Every public function in the statistics module must be introspectable.

    Ruff does not flag an undefined name in an annotation here because F821
    is suppressed for this project, so a missing import only shows up when
    something resolves the annotations.
    """
    import inspect
    import typing

    from polyzymd.analyses.shared import inferential_statistics

    functions = [
        obj
        for name, obj in vars(inferential_statistics).items()
        if not name.startswith("_")
        and inspect.isfunction(obj)
        and obj.__module__ == inferential_statistics.__name__
    ]

    assert functions
    for function in functions:
        typing.get_type_hints(function)
