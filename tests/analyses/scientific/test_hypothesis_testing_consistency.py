"""Hypothesis testing must behave the same way in every comparison plugin.

Each plugin that overrides ``compare`` runs its own pairwise tests. These
tests pin the three properties that must not depend on which plugin was
asked:

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


def _settings_fingerprint_for(settings: Any) -> str:
    """Return the cache tag the plugin that owns these settings would write."""

    from polyzymd.analyses.rmsd import RMSDSettings

    if isinstance(settings, RMSDSettings):
        from polyzymd.analyses.rmsd import RMSDAnalysis

        return RMSDAnalysis._make_settings_cache_tag(settings)
    return settings_fingerprint(settings)


def _base_metadata(settings: Any) -> dict[str, Any]:
    """Return artifact metadata accepted by the aggregate validators."""
    return {
        "settings_fingerprint": _settings_fingerprint_for(settings),
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


def _rg_case(tmp_path: Path, ttest_method: str) -> tuple[Any, ComparisonContext]:
    from polyzymd.analyses.rg import RgAnalysis, RgRunSettings, RgSettings

    settings = RgSettings(runs=[RgRunSettings(label="protein_rg", selection="protein")])

    def artifact(label: str, values: tuple[float, ...]) -> ConditionArtifact:
        return ConditionArtifact(
            analysis_name="rg",
            condition_label=label,
            replicates=[1, 2, 3],
            payload={
                "runs": [
                    {
                        "run_label": "protein_rg",
                        "selection": "protein",
                        "replicates": [1, 2, 3],
                        "n_replicates": 3,
                        "overall_mean": float(np.mean(values)),
                        "overall_sem": 0.05,
                        "per_replicate_means": list(values),
                        "per_replicate_stds": [0.2, 0.2, 0.2],
                        "per_replicate_medians": list(values),
                        "calculation_mode": "selection",
                        "fragment_weighting": "equal",
                    }
                ],
                "metrics": {},
                "replicate_metrics": {},
                "n_replicates": 3,
            },
            metadata={**_base_metadata(settings), "selection_string": "protein"},
            provenance={"frame_selection": {"equilibration": "10ns"}},
        )

    aggregated = {
        "Control": artifact("Control", LOW_VARIANCE),
        "Treated": artifact("Treated", HIGH_VARIANCE),
    }
    return RgAnalysis(), _context(tmp_path, settings, aggregated, ttest_method)


def _sasa_case(tmp_path: Path, ttest_method: str) -> tuple[Any, ComparisonContext]:
    from polyzymd.analyses.sasa import SASAAnalysis, SASARunSettings, SASASettings

    settings = SASASettings(runs=[SASARunSettings(label="protein", target_selection="chainid A")])

    def artifact(label: str, values: tuple[float, ...]) -> ConditionArtifact:
        return ConditionArtifact(
            analysis_name="sasa",
            condition_label=label,
            replicates=[1, 2, 3],
            payload={
                "run_results": [
                    {
                        "run_label": "protein",
                        "target_selection": "chainid A",
                        "context_selection": "chainid A",
                        "replicates": [1, 2, 3],
                        "n_replicates": 3,
                        "overall_mean": float(np.mean(values)),
                        "overall_sem": 0.05,
                        "per_replicate_means": list(values),
                        "zero_atom_selection": False,
                    }
                ],
                "metrics": {},
                "replicate_metrics": {},
                "n_replicates": 3,
            },
            metadata=_base_metadata(settings),
            provenance={"frame_selection": {"equilibration": "10ns"}},
        )

    aggregated = {
        "Control": artifact("Control", LOW_VARIANCE),
        "Treated": artifact("Treated", HIGH_VARIANCE),
    }
    return SASAAnalysis(), _context(tmp_path, settings, aggregated, ttest_method)


def _distances_case(tmp_path: Path, ttest_method: str) -> tuple[Any, ComparisonContext]:
    from polyzymd.analyses.distances import (
        DistancePairSettings,
        DistancesAnalysis,
        DistancesSettings,
    )

    settings = DistancesSettings(
        pairs=[
            DistancePairSettings(
                label="catalytic_pair",
                selection_a="resid 1 and name CA",
                selection_b="resid 2 and name CA",
            )
        ]
    )

    def artifact(label: str, values: tuple[float, ...]) -> ConditionArtifact:
        pair_results = [
            {
                "pair_label": "catalytic_pair",
                "selection1": "resid 1 and name CA",
                "selection2": "resid 2 and name CA",
                "threshold": None,
                "overall_mean": float(np.mean(values)),
                "overall_sem": 0.05,
                "overall_fraction_below": None,
                "sem_fraction_below": None,
                "per_replicate_means": list(values),
                "per_replicate_fractions_below": [],
                "replicates": [1, 2, 3],
                "n_replicates": 3,
            }
        ]
        return ConditionArtifact(
            analysis_name="distances",
            condition_label=label,
            replicates=[1, 2, 3],
            payload={
                "pair_results": pair_results,
                "pairs": pair_results,
                "n_replicates": 3,
            },
            metadata=_base_metadata(settings),
            provenance={"frame_selection": {"equilibration": "10ns"}},
        )

    aggregated = {
        "Control": artifact("Control", LOW_VARIANCE),
        "Treated": artifact("Treated", HIGH_VARIANCE),
    }
    return DistancesAnalysis(), _context(tmp_path, settings, aggregated, ttest_method)


def _contacts_case(tmp_path: Path, ttest_method: str) -> tuple[Any, ComparisonContext]:
    from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
    from polyzymd.analyses.contacts._identity import contacts_detection_fingerprint

    settings = ContactsSettings()

    def artifact(
        label: str, values: tuple[float, ...], trace: tuple[float, ...]
    ) -> ConditionArtifact:
        # The first residue carries the variance contrast. The second one
        # contributes a negligible amount of contact but drops out in one
        # replicate, so coverage is not degenerate.
        per_residue = [
            [value / 100.0 for value in values],
            list(trace),
        ]
        rows = [
            {
                "protein_resid": index + 1,
                "protein_resname": "ALA",
                "protein_chain_id": "A",
                "protein_group": "nonpolar",
                "contact_fraction_mean": float(np.mean(fractions)),
                "contact_fraction_per_replicate": list(fractions),
            }
            for index, fractions in enumerate(per_residue)
        ]
        coverage = [
            sum(1 for fractions in per_residue if fractions[index] > 0.0) / len(per_residue)
            for index in range(3)
        ]
        contact = [
            float(np.mean([fractions[index] for fractions in per_residue])) for index in range(3)
        ]
        return ConditionArtifact(
            analysis_name="contacts",
            condition_label=label,
            replicates=[1, 2, 3],
            payload={
                "metrics": {
                    "coverage": {
                        "name": "coverage",
                        "values": coverage,
                        "mean": float(np.mean(coverage)),
                        "sem": float(np.std(coverage, ddof=1) / np.sqrt(3)),
                        "std": float(np.std(coverage, ddof=1)),
                        "n": 3,
                    },
                    "mean_contact_fraction": {
                        "name": "mean_contact_fraction",
                        "values": contact,
                        "mean": float(np.mean(contact)),
                        "sem": float(np.std(contact, ddof=1) / np.sqrt(3)),
                        "std": float(np.std(contact, ddof=1)),
                        "n": 3,
                    },
                },
                "replicate_metrics": {
                    str(replicate): {
                        "coverage": coverage[index],
                        "mean_contact_fraction": contact[index],
                    }
                    for index, replicate in enumerate((1, 2, 3))
                },
                "n_replicates": 3,
                "n_residues": len(rows),
                "total_frames_per_replicate": [1000, 1000, 1000],
                "criteria_cutoff": float(settings.cutoff),
                "residue_stats": rows,
                "residence_time_by_polymer_type": {},
            },
            metadata={
                "contacts_detection_fingerprint": contacts_detection_fingerprint(settings),
                "compute_residence_times": bool(settings.compute_residence_times),
                "equilibration": "10ns",
            },
            provenance={
                "source": "hypothesis_testing_consistency",
                "frame_selection": {"equilibration": "10ns"},
            },
        )

    aggregated = {
        "Control": artifact("Control", LOW_VARIANCE, (1e-6, 0.0, 1e-6)),
        "Treated": artifact("Treated", HIGH_VARIANCE, (0.0, 1e-6, 1e-6)),
    }
    return ContactsAnalysis(), _context(tmp_path, settings, aggregated, ttest_method)


CASES = {
    "rmsd": _rmsd_case,
    "rg": _rg_case,
    "sasa": _sasa_case,
    "distances": _distances_case,
    "contacts": _contacts_case,
}


def _pairwise_p_values(plugin: str, result: Any) -> list[tuple[float, float | None]]:
    """Return ``(p_value, p_value_adjusted)`` for every testable pairwise test."""
    pairs: list[tuple[float, float | None]] = []
    if plugin == "rmsd":
        for comparison in result.payload["comparisons"]:
            pairs.append((comparison["p_value"], comparison["p_adjusted"]))
    elif plugin == "distances":
        for comparison in result.pairwise_comparisons:
            pairs.append((comparison.distance_p_value, comparison.distance_p_value_adjusted))
            if comparison.fraction_p_value is not None:
                pairs.append((comparison.fraction_p_value, comparison.fraction_p_value_adjusted))
    elif plugin == "contacts":
        for comparison in result.pairwise_comparisons:
            for aggregate in comparison.aggregate_comparisons:
                if aggregate.metric != "mean_contact_fraction":
                    continue
                pairs.append((aggregate.p_value, aggregate.p_value_adjusted))
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


@pytest.mark.parametrize("plugin", ["distances", "rg", "rmsd", "sasa"])
def test_single_comparison_leaves_p_value_unchanged(plugin: str, tmp_path: Path) -> None:
    """A family of one test must have an adjusted p-value equal to the raw one."""
    analysis, ctx = CASES[plugin](tmp_path, "welch")
    pairs = _pairwise_p_values(plugin, analysis.compare(ctx))

    assert len(pairs) == 1
    raw_p, adjusted_p = pairs[0]
    assert adjusted_p == pytest.approx(raw_p)


@pytest.mark.parametrize("plugin", sorted(set(CASES) - {"rmsd"}))
def test_direction_labels_require_significance(plugin: str, tmp_path: Path) -> None:
    """Direction labels must not claim a change that the test did not find.

    The ported plugins are left out because the contract reports a signed delta
    and a significance flag rather than a direction word.
    """
    from polyzymd.analyses.shared.inferential_statistics import NO_SIGNIFICANT_CHANGE

    analysis, ctx = CASES[plugin](tmp_path, "welch")
    result = analysis.compare(ctx)

    if plugin == "distances":
        checked = [
            (c.distance_significant, c.distance_direction) for c in result.pairwise_comparisons
        ]
    elif plugin == "contacts":
        checked = [
            (aggregate.significant, aggregate.direction)
            for comparison in result.pairwise_comparisons
            for aggregate in comparison.aggregate_comparisons
        ]
    else:
        checked = [(c.significant, c.direction) for c in result.pairwise_comparisons]

    assert checked
    for significant, direction in checked:
        if not significant:
            assert direction == NO_SIGNIFICANT_CHANGE


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


def test_correction_family_spans_every_metric(monkeypatch: pytest.MonkeyPatch) -> None:
    """Three conditions and two metrics make one family of six tests."""
    from polyzymd.analyses.base import MetricValue
    from polyzymd.analyses.shared import inferential_statistics
    from polyzymd.analyses.stats import default_scalar_comparison

    family_sizes: list[int] = []
    real_benjamini_hochberg = inferential_statistics.benjamini_hochberg

    def _recording_benjamini_hochberg(p_values, alpha=0.05):
        family_sizes.append(len(p_values))
        return real_benjamini_hochberg(p_values, alpha=alpha)

    monkeypatch.setattr(inferential_statistics, "benjamini_hochberg", _recording_benjamini_hochberg)

    def _metrics(offset: float) -> dict[str, MetricValue]:
        first = [1.0 + offset, 1.1 + offset, 0.9 + offset]
        second = [5.0 + offset, 5.2 + offset, 4.8 + offset]
        return {
            "first": MetricValue(
                name="first",
                mean=float(np.mean(first)),
                sem=0.05,
                replicate_values=first,
            ),
            "second": MetricValue(
                name="second",
                mean=float(np.mean(second)),
                sem=0.05,
                replicate_values=second,
            ),
        }

    result = default_scalar_comparison(
        analysis_name="family",
        project_name="family",
        metrics_by_condition={
            "A": _metrics(0.0),
            "B": _metrics(1.0),
            "C": _metrics(2.0),
        },
        control_label=None,
    )

    # Three conditions give three pairs, and both metrics join the same family.
    assert len(result.pairwise_comparisons) == 6
    assert family_sizes == [6]
    assert all(comp.p_value_adjusted is not None for comp in result.pairwise_comparisons)

    # The omnibus ANOVA is outside the family and stays uncorrected.
    assert result.anova is not None
    assert len(result.anova) == 2
    assert all(anova.p_value_adjusted is None for anova in result.anova)
