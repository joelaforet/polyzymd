"""Tests for MDAnalysis condition-artifact comparison."""

from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar

import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import Analysis, ComparisonContext, Condition, MetricValue
from polyzymd.analyses.mda.artifacts import ComparisonArtifact, ConditionArtifact
from polyzymd.analyses.mda.comparison import (
    MDAComparisonContext,
    MDAComparisonError,
    compare_condition_artifacts,
)


class _Settings(BaseModel):
    """Small settings stand-in for comparison context tests."""

    scale: float = 1.0


class _ArtifactDefaultCompareAnalysis(Analysis):
    """Analysis using default comparison dispatch for tests."""

    name: ClassVar[str] = "artifact_default_compare"
    Settings: ClassVar[type] = _Settings
    has_compute_stage: ClassVar[bool] = False
    has_aggregate_stage: ClassVar[bool] = False

    def aggregate_settings_fingerprint(self, settings: Any) -> str | None:
        """Return a deterministic settings fingerprint for artifacts."""

        del settings
        return "settings-1"

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract legacy metrics and fail if MDA dispatch should have been used."""

        if isinstance(summary, ConditionArtifact):
            raise AssertionError("MDA condition artifacts must not call extract_metrics")
        return {
            "mean_value": MetricValue(
                name="mean_value",
                mean=summary["mean"],
                sem=summary["sem"],
                replicate_values=list(summary["values"]),
            )
        }


def _condition(label: str, tmp_path: Path, replicates: tuple[int, ...] = (1, 2, 3)) -> Condition:
    """Build a synthetic condition for comparison tests."""

    return Condition(label, tmp_path / f"{label}.yaml", replicates, SimpleNamespace())


def _artifact(
    label: str,
    values: list[float],
    *,
    analysis_name: str = "rmsd",
    metric_name: str = "mean_rmsd",
    replicates: list[int] | None = None,
    settings_fingerprint: str = "settings-1",
    metric_metadata: dict[str, object] | None = None,
    skipped_replicates: list[dict[str, object]] | None = None,
    stored_mean: float | None = None,
    stored_std: float | None = None,
    stored_sem: float | None = None,
) -> ConditionArtifact:
    """Build a condition artifact with one AggregatedMetric payload."""

    replicate_ids = replicates if replicates is not None else list(range(1, len(values) + 1))
    mean = sum(values) / len(values)
    if len(values) == 1:
        std = 0.0
    else:
        std = math.sqrt(sum((value - mean) ** 2 for value in values) / (len(values) - 1))
    sem = std / math.sqrt(len(values))
    payload: dict[str, object] = {
        "metrics": {
            metric_name: {
                "name": metric_name,
                "values": values,
                "mean": mean if stored_mean is None else stored_mean,
                "sem": sem if stored_sem is None else stored_sem,
                "std": std if stored_std is None else stored_std,
                "n": len(values),
            }
        }
    }
    if metric_metadata is not None:
        payload["metric_metadata"] = {metric_name: metric_metadata}
    return ConditionArtifact(
        analysis_name=analysis_name,
        condition_label=label,
        replicates=replicate_ids,
        payload=payload,
        metadata={"settings_fingerprint": settings_fingerprint},
        source_replicates=[{"replicate": replicate} for replicate in replicate_ids],
        skipped_replicates=skipped_replicates or [],
        warnings=[f"warning for {label}"],
    )


def _context(
    *,
    analysis_name: str = "rmsd",
    expected_condition_labels: tuple[str, ...] = ("Control", "PEG"),
    expected_replicates_by_condition: dict[str, tuple[int, ...]] | None = None,
    effective_control: str | None = "Control",
    min_replicates: int = 1,
) -> MDAComparisonContext:
    """Build a deterministic MDA comparison context."""

    return MDAComparisonContext(
        analysis_name=analysis_name,
        project_name="project",
        expected_condition_labels=expected_condition_labels,
        expected_replicates_by_condition=expected_replicates_by_condition,
        control_label="Control",
        effective_control=effective_control,
        equilibration="10ns",
        settings_fingerprint="settings-1",
        min_replicates=min_replicates,
        fdr_alpha=0.1,
        ttest_method="welch",
        posthoc_method="ttest_bh",
    )


def test_compare_condition_artifacts_returns_stable_comparison_artifact() -> None:
    """Condition artifacts should feed the shared scalar comparison pipeline."""

    artifacts = [
        _artifact(
            "Control",
            [1.0, 1.1, 0.9],
            metric_metadata={"higher_is_better": False, "unit": "nm"},
        ),
        _artifact(
            "PEG",
            [1.4, 1.5, 1.6],
            metric_metadata={"higher_is_better": False, "unit": "nm"},
        ),
    ]

    result = compare_condition_artifacts(artifacts, _context())

    assert isinstance(result, ComparisonArtifact)
    assert result.analysis_name == "rmsd"
    assert result.conditions == ["Control", "PEG"]
    assert result.control_label == "Control"
    assert result.effective_control == "Control"
    assert result.payload["statistical_parameters"] == {
        "fdr_alpha": 0.1,
        "ttest_method": "welch",
        "posthoc_method": "ttest_bh",
        "control_label": "Control",
        "effective_control": "Control",
        "equilibration": "10ns",
    }
    assert result.payload["condition_summaries"][0]["mean_rmsd_mean"] == pytest.approx(1.0)
    assert result.payload["pairwise_comparisons"][0]["metric"] == "mean_rmsd"
    assert result.payload["ranking"] == ["Control", "PEG"]
    assert result.payload["metric_metadata"]["mean_rmsd"]["higher_is_better"] is False
    assert result.payload["metric_metadata"]["mean_rmsd"]["unit"] == "nm"
    assert result.provenance["source_condition_labels"] == ["Control", "PEG"]
    assert result.provenance["source_replicates"] == {"Control": [1, 2, 3], "PEG": [1, 2, 3]}
    assert result.warnings == ["warning for Control", "warning for PEG"]


def test_compare_condition_artifacts_uses_shared_statistics_pipeline() -> None:
    """MDA comparison should preserve t-test, FDR, ANOVA, and effect-size behavior."""

    from polyzymd.analyses.shared.inferential_statistics import (
        benjamini_hochberg,
        cohens_d,
        independent_ttest,
        one_way_anova,
        percent_change,
    )

    values_by_label = {
        "Control": [1.0, 2.0, 3.0],
        "PEG": [2.0, 4.0, 8.0],
        "SBMA": [5.0, 6.0, 7.0],
    }
    artifacts = [
        _artifact(label, values, metric_metadata={"higher_is_better": True})
        for label, values in values_by_label.items()
    ]
    ctx = MDAComparisonContext(
        analysis_name="rmsd",
        project_name="project",
        expected_condition_labels=("Control", "PEG", "SBMA"),
        control_label=None,
        effective_control=None,
        equilibration="10ns",
        settings_fingerprint="settings-1",
        fdr_alpha=0.2,
        ttest_method="welch",
        posthoc_method="ttest_bh",
    )

    result = compare_condition_artifacts(artifacts, ctx)
    pairwise = result.payload["pairwise_comparisons"]
    pair_by_labels = {(entry["condition_a"], entry["condition_b"]): entry for entry in pairwise}

    control_summary = result.payload["condition_summaries"][0]
    assert control_summary["mean_rmsd_mean"] == pytest.approx(2.0)
    assert control_summary["mean_rmsd_sem"] == pytest.approx(1.0 / math.sqrt(3.0))
    assert control_summary["mean_rmsd_replicate_values"] == values_by_label["Control"]

    control_vs_peg = pair_by_labels[("Control", "PEG")]
    expected_welch = independent_ttest(
        values_by_label["Control"],
        values_by_label["PEG"],
        method="welch",
    )
    expected_effect = cohens_d(values_by_label["Control"], values_by_label["PEG"])
    assert control_vs_peg["posthoc_method"] == "ttest_bh"
    assert control_vs_peg["t_statistic"] == pytest.approx(expected_welch.t_statistic)
    assert control_vs_peg["p_value"] == pytest.approx(expected_welch.p_value)
    assert control_vs_peg["cohens_d"] == pytest.approx(expected_effect.cohens_d)
    assert control_vs_peg["effect_size_interpretation"] == expected_effect.interpretation
    assert control_vs_peg["percent_change"] == pytest.approx(
        percent_change(
            sum(values_by_label["Control"]) / 3.0,
            sum(values_by_label["PEG"]) / 3.0,
        )
    )

    expected_bh = benjamini_hochberg([entry["p_value"] for entry in pairwise], alpha=0.2)
    for entry, expected in zip(pairwise, expected_bh, strict=True):
        assert entry["p_value_adjusted"] == pytest.approx(expected.adjusted_p_value)
        assert entry["significant"] is expected.significant

    expected_anova = one_way_anova(*values_by_label.values())
    assert result.payload["anova"][0]["metric"] == "mean_rmsd"
    assert result.payload["anova"][0]["f_statistic"] == pytest.approx(expected_anova.f_statistic)
    assert result.payload["anova"][0]["p_value"] == pytest.approx(expected_anova.p_value)
    assert result.payload["ranking"] == ["SBMA", "PEG", "Control"]


@pytest.mark.parametrize(
    ("override", "match"),
    [
        ({"stored_mean": 99.0}, "stored mean"),
        ({"stored_std": 99.0}, "stored std"),
        ({"stored_sem": 99.0}, "stored sem"),
    ],
)
def test_compare_condition_artifacts_rejects_stale_aggregated_statistics(
    override: dict[str, float],
    match: str,
) -> None:
    """Stored aggregate statistics must agree with replicate-level values."""

    artifacts = [
        _artifact("Control", [1.0, 2.0, 3.0], **override),
        _artifact("PEG", [2.0, 3.0, 4.0]),
    ]

    with pytest.raises(MDAComparisonError, match=match):
        compare_condition_artifacts(artifacts, _context())


@pytest.mark.parametrize(
    ("artifact", "ctx", "match"),
    [
        (
            _artifact("Control", [1.0], replicates=[]),
            _context(expected_replicates_by_condition={"Control": (1, 2), "PEG": (1, 2)}),
            "no replicates",
        ),
        (
            _artifact("Control", [1.0], replicates=[1]),
            _context(
                expected_replicates_by_condition={"Control": (1, 2), "PEG": (1, 2)},
                min_replicates=2,
            ),
            "below required minimum 2",
        ),
        (
            _artifact("Control", [1.0, 2.0], replicates=[2, 4]),
            _context(expected_replicates_by_condition={"Control": (1, 2, 3), "PEG": (1, 2)}),
            r"unexpected replicate IDs \[4\]",
        ),
    ],
)
def test_compare_condition_artifacts_rejects_bad_replicate_identity(
    artifact: ConditionArtifact,
    ctx: MDAComparisonContext,
    match: str,
) -> None:
    """MDA comparison should reject stale or insufficient replicate sets."""

    artifacts = [artifact, _artifact("PEG", [2.0, 3.0], replicates=[1, 2])]

    with pytest.raises(MDAComparisonError, match=match):
        compare_condition_artifacts(artifacts, ctx)


def test_compare_condition_artifacts_accepts_expected_replicate_subset() -> None:
    """MDA comparison should accept valid partial aggregate replicate subsets."""

    artifacts = [
        _artifact("Control", [1.0, 3.0], replicates=[1, 3]),
        _artifact("PEG", [2.0, 4.0], replicates=[1, 2]),
    ]
    ctx = _context(
        expected_replicates_by_condition={"Control": (1, 2, 3), "PEG": (1, 2, 3)},
        min_replicates=2,
    )

    result = compare_condition_artifacts(artifacts, ctx)

    assert result.provenance["source_replicates"] == {"Control": [1, 3], "PEG": [1, 2]}


@pytest.mark.parametrize(
    ("artifacts", "match"),
    [
        ([_artifact("Control", [1.0]), _artifact("Control", [1.1])], "duplicate"),
        (
            [_artifact("Control", [1.0], analysis_name="rg"), _artifact("PEG", [1.1])],
            "analysis mismatch",
        ),
        (
            [_artifact("Control", [1.0], settings_fingerprint="old"), _artifact("PEG", [1.1])],
            "settings fingerprint mismatch",
        ),
        (
            [
                _artifact("Control", [1.0], metric_name="a"),
                _artifact("PEG", [1.1], metric_name="b"),
            ],
            "inconsistent metric keys",
        ),
        ([_artifact("Control", [float("nan")]), _artifact("PEG", [1.1])], "non-finite scalar"),
        (
            [_artifact("Control", [1.0], replicates=[1, 2]), _artifact("PEG", [1.1])],
            r"n == len\(values\) == len\(replicates\)",
        ),
    ],
)
def test_compare_condition_artifacts_rejects_invalid_inputs(
    artifacts: list[ConditionArtifact],
    match: str,
) -> None:
    """Comparison should reject stale, mismatched, and invalid aggregate artifacts."""

    with pytest.raises(MDAComparisonError, match=match):
        compare_condition_artifacts(artifacts, _context(expected_condition_labels=None))


def test_compare_condition_artifacts_validates_expected_condition_labels() -> None:
    """Comparison should fail when required condition artifacts are missing."""

    with pytest.raises(MDAComparisonError, match=r"missing=\['PEG'\]"):
        compare_condition_artifacts([_artifact("Control", [1.0])], _context())


def test_default_compare_dispatches_all_mda_artifacts(tmp_path: Path) -> None:
    """Default Analysis.compare should use MDA artifact comparison when all summaries match."""

    analysis = _ArtifactDefaultCompareAnalysis()
    conditions = [_condition("Control", tmp_path), _condition("PEG", tmp_path)]
    ctx = ComparisonContext(
        name="project",
        conditions=conditions,
        excluded_conditions=[],
        control_label="Control",
        analysis_dirs={label: tmp_path / label for label in ("Control", "PEG")},
        results_dir=tmp_path / "comparison",
        equilibration="10ns",
        settings=_Settings(),
        fdr_alpha=0.2,
        ttest_method="welch",
        posthoc_method="ttest_bh",
        aggregated_results={
            "Control": _artifact("Control", [1.0, 1.1, 0.9], analysis_name=analysis.name),
            "PEG": _artifact("PEG", [1.4, 1.5, 1.6], analysis_name=analysis.name),
        },
    )

    result = analysis.compare(ctx)

    assert isinstance(result, ComparisonArtifact)
    assert result.payload["statistical_parameters"]["fdr_alpha"] == 0.2
    assert result.payload["statistical_parameters"]["ttest_method"] == "welch"


def test_default_compare_rejects_mda_artifact_replicates_outside_active_condition(
    tmp_path: Path,
) -> None:
    """Default MDA compare should validate artifacts against active condition replicates."""

    analysis = _ArtifactDefaultCompareAnalysis()
    conditions = [_condition("Control", tmp_path), _condition("PEG", tmp_path)]
    ctx = ComparisonContext(
        name="project",
        conditions=conditions,
        excluded_conditions=[],
        control_label="Control",
        analysis_dirs={label: tmp_path / label for label in ("Control", "PEG")},
        results_dir=tmp_path / "comparison",
        equilibration="10ns",
        settings=_Settings(),
        aggregated_results={
            "Control": _artifact(
                "Control",
                [1.0, 1.1, 0.9],
                analysis_name=analysis.name,
                replicates=[4, 5, 6],
            ),
            "PEG": _artifact("PEG", [1.4, 1.5, 1.6], analysis_name=analysis.name),
        },
    )

    with pytest.raises(MDAComparisonError, match=r"active replicates are \[1, 2, 3\]"):
        analysis.compare(ctx)


def test_default_compare_keeps_legacy_extract_metrics_path(tmp_path: Path) -> None:
    """Default Analysis.compare should keep legacy summaries on the existing path."""

    analysis = _ArtifactDefaultCompareAnalysis()
    conditions = [_condition("Control", tmp_path), _condition("PEG", tmp_path)]
    ctx = ComparisonContext(
        name="project",
        conditions=conditions,
        excluded_conditions=[],
        control_label="Control",
        analysis_dirs={label: tmp_path / label for label in ("Control", "PEG")},
        results_dir=tmp_path / "comparison",
        equilibration="10ns",
        settings=_Settings(),
        aggregated_results={
            "Control": {
                "mean": 1.0,
                "sem": 0.1,
                "values": [0.9, 1.0, 1.1],
                "settings_fingerprint": "settings-1",
                "replicates": [1, 2, 3],
                "n_replicates": 3,
            },
            "PEG": {
                "mean": 1.5,
                "sem": 0.1,
                "values": [1.4, 1.5, 1.6],
                "settings_fingerprint": "settings-1",
                "replicates": [1, 2, 3],
                "n_replicates": 3,
            },
        },
    )

    result = analysis.compare(ctx)

    assert result is not None
    assert result.analysis_type == analysis.name
    assert result.conditions[0].label == "Control"


def test_default_compare_rejects_mixed_legacy_and_mda_summaries(tmp_path: Path) -> None:
    """Default compare should not silently mix old aggregate JSON with MDA artifacts."""

    analysis = _ArtifactDefaultCompareAnalysis()
    conditions = [_condition("Control", tmp_path), _condition("PEG", tmp_path)]
    ctx = ComparisonContext(
        name="project",
        conditions=conditions,
        excluded_conditions=[],
        control_label="Control",
        analysis_dirs={label: tmp_path / label for label in ("Control", "PEG")},
        results_dir=tmp_path / "comparison",
        equilibration="10ns",
        settings=_Settings(),
        aggregated_results={
            "Control": _artifact("Control", [1.0], analysis_name=analysis.name),
            "PEG": {"mean": 1.5, "sem": 0.1, "values": [1.5]},
        },
    )

    with pytest.raises(MDAComparisonError, match="clear stale aggregate result.json"):
        analysis.compare(ctx)
