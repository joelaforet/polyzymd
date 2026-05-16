"""Tests for MDAnalysis condition-artifact comparison."""

from __future__ import annotations

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
) -> ConditionArtifact:
    """Build a condition artifact with one AggregatedMetric payload."""

    replicate_ids = replicates or list(range(1, len(values) + 1))
    mean = sum(values) / len(values)
    payload: dict[str, object] = {
        "metrics": {
            metric_name: {
                "name": metric_name,
                "values": values,
                "mean": mean,
                "sem": 0.1,
                "std": 0.2,
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
    effective_control: str | None = "Control",
) -> MDAComparisonContext:
    """Build a deterministic MDA comparison context."""

    return MDAComparisonContext(
        analysis_name=analysis_name,
        project_name="project",
        expected_condition_labels=expected_condition_labels,
        control_label="Control",
        effective_control=effective_control,
        equilibration="10ns",
        settings_fingerprint="settings-1",
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
