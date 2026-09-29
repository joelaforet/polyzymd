"""Aggregation must refuse artifacts produced by the superseded estimator.

The settings fingerprint does not change when the correlation estimator
changes, so a re-run of the aggregate stage alone would average old numbers
with new ones. The estimator version key closes that gap and is checked here.
"""

from __future__ import annotations

import pytest

from polyzymd.analyses.mda import (
    MDAAggregationError,
    ReplicateArtifact,
    validate_autocorrelation_estimator_version,
)
from polyzymd.analyses.shared.autocorrelation import AUTOCORRELATION_ESTIMATOR_VERSION


def _replicate_artifact(*, estimator_version: str | None) -> ReplicateArtifact:
    """Return a replicate artifact stamped with a given estimator version."""

    metadata: dict[str, object] = {"settings_fingerprint": "fingerprint"}
    if estimator_version is not None:
        metadata["autocorrelation_estimator_version"] = estimator_version
    return ReplicateArtifact(
        analysis_name="demo",
        condition_label="Control",
        replicate=1,
        payload={},
        metadata=metadata,
    )


def test_aggregation_rejects_missing_estimator_version() -> None:
    """An artifact from before the fix carries no estimator version."""

    artifact = _replicate_artifact(estimator_version=None)

    with pytest.raises(MDAAggregationError, match="no autocorrelation estimator version"):
        validate_autocorrelation_estimator_version(artifact, analysis_label="Demo")


def test_aggregation_rejects_superseded_estimator_version() -> None:
    """An artifact stamped with the old estimator must not be aggregated."""

    artifact = _replicate_artifact(estimator_version="1")

    with pytest.raises(MDAAggregationError, match="clear stale caches"):
        validate_autocorrelation_estimator_version(artifact, analysis_label="Demo")


def test_aggregation_accepts_current_estimator_version() -> None:
    """The current version passes the estimator gate."""

    artifact = _replicate_artifact(estimator_version=AUTOCORRELATION_ESTIMATOR_VERSION)

    validate_autocorrelation_estimator_version(artifact, analysis_label="Demo")
