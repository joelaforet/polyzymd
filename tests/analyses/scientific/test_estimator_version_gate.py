"""Aggregation must refuse artifacts produced by the superseded estimator.

The settings fingerprint does not change when the correlation estimator
changes, so a re-run of the aggregate stage alone would happily average old
numbers with new ones. The estimator version key closes that gap and is checked
here. Every plugin that used to call this gate from its own aggregation helper
now runs on the observable contract, where ``plugin_code_hash`` gives the same
protection, so the gate itself is tested directly.
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
        analysis_name="probe",
        condition_label="Control",
        replicate=1,
        payload={},
        metadata=metadata,
    )


def test_aggregation_rejects_missing_estimator_version() -> None:
    """An artifact from before the fix carries no estimator version."""

    with pytest.raises(MDAAggregationError, match="no autocorrelation estimator version"):
        validate_autocorrelation_estimator_version(
            _replicate_artifact(estimator_version=None), analysis_label="probe"
        )


def test_aggregation_rejects_superseded_estimator_version() -> None:
    """An artifact stamped with the old estimator must not be aggregated."""

    with pytest.raises(MDAAggregationError, match="clear stale caches"):
        validate_autocorrelation_estimator_version(
            _replicate_artifact(estimator_version="1"), analysis_label="probe"
        )


def test_aggregation_accepts_current_estimator_version() -> None:
    """The current version passes the estimator gate."""

    validate_autocorrelation_estimator_version(
        _replicate_artifact(estimator_version=AUTOCORRELATION_ESTIMATOR_VERSION),
        analysis_label="probe",
    )
