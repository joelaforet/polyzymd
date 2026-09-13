"""Aggregation must refuse artifacts produced by the superseded estimator.

The settings fingerprint does not change when the correlation estimator
changes, so a re-run of the aggregate stage alone would happily average old
numbers with new ones. The estimator version key closes that gap and is checked
here. Contract plugins get the same protection from the ``plugin_code_hash``
field of the identity block, which ``tests/analyses/test_contract.py`` covers.
"""

from __future__ import annotations

import pytest

from polyzymd.analyses.catalytic_triad import CatalyticTriadSettings, TriadPairSettings
from polyzymd.analyses.distances import DistancesSettings
from polyzymd.analyses.mda import MDAAggregationError, ReplicateArtifact
from polyzymd.analyses.shared.autocorrelation import AUTOCORRELATION_ESTIMATOR_VERSION

FINGERPRINT = "fingerprint"


def _triad_settings() -> CatalyticTriadSettings:
    """Return minimal catalytic-triad settings."""

    return CatalyticTriadSettings(
        name="triad",
        pairs=[
            TriadPairSettings(
                label="A-B",
                selection_a="resid 1 and name OD1",
                selection_b="resid 2 and name ND1",
            )
        ],
    )


def _replicate_artifact(analysis_name: str, *, estimator_version: str | None) -> ReplicateArtifact:
    """Return a replicate artifact stamped with a given estimator version."""

    metadata: dict[str, object] = {"settings_fingerprint": FINGERPRINT}
    if estimator_version is not None:
        metadata["autocorrelation_estimator_version"] = estimator_version
    return ReplicateArtifact(
        analysis_name=analysis_name,
        condition_label="Control",
        replicate=1,
        payload={},
        metadata=metadata,
    )


def _validator_call(analysis_name: str, artifact: ReplicateArtifact, tmp_path):
    """Call one plugin's replicate-artifact validator with the given artifact."""

    common = {
        "condition_label": "Control",
        "expected_replicates": [1],
        "settings_fingerprint": FINGERPRINT,
        "artifacts": [artifact],
    }
    if analysis_name == "rg":
        from polyzymd.analyses.rg._mda import _validate_and_order_artifacts as validate

        return validate(run_labels=["default"], analysis_dir=tmp_path, **common)
    if analysis_name == "distances":
        from polyzymd.analyses.distances._mda import _validate_and_order_artifacts as validate

        return validate(settings=DistancesSettings(), analysis_dir=tmp_path, **common)
    if analysis_name == "catalytic_triad":
        from polyzymd.analyses.catalytic_triad._mda import (
            _validate_and_order_artifacts as validate,
        )

        return validate(settings=_triad_settings(), analysis_dir=tmp_path, **common)
    raise AssertionError(f"unhandled analysis {analysis_name}")


PLUGINS = ["rg", "distances", "catalytic_triad"]


@pytest.mark.parametrize("analysis_name", PLUGINS)
def test_aggregation_rejects_missing_estimator_version(analysis_name: str, tmp_path) -> None:
    """An artifact from before the fix carries no estimator version."""

    artifact = _replicate_artifact(analysis_name, estimator_version=None)

    with pytest.raises(MDAAggregationError, match="no autocorrelation estimator version"):
        _validator_call(analysis_name, artifact, tmp_path)


@pytest.mark.parametrize("analysis_name", PLUGINS)
def test_aggregation_rejects_superseded_estimator_version(analysis_name: str, tmp_path) -> None:
    """An artifact stamped with the old estimator must not be aggregated."""

    artifact = _replicate_artifact(analysis_name, estimator_version="1")

    with pytest.raises(MDAAggregationError, match="clear stale caches"):
        _validator_call(analysis_name, artifact, tmp_path)


@pytest.mark.parametrize("analysis_name", PLUGINS)
def test_aggregation_accepts_current_estimator_version(analysis_name: str, tmp_path) -> None:
    """The current version passes the estimator gate, whatever happens after it."""

    artifact = _replicate_artifact(
        analysis_name, estimator_version=AUTOCORRELATION_ESTIMATOR_VERSION
    )

    try:
        _validator_call(analysis_name, artifact, tmp_path)
    except Exception as error:  # payload checks downstream of the gate may still fire
        assert "estimator version" not in str(error)
