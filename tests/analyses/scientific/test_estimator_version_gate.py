"""Aggregation must refuse artifacts produced by the superseded estimator.

The settings fingerprint does not change when the correlation estimator
changes, and neither does it change when RMSF stops thinning frames, so a
re-run of the aggregate stage alone would happily average old numbers with
new ones. Two version keys close that gap and are checked here.
"""

from __future__ import annotations

import pytest

from polyzymd.analyses.catalytic_triad import CatalyticTriadSettings, TriadPairSettings
from polyzymd.analyses.distances import DistancesSettings
from polyzymd.analyses.mda import MDAAggregationError, ReplicateArtifact
from polyzymd.analyses.rmsf._mda import RMSF_PROFILE_VERSION
from polyzymd.analyses.rmsf._mda import (
    _validate_and_order_artifacts as _validate_rmsf_artifacts,
)
from polyzymd.analyses.sasa import SASASettings
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
    if analysis_name == "rmsd":
        from polyzymd.analyses.rmsd._mda import _validate_and_order_artifacts as validate

        return validate(run_labels=["default"], **common)
    if analysis_name == "sasa":
        from polyzymd.analyses.sasa._mda import _validate_and_order_artifacts as validate

        return validate(settings=SASASettings(), **common)
    if analysis_name == "distances":
        from polyzymd.analyses.distances._mda import _validate_and_order_artifacts as validate

        return validate(settings=DistancesSettings(), analysis_dir=tmp_path, **common)
    if analysis_name == "catalytic_triad":
        from polyzymd.analyses.catalytic_triad._mda import (
            _validate_and_order_artifacts as validate,
        )

        return validate(settings=_triad_settings(), analysis_dir=tmp_path, **common)
    raise AssertionError(f"unhandled analysis {analysis_name}")


PLUGINS = ["rmsd", "sasa", "distances", "catalytic_triad"]


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


def _rmsf_artifact(profile_version: str | None) -> ReplicateArtifact:
    """Return an RMSF replicate artifact at a given profile version."""

    metadata: dict[str, object] = {
        "settings_fingerprint": FINGERPRINT,
        "selection_string": "protein and name CA",
    }
    if profile_version is not None:
        metadata["rmsf_profile_version"] = profile_version
    return ReplicateArtifact(
        analysis_name="rmsf",
        condition_label="Control",
        replicate=1,
        payload={},
        metadata=metadata,
    )


@pytest.mark.parametrize("stored_version", [None, "1"])
def test_rmsf_aggregation_rejects_old_profile_version(stored_version, tmp_path) -> None:
    """Version 1 profiles were subsampled and must not be averaged with version 2."""

    from polyzymd.analyses.rmsf import RMSFSettings

    with pytest.raises(MDAAggregationError, match="clear stale caches"):
        _validate_rmsf_artifacts(
            condition_label="Control",
            expected_replicates=[1],
            settings=RMSFSettings(),
            settings_fingerprint=FINGERPRINT,
            artifacts=[_rmsf_artifact(stored_version)],
            analysis_dir=tmp_path,
        )


def test_rmsf_aggregation_rejects_mixed_profile_versions(tmp_path) -> None:
    """One stale replicate is enough to stop the whole condition."""

    from polyzymd.analyses.rmsf import RMSFSettings

    current = _rmsf_artifact(RMSF_PROFILE_VERSION)
    stale = _rmsf_artifact("1")
    stale = stale.model_copy(update={"replicate": 2})

    with pytest.raises(MDAAggregationError, match="RMSF replicate 2"):
        _validate_rmsf_artifacts(
            condition_label="Control",
            expected_replicates=[1, 2],
            settings=RMSFSettings(),
            settings_fingerprint=FINGERPRINT,
            artifacts=[stale, current],
            analysis_dir=tmp_path,
        )
