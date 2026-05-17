"""Compatibility smoke tests for shared SASA helper imports."""

from __future__ import annotations

from polyzymd.analyses.sasa import _artifacts
from polyzymd.analyses.shared import sasa


def test_shared_sasa_facade_reexports_plugin_artifact_helpers() -> None:
    """Legacy shared imports should point at SASA plugin artifact helpers."""
    assert sasa.SASAComputationResult is _artifacts.SASAComputationResult
    assert (
        sasa.compute_sasa_artifact_compatibility_hash
        is _artifacts.compute_sasa_artifact_compatibility_hash
    )
    assert sasa.load_sasa_artifacts is _artifacts.load_sasa_artifacts
    assert sasa.save_sasa_artifacts is _artifacts.save_sasa_artifacts


def test_shared_sasa_facade_star_import_includes_compatibility_hash() -> None:
    """Legacy star imports should include the SASA compatibility hash helper."""

    namespace: dict[str, object] = {}

    exec("from polyzymd.analyses.shared.sasa import *", namespace)

    assert "compute_sasa_artifact_compatibility_hash" in sasa.__all__
    assert (
        namespace["compute_sasa_artifact_compatibility_hash"]
        is _artifacts.compute_sasa_artifact_compatibility_hash
    )
