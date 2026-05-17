"""Compatibility smoke tests for shared SASA helper imports."""

from __future__ import annotations

from polyzymd.analyses.sasa import _artifacts
from polyzymd.analyses.shared import sasa


def test_shared_sasa_facade_reexports_plugin_artifact_helpers() -> None:
    """Legacy shared imports should point at SASA plugin artifact helpers."""
    assert sasa.SASAComputationResult is _artifacts.SASAComputationResult
    assert sasa.load_sasa_artifacts is _artifacts.load_sasa_artifacts
    assert sasa.save_sasa_artifacts is _artifacts.save_sasa_artifacts
