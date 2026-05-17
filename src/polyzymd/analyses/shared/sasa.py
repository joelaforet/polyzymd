"""Compatibility facade for SASA plugin artifact helpers.

The implementation moved to :mod:`polyzymd.analyses.sasa._artifacts`. This
module remains temporarily available for existing imports.
"""

from polyzymd.analyses.sasa._artifacts import (
    NM2_TO_A2,
    SASA_ARTIFACT_COMPATIBILITY_VERSION,
    SASA_ARTIFACT_SCHEMA_NAME,
    SASA_ARTIFACT_SCHEMA_VERSION,
    SASA_COMPAT_PROBE_RADIUS_ABS_TOL,
    SASAArtifactCompatibility,
    SASAArtifactCompatibilityQuery,
    SASAArtifactContract,
    SASAComputationResult,
    SASASiblingArtifactMatch,
    build_sasa_artifact_contract,
    build_sasa_artifact_metadata,
    check_sasa_artifact_compatibility,
    compute_sasa,
    compute_sasa_artifact_compatibility_hash,
    find_sibling_sasa_artifacts,
    load_sasa_artifacts,
    resolve_selection_indices,
    save_sasa_artifacts,
    validate_target_subset,
)

__all__ = [
    "NM2_TO_A2",
    "SASA_ARTIFACT_COMPATIBILITY_VERSION",
    "SASA_ARTIFACT_SCHEMA_NAME",
    "SASA_ARTIFACT_SCHEMA_VERSION",
    "SASA_COMPAT_PROBE_RADIUS_ABS_TOL",
    "SASAArtifactCompatibility",
    "SASAArtifactCompatibilityQuery",
    "SASAArtifactContract",
    "SASAComputationResult",
    "SASASiblingArtifactMatch",
    "build_sasa_artifact_contract",
    "build_sasa_artifact_metadata",
    "check_sasa_artifact_compatibility",
    "compute_sasa",
    "compute_sasa_artifact_compatibility_hash",
    "find_sibling_sasa_artifacts",
    "load_sasa_artifacts",
    "resolve_selection_indices",
    "save_sasa_artifacts",
    "validate_target_subset",
]
