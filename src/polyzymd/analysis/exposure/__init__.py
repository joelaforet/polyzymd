"""Backward-compatibility re-export shim for exposure analysis module.

Canonical locations are now under ``polyzymd.analyses._exposure_*``.
This package re-exports all public names so that existing
``from polyzymd.analysis.exposure import ...`` statements continue to work.
Will be removed in Phase 7.
"""

from polyzymd.analyses._exposure_chaperone import (
    ChaperoneDetectionResult,
    ChaperoneEvent,
    UnassistedRefoldingEvent,
    detect_events,
    detect_events_for_residue,
)
from polyzymd.analyses._exposure_classification import (
    classify_all_residues,
    classify_residue_stability,
)
from polyzymd.analyses._exposure_config import ExposureConfig
from polyzymd.analyses._exposure_dynamics import (
    ExposureDynamicsResult,
    ResidueExposureSummary,
    analyze_exposure_dynamics,
)
from polyzymd.analyses._exposure_enrichment import (
    ChaperoneEnrichmentResult,
    GroupEnrichmentEntry,
    compute_chaperone_enrichment,
)

__all__ = [
    # Config
    "ExposureConfig",
    # Classification
    "classify_residue_stability",
    "classify_all_residues",
    # Chaperone event detection
    "ChaperoneEvent",
    "UnassistedRefoldingEvent",
    "ChaperoneDetectionResult",
    "detect_events",
    "detect_events_for_residue",
    # Dynamics
    "ResidueExposureSummary",
    "ExposureDynamicsResult",
    "analyze_exposure_dynamics",
    # Enrichment
    "GroupEnrichmentEntry",
    "ChaperoneEnrichmentResult",
    "compute_chaperone_enrichment",
]
