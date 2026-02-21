"""Exposure dynamics analysis module.

Provides tools to classify protein residue exposure stability from MD
trajectories and detect chaperone-like interaction events between polymer
and transiently exposed protein residues.

Key classes
-----------
ExposureConfig
    Configuration (thresholds, min event length).
ExposureDynamicsResult
    Per-residue exposure dynamics result (serialisable to JSON).
ResidueExposureSummary
    Per-residue summary: stability, event counts, polymer breakdown.
ChaperoneEnrichmentResult
    Dynamic enrichment of polymer types toward AA groups.

Key functions
-------------
analyze_exposure_dynamics
    Top-level orchestrator: SASA + contacts → ExposureDynamicsResult.
compute_chaperone_enrichment
    Compute enrichment scores for a single trajectory.
classify_residue_stability
    Classify a single residue as stably-exposed, stably-buried, or transient.
detect_events
    Detect chaperone and unassisted events for all residues.
"""

from polyzymd.analysis.exposure.chaperone import (
    ChaperoneDetectionResult,
    ChaperoneEvent,
    UnassistedRefoldingEvent,
    detect_events,
    detect_events_for_residue,
)
from polyzymd.analysis.exposure.classification import (
    classify_all_residues,
    classify_residue_stability,
)
from polyzymd.analysis.exposure.config import ExposureConfig
from polyzymd.analysis.exposure.dynamics import (
    ExposureDynamicsResult,
    ResidueExposureSummary,
    analyze_exposure_dynamics,
)
from polyzymd.analysis.exposure.enrichment import (
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
