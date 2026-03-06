"""Exposure dynamics and chaperone analysis module.

Provides tools to classify protein residue exposure stability from MD
trajectories, detect chaperone-like interaction events between polymer
and transiently exposed protein residues, and quantify chaperoning
effectiveness via refolding acceleration ratios and selectivity free
energies.

Key classes
-----------
ExposureConfig
    Configuration (thresholds, min event length).
ExposureDynamicsResult
    Per-residue exposure dynamics result (serialisable to JSON).
ResidueExposureSummary
    Per-residue summary: stability, event counts, polymer breakdown.
ChaperoneEventsResult
    Cached raw chaperone events with per-event frame-level data.
ChaperoneKineticsResult
    Refolding acceleration ratios rho(P, G) per (polymer, AA group).
ChaperoneSelectivityResult
    Chaperone selectivity free energy DeltaG_sel^chap(P, G).

Key functions
-------------
analyze_exposure_dynamics
    Top-level orchestrator: SASA + contacts → ExposureDynamicsResult.
compute_chaperone_kinetics
    Compute rho(P, G) from raw chaperone detection results.
compute_chaperone_selectivity
    Compute DeltaG_sel^chap(P, G) from events + SASA + contacts.
classify_residue_stability
    Classify a single residue as stably-exposed, stably-buried, or transient.
detect_events
    Detect chaperone and unassisted events for all residues.
"""

from polyzymd.analysis.exposure.chaperone import (
    ChaperoneDetectionResult,
    ChaperoneEvent,
    ChaperoneEventsResult,
    UnassistedRefoldingEvent,
    detect_events,
    detect_events_for_residue,
)
from polyzymd.analysis.exposure.chaperone_kinetics import (
    AccelerationRatioEntry,
    ChaperoneKineticsResult,
    compute_chaperone_kinetics,
)
from polyzymd.analysis.exposure.chaperone_selectivity import (
    ChaperoneSelectivityEntry,
    ChaperoneSelectivityResult,
    compute_chaperone_selectivity,
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
    "ChaperoneEventsResult",
    "detect_events",
    "detect_events_for_residue",
    # Dynamics
    "ResidueExposureSummary",
    "ExposureDynamicsResult",
    "analyze_exposure_dynamics",
    # Chaperone kinetics (rho)
    "AccelerationRatioEntry",
    "ChaperoneKineticsResult",
    "compute_chaperone_kinetics",
    # Chaperone selectivity (DeltaG_sel^chap)
    "ChaperoneSelectivityEntry",
    "ChaperoneSelectivityResult",
    "compute_chaperone_selectivity",
]
