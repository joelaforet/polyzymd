"""Backward-compatibility re-export shim.

Canonical location: ``polyzymd.analyses._exposure_chaperone``

This module re-exports chaperone event detection classes and functions from
their new home so that existing imports continue to work during the migration
period.  Will be removed in Phase 7.
"""

from polyzymd.analyses._exposure_chaperone import (
    ChaperoneDetectionResult,
    ChaperoneEvent,
    UnassistedRefoldingEvent,
    detect_events,
    detect_events_for_residue,
)

__all__ = [
    "ChaperoneEvent",
    "UnassistedRefoldingEvent",
    "ChaperoneDetectionResult",
    "detect_events_for_residue",
    "detect_events",
]
