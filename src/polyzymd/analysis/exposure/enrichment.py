"""Backward-compatibility re-export shim.

Canonical location: ``polyzymd.analyses._exposure_enrichment``

This module re-exports enrichment classes and functions from their new home
so that existing imports continue to work during the migration period.
Will be removed in Phase 7.
"""

from polyzymd.analyses._exposure_enrichment import (
    ChaperoneEnrichmentResult,
    GroupEnrichmentEntry,
    compute_chaperone_enrichment,
)

__all__ = [
    "GroupEnrichmentEntry",
    "ChaperoneEnrichmentResult",
    "compute_chaperone_enrichment",
]
