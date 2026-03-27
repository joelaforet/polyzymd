"""Backward-compatibility re-export shim.

Canonical location: ``polyzymd.analyses._exposure_classification``

This module re-exports classification utilities from their new home so that
existing imports continue to work during the migration period.  Will be
removed in Phase 7.
"""

from polyzymd.analyses._exposure_classification import (
    ResidueStability,
    classify_all_residues,
    classify_residue_stability,
)

__all__ = [
    "ResidueStability",
    "classify_residue_stability",
    "classify_all_residues",
]
