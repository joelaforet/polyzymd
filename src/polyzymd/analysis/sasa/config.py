"""Backward-compatibility re-export shim.

Canonical location: ``polyzymd.analyses._sasa_config``

This module re-exports ``SASAConfig`` from its new home so that existing
``from polyzymd.analysis.sasa.config import SASAConfig`` statements continue
to work during the migration period.  Will be removed in Phase 7.
"""

from polyzymd.analyses._sasa_config import SASAConfig

__all__ = ["SASAConfig"]
