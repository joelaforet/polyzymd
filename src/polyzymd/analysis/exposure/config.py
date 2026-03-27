"""Backward-compatibility re-export shim.

Canonical location: ``polyzymd.analyses._exposure_config``

This module re-exports ``ExposureConfig`` from its new home so that existing
``from polyzymd.analysis.exposure.config import ExposureConfig`` statements
continue to work during the migration period.  Will be removed in Phase 7.
"""

from polyzymd.analyses._exposure_config import ExposureConfig

__all__ = ["ExposureConfig"]
