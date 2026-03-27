"""Backward-compatibility shim — canonical location is ``polyzymd.analyses._contacts_surface_exposure``.

All public names are re-exported so that existing ``from polyzymd.analysis.contacts.surface_exposure
import …`` statements continue to work during the migration period.
"""

from polyzymd.analyses._contacts_surface_exposure import (  # noqa: F401
    ResidueExposure,
    SurfaceExposureFilter,
    SurfaceExposureResult,
)

__all__ = [
    "ResidueExposure",
    "SurfaceExposureFilter",
    "SurfaceExposureResult",
]
