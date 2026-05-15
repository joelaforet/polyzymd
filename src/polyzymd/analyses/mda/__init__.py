"""Public MDAnalysis extension-layer primitives.

This package exposes import-light protocol and typing helpers for integrations
that wrap MDAnalysis ``AnalysisBase`` objects without importing MDAnalysis at
module import time.
"""

from polyzymd.analyses.mda.base import (
    MDA_EXTENSION_API_VERSION,
    AnalysisBaseLike,
    MDAnalysisExtensionError,
    MDARunKwargs,
)
from polyzymd.analyses.mda.frame_selection import FrameSelection
from polyzymd.analyses.mda.universe import FileIdentity, UniverseProvenance, UniverseProvider

__all__ = [
    "MDA_EXTENSION_API_VERSION",
    "AnalysisBaseLike",
    "MDAnalysisExtensionError",
    "MDARunKwargs",
    "FrameSelection",
    "FileIdentity",
    "UniverseProvider",
    "UniverseProvenance",
]
