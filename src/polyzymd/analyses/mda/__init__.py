"""Public MDAnalysis extension-layer primitives.

This package exposes import-light protocol and typing helpers for integrations
that wrap MDAnalysis ``AnalysisBase`` objects without importing MDAnalysis at
module import time.
"""

from polyzymd.analyses.mda.artifacts import (
    MDA_ARTIFACT_SCHEMA_VERSION,
    ArtifactEnvelope,
    ArtifactManifest,
    ArtifactSidecarRef,
    ComparisonArtifact,
    ConditionArtifact,
    ReplicateArtifact,
)
from polyzymd.analyses.mda.base import (
    MDA_EXTENSION_API_VERSION,
    AnalysisBaseLike,
    MDAnalysisExtensionError,
    MDARunKwargs,
)
from polyzymd.analyses.mda.frame_selection import FrameSelection, frame_selection_payload
from polyzymd.analyses.mda.pair_distance import PairSelection, pair_distance_matrix
from polyzymd.analyses.mda.store import ArtifactStore, ArtifactStoreError
from polyzymd.analyses.mda.universe import FileIdentity, UniverseProvenance, UniverseProvider

__all__ = [
    "MDA_EXTENSION_API_VERSION",
    "MDA_ARTIFACT_SCHEMA_VERSION",
    "AnalysisBaseLike",
    "MDAnalysisExtensionError",
    "MDARunKwargs",
    "ArtifactEnvelope",
    "ArtifactManifest",
    "ArtifactSidecarRef",
    "ComparisonArtifact",
    "ConditionArtifact",
    "ReplicateArtifact",
    "ArtifactStore",
    "ArtifactStoreError",
    "FrameSelection",
    "PairSelection",
    "pair_distance_matrix",
    "frame_selection_payload",
    "FileIdentity",
    "UniverseProvider",
    "UniverseProvenance",
]
