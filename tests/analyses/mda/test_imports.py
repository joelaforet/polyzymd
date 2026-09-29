"""Tests for MDAnalysis extension-layer public imports."""

from __future__ import annotations

import importlib
import sys
from typing import Any


def test_public_facade_reexports_primitives() -> None:
    """The package facade should expose the stable extension-layer primitives."""

    from polyzymd.analyses import mda
    from polyzymd.analyses.mda.aggregation import (
        AggregatedMetric,
        ExplicitReplicateMetricPolicy,
        MDAAggregationContext,
        MDAAggregationError,
        ReplicateMetricPolicy,
        aggregate_replicate_artifacts,
        aggregate_replicate_artifacts_from_disk,
    )
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
    from polyzymd.analyses.mda.comparison import (
        MDAComparisonContext,
        MDAComparisonError,
        compare_condition_artifacts,
    )
    from polyzymd.analyses.mda.frame_selection import FrameSelection
    from polyzymd.analyses.mda.job import (
        MDAAnalysisJob,
        MDAAnalysisJobError,
        MDABackendPolicy,
        MDAFunctionAdapter,
        MDAJobResult,
        MDAUniversePolicy,
    )
    from polyzymd.analyses.mda.lifecycle import MDAReplicateJobContext
    from polyzymd.analyses.mda.plugin import (
        MDAArtifactCollector,
        MDACollectorContext,
        StrictJSONMDAResultCollector,
        frame_selection_payload,
        strict_json_payload,
    )
    from polyzymd.analyses.mda.store import ArtifactStore, ArtifactStoreError
    from polyzymd.analyses.mda.universe import FileIdentity, UniverseProvenance, UniverseProvider

    assert mda.MDA_EXTENSION_API_VERSION == MDA_EXTENSION_API_VERSION == "1"
    assert mda.MDA_ARTIFACT_SCHEMA_VERSION == MDA_ARTIFACT_SCHEMA_VERSION == "1"
    assert mda.AnalysisBaseLike is AnalysisBaseLike
    assert mda.MDAnalysisExtensionError is MDAnalysisExtensionError
    assert mda.MDARunKwargs is MDARunKwargs
    assert mda.AggregatedMetric is AggregatedMetric
    assert mda.ExplicitReplicateMetricPolicy is ExplicitReplicateMetricPolicy
    assert mda.MDAAggregationContext is MDAAggregationContext
    assert mda.MDAAggregationError is MDAAggregationError
    assert mda.ReplicateMetricPolicy is ReplicateMetricPolicy
    assert mda.aggregate_replicate_artifacts is aggregate_replicate_artifacts
    assert mda.aggregate_replicate_artifacts_from_disk is aggregate_replicate_artifacts_from_disk
    assert mda.ArtifactEnvelope is ArtifactEnvelope
    assert mda.ArtifactManifest is ArtifactManifest
    assert mda.ArtifactSidecarRef is ArtifactSidecarRef
    assert mda.ComparisonArtifact is ComparisonArtifact
    assert mda.ConditionArtifact is ConditionArtifact
    assert mda.ReplicateArtifact is ReplicateArtifact
    assert mda.ArtifactStore is ArtifactStore
    assert mda.ArtifactStoreError is ArtifactStoreError
    assert mda.MDAComparisonContext is MDAComparisonContext
    assert mda.MDAComparisonError is MDAComparisonError
    assert mda.compare_condition_artifacts is compare_condition_artifacts
    assert mda.FrameSelection is FrameSelection
    assert mda.MDAAnalysisJob is MDAAnalysisJob
    assert mda.MDAAnalysisJobError is MDAAnalysisJobError
    assert mda.MDABackendPolicy is MDABackendPolicy
    assert mda.MDAFunctionAdapter is MDAFunctionAdapter
    assert mda.MDAJobResult is MDAJobResult
    assert mda.MDAUniversePolicy is MDAUniversePolicy
    assert mda.MDAReplicateJobContext is MDAReplicateJobContext
    assert mda.MDAArtifactCollector is MDAArtifactCollector
    assert mda.MDACollectorContext is MDACollectorContext
    assert mda.StrictJSONMDAResultCollector is StrictJSONMDAResultCollector
    assert mda.frame_selection_payload is frame_selection_payload
    assert mda.strict_json_payload is strict_json_payload
    assert mda.FileIdentity is FileIdentity
    assert mda.UniverseProvider is UniverseProvider
    assert mda.UniverseProvenance is UniverseProvenance
    assert set(mda.__all__) == {
        "MDA_EXTENSION_API_VERSION",
        "MDA_ARTIFACT_SCHEMA_VERSION",
        "AnalysisBaseLike",
        "MDAnalysisExtensionError",
        "MDARunKwargs",
        "AggregatedMetric",
        "ExplicitReplicateMetricPolicy",
        "MDAAggregationContext",
        "MDAAggregationError",
        "ReplicateMetricPolicy",
        "aggregate_replicate_artifacts",
        "aggregate_replicate_artifacts_from_disk",
        "ArtifactEnvelope",
        "ArtifactManifest",
        "ArtifactSidecarRef",
        "ComparisonArtifact",
        "ConditionArtifact",
        "ReplicateArtifact",
        "ArtifactStore",
        "ArtifactStoreError",
        "MDAComparisonContext",
        "MDAComparisonError",
        "compare_condition_artifacts",
        "FrameSelection",
        "MDAAnalysisJob",
        "MDAAnalysisJobError",
        "MDABackendPolicy",
        "MDAFunctionAdapter",
        "MDAJobResult",
        "MDAUniversePolicy",
        "MDAReplicateJobContext",
        "MDAArtifactCollector",
        "MDACollectorContext",
        "StrictJSONMDAResultCollector",
        "frame_selection_payload",
        "strict_json_payload",
        "FileIdentity",
        "UniverseProvider",
        "UniverseProvenance",
    }
