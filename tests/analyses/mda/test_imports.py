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


def test_import_does_not_load_heavy_simulation_modules() -> None:
    """Importing extension primitives should not import heavy optional modules."""

    heavy_modules = ("MDAnalysis", "numpy", "openmm", "openff", "parmed", "pdbfixer")
    initially_loaded = {name for name in heavy_modules if name in sys.modules}

    importlib.import_module("polyzymd.analyses.mda")
    importlib.import_module("polyzymd.analyses.mda.aggregation")
    importlib.import_module("polyzymd.analyses.mda.artifacts")
    importlib.import_module("polyzymd.analyses.mda.base")
    importlib.import_module("polyzymd.analyses.mda.frame_selection")
    importlib.import_module("polyzymd.analyses.mda.job")
    importlib.import_module("polyzymd.analyses.mda.lifecycle")
    importlib.import_module("polyzymd.analyses.mda.plugin")
    importlib.import_module("polyzymd.analyses.mda.store")
    importlib.import_module("polyzymd.analyses.mda.universe")

    for module_name in heavy_modules:
        if module_name not in initially_loaded:
            assert module_name not in sys.modules


def test_analysis_base_like_protocol_accepts_run_results_object() -> None:
    """The protocol should match objects with ``results`` and chainable ``run``."""

    from polyzymd.analyses.mda import AnalysisBaseLike, MDARunKwargs

    class FakeAnalysisBase:
        """Minimal structural stand-in for an MDAnalysis analysis object."""

        results: dict[str, Any]

        def __init__(self) -> None:
            self.results = {}
            self.run_kwargs: dict[str, Any] = {}

        def run(self, **kwargs: Any) -> "FakeAnalysisBase":
            """Store run keyword arguments and return ``self``."""

            self.run_kwargs = dict(kwargs)
            return self

    kwargs = MDARunKwargs(start=1, stop=10, step=2, verbose=False)
    analysis = FakeAnalysisBase()

    assert isinstance(analysis, AnalysisBaseLike)
    assert analysis.run(**kwargs) is analysis
    assert analysis.run_kwargs == kwargs


def test_extension_error_is_runtime_error() -> None:
    """Extension failures should use a RuntimeError subclass."""

    from polyzymd.analyses.mda import MDAnalysisExtensionError

    assert issubclass(MDAnalysisExtensionError, RuntimeError)
