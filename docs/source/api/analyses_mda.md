# MDAnalysis Extension-Layer API

`polyzymd.analyses.mda` is the public MDAnalysis extension-layer API for
analyses that run as MDAnalysis-compatible jobs and persist PolyzyMD artifacts.
The primary contributor surface is the job, frame-selection, collector,
artifact, and artifact-store API. Use these pages for lookup. For the contributor
workflow, see {doc}`../contributor_guide/extending_analyses`.

The package is import-light: importing `polyzymd.analyses.mda` should not import
MDAnalysis or other heavy simulation dependencies. Individual jobs and helpers
perform heavy imports lazily when they actually build or run an analysis.

## API overview

| Area | Public objects | Details |
|------|----------------|---------|
| Public facade | Re-exported public MDAnalysis-layer symbols | {doc}`analyses_mda/facade` |
| Jobs | `MDAAnalysisJob`, `MDAFunctionAdapter`, `MDAJobResult`, `MDABackendPolicy`, `MDAUniversePolicy` | {doc}`analyses_mda/jobs` |
| Frame selection | `FrameSelection` | {doc}`analyses_mda/frame_selection` |
| Plugin lifecycle | `MDAReplicateJobContext`, `MDACollectorContext`, `MDAArtifactCollector`, `StrictJSONMDAResultCollector` | {doc}`analyses_mda/lifecycle` |
| Artifacts | `ReplicateArtifact`, `ConditionArtifact`, `ComparisonArtifact`, `ArtifactManifest`, `ArtifactSidecarRef` | {doc}`analyses_mda/artifacts_store` |
| Storage | `ArtifactStore`, `ArtifactStoreError` | {doc}`analyses_mda/artifacts_store` |
| Aggregation | `aggregate_replicate_artifacts`, `aggregate_replicate_artifacts_from_disk`, `ExplicitReplicateMetricPolicy` | {doc}`analyses_mda/aggregation_comparison` |
| Comparison | `compare_condition_artifacts`, `MDAComparisonContext` | {doc}`analyses_mda/aggregation_comparison` |
| Universe provenance | `UniverseProvider`, `UniverseProvenance`, `FileIdentity` | {doc}`analyses_mda/universe_primitives` |
| Shared primitives | `AnalysisBaseLike`, `MDARunKwargs`, `PairDistanceSpec`, `build_pair_distance_analysis` | {doc}`analyses_mda/universe_primitives` |

```{toctree}
:maxdepth: 1
:caption: MDAnalysis Extension Layer

analyses_mda/facade
analyses_mda/jobs
analyses_mda/frame_selection
analyses_mda/lifecycle
analyses_mda/artifacts_store
analyses_mda/aggregation_comparison
analyses_mda/universe_primitives
```
