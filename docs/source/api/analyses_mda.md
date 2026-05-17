# MDAnalysis Extension-Layer API

`polyzymd.analyses.mda` is the public extension-layer API for analyses that run
as MDAnalysis-compatible jobs and persist PolyzyMD artifacts. Use this page for
lookup. For the contributor workflow, see
{doc}`../contributor_guide/extending_analyses`.

## API overview

| Area | Public objects |
|------|----------------|
| Jobs | `MDAAnalysisJob`, `MDAFunctionAdapter`, `MDAJobResult`, `MDABackendPolicy`, `MDAUniversePolicy` |
| Frame selection | `FrameSelection` |
| Plugin lifecycle | `MDAReplicateJobContext`, `MDACollectorContext`, `MDAArtifactCollector`, `StrictJSONMDAResultCollector` |
| Artifacts | `ReplicateArtifact`, `ConditionArtifact`, `ComparisonArtifact`, `ArtifactManifest`, `ArtifactSidecarRef` |
| Storage | `ArtifactStore`, `ArtifactStoreError` |
| Aggregation | `aggregate_replicate_artifacts`, `aggregate_replicate_artifacts_from_disk`, `ExplicitReplicateMetricPolicy` |
| Comparison | `compare_condition_artifacts`, `MDAComparisonContext` |
| Universe provenance | `UniverseProvider`, `UniverseProvenance`, `FileIdentity` |
| Shared primitives | `PairDistanceSpec`, `build_pair_distance_analysis`, `AnalysisBaseLike`, `MDARunKwargs` |

The package is import-light: importing `polyzymd.analyses.mda` should not import
MDAnalysis or other heavy simulation dependencies. Individual jobs and helpers
perform heavy imports lazily when they actually build or run an analysis.

## Public facade

```{eval-rst}
.. automodule:: polyzymd.analyses.mda
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Job execution

`MDAAnalysisJob` represents one named MDAnalysis-compatible analysis run on one
replicate. It forwards `FrameSelection` kwargs to the wrapped object's `run()`
method and returns an `MDAJobResult` with the completed analysis and its results.

`MDAFunctionAdapter` is the simple-function path used by the default scaffold.
It calls a function once with a loaded universe and frame-selection kwargs; it
does not implement a custom frame loop. Functions that need per-frame iteration
should use MDAnalysis primitives internally or be replaced by an `AnalysisBase`
subclass.

`MDABackendPolicy` controls optional MDAnalysis internal backends. The default
policy forwards no backend kwargs so PolyzyMD-level parallelism remains the
default. In `comparison.yaml`, the top-level `mda_backend_policy` section maps
to this object and is intentionally opt-in:

```yaml
mda_backend_policy:
  backend: "multiprocessing"
  n_workers: 2
  n_parts: 2
```

Leave the section empty or omit it to forward no backend kwargs. Function-adapter
jobs reject non-default backend policies; use an `AnalysisBase`-compatible job
when opting into MDAnalysis internal parallelism.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.job
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Frame selection

`FrameSelection` maps PolyzyMD equilibration/window semantics to MDAnalysis
`AnalysisBase.run()` kwargs. It supports either `start`/`stop`/`step` or an
explicit `frames` selector, but not both.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.frame_selection
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Plugin lifecycle and collectors

`MDAReplicateJobContext` is passed to `Analysis.build_mda_jobs()`. It carries the
loaded universe, frame selection, universe policy, settings, replicate identity,
and artifact output location.

`MDACollectorContext` is passed to `Analysis.build_mda_collector()`. A collector
maps completed `MDAJobResult` objects into one `ReplicateArtifact`. Collectors
must convert raw MDAnalysis `Results` objects into JSON-compatible payloads or
sidecars before returning the artifact.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.lifecycle
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.mda.plugin
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Artifact envelopes and store

The artifact models are Pydantic envelopes for stable cache files:

- `ReplicateArtifact`: output from one replicate;
- `ConditionArtifact`: aggregate over replicates for one condition;
- `ComparisonArtifact`: cross-condition statistical result;
- `ArtifactSidecarRef`: validated reference to an artifact-owned sidecar.

`ArtifactStore` reads and writes canonical `result.json`, manifest files, and
sidecars. Sidecar references are relative paths with recorded size and SHA-256
hashes.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.artifacts
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.mda.store
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Aggregation and comparison

The default aggregation policy reads one finite scalar per metric from
`payload["metrics"]` or `payload["replicate_metrics"]` in each
`ReplicateArtifact`. It computes condition-level `mean`, `std`, `sem`, `n`, and
replicate `values` without loading trajectories.

The MDA comparison helper consumes `ConditionArtifact` objects, validates metric
keys, replicate identity, settings fingerprints, and aggregate-statistic
consistency, then delegates scalar statistics to the shared comparison engine.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.aggregation
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.mda.comparison
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Universe provenance and shared MDAnalysis primitives

`UniverseProvider` wraps the existing trajectory loader and records topology and
trajectory identity for provenance. Shared primitives such as the pair-distance
builder are internal building blocks used by built-in analyses and available for
contributors who need the same MDAnalysis-native behavior.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.universe
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.mda.pair_distance
   :members: PairDistanceSpec, build_pair_distance_analysis, pair_distance_version
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.mda.base
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```
