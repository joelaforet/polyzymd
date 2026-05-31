# How PolyzyMD analysis plugins work

PolyzyMD analysis plugins turn trajectory calculations into reusable ensemble
evidence. A plugin is not only a loop over frames; it is a small participant in
a larger workflow that must compare conditions, preserve provenance, reuse
cached work, render figures, and report results consistently across replicate
simulations.

That is why the current analysis architecture separates **MDAnalysis-native
trajectory work** from **PolyzyMD ensemble workflow ownership**. Contributors
write an `Analysis` subclass and describe the per-replicate trajectory jobs it
needs. PolyzyMD then handles how those per-replicate outputs become condition,
comparison, plotting, and formatting artifacts.

## The central idea

The conceptual lifecycle is:

```text
Analysis subclass
  -> build_mda_jobs()
  -> MDAAnalysisJob
  -> collector
  -> ReplicateArtifact
  -> ConditionArtifact
  -> ComparisonArtifact or custom comparison output
  -> plot() / format()
```

Each stage has a narrow responsibility. The plugin defines scientific work and
how to interpret its results. The framework supplies consistent orchestration,
storage, and cross-condition behavior.

## Why plugins build MDAAnalysisJob objects

Direct loops over frames inside a plugin-specific replicate hook are easy to
imagine, but they make every plugin responsible for concerns that are not
specific to its science: frame selection, job identity, cache behavior, artifact
shape, provenance, aggregation inputs, parallel execution, and safe plotting
inputs.

PolyzyMD instead asks compute-stage plugins to build `MDAAnalysisJob` objects
and collectors. This keeps the frame-level calculation close to MDAnalysis while
letting PolyzyMD keep ownership of the ensemble workflow.

| Owner | Responsibility |
| --- | --- |
| **MDAnalysis** | Loading trajectories, selecting atoms, iterating selected frames, and running trajectory-native calculations. |
| **Plugin code** | Choosing selections, defining settings, building `MDAAnalysisJob` objects, and collecting completed job output into meaningful plugin artifacts. |
| **PolyzyMD** | Cache identity, artifact storage, replicate-to-condition aggregation, cross-condition comparison, provenance, plotting lifecycle, and CLI formatting lifecycle. |

This split makes plugins easier to review. A maintainer can ask whether the
trajectory job computes the right quantity, whether the collector returns the
right artifact, and whether later stages read artifacts instead of rerunning
trajectory work.

## Public API boundaries

Contributor-facing plugins should use stable public facades:

- `polyzymd.analyses.base` for `Analysis`, lifecycle contexts, metrics, and
  comparison result models.
- `polyzymd.analyses.mda` for `MDAAnalysisJob`, frame selection, artifacts,
  artifact stores, and MDAnalysis integration concepts.
- Documented utilities from `polyzymd.analyses.shared` when a shared helper
  already fits the problem.

`polyzymd.analyses._framework` is private/internal. It exists behind the public
facades and should not appear in contributor imports or examples except when a
document explicitly labels it as internal implementation detail.

## What artifacts mean

Artifacts are the durable boundary between lifecycle stages. They let PolyzyMD
rerun, aggregate, compare, and plot without depending on transient Python
objects created during a trajectory pass.

| Artifact concept | Meaning |
| --- | --- |
| `ReplicateArtifact` | The validated result of one plugin for one replicate. It is produced by a collector after the `MDAAnalysisJob` finishes. |
| `ConditionArtifact` | The condition-level summary assembled from replicate artifacts and any referenced sidecars. |
| `ComparisonArtifact` | The cross-condition comparison output used by default comparison workflows, or the canonical artifact form of a plugin comparison result. |
| Custom comparison output | Some established plugins intentionally return custom comparison models. They should still follow the active output contract and remain saveable/loadable. |

Collectors should not serialize raw MDAnalysis Results objects. Those objects are
useful during compute, but they are not the storage contract for PolyzyMD
analysis results. The collector translates compute output into a
`ReplicateArtifact` with a stable payload and provenance.

## Sidecars keep artifacts practical

Artifacts should contain compact, validated metadata and summary values. Large
arrays, per-frame tables, event streams, or dense per-residue outputs belong in
sidecar files referenced by the artifact.

This design keeps artifact JSON readable and versionable while preserving access
to richer data for aggregation, comparison, and plots. A sidecar is not an
arbitrary cache file that a later stage discovers by filename. It is part of the
artifact contract: the artifact records what the sidecar is, where it is, and
why it belongs to that result.

## Plotting is artifact-only

`plot()` is a rendering stage, not a hidden compute stage. It should load cached
artifacts and referenced sidecars, then render figures from those durable
outputs. This artifact-only plotting rule prevents plots from silently changing
scientific results, rerunning expensive trajectory work, or depending on data
that was not part of the saved analysis result.

If a plot needs data that is too large for the artifact payload, the compute or
aggregation stage should write a validated sidecar and register it on the
artifact. The plot can then load that sidecar through the artifact record.

## Simple and advanced plugin shapes

PolyzyMD supports small plugins and larger package-style plugins. The choice is
about maintainability, not status.

| Choose this shape | When it fits | Typical structure |
| --- | --- | --- |
| **Simple plugin module** | The analysis has one clear metric, limited plotting, and straightforward settings. | One module with the `Analysis` subclass, settings, job construction, collection, metrics, and perhaps a small plot method. |
| **Advanced plugin package** | The analysis has multiple jobs, richer results, several plot types, substantial sidecar handling, or custom comparison output. | A package with a public `__init__.py` for the `Analysis` subclass and private helper modules such as `_mda.py`, `_plotters.py`, `_results.py`, or `_formatters.py`. |

The same lifecycle applies in both cases. A package layout only moves complexity
into smaller files; it does not change the public import boundary for
contributors or the artifact contract with PolyzyMD.

## How to think about plugin responsibility

A healthy analysis plugin answers three questions clearly:

1. **What trajectory-native work is needed?** That belongs in MDAnalysis jobs
   built by `build_mda_jobs()`.
2. **What durable result should represent one replicate?** That belongs in the
   collector that returns a `ReplicateArtifact` and registers sidecars when
   needed.
3. **What comparison or presentation should users see?** That belongs in
   aggregation, comparison, `plot()`, and `format()` stages that consume saved
   artifacts rather than raw trajectory state.

Keeping these questions separate is what lets PolyzyMD add new analyses without
changing core orchestration code.

## What to read next

- {doc}`../extending_analyses` for the current full implementation guide.
- {doc}`../../api/analyses` for the analysis plugin API overview.
- {doc}`../../api/analyses_base` for the public `polyzymd.analyses.base`
  facade.
- {doc}`../../api/analyses_mda` for the public MDAnalysis integration API.
- {doc}`../../api/analyses_shared` for documented shared analysis utilities.
