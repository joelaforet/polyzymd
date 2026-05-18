# Analyses Plugin System API

This reference page summarizes the public `polyzymd.analyses` API surface for
plugin discovery, orchestration, statistics, the base facade, the public
MDAnalysis layer, and built-in analysis plugin packages.

PolyzyMD analyses run trajectory-native work at the replicate level and lift the
outputs into condition and comparison artifacts. The stable contributor import
surfaces are `polyzymd.analyses.base`, `polyzymd.analyses.mda`, selected
`polyzymd.analyses.shared` utilities, and the built-in plugin packages listed
below.

```{eval-rst}
.. currentmodule:: polyzymd.analyses
```

## Public package facade

The package root exposes discovery helpers and selected convenience imports. Use
the narrower modules below for detailed autodoc reference; the package facade is
summarized here to avoid duplicating class entries from the dedicated pages.

## Discovery

Discovery is package-based and uses importable modules under
`polyzymd.analyses`. A plugin can be a single module or a package with private
helper modules for plotting, result models, or MDAnalysis job construction.

| Function | Purpose |
|----------|---------|
| `list_analyses()` | Return discovered plugin classes keyed by canonical name. |
| `list_all_names()` | Return canonical names and aliases. |
| `get_analysis(name)` | Resolve a plugin class by canonical name or alias. |
| `clear_cache()` | Clear the discovery cache, primarily for tests and dynamic plugin development. |

```{eval-rst}
.. automodule:: polyzymd.analyses.discovery
   :members:
   :undoc-members:
   :show-inheritance:
```

## Orchestration

The orchestrator coordinates compute, aggregation, comparison, plotting, cache
identity, and result persistence for one analysis or a set of analyses defined
by `ComparisonConfig`.

| Function | Scope |
|----------|-------|
| `run_analysis()` | Run compute and aggregation for one analysis on one condition. |
| `run_comparison()` | Run the full lifecycle for one analysis across all configured conditions. |
| `run_all_comparisons()` | Run selected or discovered analyses from one comparison configuration. |

```{eval-rst}
.. automodule:: polyzymd.analyses.orchestrator
   :members:
   :undoc-members:
   :show-inheritance:
```

## Public base facade

`polyzymd.analyses.base` is the stable public facade for contributor imports. It
re-exports the `Analysis` base class, lifecycle context objects, scalar metric
descriptors, and comparison result models from implementation modules.

For the complete class and context reference, see {doc}`analyses_base`.

At a high level, compute-stage plugins implement `build_mda_jobs()` and
`build_mda_collector()`. Collectors produce `ReplicateArtifact` objects;
aggregation combines those into `ConditionArtifact` objects; comparison produces
`ComparisonArtifact` outputs or an active custom comparison contract for plugins
that still need specialized comparison models.

The detailed autodoc for this facade lives on {doc}`analyses_base`.

## Public MDAnalysis layer

`polyzymd.analyses.mda` is the public MDAnalysis extension layer for jobs, frame
selection, collectors, artifact envelopes, artifact storage, default aggregation,
and artifact-based comparison. The primary contributor surface is documented in
{doc}`analyses_mda`.

The detailed autodoc for this layer lives on {doc}`analyses_mda`.

## Statistics helpers

`polyzymd.analyses.stats` contains reusable scalar comparison helpers used by
the default `Analysis.compare()` path and by plugin-specific comparison code.

Key public helpers include `default_scalar_comparison()` and
`format_scalar_comparison()`.

## Shared utilities

Reusable plugin utilities live in `polyzymd.analyses.shared`. They are documented
separately on {doc}`analyses_shared`; this overview intentionally omits detailed
shared-utility autodoc blocks.

## Built-in plugin packages

| Plugin name | Public package | Primary output contract | Comparison style |
|-------------|----------------|-------------------------|------------------|
| `contacts` | `polyzymd.analyses.contacts` | Contact-event artifacts and sidecars | Custom comparison |
| `distances` | `polyzymd.analyses.distances` | Pair-distance artifacts | Custom comparison |
| `hydrogen_bonds` | `polyzymd.analyses.hydrogen_bonds` | Hydrogen-bond event artifacts | Custom comparison |
| `rmsd` | `polyzymd.analyses.rmsd` | Per-run RMSD artifacts | Custom multi-run comparison |
| `rg` | `polyzymd.analyses.rg` | Per-run radius-of-gyration artifacts | Custom multi-run comparison |
| `rmsf` | `polyzymd.analyses.rmsf` | Per-residue profile artifacts | Default scalar comparison |
| `sasa` | `polyzymd.analyses.sasa` | SASA artifacts and sidecars | Custom multi-run comparison |
| `catalytic_triad` | `polyzymd.analyses.catalytic_triad` | Pair-distance-derived artifacts | Default scalar comparison |
| `secondary_structure` | `polyzymd.analyses.secondary_structure` | Secondary-structure matrix artifacts | Default scalar comparison |

Built-in plugin packages expose their public `Analysis` subclass and supported
settings/result contracts from the package root. Helper modules with leading
underscores inside those packages are implementation details unless a page
explicitly labels them as internal developer reference.

## Private framework internals

`polyzymd.analyses._framework` is private/internal infrastructure for lifecycle
contexts, artifact I/O, comparison models, and plugin contract enforcement. It
is not a contributor import surface; contributor plugins should import public
symbols from `polyzymd.analyses.base` and `polyzymd.analyses.mda`.

## Related API pages

- {doc}`analyses_base` — base class, contexts, metrics, and comparison models
- {doc}`analyses_mda` — public MDAnalysis job/artifact layer
- {doc}`analyses_shared` — reusable shared utility modules
- {doc}`overview` — package-level API overview
