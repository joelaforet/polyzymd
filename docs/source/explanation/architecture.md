# Architecture

This page explains how PolyzyMD is organized, why the major subsystems are
separated, and where contributors should look when they need to extend a
workflow.

## The high-level shape of the project

PolyzyMD is organized around a simulation lifecycle:

1. load and validate configuration
2. build a molecular system
3. run simulation workflows locally or through SLURM
4. analyze trajectories
5. compare conditions and plot results

That lifecycle is reflected in the package layout:

```text
src/polyzymd/
├── analysis/      # legacy/specialized analysis support
├── analyses/      # artifact-native analysis plugin system
├── builders/      # molecular system construction
├── cli/           # command-line entry points
├── compare/       # legacy comparison components
├── config/        # simulation and comparison configuration
├── core/          # shared domain types
├── data/          # bundled package data
├── engines/       # engine-specific integration layer
├── exporters/     # output format exporters
├── simulation/    # local simulation execution
├── templates/     # packaged example/scaffold templates
├── utils/         # general utilities
└── workflow/      # orchestration and SLURM support
```

## What each area is responsible for

### `cli/`

Defines the command-line interface and maps user commands onto the lower-level
workflow code.

### `analysis/`

Contains older or specialized analysis support that has not become the primary
contributor extension surface. New user-facing analysis plugins should generally
target `analyses/` unless maintainers point to a compatibility path here.

### `compare/`

Holds legacy comparison, plotting, and result components used by parts of the
older comparison stack. The current plugin lifecycle concentrates new analysis
comparison behavior under `analyses/`.

### `config/`

Holds the schema and loading logic for YAML configuration. If a user-facing
setting needs validation, this is usually the first place to inspect.

### `builders/`

Turns input structures into a simulation-ready system by assembling enzyme,
substrate, polymer, and solvent components.

### `simulation/`

Runs minimization, equilibration, continuation, checkpoints, and production
segments.

### `workflow/`

Handles orchestration around the simulation engine, especially SLURM job
generation, resubmission, and recovery flows.

### `engines/`

Provides engine-specific integration points, including OpenMM and GROMACS
support modules. It is a boundary for engine details rather than the place where
high-level workflows are defined.

### `data/`

Stores bundled package data such as force-field or template resources that need
to ship with PolyzyMD. It is supporting material for code paths elsewhere, not a
user results directory.

### `templates/`

Contains packaged templates and examples used by scaffolding and setup flows.
These are starting points for generated files, not the authoritative runtime
configuration schema.

### `analyses/`

The **plugin system** — the primary extension point for contributors. The
default contributor path is artifact-native and MDAnalysis-backed: a plugin
subclasses `Analysis`, builds `MDAAnalysisJob` objects, and provides collectors
that convert completed jobs into `ReplicateArtifact` objects. PolyzyMD then owns
Universe loading, cache identity, default artifact aggregation, cross-condition
comparison, statistics, and CLI output. Each analysis plugin still participates
in a unified lifecycle, but not every plugin uses every stage:
compute → aggregate → compare → plot → format.

The analysis package itself is split by responsibility:

```text
src/polyzymd/analyses/
├── base.py          # stable contributor facade: Analysis, contexts, metrics
├── discovery.py     # plugin auto-discovery
├── orchestrator.py  # comparison workflow orchestration facade
├── stats.py         # default scalar comparison pipeline
├── _framework/      # private/internal lifecycle, I/O, and contracts
├── mda/             # public MDAnalysis job, frame-selection, artifact layer
├── shared/          # contributor-shared reusable utilities
└── <plugin>/        # built-in and contributed analysis plugins
```

The private `_framework/` modules support the public facade; they are not the
contributor import surface. Contributor-facing docs and plugins should use
`polyzymd.analyses.base`, `polyzymd.analyses.mda`, and documented shared
utilities instead.

Within that lifecycle, PolyzyMD now draws a sharper boundary between
trajectory-level work and ensemble-level work:

- **MDAnalysis owns per-trajectory analysis** through function-adapter jobs or
  AnalysisBase-compatible ``run(...)`` objects
- **PolyzyMD owns ensemble/comparison workflow** including replicate discovery,
  caching, default artifact aggregation, cross-condition statistics, and CLI
  output
- **Collectors turn completed MDAnalysis jobs into `ReplicateArtifact` objects**
  with payload, metadata, provenance, warnings, and sidecar references
- **Scalar metrics for default comparison live in `payload["metrics"]`** so the
  artifact aggregator can preserve replicate-level values without reducing raw
  arrays or event tables implicitly
- **Default aggregation combines replicate artifacts into a `ConditionArtifact`**
  for each condition before cross-condition comparison
- **Comparison produces a `ComparisonArtifact`** on the canonical path, or an
  explicitly active custom output contract for plugins that still expose one
- **Composition is preferred over mixins or deep inheritance** so trajectory-
  native plugins can provide small hooks instead of re-implementing the full
  framework lifecycle
- **Advanced trajectory-native plugins keep per-trajectory logic in
  MDAnalysis-compatible job modules**
- **Derived analyses stay outside `Analysis` unless they truly process
  trajectories**; post-processing of already-aggregated outputs should remain a
  higher-level PolyzyMD concern

To add a new analysis, run `polyzymd new-analysis <name>` to scaffold the
single-file artifact-native pattern. The public default path:

- subclasses `Analysis`
- builds `MDAAnalysisJob` objects, often with `MDAAnalysisJob.from_function(...)`
- implements `build_mda_collector(...)` so completed jobs become
  `ReplicateArtifact` objects
- places scalar metrics for the default comparison pipeline in
  `payload["metrics"]`
- stores large arrays, event tables, or other rich outputs as artifact payloads
  or validated sidecars that plots can read later

Advanced packages use the explicit MDAnalysis job lifecycle:

- implement `build_mda_jobs(...)` and, when needed, `build_mda_collector(...)`
  for trajectory-native replicate work
- set `has_compute_stage = False` for compare-only plugins
- implement `aggregate()` only when `has_aggregate_stage = True`

`extract_metrics()` remains useful for compatibility with non-artifact aggregate
results and older default-comparison plugins, but it is not the preferred place
to put scalar metrics for new artifact-native analyses. New plugins should make
their default-comparison metrics explicit in artifact payloads and let PolyzyMD
perform the standard aggregation, statistics, ranking, formatting, and CLI
presentation.

Plots have the same boundary: they should read cached artifacts and sidecars,
not reload trajectories or rerun compute-stage analyses.

See {doc}`../contributor_guide/extending_analyses` for the full guide.

### Comparison infrastructure (distributed)

Comparison functionality is split across focused modules:

- `config/comparison.py` for comparison config and plotting settings
- `cli/compare.py` for `polyzymd compare` subcommands
- `analyses/stats.py` for the default scalar comparison pipeline, condition
  summaries, rankings, pairwise comparisons, and human-readable formatting
- `analyses/shared/inferential_statistics.py` for lower-level statistical
  primitives such as t-tests, ANOVA, and effect sizes
- `analyses/mda/` for MDAnalysis jobs, frame selection, artifacts, artifact
  storage, aggregation, and comparison inputs
- `analyses/shared/paths.py` for label/path helpers such as `sanitize_label()`

Established analysis package plugins often delegate plotting to `_plotters.py`
modules, but that is optional organization for larger plugins. Simple
single-file plugins can keep plotting in `plot()` when they need figures at all.

### `core/` and `utils/`

Provide shared infrastructure such as common types, experimental workflow
labeling, and helper functionality that should not be duplicated across the
package.

### `exporters/`

Contains format-export support for moving PolyzyMD outputs into other molecular
simulation ecosystems. Exporters sit at the edge of the package rather than in
the core build or simulation lifecycle.

## How data moves through the system

At a conceptual level, the flow looks like this:

```text
config.yaml
  -> config schema
  -> system builders
  -> OpenMM-ready simulation objects
  -> local or SLURM execution
  -> ReplicateArtifact files and sidecars
  -> ConditionArtifact files
  -> ComparisonArtifact or active custom comparison output
  -> plots and reports
```

This separation is intentional:

- users can stop after building or running
- analysis can be repeated without rebuilding simulations
- comparison workflows can reuse cached analysis outputs
- plotting can be rerun without recomputing the underlying statistics

## Design patterns you will encounter

### Lazy imports for heavy dependencies

Modules that depend on OpenMM or MDAnalysis often import those packages inside
functions or methods instead of at module import time. This keeps lightweight
CLI operations usable even when optional heavy dependencies are absent.

### Plugin-based extension points

Analysis is the primary extensibility axis. New analyses are usually single
files or packages that subclass `Analysis`. The framework discovers both shapes
automatically via `pkgutil` — no registries, no decorators, no imports needed.
Use `polyzymd new-analysis <name>` for the default single-file MDAnalysis-native
scaffold, or `--advanced` for a package scaffold with lazy `_mda.py` helpers.

### Separation between per-condition and cross-condition work

The unified `analyses/` lifecycle still handles both scopes in one plugin
contract. Public plugins usually:

- use `MDAAnalysisJob.from_function(...)` when one function-adapter job is
  sufficient
- collect completed jobs into `ReplicateArtifact` objects with scalar default
  comparison metrics stored under `payload["metrics"]`

Advanced plugins can also:

- use the MDAnalysis job path with `build_mda_jobs(...)` and
  `build_mda_collector(...)`
- skip compute entirely with `has_compute_stage = False`

Per-condition aggregation is likewise optional and is required only when
`has_aggregate_stage = True`. Cross-condition work still happens through
`compare()`, `plot()`, and `format()`.

The important migration rule is that simple contributors should not reimplement
the lifecycle. The default scaffold creates a single-file `Analysis` plugin
that uses `MDAAnalysisJob.from_function(...)`, returns a `ReplicateArtifact`,
stores default scalar metrics in the artifact payload, and lets PolyzyMD own
replicate discovery, caching, artifact aggregation, cross-condition statistics,
formatting, and CLI output. Advanced package plugins use
`AnalysisBase`-compatible jobs when MDAnalysis should own per-frame trajectory
iteration directly.

## Where contributors usually need to look

- **Add or validate config fields:** start in `src/polyzymd/config/`.
- **Change build behavior:** start in `src/polyzymd/builders/`.
- **Change run or restart behavior:** start in `src/polyzymd/simulation/` and
  `src/polyzymd/workflow/`.
- **Add an analysis type:** start in `src/polyzymd/analyses/`. New plugins
  usually subclass `Analysis`, build `MDAAnalysisJob` objects, return
  `ReplicateArtifact` objects through collectors, and let the artifact layer
  handle default aggregation and comparison.
- **Add comparison statistics:** use `src/polyzymd/analyses/stats.py` for
  default scalar comparison behavior and
  `src/polyzymd/analyses/shared/inferential_statistics.py` for lower-level
  statistical tests and effect-size primitives.
- **Add or change CLI commands:** start in `src/polyzymd/cli/`.

## A practical mental model

If you are new to the codebase, it helps to think in layers:

- `config` describes what should happen
- `builders` and `simulation` make it happen for one system
- `workflow` makes it practical on clusters
- `engines` isolates engine-specific details where possible
- `analyses` plugins measure and compare what happened through artifacts
- comparison workflows interpret differences across studies

That mental model is usually enough to find the right subsystem before you dive
into module-level details or API reference pages.

## Related pages

- contributor workflows: {doc}`../contributor_guide/contributing`
- extending analyses: {doc}`../contributor_guide/extending_analyses`
- SLURM usage: {doc}`../how_to/hpc_slurm`
- API details: {doc}`../api/overview`

<!-- IMAGE OPPORTUNITY: Add a left-to-right architecture diagram showing
`config -> builders -> simulation/workflow -> analysis -> analyses -> comparison workflows -> plots`,
with extension points called out at `analyses` and `workflow`. -->
