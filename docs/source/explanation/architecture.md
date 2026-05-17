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
|- cli/
|- config/
|- builders/
|- simulation/
|- workflow/
|- analyses/     # ★ plugin system — unified analysis lifecycle
|- config/comparison.py  # comparison config and plot settings
|- exporters/
|- core/
`- utils/
```

## What each area is responsible for

### `cli/`

Defines the command-line interface and maps user commands onto the lower-level
workflow code.

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

### `analyses/`

The **plugin system** — the primary extension point for contributors. The
default contributor path is MDAnalysis-native: a generated plugin builds an
`MDAAnalysisJob`, returns a canonical artifact with explicit metrics, and lets
PolyzyMD handle Universe loading, replicate cache identity, aggregation,
comparison, and CLI integration. Each analysis plugin still participates in a
unified lifecycle, but not every plugin uses every stage:
compute → aggregate → compare → plot → format.

Within that lifecycle, PolyzyMD now draws a sharper boundary between
trajectory-level work and ensemble-level work:

- **MDAnalysis owns per-trajectory analysis** through function-adapter jobs or
  AnalysisBase-compatible ``run(...)`` objects
- **PolyzyMD owns ensemble/comparison workflow** including replicate discovery,
  caching, aggregation, cross-condition statistics, plotting, and CLI output
- **Composition is preferred over mixins or deep inheritance** so trajectory-
  native plugins can provide small hooks instead of re-implementing the full
  framework lifecycle
- **Advanced trajectory-native plugins keep per-trajectory logic in
  MDAnalysis-compatible job modules**
- **Derived analyses stay outside `Analysis` unless they truly process
  trajectories**; post-processing of already-aggregated outputs should remain a
  higher-level PolyzyMD concern

To add a new analysis, run `polyzymd new-analysis <name>` to scaffold the
single-file MDAnalysis-native pattern. The public default path:

- subclasses `Analysis`
- implements `build_mda_jobs(...)` with `MDAAnalysisJob.from_function(...)`
- implements `build_mda_collector(...)` to return a `ReplicateArtifact`
- implements `extract_metrics(...)` for default scalar comparison

Advanced packages use the explicit MDAnalysis job lifecycle:

- implement `build_mda_jobs(...)` and, when needed, `build_mda_collector(...)`
  for trajectory-native replicate work
- set `has_compute_stage = False` for compare-only plugins
- implement `aggregate()` only when `has_aggregate_stage = True`

See {doc}`../contributor_guide/extending_analyses` for the full guide.

### Comparison infrastructure (distributed)

Comparison functionality is split across focused modules:

- `config/comparison.py` for comparison config and plotting settings
- `cli/compare.py` for `polyzymd compare` subcommands
- `analyses/shared/inferential_statistics.py` for t-tests, ANOVA, and effect sizes
- `analyses/shared/result_io.py` for result discovery and loading
- `analyses/shared/paths.py` for label/path helpers such as `sanitize_label()`

Established analysis package plugins often delegate plotting to `_plotters.py`
modules, but that is optional organization for larger plugins. Simple
single-file plugins can keep plotting in `plot()` when they need figures at all.

### `core/` and `utils/`

Provide shared infrastructure such as common types, experimental workflow
labeling, and helper functionality that should not be duplicated across the
package.

## How data moves through the system

At a conceptual level, the flow looks like this:

```text
config.yaml
  -> config schema
  -> system builders
  -> OpenMM-ready simulation objects
  -> local or SLURM execution
  -> analysis results on disk
  -> cross-condition comparisons
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

- use `MDAAnalysisJob.from_function(...)` when one function-adapter job is sufficient

Advanced plugins can also:

- use the MDAnalysis job path with `build_mda_jobs(...)` and
  `build_mda_collector(...)`
- skip compute entirely with `has_compute_stage = False`

Per-condition aggregation is likewise optional and is required only when
`has_aggregate_stage = True`. Cross-condition work still happens through
`compare()`, `plot()`, and `format()`.

The important migration rule is that simple contributors should not reimplement
the lifecycle. Measurement plugins provide `measure()` only; advanced
MDAnalysis-native plugins let MDAnalysis own per-trajectory iteration while
PolyzyMD owns replicate discovery, caching, ensemble aggregation,
cross-condition statistics, plotting, and CLI output.

## Where contributors usually need to look

| Goal | Start here |
|------|------------|
| add or validate config fields | `src/polyzymd/config/` |
| change build behavior | `src/polyzymd/builders/` |
| change run or restart behavior | `src/polyzymd/simulation/` and `src/polyzymd/workflow/` |
| add an analysis type | `src/polyzymd/analyses/` (default: single-file `ScalarMeasurementAnalysis`; advanced: package `Analysis` plugin) |
| add comparison statistics | `src/polyzymd/analyses/shared/inferential_statistics.py` |
| add or change CLI commands | `src/polyzymd/cli/` |

## A practical mental model

If you are new to the codebase, it helps to think in layers:

- `config` describes what should happen
- `builders` and `simulation` make it happen for one system
- `workflow` makes it practical on clusters
- `analyses` plugins measure and compare what happened
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
