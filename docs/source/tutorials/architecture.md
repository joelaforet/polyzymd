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
|- analysis/
|- compare/
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

### `analysis/`

Computes post-simulation metrics for individual conditions or trajectories.

### `compare/`

Aggregates results across conditions, applies statistics, and generates
comparison-ready outputs and plots.

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

### Registry-based extension points

Analysis, comparison, and plotting workflows are designed to be extensible.
Instead of modifying a central switch statement, new implementations are usually
registered through the existing registry patterns.

### Separation between per-condition and cross-condition work

The `analysis/` package answers questions about one simulation condition. The
`compare/` package answers questions across multiple conditions. Keeping those
roles separate helps maintain both code clarity and scientific interpretation.

## Where contributors usually need to look

| Goal | Start here |
|------|------------|
| add or validate config fields | `src/polyzymd/config/` |
| change build behavior | `src/polyzymd/builders/` |
| change run or restart behavior | `src/polyzymd/simulation/` and `src/polyzymd/workflow/` |
| add an analysis type | `src/polyzymd/analysis/` |
| add a comparison workflow | `src/polyzymd/compare/` |
| add or change CLI commands | `src/polyzymd/cli/` |

## A practical mental model

If you are new to the codebase, it helps to think in layers:

- `config` describes what should happen
- `builders` and `simulation` make it happen for one system
- `workflow` makes it practical on clusters
- `analysis` measures what happened
- `compare` interprets differences across studies

That mental model is usually enough to find the right subsystem before you dive
into module-level details or API reference pages.

## Related pages

- contributor workflows: {doc}`contributing`
- SLURM usage: {doc}`hpc_slurm`
- API details: {doc}`../api/overview`

<!-- IMAGE OPPORTUNITY: Add a left-to-right architecture diagram showing
`config -> builders -> simulation/workflow -> analysis -> compare -> plots`,
with extension points called out at `analysis`, `compare`, and `workflow`. -->
