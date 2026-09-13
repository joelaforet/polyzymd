# Architecture

This page explains how PolyzyMD is organized, why the major subsystems are
separated, and which boundaries matter when extending the project. It is a
conceptual map, not a step-by-step contributor tutorial.

## The high-level shape of the project

PolyzyMD follows the lifecycle of an enzyme-polymer molecular dynamics study:

1. load and validate configuration
2. build a molecular system
3. run simulation workflows locally or through SLURM
4. analyze trajectories into durable artifacts
5. compare conditions and create plots or reports

That lifecycle is reflected in the active package layout:

```text
src/polyzymd/
├── analyses/      # artifact-native analysis and comparison plugin system
├── builders/      # molecular system construction
├── cli/           # command-line entry points
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

Older `analysis/` and `compare/` directories may still exist in the source tree,
but they are not the primary architecture for new analysis or comparison work.
Current analysis and comparison behavior is concentrated in `analyses/`,
`config/comparison.py`, and `cli/compare.py`.

## Why the code is split this way

The main boundary is between defining a study, running simulations, and
interpreting results. Keeping these responsibilities separate lets users and
contributors change one phase without accidentally coupling it to another.

### Configuration describes intent

`config/` holds schema and loading logic for YAML configuration, including
comparison configuration. It validates what a study should do before lower-level
builders or analysis plugins act on it.

### Builders create simulation-ready systems

`builders/` turns input structures into simulation-ready molecular systems by
assembling enzyme, substrate, polymer, solvent, and related components. The
builder layer stays focused on construction; it does not own long-running job
or analysis policy.

### Simulation and workflow execute the study

`simulation/` runs local minimization, equilibration, checkpointing,
continuation, and production segments. `workflow/` handles orchestration around
those runs, especially SLURM job generation, resubmission, and recovery flows.

`engines/` isolates engine-specific integration details such as OpenMM or
GROMACS support. This keeps high-level workflows from depending directly on one
engine's file formats or object model.

### Analyses interpret completed trajectories

`analyses/` is the current analysis and comparison architecture. It is both a
plugin system and an artifact lifecycle. Plugins measure trajectories, produce
replicate artifacts, aggregate those artifacts per condition, compare conditions,
and optionally plot or format results.

This design keeps trajectory processing separate from ensemble interpretation:
MDAnalysis handles per-trajectory analysis idioms, while PolyzyMD handles study
structure, artifact identity, aggregation, comparison, and CLI integration.

## The current `analyses/` boundary

The analysis package is split by public surface and private implementation:

```text
src/polyzymd/analyses/
├── contract.py        # Observable, the kinds, aggregation, testing, the factory
├── base.py            # the Analysis lifecycle and the framework contexts
├── contract_plots.py  # one figure per observable kind
├── discovery.py       # plugin auto-discovery
├── orchestrator.py    # comparison workflow orchestration facade
├── protocols.py       # analyze() and ProtocolReport, the agent surface
├── stats.py           # the two direction-wording helpers the report uses
├── mda/               # public MDAnalysis extension layer
├── shared/            # reusable utilities shared by plugins
├── _framework/        # private lifecycle, I/O and validation
└── <plugin>.py        # one module per analysis
```

The important public/private boundary is intentional:

- `polyzymd.analyses.contract` is what a plugin author imports: `Observable`,
  `iter_frames` and `contract_analysis`.
- `polyzymd.analyses.base` holds `Analysis` and the four lifecycle contexts a
  plugin receives rather than builds.
- `polyzymd.analyses.mda` is the public MDAnalysis extension layer for frame
  selection, the replicate lifecycle, artifacts, artifact storage and Universe
  handling.
- `_framework/` modules are private implementation details behind the public
  facade.
- Plugin helper modules named `_*.py` are private to their plugin package unless
  that plugin explicitly documents them as public.

Contributors should not import from `_framework/` or rely on another plugin's
private helper modules. The stable import surface is deliberately narrower than
the internal package layout so PolyzyMD can evolve lifecycle internals without
breaking plugins.

## How the MDAnalysis lifecycle is divided

PolyzyMD and MDAnalysis share responsibility during trajectory analysis, but not
at the same layer.

PolyzyMD resolves topology and trajectory paths from the study context, applies
frame-selection policy, and provides or caches loaded MDAnalysis `Universe`
objects through the MDA lifecycle. It also owns replicate discovery, cache
identity, artifact storage, aggregation, cross-condition comparison, and CLI
output.

MDAnalysis owns the per-trajectory analysis idioms: selecting atoms and
iterating frames. A plugin's `compute()` uses those and returns `Observable`
objects holding raw per-frame or per-index numbers. The framework turns them
into project-level artifacts.

Conceptually, the current analysis flow is:

```text
config -> builders -> simulation/workflow -> analyses/artifacts -> comparison -> plots
```

Within `analyses/`, that becomes:

```text
plugin.compute()
  -> Observable objects, one per reported quantity
  -> ReplicateArtifact, with the per-frame series in an NPZ sidecar
  -> ConditionArtifact, one aggregate per observable across replicates
  -> ComparisonArtifact, one test per observable across conditions
  -> plots and formatted CLI output
```

There is one path. A plugin reports observables; the framework reduces each one
according to its `kind`, aggregates across replicates, tests across conditions
under one Benjamini-Hochberg family, and draws the figure that kind calls for.
Plots read cached artifacts and sidecars rather than reloading trajectories or
rerunning the compute stage.

For concrete commands and code examples, see
{doc}`../contributor_guide/analysis_plugins/index`.

## Comparison infrastructure is distributed

There is no separate active comparison stack that new plugins should target.
Comparison behavior is distributed across focused modules:

- `config/comparison.py` describes comparison and plotting settings.
- `cli/compare.py` exposes the `polyzymd compare` command group.
- `analyses/contract.py` implements aggregation across replicates and the
  cross-condition tests, with one correction family per run.
- `analyses/contract_plots.py` draws the figure each observable kind calls for.
- `analyses/shared/inferential_statistics.py` provides lower-level statistical
  primitives such as t-tests, ANOVA, and effect sizes.
- `analyses/mda/` provides the artifact layer that carries replicate,
  condition, and comparison data through the lifecycle.

This distribution keeps comparison close to the artifacts it consumes while
still allowing the CLI and configuration layers to remain stable entry points.

## Supporting packages

### `core/` and `utils/`

`core/` and `utils/` provide shared infrastructure such as common types,
experimental workflow labeling, and helper functionality that should not be
duplicated across the package.

### `data/` and `templates/`

`data/` stores bundled package resources such as force-field or template data
that need to ship with PolyzyMD. It is not a user results directory.

`templates/` contains packaged templates and examples used by scaffolding and
setup flows. These are starting points for generated files, not the
authoritative runtime schema.

### `exporters/`

`exporters/` contains format-export support for moving PolyzyMD outputs into
other molecular simulation ecosystems. Exporters sit at the edge of the package
rather than in the core build or simulation lifecycle.

## How data moves through the system

At a conceptual level, data moves from declared intent to generated evidence:

```text
config.yaml
  -> validated config objects
  -> system builders
  -> simulation objects and run directories
  -> local or SLURM execution
  -> ReplicateArtifact files and sidecars
  -> condition artifacts
  -> comparison artifacts or documented custom outputs
  -> plots and reports
```

This separation is intentional:

- users can stop after building or running
- analysis can be repeated without rebuilding simulations
- comparison workflows can reuse cached analysis outputs
- plotting can be rerun without recomputing statistics
- plugin internals can change while public artifact and facade contracts remain
  stable

## Design patterns you will encounter

### Lazy imports for heavy dependencies

Modules that depend on OpenMM, OpenFF, MDAnalysis, or other heavy scientific
packages often import those packages inside functions or methods instead of at
module import time. This keeps lightweight CLI and documentation operations
usable even when optional heavy dependencies are absent.

### Plugin-based extension points

Analysis is the primary extensibility axis. A plugin is one module under
`analyses/` holding a settings model, a `compute()` and one call to
`contract_analysis()`. The framework discovers it through `pkgutil`, so
contributors do not need registries, decorators or core imports to make a plugin
available.

The reason for this design is the open-closed principle: new analyses should be
added by extension, not by modifying the orchestrator or CLI every time a metric
is introduced.

### Public facade over private lifecycle internals

The analysis framework uses private modules internally because lifecycle code,
artifact I/O, and validation contracts are implementation details. The public
facade keeps contributor imports stable while giving maintainers room to improve
internals.

This is why documentation points contributors to `polyzymd.analyses.contract`,
`polyzymd.analyses.base` and `polyzymd.analyses.mda`, not to `_framework/`.

## Where contributors usually need to look

- **Configuration behavior:** `src/polyzymd/config/`
- **Build behavior:** `src/polyzymd/builders/`
- **Run, restart, or cluster behavior:** `src/polyzymd/simulation/` and
  `src/polyzymd/workflow/`
- **Analysis and comparison plugins:** `src/polyzymd/analyses/`
- **Comparison configuration and CLI entry points:** `config/comparison.py` and
  `cli/compare.py`
- **CLI commands:** `src/polyzymd/cli/`

For the chain-ID convention used by selections and interpretation, see
{doc}`residue_assignment`.

## A practical mental model

If you are new to the codebase, think in layers:

- `config` describes what should happen
- `builders` and `simulation` make it happen for one system
- `workflow` makes it practical on clusters
- `engines` isolates engine-specific details where possible
- `analyses` plugins measure trajectories and preserve evidence as artifacts
- comparison workflows interpret differences across study conditions

That mental model is usually enough to find the right subsystem before diving
into module-level details or API reference pages.

## Related pages

- contributor workflows: {doc}`../contributor_guide/contributing`
- writing an analysis: {doc}`../contributor_guide/analysis_plugins/index`
- chain conventions: {doc}`residue_assignment`
- SLURM usage: {doc}`../how_to/hpc_slurm`
- API reference: {doc}`../api/index`

<!-- IMAGE OPPORTUNITY: Add a left-to-right architecture diagram showing
`config -> builders -> simulation/workflow -> analyses/artifacts -> comparison -> plots`,
with extension points called out at `analyses` and `workflow`. -->
