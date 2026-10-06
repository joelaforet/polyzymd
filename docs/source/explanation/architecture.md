# Architecture

This page explains how PolyzyMD is organized, why the major subsystems are
separated, and which boundaries matter when extending the project. It is a
conceptual map, not a step-by-step contributor tutorial.

## The high-level shape of the project

PolyzyMD follows the lifecycle of an enzyme-polymer molecular dynamics study:

1. load and validate configuration
2. build a molecular system
3. run simulation workflows locally or through SLURM
4. analyze trajectories into stored per-replicate results
5. compare conditions and create plots or reports

That lifecycle is reflected in the active package layout:

```text
src/polyzymd/
├── analyses/      # study API, analysis functions and the analyze protocol
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

Analysis and comparison behavior is concentrated in `analyses/` and the
`polyzymd analyze` command in `cli/analyze.py`.

## Why the code is split this way

The main boundary is between defining a study, running simulations, and
interpreting results. Keeping these responsibilities separate lets users and
contributors change one phase without accidentally coupling it to another.

### Configuration describes intent

`config/` holds schema and loading logic for YAML configuration, including
comparison configuration. It validates what a study should do before lower-level
builders or analyses act on it.

### Builders create simulation-ready systems

`builders/` turns input structures into simulation-ready molecular systems by
assembling enzyme, substrate, polymer, solvent, and related components. The
builder layer stays focused on construction; it does not own long-running job
or analysis policy.

### Simulation and workflow execute the study

`simulation/` runs local minimization, equilibration, checkpointing,
continuation, and production segments. `workflow/` handles orchestration around
those runs, especially SLURM job generation, resubmission, and recovery flows,
and writes the SLURM jobs of `polyzymd analyze --submit`.

`engines/` isolates engine-specific integration details such as OpenMM or
GROMACS support. This keeps high-level workflows from depending directly on one
engine's file formats or object model.

### Analyses interpret completed trajectories

`analyses/` turns completed trajectories into comparisons. A `Study` loads every
replicate of every condition from its simulation `config.yaml` as an
MDAnalysis `Universe` and removes the equilibration window. An analysis is a
plain function of atom groups or a `Universe`, which the study runs on every
production frame (`study.timeseries`) or once per replicate
(`study.per_replicate`). The study stores each replicate's result with a record
of how it was made, reduces it to one value per replicate, and computes
intervals and tests across conditions with the replicate as the sampling unit.

This design keeps trajectory processing separate from ensemble interpretation:
MDAnalysis handles per-trajectory analysis idioms, while PolyzyMD handles study
structure, result identity, aggregation, comparison, and CLI integration.

## The current `analyses/` boundary

```text
src/polyzymd/analyses/
├── study.py         # Study, Condition, Replicate: loading and the equilibration window
├── timeseries.py    # study.timeseries / per_replicate, stored results, summaries and tests
├── functions.py     # the shipped analysis functions
├── reference.py     # reference structures for RMSD, RMSF and native contacts
├── figures.py       # figures drawn from stored values
├── protocols.py     # polyzymd analyze: runs a shipped analysis and builds the report
├── universe.py      # UniverseProvider: loads each replicate, records its input files
├── identity.py      # compute_config_hash, recorded by every stored result
└── shared/          # selections, statistics, plotting and loader utilities
```

The public surface for analysis code is `polyzymd.Study` (also
`polyzymd.analyses.study.Study`), `polyzymd.analyses.functions` and
`polyzymd.analyses.analyze`. Modules and functions whose names start with `_`
are private.

## How the MDAnalysis lifecycle is divided

PolyzyMD and MDAnalysis share responsibility during trajectory analysis, but not
at the same layer.

PolyzyMD resolves topology and trajectory paths from each simulation config,
joins the production segments in order, applies the equilibration window and
stride, and loads each replicate's `Universe`. It also owns replicate
discovery, the record that decides whether a stored result is reused, storage,
aggregation, cross-condition comparison, and CLI output.

MDAnalysis owns the per-trajectory analysis idioms: selecting atoms, iterating
frames, and running `AnalysisBase`-compatible work such as
`HydrogenBondAnalysis`. A per-frame function runs through MDAnalysis
`AnalysisFromFunction`.

Conceptually, the analysis flow is:

```text
config -> builders -> simulation/workflow -> analyses -> comparison -> plots
```

Within `analyses/`, that becomes:

```text
Study.from_configs
  -> one Universe per replicate, production frames only
  -> study.timeseries or study.per_replicate
  -> polyzymd_results/<name>/<condition>/replicate_<n>/ with record.json
  -> one value per replicate
  -> summary() or compare(): intervals, tests, ProtocolReport
  -> figures drawn from the stored values
```

For code examples, see {doc}`../how_to/study_api`.

## Comparison infrastructure

- `analyses/timeseries.py` holds the per-replicate values and their
  `summary()` and `compare()`, which give Student t intervals, Welch's or
  Student's t tests, Benjamini-Hochberg correction and effect sizes.
- `analyses/shared/inferential_statistics.py` provides the statistical
  primitives such as t-tests and effect sizes.
- `analyses/protocols.py` defines the `ProtocolReport` that both
  `polyzymd analyze` and Python comparisons return.

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
  -> per-replicate results with their records
  -> condition summaries and comparisons
  -> plots and reports
```

This separation is intentional:

- users can stop after building or running
- analysis can be repeated without rebuilding simulations
- comparisons reuse stored per-replicate results whose records still match
- plotting can be rerun without remeasuring trajectories

## Design patterns you will encounter

### Lazy imports for heavy dependencies

Modules that depend on OpenMM, OpenFF, MDAnalysis, or other heavy scientific
packages often import those packages inside functions or methods instead of at
module import time. This keeps lightweight CLI and documentation operations
usable even when optional heavy dependencies are absent.

### Functions as extension points

Analysis is the primary extensibility axis. A new analysis is a function of
MDAnalysis atom groups or a `Universe` that returns a number, or one number per
label, and the study API runs it on every replicate. Nothing has to be
registered: the function is passed to `study.timeseries` or
`study.per_replicate`, and the shipped analyses in `analyses/functions.py` are
written the same way.

## Where contributors usually need to look

- **Configuration behavior:** `src/polyzymd/config/`
- **Build behavior:** `src/polyzymd/builders/`
- **Run, restart, or cluster behavior:** `src/polyzymd/simulation/` and
  `src/polyzymd/workflow/`
- **Analyses and comparisons:** `src/polyzymd/analyses/`
- **The analyze command:** `cli/analyze.py`
- **CLI commands:** `src/polyzymd/cli/`

For the chain-ID convention used by selections and interpretation, see
{doc}`residue_assignment`.

## A practical mental model

If you are new to the codebase, think in layers:

- `config` describes what should happen
- `builders` and `simulation` make it happen for one system
- `workflow` makes it practical on clusters
- `engines` isolates engine-specific details where possible
- `analyses` functions measure trajectories, and the study stores every
  per-replicate result with its record
- comparisons interpret differences across study conditions

That mental model is usually enough to find the right subsystem before diving
into module-level details or API reference pages.

## Related pages

- contributor workflows: {doc}`../contributor_guide/contributing`
- adding an analysis: {doc}`../contributor_guide/adding_an_analysis`
- chain conventions: {doc}`residue_assignment`
- SLURM usage: {doc}`../how_to/hpc_slurm`
- API reference: {doc}`../api/index`

<!-- IMAGE OPPORTUNITY: Add a left-to-right architecture diagram showing
`config -> builders -> simulation/workflow -> analyses -> comparison -> plots`,
with extension points called out at `analyses` and `workflow`. -->
