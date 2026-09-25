# Which analysis entry point should I use

PolyzyMD exposes several ways into the analysis machinery, and more than one of
them can produce the same number. They differ in how much you have to build
yourself and in how much the framework guarantees about the result. This page
says which one to reach for, and why the others exist.

## The short answer

| What you want | Use |
|---|---|
| A number, with its unit and uncertainty, from the shell | `polyzymd analyze` |
| The same from Python | `polyzymd.analyses.protocols.analyze` |
| The full comparison workflow, including plots and cached artifacts | `polyzymd compare run` |
| To drive the pipeline from a script that already holds a config | `run_analysis` and `run_comparison` |
| A `Universe` for MDAnalysis work PolyzyMD does not cover | `TrajectoryLoader.load_universe` |
| To write a new analysis plugin | `UniverseProvider`, inside the plugin only |

## Getting a number

`polyzymd analyze` and `protocols.analyze` are the same code path. They build a
comparison in memory from simulation config paths, run the plugin pipeline and
return a `ProtocolReport` in which every number states what it is. Nothing has
to be authored first, and nothing has to be looked up afterwards to learn what
an error bar means. This is the right default, and it is the only entry point
that guarantees a unit, an interval with a named method, a replicate count and
a provenance block on every result.

The cost is that it exposes one primary metric per report. A plugin that
reports several metrics names the rest in `all_metrics`, but summarises only
the first.

## Running the full workflow

`polyzymd compare run` reads a `comparison.yaml`, runs the same pipeline, writes
the cached condition and comparison artifacts, and generates the plots. Use it
when the output is a figure or a cached artifact that later commands will read,
when the conditions are worth writing down and reusing, or when the run is
being submitted to SLURM. Its `--format agent` prints the same compact report
as `polyzymd analyze`, so moving between them costs nothing.

`compare run` is also the only route that reaches the plot settings, the
semantic colour configuration and the per-plugin settings blocks in their full
form. `polyzymd analyze --set` covers the common settings, not all of them.

## Orchestrating from Python

`run_analysis` runs one analysis for one condition and `run_comparison` runs the
whole pipeline for one analysis across conditions. Both are the layer under
`protocols.analyze`, and both expect you to build the `Condition` and settings
objects, or a `ComparisonConfig`, yourself. Use them when a script already holds
those objects, for instance a SLURM worker or a notebook that loops over many
comparison projects. They return the pipeline's own result objects, which carry
more detail than `ProtocolReport` and correspondingly more shapes to handle.

## Loading trajectories yourself

`TrajectoryLoader.load_universe(replicate)` is the canonical way to get an
MDAnalysis `Universe` from a PolyzyMD config. It chains the production segments
in order, checks their lineage, and enriches elements. Use it when you need a
measurement no plugin provides. Do not use it to reimplement one that exists:
a hand-written loop will not equilibrate uniformly, will not aggregate across
replicates the way the framework does, and will produce a number with no
provenance.

`UniverseProvider.from_config(config).load_universe(replicate)` wraps
`TrajectoryLoader` and records file hashes, engine identity and warnings into
the artifact envelope. That provenance is what makes a replicate artifact
reproducible, so the framework uses it inside the plugin lifecycle. Plugins get
it handed to them through `MDAReplicateJobContext`; analysis code outside a
plugin should not construct one, because the provenance it records has nowhere
to go.

Loading a file with `mda.Universe()` directly skips element enrichment and the
lineage check. A few plugins do it for topology-only or external-reference
loads, where no trajectory is involved. It is not a general-purpose route.

## Why there is more than one

The layers grew from the bottom. `TrajectoryLoader` came first, the plugin
framework was built on it, the CLI was built on the framework, and the protocol
was added last because the framework alone still demanded too much of a caller
who only wanted a number. Each layer is still used by the one above it, and
only the top one is meant to be the starting point. The layers inside the plugin
framework are planned to be merged into a single runner in a later release;
`TrajectoryLoader` and the protocol will stay as they are.

## See also

- {doc}`../how_to/analysis_agent_protocol` for the one-command recipe.
- {doc}`../reference/analysis_protocol_report` for the report schema.
- {doc}`analysis_concepts` for the plugin lifecycle these entry points drive.
- {doc}`analysis_statistics_best_practices` for what the reported uncertainties
  mean.
