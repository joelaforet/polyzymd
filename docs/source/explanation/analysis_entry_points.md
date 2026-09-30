# Which analysis entry point should I use

PolyzyMD exposes several ways into the analysis machinery, and more than one of
them can produce the same number. They differ in how much you have to write
yourself and in how much PolyzyMD guarantees about the result. This page says
which one to reach for, and why the others exist.

## The short answer

| What you want | Use |
|---|---|
| A shipped analysis, with its unit and uncertainty, from the shell | `polyzymd analyze NAME -c A/config.yaml -c B/config.yaml` |
| The same from Python | `polyzymd.analyses.analyze` |
| A measurement of your own, or a combination of shipped functions | `pz.Study.from_configs` with `study.timeseries` or `study.per_replicate` |
| A `Universe` for MDAnalysis work the study API does not cover | `replicate.universe()` of a study, or `TrajectoryLoader.load_universe` |

## Getting a number

`polyzymd analyze` and `polyzymd.analyses.analyze` are the same code path. They
build a study from the simulation config paths, run the shipped function of the
named analysis on every replicate, and return a `ProtocolReport` in which every
number states its unit, its interval and the method of the interval, its
replicate count, and a provenance block. Nothing has to be written first, and
this is the right default for the analyses PolyzyMD ships.

One report holds one result. An analysis that measures several, such as the
per-residue and whole-protein results of `sasa`, lists them all in `all_runs`
and reports the one that `--run` names, by default the first.

## Writing your own measurement

`pz.Study.from_configs` loads the same replicates in Python.
`study.timeseries(function, ...)` runs a per-frame function on every production
frame of every replicate, and `study.per_replicate(function, ...)` runs a
function once per replicate on its production frames. Both store each
replicate's result with a record of how it was made, and their results reduce
to one value per replicate whose `summary()` and `compare()` give the same
report as `polyzymd analyze`. `Timeseries.transform` combines stored series,
for example the hydrogen bonds of a catalytic triad in
{doc}`../how_to/analysis_triad_quickstart`. Use this route when no shipped
analysis measures what you need; see {doc}`analysis_api`.

## Loading trajectories yourself

`replicate.universe()` gives the `Universe` of one replicate of a study, with
its production segments in order, and `replicate.frames` the production frame
indices after the equilibration window. Use it for work that does not fit one
function per frame or per replicate, such as a principal component basis fitted
on all replicates at once.

`TrajectoryLoader.load_universe(replicate)` gives the same `Universe` from a
simulation config without a study. It chains the production segments in order
and fills in the element attribute from the atom type if possible, otherwise
from the atom name.

Loading a file with `mda.Universe()` directly leaves the element attribute empty
when the topology file has none, and does not chain the production segments.
It is not a general-purpose route.

## See also

- {doc}`../how_to/analysis_agent_protocol` for the one-command recipe.
- {doc}`../reference/analysis_protocol_report` for the report schema.
- {doc}`analysis_concepts` for the steps these entry points run.
- {doc}`analysis_statistics_best_practices` for what the reported uncertainties
  mean.
