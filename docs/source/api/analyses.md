# Study API modules

The modules below hold the study API and `polyzymd analyze`. Every analysis
that `polyzymd analyze` offers is a function in `polyzymd.analyses.functions`,
listed in {doc}`../reference/analysis_functions`, and
{doc}`../explanation/analysis_api` shows how to use the modules together.

| Module | What it holds |
|---|---|
| `polyzymd.analyses.study` | `Study`, `Condition` and `Replicate`: every replicate of every condition as an MDAnalysis `Universe` with its equilibration window removed |
| `polyzymd.analyses.timeseries` | `Study.timeseries` and `Study.per_replicate`, the stored `Timeseries` and `ReplicateValues`, their summaries and comparisons |
| `polyzymd.analyses.figures` | Figures drawn from stored values |
| `polyzymd.analyses.reference` | Reference structures (frame, centroid, average or file) for RMSD, RMSF and native contacts |
| `polyzymd.analyses.protocols` | `polyzymd analyze` in Python (`analyze`) and the `ProtocolReport` it returns |
| `polyzymd.analyses.universe` | `UniverseProvider`, which loads a replicate's `Universe` and records the path, size and modification time of each input file |
| `polyzymd.analyses.identity` | `compute_config_hash`, the hash of a simulation config that every stored result records |

## Study

```{eval-rst}
.. automodule:: polyzymd.analyses.study
   :members:
   :show-inheritance:
   :no-index:
```

## Timeseries and per-replicate values

```{eval-rst}
.. automodule:: polyzymd.analyses.timeseries
   :members:
   :show-inheritance:
   :no-index:
```

## Figures

```{eval-rst}
.. automodule:: polyzymd.analyses.figures
   :members:
   :no-index:
```

## Reference structures

```{eval-rst}
.. automodule:: polyzymd.analyses.reference
   :members:
   :no-index:
```

## polyzymd analyze in Python

```{eval-rst}
.. automodule:: polyzymd.analyses.protocols
   :members: analyze, ProtocolReport, ConditionReport, PairwiseReport, ProtocolProvenance, FUNCTION_ANALYSES
   :show-inheritance:
   :no-index:
```

## Universe loading and input file records

```{eval-rst}
.. automodule:: polyzymd.analyses.universe
   :members:
   :show-inheritance:
   :no-index:
```

## Config hash

```{eval-rst}
.. automodule:: polyzymd.analyses.identity
   :members:
   :no-index:
```

## Related API pages

- {doc}`analyses_shared`: trajectory loading, statistics, plotting and selection helpers
- {doc}`config`: `PlotSettings` and the other plot settings models
