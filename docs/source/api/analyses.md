# Study API modules

The modules below hold the study API and `polyzymd analyze`. Every analysis
that `polyzymd analyze` offers is a function in `polyzymd.analyses.functions`,
listed in {doc}`../reference/analysis_functions`, and
{doc}`../how_to/study_api` shows how to use the modules together.

| Module | What it holds |
|---|---|
| `polyzymd.analyses.study` | `Study`, `Condition` and `Replicate`: every replicate of every condition as an MDAnalysis `Universe` with its equilibration window removed |
| `polyzymd.analyses.project` | `Project`: the studies of `project.yaml`, and their stored results in one table |
| `polyzymd.analyses.timeseries` | `Study.timeseries` and `Study.per_replicate`, the stored `Timeseries` and `ReplicateValues`, their summaries and comparisons |
| `polyzymd.analyses.figures` | Figures drawn from stored values |
| `polyzymd.analyses.reference` | Reference structures (frame, centroid, average or file) for RMSD, RMSF and native contacts |
| `polyzymd.analyses.protocols` | `polyzymd analyze` in Python (`analyze`) and the `ProtocolReport` it returns |
| `polyzymd.analyses.universe` | `UniverseProvider`, which loads a replicate's `Universe` |
| `polyzymd.analyses.study_file`, `results`, `user_functions`, `study_scaffold`, `study_git`, `study_metadata`, `study_freeze`, `study_upload_guide` | `study.yaml` and `data.local.yaml`, reading stored results without trajectories, running a study's own functions, creating a study folder, recording its git state, freezing it for publication, and preparing its upload |
| `polyzymd.analyses.identity` | `compute_config_hash`, the hash of a simulation config that every stored result records |

## Study

```{eval-rst}
.. automodule:: polyzymd.analyses.study
   :members:
   :show-inheritance:
   :no-index:
```

## Project

```{eval-rst}
.. automodule:: polyzymd.analyses.project
   :members:
   :no-index:
```

## Study files and stored results

```{eval-rst}
.. automodule:: polyzymd.analyses.study_file
   :members:
   :no-index:

.. automodule:: polyzymd.analyses.results
   :members:
   :no-index:

.. automodule:: polyzymd.analyses.user_functions
   :members:
   :no-index:

.. automodule:: polyzymd.analyses.study_scaffold
   :members:
   :no-index:

.. automodule:: polyzymd.analyses.study_git
   :members:
   :no-index:

.. automodule:: polyzymd.analyses.study_metadata
   :members:
   :no-index:

.. automodule:: polyzymd.analyses.study_freeze
   :members:
   :no-index:

.. automodule:: polyzymd.analyses.study_upload_guide
   :members:
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
   :members: analyze, ProtocolReport, ConditionReport, PairwiseReport, ProtocolProvenance, ANALYSES, ShippedAnalysis
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
