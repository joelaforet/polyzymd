# Universe provenance and shared MDAnalysis primitives

`UniverseProvider` wraps the existing trajectory loader and records topology and
trajectory identity for provenance.

`AnalysisBaseLike` and `MDARunKwargs` describe the lightweight protocol used by
jobs and adapters. The pair-distance helper is a public but specialized
primitive for distance-family analyses that need the same MDAnalysis-native
pair-distance behavior as built-in plugins; most contributors only need the job,
frame-selection, artifact, and store objects above.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.universe
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.mda.pair_distance
   :members: PairDistanceSpec, build_pair_distance_analysis, pair_distance_version
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.mda.base
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```
