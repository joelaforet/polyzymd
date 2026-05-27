# Aggregation and comparison

The default aggregation policy reads one finite scalar per metric from
`payload["metrics"]` or `payload["replicate_metrics"]` in each
`ReplicateArtifact`. It computes condition-level `mean`, `std`, `sem`, `n`, and
replicate `values` without loading trajectories.

The MDA comparison helper consumes `ConditionArtifact` objects, validates metric
keys, replicate identity, settings fingerprints, and aggregate-statistic
consistency, then delegates scalar statistics to the shared comparison engine.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.aggregation
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.mda.comparison
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```
