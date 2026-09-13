# Aggregation

The default aggregation policy reads one finite scalar per metric from
`payload["metrics"]` or `payload["replicate_metrics"]` in each
`ReplicateArtifact`. It computes condition-level `mean`, `std`, `sem`, `n` and
the replicate `values` without loading trajectories.

Cross-condition comparison is not here. A contract plugin compares through
`polyzymd.analyses.contract.compare_observables`, which tests one observable at
a time across conditions under a single Benjamini-Hochberg family.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.aggregation
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```
