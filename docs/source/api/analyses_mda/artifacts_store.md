# Artifact envelopes and store

The artifact models are Pydantic envelopes for stable cache files:

- `ReplicateArtifact`: output from one replicate;
- `ConditionArtifact`: aggregate over replicates for one condition;
- `ComparisonArtifact`: cross-condition statistical result;
- `ArtifactSidecarRef`: validated reference to an artifact-owned sidecar.

`ArtifactStore` reads and writes canonical `result.json`, manifest files, and
sidecars. Sidecar references are relative paths with recorded size and SHA-256
hashes.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.artifacts
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.mda.store
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```
