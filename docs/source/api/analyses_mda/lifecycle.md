# Plugin lifecycle and collectors

`MDAReplicateJobContext` is passed to `Analysis.build_mda_jobs()`. It carries the
loaded universe, frame selection, universe policy, settings, replicate identity,
and artifact output location.

`MDACollectorContext` is passed to `Analysis.build_mda_collector()`. A collector
maps completed `MDAJobResult` objects into one `ReplicateArtifact`. Collectors
must convert raw MDAnalysis `Results` objects into JSON-compatible payloads or
sidecars before returning the artifact.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.lifecycle
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.mda.plugin
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```
