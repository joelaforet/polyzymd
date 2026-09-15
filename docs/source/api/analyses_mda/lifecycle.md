# Replicate lifecycle

`run_replicate` runs one replicate of a contract plugin. It loads the universe,
resolves the production window, builds the artifact store, calls the plugin, and
checks the artifact that comes back. `MDAReplicateJobContext` is what the plugin
receives: the loaded universe, the frame selection, the universe policy, the
artifact store, and the warnings the framework already holds.

An artifact holding a raw MDAnalysis `Results` object anywhere is rejected
rather than written, because serializing one pickles a whole universe into the
JSON.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.lifecycle
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```
