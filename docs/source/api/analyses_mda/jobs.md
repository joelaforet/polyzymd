# Job execution

`MDAAnalysisJob` represents one named MDAnalysis-compatible analysis run on one
replicate. It forwards `FrameSelection` kwargs to the wrapped object's `run()`
method and returns an `MDAJobResult` with the completed analysis and its results.

`MDAFunctionAdapter` is the simple-function path used by the default scaffold.
It calls a function once with a loaded universe and frame-selection kwargs; it
does not implement a custom frame loop. Functions that need per-frame iteration
should use MDAnalysis primitives internally or be replaced by an `AnalysisBase`
subclass.

`MDABackendPolicy` controls optional MDAnalysis internal backends. The default
policy forwards no backend kwargs so PolyzyMD-level parallelism remains the
default. In `comparison.yaml`, the top-level `mda_backend_policy` section maps
to this object and is intentionally opt-in:

```yaml
mda_backend_policy:
  backend: "multiprocessing"
  n_workers: 2
  n_parts: 2
```

Leave the section empty or omit it to forward no backend kwargs. Function-adapter
jobs reject non-default backend policies; use an `AnalysisBase`-compatible job
when opting into MDAnalysis internal parallelism.

```{eval-rst}
.. automodule:: polyzymd.analyses.mda.job
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```
