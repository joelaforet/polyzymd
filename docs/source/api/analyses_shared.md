# Analysis Shared Utilities

This reference page documents contributor-facing utilities in
`polyzymd.analyses.shared`. These modules provide reusable building blocks for
analysis plugins; framework internals and plugin-private helpers are documented
with their owning packages.

The package root re-exports common helpers for convenience. Import specialized
selectors, grouping classes, and module-specific helpers from their submodules.

## Trajectory loading and windows

Use these modules to locate trajectories, parse time values, and resolve the
trajectory window passed into MDAnalysis job lifecycles.

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.loader
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.window
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Alignment and representative frames

Alignment helpers in `polyzymd.analyses.shared.alignment` standardize
reference-mode handling. Centroid helpers support plugins that need
representative frames or structures.

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.centroid
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Time-series statistics and convergence

These modules provide statistical summaries, autocorrelation-aware estimates,
inferential tests, and convergence diagnostics used by built-in and contributor
plugins.

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.autocorrelation
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.statistics
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.inferential_statistics
   :members:
   :exclude-members: cohens_d
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.convergence
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Plotting

Plotting helpers centralize figure themes, output paths, axis styling, legends,
grouped bars, and matrix annotations.

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.plotting
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Selections, selectors, and grouping

Selection helpers extend MDAnalysis selections. Selector and grouping packages
provide reusable abstractions for selecting molecules or classifying residues in
plugin settings and analysis code.

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.selections
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.selectors
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.selectors.base
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.selectors.protein
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.selectors.polymer
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.selectors.solvent
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.groupings
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.groupings.base
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.aa_classification
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Diagnostics and path helpers

Diagnostics helpers validate selections and analysis inputs. The module is
`polyzymd.analyses.shared.diagnostics`. Path helpers standardize
artifact-oriented file locations used by analysis plugins.

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.paths
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Multi-run comparison and formatting

Multi-run helpers support plugins that compare several named runs or entities
per condition, such as RMSD, radius of gyration, and SASA analyses.

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.multi_run_comparison
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.multi_run_formatting
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Compatibility facades

`polyzymd.analyses.shared.config_hash` and `polyzymd.analyses.shared.sasa`
remain importable as temporary compatibility facades for older callers. They are
not preferred contributor extension points, and this page intentionally does not
include autodoc blocks for them.
