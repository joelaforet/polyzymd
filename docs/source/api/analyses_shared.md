# Analysis Shared Utilities

This reference page documents contributor-facing utilities in
`polyzymd.analyses.shared`. These modules provide reusable building blocks for
the shipped analysis functions and your own.

The package root re-exports common helpers for convenience. Import grouping
classes and module-specific helpers from their submodules.

## Trajectory loading and windows

Use these modules to locate trajectories, parse time values, and resolve the
trajectory window of each replicate.

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.loader
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

.. automodule:: polyzymd.analyses.shared.gromacs
   :members:
   :no-index:

.. automodule:: polyzymd.analyses.shared.window
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:
```

## Alignment and representative frames

Alignment helpers in `polyzymd.analyses.shared.alignment` standardize
reference-mode handling. Centroid helpers support analyses that need
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
and inferential tests used by the study API's summaries and comparisons and by
the analysis functions.

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

## Selections and residue grouping

Selection helpers extend MDAnalysis selections with midpoints and centres of
mass. The grouping and amino-acid classification modules classify protein
residues (for example aromatic, polar, charged) and give the maximum solvent
accessible surface area of each residue type.

```{eval-rst}
.. automodule:: polyzymd.analyses.shared.selections
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

## Diagnostics

Diagnostics helpers validate selections and the equilibration window. The
module is `polyzymd.analyses.shared.diagnostics`.
