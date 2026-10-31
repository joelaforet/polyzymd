# Explanation

Explanation pages give the reasons behind PolyzyMD's design and the
assumptions behind each analysis. For the steps of a task, see
{doc}`../how_to/index`. To look up a key or an option, see
{doc}`../reference/index`.

Read these pages in this order to understand an analysis result:

1. {doc}`analysis_concepts`: the steps of an analysis and what PolyzyMD
   stores.
2. {doc}`analysis_statistics_best_practices`: why PolyzyMD uses one value per
   replicate, and how it compares conditions.
3. {doc}`convergence_detection`: how PolyzyMD checks that a series has
   reached equilibrium.
4. The page of the metric you report, under "Metrics" below.

## Studies and projects

```{toctree}
:caption: Studies and projects
:maxdepth: 1

projects
study_folders
```

## Analysis

```{toctree}
:caption: Analysis
:maxdepth: 1

analysis_concepts
analysis_api
analysis_entry_points
analysis_statistics_best_practices
convergence_detection
```

## Metrics

```{toctree}
:caption: Metrics
:maxdepth: 1

analysis_rmsd_best_practices
analysis_rg_best_practices
analysis_rmsf_best_practices
analysis_reference_selection
analysis_contact_lifetimes
analysis_triad_best_practices
```

## Simulation design

```{toctree}
:caption: Simulation design
:maxdepth: 1

simulation_safeguards
physics_audit
modified_residues
residue_assignment
gromacs_parallelism
```
