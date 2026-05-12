# Explanation

Explanation pages help you understand *why* PolyzyMD works the way it does, how
to interpret outputs, and what assumptions or tradeoffs matter.

```{tip}
Looking for step-by-step instructions? Go to {doc}`../how_to/index`.
Need to look up a setting or option? Go to {doc}`../reference/index`.
```

::::{grid} 2
:gutter: 3

:::{grid-item-card} Concepts and Design
:link: analysis_concepts
:link-type: doc

How the analysis pipeline, chain conventions, and plugin architecture fit
together.
:::

:::{grid-item-card} Protein Modifications
:link: protein_modification_architecture
:link-type: doc

How covalent residue modifications and polymer conjugation are represented.
:::

:::{grid-item-card} Interpretation and Best Practices
:link: analysis_statistics_best_practices
:link-type: doc

Replication, convergence, multiple testing, and per-plugin best practices.
:::

:::{grid-item-card} Experimental Analysis Concepts
:link: polymer_bridging_interpretation
:link-type: doc

Interpretation guidance for experimental analysis workflows.
:::

::::

## Concepts and Design

```{toctree}
:maxdepth: 1

Analysis System Concepts <analysis_concepts>
Architecture <architecture>
Protein Modification Architecture <protein_modification_architecture>
Residue Assignment and Chain Conventions <residue_assignment>
Colored Logging <colored_logging>
```

## Interpretation and Best Practices

```{toctree}
:maxdepth: 1

Statistical Best Practices for Analysis <analysis_statistics_best_practices>
RMSD Best Practices <analysis_rmsd_best_practices>
Convergence Detection <convergence_detection>
Rg Best Practices <analysis_rg_best_practices>
RMSF Best Practices <analysis_rmsf_best_practices>
RMSF Reference Selection <analysis_reference_selection>
RMSF Verification <analysis_rmsf_verification>
Catalytic Triad Best Practices <analysis_triad_best_practices>
```

## Experimental Analysis Concepts

These workflows remain available in PolyzyMD, but their interpretation is still
evolving and they should be treated as experimental. See also the how-to guides
for {doc}`../how_to/analysis_binding_preference`,
{doc}`../how_to/analysis_binding_free_energy`,
{doc}`../how_to/analysis_polymer_affinity`, and
{doc}`../how_to/analysis_exposure_dynamics`.

```{toctree}
:maxdepth: 1

Polymer Bridging Interpretation <polymer_bridging_interpretation>
```
