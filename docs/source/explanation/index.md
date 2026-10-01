# Explanation

Use Explanation when you need the reasoning behind PolyzyMD: design choices,
analysis assumptions, statistical interpretation, and metric-specific caveats.
If you want commands for a task, start with {doc}`../how_to/index`. If you need
to look up options, file layouts, or schemas, use {doc}`../reference/index`.

```{tip}
New to PolyzyMD analysis? Start with {doc}`analysis_concepts`, then
read the statistical and reference-selection pages listed in the contributor
pathway below.
```

::::{grid} 2
:gutter: 3

:::{grid-item-card} New contributor path
:link: analysis_concepts
:link-type: doc

Read the concepts, statistics, convergence, and reference-selection pages that
shape how an analysis is designed.
:::

:::{grid-item-card} Statistical interpretation
:link: analysis_statistics_best_practices
:link-type: doc

Understand replication, uncertainty, multiple comparisons, and convergence
before comparing simulation conditions.
:::

:::{grid-item-card} RMSF interpretation
:link: analysis_rmsf_best_practices
:link-type: doc

Learn why RMSF depends on alignment, residue mapping, replicate treatment, and
the selected structural reference.
:::

:::{grid-item-card} Architecture and conventions
:link: architecture
:link-type: doc

See the design rationale and implementation constraints that keep PolyzyMD
analyses consistent. Chain conventions are covered in the residue-assignment
page.
:::

::::

## New contributor pathway

If you are adding or reviewing an analysis, read these pages first:

1. {doc}`analysis_concepts` for the analysis steps and what is stored.
2. {doc}`analysis_statistics_best_practices` for replicate-level interpretation
   and statistical expectations.
3. {doc}`convergence_detection` for deciding whether trajectory summaries are
   interpretable.
4. {doc}`analysis_reference_selection` for how structural references affect
   RMSF and related fluctuation metrics.

## Concepts and Design

```{toctree}
:maxdepth: 1

Analysing a set of simulations <analysis_api>
Study folders: publishing a reproducible MD study <study_folders>
Analysis system concepts <analysis_concepts>
Which analysis entry point should I use <analysis_entry_points>
Architecture and design rationale <architecture>
Residue assignment and chain conventions <residue_assignment>
Why PolyzyMD guards the prepared structure <simulation_safeguards>
Why PolyzyMD uses colored logging <colored_logging>
```

## Interpretation and Best Practices

Start with the cross-cutting interpretation pages, then use the metric-specific
pages for caveats tied to particular metrics.

### Foundations

```{toctree}
:maxdepth: 1

Statistics best practices for MD analysis <analysis_statistics_best_practices>
Establishing convergence in MD simulations <convergence_detection>
Methods and references <references>
```

### Metric caveats

```{toctree}
:maxdepth: 1

RMSD interpretation: use, limits, and cautions <analysis_rmsd_best_practices>
Rg analysis: best practices <analysis_rg_best_practices>
RMSF analysis: statistical best practices <analysis_rmsf_best_practices>
How RMSF reference selection changes interpretation <analysis_reference_selection>
RMSF implementation verification <analysis_rmsf_verification>
SASA implementation verification <analysis_sasa_verification>
Contacts implementation verification <analysis_contacts_verification>
How long polymer contacts last <analysis_contact_lifetimes>
Native contacts implementation verification <analysis_native_contacts_verification>
Hydrogen bonds implementation verification <analysis_hydrogen_bonds_verification>
Catalytic triad: interpretation and best practices <analysis_triad_best_practices>
```
