# Explanation

Use Explanation when you need the reasoning behind PolyzyMD: design choices,
analysis assumptions, statistical interpretation, and plugin-specific caveats.
If you want commands for a task, start with {doc}`../how_to/index`. If you need
to look up options, file layouts, or schemas, use {doc}`../reference/index`.

```{tip}
New to the analysis plugin system? Start with {doc}`analysis_concepts`, then
read the statistical and reference-selection pages listed in the contributor
pathway below.
```

::::{grid} 2
:gutter: 3

:::{grid-item-card} New contributor path
:link: analysis_concepts
:link-type: doc

Read the concepts, statistics, convergence, and reference-selection pages that
shape analysis-plugin design.
:::

:::{grid-item-card} Protein modifications
:link: protein_modification_architecture
:link-type: doc

How covalent residue modifications and polymer conjugation are represented.
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

If you are adding or reviewing an analysis plugin, read these pages first:

1. {doc}`analysis_concepts` for the analysis lifecycle and plugin boundaries.
2. {doc}`analysis_statistics_best_practices` for replicate-level interpretation
   and statistical expectations.
3. {doc}`convergence_detection` for deciding whether trajectory summaries are
   interpretable.
4. {doc}`analysis_reference_selection` for how structural references affect
   RMSF and related fluctuation metrics.

## Concepts and Design

```{toctree}
:maxdepth: 1

Analysis system concepts <analysis_concepts>
Architecture and design rationale <architecture>
Protein modification architecture <protein_modification_architecture>
Residue assignment and chain conventions <residue_assignment>
Why PolyzyMD uses colored logging <colored_logging>
```

## Interpretation and Best Practices

Start with the cross-cutting interpretation pages, then use the plugin-specific
pages for caveats tied to particular metrics.

### Foundations

```{toctree}
:maxdepth: 1

Statistics best practices for MD analysis <analysis_statistics_best_practices>
Establishing convergence in MD simulations <convergence_detection>
```

### Metric and plugin caveats

```{toctree}
:maxdepth: 1

RMSD interpretation: use, limits, and cautions <analysis_rmsd_best_practices>
Rg analysis: best practices <analysis_rg_best_practices>
RMSF analysis: statistical best practices <analysis_rmsf_best_practices>
How RMSF reference selection changes interpretation <analysis_reference_selection>
RMSF implementation verification <analysis_rmsf_verification>
Catalytic triad analysis: interpretation and best practices <analysis_triad_best_practices>
```
