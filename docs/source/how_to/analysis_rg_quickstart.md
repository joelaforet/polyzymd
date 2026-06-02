# Rg Analysis: Quick Start

Compute Radius of Gyration timeseries to track structural compactness for
protein and polymer selections with autocorrelation-aware uncertainty.

```{versionadded} 1.3.0
The Rg analysis plugin was added in PolyzyMD 1.3.0.
```

```{note}
This page focuses on getting results quickly. For full field-level settings,
output schema details, plot variants, and troubleshooting lookup, see
{doc}`../reference/analysis_rg_reference`.
```

:::{admonition} Environment Setup
:class: tip

All commands below assume you have activated the PolyzyMD pixi environment:

```bash
pixi shell -e build
```

Alternatively, prefix each command with `pixi run -e build`.
:::

## TL;DR

```bash
# Configure Rg runs in comparison.yaml, then run:
polyzymd compare run rg -f comparison.yaml --eq-time 10ns

# Run all enabled analyses in the same workflow
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Force recompute and machine-readable output
polyzymd compare run rg -f comparison.yaml --eq-time 10ns --recompute --format json
```

## Prerequisites

Before running Rg analysis, you need:

1. Completed production simulation data for at least one replicate
2. A `comparison.yaml` file with conditions and plugin settings
3. Trajectory files in the scratch location used by each condition config

Verify your setup:

```bash
ls $(polyzymd info -c config.yaml --scratch-dir)/production_*/
```

## What Rg Analysis Provides

The Rg plugin computes per-run compactness statistics and comparison outputs:

| Feature | Description |
|---------|-------------|
| **Mean Rg** | Average radius of gyration (Å) for the selected atom group |
| **SEM** | Autocorrelation-corrected standard error |
| **Median / Min / Max / Final Rg** | Robust center, range, and endpoint diagnostics |
| **Timeseries** | Full per-frame Rg stored in NPZ sidecars |
| **Multi-run support** | Multiple named selections in one plugin section |

```{tip}
Rg complements RMSD and RMSF:

- **Rg** answers compactness questions
- **RMSD** answers reference-deviation questions
- **RMSF** answers per-residue flexibility questions

Rg is translation and rotation invariant, so it does not require alignment or
reference structures.
```

## Basic Usage

`````{tab-set}

````{tab-item} YAML (Recommended)
Define Rg runs in `comparison.yaml`:

```yaml
name: "rg_quickstart"
control: "no_polymer"

conditions:
  - label: "no_polymer"
    config: "configs/no_polymer.yaml"
    replicates: [1, 2, 3]
  - label: "with_polymer"
    config: "configs/with_polymer.yaml"
    replicates: [1, 2, 3]

plugins:
  rg:
    runs:
      - label: "Whole Protein"
        selection: "protein"
      - label: "Protein Backbone"
        selection: "protein and name CA"
```

Run analysis:

```bash
polyzymd compare run rg -f comparison.yaml --eq-time 10ns
polyzymd compare run-all -f comparison.yaml --eq-time 10ns
polyzymd compare run rg -f comparison.yaml --eq-time 10ns --recompute
```
````

````{tab-item} CLI
Single plugin run:

```bash
polyzymd compare run rg -f comparison.yaml --eq-time 10ns
```

Run all enabled plugins:

```bash
polyzymd compare run-all -f comparison.yaml --eq-time 10ns
```
````

`````

## Multi-Run Configuration

Rg uses a `runs` list. Each run defines a label and selection, and each run is
computed independently for every replicate.

```yaml
plugins:
  rg:
    runs:
      - label: "Whole Protein"
        selection: "protein"

      - label: "Protein Backbone"
        selection: "protein and name CA"

      - label: "Core Region"
        selection: "protein and name CA and resid 20:250"

      - label: "Polymer"
        selection: "chainid C"
```

```{important}
Runs are not replicates.

- A **run** is a named metric definition within the plugin
- A **replicate** is an independent simulation repeat (`run_1`, `run_2`, ...)

All configured runs are evaluated for each available replicate.
```

## Key Mode: Fragment-Aware Rg

```{versionadded} 1.3.0
Fragment-aware Rg calculation was added in PolyzyMD 1.3.0.
```

Use fragment mode when a selection contains many disconnected molecules (for
example, many polymer chains) and you want average fragment compactness rather
than whole-cloud compactness.

```yaml
plugins:
  rg:
    runs:
      - label: "protein_rg"
        selection: "protein"

      - label: "polymer_blob_rg"
        selection: "resname SBM or resname EGM or resname EGP"
        calculation_mode: "fragments"
        fragment_weighting: "equal"
```

Fragment mode details, weighting behavior, and related output fields are
documented in {doc}`../reference/analysis_rg_reference`.

## Comparing Rg Across Conditions

Run condition comparisons with the standard compare command:

```bash
polyzymd compare run rg -f comparison.yaml --eq-time 10ns
```

Per run, PolyzyMD reports:

- Ranking by mean Rg (lower means more compact)
- Pairwise tests with p-values and effect sizes
- Direction labels (`compaction`, `expansion`, `unchanged`)
- ANOVA when 3+ conditions are present

Example output:

```text
Rg Comparison — Whole Protein
===============================
Ranking: With Polymer > No Polymer (lower Rg = more compact)

No Polymer:   18.256 ± 0.044 Å
With Polymer: 17.812 ± 0.038 Å

With Polymer vs No Polymer:
  Change: -2.4% (compaction)
  p-value: 0.0123 *
  Cohen's d: 1.87 (large)
```

For full comparison workflow context, see {doc}`analysis_compare_conditions`.

## Reference and Troubleshooting

For complete lookup material, see {doc}`../reference/analysis_rg_reference`,
including:

- Full `RgRunSettings` and `RgSettings` field tables
- Output directory layout and JSON/NPZ structures
- Plot type details and `plot_settings.rg` options
- CLI option lookup
- Troubleshooting cases and fixes
- Rg vs RMSD and Rg vs RMSF comparison tables

For interpretation guidance, see
{doc}`../explanation/analysis_rg_best_practices` and
{doc}`../explanation/analysis_statistics_best_practices`.

## Next Steps

- {doc}`../reference/analysis_rg_reference`
- {doc}`analysis_compare_conditions`
- {doc}`analysis_rmsd_quickstart`
- {doc}`analysis_rmsf_quickstart`
- {doc}`analysis_distances_quickstart`
- {doc}`analysis_contacts_quickstart`
