# Rg Analysis: Quick Start

Compute radius of gyration timeseries to track structural compactness for
protein and polymer selections, with uncertainty taken across replicates.

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

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
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

The Rg plugin reports one observable per run, and the framework turns it into
condition-level statistics:

| Feature | Description |
|---------|-------------|
| Mean Rg | Mean over replicates of each replicate's mean radius of gyration (Å) |
| SEM and 95 percent interval | Taken across replicates, with the replicate as the sampling unit |
| Timeseries | Full per-frame Rg stored in the `observables.npz` sidecar |
| Fragment profile | In fragment mode, the mean Rg of each bonded fragment and their distribution |
| Multi-run support | Multiple named selections in one plugin section |

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

Fragment mode needs a topology that carries bonds, because fragments are
connected components of the bond graph. If the selected atoms have no bonds, the
run fails with `TopologyBondsMissingError` instead of quietly measuring the
whole selection as one fragment. This is common for solvated systems above
99999 atoms, where MDAnalysis skips the CONECT records that OpenMM writes in
hexadecimal. Load a topology that carries bonds, or guess bonds for the protein
and polymer selection. If whole-selection Rg is genuinely what you want, set
`allow_single_fragment_fallback: true` on the run and say so when you report the
number.

Fragment mode details, weighting behavior, and related output fields are
documented in {doc}`../reference/analysis_rg_reference`.

## Comparing Rg Across Conditions

Run condition comparisons with the standard compare command:

```bash
polyzymd compare run rg -f comparison.yaml --eq-time 10ns
```

Per observable, PolyzyMD reports the mean and 95 percent interval of each
condition and a test of every condition against the control, with one
Benjamini-Hochberg family over the whole run.

Example output:

```text
# rg  eq 10ns
No Polymer    rg_whole_protein  mean_of_timeseries  mean 18.26 A  sem 0.044  ci95 18.07 to 18.45  n 3
With Polymer  rg_whole_protein  mean_of_timeseries  mean 17.81 A  sem 0.038  ci95 17.65 to 17.98  n 3
No Polymer vs With Polymer  rg_whole_protein  delta -0.444  p_adj 0.0123  test student_t  correction benjamini_hochberg  significant
```

For full comparison workflow context, see {doc}`analysis_compare_conditions`.

## Reference and Troubleshooting

For complete lookup material, see {doc}`../reference/analysis_rg_reference`,
including:

- Full `RgRunSettings` and `RgSettings` field tables
- The observables each run reports, with kinds and units
- Output directory layout and artifact payloads
- Error messages and their fixes
- Rg against RMSD and RMSF

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
