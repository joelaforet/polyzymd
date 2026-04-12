# Tutorial: Measure Surface Accessibility with the SASA Plugin

This tutorial walks you through using the SASA (Solvent Accessible Surface
Area) analysis plugin to measure how polymer conjugation affects enzyme surface
exposure. By the end, you will have configured a multi-run SASA analysis,
executed it, and interpreted the comparison results.

## What You Will Learn

- How to configure multiple SASA "runs" with different target and context
  selections
- How the target/context model lets you isolate polymer shielding effects
- How to run SASA analysis locally and on HPC
- How to read the comparison output and identify shielding signals in the plots

## Prerequisites

Before starting, make sure you have:

- A working pixi environment (`pixi install -e build`)
- Completed simulation trajectories for at least two conditions
- A `comparison.yaml` defining your conditions (see
  {doc}`../how_to/analysis_compare_conditions`)
- Familiarity with the PolyzyMD chain convention (A=protein, B=substrate,
  C=polymer)

If you have not run a basic analysis yet, complete {doc}`first_analysis` first.

## How the SASA Plugin Works

The SASA plugin computes Solvent Accessible Surface Area using the
Shrake-Rupley algorithm via MDAnalysis. It reports per-frame total SASA and
per-residue SASA profiles.

The key feature is the **multi-run model**. Instead of computing one number,
you define multiple "runs" — each with its own **target** and **context**
selections:

| Selection | What it controls |
|-----------|-----------------|
| `target_selection` | Atoms whose SASA is **reported** |
| `context_selection` | Atoms that are **considered as blocking surface** during the calculation |

For example, setting `target_selection: "protein"` and
`context_selection: "protein or chainID C"` reports the protein's SASA while
treating the polymer (chain C) as an obstruction. Comparing this to a run
where `context_selection: "protein"` (polymer ignored) tells you how much
surface area the polymer covers.

## Step 1: Configure SASA Runs

Add a `sasa` section to the `plugins:` block in your `comparison.yaml`. Each
run defines one SASA calculation with its own target and context.

Here is a four-run configuration that fully characterizes polymer shielding:

```yaml
plugins:
  sasa:
    runs:
      - label: "protein_isolated"
        target_selection: "protein"
        context_selection: "protein"

      - label: "protein_with_polymer"
        target_selection: "protein"
        context_selection: "protein or chainID C"

      - label: "active_site_isolated"
        target_selection: "protein and (resid 77 or resid 156 or resid 262)"
        context_selection: "protein"

      - label: "active_site_with_polymer"
        target_selection: "protein and (resid 77 or resid 156 or resid 262)"
        context_selection: "protein or chainID C"

    probe_radius_nm: 0.14
    n_sphere_points: 960
```

### What each run tells you

| Run | Question it answers |
|-----|---------------------|
| `protein_isolated` | What is the total protein SASA when only protein atoms block solvent? (baseline) |
| `protein_with_polymer` | What is the protein SASA when the polymer can also block solvent? (shielded) |
| `active_site_isolated` | How exposed is the active site without polymer effects? |
| `active_site_with_polymer` | How exposed is the active site when the polymer is present? |

The difference between `protein_isolated` and `protein_with_polymer` is the
polymer's shielding effect. A decrease in SASA indicates that the polymer
covers part of the protein surface.

:::{tip}
For the "No Polymer" control condition, the `protein_with_polymer` run will
produce the same result as `protein_isolated` because there is no chain C. This
is expected and provides a useful internal consistency check.
:::

### SASA algorithm parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `probe_radius_nm` | `0.14` | Shrake-Rupley probe radius in nm (standard water-sized probe) |
| `n_sphere_points` | `960` | Number of points on each atom's test sphere (higher = more accurate, slower) |
| `chunk_size` | `100` | Frames processed per chunk (lower = less memory, slower) |

The defaults are suitable for most enzyme-polymer systems. Increase
`n_sphere_points` to 1500+ only if you need very high precision for
small SASA differences.

## Step 2: The Full comparison.yaml

Here is a complete example for a CALB enzyme study:

```yaml
name: "calb_sasa_study"
description: "SASA analysis for CALB with polymer conjugates"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../noPoly_CALB_pNPB/config.yaml"
    replicates: [1, 2, 3]

  - label: "SBMA-100"
    config: "../SBMA_100_CALB_pNPB/config.yaml"
    replicates: [1, 2, 3]

  - label: "EGMA-100"
    config: "../EGMA_100_CALB_pNPB/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  sasa:
    runs:
      - label: "protein_isolated"
        target_selection: "protein"
        context_selection: "protein"
      - label: "protein_with_polymer"
        target_selection: "protein"
        context_selection: "protein or chainID C"
      - label: "active_site_isolated"
        target_selection: "protein and (resid 77 or resid 156 or resid 262)"
        context_selection: "protein"
      - label: "active_site_with_polymer"
        target_selection: "protein and (resid 77 or resid 156 or resid 262)"
        context_selection: "protein or chainID C"
    probe_radius_nm: 0.14
    n_sphere_points: 960

plot_settings:
  format: "png"
  dpi: 300
  style: "publication"
```

## Step 3: Run Locally

For small systems or interactive exploration, run the SASA analysis locally:

```bash
pixi run -e build polyzymd compare run sasa \
    --comparison-yaml comparison.yaml
```

This runs the full pipeline sequentially: `compute_replicate` for every
replicate, `aggregate` for every condition, then `compare` and `plot`. Expect
this to take several minutes per replicate depending on trajectory length and
system size.

:::{note}
SASA computation is CPU-intensive and memory-intensive because the
Shrake-Rupley algorithm iterates over every atom in the context selection for
every frame. The `chunk_size` parameter controls how many frames are loaded
into memory at once.
:::

## Step 4: Run on HPC

For large systems or many replicates, submit the analysis as SLURM jobs:

```bash
pixi run -e build polyzymd compare submit sasa \
    --comparison-yaml comparison.yaml \
    --partition aa100 \
    --mem 8G \
    --time 02:00:00
```

This submits a DAG of SLURM jobs that process replicates in parallel. See
{doc}`../how_to/hpc_execution` for the full HPC guide, including dry runs, monitoring,
and troubleshooting.

:::{tip}
SASA jobs are marked with `execution_cost_hint = "high"` in the plugin.
Allocate at least 8 GB of memory and 1–2 hours of wall time per replicate
for systems with 50,000+ atoms.
:::

## Step 5: Interpret the Results

### Comparison JSON

The comparison result is saved to
`comparison/sasa/result.json`. It contains:

- **Per-condition summaries** with mean SASA ± SEM for each run
- **Pairwise comparisons** between conditions for each run (t-test, Cohen's d,
  percent change)
- **ANOVA results** for each run when three or more conditions are compared
- **Rankings** of conditions by SASA for each run
- **Direction labels**: `"shielding"` (SASA decreased >1%), `"exposure"` (SASA
  increased >1%), or `"unchanged"`

Here is a simplified excerpt showing how to read the key fields:

```json
{
  "run_labels": ["protein_isolated", "protein_with_polymer"],
  "conditions": [
    {
      "label": "No Polymer",
      "run_summaries": [
        {
          "label": "protein_isolated",
          "mean_sasa": 12450.3,
          "sem_sasa": 85.2
        },
        {
          "label": "protein_with_polymer",
          "mean_sasa": 12448.1,
          "sem_sasa": 84.9
        }
      ]
    },
    {
      "label": "SBMA-100",
      "run_summaries": [
        {
          "label": "protein_isolated",
          "mean_sasa": 12380.5,
          "sem_sasa": 91.0
        },
        {
          "label": "protein_with_polymer",
          "mean_sasa": 11200.7,
          "sem_sasa": 102.3
        }
      ]
    }
  ],
  "pairwise_comparisons": [
    {
      "run_label": "protein_with_polymer",
      "condition_a": "No Polymer",
      "condition_b": "SBMA-100",
      "p_value": 0.003,
      "cohens_d": -1.82,
      "direction": "shielding",
      "significant": true,
      "percent_change": -10.0
    }
  ]
}
```

### How to read the results

1. **Compare `protein_isolated` across conditions.** This measures intrinsic
   protein compactness differences (without polymer effects). Small differences
   here indicate that the protein folds similarly across conditions.

2. **Compare `protein_with_polymer` across conditions.** The "No Polymer"
   condition serves as the baseline. Conditions with polymer show lower SASA
   in this run if the polymer shields the surface.

3. **Calculate the shielding effect.** For each polymer condition, subtract
   `protein_with_polymer` from `protein_isolated`. The larger the difference,
   the more surface the polymer covers.

4. **Check active site runs.** If `active_site_with_polymer` is significantly
   lower than `active_site_isolated` for a polymer condition, the polymer may
   be blocking substrate access — a concern for enzyme activity.

### Plots

The SASA plugin generates three types of plots:

| Plot | File pattern | What it shows |
|------|-------------|---------------|
| **Comparison bars** | `sasa_comparison_<run>.png` | Mean SASA ± SEM per condition, with replicate scatter points |
| **Time series** | `sasa_timeseries_<run>.png` | Per-frame SASA traces overlaid for each condition |
| **Residue profiles** | `sasa_profile_<run>.png` | Per-residue mean SASA across conditions |

The bar plots are the most informative for quick assessment. Look for
conditions where the `protein_with_polymer` bar is significantly lower than
the `protein_isolated` bar — this is the polymer shielding signal.

## Common Configurations

:::{admonition} Recipe collection (how-to mode)
:class: note

The configurations below are **task-oriented recipes** rather than step-by-step
tutorial content. Use them as starting points for your own SASA analysis.
:::

### Minimal: whole-protein SASA only

```yaml
plugins:
  sasa:
    runs:
      - label: "protein_total"
        target_selection: "protein"
```

When `context_selection` is omitted, it defaults to match `target_selection`.
This measures the protein's self-SASA without considering any other molecules.

### Two-run: basic shielding comparison

```yaml
plugins:
  sasa:
    runs:
      - label: "protein_only"
        target_selection: "protein"
        context_selection: "protein"
      - label: "protein_full_context"
        target_selection: "protein"
        context_selection: "all"
```

Using `"all"` as the context includes everything (protein, polymer, substrate,
solvent). This gives the "true" SASA but makes comparisons harder because
solvent box size differences between conditions can affect the result.

### Monomer-specific shielding

If your polymer contains specific monomer types, you can test which monomer
contributes more to shielding:

```yaml
plugins:
  sasa:
    runs:
      - label: "protein_isolated"
        target_selection: "protein"
        context_selection: "protein"
      - label: "protein_with_sbma"
        target_selection: "protein"
        context_selection: "protein or resname SBMA"
      - label: "protein_with_egma"
        target_selection: "protein"
        context_selection: "protein or resname EGMA"
```

:::{warning}
Ensure your monomer residue names (`resname SBMA`, `resname EGMA`) match
the actual residue names in your topology. Check with:

```bash
pixi run -e build python -c "
import MDAnalysis as mda
u = mda.Universe('solvated_system.pdb')
print(set(u.select_atoms('chainID C').residues.resnames))
"
```
:::

### Active site focus with specific residues

```yaml
plugins:
  sasa:
    runs:
      - label: "catalytic_triad"
        target_selection: "protein and (resid 77 or resid 156 or resid 262)"
        context_selection: "protein or chainID C"
      - label: "binding_pocket"
        target_selection: "protein and (resid 77 or resid 156 or resid 262 or resid 80 or resid 155)"
        context_selection: "protein or chainID C"
```

### Stride for long trajectories

If your trajectory has many frames and SASA computation is slow, increase
the stride to analyze every Nth frame:

```yaml
plugins:
  sasa:
    runs:
      - label: "protein_with_polymer"
        target_selection: "protein"
        context_selection: "protein or chainID C"
        stride: 5
```

:::{note}
A stride of 5 uses every 5th frame, reducing computation time by ~5x. The
plugin still computes autocorrelation-corrected SEM on the subsampled data.
:::

## Result Models

:::{admonition} Reference material
:class: note

This section is **reference-style** content for plugin developers and advanced
users who need to inspect the data models programmatically.
:::

For plugin developers or advanced users, the SASA result hierarchy is:

| Model | Level | Key fields |
|-------|-------|------------|
| `SASARunResult` | Per-replicate, per-run | `mean_sasa`, `std_sasa`, `sem_sasa`, `n_frames_used`, `n_target_atoms` |
| `SASAResult` | Per-replicate (all runs) | `run_results: list[SASARunResult]` |
| `SASARunAggregatedResult` | Per-condition, per-run | `overall_mean`, `overall_sem`, `per_replicate_means`, `per_residue_mean_sasa` |
| `SASAAggregatedResult` | Per-condition (all runs) | `run_results: list[SASARunAggregatedResult]` |
| `SASAComparisonResult` | Cross-condition | `pairwise_comparisons`, `ranking_by_run`, `anova_by_run` |

Raw per-frame and per-residue SASA data is stored as NPZ sidecars alongside
the JSON result files.

## What You Have Now

After following this tutorial, you have:

- a multi-run SASA configuration that isolates the polymer shielding effect
- comparison results with per-run statistical tests across conditions
- bar, time series, and residue-profile plots for visual assessment
- the understanding to design custom SASA run configurations for your system

## See Also

- {doc}`../how_to/hpc_execution` — Submitting analysis jobs to SLURM
- {doc}`../how_to/analysis_compare_conditions` — Setting up comparison.yaml
- {doc}`../contributor_guide/extending_analyses` — Writing your own analysis plugin
- {doc}`../explanation/analysis_statistics_best_practices` — Autocorrelation and uncertainty
