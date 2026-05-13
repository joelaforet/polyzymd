# RMSD Analysis: Quick Start

Compute RMSD timeseries for protein and polymer structures with proper
statistical handling in under 5 minutes.

```{versionadded} 1.3.0
The RMSD analysis plugin was added in PolyzyMD 1.3.0.
```

```{note}
**Want to understand the statistics?** This guide focuses on getting results
quickly. For proper uncertainty quantification (autocorrelation correction,
SEM vs. SD) and interpretation of RMSD curves, see the
{doc}`../explanation/analysis_rmsd_best_practices`.
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
# Configure RMSD runs in comparison.yaml, then run:
polyzymd compare run rmsd -f comparison.yaml --eq-time 10ns

# Run all enabled analyses in the same workflow
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Force recompute and machine-readable output
polyzymd compare run rmsd -f comparison.yaml --eq-time 10ns --recompute --format json
```

## Prerequisites

Before running RMSD analysis, you need:

1. **Completed production simulation(s)** — at least one replicate
2. **Comparison config** — `comparison.yaml` with conditions and plugin settings
3. **Trajectory files** — in the scratch directory specified in config

Verify your setup:

```bash
# Check that trajectories exist
ls $(polyzymd info -c config.yaml --scratch-dir)/production_*/
```

## What RMSD Analysis Provides

The RMSD analysis module computes:

| Feature | Description |
|---------|-------------|
| **Mean RMSD** | Average deviation from reference structure (Å) |
| **SEM** | Autocorrelation-corrected standard error of the mean |
| **Median RMSD** | Robust central tendency measure |
| **Min / Max RMSD** | Extremes of conformational deviation |
| **Final RMSD** | Last-frame RMSD (convergence diagnostic) |
| **Timeseries** | Full per-frame RMSD saved as NPZ sidecar |
| **Multi-run** | Multiple named selections in a single analysis |
| **Convergence Detection** | Sliding-window slope diagnostic; detects when RMSD has plateaued |

```{tip}
**RMSD vs RMSF vs Distances — when to use which:**
- **RMSD**: Global structural deviation over time — "is the protein drifting?"
- **RMSF**: Per-residue fluctuation around average — "which residues are flexible?"
- **Distances**: Specific atom-pair distances — "is this H-bond intact?"
```

## Basic Usage

`````{tab-set}

````{tab-item} YAML (Recommended)
For reproducible analysis, define RMSD runs in `comparison.yaml`:

```yaml
# comparison.yaml
name: "rmsd_quickstart"
control: "no_polymer"

conditions:
  - label: "no_polymer"
    config: "configs/no_polymer.yaml"
    replicates: [1, 2, 3]
  - label: "with_polymer"
    config: "configs/with_polymer.yaml"
    replicates: [1, 2, 3]

plugins:
  rmsd:
    runs:
      - label: "Protein Backbone"
        selection: "protein and name CA"
        alignment_selection: "protein and name CA"
        reference_mode: "centroid"
```

Then run:

```bash
# Run RMSD analysis only
polyzymd compare run rmsd -f comparison.yaml --eq-time 10ns

# Run all enabled plugins in comparison.yaml
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Force recompute
polyzymd compare run rmsd -f comparison.yaml --eq-time 10ns --recompute
```

**Benefits:**
- Version-controlled, reproducible
- Self-documenting experiment setup
- Easy to re-run with different parameters
````

````{tab-item} CLI
### Single analysis run

```bash
polyzymd compare run rmsd -f comparison.yaml --eq-time 10ns
```

**Expected behavior:**

```text
Loading comparison config from: comparison.yaml
Running plugin: rmsd
  Equilibration: 10ns
  Conditions: no_polymer, with_polymer
  Runs: Protein Backbone
RMSD comparison complete
```

### All enabled analyses

Run RMSD plus any other enabled plugins:

```bash
polyzymd compare run-all -f comparison.yaml --eq-time 10ns
```
````

`````

## Multi-Run Configuration

The RMSD plugin uses a **runs** list, where each run defines a named RMSD
calculation with its own selection, alignment, and reference settings. This
lets you track multiple structural metrics in a single analysis pass.

```yaml
plugins:
  rmsd:
    runs:
      - label: "Protein Backbone"
        selection: "protein and name CA"
        alignment_selection: "protein and name CA"
        reference_mode: "centroid"

      - label: "Active Site"
        selection: "protein and (resid 77 or resid 133 or resid 156) and name CA"
        alignment_selection: "protein and name CA"
        reference_mode: "centroid"

      - label: "Polymer Core"
        selection: "chainID C and not name H*"
        alignment_selection: "protein and name CA"
        reference_mode: "average"
```

Each run produces an independent RMSD timeseries. During comparison, each run
is ranked and tested separately — averaging RMSD from different selections is
not meaningful.

```{important}
**Runs ≠ Replicates.** A "run" is a named RMSD selection (e.g., "Protein
Backbone" vs "Active Site"). A "replicate" is an independent simulation repeat
(run_1, run_2, run_3). All configured runs are computed for every replicate.
```

## External Reference Mode

Use `reference_mode: "external"` when you want to measure deviation from a
specific known structure, such as a crystal structure representing the
catalytically competent geometry.

```yaml
plugins:
  rmsd:
    runs:
      - label: "Crystal Deviation"
        selection: "protein and name CA"
        alignment_selection: "protein and name CA"
        reference_mode: "external"
        reference_file: "/path/to/crystal_structure.pdb"
```

```{note}
When using `external` reference mode, the external PDB must contain atoms
matching the `selection` string. PolyzyMD validates that atom counts match
between the trajectory and external reference and raises an error on mismatch.
```

**When to use external reference:**

| Mode | Question Answered |
|------|-------------------|
| `centroid` (default) | How much does the structure deviate from its most populated conformation? |
| `average` | How much does the structure deviate from its time-averaged conformation? |
| `frame` | How much does the structure deviate from a specific frame? |
| `external` | How much does the structure deviate from a known functional geometry? |

```{tip}
For enzyme studies, consider running **two RMSD runs**: one with `centroid`
mode (overall stability) and one with `external` mode pointing to a crystal
structure (catalytic competence). These answer complementary questions.
```

## Comparing RMSD Across Conditions

To statistically compare RMSD across multiple simulation conditions (e.g.,
different polymer compositions), use the `compare run rmsd` command:

```bash
# Add rmsd section to comparison.yaml, then:
polyzymd compare run rmsd -f comparison.yaml --eq-time 10ns
```

This provides **per-run**:
- **Ranking**: Conditions sorted by mean RMSD (lowest = most stable)
- **Pairwise t-tests**: With p-values, Cohen's d, percent change
- **Direction labels**: `stabilizing` (lower RMSD), `destabilizing` (higher), or `unchanged`
- **ANOVA**: Omnibus test when 3+ conditions are present

**Example output:**

```text
RMSD Comparison — Protein Backbone
===================================
Ranking: With Polymer > No Polymer (lower RMSD = more stable)

No Polymer:   1.856 ± 0.034 Å
With Polymer: 1.612 ± 0.028 Å

With Polymer vs No Polymer:
  Change: -13.1% (stabilizing)
  p-value: 0.0089 *
  Cohen's d: 2.41 (large)
```

See {doc}`analysis_compare_conditions` for the full multi-plugin comparison
workflow.

## Reference and Troubleshooting

For the full list of configuration fields, default values, output file
structure, plotting options, convergence details, CLI options, and
troubleshooting fixes, see {doc}`../reference/analysis_rmsd_reference`.

For deeper interpretation guidance, see {doc}`../explanation/analysis_rmsd_best_practices`
and {doc}`../explanation/convergence_detection`.

## Next Steps

- **Understand RMSD interpretation**: {doc}`../explanation/analysis_rmsd_best_practices`
- **Compare conditions**: {doc}`analysis_compare_conditions`
- **RMSF analysis**: {doc}`analysis_rmsf_quickstart`
- **Understand statistics**: {doc}`../explanation/analysis_statistics_best_practices`
- **Distance analysis**: {doc}`analysis_distances_quickstart`
- **Contact analysis**: {doc}`analysis_contacts_quickstart`
