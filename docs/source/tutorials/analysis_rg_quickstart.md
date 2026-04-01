# Rg Analysis: Quick Start

Compute Radius of Gyration timeseries measuring structural compactness for
protein and polymer selections with proper statistical handling.

```{versionadded} 1.3.0
The Rg analysis plugin was added in PolyzyMD 1.3.0.
```

```{note}
**Want to understand the statistics?** This guide focuses on getting results
quickly. For proper uncertainty quantification (autocorrelation correction,
SEM vs. SD) and interpretation of Rg curves, see the
[Rg Best Practices Guide](analysis_rg_best_practices.md).
```

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

1. **Completed production simulation(s)** — at least one replicate
2. **Comparison config** — `comparison.yaml` with conditions and plugin settings
3. **Trajectory files** — in the scratch directory specified in config

Verify your setup:

```bash
# Check that trajectories exist
ls $(polyzymd info -c config.yaml --scratch-dir)/production_*/
```

## What Rg Analysis Provides

The Rg analysis module computes:

| Feature | Description |
|---------|-------------|
| **Mean Rg** | Average radius of gyration (Å) — measures structural compactness |
| **SEM** | Autocorrelation-corrected standard error of the mean |
| **Median Rg** | Robust central tendency measure |
| **Min / Max Rg** | Extremes of compactness/expansion |
| **Final Rg** | Last-frame Rg (convergence diagnostic) |
| **Timeseries** | Full per-frame Rg saved as NPZ sidecar |
| **Multi-run** | Multiple named selections in a single analysis |

```{tip}
**Rg vs RMSD vs RMSF — when to use which:**
- **Rg**: Overall structural compactness — "is the protein expanding or compacting?"
- **RMSD**: Global structural deviation from a reference — "is the protein drifting?"
- **RMSF**: Per-residue fluctuation around average — "which residues are flexible?"

Rg does **not** require alignment or a reference structure — it is inherently
translation and rotation invariant. This makes it a complementary metric to
RMSD for assessing conformational stability.
```

## Basic Usage

`````{tab-set}

````{tab-item} YAML (Recommended)
For reproducible analysis, define Rg runs in `comparison.yaml`:

```yaml
# comparison.yaml
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
    enabled: true
    runs:
      - label: "Whole Protein"
        selection: "protein"
      - label: "Protein Backbone"
        selection: "protein and name CA"
```

Then run:

```bash
# Run Rg analysis only
polyzymd compare run rg -f comparison.yaml --eq-time 10ns

# Run all enabled plugins in comparison.yaml
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Force recompute
polyzymd compare run rg -f comparison.yaml --eq-time 10ns --recompute
```

**Benefits:**
- Version-controlled, reproducible
- Self-documenting experiment setup
- Easy to re-run with different parameters
````

````{tab-item} CLI
### Single analysis run

```bash
polyzymd compare run rg -f comparison.yaml --eq-time 10ns
```

**Expected behavior:**

```
Loading comparison config from: comparison.yaml
Running plugin: rg
  Equilibration: 10ns
  Conditions: no_polymer, with_polymer
  Runs: Whole Protein, Protein Backbone
Rg comparison complete
```

### All enabled analyses

Run Rg plus any other enabled plugins:

```bash
polyzymd compare run-all -f comparison.yaml --eq-time 10ns
```
````

````{tab-item} Python
Use the analysis plugin workflow for reproducible results.

```python
from pathlib import Path
import subprocess

# Rg analysis is configured in comparison.yaml and executed
# through the compare workflow
comparison_file = Path("comparison.yaml")
subprocess.run(
    [
        "polyzymd",
        "compare",
        "run",
        "rg",
        "-f",
        str(comparison_file),
        "--eq-time",
        "10ns",
    ],
    check=True,
)
```

For programmatic access to results:

```python
from polyzymd.analyses.discovery import get_analysis
from polyzymd.analyses.orchestrator import run_comparison
from polyzymd.compare.config import ComparisonConfig

config = ComparisonConfig.from_yaml("comparison.yaml")
analysis = get_analysis("rg")()
pipeline_result = run_comparison(analysis, config, equilibration="10ns")
result = pipeline_result["comparison"]

# Per-run rankings (lowest Rg = most compact)
for run_label, ranking in result.ranking_by_run.items():
    print(f"{run_label}: {ranking}")
```
````

`````

## Multi-Run Configuration

The Rg plugin uses a **runs** list, where each run defines a named Rg
calculation with its own selection. Unlike RMSD, Rg runs only need a `label`
and `selection` — there is no alignment or reference structure to configure
because Rg is inherently rotation and translation invariant.

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
        selection: "chainID C"
```

Each run produces an independent Rg timeseries. During comparison, each run
is ranked and tested separately — averaging Rg from different selections is
not meaningful.

```{important}
**Runs ≠ Replicates.** A "run" is a named Rg selection (e.g., "Whole
Protein" vs "Polymer"). A "replicate" is an independent simulation repeat
(run_1, run_2, run_3). All configured runs are computed for every replicate.
```

## Why No Alignment or Reference?

Unlike RMSD, Rg does not require trajectory alignment or a reference structure.
This is because Rg depends only on the mass-weighted distances of atoms from
their common center of mass — a quantity that is invariant under translations
and rotations.

This means:
- **No `alignment_selection`** field — alignment is unnecessary
- **No `reference_mode`** field — there is no reference structure
- **No `reference_file`** or `reference_frame`** — not applicable
- **No `centroid_selection`** — not applicable

Each run needs only `label` and `selection`. This makes Rg configuration
simpler than RMSD while providing complementary information about structural
compactness.

```{tip}
Because Rg is alignment-free, it avoids artifacts that can arise from
imperfect alignment in RMSD calculations. If you see unexpected RMSD behavior,
Rg can serve as an independent check — both should agree on whether a protein
is undergoing conformational changes.
```

## Configuration Reference

All fields for `RgRunSettings`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | `str` | *required* | Human-readable run label (must be unique) |
| `selection` | `str` | *required* | MDAnalysis selection for Rg calculation |

Top-level `RgSettings` contains a single field:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `runs` | `list[RgRunSettings]` | *required* | One or more named Rg runs (at least one required) |

```{note}
Run labels must be unique within a single `comparison.yaml`. Duplicate labels
raise a validation error.
```

```{warning}
Unlike RMSD (which defaults `selection` to `"protein and name CA"`), Rg has
**no default selection** — you must always specify `selection` explicitly for
each run. Omitting it raises a validation error.
```

## Output Files

Results are saved in your project's analysis directory:

```text
<projects_directory>/
└── analysis/
    └── rg/
        ├── run_1/
        │   ├── rg_eq10.00ns.json
        │   ├── rg_Whole Protein_timeseries.npz
        │   └── rg_Protein Backbone_timeseries.npz
        ├── run_2/
        │   ├── rg_eq10.00ns.json
        │   ├── rg_Whole Protein_timeseries.npz
        │   └── rg_Protein Backbone_timeseries.npz
        ├── run_3/
        │   └── ...
        └── aggregated/
            └── rg_reps1-3_eq10.00ns.json
```

Each replicate directory contains:
- **JSON result** — summary statistics for all configured runs
- **NPZ sidecar(s)** — raw per-frame Rg timeseries (one per run)

### JSON Result Structure

Per-replicate result (`RgResult`):

```python
{
    "config_hash": "abc123...",
    "replicate": 1,
    "equilibration_time": 10.0,
    "equilibration_unit": "ns",
    "selection_string": "protein; protein and name CA",
    "n_frames_total": 10000,
    "n_frames_used": 9000,
    "trajectory_files": ["..."],
    "run_results": [
        {
            "run_label": "Whole Protein",
            "selection": "protein",
            "mean_rg": 18.234,
            "std_rg": 0.412,
            "median_rg": 18.191,
            "min_rg": 17.087,
            "max_rg": 19.604,
            "final_rg": 18.356,
            "sem_rg": 0.098,
            "correlation_time": 3821.3,
            "correlation_time_unit": "ps",
            "n_independent_frames": 19,
            "statistical_inefficiency": 473.7,
            "n_frames_total": 10000,
            "n_frames_used": 9000,
            "npz_path": ".../rg_Whole Protein_timeseries.npz",
            "time_unit": "ns",
            "timestep_ps": 10.0
        }
    ]
}
```

Aggregated result (`RgAggregatedResult`):

```python
{
    "replicates": [1, 2, 3],
    "n_replicates": 3,
    "run_results": [
        {
            "run_label": "Whole Protein",
            "selection": "protein",
            "overall_mean": 18.256,
            "overall_sem": 0.044,
            "overall_median": 18.223,
            "per_replicate_means": [18.234, 18.291, 18.244],
            "per_replicate_stds": [0.412, 0.398, 0.424],
            "per_replicate_medians": [18.191, 18.262, 18.216]
        }
    ]
}
```

## Visualization

`````{tab-set}

````{tab-item} CLI
Generate plugin plots after running comparisons:

```bash
polyzymd compare run rg -f comparison.yaml --eq-time 10ns
polyzymd compare plot-all -f comparison.yaml
```

Plots are saved to `<projects_directory>/plots/rg/`.
````

````{tab-item} Python
Plotting is handled by the Rg plugin `plot()` method, which is invoked
through the compare plotting workflow:

```python
import subprocess

subprocess.run(
    ["polyzymd", "compare", "plot-all", "-f", "comparison.yaml"],
    check=True,
)
```
````

`````

### Available Plot Types

The Rg plugin generates figures through `polyzymd compare plot-all`:

| Plot output | Description |
|-------------|-------------|
| `rg_timeseries_<run>.png` | Mean Rg vs time with SEM shading, one per run |
| `rg_comparison_<run>.png` | Grouped bar chart of mean Rg across conditions, one per run |

**Timeseries plot features:**
- Mean Rg curve per condition with SEM shading
- Legend placed outside the plot area (`bbox_to_anchor=(1.02, 0.5)`)
- Optional per-replicate traces via `show_per_replicate: true`

### Plot Settings

Rg plot behavior can be customized in `comparison.yaml`:

```yaml
plot_settings:
  rg:
    show_per_replicate: false    # Overlay individual replicate traces
    figsize: [10, 6]             # Default figure size (bar charts)
    timeseries_figsize: [12, 5]  # Timeseries figure size (wider)
```

## Common Options

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --file` | `comparison.yaml` | Path to comparison configuration |
| `--eq-time` | `0ns` | Equilibration time to skip |
| `--recompute` | off | Ignore cached results and recompute |
| `--format` | `table` | Output format (`table` or `json`) |
| `-o, --output` | (none) | Save formatted output to file |
| `-q, --quiet` | off | Suppress INFO messages |
| `--debug` | off | Enable DEBUG logging |

### Replicate Specification

Replicates are configured per condition in `comparison.yaml`:

```yaml
conditions:
  - label: "no_polymer"
    config: "configs/no_polymer.yaml"
    replicates: [1, 3, 5]
```

### Equilibration Time

Always skip the equilibration period:

```bash
# Skip first 10 ns
polyzymd compare run rg -f comparison.yaml --eq-time 10ns

# Skip first 100 ns (longer equilibration)
polyzymd compare run rg -f comparison.yaml --eq-time 100ns
```

```{tip}
A good rule of thumb: skip at least 5–10% of your total simulation time,
or until Rg has plateaued. Use the Rg timeseries (or RMSD timeseries) as
the diagnostic for choosing an appropriate equilibration cutoff.
```

## Comparing Rg Across Conditions

To statistically compare Rg across multiple simulation conditions (e.g.,
different polymer compositions), use the `compare run rg` command:

```bash
# Add rg section to comparison.yaml, then:
polyzymd compare run rg -f comparison.yaml --eq-time 10ns
```

This provides **per-run**:
- **Ranking**: Conditions sorted by mean Rg (lowest = most compact)
- **Pairwise t-tests**: With p-values, Cohen's d, percent change
- **Direction labels**: `compaction` (lower Rg), `expansion` (higher Rg), or `unchanged`
- **ANOVA**: Omnibus test when 3+ conditions are present

**Example output:**

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

See [Comparing Conditions Guide](analysis_compare_conditions.md) for the full
multi-plugin comparison documentation.

## Troubleshooting

### "Selection matched no atoms"

**Cause:** MDAnalysis selection doesn't match any atoms in your topology.

**Fix:**
- Check residue numbering in your PDB vs. MDAnalysis (0-indexed vs 1-indexed)
- Verify atom names match your topology
- Use `polyzymd --debug compare run rg -f comparison.yaml ...` for detailed
  diagnostics

### "At least one Rg run must be defined"

**Cause:** The `runs` list in `plugins.rg` is empty or missing.

**Fix:** Add at least one run entry with `label` and `selection` fields:

```yaml
plugins:
  rg:
    runs:
      - label: "Whole Protein"
        selection: "protein"
```

### "Rg run labels must be unique"

**Cause:** Two or more runs share the same `label` value.

**Fix:** Ensure each run has a distinct label:

```yaml
plugins:
  rg:
    runs:
      - label: "Whole Protein"       # OK — unique
        selection: "protein"
      - label: "Protein Backbone"    # OK — unique
        selection: "protein and name CA"
```

### "Equilibration removed all frames"

**Cause:** `--eq-time` is longer than the trajectory.

**Fix:** Reduce the equilibration time or check that simulations completed
successfully:

```bash
# If your simulation is 50 ns, don't skip 100 ns:
polyzymd compare run rg -f comparison.yaml --eq-time 10ns  # correct
```

### Very large Rg fluctuations (> 5 Å std)

**Cause:** Usually indicates partial unfolding or a very flexible selection.

**Fix:**
- Check that your selection is appropriate (whole protein vs backbone)
- Verify trajectory files are complete
- Inspect the timeseries for jumps indicating conformational transitions
- Consider using a more specific selection (e.g., core residues only)

### "Low statistical reliability" warning

**Cause:** Long correlation time relative to trajectory length.

**This is informational, not an error.** Results are still valid but
uncertainties may be underestimated.

**Mitigation:**
- Use multiple replicates (aggregated SEM is more reliable)
- Run longer simulations
- Results are still useful for qualitative comparisons

### Missing replicate data

**Message:** `Skipping replicate N: trajectory data not found`

**Cause:** The requested replicate hasn't completed or path is incorrect.

**Fix:** This is informational — analysis continues with available replicates.
Check simulation status if unexpected.

## Rg vs RMSD Comparison

| Feature | Rg | RMSD |
|---------|-----|------|
| **Measures** | Structural compactness (mass-weighted size) | Deviation from a reference structure |
| **Output** | One value per frame (timeseries) | One value per frame (timeseries) |
| **Reference** | None needed (translation/rotation invariant) | Required (centroid/average/external) |
| **Alignment** | Not needed | Required before calculation |
| **Configuration** | `label` + `selection` only | `label` + `selection` + alignment + reference |
| **Direction labels** | `compaction` / `expansion` / `unchanged` | `stabilizing` / `destabilizing` / `unchanged` |
| **Detects** | Unfolding, compaction, swelling | Conformational drift, structural deviation |
| **Typical range** | 12–25 Å (whole protein) | 0.5–5 Å (Cα atoms) |
| **Best for** | Folding/unfolding monitoring, compactness changes | Equilibration assessment, stability comparison |

```{tip}
Use both Rg and RMSD together for a comprehensive picture: RMSD tells you how
much the structure deviates from a reference, while Rg tells you whether it is
expanding or compacting. A protein can have low RMSD (staying near its
reference) but changing Rg (gradual expansion), or vice versa.
```

## Next Steps

- **Understand Rg interpretation**: [Rg Best Practices](analysis_rg_best_practices.md)
- **Compare conditions**: [Comparing Conditions Guide](analysis_compare_conditions.md)
- **RMSD analysis**: [RMSD Quick Start](analysis_rmsd_quickstart.md)
- **RMSF analysis**: [RMSF Quick Start](analysis_rmsf_quickstart.md)
- **Understand statistics**: [Statistical Best Practices](analysis_statistics_best_practices.md)
- **Distance analysis**: [Distance Quick Start](analysis_distances_quickstart.md)
- **Contact analysis**: [Contacts Quick Start](analysis_contacts_quickstart.md)
