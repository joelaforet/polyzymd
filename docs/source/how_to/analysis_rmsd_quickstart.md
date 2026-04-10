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
[RMSD Best Practices Guide](../explanation/analysis_rmsd_best_practices.md).
```

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
    enabled: true
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

```
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

````{tab-item} Python
Use the analysis plugin workflow for reproducible results.

```python
from pathlib import Path
import subprocess

# RMSD analysis is configured in comparison.yaml and executed
# through the compare workflow
comparison_file = Path("comparison.yaml")
subprocess.run(
    [
        "polyzymd",
        "compare",
        "run",
        "rmsd",
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
from polyzymd.config.comparison import ComparisonConfig

config = ComparisonConfig.from_yaml("comparison.yaml")
analysis = get_analysis("rmsd")()
pipeline_result = run_comparison(analysis, config, equilibration="10ns")
result = pipeline_result["comparison"]

# Per-run rankings (lowest RMSD = most stable)
for run_label, ranking in result.ranking_by_run.items():
    print(f"{run_label}: {ranking}")
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

## Configuration Reference

All fields for `RMSDRunSettings`:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `label` | `str` | *required* | Human-readable run label (must be unique) |
| `selection` | `str` | `"protein and name CA"` | MDAnalysis selection for RMSD calculation |
| `alignment_selection` | `str` | `"protein and name CA"` | MDAnalysis selection for trajectory alignment |
| `reference_mode` | `str` | `"centroid"` | Reference mode: `centroid`, `average`, `frame`, or `external` |
| `reference_frame` | `int` | `0` | 0-indexed frame for `reference_mode: frame` |
| `reference_file` | `str \| None` | `null` | Path to external PDB for `reference_mode: external` |
| `centroid_selection` | `str \| None` | `null` | Selection for centroid finding; defaults to `alignment_selection` |
| `convergence_window_size_ns` | `float` | `15.0` | Sliding window size for convergence detection (ns) |
| `convergence_step_size_ns` | `float` | `5.0` | Step between successive windows (ns) |
| `convergence_slope_threshold` | `float` | `0.0005` | Max absolute slope to qualify as "flat" (Å/ns) |
| `convergence_sustained_for_ns` | `float` | `15.0` | Required sustained duration below threshold (ns) |

Top-level `RMSDSettings` contains a single field:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `runs` | `list[RMSDRunSettings]` | *required* | One or more named RMSD runs (at least one required) |

```{note}
Run labels must be unique within a single `comparison.yaml`. Duplicate labels
raise a validation error.
```

## Output Files

Results are saved in your project's analysis directory:

```text
<projects_directory>/
└── analysis/
    └── rmsd/
        ├── run_1/
        │   ├── rmsd_eq10.00ns.json
        │   ├── rmsd_Protein Backbone_timeseries.npz
        │   └── rmsd_Active Site_timeseries.npz
        ├── run_2/
        │   ├── rmsd_eq10.00ns.json
        │   ├── rmsd_Protein Backbone_timeseries.npz
        │   └── rmsd_Active Site_timeseries.npz
        ├── run_3/
        │   └── ...
        └── aggregated/
            └── rmsd_reps1-3_eq10.00ns.json
```

Each replicate directory contains:
- **JSON result** — summary statistics for all configured runs
- **NPZ sidecar(s)** — raw per-frame RMSD timeseries (one per run)

### JSON Result Structure

Per-replicate result (`RMSDResult`):

```python
{
    "config_hash": "abc123...",
    "replicate": 1,
    "equilibration_time": 10.0,
    "equilibration_unit": "ns",
    "selection_string": "protein and name CA; ...",
    "n_frames_total": 10000,
    "n_frames_used": 9000,
    "trajectory_files": ["..."],
    "run_results": [
        {
            "run_label": "Protein Backbone",
            "selection": "protein and name CA",
            "alignment_selection": "protein and name CA",
            "reference_mode": "centroid",
            "mean_rmsd": 1.823,
            "std_rmsd": 0.312,
            "median_rmsd": 1.791,
            "min_rmsd": 0.987,
            "max_rmsd": 3.104,
            "final_rmsd": 1.956,
            "sem_rmsd": 0.078,
            "correlation_time": 4521.3,
            "correlation_time_unit": "ps",
            "n_independent_frames": 16,
            "statistical_inefficiency": 562.7,
            "n_frames_total": 10000,
            "n_frames_used": 9000,
            "npz_path": ".../rmsd_Protein Backbone_timeseries.npz",
            "time_unit": "ns",
            "timestep_ps": 10.0,
            "converged": true,
            "convergence_assessable": true,
            "convergence_time_ns": 12.5,
            "convergence_message": "Converged at 12.500 ns"
        }
    ]
}
```

Aggregated result (`RMSDAggregatedResult`):

```python
{
    "replicates": [1, 2, 3],
    "n_replicates": 3,
    "run_results": [
        {
            "run_label": "Protein Backbone",
            "selection": "protein and name CA",
            "overall_mean": 1.856,
            "overall_sem": 0.034,
            "overall_median": 1.823,
            "per_replicate_means": [1.823, 1.891, 1.854],
            "per_replicate_stds": [0.312, 0.298, 0.324],
            "per_replicate_medians": [1.791, 1.862, 1.816],
            "n_converged_replicates": 3,
            "convergence_fraction": 1.0,
            "mean_convergence_time_ns": 13.2,
            "median_convergence_time_ns": 12.5
        }
    ]
}
```

## Visualization

`````{tab-set}

````{tab-item} CLI
Generate plugin plots after running comparisons:

```bash
polyzymd compare run rmsd -f comparison.yaml --eq-time 10ns
polyzymd compare plot-all -f comparison.yaml
```

Plots are saved to `<projects_directory>/plots/rmsd/`.
````

````{tab-item} Python
Plotting is handled by the RMSD plugin `plot()` method, which is invoked
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

The RMSD plugin generates figures through `polyzymd compare plot-all`:

| Plot output | Description |
|-------------|-------------|
| `rmsd_timeseries_<run>.png` | Mean RMSD vs time with SEM shading, one per run |
| `rmsd_comparison_<run>.png` | Grouped bar chart of mean RMSD across conditions, one per run |
| `rmsd_convergence_<condition>_<run>.png` | Dual-axis plot: RMSD timeseries with sliding-window slope and convergence marker (requires `show_convergence_plots: true`) |

**Timeseries plot features:**
- Mean RMSD curve per condition with SEM shading
- Legend placed outside the plot area (`bbox_to_anchor=(1.02, 0.5)`)
- Optional per-replicate traces via `show_per_replicate: true`

### Plot Settings

RMSD plot behavior can be customized in `comparison.yaml`:

```yaml
plot_settings:
  rmsd:
    show_per_replicate: false    # Overlay individual replicate traces
    figsize: [10, 6]             # Default figure size (bar charts)
    timeseries_figsize: [12, 5]  # Timeseries figure size (wider)
    show_convergence_plots: false  # Generate per-replicate convergence diagnostics
    convergence_figsize: [12, 5]   # Convergence panel figure size
```

## Convergence Detection

Convergence detection is always on — every RMSD run automatically applies a
sliding-window slope heuristic to determine whether the RMSD timeseries has
plateaued. This is a purely additive diagnostic: it does not affect ranking,
statistical tests, or any other comparison output. Convergence results appear
as additional fields in per-replicate and aggregated JSON files, and optional
convergence plots can be enabled via `show_convergence_plots: true`.

For a conceptual explanation of the algorithm, its parameters, and its
limitations, see {doc}`/explanation/convergence_detection`.

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
polyzymd compare run rmsd -f comparison.yaml --eq-time 10ns

# Skip first 100 ns (longer equilibration)
polyzymd compare run rmsd -f comparison.yaml --eq-time 100ns
```

```{tip}
A good rule of thumb: skip at least 5-10% of your total simulation time,
or until RMSD has plateaued. The RMSD timeseries itself is the best diagnostic
for choosing an appropriate equilibration cutoff.
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

See [Comparing Conditions Guide](analysis_compare_conditions.md) for the full
multi-plugin comparison documentation.

## Troubleshooting

### "Selection matched no atoms"

**Cause:** MDAnalysis selection doesn't match any atoms in your topology.

**Fix:**
- Check residue numbering in your PDB vs. MDAnalysis (0-indexed vs 1-indexed)
- Verify atom names match your topology
- Use `polyzymd --debug compare run rmsd -f comparison.yaml ...` for detailed
  diagnostics

### "At least one RMSD run must be defined"

**Cause:** The `runs` list in `plugins.rmsd` is empty or missing.

**Fix:** Add at least one run entry with a `label` field:

```yaml
plugins:
  rmsd:
    runs:
      - label: "Protein Backbone"
```

### "reference_file does not exist"

**Cause:** Using `reference_mode: external` but the PDB path is invalid.

**Fix:** Provide an absolute path or a path relative to the working directory:

```yaml
reference_mode: "external"
reference_file: "/absolute/path/to/crystal.pdb"
```

### "atom count mismatch between trajectory and external PDB"

**Cause:** The `selection` string matches different numbers of atoms in the
trajectory vs. the external reference PDB.

**Fix:**
- Ensure both systems use the same atom naming convention
- Check that the external PDB contains the same residues as your simulation
- Use a more specific selection if topologies differ

### Very high RMSD values (> 10 Å)

**Cause:** Usually indicates alignment issues, wrong selection, or unfolding.

**Fix:**
- Check that `alignment_selection` matches atoms in your system
- Try `reference_mode: "average"` to compare
- Verify trajectory files are complete
- Check for protein unfolding or large conformational changes

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

## RMSD vs RMSF Comparison

| Feature | RMSD | RMSF |
|---------|------|------|
| **Measures** | Global deviation from reference | Per-residue fluctuation |
| **Output** | One value per frame (timeseries) | One value per residue (profile) |
| **Reference** | Fixed structure (centroid/average/external) | Time-averaged position |
| **Detects** | Conformational drift, unfolding | Flexible loops, rigid core |
| **Multi-run** | Yes (`runs` list with different selections) | Single selection |
| **Best for** | Equilibration assessment, stability comparison | Flexibility mapping |

```{tip}
Use RMSD first to assess overall stability and choose equilibration time,
then use RMSF to identify which regions drive flexibility differences.
```

## Next Steps

- **Understand RMSD interpretation**: [RMSD Best Practices](../explanation/analysis_rmsd_best_practices.md)
- **Compare conditions**: [Comparing Conditions Guide](analysis_compare_conditions.md)
- **RMSF analysis**: [RMSF Quick Start](analysis_rmsf_quickstart.md)
- **Understand statistics**: [Statistical Best Practices](../explanation/analysis_statistics_best_practices.md)
- **Distance analysis**: [Distance Quick Start](analysis_distances_quickstart.md)
- **Contact analysis**: [Contacts Quick Start](analysis_contacts_quickstart.md)
