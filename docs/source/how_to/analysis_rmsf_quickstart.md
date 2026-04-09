# RMSF Analysis: Quick Start

Get RMSF (Root Mean Square Fluctuation) results in under 5 minutes.

```{note}
**Want to understand the statistics?** This guide focuses on getting results quickly.
For proper uncertainty quantification and statistical best practices, see the
[Best Practices Guide](../explanation/analysis_rmsf_best_practices.md).
```

## TL;DR

```bash
# Single analysis (replicates defined in comparison.yaml)
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns

# All enabled analyses (including RMSF)
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Save formatted output
polyzymd compare run rmsf -f comparison.yaml --format table -o rmsf_summary.txt

# Force recompute (ignore cache)
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns --recompute
```

## Prerequisites

Before running RMSF analysis, you need:

1. **Completed production simulation(s)** - at least one replicate
2. **Config file** - the same `config.yaml` used for the simulation
3. **Trajectory files** - in the scratch directory specified in config

Verify your setup:

```bash
# Check that trajectories exist
ls $(polyzymd info -c config.yaml --scratch-dir)/production_*/
```

## Basic Usage

### Single Replicate

```bash
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns
```

**Expected output:**

```
Loading configuration from: config.yaml
RMSF Analysis: MySimulation
  Replicates: 1
  Equilibration: 10ns
  Selection: protein and name CA
  Alignment: centroid

RMSF Analysis Complete
  Mean RMSF: 0.621 ± 0.015 Å
  Min RMSF:  0.248 Å
  Max RMSF:  3.160 Å
```

### Understanding the Output

| Field | Meaning |
|-------|---------|
| Mean RMSF | Average fluctuation across all residues |
| Min RMSF | Most rigid residue (usually buried core) |
| Max RMSF | Most flexible residue (usually loops/termini) |
| ± value | Standard error of the mean |

## Multiple Replicates (Recommended)

Running multiple replicates provides more reliable statistics:

```bash
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns
```

**Replicate specification formats:**

| Format | Meaning |
|--------|---------|
| `-r 1` | Single replicate |
| `-r 1-5` | Replicates 1 through 5 |
| `-r 1,3,5` | Specific replicates |

**Aggregated output:**

```
RMSF Analysis Complete (Aggregated)
  Replicates: 1-3
  Mean RMSF: 0.715 ± 0.020 Å
  Min RMSF:  0.304 Å
  Max RMSF:  4.206 Å
```

The aggregated result combines per-residue RMSF values across replicates,
reporting the mean and standard error.

## Using comparison.yaml

For reproducible, version-controlled analysis configuration, use
`comparison.yaml` and put RMSF settings under `plugins.rmsf`.

### Create comparison.yaml

```yaml
# comparison.yaml
name: "rmsf_quickstart"
description: "Single-condition RMSF quick start"
control: "MySimulation"

conditions:
  - label: "MySimulation"
    config: "config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  rmsf:
    enabled: true
    selection: "protein and name CA"    # MDAnalysis selection
    reference_mode: "centroid"          # centroid | average | frame | external
    reference_frame: 0                   # Used if reference_mode: frame
    # reference_file: "/path/to/crystal.pdb"  # Used if reference_mode: external
```

### Run RMSF

```bash
# Run RMSF analysis
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns

# Run all enabled analyses in one command
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Run all enabled analyses and generate plots
polyzymd compare run-all -f comparison.yaml --eq-time 10ns --plot

# Force recompute
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns --recompute
```

### Configuration Options

| Field (`plugins.rmsf`) | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `true` | Whether to run RMSF analysis |
| `selection` | str | `"protein and name CA"` | MDAnalysis selection for RMSF calculation |
| `reference_mode` | str | `"centroid"` | Alignment reference: `centroid`, `average`, `frame`, or `external` |
| `reference_frame` | int | `0` | Frame number when `reference_mode: frame` |
| `reference_file` | str | `null` | Path to external PDB when `reference_mode: external` |

```{tip}
**When to use YAML vs CLI:** Use `comparison.yaml` for standard,
reproducible workflows. Use `polyzymd compare run rmsf` for quick reruns,
format changes, or `--recompute`.
```

## Comparing Two Conditions

To compare conditions (e.g., with vs. without polymer) with proper statistical 
analysis, use one of these approaches:

`````{tab-set}

````{tab-item} YAML (Recommended)
Create a `comparison.yaml` file to define your conditions, then run the 
comparison with a single command:

```yaml
# comparison.yaml
name: "polymer_stabilization"
description: "Effect of polymer on enzyme stability"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]

  - label: "With Polymer"
    config: "../with_polymer/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  rmsf:
    selection: "protein and name CA"
    reference_mode: "centroid"
    # For external reference mode, use:
    # reference_mode: "external"
    # reference_file: "/path/to/crystal_structure.pdb"
```

```bash
# Run comparison with automatic t-tests, effect sizes, and ranking
polyzymd compare run rmsf -f comparison.yaml

# Output formats
polyzymd compare run rmsf -f comparison.yaml --format table     # Human-readable
polyzymd compare run rmsf -f comparison.yaml --format json      # Machine-readable
```

**Output includes:**
- Mean RMSF ± SEM for each condition
- \% change relative to control
- p-value (two-sample t-test)
- Cohen's d effect size
- Ranking (lowest RMSF = most stable)

See [Comparing Conditions](analysis_compare_conditions.md) for the full guide.
````

````{tab-item} CLI
Run the RMSF comparison pipeline directly:

```bash
# Step 1: Create comparison.yaml (or use polyzymd compare init)
polyzymd compare init my_comparison
# Edit my_comparison/comparison.yaml to point to your configs

# Step 2: Run comparison (uses cached RMSF results when available)
cd my_comparison
polyzymd compare run rmsf --eq-time 10ns
```

The comparison command automatically loads cached results if available, 
so you don't recompute RMSF.
````

````{tab-item} Python
Use the analysis orchestrator for programmatic comparison:

```python
from polyzymd.analyses.discovery import get_analysis
from polyzymd.analyses.orchestrator import run_comparison
from polyzymd.compare.config import ComparisonConfig

# Load comparison configuration (must have plugins.rmsf section)
config = ComparisonConfig.from_yaml("comparison.yaml")

# Run comparison (computes RMSF if not cached)
analysis = get_analysis("rmsf")()
pipeline_result = run_comparison(analysis, config, equilibration="10ns")
result = pipeline_result["comparison"]

# Access results
print(f"Ranking (most stable first): {result.ranking}")

for cond in result.conditions:
    print(f"{cond.label}: {cond.mean_rmsf:.3f} ± {cond.sem_rmsf:.3f} Å")

# Statistical comparisons
for comp in result.pairwise_comparisons:
    sig = "*" if comp.significant else ""
    print(f"{comp.condition_b} vs {comp.condition_a}: "
          f"{comp.percent_change:+.1f}%, p={comp.p_value:.4f}{sig}, "
          f"d={comp.cohens_d:.2f}")

# Canonical cache path written by the orchestrator
print(pipeline_result["comparison_path"])
```

**Example output:**
```
Ranking (most stable first): ['With Polymer', 'No Polymer']
No Polymer: 0.715 ± 0.020 Å
With Polymer: 0.551 ± 0.034 Å
With Polymer vs No Polymer: -22.9%, p=0.0211*, d=4.06
```
````

`````

```{tip}
For proper statistical interpretation (understanding p-values with small N, 
effect sizes, ANOVA for 3+ conditions), see the 
[Best Practices Guide](../explanation/analysis_rmsf_best_practices.md#comparing-conditions).
```

## Output Files

Results are saved in your project's analysis directory:

```
<projects_directory>/
└── analysis/
    └── rmsf/
        ├── run_1/
        │   └── rmsf_eq10ns.json      # Single replicate result
        ├── run_2/
        │   └── rmsf_eq10ns.json
        ├── run_3/
        │   └── rmsf_eq10ns.json
        └── aggregated/
            └── rmsf_reps1-3_eq10ns.json  # Combined results
```

### JSON Result Structure

```python
{
    "overall_mean_rmsf": 0.715,
    "overall_sem_rmsf": 0.020,
    "overall_min_rmsf": 0.304,
    "overall_max_rmsf": 4.206,
    "per_replicate_mean_rmsf": [0.755, 0.693, 0.696],
    "mean_rmsf_per_residue": [0.45, 0.52, ...],  # Per-residue values
    "sem_rmsf_per_residue": [0.02, 0.03, ...],
    "residue_ids": [1, 2, 3, ...],
    "residue_names": ["MET", "ALA", "SER", ...],
    "correlation_time": 15394.5,  # ps
    "n_independent_frames": 6,
    # ... additional metadata
}
```

## Common Options

| Option | Default | Description |
|--------|---------|-------------|
| `-f, --file` | `comparison.yaml` | Comparison config file |
| `--eq-time` | from YAML (or `0ns`) | Equilibration time to skip |
| `--recompute` | off | Ignore cached results |
| `--format` | `table` | Output format (`table` or `json`) |
| `-o, --output` | (none) | Save formatted output to file |
| `-q, --quiet` | off | Suppress INFO messages |
| `--debug` | off | Enable DEBUG logging |

RMSF-specific controls (`selection`, `reference_mode`, `reference_frame`,
`reference_file`) are configured in `comparison.yaml` under `plugins.rmsf`.

### Equilibration Time

Always skip the equilibration period:

```bash
# Skip first 10 ns
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns

# Skip first 100 ns (longer equilibration)
polyzymd compare run rmsf -f comparison.yaml --eq-time 100ns
```

```{tip}
A good rule of thumb: skip at least 5-10% of your total simulation time,
or until RMSD has plateaued.
```

### Selection String

Change which atoms are analyzed by editing `plugins.rmsf.selection` in
`comparison.yaml`:

```yaml
plugins:
  rmsf:
    selection: "protein and backbone"
```

```yaml
plugins:
  rmsf:
    selection: "protein and name CA and resid 50-100"
```

```yaml
plugins:
  rmsf:
    selection: "(protein and name CA) or resname LIG"
```

```{tip}
**Trimming flexible termini:** N- and C-terminal loops often have very high
RMSF that dominates summary statistics and obscures active-site signals. Use a
residue range selection to exclude them:
`selection: "protein and name CA and resid 5:175"`
This is especially useful with `external` reference mode for catalytic
competence analysis.
```

### Reference Mode

Choose how the trajectory is aligned before RMSF calculation:

```yaml
# Default: align to most populated conformation
plugins:
  rmsf:
    reference_mode: "centroid"
```

```yaml
# Align to average structure
plugins:
  rmsf:
    reference_mode: "average"
```

```yaml
# Align to specific frame (e.g., catalytically competent)
plugins:
  rmsf:
    reference_mode: "frame"
    reference_frame: 500
```

```yaml
# Align to external crystal structure (condition-independent reference)
plugins:
  rmsf:
    reference_mode: "external"
    reference_file: "/path/to/crystal_structure.pdb"
```

See [Reference Structure Selection](../explanation/analysis_reference_selection.md) for details
on when to use each mode.

## Troubleshooting

### "Working directory not found"

**Cause**: Config points to wrong scratch directory

**Fix**: Update `output.scratch_directory` in config.yaml to match where
your trajectories are stored.

### Very high RMSF values (> 10 Å)

**Cause**: Usually indicates alignment issues or wrong selection

**Fix**: 
- Check that your selection string matches atoms in your system
- Try `reference_mode: "average"` to compare
- Verify trajectory files are complete

### "Low statistical reliability" warning

**Cause**: Correlation time is long relative to simulation length

**This is informational, not an error.** Your results are still valid, but
uncertainties may be underestimated. See [Best Practices Guide](../explanation/analysis_rmsf_best_practices.md#understanding-the-warnings) for details.

**Quick fixes**:
- Use multiple replicates (recommended)
- Results are still useful for qualitative comparisons

### No output / silent failure

**Fix**: Add debug flag to see detailed logging:

```bash
polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns --debug
```

### Missing Replicate Warning

**Message**: `Skipping replicate N: trajectory data not found`

**Cause**: The requested replicate hasn't been simulated yet or the path is incorrect

**Fix**: This is informational - analysis continues with available replicates.
If this is unexpected, check that the simulation completed and paths are correct.
See [Handling Incomplete Data](../explanation/analysis_rmsf_best_practices.md#handling-incomplete-data)
for details.

## Next Steps

- **Understand the statistics**: [Best Practices Guide](../explanation/analysis_rmsf_best_practices.md)
- **Choose reference structures**: [Reference Selection Guide](../explanation/analysis_reference_selection.md)
- **Analyze distances**: `polyzymd compare run distances --help`
