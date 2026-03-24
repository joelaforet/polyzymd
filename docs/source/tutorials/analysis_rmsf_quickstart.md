# RMSF Analysis: Quick Start

Get RMSF (Root Mean Square Fluctuation) results in under 5 minutes.

```{note}
**Want to understand the statistics?** This guide focuses on getting results quickly.
For proper uncertainty quantification and statistical best practices, see the
[Best Practices Guide](analysis_rmsf_best_practices.md).
```

## TL;DR

```bash
# Single replicate (basic)
polyzymd analyze rmsf -c config.yaml -r 1 --eq-time 10ns

# Multiple replicates (recommended)
polyzymd analyze rmsf -c config.yaml -r 1-3 --eq-time 10ns

# With plot
polyzymd analyze rmsf -c config.yaml -r 1-3 --eq-time 10ns --plot

# Force recompute (ignore cache)
polyzymd analyze rmsf -c config.yaml -r 1-3 --eq-time 10ns --recompute
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
polyzymd analyze rmsf -c config.yaml -r 1 --eq-time 10ns
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
polyzymd analyze rmsf -c config.yaml -r 1-3 --eq-time 10ns
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

## Using analysis.yaml

For reproducible, version-controlled analysis configuration, use `analysis.yaml` 
instead of CLI flags. Place this file alongside your `config.yaml`.

### Create analysis.yaml

```yaml
# analysis.yaml
replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

rmsf:
  enabled: true
  selection: "protein and name CA"    # MDAnalysis selection
  reference_mode: "centroid"          # centroid | average | frame | external
  reference_frame: null               # Required if reference_mode: frame
  # reference_file: "/path/to/crystal.pdb"  # Required if reference_mode: external
```

### Run All Enabled Analyses

```bash
# Initialize a template (if starting fresh)
polyzymd analyze init

# Run all analyses defined in analysis.yaml
polyzymd analyze run

# Force recompute
polyzymd analyze run --recompute
```

### Configuration Options

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `true` | Whether to run RMSF analysis |
| `selection` | str | `"protein and name CA"` | MDAnalysis selection for RMSF calculation |
| `reference_mode` | str | `"centroid"` | Alignment reference: `centroid`, `average`, `frame`, or `external` |
| `reference_frame` | int | `null` | Frame number when `reference_mode: frame` (1-indexed) |
| `reference_file` | str | `null` | Path to external PDB when `reference_mode: external` |

```{tip}
**When to use analysis.yaml vs CLI:** Use `analysis.yaml` for standard, 
reproducible workflows. Use CLI flags (`polyzymd analyze rmsf ...`) for 
one-off analyses or when exploring different parameters.
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

# RMSF-specific settings (required for polyzymd compare rmsf)
rmsf:
  selection: "protein and name CA"
  reference_mode: "centroid"
  # For external reference mode, use:
  # reference_mode: "external"
  # reference_file: "/path/to/crystal_structure.pdb"
```

```bash
# Run comparison with automatic t-tests, effect sizes, and ranking
polyzymd compare rmsf -f comparison.yaml

# Output formats
polyzymd compare rmsf -f comparison.yaml --format markdown  # For docs
polyzymd compare rmsf -f comparison.yaml --format json      # Machine-readable
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
Run analysis on each condition separately, then use `polyzymd compare show`:

```bash
# Step 1: Analyze each condition
polyzymd analyze rmsf -c no_polymer/config.yaml -r 1-3 --eq-time 10ns
polyzymd analyze rmsf -c with_polymer/config.yaml -r 1-3 --eq-time 10ns

# Step 2: Create comparison.yaml (or use polyzymd compare init)
polyzymd compare init my_comparison
# Edit my_comparison/comparison.yaml to point to your configs

# Step 3: Run comparison (uses cached RMSF results)
cd my_comparison
polyzymd compare rmsf
```

The comparison command automatically loads cached results if available, 
so you don't recompute RMSF.
````

````{tab-item} Python
Use `RMSFComparator` for programmatic comparison with full statistical output:

```python
from polyzymd.compare import ComparisonConfig, RMSFComparator

# Load comparison configuration (must have rmsf: section)
config = ComparisonConfig.from_yaml("comparison.yaml")

# Run comparison (computes RMSF if not cached)
comparator = RMSFComparator(
    config=config,
    rmsf_config=config.rmsf,
    equilibration="10ns",
)
result = comparator.compare()

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

# Save result for later
result.save("results/rmsf_comparison.json")
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
[Best Practices Guide](analysis_rmsf_best_practices.md#comparing-conditions).
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
| `-c, --config` | (required) | Path to config.yaml |
| `-r, --replicates` | (required) | Which replicates to analyze |
| `--eq-time` | `0ns` | Equilibration time to skip |
| `--selection` | `protein and name CA` | Atoms for RMSF calculation |
| `--reference-mode` | `centroid` | Reference structure type (`centroid`, `average`, `frame`, `external`) |
| `--reference-file` | (none) | Path to external PDB (required when `--reference-mode external`) |
| `--plot` | off | Generate matplotlib figure |
| `--recompute` | off | Ignore cached results |
| `-o, --output-dir` | (auto) | Custom output location |

### Equilibration Time

Always skip the equilibration period:

```bash
# Skip first 10 ns
polyzymd analyze rmsf -c config.yaml -r 1 --eq-time 10ns

# Skip first 100 ns (longer equilibration)
polyzymd analyze rmsf -c config.yaml -r 1 --eq-time 100ns
```

```{tip}
A good rule of thumb: skip at least 5-10% of your total simulation time,
or until RMSD has plateaued.
```

### Selection String

Change which atoms are analyzed:

```bash
# All backbone atoms (not just CA)
polyzymd analyze rmsf -c config.yaml -r 1 --selection "protein and backbone"

# Specific residue range
polyzymd analyze rmsf -c config.yaml -r 1 --selection "protein and name CA and resid 50-100"

# Include ligand
polyzymd analyze rmsf -c config.yaml -r 1 --selection "(protein and name CA) or resname LIG"
```

```{tip}
**Trimming flexible termini:** N- and C-terminal loops often have very high
RMSF that dominates summary statistics and obscures active-site signals. Use a
residue range selection to exclude them:
`--selection "protein and name CA and resid 5:175"`
This is especially useful with `external` reference mode for catalytic
competence analysis.
```

### Reference Mode

Choose how the trajectory is aligned before RMSF calculation:

```bash
# Default: align to most populated conformation
polyzymd analyze rmsf -c config.yaml -r 1 --reference-mode centroid

# Align to average structure
polyzymd analyze rmsf -c config.yaml -r 1 --reference-mode average

# Align to specific frame (e.g., catalytically competent)
polyzymd analyze rmsf -c config.yaml -r 1 --reference-mode frame --reference-frame 500

# Align to external crystal structure (condition-independent reference)
polyzymd analyze rmsf -c config.yaml -r 1 --reference-mode external \
    --reference-file /path/to/crystal_structure.pdb
```

See [Reference Structure Selection](analysis_reference_selection.md) for details
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
- Try `--reference-mode average` to compare
- Verify trajectory files are complete

### "Low statistical reliability" warning

**Cause**: Correlation time is long relative to simulation length

**This is informational, not an error.** Your results are still valid, but
uncertainties may be underestimated. See [Best Practices Guide](analysis_rmsf_best_practices.md#understanding-the-warnings) for details.

**Quick fixes**:
- Use multiple replicates (recommended)
- Results are still useful for qualitative comparisons

### No output / silent failure

**Fix**: Add debug flag to see detailed logging:

```bash
polyzymd --debug analyze rmsf -c config.yaml -r 1 --eq-time 10ns
```

### Missing Replicate Warning

**Message**: `Skipping replicate N: trajectory data not found`

**Cause**: The requested replicate hasn't been simulated yet or the path is incorrect

**Fix**: This is informational - analysis continues with available replicates.
If this is unexpected, check that the simulation completed and paths are correct.
See [Handling Incomplete Data](analysis_rmsf_best_practices.md#handling-incomplete-data)
for details.

## Next Steps

- **Understand the statistics**: [Best Practices Guide](analysis_rmsf_best_practices.md)
- **Choose reference structures**: [Reference Selection Guide](analysis_reference_selection.md)
- **Analyze distances**: `polyzymd analyze distance --help`
