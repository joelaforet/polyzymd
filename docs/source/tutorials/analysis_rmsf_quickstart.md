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

## Comparing Two Conditions

To compare conditions (e.g., with vs. without polymer):

### Step 1: Run analysis on both conditions

```bash
# Condition A: No polymer
polyzymd analyze rmsf -c no_polymer/config.yaml -r 1-3 --eq-time 10ns

# Condition B: With polymer
polyzymd analyze rmsf -c with_polymer/config.yaml -r 1-3 --eq-time 10ns
```

### Step 2: Load and compare results

```python
import json
import numpy as np

# Load aggregated results
with open("no_polymer/analysis/rmsf/aggregated/rmsf_reps1-3_eq10ns.json") as f:
    no_poly = json.load(f)

with open("with_polymer/analysis/rmsf/aggregated/rmsf_reps1-3_eq10ns.json") as f:
    with_poly = json.load(f)

# Compare
print(f"No Polymer:   {no_poly['overall_mean_rmsf']:.3f} ± {no_poly['overall_sem_rmsf']:.3f} Å")
print(f"With Polymer: {with_poly['overall_mean_rmsf']:.3f} ± {with_poly['overall_sem_rmsf']:.3f} Å")

# Difference
diff = no_poly['overall_mean_rmsf'] - with_poly['overall_mean_rmsf']
print(f"Difference:   {diff:.3f} Å ({100*diff/no_poly['overall_mean_rmsf']:.1f}%)")
```

```{tip}
For proper statistical testing (t-tests, effect sizes), see the
[Comparing Conditions](analysis_rmsf_best_practices.md#comparing-conditions)
section in the Best Practices Guide.
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
| `--reference-mode` | `centroid` | Reference structure type |
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

### Reference Mode

Choose how the trajectory is aligned before RMSF calculation:

```bash
# Default: align to most populated conformation
polyzymd analyze rmsf -c config.yaml -r 1 --reference-mode centroid

# Align to average structure
polyzymd analyze rmsf -c config.yaml -r 1 --reference-mode average

# Align to specific frame (e.g., catalytically competent)
polyzymd analyze rmsf -c config.yaml -r 1 --reference-mode frame --reference-frame 500
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

**Fix**: Add verbose flag to see detailed logging:

```bash
polyzymd -v analyze rmsf -c config.yaml -r 1 --eq-time 10ns
```

## Next Steps

- **Understand the statistics**: [Best Practices Guide](analysis_rmsf_best_practices.md)
- **Choose reference structures**: [Reference Selection Guide](analysis_reference_selection.md)
- **Analyze distances**: `polyzymd analyze distance --help`
