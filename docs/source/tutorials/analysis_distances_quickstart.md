# Distance Analysis: Quick Start

Compute inter-atomic distances with proper statistical handling in under 5 minutes.

```{note}
**Want to understand the statistics?** This guide focuses on getting results quickly.
For proper uncertainty quantification (autocorrelation correction, SEM vs. SD),
see the [Statistical Best Practices Guide](analysis_rmsf_best_practices.md).
```

## TL;DR

```bash
# Single distance pair, single replicate
polyzymd analyze distances -c config.yaml -r 1 --eq-time 10ns \
    --pair "resid 77 and name OG : resid 133 and name NE2"

# Multiple pairs with contact threshold
polyzymd analyze distances -c config.yaml -r 1-3 --eq-time 10ns \
    --pair "resid 77 and name OG : resid 133 and name NE2" \
    --pair "resid 133 and name NE2 : resid 156 and name OD1" \
    --threshold 3.5

# With plots
polyzymd analyze distances -c config.yaml -r 1 --eq-time 10ns \
    --pair "resid 77 and name OG : resid 133 and name NE2" --plot
```

## Prerequisites

Before running distance analysis, you need:

1. **Completed production simulation(s)** - at least one replicate
2. **Config file** - the same `config.yaml` used for the simulation
3. **Trajectory files** - in the scratch directory specified in config

Verify your setup:

```bash
# Check that trajectories exist
ls $(polyzymd info -c config.yaml --scratch-dir)/production_*/
```

## What Distance Analysis Provides

The distance analysis module computes:

| Feature | Description |
|---------|-------------|
| **Mean distance** | Average distance over trajectory (equilibrated portion) |
| **SEM** | Autocorrelation-corrected standard error of the mean |
| **Mode (KDE peak)** | Most probable distance from kernel density estimation |
| **Contact fraction** | \% of frames below a distance threshold |
| **Distribution** | Full histogram and KDE for visualization |

```{tip}
**When to use distances vs. contacts vs. triad:**
- **Distances**: Specific atom pairs with continuous distance values
- **Contacts**: All residue-residue contacts at an interface (binary count)
- **Triad**: Pre-defined catalytic geometry with simultaneous contact analysis
```

## Basic Usage

`````{tab-set}

````{tab-item} YAML (Recommended)
For reproducible analysis, define distance pairs in `analysis.yaml`:

```yaml
# analysis.yaml (alongside config.yaml)
replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

distances:
  enabled: true
  pairs:
    - label: "Ser77-His156"
      selection_a: "resid 77 and name OG"
      selection_b: "resid 156 and name NE2"
    - label: "His156-Asp133"
      selection_a: "resid 156 and name ND1"
      selection_b: "midpoint(resid 133 and name OD1 OD2)"
```

Then run:

```bash
# Initialize template (if starting fresh)
polyzymd analyze init

# Run all enabled analyses (uses analysis.yaml)
polyzymd analyze run

# Force recompute
polyzymd analyze run --recompute
```

**Benefits:**
- Version-controlled, reproducible
- Self-documenting experiment setup
- Easy to re-run with different parameters
````

````{tab-item} CLI
### Single Distance Pair

```bash
polyzymd analyze distances -c config.yaml -r 1 --eq-time 10ns \
    --pair "resid 77 and name OG : resid 133 and name NE2"
```

**Expected output:**

```
Loading configuration from: config.yaml
Distance Analysis: MySimulation
  Replicates: 1
  Equilibration: 10ns
  Distance pairs: 1
    1. resid 77 and name OG <-> resid 133 and name NE2

Distance Analysis Complete
  resid77_OG-resid133_NE2:
    Mean: 3.42 ± 0.15 Å
    Min:  2.61 Å
    Max:  5.87 Å
```

### Multiple Pairs

Specify `--pair` multiple times:

```bash
polyzymd analyze distances -c config.yaml -r 1 --eq-time 10ns \
    --pair "resid 77 and name OG : resid 133 and name NE2" \
    --pair "resid 133 and name NE2 : resid 156 and name OD1"
```

### Multiple Replicates

Aggregates results with SEM across replicates:

```bash
polyzymd analyze distances -c config.yaml -r 1-3 --eq-time 10ns \
    --pair "resid 77 and name OG : resid 133 and name NE2"
```

**Output:**

```
Distance Analysis Complete (Aggregated)
  Replicates: 1-3
  resid77_OG-resid133_NE2:
    Mean: 3.38 ± 0.08 Å (SEM across 3 replicates)
```
````

````{tab-item} Python
Use `DistanceCalculator` for programmatic analysis:

```python
from polyzymd.config.schema import SimulationConfig
from polyzymd.analysis import DistanceCalculator

# Load configuration
config = SimulationConfig.from_yaml("config.yaml")

# Define distance pairs
pairs = [
    ("resid 77 and name OG", "resid 156 and name NE2"),
    ("resid 156 and name ND1", "midpoint(resid 133 and name OD1 OD2)"),
]

# Create calculator
calc = DistanceCalculator(
    config=config,
    pairs=pairs,
    equilibration="10ns",
    threshold=3.5,  # Optional: contact analysis
)

# Single replicate
result = calc.compute(replicate=1)

for pr in result.pair_results:
    print(f"{pr.pair_label}: {pr.mean_distance:.2f} ± {pr.sem_distance:.2f} Å")
    if pr.fraction_below_threshold is not None:
        print(f"  Contact fraction: {pr.fraction_below_threshold:.1%}")

# Multiple replicates (aggregated with SEM)
agg_result = calc.compute_aggregated(replicates=[1, 2, 3])

for pr in agg_result.pair_results:
    print(f"{pr.pair_label}: {pr.overall_mean:.2f} ± {pr.overall_sem:.2f} Å")
```
````

`````

## Contact Threshold Analysis

Add `--threshold` to compute the fraction of frames where the distance is below
a cutoff (useful for hydrogen bond analysis, active site proximity, etc.):

`````{tab-set}

````{tab-item} YAML (Recommended)
```yaml
# analysis.yaml
distances:
  enabled: true
  threshold: 3.5  # Angstroms (H-bond cutoff)
  pairs:
    - label: "Ser77-His156"
      selection_a: "resid 77 and name OG"
      selection_b: "resid 156 and name NE2"
```

```bash
polyzymd analyze run
```
````

````{tab-item} CLI
```bash
polyzymd analyze distances -c config.yaml -r 1-3 --eq-time 10ns \
    --pair "resid 77 and name OG : resid 133 and name NE2" \
    --threshold 3.5
```

**Output:**

```
Distance Analysis Complete
  resid77_OG-resid133_NE2:
    Mean: 3.42 ± 0.15 Å
    Min:  2.61 Å
    Max:  5.87 Å
    Contact fraction (<3.5Å): 62.4%
```
````

````{tab-item} Python
```python
calc = DistanceCalculator(
    config=config,
    pairs=pairs,
    equilibration="10ns",
    threshold=3.5,  # Contact cutoff in Angstroms
)

result = calc.compute(replicate=1)
for pr in result.pair_results:
    if pr.fraction_below_threshold is not None:
        print(f"{pr.pair_label}: {pr.fraction_below_threshold:.1%} below {pr.threshold} Å")
```
````

`````

## Special Selection Syntax

PolyzyMD extends MDAnalysis selections with two special position modes:

| Syntax | Description | Use Case |
|--------|-------------|----------|
| `midpoint(selection)` | Geometric midpoint of selected atoms | Carboxylate groups (Asp, Glu) |
| `com(selection)` | Center of mass of selected atoms | Entire residues, aromatic rings |

### Examples

```yaml
# Midpoint of Asp carboxylate oxygens
selection_a: "midpoint(resid 133 and name OD1 OD2)"

# Center of mass of entire ligand
selection_b: "com(resname LIG)"

# Standard single atom
selection_a: "resid 77 and name OG"
```

```{tip}
Use `midpoint()` for carboxylate groups (Asp, Glu) where either oxygen can
accept a hydrogen bond. This gives a single representative point instead of
choosing arbitrarily between OD1/OD2 or OE1/OE2.
```

### CLI Syntax

On the command line, use quotes to protect the special syntax:

```bash
polyzymd analyze distances -c config.yaml -r 1 --eq-time 10ns \
    --pair "resid 156 and name ND1 : midpoint(resid 133 and name OD1 OD2)"
```

## Output Files

Results are saved in your project's analysis directory:

```
<projects_directory>/
└── analysis/
    └── distances/
        ├── run_1/
        │   └── distances_resid77_OG-resid133_NE2_eq10ns.json
        ├── run_2/
        │   └── distances_resid77_OG-resid133_NE2_eq10ns.json
        ├── run_3/
        │   └── distances_resid77_OG-resid133_NE2_eq10ns.json
        └── aggregated/
            └── distances_reps1-3_eq10ns.json
```

### JSON Result Structure

```python
{
    "pair_results": [
        {
            "pair_label": "resid77_OG-resid133_NE2",
            "selection1": "resid 77 and name OG",
            "selection2": "resid 133 and name NE2",
            "mean_distance": 3.42,
            "std_distance": 0.87,
            "sem_distance": 0.15,  # Autocorrelation-corrected
            "median_distance": 3.31,
            "min_distance": 2.61,
            "max_distance": 5.87,
            "kde_peak": 3.18,  # Mode from KDE
            "threshold": 3.5,
            "fraction_below_threshold": 0.624,
            "correlation_time": 245.3,  # ps
            "n_independent_frames": 34,
            "histogram_edges": [...],
            "histogram_counts": [...],
            "kde_x": [...],
            "kde_y": [...]
        }
    ],
    "n_frames_total": 10000,
    "n_frames_used": 9000,
    "equilibration_time": 10.0,
    "equilibration_unit": "ns",
    # ... additional metadata
}
```

## Visualization

`````{tab-set}

````{tab-item} CLI
Generate plots automatically with `--plot`:

```bash
polyzymd analyze distances -c config.yaml -r 1 --eq-time 10ns \
    --pair "resid 77 and name OG : resid 133 and name NE2" \
    --plot
```

Plots are saved to `<projects_directory>/plots/distances/`.
````

````{tab-item} Python
Use the plotting functions for custom figures:

```python
from polyzymd.analysis.distances import (
    plot_distance_histogram,
    plot_distance_timeseries,
    plot_distance_comparison,
    plot_contact_fraction_bar,
)

# Single distribution
result = calc.compute(replicate=1)
fig, ax = plot_distance_histogram(result.pair_results[0])
fig.savefig("distance_histogram.png")

# Time series (requires store_distributions=True)
fig, ax = plot_distance_timeseries(result.pair_results[0])
fig.savefig("distance_timeseries.png")

# Compare multiple conditions
results_no_poly = calc_no_poly.compute(replicate=1)
results_with_poly = calc_with_poly.compute(replicate=1)

fig, ax = plot_distance_comparison(
    [results_no_poly.pair_results[0], results_with_poly.pair_results[0]],
    labels=["No Polymer", "With Polymer"],
)
fig.savefig("distance_comparison.png")
```
````

`````

### Available Plot Types

The CLI `--plot` flag generates histograms automatically. For other plot types,
use the Python API:

| Function | Description | CLI | Python |
|----------|-------------|:---:|:------:|
| `plot_distance_histogram` | Distribution with optional threshold line | ✓ | ✓ |
| `plot_distance_timeseries` | Distance over frame number | | ✓ |
| `plot_distance_comparison` | Overlay multiple conditions | | ✓ |
| `plot_contact_fraction_bar` | Bar chart of contact fractions | ✓* | ✓ |

*\*Only generated when `--threshold` is specified with multiple replicates.*

```{note}
**Want more CLI plot options?** See [Issue #27](https://github.com/joelaforet/polyzymd/issues/27) 
and [Issue #28](https://github.com/joelaforet/polyzymd/issues/28) for planned enhancements
to automatic plot generation.
```

## Common Options

| Option | Default | Description |
|--------|---------|-------------|
| `-c, --config` | (required) | Path to config.yaml |
| `-r, --replicates` | `1` | Which replicates to analyze |
| `--eq-time` | `0ns` | Equilibration time to skip |
| `--pair` | (required) | Distance pair as `selection1 : selection2` |
| `--threshold` | (none) | Contact cutoff in Angstroms |
| `--plot` | off | Generate matplotlib figures |
| `--recompute` | off | Ignore cached results |
| `-o, --output-dir` | (auto) | Custom output location |

### Replicate Specification

| Format | Meaning |
|--------|---------|
| `-r 1` | Single replicate |
| `-r 1-5` | Replicates 1 through 5 |
| `-r 1,3,5` | Specific replicates |

## Troubleshooting

### "Selection matched no atoms"

**Cause:** MDAnalysis selection doesn't match any atoms in your topology.

**Fix:**
- Check residue numbering in your PDB vs. MDAnalysis (0-indexed vs 1-indexed)
- Verify atom names match your topology: `protein and resid 77` to see available atoms
- Use `polyzymd -v analyze distances ...` for detailed selection diagnostics

### Very wide distance distribution

**Cause:** The selected atoms may be flexible or the selection is too broad.

**Fix:**
- Ensure selections resolve to single atoms (or use `midpoint()`/`com()`)
- Check that `selection1` and `selection2` are correctly specified
- Visualize the selections in a molecular viewer

### "Low statistical reliability" warning

**Cause:** Long correlation time relative to trajectory length.

**This is informational, not an error.** Results are still valid but uncertainties
may be underestimated.

**Mitigation:**
- Use multiple replicates (aggregated SEM is more reliable)
- Run longer simulations
- Results are still useful for qualitative comparisons

### Missing replicate data

**Message:** `Skipping replicate N: trajectory data not found`

**Cause:** The requested replicate hasn't completed or path is incorrect.

**Fix:** This is informational—analysis continues with available replicates.
Check simulation status if unexpected.

## Comparison with Catalytic Triad Analysis

Distance analysis and [Catalytic Triad Analysis](analysis_triad_quickstart.md)
both measure atom-pair distances, but serve different purposes:

| Feature | Distances | Catalytic Triad |
|---------|-----------|-----------------|
| **Focus** | Any atom pairs | Pre-defined catalytic geometry |
| **Configuration** | `analysis.yaml` or CLI | `comparison.yaml` with conditions |
| **Multi-condition** | Run separately per condition | Built-in condition comparison |
| **Simultaneous contacts** | Not computed | Key metric (all pairs < threshold) |
| **Use case** | Ad-hoc distance measurements | Structured enzyme comparisons |

```{tip}
Use **distances** for exploratory analysis of specific interactions.
Use **catalytic triad** when comparing enzyme integrity across conditions.
```

## Next Steps

- **Catalytic triad analysis**: [Triad Quick Start](analysis_triad_quickstart.md)
- **Understand statistics**: [Statistical Best Practices](analysis_rmsf_best_practices.md)
- **Compare conditions**: [Comparing Conditions Guide](analysis_compare_conditions.md)
- **Contact analysis**: [Contacts Quick Start](analysis_contacts_quickstart.md)
