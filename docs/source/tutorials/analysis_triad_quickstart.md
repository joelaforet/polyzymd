# Catalytic Triad Analysis: Quick Start

Analyze catalytic triad geometry and integrity in under 5 minutes.

```{note}
**Want to understand the statistics?** This guide focuses on getting results quickly.
For proper interpretation and statistical best practices, see the
[Best Practices Guide](analysis_triad_best_practices.md).
```

## TL;DR

```bash
# Single replicate
polyzymd analyze triad -c comparison.yaml --condition "No Polymer" -r 1 --eq-time 100ns

# Multiple replicates (recommended)
polyzymd analyze triad -c comparison.yaml --condition "No Polymer" --eq-time 100ns

# All conditions in comparison.yaml
polyzymd analyze triad -c comparison.yaml --eq-time 100ns

# Force recompute (ignore cache)
polyzymd analyze triad -c comparison.yaml --condition "No Polymer" --recompute
```

## Prerequisites

Before running catalytic triad analysis, you need:

1. **Completed production simulation(s)** - at least one replicate
2. **A comparison.yaml file** - with a `catalytic_triad` section defining your triad
3. **Trajectory files** - in the scratch directory specified in your simulation configs

## Setting Up comparison.yaml

The triad analysis uses a `comparison.yaml` file to define both your simulation
conditions and the catalytic triad geometry to analyze.

### Basic Structure

```yaml
name: "my_enzyme_study"
description: "Effect of polymer on enzyme stability"

conditions:
  - label: "No Polymer"
    config: "../projects/noPoly/config.yaml"
    replicates: [1, 2, 3]

  - label: "With Polymer"
    config: "../projects/withPoly/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "100ns"

# Catalytic triad configuration
catalytic_triad:
  name: "MyEnzyme_catalytic_triad"
  description: "Ser-His-Asp catalytic triad"
  threshold: 3.5  # Angstroms - H-bond distance cutoff
  pairs:
    - label: "Asp-His"
      selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"
      selection_b: "protein and resid 156 and name ND1"
    - label: "His-Ser"
      selection_a: "protein and resid 156 and name NE2"
      selection_b: "protein and resid 77 and name OG"
```

### Example: LipA (Lipase A)

For the Ser-His-Asp catalytic triad of Bacillus subtilis Lipase A:

```yaml
catalytic_triad:
  name: "LipA_catalytic_triad"
  description: "Ser-His-Asp catalytic triad of Lipase A (Bacillus subtilis)"
  threshold: 3.5
  pairs:
    - label: "Asp133-His156"
      selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"
      selection_b: "protein and resid 156 and name ND1"
    - label: "His156-Ser77"
      selection_a: "protein and resid 156 and name NE2"
      selection_b: "protein and resid 77 and name OG"
```

### Selection Syntax

PolyzyMD supports three types of atom selections:

| Syntax | Description | Example |
|--------|-------------|---------|
| Standard | MDAnalysis selection | `protein and resid 77 and name OG` |
| `midpoint()` | Geometric midpoint of selected atoms | `midpoint(protein and resid 133 and name OD1 OD2)` |
| `com()` | Center of mass of selected atoms | `com(protein and resid 133 and name OD1 OD2)` |

:::{warning}
**Chain-Aware Selections Required**

Residue numbers restart at 1 for each chain in PolyzyMD systems. A selection like
`resid 77` will match residues from **all chains** (protein, polymer, and water).

For protein residues, always use `protein and resid X`:

```yaml
# INCORRECT - may include polymer/water atoms with same residue number
selection_a: "resid 77 and name OG"

# CORRECT - restricts to protein chain only
selection_a: "protein and resid 77 and name OG"
```

PolyzyMD will emit a runtime warning if your selection spans multiple chains,
but it's best to write correct selections from the start.
:::

```{tip}
Use `midpoint()` for carboxylate groups (Asp, Glu) where either oxygen can
accept an H-bond. This gives a single representative point for the acceptor.
```

## Basic Usage

`````{tab-set}

````{tab-item} YAML (Recommended)
For reproducible, version-controlled analysis, define triad analysis settings
in a `comparison.yaml` file (see [Setting Up comparison.yaml](#setting-up-comparisonyaml)):

```yaml
# comparison.yaml
name: "my_enzyme_study"

conditions:
  - label: "No Polymer"
    config: "../projects/noPoly/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "100ns"

catalytic_triad:
  name: "LipA_catalytic_triad"
  description: "Ser-His-Asp catalytic triad"
  threshold: 3.5
  pairs:
    - label: "Asp133-His156"
      selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"
      selection_b: "protein and resid 156 and name ND1"
    - label: "His156-Ser77"
      selection_a: "protein and resid 156 and name NE2"
      selection_b: "protein and resid 77 and name OG"
```

Then run with minimal CLI arguments:

```bash
# Analyze all conditions using settings from comparison.yaml
polyzymd analyze triad -c comparison.yaml

# Analyze specific condition
polyzymd analyze triad -c comparison.yaml --condition "No Polymer"

# Force recompute (ignore cache)
polyzymd analyze triad -c comparison.yaml --recompute
```

**Benefits of YAML configuration:**
- All analysis parameters are version-controlled
- Reproducible across team members and machines
- Self-documenting experiment setup
````

````{tab-item} CLI
### Single Replicate

```bash
polyzymd analyze triad -c comparison.yaml --condition "No Polymer" -r 1 --eq-time 100ns
```

**Expected output:**

```
Loading comparison config from: comparison.yaml
Triad: LipA_catalytic_triad
  Description: Ser-His-Asp catalytic triad of Lipase A (Bacillus subtilis)
  Pairs: 2
    - Asp133-His156
    - His156-Ser77
  Threshold: 3.5 A
  Equilibration: 100ns

=== No Polymer ===
  Replicates: [1]

Triad Analysis Complete
  Asp133-His156: 2.91 A (93.2% below threshold)
  His156-Ser77: 3.12 A (78.4% below threshold)

  Simultaneous contact: 74.1%
    (SEM: +/-8.2%)
```

### Multiple Replicates (Recommended)

Omit the `-r` flag to analyze all replicates defined for that condition:

```bash
polyzymd analyze triad -c comparison.yaml --condition "No Polymer" --eq-time 100ns
```

### Analyzing All Conditions

To analyze every condition in your comparison.yaml:

```bash
polyzymd analyze triad -c comparison.yaml --eq-time 100ns
```

This loops through all conditions and reports results for each.
````

````{tab-item} Python
Use the Python API for programmatic analysis and integration with custom
workflows:

```python
from polyzymd.compare.config import ComparisonConfig
from polyzymd.analysis.triad import CatalyticTriadAnalyzer
from polyzymd.config.loader import load_config

# Load comparison configuration
comp_config = ComparisonConfig.from_yaml("comparison.yaml")

# Get the first condition's simulation config
condition = comp_config.conditions[0]  # "No Polymer"
sim_config = load_config(condition.config)

# Create analyzer
analyzer = CatalyticTriadAnalyzer(
    config=sim_config,
    triad_config=comp_config.catalytic_triad,
    equilibration="100ns",
)

# Single replicate analysis
result = analyzer.compute(replicate=1)
print(f"Simultaneous contact: {result.simultaneous_contact_fraction * 100:.1f}%")

# Per-pair results
for pair in result.pair_results:
    print(f"  {pair.pair_label}: {pair.mean_distance:.2f} A "
          f"({pair.fraction_below_threshold * 100:.1f}% below threshold)")

# Multi-replicate aggregation
agg_result = analyzer.compute_aggregated(replicates=[1, 2, 3])
print(f"\nAggregated: {agg_result.overall_simultaneous_contact * 100:.1f} "
      f"+/- {agg_result.sem_simultaneous_contact * 100:.1f}%")

# Save results
result.save("triad_rep1.json")
agg_result.save("triad_aggregated.json")
```

**When to use Python:**
- Integrating triad analysis into larger pipelines
- Custom post-processing or visualization
- Programmatic iteration over many conditions
- Combining with other analysis modules
````

`````

### Understanding the Output

| Field | Meaning |
|-------|---------|
| Per-pair distance | Mean distance between the two selections |
| \% below threshold | Fraction of frames where that pair is in contact |
| Simultaneous contact | Fraction where **ALL** pairs are in contact at once |
| SEM | Standard error (autocorrelation-corrected) |

The **simultaneous contact fraction** is the key metric - it tells you what
percentage of simulation time the catalytic triad maintains proper geometry
for catalysis.

### Aggregated Output

When analyzing multiple replicates, you get aggregated statistics:

```
Triad Analysis Complete (Aggregated)
  Asp133-His156: 3.09 +/- 0.21 A (88.7 +/- 5.2% below)
  His156-Ser77: 4.03 +/- 1.07 A (65.3 +/- 18.4% below)

  Simultaneous contact: 49.9 +/- 27.3%
  Per-replicate:
    Rep 1: 74.1%
    Rep 2: 51.2%
    Rep 3: 24.5%
```

The aggregated result shows:
- Mean and SEM across replicates for each pair
- Overall simultaneous contact with SEM
- Per-replicate breakdown to assess variability

## Comparing Conditions

To compare triad geometry across conditions (e.g., with vs. without polymer)
with proper statistical analysis, use one of these approaches:

`````{tab-set}

````{tab-item} YAML (Recommended)
Create a `comparison.yaml` file with your conditions, then run the comparison
with a single command:

```yaml
# comparison.yaml
name: "polymer_triad_study"
description: "Effect of polymer on catalytic triad integrity"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]

  - label: "With Polymer"
    config: "../with_polymer/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "100ns"

catalytic_triad:
  name: "LipA_catalytic_triad"
  threshold: 3.5
  pairs:
    - label: "Asp133-His156"
      selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"
      selection_b: "protein and resid 156 and name ND1"
    - label: "His156-Ser77"
      selection_a: "protein and resid 156 and name NE2"
      selection_b: "protein and resid 77 and name OG"
```

```bash
# Run comparison with automatic t-tests, effect sizes, and ranking
polyzymd compare triad -f comparison.yaml

# Output formats
polyzymd compare triad -f comparison.yaml --format markdown  # For docs
polyzymd compare triad -f comparison.yaml --format json      # Machine-readable
```

**Output includes:**
- Simultaneous contact \% ± SEM for each condition
- \% change relative to control
- p-value (two-sample t-test)
- Cohen's d effect size
- Ranking (highest contact = best triad integrity)

See [Comparing Conditions](analysis_compare_conditions.md) for the full guide.
````

````{tab-item} CLI
Run analysis on each condition separately, then use `polyzymd compare show`:

```bash
# Step 1: Analyze each condition
polyzymd analyze triad -c comparison.yaml --condition "No Polymer"
polyzymd analyze triad -c comparison.yaml --condition "With Polymer"

# Step 2: Run comparison (uses cached triad results)
polyzymd compare triad -f comparison.yaml
```

The comparison command automatically loads cached results if available,
so you don't recompute triad analysis.
````

````{tab-item} Python
Use `TriadComparator` for programmatic comparison with full statistical output:

```python
from polyzymd.compare import ComparisonConfig, TriadComparator

# Load comparison configuration (must have catalytic_triad: section)
config = ComparisonConfig.from_yaml("comparison.yaml")

# Run comparison (computes triad analysis if not cached)
comparator = TriadComparator(config, equilibration="100ns")
result = comparator.compare()

# Access results
print(f"Ranking (best triad first): {result.ranking}")

for cond in result.conditions:
    contact_pct = cond.mean_simultaneous_contact * 100
    sem_pct = cond.sem_simultaneous_contact * 100
    print(f"{cond.label}: {contact_pct:.1f} ± {sem_pct:.1f}%")

# Statistical comparisons
for comp in result.pairwise_comparisons:
    sig = "*" if comp.significant else ""
    print(f"{comp.condition_b} vs {comp.condition_a}: "
          f"{comp.percent_change:+.1f}%, p={comp.p_value:.4f}{sig}, "
          f"d={comp.cohens_d:.2f}")

# Save result for later
result.save("results/triad_comparison.json")
```

**Example output:**
```
Ranking (best triad first): ['With Polymer', 'No Polymer']
No Polymer: 49.9 ± 27.3%
With Polymer: 87.3 ± 2.2%
With Polymer vs No Polymer: +74.9%, p=0.0892, d=1.93
```
````

`````

```{tip}
For proper statistical interpretation (understanding p-values with small N,
effect sizes, ANOVA for 3+ conditions), see the
[Best Practices Guide](analysis_triad_best_practices.md#comparing-across-conditions).
```

## Output Files

Results are saved in your project's analysis directory:

```
<projects_directory>/
└── analysis/
    └── triad/
        ├── run_1/
        │   └── triad_LipA_catalytic_triad_eq100ns.json
        ├── run_2/
        │   └── triad_LipA_catalytic_triad_eq100ns.json
        ├── run_3/
        │   └── triad_LipA_catalytic_triad_eq100ns.json
        └── aggregated/
            └── triad_LipA_catalytic_triad_reps1-3_eq100ns.json
```

### JSON Result Structure (Single Replicate)

```python
{
    "analysis_type": "catalytic_triad",
    "triad_name": "LipA_catalytic_triad",
    "triad_description": "Ser-His-Asp catalytic triad...",
    "threshold": 3.5,
    "simultaneous_contact_fraction": 0.741,
    "n_frames_simultaneous": 1482,
    "sim_contact_sem": 0.082,
    "sim_contact_correlation_time": 8542.5,
    "sim_contact_n_independent": 12,
    "pair_results": [
        {
            "pair_label": "Asp133-His156",
            "mean_distance": 2.91,
            "std_distance": 0.45,
            "fraction_below_threshold": 0.932,
            ...
        },
        ...
    ],
    "n_frames_used": 2000,
    "n_frames_total": 3000,
    ...
}
```

### JSON Result Structure (Aggregated)

```python
{
    "analysis_type": "catalytic_triad_aggregated",
    "n_replicates": 3,
    "replicates": [1, 2, 3],
    "overall_simultaneous_contact": 0.499,
    "sem_simultaneous_contact": 0.273,
    "per_replicate_simultaneous": [0.741, 0.512, 0.245],
    "pair_results": [
        {
            "pair_label": "Asp133-His156",
            "overall_mean": 3.09,
            "overall_sem": 0.21,
            "per_replicate_means": [2.91, 3.05, 3.31],
            ...
        },
        ...
    ],
    ...
}
```

## Common Options

| Option | Default | Description |
|--------|---------|-------------|
| `-c, --comparison` | (required) | Path to comparison.yaml |
| `--condition` | (all) | Specific condition label to analyze |
| `-r, --replicates` | (from yaml) | Override replicate specification |
| `--eq-time` | (from yaml) | Equilibration time to skip |
| `--recompute` | off | Ignore cached results |
| `-o, --output-dir` | (auto) | Custom output location |

### Replicate Specification Formats

| Format | Meaning |
|--------|---------|
| `-r 1` | Single replicate |
| `-r 1-5` | Replicates 1 through 5 |
| `-r 1,3,5` | Specific replicates |

### Equilibration Time

Always skip the equilibration period to ensure you're analyzing equilibrated
conformations:

```bash
# Skip first 100 ns
polyzymd analyze triad -c comparison.yaml --eq-time 100ns

# Skip first 10 ns (shorter simulations)
polyzymd analyze triad -c comparison.yaml --eq-time 10ns
```

## Troubleshooting

### "No catalytic_triad section found"

**Cause**: comparison.yaml doesn't have a `catalytic_triad` section

**Fix**: Add the `catalytic_triad` section to your comparison.yaml. See the
[Setting Up comparison.yaml](#setting-up-comparisonyaml) section above.

### Very High Distances (> 10 A)

**Cause**: Usually indicates wrong atom selections or residue numbering

**Fix**:
- Check that residue numbers match your PDB file
- Verify atom names (OD1/OD2 for Asp, ND1/NE2 for His, OG for Ser)
- Load trajectory in VMD/PyMOL to visually verify selections

### Very Low Simultaneous Contact (< 10%)

**Cause**: Could indicate:
- Triad disruption (real result)
- Wrong threshold
- Wrong selections

**Fix**:
- Check individual pair percentages - which pair is failing?
- Try increasing threshold slightly (e.g., 4.0 A)
- Visualize trajectory to confirm triad state

### Selection Errors

**"Selection 'resid X and name Y' returned 0 atoms"**

**Cause**: Atom name or residue number doesn't exist

**Fix**:
- Check residue numbering in your topology
- Verify atom names match force field conventions
- For `midpoint()`, ensure all atoms in the selection exist

### Low Statistical Reliability Warning

**Cause**: Correlation time is long relative to simulation length

**Fix**: This is informational. Use multiple replicates for robust statistics.
See [Best Practices Guide](analysis_triad_best_practices.md#autocorrelation-analysis).

### Missing Replicate Warning

**Message**: `Skipping replicate N: trajectory data not found`

**Cause**: The requested replicate hasn't been simulated yet or the path is incorrect

**Fix**: This is informational - analysis continues with available replicates.
If this is unexpected, check that the simulation completed and paths are correct.
See [Handling Incomplete Data](analysis_triad_best_practices.md#handling-incomplete-data)
for details.

## Next Steps

- **Understand the statistics**: [Best Practices Guide](analysis_triad_best_practices.md)
- **Compare conditions**: [Comparing Conditions Guide](analysis_compare_conditions.md)
- **Analyze flexibility**: `polyzymd analyze rmsf --help`
