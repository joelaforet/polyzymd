# Comparing Conditions Across Simulations

Statistically compare RMSF (and other metrics) across multiple simulation conditions
with automated t-tests, effect size calculations, and ranking.

```{note}
**New to RMSF analysis?** Start with the [Quick Start Guide](analysis_rmsf_quickstart.md)
to run individual analyses, then return here to compare conditions.
```

## TL;DR

```bash
# 1. Create comparison project
polyzymd compare init my_polymer_study

# 2. Edit the template (add your conditions)
vim my_polymer_study/comparison.yaml

# 3. Run comparison
cd my_polymer_study
polyzymd compare rmsf --eq-time 10ns

# Output formats
polyzymd compare rmsf --format table     # Console table (default)
polyzymd compare rmsf --format markdown  # For documentation
polyzymd compare rmsf --format json      # Machine-readable
```

## Overview

The `polyzymd compare` module provides a streamlined workflow for comparing
simulation conditions with proper statistical analysis. Instead of manually
loading JSON files and running t-tests, you:

1. **Define conditions** in a YAML file
2. **Run one command** to get statistics
3. **Get publication-ready output** with p-values, effect sizes, and rankings

### When to Use

| Use `polyzymd compare` | Use manual analysis |
|------------------------|---------------------|
| Comparing 2+ conditions | Single condition analysis |
| Need statistical tests | Just need RMSF values |
| Want consistent formatting | Custom visualization |
| Routine comparisons | One-off exploration |

## Quick Start Workflow

### Step 1: Initialize a Comparison Project

```bash
polyzymd compare init polymer_stability_study
```

This creates:

```
polymer_stability_study/
├── comparison.yaml    # Configuration template to edit
├── results/           # Output JSON files (auto-populated)
└── figures/           # Comparison plots (from polyzymd compare plot)
```

### Step 2: Edit comparison.yaml

Open `comparison.yaml` and define your conditions:

```yaml
name: "polymer_stability_study"
description: "Effect of polymer composition on enzyme stability"

# Control condition for relative comparisons (optional)
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../projects/noPoly_LipA_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% SBMA"
    config: "../projects/SBMA_100/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% EGMA"
    config: "../projects/EGMA_100/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"
  selection: "protein and name CA"
```

### Step 3: Validate Configuration

Before running, check that all paths are correct:

```bash
polyzymd compare validate
```

Expected output:

```
Validating: comparison.yaml

Validation PASSED

  Name: polymer_stability_study
  Conditions: 3
    - No Polymer (control): 3 replicates
    - 100% SBMA: 3 replicates
    - 100% EGMA: 3 replicates
  Equilibration: 10ns
  Selection: protein and name CA
```

### Step 4: Run Comparison

```bash
polyzymd compare rmsf
```

The module will:
1. Load (or compute) RMSF for each condition
2. Run statistical comparisons
3. Display results and save to `results/`

## The comparison.yaml Configuration

### Full Schema Reference

```yaml
# ============================================================================
# Required Fields
# ============================================================================

# Project name (used in output filenames)
name: "my_comparison"

# At least 2 conditions required
conditions:
  - label: "Condition A"           # Display name (required)
    config: "path/to/config.yaml"  # Simulation config (required)
    replicates: [1, 2, 3]          # Which replicates (required)

  - label: "Condition B"
    config: "path/to/config.yaml"
    replicates: [1, 2, 3]

# ============================================================================
# Optional Fields
# ============================================================================

# Description for documentation
description: "What you're comparing"

# Control condition for relative comparisons
# - If set: all conditions compared vs control
# - If null: all pairwise comparisons
control: "Condition A"  # or null

# Default analysis parameters (can be overridden on CLI)
defaults:
  equilibration_time: "10ns"        # Time to skip
  selection: "protein and name CA"  # Atoms for RMSF
  reference_mode: "centroid"        # Alignment reference
```

### Path Resolution

Paths in `config:` can be:

| Path Type | Example | Resolution |
|-----------|---------|------------|
| Relative | `../projects/foo/config.yaml` | Relative to comparison.yaml |
| Absolute | `/home/user/sims/config.yaml` | Used as-is |

### Replicate Specification

```yaml
# Single replicate
replicates: [1]

# Range
replicates: [1, 2, 3]

# Non-contiguous
replicates: [1, 3, 5]
```

## Understanding the Output

### Console Table Format

```
RMSF Comparison: polymer_stability_study
============================================================
Equilibration: 10ns
Selection: protein and name CA
Control: No Polymer

Condition Summary (ranked by RMSF, lowest first)
------------------------------------------------------------
Rank  Condition            Mean RMSF    SEM        N
------------------------------------------------------------
1     100% SBMA               0.551 A    0.0344  2
2     100% EGMA               0.597 A    0.0725  3
3     No Polymer              0.715 A    0.0203  3   *
4     50/50 Mix               0.728 A    0.0336  3
------------------------------------------------------------
* = control condition

Pairwise Comparisons
--------------------------------------------------------------------------------
Comparison                     % Change   p-value      Cohen's d  Effect
--------------------------------------------------------------------------------
100% SBMA vs No Polymer        -22.9%     0.0211*      4.06       large
100% EGMA vs No Polymer        -16.4%     0.1944       1.27       large
50/50 Mix vs No Polymer        +1.9%      0.7445       -0.29      small
--------------------------------------------------------------------------------
* p < 0.05
```

### Output Fields Explained

| Field | Meaning |
|-------|---------|
| **Mean RMSF** | Average RMSF across all replicates (Angstroms) |
| **SEM** | Standard error of the mean |
| **N** | Number of replicates |
| **% Change** | Relative to control (negative = more stable) |
| **p-value** | Two-sample t-test p-value |
| **Cohen's d** | Effect size (positive = control higher) |
| **Effect** | Effect size interpretation |

### Ranking

Conditions are ranked by mean RMSF:
- **Rank 1** = Lowest RMSF = Most stable
- Lower RMSF indicates less flexibility/movement

## Statistical Analysis

### T-Tests

Each condition is compared to the control (or all pairs if no control) using
an independent two-sample t-test:

- **p < 0.05**: Statistically significant difference (marked with *)
- **p > 0.05**: Not statistically significant

```{tip}
With small sample sizes (N=3), you can only detect large effects.
Non-significant p-values don't mean "no effect" -- just insufficient evidence.
See [Best Practices Guide](analysis_rmsf_best_practices.md#comparing-conditions)
for interpretation guidelines.
```

### Effect Size (Cohen's d)

Cohen's d quantifies the magnitude of the difference:

| Cohen's d | Interpretation |
|-----------|---------------|
| < 0.2 | Negligible |
| 0.2 - 0.5 | Small |
| 0.5 - 0.8 | Medium |
| > 0.8 | Large |

**Why report effect size?** A large effect (d > 0.8) suggests a meaningful
difference even if p > 0.05 due to small sample size.

### ANOVA

For 3+ conditions, one-way ANOVA tests whether **any** condition differs:

```
One-way ANOVA
----------------------------------------
F-statistic: 3.151
p-value:     0.0955
Significant: No (alpha=0.05)
```

ANOVA p < 0.05 indicates at least one condition differs from the others.
Use pairwise comparisons to identify which ones.

## Output Formats

### Table (Default)

Console-friendly ASCII table:

```bash
polyzymd compare rmsf --format table
```

### Markdown

Publication-ready tables for documentation:

```bash
polyzymd compare rmsf --format markdown -o report.md
```

### JSON

Machine-readable for further analysis:

```bash
polyzymd compare rmsf --format json -o results.json
```

## Working Example

### Comparing Polymer Stabilization

Setup for a study comparing enzyme stability with different polymer coatings:

**comparison.yaml:**

```yaml
name: "SBMA_EGMA_stabilization"
description: "Does SBMA stabilize LipA better than EGMA?"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../noPoly_LipA_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% SBMA"
    config: "../SBMA_100_DMSO/config.yaml"
    replicates: [1, 2]

  - label: "100% EGMA"
    config: "../EGMA_100_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "50/50 SBMA:EGMA"
    config: "../SBMA_EGMA_50_50_DMSO/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"
```

**Run and interpret:**

```bash
polyzymd compare rmsf
```

**Key findings from output:**

| Condition | Mean RMSF | vs Control | Significant? |
|-----------|-----------|------------|--------------|
| 100% SBMA | 0.551 A | -22.9% | Yes (p=0.021) |
| 100% EGMA | 0.597 A | -16.4% | No (p=0.194) |
| 50/50 Mix | 0.728 A | +1.9% | No (p=0.745) |

**Interpretation:**
- 100% SBMA significantly stabilizes the enzyme (22.9% lower RMSF, p < 0.05)
- 100% EGMA shows a large effect (d=1.27) but isn't statistically significant
- 50/50 Mix shows no benefit over the control

## Saved Results

Results are automatically saved to `results/`:

```
my_study/
├── comparison.yaml
├── results/
│   └── rmsf_comparison_my_study.json  # Full results
└── figures/
```

### Result JSON Structure

```python
{
    "metric": "rmsf",
    "name": "my_study",
    "control_label": "No Polymer",
    "conditions": [
        {
            "label": "No Polymer",
            "mean_rmsf": 0.715,
            "sem_rmsf": 0.020,
            "n_replicates": 3,
            "replicate_values": [0.755, 0.693, 0.696]
        },
        ...
    ],
    "pairwise_comparisons": [
        {
            "condition_a": "No Polymer",
            "condition_b": "100% SBMA",
            "percent_change": -22.9,
            "p_value": 0.0211,
            "cohens_d": 4.06,
            "significant": true
        },
        ...
    ],
    "ranking": ["100% SBMA", "100% EGMA", "No Polymer", "50/50 Mix"],
    ...
}
```

### Viewing Saved Results

```bash
# Display a saved result
polyzymd compare show results/rmsf_comparison_my_study.json

# Different format
polyzymd compare show results/rmsf_comparison_my_study.json --format markdown
```

## CLI Reference

### polyzymd compare init

Create a new comparison project:

```bash
polyzymd compare init NAME [OPTIONS]

Arguments:
  NAME                   Project name (creates directory)

Options:
  --eq-time TEXT         Default equilibration time [default: 10ns]
  -o, --output-dir PATH  Parent directory [default: current]
```

### polyzymd compare rmsf

Run RMSF comparison:

```bash
polyzymd compare rmsf [OPTIONS]

Options:
  -f, --file PATH                 Config file [default: comparison.yaml]
  --eq-time TEXT                  Override equilibration time
  --selection TEXT                Override atom selection
  --recompute                     Force recompute RMSF
  --format [table|markdown|json]  Output format [default: table]
  -o, --output PATH               Save formatted output to file
  -v, --verbose                   Show detailed logging
```

### polyzymd compare validate

Check configuration:

```bash
polyzymd compare validate [OPTIONS]

Options:
  -f, --file PATH    Config file [default: comparison.yaml]
```

### polyzymd compare show

Display saved results:

```bash
polyzymd compare show RESULT_FILE [OPTIONS]

Arguments:
  RESULT_FILE                     Path to saved JSON

Options:
  --format [table|markdown|json]  Output format [default: table]
```

### polyzymd compare plot

Generate publication-ready plots from saved results:

```bash
polyzymd compare plot RESULT_FILE [OPTIONS]

Arguments:
  RESULT_FILE                     Path to saved comparison JSON

Options:
  -o, --output-dir PATH           Output directory [default: figures/]
  --format [png|pdf|svg]          Image format [default: png]
  --dpi INTEGER                   Resolution for PNG [default: 150]
  --summary / --no-summary        Generate summary panel [default: yes]
  --show / --no-show              Display interactively [default: no]
```

## Generating Plots

The `polyzymd compare plot` command creates publication-ready figures from
comparison results.

### Quick Start

```bash
# Generate all plots
polyzymd compare plot results/rmsf_comparison_my_study.json

# High resolution for publication
polyzymd compare plot results/rmsf_comparison_my_study.json --dpi 300

# PDF format (vector graphics)
polyzymd compare plot results/rmsf_comparison_my_study.json --format pdf

# Preview interactively
polyzymd compare plot results/rmsf_comparison_my_study.json --show
```

### Generated Plots

| File | Description |
|------|-------------|
| `rmsf_comparison.png` | Bar chart of mean RMSF by condition, sorted by stability |
| `percent_change.png` | Bar chart of % change vs control |
| `effect_sizes.png` | Forest plot of Cohen's d effect sizes |
| `summary_panel.png` | Combined 3-panel figure for presentations |

### Color Coding

Plots use consistent color coding to indicate statistical significance:

| Color | Meaning |
|-------|---------|
| **Green** | Significant improvement (p < 0.05, RMSF decreased) |
| **Blue** | Large effect (Cohen's d > 0.8) but not statistically significant |
| **Gray** | Control condition, or no meaningful effect |
| **Red** | Worse than control (RMSF increased) |

### Example Output

After running:

```bash
polyzymd compare plot results/rmsf_comparison_polymer_study.json -o figures/
```

You get:

```
Generated plots:
  - figures/rmsf_comparison.png
  - figures/percent_change.png
  - figures/effect_sizes.png
  - figures/summary_panel.png
```

The **summary panel** is ideal for PowerPoint presentations, combining all
three visualizations in a single figure with labeled panels (A, B, C).

## Python API

### Programmatic Comparison

```python
from polyzymd.compare import ComparisonConfig, RMSFComparator

# Load configuration
config = ComparisonConfig.from_yaml("comparison.yaml")

# Run comparison
comparator = RMSFComparator(config, equilibration="10ns")
result = comparator.compare()

# Access results
print(f"Most stable: {result.ranking[0]}")
for cond in result.conditions:
    print(f"{cond.label}: {cond.mean_rmsf:.3f} +/- {cond.sem_rmsf:.3f} A")

# Statistical comparisons
for comp in result.pairwise_comparisons:
    if comp.significant:
        print(f"{comp.condition_b} vs {comp.condition_a}: "
              f"p={comp.p_value:.4f}, d={comp.cohens_d:.2f}")
```

### Formatting Results

```python
from polyzymd.compare import format_markdown, format_console_table

# Get markdown output
md_text = format_markdown(result)
with open("report.md", "w") as f:
    f.write(md_text)

# Console output
print(format_console_table(result))
```

### Loading Saved Results

```python
from polyzymd.compare import ComparisonResult

# Load from JSON
result = ComparisonResult.load("results/rmsf_comparison_my_study.json")

# Access data
control = result.get_condition("No Polymer")
print(f"Control RMSF: {control.mean_rmsf:.3f} A")

# Get specific comparison
comp = result.get_comparison("100% SBMA")
print(f"SBMA vs control: {comp.percent_change:+.1f}%, p={comp.p_value:.4f}")
```

### Plotting with Python

```python
from polyzymd.compare import (
    ComparisonResult,
    plot_rmsf_comparison,
    plot_percent_change,
    plot_effect_sizes,
    plot_summary_panel,
)
import matplotlib.pyplot as plt

# Load result
result = ComparisonResult.load("results/rmsf_comparison_my_study.json")

# Generate individual plots
fig1 = plot_rmsf_comparison(result, save_path="figures/rmsf.png", dpi=200)
fig2 = plot_percent_change(result, save_path="figures/pct_change.png")
fig3 = plot_effect_sizes(result, save_path="figures/effects.png")

# Generate combined summary panel
fig_summary = plot_summary_panel(
    result,
    title="Polymer Stabilization Study",
    save_path="figures/summary.png",
    dpi=300,
)

# Show interactively
plt.show()
```

### Plot Customization

Each plotting function accepts customization parameters:

```python
# Custom figure size and title
fig = plot_rmsf_comparison(
    result,
    figsize=(12, 8),
    title="Effect of Polymer Coating on Enzyme Flexibility",
    horizontal=True,      # Horizontal bars (default)
    sort_by_rmsf=True,    # Sort by RMSF (default)
    color_by_effect=True, # Color by statistical effect (default)
    show_significance=True,  # Show * for p<0.05 (default)
)

# Vertical bars for fewer conditions
fig = plot_rmsf_comparison(result, horizontal=False)
```

## Troubleshooting

### "Config not found for 'Condition'"

**Cause:** Path in `config:` is incorrect

**Fix:** Check that the path exists relative to comparison.yaml:

```bash
ls ../projects/my_condition/config.yaml
```

### "Need at least 2 conditions"

**Cause:** Only one condition defined in comparison.yaml

**Fix:** Add at least one more condition to compare

### High p-value Despite Large Effect

**Cause:** Small sample size (N=2-3)

**This is expected.** With few replicates, you can only detect very large
effects. The large Cohen's d suggests a real difference exists; run more
replicates to achieve statistical significance.

### Results Don't Match Manual Calculation

**Cause:** Different equilibration time or selection string

**Fix:** Ensure `--eq-time` and `--selection` match what you used for
individual RMSF calculations. Check `defaults:` in comparison.yaml.

## See Also

- [RMSF Quick Start](analysis_rmsf_quickstart.md) -- Run individual RMSF analysis
- [Statistical Best Practices](analysis_rmsf_best_practices.md) -- Understanding statistics
- [Reference Selection](analysis_reference_selection.md) -- Alignment options
