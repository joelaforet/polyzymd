# Comparing Conditions Across Simulations

Statistically compare RMSF, contacts, distances, and other metrics across multiple simulation conditions
with automated t-tests, effect size calculations, and ranking.

```{note}
**New to analysis?** Start with the individual quick start guides:
- [RMSF Quick Start](analysis_rmsf_quickstart.md) for flexibility analysis
- [Contacts Quick Start](analysis_contacts_quickstart.md) for polymer-protein contacts
- [Distance Analysis Quick Start](analysis_distances_quickstart.md) for inter-atomic distances

Then return here to compare conditions.
```

## TL;DR

```bash
# 1. Create comparison project
polyzymd compare init my_polymer_study

# 2. Edit the template (add your conditions)
vim my_polymer_study/comparison.yaml

# 3. Run comparisons
cd my_polymer_study
polyzymd compare run rmsf --eq-time 10ns      # Compare flexibility
polyzymd compare run triad --eq-time 10ns     # Compare triad geometry (if defined)
polyzymd compare run contacts --eq-time 10ns  # Compare polymer-protein contacts
polyzymd compare run distances --eq-time 10ns # Compare inter-atomic distances

# Output formats
polyzymd compare run rmsf --format table     # Console table (default)
polyzymd compare run rmsf --format markdown  # For documentation
polyzymd compare run rmsf --format json      # Machine-readable
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

# Analysis settings: WHAT to analyze (shared across conditions)
analysis_settings:
  rmsf:
    selection: "protein and name CA"

# Comparison settings: HOW to compare (must have entry for each analysis)
comparison_settings:
  rmsf: {}  # Empty is OK, but must be present
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

`````{tab-set}
````{tab-item} CLI (Recommended)
```bash
polyzymd compare rmsf
```

The module will:
1. Load (or compute) RMSF for each condition
2. Run statistical comparisons
3. Display results and save to `results/`
````

````{tab-item} Python
```python
from polyzymd.compare import ComparisonConfig, RMSFComparator

# Load configuration
config = ComparisonConfig.from_yaml("comparison.yaml")

# Get RMSF settings from analysis_settings
rmsf_settings = config.analysis_settings.get("rmsf")

# Run RMSF comparison
comparator = RMSFComparator(
    config=config,
    rmsf_settings=rmsf_settings,
    equilibration="10ns",
)
result = comparator.compare()

# Display results
print(f"Most stable: {result.ranking[0]}")
for cond in result.conditions:
    print(f"{cond.label}: {cond.mean_rmsf:.3f} ± {cond.sem_rmsf:.3f} Å")

# Save to JSON
result.save("results/rmsf_comparison.json")
```
````
`````

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

# Default analysis parameters (shared across all analyses)
defaults:
  equilibration_time: "10ns"        # Time to skip (used by all analyses)

# ============================================================================
# Analysis Settings (WHAT to analyze - applied to all conditions)
# ============================================================================
# Each key enables that analysis type. Presence of a section enables it.

analysis_settings:
  # RMSF Analysis
  rmsf:
    selection: "protein and name CA"  # Atoms for RMSF calculation
    reference_mode: "centroid"        # centroid, average, or frame
    # reference_frame: 500            # Required if reference_mode is "frame"

  # Catalytic triad comparison (required for `polyzymd compare triad`)
  # catalytic_triad:
  #   name: "enzyme_catalytic_triad"
  #   threshold: 3.5
  #   pairs:
  #     - label: "Asp-His"
  #       selection_a: "midpoint(resid 133 and name OD1 OD2)"
  #       selection_b: "resid 156 and name ND1"

  # Contacts comparison (required for `polyzymd compare contacts`)
  # contacts:
  #   polymer_selection: "resname SBM EGM"
  #   cutoff: 4.5

  # Distances comparison (required for `polyzymd compare distances`)
  # distances:
  #   threshold: 3.5  # Global default threshold (Angstroms, optional)
  #   pairs:
  #     - label: "Catalytic H-bond"
  #       selection_a: "resid 77 and name OG"
  #       selection_b: "resid 133 and name NE2"
  #       threshold: 3.5  # Per-pair threshold (optional, overrides global)
  #     - label: "Lid Opening"
  #       selection_a: "com(resid 141-148)"
  #       selection_b: "com(resid 281-289)"
  #       threshold: 15.0  # Different threshold for large-scale motion

# ============================================================================
# Comparison Settings (HOW to compare - statistical parameters)
# ============================================================================
# Each analysis in analysis_settings MUST have a corresponding entry here.
# Use empty {} for analyses with no comparison-specific parameters.

comparison_settings:
  rmsf: {}  # No comparison-specific parameters

  # catalytic_triad: {}

  # contacts:
  #   fdr_alpha: 0.05           # FDR for Benjamini-Hochberg correction
  #   min_effect_size: 0.5      # Cohen's d threshold
  #   top_residues: 10          # Top residues to show in console

  # distances: {}  # No comparison-specific parameters
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

`````{tab-set}
````{tab-item} YAML + CLI (Recommended)
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

analysis_settings:
  rmsf:
    selection: "protein and name CA"

comparison_settings:
  rmsf: {}
```

**Run:**

```bash
polyzymd compare rmsf
```
````

````{tab-item} Python
```python
from polyzymd.compare import ComparisonConfig, RMSFComparator

# Load comparison configuration
config = ComparisonConfig.from_yaml("comparison.yaml")

# Get RMSF settings from analysis_settings
rmsf_settings = config.analysis_settings.get("rmsf")

# Run RMSF comparison
comparator = RMSFComparator(
    config=config,
    rmsf_settings=rmsf_settings,
    equilibration="10ns",
)
result = comparator.compare()

# Print summary table
print(f"RMSF Comparison: {result.name}")
print(f"Control: {result.control_label}")
print()

# Ranked conditions
print("Conditions (ranked by RMSF, lowest first):")
for i, label in enumerate(result.ranking, 1):
    cond = result.get_condition(label)
    marker = " *" if label == result.control_label else ""
    print(f"  {i}. {label}: {cond.mean_rmsf:.3f} ± {cond.sem_rmsf:.3f} Å{marker}")

print()

# Pairwise comparisons
print("Pairwise comparisons vs control:")
for comp in result.pairwise_comparisons:
    sig = "*" if comp.significant else ""
    print(f"  {comp.condition_b} vs {comp.condition_a}: "
          f"{comp.percent_change:+.1f}%, p={comp.p_value:.3f}{sig}, d={comp.cohens_d:.2f}")

# Save results
result.save("results/rmsf_comparison.json")
```
````
`````

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
  -q, --quiet                     Suppress INFO messages
  --debug                         Enable DEBUG logging
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

`````{tab-set}
````{tab-item} CLI (Recommended)
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
````

````{tab-item} Python
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
fig1 = plot_rmsf_comparison(result, save_path="figures/rmsf.png", dpi=300)
fig2 = plot_percent_change(result, save_path="figures/pct_change.png", dpi=300)
fig3 = plot_effect_sizes(result, save_path="figures/effects.png", dpi=300)

# Generate combined summary panel
fig_summary = plot_summary_panel(
    result,
    title="Polymer Stabilization Study",
    save_path="figures/summary.png",
    dpi=300,
)

# Preview interactively
plt.show()
```
````
`````

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

# Load configuration (must have analysis_settings.rmsf section)
config = ComparisonConfig.from_yaml("comparison.yaml")

# Get RMSF settings
rmsf_settings = config.analysis_settings.get("rmsf")

# Run comparison
comparator = RMSFComparator(
    config=config,
    rmsf_settings=rmsf_settings,
    equilibration="10ns",
)
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

## Comparing Catalytic Triad Geometry

In addition to RMSF (global flexibility), you can compare **catalytic triad integrity**
across conditions. This is useful for enzymes where active site geometry is crucial
for catalytic function.

### What is Simultaneous Contact Fraction?

The catalytic triad comparison analyzes the **simultaneous contact fraction** -- the
percentage of simulation frames where ALL distance pairs in your triad are below
the contact threshold at the same time. Higher values indicate better triad integrity.

For example, a Ser-His-Asp catalytic triad:
- 95% simultaneous contact = triad is intact most of the time
- 50% simultaneous contact = triad frequently disrupted
- 10% simultaneous contact = triad rarely intact

### Adding Catalytic Triad to comparison.yaml

Add a `catalytic_triad` section to your `analysis_settings`:

```yaml
name: "polymer_stability_study"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../noPoly_LipA_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% SBMA"
    config: "../SBMA_100_DMSO/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

# Define your enzyme's catalytic triad in analysis_settings
analysis_settings:
  catalytic_triad:
    name: "LipA_Ser-His-Asp"
    description: "Lipase A catalytic triad"
    threshold: 3.5  # Angstroms (H-bond cutoff)
    pairs:
      - label: "Asp133-His156"
        selection_a: "midpoint(resid 133 and name OD1 OD2)"
        selection_b: "resid 156 and name ND1"
      - label: "His156-Ser77"
        selection_a: "resid 156 and name NE2"
        selection_b: "resid 77 and name OG"

# Must have corresponding entry in comparison_settings
comparison_settings:
  catalytic_triad: {}
```

### Running Triad Comparison

`````{tab-set}
````{tab-item} CLI (Recommended)
```bash
# From your comparison project directory
polyzymd compare triad

# With options
polyzymd compare triad --eq-time 10ns --format markdown
polyzymd compare triad -o triad_report.md
```
````

````{tab-item} Python
```python
from polyzymd.compare import ComparisonConfig, TriadComparator

# Load configuration (must have analysis_settings.catalytic_triad section)
config = ComparisonConfig.from_yaml("comparison.yaml")

# Get triad settings from analysis_settings
triad_settings = config.analysis_settings.get("catalytic_triad")

# Run triad comparison
comparator = TriadComparator(
    config=config,
    triad_settings=triad_settings,
    equilibration="10ns",
)
result = comparator.compare()

# Print results
print(f"Best triad integrity: {result.ranking[0]}")
for cond in result.conditions:
    contact_pct = cond.mean_simultaneous_contact * 100
    sem_pct = cond.sem_simultaneous_contact * 100
    print(f"{cond.label}: {contact_pct:.1f}% ± {sem_pct:.1f}%")

# Save to JSON
result.save("results/triad_comparison.json")
```
````
`````

### Example Output

```
Catalytic Triad Comparison: polymer_stability_study
======================================================================
Triad: LipA_Ser-His-Asp
Description: Lipase A catalytic triad
Pairs: Asp133-His156, His156-Ser77
Contact threshold: 3.5 A
Equilibration: 10ns
Control: No Polymer

Condition Summary (ranked by simultaneous contact, highest first)
----------------------------------------------------------------------
Rank  Condition            Contact %    SEM        N
----------------------------------------------------------------------
1     100% SBMA               87.3%      2.15%  3
2     No Polymer              72.1%      3.42%  3   *
3     50/50 Mix               68.5%      4.21%  3
----------------------------------------------------------------------
* = control condition

Pairwise Comparisons
-------------------------------------------------------------------------------------
Comparison                     % Change   p-value      Cohen's d  Effect
-------------------------------------------------------------------------------------
100% SBMA vs No Polymer        +21.1%     0.0156*      2.89       large
50/50 Mix vs No Polymer        -5.0%      0.5234       -0.52      medium
-------------------------------------------------------------------------------------
* p < 0.05
Positive % change = improved triad contact

Interpretation
----------------------------------------------------------------------
Best triad integrity: 100% SBMA (87.3% simultaneous contact)
  -> 21.1% higher than control (No Polymer)
  -> Statistically significant (p=0.0156, d=2.89 [large])
```

### Per-Pair Distance Table

The output also includes a per-pair distance summary showing how each
individual H-bond distance compares across conditions:

```
Per-Pair Distances (Mean ± SEM across replicates)
------------------------------------------------------------------------------------------
Condition            Asp133-His156    His156-Ser77
------------------------------------------------------------------------------------------
100% SBMA            2.81±0.08        2.74±0.05
No Polymer           3.12±0.11        2.89±0.09
50/50 Mix            3.28±0.15        2.95±0.12
------------------------------------------------------------------------------------------
```

### Interpreting Results

| Contact % | Interpretation |
|-----------|---------------|
| > 90% | Excellent triad integrity |
| 70-90% | Good triad integrity |
| 50-70% | Moderate disruption |
| < 50% | Significant triad disruption |

**Key metrics:**
- **% Change**: Positive = more triad contact = better
- **p-value**: < 0.05 indicates statistically significant difference
- **Cohen's d**: Effect size magnitude

### CLI Reference for Triad

```bash
polyzymd compare triad [OPTIONS]

Options:
  -f, --file PATH                 Config file [default: comparison.yaml]
  --eq-time TEXT                  Override equilibration time
  --recompute                     Force recompute triad analysis
  --format [table|markdown|json]  Output format [default: table]
  -o, --output PATH               Save formatted output to file
  -q, --quiet                     Suppress INFO messages
  --debug                         Enable DEBUG logging
```

### Python API for Triad Comparison

```python
from polyzymd.compare import ComparisonConfig, TriadComparator, format_triad_result

# Load configuration
config = ComparisonConfig.from_yaml("comparison.yaml")

# Get triad settings from analysis_settings
triad_settings = config.analysis_settings.get("catalytic_triad")

# Run triad comparison
comparator = TriadComparator(
    config=config,
    triad_settings=triad_settings,
    equilibration="10ns",
)
result = comparator.compare()

# Access results
print(f"Best triad integrity: {result.ranking[0]}")
for cond in result.conditions:
    contact_pct = cond.mean_simultaneous_contact * 100
    print(f"{cond.label}: {contact_pct:.1f}% ± {cond.sem_simultaneous_contact*100:.1f}%")

# Format output
print(format_triad_result(result, format="markdown"))

# Save result
result.save("results/triad_comparison.json")
```

### Loading Saved Triad Results

```python
from polyzymd.compare import TriadComparisonResult

# Load from JSON
result = TriadComparisonResult.load("results/triad_comparison_my_study.json")

# Access condition data
control = result.get_condition("No Polymer")
print(f"Control contact: {control.mean_simultaneous_contact * 100:.1f}%")

# Get pairwise comparison
comp = result.get_comparison("100% SBMA")
if comp and comp.significant:
    print(f"SBMA significantly improves triad contact (p={comp.p_value:.4f})")
```

## Comparing Polymer-Protein Contacts

Compare **polymer-protein contact statistics** across conditions to understand how
different polymer compositions affect protein-polymer interactions.

```{note}
**New to contacts analysis?** Start with the [Contacts Quick Start](analysis_contacts_quickstart.md)
to run individual analyses, then return here to compare conditions.
```

### Key Metrics

The contacts comparison analyzes two aggregate metrics:

| Metric | Description | Higher means... |
|--------|-------------|-----------------|
| **Coverage** | % of protein residues contacted by polymer | More extensive binding |
| **Mean Contact Fraction** | Average % of frames each residue is in contact | Stronger/more persistent binding |

Additionally, **residence times by polymer type** show how long each polymer type
(e.g., SBMA vs EGMA) maintains contacts, revealing selectivity differences.

### Running Contacts Comparison

`````{tab-set}
````{tab-item} CLI (Recommended)
```bash
# Basic comparison
polyzymd compare contacts

# With custom equilibration time
polyzymd compare contacts --eq-time 10ns

# Override polymer selection (only SBMA monomers)
polyzymd compare contacts --polymer-selection "resname SBM"

# Different output formats
polyzymd compare contacts --format markdown -o contacts_report.md
polyzymd compare contacts --format json
```
````

````{tab-item} Python
```python
from polyzymd.compare import ComparisonConfig, ContactsComparator

# Load configuration
config = ComparisonConfig.from_yaml("comparison.yaml")

# Get contacts settings from analysis_settings and comparison_settings
analysis_settings = config.analysis_settings.get("contacts")
comparison_settings = config.comparison_settings.get("contacts")

# Run comparison
comparator = ContactsComparator(
    config=config,
    analysis_settings=analysis_settings,
    comparison_settings=comparison_settings,
    equilibration="10ns",
)
result = comparator.compare()

# Access results
print(f"Highest coverage: {result.ranking_by_coverage[0]}")
print(f"Highest contact: {result.ranking_by_contact_fraction[0]}")

for cond in result.conditions:
    print(f"{cond.label}: {cond.coverage_mean*100:.1f}% coverage, "
          f"{cond.contact_fraction_mean*100:.1f}% contact")

# Save result
result.save("results/contacts_comparison.json")
```
````
`````

### Optional: Contacts Configuration in comparison.yaml

Add a `contacts` section to `analysis_settings` and `comparison_settings`:

```yaml
name: "polymer_stability_study"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../noPoly_LipA_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% SBMA"
    config: "../SBMA_100_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% EGMA"
    config: "../EGMA_100_DMSO/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

# Configure contacts analysis in analysis_settings
analysis_settings:
  contacts:
    polymer_selection: "resname SBM EGM"  # MDAnalysis selection
    protein_selection: "protein"
    cutoff: 4.5                            # Angstroms
    contact_criteria: "heavy_atom"

# Comparison-specific parameters in comparison_settings
comparison_settings:
  contacts:
    fdr_alpha: 0.05                        # Benjamini-Hochberg FDR
    min_effect_size: 0.5                   # Cohen's d threshold
    top_residues: 10                       # Top residues to show
```

If no `contacts` section is provided, defaults are used.

### Handling Conditions Without Polymer

Conditions without polymer atoms (e.g., "No Polymer" controls) are **automatically
excluded** from contacts analysis since there's nothing to measure. The command
will warn you:

```
Note: 1 condition(s) auto-excluded (no polymer atoms): No Polymer
```

This is expected behavior. Comparisons will be made between polymer-containing
conditions only.

### Example Output

```
Polymer-Protein Contacts Comparison: polymer_ratio_study
================================================================================
Analysis: polymer_protein_contacts
Polymer selection: resname SBM EGM
Contact cutoff: 4.5 A
Contact criteria: heavy_atom
Equilibration: 10ns
Auto-excluded (no polymer): No Polymer

Condition Summary - Coverage (ranked, highest first)
--------------------------------------------------------------------------------
Rank  Condition                 Coverage     SEM        N   
--------------------------------------------------------------------------------
1     100% EGMA                     88.4%       0.55%  3    
2     25% SBMA / 75% EGMA           86.9%       1.44%  3    
3     50% SBMA / 50% EGMA           82.7%       0.66%  3    
4     75% SBMA / 25% EGMA           82.7%       1.57%  3    
5     100% SBMA                     74.9%       0.28%  2    
--------------------------------------------------------------------------------

Condition Summary - Mean Contact Fraction (ranked, highest first)
--------------------------------------------------------------------------------
Rank  Condition                 Contact %    SEM        N   
--------------------------------------------------------------------------------
1     75% SBMA / 25% EGMA           30.2%       0.50%  3    
2     25% SBMA / 75% EGMA           29.1%       5.35%  3    
3     100% EGMA                     25.3%       2.64%  3    
4     100% SBMA                     22.9%       1.47%  2    
5     50% SBMA / 50% EGMA           22.9%       2.18%  3    
--------------------------------------------------------------------------------

Residence Time by Polymer Type (frames)
--------------------------------------------------------------------------------
Condition                          EGM          SBM
--------------------------------------------------------------------------------
100% SBMA                           --  10.0±0.2 
75% SBMA / 25% EGMA         7.8±0.2    9.3±0.0 
50% SBMA / 50% EGMA         7.3±0.3    8.6±0.5 
25% SBMA / 75% EGMA         7.1±0.7   10.5±0.4 
100% EGMA                   7.1±0.4            --
--------------------------------------------------------------------------------

Aggregate Comparisons
-----------------------------------------------------------------------------------------------
Comparison                     Metric          % Change   p-value      Cohen d    Effect      
-----------------------------------------------------------------------------------------------
100% EGMA vs 100% SBMA         coverage        +18.1%     0.0004*      -16.64     large       
100% EGMA vs 100% SBMA         mean contact f  +10.7%     0.5445       -0.62      medium      
...
-----------------------------------------------------------------------------------------------
* p < 0.05; positive % change = more contact in treatment

One-way ANOVA
------------------------------------------------------------
Metric                    F-stat       p-value      Significant 
------------------------------------------------------------
coverage                  18.323       0.0002       Yes*        
mean contact fraction     1.200        0.3748       No          
------------------------------------------------------------
```

### Interpreting Results

**Coverage rankings:**
- Higher coverage = polymer interacts with more of the protein surface
- 100% EGMA shows highest coverage (88.4%) - broader but possibly weaker binding

**Contact fraction rankings:**
- Higher mean contact = more persistent interactions per residue
- 75% SBMA / 25% EGMA shows highest contact fraction (30.2%) - more stable binding

**Residence time by polymer type:**
- SBMA (SBM) tends to have longer residence times than EGMA (EGM)
- This suggests SBMA forms more persistent interactions
- Useful for understanding polymer selectivity

**Statistical tests:**
- ANOVA tests whether any condition differs overall
- Pairwise comparisons with Benjamini-Hochberg FDR correction
- Cohen's d quantifies effect magnitude independent of sample size

### CLI Reference for Contacts

```bash
polyzymd compare contacts [OPTIONS]

Options:
  -f, --file PATH                 Config file [default: comparison.yaml]
  --eq-time TEXT                  Override equilibration time
  --polymer-selection TEXT        Override polymer selection (MDAnalysis syntax)
  --cutoff FLOAT                  Override contact cutoff (Angstroms)
  --fdr-alpha FLOAT               FDR alpha for multiple testing correction
  --recompute                     Force recompute contacts analysis
  --format [table|markdown|json]  Output format [default: table]
  -o, --output PATH               Save formatted output to file
  -q, --quiet                     Suppress INFO messages
  --debug                         Enable DEBUG logging
```

### Python API for Contacts Comparison

```python
from polyzymd.compare import (
    ComparisonConfig,
    ContactsComparator,
    format_contacts_result,
)

# Load configuration
config = ComparisonConfig.from_yaml("comparison.yaml")

# Get contacts settings from analysis_settings and comparison_settings
analysis_settings = config.analysis_settings.get("contacts")
comparison_settings = config.comparison_settings.get("contacts")

# Run comparison
comparator = ContactsComparator(
    config=config,
    analysis_settings=analysis_settings,
    comparison_settings=comparison_settings,
    equilibration="10ns",
)
result = comparator.compare()

# Access results
print(f"Highest coverage: {result.ranking_by_coverage[0]}")
print(f"Highest contact: {result.ranking_by_contact_fraction[0]}")

for cond in result.conditions:
    print(f"{cond.label}: {cond.coverage_mean*100:.1f}% coverage, "
          f"{cond.contact_fraction_mean*100:.1f}% contact")
    
    # Residence time by polymer type
    for poly_type, (mean, sem) in cond.residence_time_by_polymer_type.items():
        print(f"  {poly_type}: {mean:.1f} ± {sem:.1f} frames")

# Format output
print(format_contacts_result(result, format="markdown"))

# Save result
result.save("results/contacts_comparison.json")
```

### Loading Saved Contacts Results

```python
from polyzymd.compare import ContactsComparisonResult

# Load from JSON
result = ContactsComparisonResult.load("results/contacts_comparison_my_study.json")

# Access condition data
for cond in result.conditions:
    print(f"{cond.label}: coverage={cond.coverage_mean*100:.1f}%")

# Get aggregate comparisons
for comp in result.aggregate_comparisons:
    if comp.significant:
        print(f"{comp.condition_a} vs {comp.condition_b} ({comp.metric}): "
              f"p={comp.p_value:.4f}")
```

### Contacts vs RMSF: Complementary Analyses

| Analysis | Question Answered |
|----------|-------------------|
| **RMSF** | Does polymer stabilize the enzyme (reduce flexibility)? |
| **Contacts** | Where and how strongly does polymer bind? |
| **Combined** | Do contact hotspots correlate with stabilization? |

For mechanistic insights correlating contacts with flexibility changes, see
`polyzymd compare report` (coming soon).

## Comparing Distances Across Conditions

Compare **inter-atomic distances** across conditions with statistical analysis.
This is useful for tracking specific interactions (e.g., substrate proximity,
hydrogen bond distances) that may change with different polymer environments.

```{note}
**New to distance analysis?** Start with the [Distance Analysis Quick Start](analysis_distances_quickstart.md)
to understand distance pair definitions and selection syntax, then return here
to compare conditions.
```

### Key Metrics

The distances comparison provides **dual-metric ranking**:

| Metric | Description | Ranking |
|--------|-------------|---------|
| **Mean Distance** | Average distance across trajectory (primary) | Lowest first (closer = better) |
| **Fraction Below Threshold** | % of frames below contact threshold (secondary) | Highest first (more contact = better) |

This dual approach captures both the typical distance AND the frequency of close contacts.

### Running Distances Comparison

`````{tab-set}
````{tab-item} CLI (Recommended)
```bash
# Basic comparison
polyzymd compare run distances -f comparison.yaml

# With custom equilibration time
polyzymd compare run distances -f comparison.yaml --eq-time 10ns

# Different output formats
polyzymd compare run distances -f comparison.yaml --format markdown
polyzymd compare run distances -f comparison.yaml --format json -o distances_report.json
```
````

````{tab-item} Python
```python
from polyzymd.compare import ComparisonConfig
from polyzymd.compare.comparators import DistancesComparator

# Load configuration
config = ComparisonConfig.from_yaml("comparison.yaml")

# Get distances settings from analysis_settings
distances_settings = config.analysis_settings.get("distances")

# Run comparison
comparator = DistancesComparator(
    config=config,
    analysis_settings=distances_settings,
    equilibration="10ns",
)
result = comparator.compare()

# Access results
print(f"Closest mean distance: {result.ranking[0]}")
if result.ranking_by_fraction:
    print(f"Highest contact fraction: {result.ranking_by_fraction[0]}")

for cond in result.conditions:
    print(f"{cond.label}: {cond.overall_mean_distance:.2f} A")
    if cond.overall_fraction_below is not None:
        print(f"  Contact fraction: {cond.overall_fraction_below*100:.1f}%")

# Save result
result.save("results/distances_comparison.json")
```
````
`````

### Adding Distances to comparison.yaml

Add a `distances` section to your `analysis_settings`:

```yaml
name: "substrate_proximity_study"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../noPoly_LipA_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% SBMA"
    config: "../SBMA_100_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% EGMA"
    config: "../EGMA_100_DMSO/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

# Define distance pairs in analysis_settings
analysis_settings:
  distances:
    threshold: 3.5  # Global default threshold (Angstroms, optional)
    pairs:
      - label: "Catalytic H-bond"
        selection_a: "resid 77 and name OG"
        selection_b: "resid 133 and name NE2"
        threshold: 3.5  # Per-pair threshold (overrides global)
      - label: "Lid Domain Opening"
        selection_a: "com(resid 141-148)"
        selection_b: "com(resid 281-289)"
        threshold: 15.0  # Different threshold for this pair

# Must have corresponding entry in comparison_settings
comparison_settings:
  distances: {}
```

### Per-Pair Thresholds

Different distance pairs often have different biologically relevant thresholds:

| Type of Distance | Typical Threshold |
|------------------|-------------------|
| Hydrogen bond | 3.0 - 3.5 A |
| Salt bridge | 4.0 - 4.5 A |
| Aromatic stacking | 4.0 - 5.0 A |
| Domain separation | 10 - 20 A |
| Lid opening | 15 - 25 A |

**Threshold resolution order:**
1. Per-pair `threshold` in the pair definition (highest priority)
2. Global `threshold` in the `distances` section (fallback)
3. No threshold (fraction below not computed)

**Example with mixed thresholds:**

```yaml
analysis_settings:
  distances:
    threshold: 4.0  # Default for pairs without explicit threshold
    pairs:
      - label: "Ser77-His133 H-bond"
        selection_a: "resid 77 and name OG"
        selection_b: "resid 133 and name NE2"
        threshold: 3.5  # H-bond cutoff

      - label: "Asp156-His133 H-bond"
        selection_a: "midpoint(resid 156 and name OD1 OD2)"
        selection_b: "resid 133 and name ND1"
        # Uses global threshold: 4.0 A

      - label: "Lid-to-Core Distance"
        selection_a: "com(resid 141-148)"
        selection_b: "com(resid 281-289)"
        threshold: 15.0  # Large-scale motion threshold
```

**Threshold cache invalidation:** When you change a threshold value, PolyzyMD
automatically detects the mismatch and recomputes contact fractions from the
stored per-replicate distance data. This avoids expensive trajectory
reprocessing - only the statistical aggregation is recalculated.

### Example Output

```
Distance Comparison: substrate_proximity_study
================================================================================
Pairs analyzed: 2
Pair labels: Ser77-Substrate, His156-Substrate
Contact threshold: 3.5 A
Equilibration: 10ns
Control: No Polymer

Condition Summary (ranked by mean distance, lowest first)
--------------------------------------------------------------------------------
Rank  Condition                 Mean Dist    SEM        % Below    N   
--------------------------------------------------------------------------------
1     100% SBMA                 7.46 A       0.421         0.4%   3    
2     100% EGMA                 7.66 A       0.270         1.6%   3    
3     No Polymer                8.02 A       0.315         0.0%   3   *
--------------------------------------------------------------------------------
* = control condition

Secondary Ranking (by % below threshold, highest first)
------------------------------------------------------------
1     100% EGMA                 1.6% (SEM: 0.31%)
2     100% SBMA                 0.4% (SEM: 0.31%)
3     No Polymer                0.0% (SEM: 0.00%)

Per-Pair Distances (Mean +/- SEM across replicates)
------------------------------------------------------------------------------------------
Condition                 Ser77-Substrate His156-Substrate
------------------------------------------------------------------------------------------
100% SBMA                 7.46+/-0.42     6.12+/-0.25    
100% EGMA                 7.66+/-0.27     6.45+/-0.18    
No Polymer                8.02+/-0.31     6.89+/-0.22    
------------------------------------------------------------------------------------------

Pairwise Comparisons (Distance Metric)
------------------------------------------------------------------------------------------
Comparison                     % Change   p-value      Cohen d    Effect       Direction 
------------------------------------------------------------------------------------------
100% SBMA vs No Polymer        -7.0%      0.0451*      0.87       large        closer    
100% EGMA vs No Polymer        -4.4%      0.1234       0.70       medium       closer    
------------------------------------------------------------------------------------------
* p < 0.05
Negative % change = lower distance (closer)

Pairwise Comparisons (Fraction Below Threshold)
------------------------------------------------------------------------------------------
Comparison                     % Change   p-value      Cohen d    Effect       Direction   
------------------------------------------------------------------------------------------
100% SBMA vs No Polymer        +0.4%      0.2161       -1.20      large        more_contact
100% EGMA vs No Polymer        +1.6%      0.0065*      -4.25      large        more_contact
------------------------------------------------------------------------------------------
* p < 0.05
Positive % change = more frames below threshold (more contact)

One-way ANOVA
--------------------------------------------------
Distance metric:
  F-statistic: 3.901
  p-value:     0.0512
  Significant: No (alpha=0.05)
Fraction metric:
  F-statistic: 5.880
  p-value:     0.0171
  Significant: Yes (alpha=0.05)

Interpretation
--------------------------------------------------------------------------------
Closest mean distance: 100% SBMA (7.46 A)
  -> 7.0% closer than control (No Polymer)
  -> Statistically significant (p=0.0451, d=0.87 [large])

Highest contact fraction: 100% EGMA (1.6% below threshold)

Analysis completed: 2026-02-16 21:03:40
PolyzyMD version: 1.0.0
```

### Interpreting Results

**Primary metric (Mean Distance):**
- Lower distance = closer = atoms more frequently in proximity
- Ranking from lowest to highest (Rank 1 = closest)
- Negative % change vs control = improvement (closer)

**Secondary metric (Fraction Below Threshold):**
- Only computed if `threshold` is specified in config
- Higher fraction = more frames with close contact
- Ranking from highest to lowest (Rank 1 = most contact)
- Positive % change vs control = improvement (more contact)

**Why dual metrics?**
- Mean distance captures typical behavior
- Fraction below threshold captures extreme events (e.g., catalytic encounters)
- A condition might have similar mean distance but more frequent close approaches

### Statistical Analysis

Both metrics undergo independent statistical testing:

| Test | Applied To | Interpretation |
|------|-----------|----------------|
| **t-test** | Each condition vs control | p < 0.05 = significant difference |
| **Cohen's d** | Each comparison | Effect magnitude (regardless of p-value) |
| **ANOVA** | All conditions | Any condition differs? (3+ conditions) |

**Effect size interpretation:**

| Cohen's d | Interpretation |
|-----------|---------------|
| < 0.2 | Negligible |
| 0.2 - 0.5 | Small |
| 0.5 - 0.8 | Medium |
| > 0.8 | Large |

### CLI Reference for Distances

```bash
polyzymd compare run distances [OPTIONS]

Options:
  -f, --file PATH                 Config file [default: comparison.yaml]
  --eq-time TEXT                  Override equilibration time
  --recompute                     Force recompute distance analysis
  --format [table|markdown|json]  Output format [default: table]
  -o, --output PATH               Save formatted output to file
  -q, --quiet                     Suppress INFO messages
  --debug                         Enable DEBUG logging
```

### Python API for Distances Comparison

```python
from polyzymd.compare import ComparisonConfig
from polyzymd.compare.comparators import DistancesComparator
from polyzymd.compare.distances_formatters import format_distances_result

# Load configuration
config = ComparisonConfig.from_yaml("comparison.yaml")

# Get distances settings
distances_settings = config.analysis_settings.get("distances")

# Run comparison
comparator = DistancesComparator(
    config=config,
    analysis_settings=distances_settings,
    equilibration="10ns",
)
result = comparator.compare()

# Access primary ranking (by mean distance)
print(f"Closest: {result.ranking[0]}")
for cond in result.conditions:
    print(f"{cond.label}: {cond.overall_mean_distance:.2f} ± {cond.overall_sem_distance:.3f} A")

# Access secondary ranking (by fraction below threshold)
if result.ranking_by_fraction:
    print(f"\nHighest contact: {result.ranking_by_fraction[0]}")
    for cond in result.conditions:
        if cond.overall_fraction_below is not None:
            print(f"{cond.label}: {cond.overall_fraction_below*100:.1f}%")

# Access per-pair details
for cond in result.conditions:
    print(f"\n{cond.label}:")
    for pair in cond.pair_summaries:
        print(f"  {pair.label}: {pair.mean_distance:.2f} ± {pair.sem_distance:.2f} A")

# Format output
print(format_distances_result(result, format="markdown"))

# Save result
result.save("results/distances_comparison.json")
```

### Loading Saved Distance Results

```python
from polyzymd.compare.results import DistanceComparisonResult

# Load from JSON
result = DistanceComparisonResult.load("results/distances_comparison_my_study.json")

# Access condition data
control = result.get_condition("No Polymer")
print(f"Control mean distance: {control.overall_mean_distance:.2f} A")

# Get pairwise comparison
comp = result.get_comparison("100% SBMA")
if comp and comp.distance_significant:
    print(f"SBMA significantly closer (p={comp.distance_p_value:.4f})")
```

### Use Cases for Distance Comparison

| Use Case | Configuration |
|----------|---------------|
| **Substrate binding** | Distance from catalytic residues to substrate atoms |
| **Active site geometry** | Similar to triad, but for non-catalytic interactions |
| **Polymer-residue proximity** | Distance from polymer termini to specific residues |
| **Conformational changes** | Distance between domains or loops |

### Distances vs Catalytic Triad Comparison

| Feature | `compare run distances` | `compare run triad` |
|---------|-------------------------|---------------------|
| **Metric** | Mean distance + fraction | Simultaneous contact fraction |
| **Pairs** | Any atom pairs | Pre-defined triad geometry |
| **Ranking** | Dual (distance + fraction) | Single (simultaneous contact) |
| **Use case** | General distance tracking | Catalytic geometry integrity |

```{tip}
Use **distances** for monitoring specific interactions with dual-metric analysis.
Use **triad** when all pairs must be in contact simultaneously (catalytic triad geometry).
```

## See Also

- [RMSF Quick Start](analysis_rmsf_quickstart.md) -- Run individual RMSF analysis
- [Contacts Quick Start](analysis_contacts_quickstart.md) -- Run individual contacts analysis
- [Distance Analysis Quick Start](analysis_distances_quickstart.md) -- Run individual distance analysis
- [Catalytic Triad Analysis](analysis_triad_quickstart.md) -- Run individual triad analysis
- [Statistical Best Practices](analysis_rmsf_best_practices.md) -- Understanding statistics
- [Reference Selection](analysis_reference_selection.md) -- Alignment options
