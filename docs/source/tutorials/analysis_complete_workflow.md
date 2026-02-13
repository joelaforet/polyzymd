# Complete Analysis Workflow: From Finished Simulations to Publication

Your simulations just finished. Now what?

This guide walks you through the **complete PolyzyMD analysis workflow**—from
raw trajectories to publication-ready figures comparing multiple conditions.
By the end, you'll have analyzed protein flexibility (RMSF), active site
geometry (catalytic triad), inter-atomic distances, and polymer-protein 
contacts across all your simulation conditions.

```{note}
**Prerequisites:**
- Completed production simulations with trajectory files (`.dcd`, `.xtc`)
- A `config.yaml` file for each simulation condition
- `solvated_system.pdb` topology in each run directory (created during build)

**Time estimate:** 30-60 minutes for a typical 3-condition study
```

## Overview: The Analysis Pipeline

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        YOUR SIMULATION DATA                              │
│  condition_A/config.yaml    condition_B/config.yaml    condition_C/...  │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STEP 1: Create analysis.yaml for each condition                        │
│  Define: replicates, equilibration time, analysis parameters            │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STEP 2: Run individual analyses                                         │
│  polyzymd analyze run  →  RMSF, triad, contacts per condition           │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STEP 3: Create comparison.yaml to compare conditions                   │
│  Define: conditions, control, comparison parameters                      │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STEP 4: Run statistical comparisons                                     │
│  polyzymd compare rmsf/triad/contacts  →  p-values, effect sizes        │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│  STEP 5: Generate publication figures                                    │
│  polyzymd compare plot  →  bar charts, forest plots, summary panels     │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Step 1: Configure Your Analyses

### 1.1 Understand Your Project Structure

A typical PolyzyMD simulation campaign looks like this:

```
my_enzyme_study/
├── noPoly_enzyme_DMSO/           # Control: no polymer
│   ├── config.yaml
│   └── scratch/
│       ├── solvated_system.pdb
│       ├── rep_1/production_1.dcd
│       ├── rep_2/production_1.dcd
│       └── rep_3/production_1.dcd
│
├── SBMA_100_enzyme_DMSO/         # Condition: 100% SBMA polymer
│   ├── config.yaml
│   └── scratch/...
│
├── EGMA_100_enzyme_DMSO/         # Condition: 100% EGMA polymer
│   ├── config.yaml
│   └── scratch/...
│
└── analysis/                      # We'll create this
    └── comparison.yaml
```

### 1.2 Create analysis.yaml for Each Condition

Each simulation condition needs an `analysis.yaml` file that lives **alongside**
its `config.yaml`. This file defines which analyses to run and their parameters.

**Navigate to your first condition:**

```bash
cd noPoly_enzyme_DMSO/
```

**Initialize a template:**

```bash
polyzymd analyze init
```

This creates an `analysis.yaml` with default values. Now customize it:

```yaml
# analysis.yaml - No Polymer Control
# ===================================

# Which replicates to analyze
replicates: [1, 2, 3]

# Shared parameters
defaults:
  equilibration_time: "10ns"    # Skip first 10ns for equilibration

# ─────────────────────────────────────────────────────────────────────────
# RMSF Analysis: Protein flexibility
# ─────────────────────────────────────────────────────────────────────────
rmsf:
  enabled: true
  selection: "protein and name CA"    # Alpha carbons only
  reference_mode: "average"           # Align to average structure

# ─────────────────────────────────────────────────────────────────────────
# Catalytic Triad Analysis: Active site geometry
# ─────────────────────────────────────────────────────────────────────────
catalytic_triad:
  enabled: true
  name: "Ser-His-Asp"
  threshold: 3.5                      # H-bond distance cutoff (Angstroms)
  pairs:
    # Ser77 OG ↔ His156 NE2 (nucleophile-histidine)
    - label: "Ser77-His156"
      selection_a: "resid 77 and name OG"
      selection_b: "resid 156 and name NE2"
    # His156 ND1 ↔ Asp133 OD1/OD2 midpoint (histidine-aspartate)
    - label: "His156-Asp133"
      selection_a: "resid 156 and name ND1"
      selection_b: "midpoint(resid 133 and name OD1 OD2)"

# ─────────────────────────────────────────────────────────────────────────
# Distance Analysis: Specific atom-pair distances
# ─────────────────────────────────────────────────────────────────────────
distances:
  enabled: true
  pairs:
    # Substrate to catalytic serine (reaction distance)
    - label: "Substrate-Ser77"
      selection_a: "resname SUB and name C1"      # Replace with your substrate
      selection_b: "resid 77 and name OG"
    # Histidine to substrate (proton shuttle)
    - label: "His156-Substrate"
      selection_a: "resid 156 and name NE2"
      selection_b: "resname SUB and name O1"

# ─────────────────────────────────────────────────────────────────────────
# Contacts Analysis: Polymer-protein interactions
# ─────────────────────────────────────────────────────────────────────────
# NOTE: For the no-polymer control, contacts analysis will be skipped
# automatically (no polymer atoms to analyze). Enable it anyway for
# consistency across conditions.
contacts:
  enabled: true
  polymer_selection: "chainID C"      # PolyzyMD chain convention
  protein_selection: "protein"
  cutoff: 4.5                         # Contact distance (Angstroms)
  compute_residence_times: true       # Track contact durations
```

```{tip}
**Finding your catalytic triad residues:**
Check your enzyme's literature or use PyMOL/VMD to identify the catalytic
residues. Common triads:
- **Lipases:** Ser-His-Asp (e.g., Ser77, His156, Asp133)
- **Proteases:** Ser-His-Asp or Cys-His-Asn
- **Esterases:** Ser-Glu-His
```

### 1.3 Copy to Other Conditions

Copy the same `analysis.yaml` to your other conditions:

```bash
cp analysis.yaml ../SBMA_100_enzyme_DMSO/
cp analysis.yaml ../EGMA_100_enzyme_DMSO/
```

```{note}
The `analysis.yaml` can be **identical** across conditions. The differences
come from each condition's `config.yaml` pointing to different trajectories.
```

---

## Step 2: Run Individual Analyses

### 2.1 Run All Analyses for Each Condition

For each condition, run the complete analysis suite:

`````{tab-set}
````{tab-item} One at a time
```bash
# Condition 1: No Polymer (control)
cd noPoly_enzyme_DMSO/
polyzymd analyze run

# Condition 2: SBMA
cd ../SBMA_100_enzyme_DMSO/
polyzymd analyze run

# Condition 3: EGMA
cd ../EGMA_100_enzyme_DMSO/
polyzymd analyze run
```
````

````{tab-item} Loop (bash)
```bash
for condition in noPoly_enzyme_DMSO SBMA_100_enzyme_DMSO EGMA_100_enzyme_DMSO; do
    echo "Analyzing $condition..."
    cd $condition
    polyzymd analyze run
    cd ..
done
```
````

````{tab-item} Parallel (GNU parallel)
```bash
# Run all conditions in parallel (if you have the cores)
parallel -j3 "cd {} && polyzymd analyze run" ::: \
    noPoly_enzyme_DMSO \
    SBMA_100_enzyme_DMSO \
    EGMA_100_enzyme_DMSO
```
````
`````

**Expected output:**

```
Loading configuration from: config.yaml
Analysis: noPoly_enzyme_DMSO

Enabled analyses: rmsf, catalytic_triad, contacts

[1/3] RMSF Analysis
  Replicates: 1, 2, 3
  Selection: protein and name CA
  Reference mode: average
  Processing replicate 1... done (181 residues, mean RMSF: 0.72 Å)
  Processing replicate 2... done (181 residues, mean RMSF: 0.69 Å)
  Processing replicate 3... done (181 residues, mean RMSF: 0.70 Å)
  Results saved: analysis/rmsf/

[2/3] Catalytic Triad Analysis
  Triad: Ser-His-Asp
  Pairs: Ser77-His156, His156-Asp133
  Threshold: 3.5 Å
  Processing replicate 1... done (72.1% simultaneous contact)
  Processing replicate 2... done (74.3% simultaneous contact)
  Processing replicate 3... done (71.8% simultaneous contact)
  Results saved: analysis/triad/

[3/3] Contacts Analysis
  Note: No polymer atoms found (chainID C). Skipping.

Analysis complete! Results in: analysis/
```

### 2.2 Verify Results

Check that results were created:

```bash
ls -la analysis/
```

```
analysis/
├── rmsf/
│   ├── rmsf_rep1.json
│   ├── rmsf_rep2.json
│   └── rmsf_rep3.json
├── triad/
│   ├── triad_rep1.json
│   ├── triad_rep2.json
│   └── triad_rep3.json
└── contacts/           # Empty for no-polymer control
```

### 2.3 Quick Inspection (Optional)

Preview individual results before comparison:

`````{tab-set}
````{tab-item} RMSF
```bash
# View RMSF summary for one replicate
polyzymd analyze rmsf -c config.yaml -r 1 --eq-time 10ns
```
````

````{tab-item} Triad
```bash
# View triad geometry for one replicate
polyzymd analyze triad -c config.yaml -r 1 --eq-time 10ns
```
````

````{tab-item} Distances
```bash
# View distances for one replicate
polyzymd analyze distances -c config.yaml -r 1 --eq-time 10ns \
    --pair "resname SUB and name C1 : resid 77 and name OG"
```
````

````{tab-item} Contacts
```bash
# View contacts for one replicate (polymer conditions only)
polyzymd analyze contacts -c config.yaml -r 1 --eq-time 10ns
```
````
`````

---

## Step 3: Set Up Cross-Condition Comparison

Now that individual analyses are complete, we'll compare conditions statistically.

### 3.1 Create a Comparison Project

```bash
# Go to your study root directory
cd my_enzyme_study/

# Create comparison project
polyzymd compare init polymer_stabilization_study
```

This creates:

```
polymer_stabilization_study/
├── comparison.yaml    # Template to edit
├── results/           # Comparison outputs (auto-populated)
└── figures/           # Publication figures (auto-populated)
```

### 3.2 Configure comparison.yaml

Edit `comparison.yaml` to define your conditions:

```yaml
# comparison.yaml - Polymer Stabilization Study
# ==============================================

name: "polymer_stabilization_study"
description: "Effect of SBMA vs EGMA polymer conjugation on enzyme stability"

# Control condition for relative comparisons
control: "No Polymer"

# ─────────────────────────────────────────────────────────────────────────
# Conditions to Compare
# ─────────────────────────────────────────────────────────────────────────
conditions:
  - label: "No Polymer"
    config: "../noPoly_enzyme_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% SBMA"
    config: "../SBMA_100_enzyme_DMSO/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% EGMA"
    config: "../EGMA_100_enzyme_DMSO/config.yaml"
    replicates: [1, 2, 3]

# ─────────────────────────────────────────────────────────────────────────
# Default Parameters
# ─────────────────────────────────────────────────────────────────────────
defaults:
  equilibration_time: "10ns"

# ─────────────────────────────────────────────────────────────────────────
# RMSF Comparison Settings
# ─────────────────────────────────────────────────────────────────────────
rmsf:
  selection: "protein and name CA"
  reference_mode: "average"

# ─────────────────────────────────────────────────────────────────────────
# Catalytic Triad Comparison Settings
# ─────────────────────────────────────────────────────────────────────────
catalytic_triad:
  name: "Ser-His-Asp"
  threshold: 3.5
  pairs:
    - label: "Ser77-His156"
      selection_a: "resid 77 and name OG"
      selection_b: "resid 156 and name NE2"
    - label: "His156-Asp133"
      selection_a: "resid 156 and name ND1"
      selection_b: "midpoint(resid 133 and name OD1 OD2)"

# ─────────────────────────────────────────────────────────────────────────
# Contacts Comparison Settings
# ─────────────────────────────────────────────────────────────────────────
contacts:
  polymer_selection: "chainID C"
  protein_selection: "protein"
  cutoff: 4.5
```

### 3.3 Validate Your Configuration

Before running comparisons, verify all paths are correct:

```bash
cd polymer_stabilization_study/
polyzymd compare validate
```

**Expected output:**

```
Validating: comparison.yaml

Validation PASSED

  Name: polymer_stabilization_study
  Description: Effect of SBMA vs EGMA polymer conjugation on enzyme stability
  Control: No Polymer
  
  Conditions: 3
    ✓ No Polymer: 3 replicates (config found)
    ✓ 100% SBMA: 3 replicates (config found)
    ✓ 100% EGMA: 3 replicates (config found)
  
  Analyses configured:
    ✓ rmsf: protein and name CA
    ✓ catalytic_triad: Ser-His-Asp (2 pairs)
    ✓ contacts: chainID C ↔ protein
```

---

## Step 4: Run Statistical Comparisons

### 4.1 Compare RMSF (Protein Flexibility)

`````{tab-set}
````{tab-item} CLI
```bash
polyzymd compare rmsf
```
````

````{tab-item} Python
```python
from polyzymd.compare import ComparisonConfig, RMSFComparator

config = ComparisonConfig.from_yaml("comparison.yaml")
comparator = RMSFComparator(config, rmsf_config=config.rmsf, equilibration="10ns")
result = comparator.compare()

print(f"Most stable: {result.ranking[0]}")
result.save("results/rmsf_comparison.json")
```
````
`````

**Example output:**

```
RMSF Comparison: polymer_stabilization_study
============================================================
Selection: protein and name CA
Equilibration: 10ns
Control: No Polymer

Condition Summary (ranked by RMSF, lowest first)
------------------------------------------------------------
Rank  Condition            Mean RMSF    SEM        N
------------------------------------------------------------
1     100% SBMA               0.551 Å    0.034  3
2     100% EGMA               0.597 Å    0.073  3
3     No Polymer              0.715 Å    0.020  3   *
------------------------------------------------------------
* = control condition

Pairwise Comparisons vs Control
--------------------------------------------------------------------------------
Comparison                     % Change   p-value      Cohen's d  Effect
--------------------------------------------------------------------------------
100% SBMA vs No Polymer        -22.9%     0.0211*      4.06       large
100% EGMA vs No Polymer        -16.4%     0.1944       1.27       large
--------------------------------------------------------------------------------
* p < 0.05

Results saved: results/rmsf_comparison_polymer_stabilization_study.json
```

**Interpretation:**
- 100% SBMA significantly reduces RMSF by 22.9% (p < 0.05)
- 100% EGMA shows a large effect (Cohen's d = 1.27) but isn't statistically
  significant with N=3 replicates

### 4.2 Compare Catalytic Triad Geometry

`````{tab-set}
````{tab-item} CLI
```bash
polyzymd compare triad
```
````

````{tab-item} Python
```python
from polyzymd.compare import ComparisonConfig, TriadComparator

config = ComparisonConfig.from_yaml("comparison.yaml")
comparator = TriadComparator(config, equilibration="10ns")
result = comparator.compare()

print(f"Best triad integrity: {result.ranking[0]}")
result.save("results/triad_comparison.json")
```
````
`````

**Example output:**

```
Catalytic Triad Comparison: polymer_stabilization_study
======================================================================
Triad: Ser-His-Asp
Pairs: Ser77-His156, His156-Asp133
Threshold: 3.5 Å
Control: No Polymer

Condition Summary (ranked by simultaneous contact, highest first)
----------------------------------------------------------------------
Rank  Condition            Contact %    SEM        N
----------------------------------------------------------------------
1     100% SBMA               87.3%      2.15%  3
2     No Polymer              72.1%      3.42%  3   *
3     100% EGMA               68.5%      4.21%  3
----------------------------------------------------------------------

Pairwise Comparisons
---------------------------------------------------------------------------
Comparison                     % Change   p-value      Cohen's d  Effect
---------------------------------------------------------------------------
100% SBMA vs No Polymer        +21.1%     0.0156*      2.89       large
100% EGMA vs No Polymer        -5.0%      0.5234       -0.52      medium
---------------------------------------------------------------------------
Positive % change = improved triad integrity

Results saved: results/triad_comparison_polymer_stabilization_study.json
```

**Interpretation:**
- 100% SBMA significantly improves triad integrity by 21.1% (p < 0.05)
- 100% EGMA slightly reduces triad integrity but not significantly

### 4.3 Analyze Distances (Single Condition)

Distance analysis measures specific atom-pair distances over time. Unlike the
comparison-focused RMSF and triad analyses, distances are typically computed
per-condition and then compared manually or via custom scripts.

`````{tab-set}
````{tab-item} CLI
```bash
# Single replicate with threshold analysis
polyzymd analyze distances -c config.yaml -r 1 --eq-time 10ns \
    --pair "resname SUB and name C1 : resid 77 and name OG" \
    --threshold 3.5

# Multiple replicates (aggregated with SEM)
polyzymd analyze distances -c config.yaml -r 1-3 --eq-time 10ns \
    --pair "resname SUB and name C1 : resid 77 and name OG" \
    --pair "resid 156 and name NE2 : resname SUB and name O1" \
    --threshold 3.5 --plot
```
````

````{tab-item} Python
```python
from polyzymd.config.schema import SimulationConfig
from polyzymd.analysis import DistanceCalculator

config = SimulationConfig.from_yaml("config.yaml")

# Define distance pairs
pairs = [
    ("resname SUB and name C1", "resid 77 and name OG"),  # Substrate-Ser
    ("resid 156 and name NE2", "resname SUB and name O1"),  # His-Substrate
]

# Create calculator with threshold for contact analysis
calc = DistanceCalculator(
    config=config,
    pairs=pairs,
    equilibration="10ns",
    threshold=3.5,  # Contact cutoff
)

# Single replicate
result = calc.compute(replicate=1)
for pr in result.pair_results:
    print(f"{pr.pair_label}: {pr.mean_distance:.2f} ± {pr.sem_distance:.2f} Å")
    if pr.fraction_below_threshold is not None:
        print(f"  Contact fraction: {pr.fraction_below_threshold:.1%}")

# Aggregate across replicates
agg = calc.compute_aggregated(replicates=[1, 2, 3])
```
````
`````

**Example output:**

```
Distance Analysis Complete (Aggregated)
  Replicates: 1-3
  resname_SUB_C1-resid77_OG:
    Mean: 4.21 ± 0.18 Å (SEM across 3 replicates)
    Contact fraction (<3.5Å): 32.4% ± 4.2%
  resid156_NE2-resname_SUB_O1:
    Mean: 3.89 ± 0.25 Å
    Contact fraction (<3.5Å): 48.7% ± 5.1%
```

```{tip}
**Distance vs. Triad analysis:** Use **triad** when you need simultaneous
contact analysis (all pairs below threshold at once). Use **distances** when
you want individual pair statistics, KDE-based mode estimation, or custom
atom selections beyond the catalytic triad.
```

For the full distance analysis guide, see [Distance Analysis Quick Start](analysis_distances_quickstart.md).

### 4.4 Compare Polymer-Protein Contacts

`````{tab-set}
````{tab-item} CLI
```bash
polyzymd compare contacts
```
````

````{tab-item} Python
```python
from polyzymd.compare import ComparisonConfig, ContactsComparator

config = ComparisonConfig.from_yaml("comparison.yaml")
comparator = ContactsComparator(config, equilibration="10ns")
result = comparator.compare()

print(f"Highest coverage: {result.ranking_by_coverage[0]}")
result.save("results/contacts_comparison.json")
```
````
`````

**Example output:**

```
Polymer-Protein Contacts Comparison: polymer_stabilization_study
================================================================================
Polymer selection: chainID C
Contact cutoff: 4.5 Å
Note: 1 condition auto-excluded (no polymer atoms): No Polymer

Condition Summary - Coverage (ranked, highest first)
--------------------------------------------------------------------------------
Rank  Condition            Coverage     SEM        N   
--------------------------------------------------------------------------------
1     100% EGMA               88.4%      0.55%  3    
2     100% SBMA               74.9%      0.28%  3    
--------------------------------------------------------------------------------

Condition Summary - Mean Contact Fraction (ranked, highest first)
--------------------------------------------------------------------------------
Rank  Condition            Contact %    SEM        N   
--------------------------------------------------------------------------------
1     100% SBMA               30.2%      0.50%  3    
2     100% EGMA               25.3%      2.64%  3    
--------------------------------------------------------------------------------

Residence Time by Polymer Type (frames)
--------------------------------------------------------------------------------
Condition            EGM          SBM
--------------------------------------------------------------------------------
100% SBMA             --    10.0±0.2 
100% EGMA       7.1±0.4           --
--------------------------------------------------------------------------------

Results saved: results/contacts_comparison_polymer_stabilization_study.json
```

**Interpretation:**
- EGMA has broader coverage (88.4% of residues contacted)
- SBMA has higher contact fraction (30.2%) and longer residence times (10.0 frames)
- This suggests SBMA forms fewer but more persistent contacts

---

## Step 5: Generate Publication Figures

### 5.1 Generate All Plots

`````{tab-set}
````{tab-item} CLI
```bash
# RMSF comparison plots
polyzymd compare plot results/rmsf_comparison_polymer_stabilization_study.json \
    --dpi 300 --format png

# Triad comparison plots  
polyzymd compare plot results/triad_comparison_polymer_stabilization_study.json \
    --dpi 300 --format png

# For vector graphics (publications)
polyzymd compare plot results/rmsf_comparison_polymer_stabilization_study.json \
    --format pdf
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

# Load RMSF results
rmsf_result = ComparisonResult.load(
    "results/rmsf_comparison_polymer_stabilization_study.json"
)

# Generate individual plots
plot_rmsf_comparison(rmsf_result, save_path="figures/rmsf_bars.png", dpi=300)
plot_percent_change(rmsf_result, save_path="figures/rmsf_pct_change.png", dpi=300)
plot_effect_sizes(rmsf_result, save_path="figures/rmsf_effects.png", dpi=300)

# Generate combined summary panel (ideal for presentations)
plot_summary_panel(
    rmsf_result,
    title="Effect of Polymer Coating on Enzyme Flexibility",
    save_path="figures/rmsf_summary.png",
    dpi=300,
)
```
````
`````

### 5.2 Generated Files

After running plot commands:

```
polymer_stabilization_study/
├── comparison.yaml
├── results/
│   ├── rmsf_comparison_polymer_stabilization_study.json
│   ├── triad_comparison_polymer_stabilization_study.json
│   └── contacts_comparison_polymer_stabilization_study.json
└── figures/
    ├── rmsf_comparison.png       # Bar chart by condition
    ├── percent_change.png        # % change vs control
    ├── effect_sizes.png          # Forest plot of Cohen's d
    └── summary_panel.png         # Combined 3-panel figure
```

### 5.3 Color Coding in Plots

All plots use consistent color coding:

| Color | Meaning |
|-------|---------|
| **Green** | Significant improvement (p < 0.05) |
| **Blue** | Large effect (Cohen's d > 0.8) but not statistically significant |
| **Gray** | Control condition, or negligible effect |
| **Red** | Worse than control |

---

## Step 6: Export Results for Publications

### 6.1 Markdown Tables

Generate publication-ready tables:

```bash
# RMSF table
polyzymd compare rmsf --format markdown -o results/rmsf_table.md

# Triad table
polyzymd compare triad --format markdown -o results/triad_table.md
```

### 6.2 JSON for Further Analysis

All results are saved as JSON for custom analysis:

```python
import json

with open("results/rmsf_comparison_polymer_stabilization_study.json") as f:
    data = json.load(f)

# Access raw data
for cond in data["conditions"]:
    print(f"{cond['label']}: {cond['mean_rmsf']:.3f} ± {cond['sem_rmsf']:.3f} Å")
    print(f"  Replicate values: {cond['replicate_values']}")
```

---

## Quick Reference: All Commands

### Individual Analysis Commands

| Command | Description |
|---------|-------------|
| `polyzymd analyze init` | Create analysis.yaml template |
| `polyzymd analyze run` | Run all enabled analyses |
| `polyzymd analyze rmsf` | Run RMSF analysis only |
| `polyzymd analyze triad` | Run catalytic triad analysis only |
| `polyzymd analyze distances` | Run distance analysis only |
| `polyzymd analyze contacts` | Run contacts analysis only |

### Comparison Commands

| Command | Description |
|---------|-------------|
| `polyzymd compare init NAME` | Create comparison project |
| `polyzymd compare validate` | Check configuration |
| `polyzymd compare rmsf` | Compare RMSF across conditions |
| `polyzymd compare triad` | Compare triad geometry |
| `polyzymd compare contacts` | Compare contact statistics |
| `polyzymd compare plot FILE` | Generate publication figures |
| `polyzymd compare show FILE` | Display saved results |

### Common Options

| Option | Description |
|--------|-------------|
| `--eq-time 10ns` | Override equilibration time |
| `--format table/markdown/json` | Output format |
| `-o FILE` | Save output to file |
| `--recompute` | Force recomputation (ignore cache) |
| `-v, --verbose` | Show detailed logging |

---

## Troubleshooting

### "Config not found" Error

**Cause:** Path in comparison.yaml is incorrect relative to the comparison directory.

**Fix:** Use relative paths from where comparison.yaml is located:
```yaml
config: "../noPoly_enzyme_DMSO/config.yaml"  # One directory up
```

### "No polymer atoms found"

**Cause:** The control condition has no polymer (expected).

**This is normal.** Contacts analysis is automatically skipped for conditions
without polymer atoms. The comparison will exclude these conditions.

### High p-values Despite Large Effects

**Cause:** Small sample size (N=2-3 replicates).

**This is expected.** With few replicates, only very large effects achieve
statistical significance. Report both p-values AND effect sizes (Cohen's d).
A large effect (d > 0.8) suggests a meaningful difference worth investigating
with more replicates.

### RMSF Values Don't Match Manual Calculation

**Cause:** Different equilibration time or reference mode.

**Fix:** Ensure `--eq-time` and reference settings match across analyses.
Check the `rmsf:` section in both analysis.yaml and comparison.yaml.

---

## Next Steps

Now that you've completed the analysis workflow:

1. **Deep dive into specific analyses:**
   - [RMSF Best Practices](analysis_rmsf_best_practices.md) — Reference modes, uncertainty
   - [Triad Best Practices](analysis_triad_best_practices.md) — H-bond definitions
   - [Contacts Cookbook](analysis_contacts_cookbook.md) — Advanced queries

2. **Understand the statistics:**
   - [Comparing Conditions Guide](analysis_compare_conditions.md) — Full statistical details

3. **Customize selections:**
   - [Reference Selection Guide](analysis_reference_selection.md) — MDAnalysis syntax

---

## See Also

- [RMSF Quick Start](analysis_rmsf_quickstart.md)
- [Catalytic Triad Quick Start](analysis_triad_quickstart.md)
- [Distance Analysis Quick Start](analysis_distances_quickstart.md)
- [Contacts Quick Start](analysis_contacts_quickstart.md)
- [Comparing Conditions](analysis_compare_conditions.md)
- [CLI Reference](cli_reference.md)
