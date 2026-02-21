# Complete Analysis Workflow: From Finished Simulations to Publication

Your simulations just finished. Now what?

This guide walks you through the **complete PolyzyMD analysis workflow** — from
raw trajectories to publication-ready figures comparing multiple conditions.
By the end, you will have analyzed protein flexibility (RMSF), active site
geometry (catalytic triad), inter-atomic distances, polymer-protein contacts,
binding preference enrichment, and chaperone-like exposure dynamics across all
your simulation conditions.

```{note}
**Prerequisites:**
- Completed production simulations with trajectory files (`.dcd`, `.xtc`)
- A `config.yaml` file for each simulation condition
- `solvated_system.pdb` topology in each run directory (created during build)

**Time estimate:** 30–90 minutes for a typical 3-condition study (longer if
computing SASA for exposure dynamics)
```

## Overview: The Analysis Pipeline

```
┌──────────────────────────────────────────────────────────────────────────┐
│                         YOUR SIMULATION DATA                              │
│  condition_A/config.yaml   condition_B/config.yaml   condition_C/...    │
└──────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  STEP 1: Create analysis.yaml for each condition                         │
│  Define: replicates, equilibration time, analysis parameters             │
└──────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  STEP 2: Run individual analyses                                          │
│  polyzymd analyze run  →  RMSF, triad, distances, contacts per condition │
└──────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  STEP 3: Create comparison.yaml to compare conditions                    │
│  Define: conditions, control, comparison parameters                       │
└──────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  STEP 4: Run statistical comparisons                                      │
│  polyzymd compare rmsf         →  protein flexibility                    │
│  polyzymd compare triad        →  active site geometry                   │
│  polyzymd compare contacts     →  polymer coverage & contact intensity   │
│  polyzymd compare contacts*    →  binding preference enrichment          │
│  polyzymd compare exposure     →  chaperone-like activity (new)          │
└──────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  STEP 5: Generate publication figures                                     │
│  polyzymd compare plot-all  →  bar charts, KDE panels, summary plots     │
└──────────────────────────────────────────────────────────────────────────┘
```

*Binding preference is enabled inside the `contacts` section of `comparison.yaml`.

---

## Analysis Types at a Glance

| Analysis | Command | What it measures | Guide |
|---|---|---|---|
| RMSF | `compare rmsf` | Per-residue flexibility | [Quick Start](analysis_rmsf_quickstart.md) |
| Catalytic Triad | `compare triad` | Active site H-bond geometry | [Quick Start](analysis_triad_quickstart.md) |
| Distances | `analyze distances` | Specific atom-pair distances | [Quick Start](analysis_distances_quickstart.md) |
| Contacts | `compare contacts` | Polymer coverage & contact fraction | [Quick Start](analysis_contacts_quickstart.md) |
| Binding Preference | `compare contacts`* | Amino acid class enrichment | [Full Guide](analysis_binding_preference.md) |
| Exposure Dynamics | `compare exposure` | Chaperone-like polymer activity | [Full Guide](analysis_exposure_dynamics.md) |

*Binding preference is part of the contacts analysis; enable `compute_binding_preference: true` in your config.

---

## Step 1: Configure Your Analyses

### 1.1 Understand Your Project Structure

A typical PolyzyMD simulation campaign looks like this:

```text
my_enzyme_study/
├── noPoly_enzyme_DMSO/           # Control: no polymer
│   ├── config.yaml
│   └── scratch/ → /scratch/user/project/noPoly/  # Symlink*
│       ├── solvated_system.pdb
│       ├── rep_1/production_1.dcd
│       ├── rep_2/production_1.dcd
│       └── rep_3/production_1.dcd
│
├── SBMA_100_enzyme_DMSO/         # Condition: 100% SBMA polymer
│   ├── config.yaml
│   └── scratch/ → ...            # Symlink*
│
├── EGMA_100_enzyme_DMSO/         # Condition: 100% EGMA polymer
│   ├── config.yaml
│   └── scratch/ → ...            # Symlink*
│
└── analysis/                      # We'll create this
    └── comparison.yaml
```

> **\*Note:** The `scratch/` directories are typically **symlinks** to a separate
> HPC scratch partition where large trajectory files are stored. The `config.yaml`
> file contains the actual path to the scratch directory, so PolyzyMD resolves
> these paths automatically during analysis.

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

This creates an `analysis.yaml` with default values. Customize it:

```yaml
# analysis.yaml - No Polymer Control
# ===================================

# Which replicates to analyze
replicates: [1, 2, 3]

# Shared parameters
defaults:
  equilibration_time: "10ns"    # Skip first 10 ns for equilibration

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
      selection_a: "protein and resid 77 and name OG"
      selection_b: "protein and resid 156 and name NE2"
    # His156 ND1 ↔ Asp133 OD1/OD2 midpoint (histidine-aspartate)
    - label: "His156-Asp133"
      selection_a: "protein and resid 156 and name ND1"
      selection_b: "midpoint(protein and resid 133 and name OD1 OD2)"

# ─────────────────────────────────────────────────────────────────────────
# Distance Analysis: Specific atom-pair distances
# ─────────────────────────────────────────────────────────────────────────
distances:
  enabled: true
  pairs:
    # Substrate to catalytic serine (reaction distance)
    - label: "Substrate-Ser77"
      selection_a: "resname SUB and name C1"      # Replace with your substrate
      selection_b: "protein and resid 77 and name OG"
    # Histidine to substrate (proton shuttle)
    - label: "His156-Substrate"
      selection_a: "protein and resid 156 and name NE2"
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

Always use `"protein and resid X"` in selections — without `protein and`, your
selection may accidentally match atoms from polymer or solvent chains.
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

Enabled analyses: rmsf, catalytic_triad, distances, contacts

[1/4] RMSF Analysis
  Replicates: 1, 2, 3
  Selection: protein and name CA
  Reference mode: average
  Processing replicate 1... done (181 residues, mean RMSF: 0.72 Å)
  Processing replicate 2... done (181 residues, mean RMSF: 0.69 Å)
  Processing replicate 3... done (181 residues, mean RMSF: 0.70 Å)
  Results saved: analysis/rmsf/

[2/4] Catalytic Triad Analysis
  Triad: Ser-His-Asp
  Pairs: Ser77-His156, His156-Asp133
  Threshold: 3.5 Å
  Processing replicate 1... done (72.1% simultaneous contact)
  Processing replicate 2... done (74.3% simultaneous contact)
  Processing replicate 3... done (71.8% simultaneous contact)
  Results saved: analysis/triad/

[3/4] Distance Analysis
  Pairs: Substrate-Ser77, His156-Substrate
  Results saved: analysis/distances/

[4/4] Contacts Analysis
  Note: No polymer atoms found (chainID C). Skipping.

Analysis complete! Results in: analysis/
```

### 2.2 Verify Results

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
├── distances/
│   ├── distances_rep1.json
│   ├── distances_rep2.json
│   └── distances_rep3.json
└── contacts/           # Empty for no-polymer control; populated for polymer conditions
    ├── contacts_rep1.json
    ├── contacts_rep2.json
    └── contacts_rep3.json
```

```{note}
The `contacts/` directory will be empty for the no-polymer control — that is
expected. The comparison step automatically excludes conditions without polymer
when running `compare contacts` and `compare exposure`.

Exposure dynamics analysis reads from the cached `contacts_repN.json` files
produced here. **Run contacts first before running `compare exposure`.**
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
polyzymd analyze distances -c config.yaml -r 1 --eq-time 10ns
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

### 3.1 Create a Comparison Project

```bash
# Go to your study root directory
cd my_enzyme_study/

# Create comparison project
polyzymd compare init -n polymer_stabilization_study
```

This creates:

```
polymer_stabilization_study/
├── comparison.yaml    # Template to edit
├── results/           # Comparison outputs (auto-populated)
├── figures/           # Publication figures (auto-populated)
└── structures/        # Place your enzyme PDB here (for binding preference)
    └── README.md
```

### 3.2 Configure comparison.yaml

Edit `comparison.yaml` to define your conditions and analyses. This single file
drives all comparison analyses:

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
# Analysis Settings (WHAT to analyze)
# ─────────────────────────────────────────────────────────────────────────
analysis_settings:

  rmsf:
    selection: "protein and name CA"
    reference_mode: "average"

  catalytic_triad:
    name: "Ser-His-Asp"
    threshold: 3.5
    pairs:
      - label: "Ser77-His156"
        selection_a: "protein and resid 77 and name OG"
        selection_b: "protein and resid 156 and name NE2"
      - label: "His156-Asp133"
        selection_a: "protein and resid 156 and name ND1"
        selection_b: "midpoint(protein and resid 133 and name OD1 OD2)"

  contacts:
    polymer_selection: "chainID C"
    protein_selection: "protein"
    cutoff: 4.5
    grouping: "aa_class"
    compute_residence_times: true

    # Binding preference: which amino acid classes does each polymer prefer?
    # Requires your enzyme's static PDB in the structures/ directory.
    compute_binding_preference: true
    surface_exposure_threshold: 0.2     # 20% relative SASA = surface exposed
    enzyme_pdb_for_sasa: "structures/enzyme.pdb"
    include_default_aa_groups: true     # aromatic, polar, charged, nonpolar

  # Exposure dynamics: chaperone-like polymer activity
  # Reads from cached contacts_repN.json — run contacts first.
  exposure:
    exposure_threshold: 0.20          # Fraction SASA defining "exposed"
    transient_lower: 0.20             # Residue must dip below this threshold
    transient_upper: 0.80             # ...and rise above this threshold
    min_event_length: 1               # Minimum frames to count as an event
    protein_chain: "A"
    protein_selection: "protein"
    polymer_selection: "chainID C"
    probe_radius_nm: 0.14
    n_sphere_points: 960

# ─────────────────────────────────────────────────────────────────────────
# Comparison Settings (HOW to compare — statistical parameters)
# ─────────────────────────────────────────────────────────────────────────
comparison_settings:
  rmsf: {}

  catalytic_triad: {}

  contacts:
    fdr_alpha: 0.05
    min_effect_size: 0.5

  exposure: {}
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
  Control: No Polymer

  Conditions: 3
    ✓ No Polymer: 3 replicates (config found)
    ✓ 100% SBMA: 3 replicates (config found)
    ✓ 100% EGMA: 3 replicates (config found)

  Analyses configured:
    ✓ rmsf
    ✓ catalytic_triad: Ser-His-Asp (2 pairs)
    ✓ contacts: chainID C ↔ protein, cutoff 4.5 Å
    ✓ exposure: threshold 0.20
```

---

## Step 4: Run Statistical Comparisons

### 4.1 Compare RMSF (Protein Flexibility)

Compares mean per-residue RMSF across conditions. Lower RMSF = more rigid protein.

```bash
polyzymd compare rmsf
```

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
```

**Interpretation:** 100% SBMA significantly reduces RMSF by 22.9% (p < 0.05).
100% EGMA shows a large effect (Cohen's d = 1.27) but does not reach significance
with N=3 replicates.

For detailed guidance: [RMSF Quick Start](analysis_rmsf_quickstart.md) · [RMSF Best Practices](analysis_rmsf_best_practices.md)

---

### 4.2 Compare Catalytic Triad Geometry

Compares the fraction of time all active-site H-bond distances are simultaneously
within threshold. Higher fraction = better-maintained active site geometry.

```bash
polyzymd compare triad
```

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
```

**Interpretation:** SBMA significantly improves active site geometry (+21.1%,
p < 0.05). EGMA slightly reduces it, but the difference is not significant.

For detailed guidance: [Catalytic Triad Quick Start](analysis_triad_quickstart.md) · [Triad Best Practices](analysis_triad_best_practices.md)

---

### 4.3 Analyze Distances

Distance analysis measures specific atom-pair distances over time. Use this
for individual atom-pair statistics, KDE-based mode estimation, or any
inter-atomic distance that falls outside the catalytic triad framework.

```bash
# Single replicate
polyzymd analyze distances -c config.yaml -r 1 --eq-time 10ns

# Multiple replicates with a contact threshold
polyzymd compare distances --eq-time 10ns
```

```{tip}
**Distances vs. triad:** Use **triad** when you need simultaneous contact
analysis (all pairs below threshold at once). Use **distances** when you
want individual pair statistics, KDE distributions, or custom atom selections.
```

**Example output:**

```
Distance Analysis (Aggregated, replicates 1–3)
  Substrate-Ser77:
    Mean: 4.21 ± 0.18 Å
    Contact fraction (<3.5 Å): 32.4% ± 4.2%
  His156-Substrate:
    Mean: 3.89 ± 0.25 Å
    Contact fraction (<3.5 Å): 48.7% ± 5.1%
```

For detailed guidance: [Distance Analysis Quick Start](analysis_distances_quickstart.md)

---

### 4.4 Compare Polymer-Protein Contacts

Contacts analysis provides three complementary views:

1. **Coverage** — What fraction of protein residues does the polymer contact?
2. **Contact fraction** — How often, on average, is each contacted residue in contact?
3. **Residence times** — How persistent are individual contacts?

```bash
polyzymd compare contacts
```

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
```

**Interpretation:** EGMA contacts more of the protein surface (88.4%) while
SBMA makes fewer but more persistent contacts (higher contact fraction and longer
residence times).

For detailed guidance: [Contacts Quick Start](analysis_contacts_quickstart.md) · [Contacts Cookbook](analysis_contacts_cookbook.md)

---

### 4.5 Analyze Binding Preference

Binding preference analysis answers: *"Does the polymer preferentially contact
certain amino acid classes, controlling for surface availability?"*

Enable `compute_binding_preference: true` in the `contacts` section of your
`comparison.yaml` (shown in Step 3.2 above), then re-run contacts:

```bash
# Copy your enzyme's reference PDB to the structures directory first
cp /path/to/enzyme.pdb structures/

polyzymd compare contacts
```

The binding preference table appears at the end of the contacts output:

```
Binding Preference - Enrichment by Amino Acid Class
--------------------------------------------------------------------------------
Surface exposure threshold: 20% relative SASA

  Polymer: EGM
  Protein Group        100% EGMA
  ----------------------------------------
  aromatic             1.90±0.06 +
  charged_negative     0.45±0.04 -
  charged_positive     0.79±0.10 -
  nonpolar             1.31±0.04 +
  polar                0.71±0.02 -

  Polymer: SBM
  Protein Group        100% SBMA
  ----------------------------------------
  aromatic             1.29±0.03 +
  charged_negative     0.90±0.05 -
  charged_positive     1.02±0.03 +
  nonpolar             0.95±0.03 -
  polar                1.00±0.00 =

  + = enriched (>1.0), - = depleted (<1.0)
```

**Interpretation:** EGMA strongly prefers aromatic (1.9×) and nonpolar (1.3×)
residues, consistent with its hydrophobic character. SBMA shows more balanced
binding. A value of 1.9 means the polymer contacts aromatic residues 1.9× more
often than expected given surface availability.

For the full formula derivation and interpretation guide:
[Binding Preference Analysis](analysis_binding_preference.md)

---

### 4.6 Compare Exposure Dynamics

Exposure dynamics analysis answers: *"Does the polymer act as a molecular
chaperone — co-localizing with transiently exposed protein residues during
their exposure episodes?"*

This is a two-metric analysis:

- **Chaperone fraction** — Of all exposure episodes for a transiently exposed
  residue, what fraction have polymer contact during the episode?
- **Transient residue fraction** — What fraction of protein residues undergo
  transient exposure events?

```{note}
Exposure dynamics reads from the cached `contacts_repN.json` files produced
by contacts analysis. Always run contacts first.
```

```bash
polyzymd compare exposure
```

**Optional flags:**

```bash
# Force recompute SASA from trajectory (slow; only needed if trajectory changed)
polyzymd compare exposure --recompute-sasa

# Force recompute exposure classification (fast)
polyzymd compare exposure --recompute-exposure

# Override exposure threshold
polyzymd compare exposure --exposure-threshold 0.25

# Override polymer residue names for enrichment
polyzymd compare exposure --polymer-resnames "SBM,EGM"
```

**Example output:**

```
Exposure Dynamics Comparison: polymer_stabilization_study
================================================================================
Metric: chaperone_fraction
Equilibration: 10ns
Note: 1 condition auto-excluded (no polymer): No Polymer

Condition Summary - Chaperone Fraction (ranked, highest first)
--------------------------------------------------------------------------------
Rank  Condition            Chaperone %    SEM        N
--------------------------------------------------------------------------------
1     100% SBMA               42.3%      3.1%   3
2     100% EGMA               28.7%      4.8%   3
--------------------------------------------------------------------------------

Condition Summary - Transient Residue Fraction (ranked, highest first)
--------------------------------------------------------------------------------
Rank  Condition            Transient %    SEM        N
--------------------------------------------------------------------------------
1     100% EGMA               18.2%      1.4%   3
2     100% SBMA               14.6%      1.9%   3
--------------------------------------------------------------------------------

Chaperone Event Counts (mean over replicates)
--------------------------------------------------------------------------------
Condition            Transient Res   Chaperone Events   Unassisted Events
--------------------------------------------------------------------------------
100% SBMA                     26.5               31.2               42.1
100% EGMA                     32.9               23.8               59.1
--------------------------------------------------------------------------------

Pairwise Statistical Comparisons (chaperone_fraction)
---------------------------------------------------------------------------
Comparison                     % Change   p-value      Cohen's d  Effect
---------------------------------------------------------------------------
100% SBMA vs 100% EGMA         +47.4%     0.0412*      2.14       large
---------------------------------------------------------------------------
* p < 0.05
```

**Interpretation:** SBMA has a significantly higher chaperone fraction (42.3%
vs 28.7%) — when protein residues transiently expose, SBMA is more likely to be
present during that exposure window. EGMA contacts a larger fraction of transiently
exposed residues but provides less persistent co-localization during individual
episodes.

For the full derivation, worked examples, and limitations:
[Exposure Dynamics Analysis](analysis_exposure_dynamics.md)

---

## Step 5: Generate Publication Figures

### 5.1 Generate All Plots

```bash
cd polymer_stabilization_study/

# Generate all configured plots
polyzymd compare plot-all

# Generate plots for a specific analysis only
polyzymd compare plot-all -a rmsf
polyzymd compare plot-all -a catalytic_triad

# List available plot types
polyzymd compare plot-all --list-available
```

### 5.2 Generated File Structure

After running comparisons and plots:

```
polymer_stabilization_study/
├── comparison.yaml
├── structures/
│   └── enzyme.pdb
├── results/
│   ├── rmsf_comparison_polymer_stabilization_study.json
│   ├── triad_comparison_polymer_stabilization_study.json
│   ├── contacts_comparison_polymer_stabilization_study.json
│   └── exposure_comparison_polymer_stabilization_study.json
└── figures/
    ├── rmsf_comparison.png
    ├── rmsf_profile.png
    ├── triad_kde_panel.png
    ├── triad_threshold_bars.png
    ├── distance_kde.png
    └── distance_threshold_bars.png
```

---

## Step 6: Export for Publication

### 6.1 Markdown Tables

Generate formatted tables for manuscripts or supplementary material:

```bash
polyzymd compare rmsf --format markdown -o results/rmsf_table.md
polyzymd compare triad --format markdown -o results/triad_table.md
polyzymd compare contacts --format markdown -o results/contacts_table.md
polyzymd compare exposure --format markdown -o results/exposure_table.md
```

### 6.2 JSON for Custom Analysis

All results are saved as JSON for downstream processing:

```python
import json

with open("results/rmsf_comparison_polymer_stabilization_study.json") as f:
    data = json.load(f)

for cond in data["conditions"]:
    print(f"{cond['label']}: {cond['mean_rmsf']:.3f} ± {cond['sem_rmsf']:.3f} Å")
```

---

## Quick Reference: All Commands

### Individual Analysis Commands

| Command | Description |
|---------|-------------|
| `polyzymd analyze init` | Create `analysis.yaml` template |
| `polyzymd analyze run` | Run all enabled analyses |
| `polyzymd analyze rmsf` | Run RMSF analysis only |
| `polyzymd analyze triad` | Run catalytic triad analysis only |
| `polyzymd analyze distances` | Run distance analysis only |
| `polyzymd analyze contacts` | Run contacts analysis only |

### Comparison Commands

| Command | Description |
|---------|-------------|
| `polyzymd compare init -n NAME` | Create comparison project |
| `polyzymd compare validate` | Check configuration |
| `polyzymd compare rmsf` | Compare RMSF across conditions |
| `polyzymd compare triad` | Compare active site geometry |
| `polyzymd compare contacts` | Compare contacts + binding preference |
| `polyzymd compare exposure` | Compare chaperone-like activity |
| `polyzymd compare run TYPE` | Generic runner (any registered type) |
| `polyzymd compare plot-all` | Generate all configured plots |

### Common Options

| Option | Description |
|--------|-------------|
| `--eq-time 10ns` | Override equilibration time |
| `--format table/markdown/json` | Output format |
| `-o FILE` | Save output to file |
| `--recompute` | Force recomputation (ignore cache) |
| `-q, --quiet` | Suppress INFO messages |
| `--debug` | Enable DEBUG logging |

### Exposure-Specific Options

| Option | Description |
|--------|-------------|
| `--recompute-sasa` | Force recompute SASA from trajectory |
| `--recompute-exposure` | Force recompute exposure classification |
| `--exposure-threshold FLOAT` | Override SASA threshold (default: 0.20) |
| `--polymer-resnames LIST` | Override polymer residue names (comma-separated) |

---

## Recommended Analysis Order

When running a full study, follow this order to satisfy dependencies:

```
1. polyzymd analyze run           (per condition — produces contacts_repN.json)
      ↓
2. polyzymd compare rmsf          (no dependencies)
   polyzymd compare triad         (no dependencies)
   polyzymd compare distances     (no dependencies)
   polyzymd compare contacts      (reads contacts_repN.json)
      ↓ (contacts must be done first)
3. polyzymd compare exposure      (reads contacts_repN.json + computes SASA)
```

Steps 2 are independent of each other and can be run in any order.
Step 3 requires contacts to be complete for all conditions.

---

## Troubleshooting

### "Config not found" Error

**Cause:** Path in `comparison.yaml` is incorrect relative to the comparison directory.

**Fix:** Use relative paths from where `comparison.yaml` is located:
```yaml
config: "../noPoly_enzyme_DMSO/config.yaml"  # One directory up
```

### "No polymer atoms found"

**Cause:** The control condition has no polymer. This is expected.

Contacts analysis is automatically skipped for conditions without polymer atoms.
Both `compare contacts` and `compare exposure` exclude these conditions automatically —
you will see a yellow notice in the output.

### "Contacts not found for rep N"

**Cause:** Running `compare exposure` before `compare contacts`.

**Fix:** Run contacts analysis for all conditions first:
```bash
for condition in SBMA_100_enzyme_DMSO EGMA_100_enzyme_DMSO; do
    cd $condition && polyzymd analyze run && cd ..
done
cd polymer_stabilization_study
polyzymd compare contacts      # Produces contacts_repN.json cache
polyzymd compare exposure      # Now has what it needs
```

### High p-values Despite Large Effects

**Cause:** Small sample size (N = 2–3 replicates).

With few replicates, only very large effects achieve statistical significance.
Report both p-values **and** effect sizes (Cohen's d). A large effect (d > 0.8)
indicates a meaningful difference worth investigating further.
See [Statistics Best Practices](analysis_statistics_best_practices.md).

### RMSF Values Differ from Manual Calculation

**Cause:** Different equilibration time or reference mode.

**Fix:** Ensure `--eq-time` and reference settings match across `analysis.yaml`
and `comparison.yaml`. Check [RMSF Best Practices](analysis_rmsf_best_practices.md).

---

## Next Steps

1. **Deep dive into specific analyses:**
   - [RMSF Best Practices](analysis_rmsf_best_practices.md) — Reference modes, uncertainty
   - [Triad Best Practices](analysis_triad_best_practices.md) — H-bond definitions
   - [Contacts Cookbook](analysis_contacts_cookbook.md) — Advanced contact queries
   - [Binding Preference Analysis](analysis_binding_preference.md) — Amino acid class enrichment
   - [Exposure Dynamics Analysis](analysis_exposure_dynamics.md) — Chaperone fraction, transient exposure

2. **Understand the statistics:**
   - [Statistics Best Practices](analysis_statistics_best_practices.md) — Autocorrelation, SEM correction
   - [Comparing Conditions Guide](analysis_compare_conditions.md) — Full statistical details

3. **Customize selections:**
   - [Reference Selection Guide](analysis_reference_selection.md) — MDAnalysis syntax

4. **Extend the framework:**
   - [Extending Comparators](extending_comparators.md) — Add new comparison types
   - [Extending Plotters](extending_plotters.md) — Add new plot types

---

## See Also

- [RMSF Quick Start](analysis_rmsf_quickstart.md)
- [Catalytic Triad Quick Start](analysis_triad_quickstart.md)
- [Distance Analysis Quick Start](analysis_distances_quickstart.md)
- [Contacts Quick Start](analysis_contacts_quickstart.md)
- [Binding Preference Analysis](analysis_binding_preference.md)
- [Exposure Dynamics Analysis](analysis_exposure_dynamics.md)
- [Statistics Best Practices](analysis_statistics_best_practices.md)
- [Comparing Conditions](analysis_compare_conditions.md)
- [CLI Reference](cli_reference.md)
