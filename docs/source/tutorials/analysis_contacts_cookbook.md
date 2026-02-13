# Contacts Analysis Cookbook

Worked examples for answering common scientific questions about polymer-protein 
interactions using PolyzyMD's contacts analysis module.

```{note}
**Prerequisites:** Complete a contacts analysis first using the 
[Contacts Quick Start](analysis_contacts_quickstart.md), then return here for 
advanced analysis techniques.
```

## Overview

This cookbook answers questions like:

| Recipe | Scientific Question |
|--------|---------------------|
| [Recipe 1](#recipe-1-do-zwitterionic-polymers-preferentially-contact-charged-residues) | Do zwitterionic polymers preferentially contact charged residues? |
| [Recipe 2](#recipe-2-which-protein-regions-does-my-polymer-bind) | Which protein regions does my polymer bind? |
| [Recipe 3](#recipe-3-comparing-residence-times-by-polymer-type) | Do zwitterionic polymers form more persistent contacts? |
| [Recipe 4](#recipe-4-active-site-contacts-vs-surface-contacts) | Does polymer bind near the active site? |
| [Recipe 5](#recipe-5-complex-queries-with-compositeselector) | Complex multi-criteria selections |

---

## Recipe 1: Do zwitterionic polymers preferentially contact charged residues?

### Scientific Question

You've synthesized copolymers with zwitterionic (SBMA) and PEG-like (EGMA) 
monomers. You want to know if the zwitterionic component shows preferential 
binding to charged protein residues compared to nonpolar residues.

### CLI Approach: Separate Analyses

Run two separate contact analyses with different polymer selections:

```bash
# Analyze SBMA contacts only
polyzymd analyze contacts -c config.yaml -r 1-3 --eq-time 10ns \
    --polymer-selection "segid C and resname SBM" \
    -o analysis/contacts_sbma/

# Analyze EGMA contacts only
polyzymd analyze contacts -c config.yaml -r 1-3 --eq-time 10ns \
    --polymer-selection "segid C and resname EGM" \
    -o analysis/contacts_egma/
```

Then compare the outputs manually or load both in Python.

### Python Approach: Interaction Matrix

The more powerful approach uses `interaction_matrix()` to get a complete 
breakdown in one analysis:

```python
from polyzymd.analysis.contacts.results import ContactResult
from polyzymd.analysis.common.groupings import CustomGrouping

# Load a contacts result (analyzed with all polymer types)
result = ContactResult.load("analysis/contacts/contacts_rep1.json")

# Define custom polymer groupings
polymer_grouping = CustomGrouping.from_groups({
    "zwitterionic": ["SBM", "SBMA", "MPC"],
    "peg_like": ["EGM", "EGMA", "OEGMA"],
})

# Get interaction matrix with custom grouping
matrix = result.interaction_matrix(
    metric="contact_fraction",
    polymer_grouping=polymer_grouping,
)

# Compare zwitterionic vs PEG-like contacts with different AA classes
print("Contact fractions by polymer type and AA class:")
print("-" * 50)
for polymer_group in ["zwitterionic", "peg_like"]:
    print(f"\n{polymer_group}:")
    for aa_class in ["charged_positive", "charged_negative", "aromatic", "nonpolar"]:
        if aa_class in matrix[polymer_group]:
            frac = matrix[polymer_group][aa_class]
            print(f"  {aa_class}: {frac:.1%}")
```

**Example output:**
```
Contact fractions by polymer type and AA class:
--------------------------------------------------

zwitterionic:
  charged_positive: 47.9%
  charged_negative: 52.3%
  aromatic: 45.4%
  nonpolar: 27.2%

peg_like:
  charged_positive: 31.5%
  charged_negative: 38.1%
  aromatic: 37.0%
  nonpolar: 33.0%
```

### Interpreting Results

In this example:
- **Zwitterionic polymers** show higher contact fractions with charged residues 
  (47.9% and 52.3%) compared to nonpolar residues (27.2%)
- **PEG-like polymers** show more uniform contact across AA classes
- The zwitterionic preference for charged residues is consistent with 
  electrostatic complementarity

```{tip}
For statistical comparison across multiple replicates/conditions, use 
`polyzymd compare contacts` with custom polymer selections.
```

---

## Recipe 2: Which protein regions does my polymer bind?

### Scientific Question

You want a breakdown of polymer contacts by amino acid type to understand 
the binding surface characteristics.

### Python Approach: Coverage by Group

```python
from polyzymd.analysis.contacts.results import ContactResult

result = ContactResult.load("analysis/contacts/contacts_rep1.json")

# Get coverage (fraction of residues contacted) by AA class
coverage = result.coverage_by_group()

print("Polymer coverage by amino acid class:")
print("-" * 40)
for group, frac in sorted(coverage.items(), key=lambda x: -x[1]):
    # Visual bar
    bar = "#" * int(frac * 20)
    print(f"{group:20s} {frac:6.1%} {bar}")
```

**Example output:**
```
Polymer coverage by amino acid class:
----------------------------------------
charged_negative     100.0% ####################
aromatic             100.0% ####################
charged_positive     100.0% ####################
polar                 93.5% ##################
nonpolar              86.2% #################
```

### Understanding the Default AA Classification

PolyzyMD uses `ProteinAAClassification` which groups amino acids as:

| Group | Residues |
|-------|----------|
| `aromatic` | TRP, PHE, TYR (optionally HIS) |
| `charged_positive` | LYS, ARG |
| `charged_negative` | ASP, GLU |
| `polar` | SER, THR, ASN, GLN, HIS, CYS |
| `nonpolar` | ALA, VAL, LEU, ILE, MET, PRO, GLY |

### Using Custom Classifications

Define your own residue groupings:

```python
from polyzymd.analysis.common.groupings import CustomGrouping

# Custom grouping for hydrophobicity
hydrophobicity_grouping = CustomGrouping.from_groups({
    "hydrophobic": ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"],
    "hydrophilic": ["SER", "THR", "ASN", "GLN", "LYS", "ARG", "ASP", "GLU", "HIS"],
    "amphipathic": ["TYR", "CYS", "GLY"],
})

# Use with interaction_matrix
matrix = result.interaction_matrix(
    metric="contact_fraction",
    protein_grouping=hydrophobicity_grouping,
)
```

---

## Recipe 3: Comparing residence times by polymer type

### Scientific Question

Do zwitterionic polymers form more persistent (longer-lasting) contacts with 
the protein compared to PEG-like polymers?

### Python Approach: Residence Time Summary

```python
from polyzymd.analysis.contacts.results import ContactResult

result = ContactResult.load("analysis/contacts/contacts_rep1.json")

# Get residence time statistics by polymer type
rt_summary = result.residence_time_summary()

print("Residence time by polymer type:")
print("-" * 50)
for ptype, stats in rt_summary.items():
    mean_frames = stats["mean_frames"]
    max_frames = stats["max_frames"]
    n_events = stats["n_events"]
    print(f"{ptype}:")
    print(f"  Mean: {mean_frames:.1f} frames")
    print(f"  Max:  {max_frames} frames")
    print(f"  Events: {n_events}")
```

**Example output:**
```
Residence time by polymer type:
--------------------------------------------------
SBM:
  Mean: 4.6 frames
  Max:  181 frames
  Events: 2847
EGM:
  Mean: 4.2 frames
  Max:  272 frames
  Events: 3102
```

### Converting to Physical Time

Multiply by your trajectory's output interval:

```python
# If trajectory saves every 100 ps
timestep_ps = 100  # ps per frame

for ptype, stats in rt_summary.items():
    mean_ns = stats["mean_frames"] * timestep_ps / 1000
    max_ns = stats["max_frames"] * timestep_ps / 1000
    print(f"{ptype}: mean={mean_ns:.2f} ns, max={max_ns:.1f} ns")
```

### Residence Time by AA Class

Get residence times broken down by both polymer type AND amino acid class:

```python
matrix = result.interaction_matrix(metric="residence_time")

print("Mean residence time (frames) by polymer type and AA class:")
for ptype in matrix.keys():
    print(f"\n{ptype}:")
    for aa_class, mean_frames in matrix[ptype].items():
        print(f"  {aa_class}: {mean_frames:.1f} frames")
```

---

## Recipe 4: Active site contacts vs surface contacts

### Scientific Question

Does your polymer bind preferentially near the enzyme's active site, or only 
on the exposed surface? This is important for understanding whether polymer 
conjugation might affect catalytic function.

### CLI Approach: Region-Specific Analysis

```bash
# Active site region (example: LipA catalytic triad vicinity)
polyzymd analyze contacts -c config.yaml -r 1-3 --eq-time 10ns \
    --protein-selection "protein and (resid 75-80 or resid 130-140 or resid 153-160)" \
    -o analysis/contacts_active_site/

# Surface residues (everything else)
polyzymd analyze contacts -c config.yaml -r 1-3 --eq-time 10ns \
    --protein-selection "protein and not (resid 75-80 or resid 130-140 or resid 153-160)" \
    -o analysis/contacts_surface/
```

### Python Approach: Proximity-Based Selection

Use `ProteinResiduesNearReference` to select residues within a cutoff of 
reference atoms (e.g., catalytic residues):

```python
from polyzymd.analysis.common.selectors import ProteinResiduesNearReference

# Select protein residues within 8 A of the catalytic triad
active_site_selector = ProteinResiduesNearReference(
    reference_selection="resid 77 133 156 and name CA",
    cutoff=8.0,
)

# Use in analysis (programmatic API)
from polyzymd.analysis.contacts import ContactAnalyzer

analyzer = ContactAnalyzer(
    universe=universe,
    protein_selector=active_site_selector,
    # ... other parameters
)
```

### Comparing Active Site vs Surface Contact Fractions

After running both analyses, compare the results:

```python
from polyzymd.analysis.contacts.results import ContactResult

active_site = ContactResult.load("analysis/contacts_active_site/contacts_rep1.json")
surface = ContactResult.load("analysis/contacts_surface/contacts_rep1.json")

print(f"Active site mean contact: {active_site.mean_contact_fraction():.1%}")
print(f"Surface mean contact:     {surface.mean_contact_fraction():.1%}")

# If active site contact is much lower, polymer avoids the catalytic region
ratio = surface.mean_contact_fraction() / active_site.mean_contact_fraction()
if ratio > 2:
    print("Polymer preferentially binds surface (good for catalysis)")
elif ratio < 0.5:
    print("Polymer preferentially binds active site (may affect catalysis)")
else:
    print("Polymer binds uniformly across protein")
```

---

## Recipe 5: Complex queries with CompositeSelector

### Scientific Question

Find contacts where BOTH conditions are met:
- Polymer is zwitterionic (SBMA/MPC)
- Protein residue is aromatic AND near the active site

This requires combining multiple selection criteria.

### Python Approach: CompositeSelector

```python
from polyzymd.analysis.common.selectors import (
    PolymerResiduesByType,
    ProteinResiduesByGroup,
    ProteinResiduesNearReference,
    CompositeSelector,
)
from polyzymd.analysis.common.groupings import ProteinAAClassification

# 1. Polymer selector: zwitterionic monomers only
polymer_selector = PolymerResiduesByType(
    residue_names=["SBM", "SBMA", "MPC"]
)

# 2. Protein selector: aromatic residues
aromatic_selector = ProteinResiduesByGroup(
    grouping=ProteinAAClassification(),
    groups=["aromatic"],
)

# 3. Protein selector: near active site
near_active_site = ProteinResiduesNearReference(
    reference_selection="resid 77 133 156",
    cutoff=10.0,
)

# 4. Combine protein selectors with AND logic
aromatic_near_active_site = CompositeSelector(
    selectors=[aromatic_selector, near_active_site],
    mode="intersection",  # AND logic
)

# 5. Use in analysis
from polyzymd.analysis.contacts import ContactAnalyzer

analyzer = ContactAnalyzer(
    universe=universe,
    polymer_selector=polymer_selector,
    protein_selector=aromatic_near_active_site,
    cutoff=4.5,
)
result = analyzer.run()
```

### CompositeSelector Modes

| Mode | Behavior |
|------|----------|
| `"union"` | OR logic - residue matches if it passes ANY selector |
| `"intersection"` | AND logic - residue matches if it passes ALL selectors |

### Example: Exclude Specific Regions

Select all protein residues EXCEPT those near the active site:

```python
from polyzymd.analysis.common.selectors import (
    ProteinResidues,
    ProteinResiduesNearReference,
    CompositeSelector,
)

# All protein residues
all_protein = ProteinResidues()

# Residues near active site (to exclude)
near_active = ProteinResiduesNearReference(
    reference_selection="resid 77 133 156",
    cutoff=8.0,
    exclude=True,  # Invert selection
)

# Combine: all protein AND (not near active site)
surface_only = CompositeSelector(
    selectors=[all_protein, near_active],
    mode="intersection",
)
```

---

## API Reference Summary

### ContactResult Methods

| Method | Returns | Use Case |
|--------|---------|----------|
| `interaction_matrix(metric)` | `dict[polymer][aa_group]` | Compare polymer types x AA classes |
| `coverage_by_group()` | `dict[aa_group] -> float` | Which AA classes are contacted |
| `residence_time_summary()` | `dict[polymer_type] -> stats` | Residence time by polymer |
| `mean_contact_fraction()` | `float` | Overall contact fraction |
| `coverage()` | `float` | Fraction of protein residues contacted |

### Selector Classes

| Selector | Purpose |
|----------|---------|
| `PolymerResiduesByType` | Filter polymer by monomer type (resname) |
| `ProteinResiduesByGroup` | Filter protein by AA classification |
| `ProteinResiduesNearReference` | Proximity to reference atoms |
| `CompositeSelector` | Combine selectors with AND/OR logic |
| `MDAnalysisSelector` | Arbitrary MDAnalysis selection string |

### Grouping Classes

| Grouping | Purpose |
|----------|---------|
| `ProteinAAClassification` | Standard AA groups (aromatic, charged, polar, nonpolar) |
| `CustomGrouping.from_groups()` | Define your own classification |

### Imports

```python
# Results
from polyzymd.analysis.contacts.results import ContactResult

# Selectors
from polyzymd.analysis.common.selectors import (
    PolymerResiduesByType,
    ProteinResiduesByGroup,
    ProteinResiduesNearReference,
    CompositeSelector,
    MDAnalysisSelector,
)

# Groupings
from polyzymd.analysis.common.groupings import (
    ProteinAAClassification,
    CustomGrouping,
)
```

---

## See Also

- [Contacts Quick Start](analysis_contacts_quickstart.md) - basic usage and CLI reference
- [Comparing Conditions](analysis_compare_conditions.md) - statistical comparison across simulations
- [RMSF Quick Start](analysis_rmsf_quickstart.md) - complementary flexibility analysis
