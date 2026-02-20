# Binding Preference Analysis

Understand **which amino acid classes your polymer preferentially binds** using
enrichment analysis. This module answers the fundamental question: *"Does my
polymer type show preferential affinity for certain types of protein residues?"*

```{note}
Binding preference analysis builds on [contacts analysis](analysis_contacts_quickstart.md).
Run contacts analysis first, then enable binding preference in your comparison
configuration.
```

## Scientific Motivation

When designing enzyme-polymer conjugates, understanding polymer binding
preferences is crucial for several reasons:

1. **Active site protection**: Does the polymer avoid or shield the catalytic
   residues?
2. **Electrostatic complementarity**: Do zwitterionic polymers preferentially
   bind charged surface residues?
3. **Hydrophobic interactions**: Do PEG-like polymers show affinity for
   nonpolar patches?
4. **Rational design**: Which polymer chemistry best matches your enzyme's
   surface properties?

Raw contact counts are misleading because protein surfaces have unequal
distributions of amino acid types. A polymer might contact 50% of aromatic
residues and 50% of polar residues, but if the surface has 5 aromatic residues
and 50 polar residues, the polymer clearly *prefers* aromatic residues.

**Binding preference analysis normalizes by surface availability** to reveal
true preferences.

## The Enrichment Formula

For each (polymer type, protein group) pair, we compute a **zero-centered
enrichment** value:

```{math}
\text{Enrichment} = \frac{\text{Contact Share}}{\text{Expected Share}} - 1
```

where:

```{math}
\text{Contact Share} = \frac{\sum_{\text{residues in group}} \text{contact frames}}{\sum_{\text{all residues}} \text{contact frames}}
```

and:

```{math}
\text{Expected Share} = \frac{\text{exposed residues in group}}{\text{total exposed residues}}
```

This **surface-availability normalization** asks: *"Given how much of the
protein surface is aromatic/charged/etc., does this polymer type contact that
surface proportionally?"*

### Interpretation (Zero-Centered Scale)

| Enrichment | Meaning |
|------------|---------|
| **> 0** | **Preferential binding** — polymer contacts this group more than expected |
| **= 0** | **Neutral** — contact frequency matches random chance |
| **< 0** | **Avoidance** — polymer contacts this group less than expected |
| **= -1** | **Complete avoidance** — no contacts observed |

**Examples:**
- `+0.45` means "45% more contacts than expected"
- `-0.30` means "30% fewer contacts than expected"
- `+1.00` means "2× as many contacts as expected" (100% more)

### Why Surface-Availability Normalization?

The correct baseline for "expected contacts" is **protein surface availability**,
not polymer composition. Here's why:

**The question we're asking**: *"Does this polymer type preferentially bind
aromatic residues?"*

**The correct comparison**: If aromatic residues comprise 10% of the protein's
surface, and this polymer contacts aromatics 15% of the time, then enrichment
= (0.15 / 0.10) - 1 = +0.50 (50% preference).

**What NOT to compare**: The fraction of the polymer that is SBMA vs EGMA is
irrelevant to this question — it tells us about polymer composition, not about
which protein groups the polymer prefers.

```{note}
Polymer composition (residue counts, heavy atom counts) is stored as **metadata**
in the results for secondary analysis, but is not used in the enrichment formula.
```

### Why Contact Frames, Not Binary Counts?

The formula uses **contact frames** (duration-weighted) rather than binary
"contacted yes/no" counts. This ensures that:

- A residue contacted for 60% of the simulation contributes 60× more than one
  contacted for 1 frame
- Transient, non-specific contacts don't dominate the signal
- Stable, functionally relevant interactions are emphasized

### Surface Exposure Filtering

Only **surface-exposed residues** are considered in the enrichment calculation.
Buried residues are excluded because they're physically inaccessible to polymer.

Surface exposure is determined using **Solvent Accessible Surface Area (SASA)**:
- Residues with relative SASA > threshold (default 20%) are considered exposed
- This is computed from a static enzyme PDB structure
- The threshold can be adjusted in configuration

## Quick Start

### Enable Binding Preference in comparison.yaml

`````{tab-set}
````{tab-item} YAML (Recommended)
Add `compute_binding_preference: true` to your contacts analysis settings:

```yaml
# comparison.yaml
name: "Polymer Binding Preference Study"
control: "No Polymer"

structures:
  enzyme_pdb: "structures/enzyme.pdb"  # For SASA calculation

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% SBMA"
    config: "../sbma_100/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% EGMA"
    config: "../egma_100/config.yaml"
    replicates: [1, 2, 3]

analysis_settings:
  contacts:
    name: "polymer_contacts"
    polymer_selection: "resname SBM EGM"
    protein_selection: "protein"
    cutoff: 4.5
    
    # Binding preference settings
    compute_binding_preference: true
    surface_exposure_threshold: 0.2      # 20% relative SASA
    enzyme_pdb_for_sasa: "structures/enzyme.pdb"
    include_default_aa_groups: true
```

Then run the comparison:

```bash
polyzymd compare contacts -f comparison.yaml
```
````

````{tab-item} CLI
The CLI automatically uses binding preference settings from `comparison.yaml`.
You cannot enable binding preference via CLI flags alone—it must be configured
in YAML.

```bash
# Run comparison (binding preference enabled in YAML)
polyzymd compare contacts -f comparison.yaml

# View results in markdown format
polyzymd compare contacts -f comparison.yaml --format markdown
```
````

````{tab-item} Python
```python
from pathlib import Path
from polyzymd.analysis.contacts import ContactResult
from polyzymd.analysis.contacts.binding_preference import compute_binding_preference
from polyzymd.analysis.contacts.surface_exposure import SurfaceExposureFilter

# Load contact results
contacts = ContactResult.load("analysis/contacts/contacts_rep1.json")

# Compute surface exposure
exposure_filter = SurfaceExposureFilter(threshold=0.2)
surface_exposure = exposure_filter.calculate("structures/enzyme.pdb")

# Define protein groups (resolved to residue IDs)
# You can use the default AA class groups or define custom ones
from polyzymd.analysis.common.aa_classification import AA_CLASS_RESIDUES

protein_groups = {
    "aromatic": {12, 45, 67, 89, 102},      # resids of aromatic residues
    "charged_positive": {23, 34, 56, 78},   # LYS, ARG resids
    "charged_negative": {15, 28, 91},       # ASP, GLU resids
    "nonpolar": {5, 10, 20, 30, 40, 50},    # hydrophobic resids
    "polar": {8, 18, 25, 35, 42, 55},       # SER, THR, etc.
}

# Extract polymer composition from trajectory (stored as metadata)
from polyzymd.analysis.contacts.binding_preference import extract_polymer_composition
polymer_composition = extract_polymer_composition(universe)  # universe from trajectory

# Compute binding preference
result = compute_binding_preference(
    contact_result=contacts,
    surface_exposure=surface_exposure,
    protein_groups=protein_groups,
    polymer_composition=polymer_composition,
)

# Access enrichment values (normalized by protein surface availability)
print("Enrichment matrix:")
for polymer_type, groups in result.enrichment_matrix().items():
    print(f"\n{polymer_type}:")
    for group, enrichment in groups.items():
        marker = "+" if enrichment > 0 else "-" if enrichment < 0 else "="
        print(f"  {group}: {enrichment:+.2f} {marker}")

# Get detailed entry for a specific pair
entry = result.get_entry("SBM", "aromatic")
print(f"\nSBM → aromatic:")
print(f"  Contact share: {entry.contact_share:.3f}")
print(f"  Expected share (surface): {entry.expected_share:.3f}")
print(f"  Enrichment: {entry.enrichment:+.2f}")
```
````
`````

### Example Output

When you run `polyzymd compare contacts` with binding preference enabled, you'll
see output like this:

```
Binding Preference - Enrichment by Amino Acid Class
-----------------------------------------------------------------------------------------------
Surface exposure threshold: 20% relative SASA
Enrichment normalized by: protein surface availability

  Polymer: EGM
  Protein Group        0% SBMA / 100%     25% SBMA / 75%     50% SBMA / 50%     
  ----------------------------------------------------------------------------------
  aromatic            +0.90±0.06 +       +0.65±0.16 +       +0.50±0.06 +        
  charged_negative    -0.55±0.04 -       -0.44±0.10 -       -0.39±0.08 -        
  charged_positive    -0.21±0.10 -       -0.11±0.07 -       -0.06±0.04 -        
  nonpolar            +0.31±0.04 +       +0.14±0.01 +       +0.13±0.06 +        
  polar               -0.29±0.02 -       -0.13±0.05 -       -0.12±0.03 -        

  Polymer: SBM
  Protein Group        100% SBMA / 0%     25% SBMA / 75%     50% SBMA / 50%     
  ----------------------------------------------------------------------------------
  aromatic            +0.29±0.03 +       +0.83±0.09 +       +0.54±0.07 +        
  charged_negative    -0.10±0.05 -       -0.24±0.18 -       -0.29±0.11 -        
  charged_positive    +0.02±0.03 +       +0.22±0.10 +       +0.10±0.02 +        
  nonpolar            -0.05±0.03 -       -0.09±0.12 -       -0.12±0.05 -        
  polar               +0.00±0.00 =       -0.17±0.04 -       +0.01±0.06 +        

  + = enriched (>0), - = depleted (<0)
```

### Interpreting These Results

From the example above, we can draw several conclusions:

1. **EGMA (EGM)** shows:
   - Strong preference for **aromatic** residues (+0.90 = 90% more contacts than expected)
   - Strong preference for **nonpolar** residues (+0.31 = 31% more contacts)
   - Avoidance of **charged** residues (-0.55 for negative, -0.21 for positive)
   - This is consistent with EGMA's hydrophobic PEG-like character

2. **SBMA (SBM)** shows:
   - Moderate preference for **aromatic** residues (+0.29)
   - Slight preference for **charged positive** residues (+0.02)
   - Near-neutral behavior for most groups
   - This reflects SBMA's zwitterionic nature (balanced charges)

3. **Composition effects**: As SBMA ratio increases (25% → 50% → 100%):
   - EGMA's aromatic preference decreases (+0.90 → +0.65 → +0.50)
   - This suggests competition for aromatic binding sites

## Configuration Options

### ContactsAnalysisSettings

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `compute_binding_preference` | bool | `false` | Enable binding preference analysis |
| `surface_exposure_threshold` | float | `0.2` | Relative SASA threshold (0.0-1.0) |
| `enzyme_pdb_for_sasa` | str | `null` | Path to enzyme PDB for SASA calculation |
| `include_default_aa_groups` | bool | `true` | Use standard AA class groupings |
| `protein_groups` | dict | `null` | Custom protein groups `{name: [resids]}` |
| `polymer_type_selections` | dict | `null` | Custom polymer type definitions (see below) |

### Polymer Type Selections

By default, polymer types are auto-detected from unique residue names in the
polymer selection. You can explicitly define polymer types with custom
selections:

```yaml
analysis_settings:
  contacts:
    polymer_selection: "chainID C"
    compute_binding_preference: true
    
    # Explicit polymer type definitions
    polymer_type_selections:
      SBMA: "chainID C and resname SBM"
      EGMA: "chainID C and resname EGM"
      OEGMA: "chainID C and resname OEG"
```

This is useful when:
- Residue names don't match your preferred labels
- You want to group multiple residue types together
- You need fine-grained control over what counts as each polymer type

### Custom Protein Groups

You can define custom protein groups to analyze specific regions of interest:

```yaml
analysis_settings:
  contacts:
    compute_binding_preference: true
    surface_exposure_threshold: 0.2
    enzyme_pdb_for_sasa: "structures/enzyme.pdb"
    
    # Include default AA class groups (aromatic, polar, etc.)
    include_default_aa_groups: true
    
    # Add custom groups (can override defaults if names conflict)
    protein_groups:
      active_site: [77, 133, 156]           # Catalytic triad
      lid_domain: [140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150]
      binding_pocket: [50, 51, 52, 80, 81, 82, 83]
```

```{tip}
If a custom group name matches a default AA class name (e.g., `aromatic`),
your custom definition takes precedence.
```

### Surface Exposure Threshold

The `surface_exposure_threshold` controls which residues are considered
"surface accessible":

| Threshold | Effect |
|-----------|--------|
| `0.1` (10%) | More permissive — includes partially buried residues |
| `0.2` (20%) | **Default** — standard definition of "exposed" |
| `0.3` (30%) | More stringent — only highly exposed residues |
| `0.5` (50%) | Very stringent — only fully exposed loops/turns |

```{note}
Lower thresholds include more residues in the analysis, which may dilute
the signal. Higher thresholds focus on the most accessible residues but
may have smaller sample sizes.
```

## Default Amino Acid Classifications

When `include_default_aa_groups: true`, the following standard groupings are used:

| Group | Amino Acids | Chemical Property |
|-------|-------------|-------------------|
| `aromatic` | TRP, PHE, TYR | π-stacking, hydrophobic |
| `charged_positive` | LYS, ARG | Cationic at pH 7 |
| `charged_negative` | ASP, GLU | Anionic at pH 7 |
| `nonpolar` | ALA, VAL, LEU, ILE, MET, PRO, GLY | Hydrophobic/aliphatic |
| `polar` | SER, THR, ASN, GLN, HIS, CYS | H-bond donors/acceptors |

```{note}
Histidine (HIS) is classified as **polar** by default because its protonation
state is pH-dependent. At low pH, it would be positively charged.
```

## Output Files

Binding preference results are saved alongside contact analysis results:

```
project/
└── condition_name/
    └── analysis/
        └── contacts/
            ├── contacts_rep1.json
            ├── contacts_rep2.json
            ├── contacts_rep3.json
            ├── binding_preference_rep1.json          # Per-replicate
            ├── binding_preference_rep2.json
            ├── binding_preference_rep3.json
            └── binding_preference_aggregated_reps1-3.json  # Aggregated
```

### JSON Structure

```json
{
  "entries": [
    {
      "polymer_type": "SBM",
      "protein_group": "aromatic",
      "total_contact_frames": 12543,
      "mean_contact_fraction": 0.45,
      "n_residues_in_group": 20,
      "n_exposed_in_group": 7,
      "n_residues_contacted": 6,
      "contact_share": 0.164,
      "expected_share": 0.086,
      "enrichment": 0.907,
      "polymer_residue_count": 50,
      "total_polymer_residues": 100,
      "polymer_heavy_atom_count": 750,
      "total_polymer_heavy_atoms": 1150
    }
  ],
  "polymer_composition": {
    "residue_counts": {"SBM": 50, "EGM": 50},
    "heavy_atom_counts": {"SBM": 750, "EGM": 400}
  },
  "n_frames": 10000,
  "total_exposed_residues": 81,
  "surface_exposure_threshold": 0.2,
  "schema_version": 3
}
```

```{note}
Schema version 3 uses **surface-availability normalization** for enrichment.
Polymer composition (residue/atom counts) is stored as metadata for secondary
analysis but is not used in the enrichment calculation.
```

## Statistical Treatment

### Aggregation Across Replicates

When multiple replicates are available, binding preference is computed
**per-replicate** and then aggregated:

1. Compute enrichment for each replicate independently
2. Report mean enrichment ± SEM across replicates
3. Store per-replicate values for downstream analysis

This approach:
- Preserves replicate independence for statistical testing
- Provides uncertainty estimates (SEM)
- Enables detection of outlier replicates

### Uncertainty in Enrichment

The uncertainty in enrichment comes from:
- Variation in contact patterns across replicates
- Different initial configurations
- Thermal fluctuations in the simulation

The SEM (Standard Error of the Mean) quantifies confidence in the mean
enrichment across replicates.

## Worked Example: Comparing Polymer Preferences

### Scientific Question

*"Do zwitterionic (SBMA) and PEG-like (EGMA) polymers show different binding
preferences? Which polymer type is better for protecting the enzyme's active
site?"*

### Analysis Workflow

```yaml
# comparison.yaml
name: "SBMA vs EGMA Binding Preferences"
control: "No Polymer"

structures:
  enzyme_pdb: "structures/lipase.pdb"

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% SBMA"
    config: "../sbma_100/config.yaml"
    replicates: [1, 2, 3]

  - label: "100% EGMA"
    config: "../egma_100/config.yaml"
    replicates: [1, 2, 3]

  - label: "50% SBMA / 50% EGMA"
    config: "../copolymer_50_50/config.yaml"
    replicates: [1, 2, 3]

analysis_settings:
  contacts:
    polymer_selection: "resname SBM EGM"
    cutoff: 4.5
    
    compute_binding_preference: true
    surface_exposure_threshold: 0.2
    enzyme_pdb_for_sasa: "structures/lipase.pdb"
    include_default_aa_groups: true
    
    # Custom groups for mechanistic insight
    protein_groups:
      catalytic_triad: [77, 133, 156]
      oxyanion_hole: [10, 77]
      lid_helix: [140, 141, 142, 143, 144, 145, 146, 147, 148, 149]
```

### Run Analysis

```bash
polyzymd compare contacts -f comparison.yaml
```

### Interpretation Framework

| Observation | Interpretation | Design Implication |
|-------------|----------------|-------------------|
| EGMA aromatic enrichment > +1.0 | Strong π-stacking with surface Trp/Phe/Tyr | May bind near aromatic-rich active sites |
| SBMA charged_positive ≈ 0 | Neutral towards Lys/Arg despite zwitterionic nature | Balanced electrostatics |
| SBMA catalytic_triad < -0.5 | Avoids active site | Good for maintaining catalytic activity |
| EGMA lid_helix > +0.5 | Preferentially binds lid domain | May affect lid dynamics |

## API Reference

### compute_binding_preference

```python
from polyzymd.analysis.contacts.binding_preference import (
    compute_binding_preference,
    extract_polymer_composition,
)

# First, extract polymer composition from trajectory
polymer_composition = extract_polymer_composition(
    universe,                    # MDAnalysis Universe with trajectory
    polymer_type_selections=None,  # Optional: custom polymer type definitions
)

# Then compute binding preference
result = compute_binding_preference(
    contact_result,      # ContactResult from trajectory analysis
    surface_exposure,    # SurfaceExposureResult from SASA calculation
    protein_groups,      # dict[str, set[int]] of group name -> resids
    polymer_composition, # PolymerComposition from extract_polymer_composition()
    polymer_types=None,  # Optional filter for polymer types
)
```

### BindingPreferenceResult Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `enrichment_matrix()` | `dict[str, dict[str, float]]` | Enrichment values (zero-centered) |
| `contact_fraction_matrix()` | `dict[str, dict[str, float]]` | Mean contact fractions |
| `contact_share_matrix()` | `dict[str, dict[str, float]]` | Contact shares |
| `get_enrichment(polymer, group)` | `float` | Single enrichment value |
| `get_entry(polymer, group)` | `BindingPreferenceEntry` | Full entry with all metrics |
| `to_dataframe()` | `pd.DataFrame` | Convert to pandas for analysis |
| `save(path)` | `None` | Save to JSON |
| `load(path)` | `BindingPreferenceResult` | Load from JSON |

### BindingPreferenceEntry Fields

| Field | Type | Description |
|-------|------|-------------|
| `polymer_type` | str | Polymer residue type (e.g., "SBM") |
| `protein_group` | str | Protein group label (e.g., "aromatic") |
| `total_contact_frames` | int | Sum of contact frames for group |
| `mean_contact_fraction` | float | Average per-residue contact fraction |
| `n_residues_in_group` | int | Total residues in group (all) |
| `n_exposed_in_group` | int | Surface-exposed residues in group |
| `n_residues_contacted` | int | Exposed residues with any contact |
| `contact_share` | float | Fraction of polymer's contacts to group |
| `expected_share` | float | Expected share (surface availability) |
| `enrichment` | float | Zero-centered enrichment |
| `polymer_residue_count` | int | Residues of this polymer type (metadata) |
| `total_polymer_residues` | int | Total polymer residues (metadata) |
| `polymer_heavy_atom_count` | int | Heavy atoms of this polymer type (metadata) |
| `total_polymer_heavy_atoms` | int | Total polymer heavy atoms (metadata) |

### PolymerComposition

```python
from polyzymd.analysis.contacts.binding_preference import (
    PolymerComposition,
    extract_polymer_composition,
)

# Auto-extract from trajectory
composition = extract_polymer_composition(universe)

# Or create manually
composition = PolymerComposition(
    residue_counts={"SBM": 50, "EGM": 50},
    heavy_atom_counts={"SBM": 750, "EGM": 400},
)

# Access properties
print(composition.total_residues)        # 100
print(composition.residue_fraction("SBM"))  # 0.5
print(composition.heavy_atom_fraction("SBM"))  # 0.652 (larger monomer)
```

## System Coverage Analysis

While binding preference answers *"What does each polymer type prefer?"*,
**system coverage** answers a complementary question: *"What does this polymer
mixture collectively cover?"*

### When to Use System Coverage

System coverage is useful when comparing **copolymer compositions** (conditions)
rather than individual polymer types:

- *"Does a 70:30 SBMA:EGMA mixture cover aromatic residues differently than
  a 30:70 mixture?"*
- *"Which copolymer ratio provides better coverage of the active site?"*
- *"How does aggregate polymer coverage change across my experimental conditions?"*

### The Coverage Formula

For each protein group, system coverage computes:

```{math}
\text{Coverage Share} = \frac{\sum_{\text{all polymers}} \text{contacts to group}}{\sum_{\text{all polymers}} \text{total contacts}}
```

```{math}
\text{Coverage Enrichment} = \frac{\text{Coverage Share}}{\text{Expected Share}} - 1
```

where Expected Share is based on protein surface availability (same as binding
preference).

### How It Differs from Binding Preference

| Aspect | Binding Preference | System Coverage |
|--------|-------------------|-----------------|
| **Question** | What does SBMA prefer? | What does this mixture cover? |
| **Scope** | Per polymer type | Collapsed across all polymers |
| **Use case** | Compare polymer chemistries | Compare copolymer ratios |
| **Comparable across** | Conditions (abundance cancels) | Conditions (aggregate behavior) |

```{important}
Both metrics are **automatically computed** when you enable binding preference.
System coverage is stored in `result.system_coverage`.
```

### Partition-Based Architecture (Schema v2)

System coverage uses a **partition-based architecture** that ensures mathematically
correct enrichment calculations. Each partition is a collection of **mutually
exclusive** protein groups that together cover the protein surface.

```{important}
**Why Partitions Must Be Mutually Exclusive**

The enrichment formula requires that both `coverage_share` and `expected_share`
sum to exactly 1.0 across all groups in a partition:

$$\sum_{\text{groups}} \text{coverage\_share}_i = 1.0$$
$$\sum_{\text{groups}} \text{expected\_share}_i = 1.0$$

If groups **overlap** (share residues), the `expected_share` sum exceeds 1.0,
causing systematically negative enrichments—a subtle but serious bug.
```

#### Built-in Partitions

System coverage automatically computes several partition types:

| Partition | Groups | Description |
|-----------|--------|-------------|
| **AA Class** | aromatic, polar, nonpolar, charged_positive, charged_negative | 5-way classification by amino acid chemistry |
| **Custom Binary** | {custom_group, rest_of_protein} | One partition per custom group in `protein_groups` |
| **Combined Custom** | {all custom groups, rest_of_protein} | Only if custom groups don't overlap |
| **User-Defined** | User-specified groups from `protein_partitions` | Configurable via YAML |

### User-Defined Partitions

For specialized analyses, you can define custom partitions that organize protein
regions into mutually exclusive groups:

`````{tab-set}
````{tab-item} YAML
```yaml
# comparison.yaml
analysis_settings:
  contacts:
    # First, define individual protein groups (can overlap)
    protein_groups:
      lid_helix_5: "resid 141:155"
      lid_helix_10: "resid 281:295"
      catalytic_triad: "resid 131 163 194"
      active_site_loop: "resid 75:85"
      
    # Then, define partitions (must be mutually exclusive within each partition)
    protein_partitions:
      lid_helices:
        - lid_helix_5
        - lid_helix_10
      catalytic_regions:
        - catalytic_triad
        - active_site_loop
```
````

````{tab-item} Python
```python
from polyzymd.compare.settings import ContactsAnalysisSettings

settings = ContactsAnalysisSettings(
    protein_groups={
        "lid_helix_5": "resid 141:155",
        "lid_helix_10": "resid 281:295",
        "catalytic_triad": "resid 131 163 194",
        "active_site_loop": "resid 75:85",
    },
    protein_partitions={
        "lid_helices": ["lid_helix_5", "lid_helix_10"],
        "catalytic_regions": ["catalytic_triad", "active_site_loop"],
    },
)
```
````
`````

```{note}
**Automatic `rest_of_protein`**: If your partition groups don't cover all
exposed protein residues, a `rest_of_protein` element is automatically added
to ensure the partition sums to 1.0.
```

#### Validation at Config Load Time

PolyzyMD validates partition definitions when the config is loaded, **failing
fast** with a clear error if groups overlap:

```python
# This will raise ValidationError at load time:
protein_partitions:
  bad_partition:
    - lid_helix_5           # resid 141:155
    - lid_domain_combined   # resid 130:160 — overlaps with lid_helix_5!
```

Error message:
```
ValidationError: Partition 'bad_partition' has overlapping groups: 
lid_helix_5 and lid_domain_combined share residues {141, 142, ..., 155}.
Partitions must contain mutually exclusive groups.
```

### Accessing System Coverage

`````{tab-set}
````{tab-item} Python
```python
from polyzymd.analysis.contacts import BindingPreferenceResult

# Load binding preference result (includes system coverage)
result = BindingPreferenceResult.load("binding_preference_rep1.json")

# Access system coverage
if result.system_coverage:
    coverage = result.system_coverage
    
    # Get AA class coverage enrichment
    print("AA Class Coverage Enrichment:")
    for entry in coverage.aa_class_coverage.entries:
        print(f"  {entry.partition_element}: {entry.coverage_enrichment:+.2f}")
    
    # Get custom group enrichment (vs rest_of_protein)
    print("\nCustom Group Coverage:")
    for group_name in coverage.custom_group_names():
        enrichment = coverage.get_custom_group_enrichment(group_name)
        print(f"  {group_name}: {enrichment:+.2f}")
    
    # Access user-defined partitions
    print("\nUser-Defined Partitions:")
    for partition_name in coverage.user_partition_names():
        partition = coverage.get_user_partition(partition_name)
        print(f"  {partition_name}:")
        for entry in partition.entries:
            print(f"    {entry.partition_element}: {entry.coverage_enrichment:+.2f}")
```
````

````{tab-item} JSON Structure (Schema v2)
System coverage uses a partition-based JSON structure:

```json
{
  "entries": [...],
  "polymer_composition": {...},
  "system_coverage": {
    "aa_class_coverage": {
      "partition_name": "aa_class",
      "entries": [
        {
          "partition_element": "aromatic",
          "total_contact_frames": 25432,
          "coverage_share": 0.164,
          "expected_share": 0.086,
          "coverage_enrichment": 0.907,
          "n_exposed_in_group": 7,
          "n_residues_in_group": 20,
          "polymer_contributions": {"SBM": 0.45, "EGM": 0.55}
        },
        {"partition_element": "polar", ...},
        {"partition_element": "nonpolar", ...},
        {"partition_element": "charged_positive", ...},
        {"partition_element": "charged_negative", ...}
      ]
    },
    "custom_group_coverages": {
      "lid_helix_5": {
        "partition_name": "lid_helix_5_vs_rest",
        "entries": [
          {"partition_element": "lid_helix_5", ...},
          {"partition_element": "rest_of_protein", ...}
        ]
      }
    },
    "user_defined_partitions": {
      "lid_helices": {
        "partition_name": "lid_helices",
        "entries": [
          {"partition_element": "lid_helix_5", ...},
          {"partition_element": "lid_helix_10", ...},
          {"partition_element": "rest_of_protein", ...}
        ]
      }
    },
    "n_frames": 10000,
    "total_contact_frames": 154892,
    "total_exposed_residues": 81,
    "polymer_types_included": ["SBM", "EGM"],
    "schema_version": 2
  },
  "schema_version": 4
}
```
````
`````

### Interpreting Polymer Contributions

The `polymer_contributions` field shows how much each polymer type contributed
to the coverage of a specific protein group. For example:

```python
# For aromatic coverage:
polymer_contributions = {"SBM": 0.45, "EGM": 0.55}
# 45% of contacts to aromatics came from SBMA
# 55% of contacts to aromatics came from EGMA
```

This helps understand:
- Which polymer is "doing the work" for each protein group
- Whether coverage is balanced across polymer types
- How polymer ratio affects coverage patterns

### Aggregation Across Replicates

System coverage is automatically aggregated when you aggregate binding preference:

```python
from polyzymd.analysis.contacts import (
    aggregate_binding_preference,
    AggregatedBindingPreferenceResult,
)

# Aggregate multiple replicates
aggregated = aggregate_binding_preference([result1, result2, result3])

# Access aggregated system coverage
if aggregated.system_coverage:
    # AA class coverage with statistics
    for entry in aggregated.system_coverage.aa_class_coverage.entries:
        print(f"{entry.partition_element}:")
        print(f"  Mean enrichment: {entry.mean_coverage_enrichment:+.2f}")
        print(f"  SEM: ±{entry.sem_coverage_enrichment:.2f}")
        print(f"  Replicates: {entry.n_replicates}")
    
    # User-defined partitions with statistics
    for partition_name in aggregated.system_coverage.user_partition_names():
        partition = aggregated.system_coverage.get_user_partition(partition_name)
        print(f"\n{partition_name}:")
        for entry in partition.entries:
            print(f"  {entry.partition_element}: {entry.mean_coverage_enrichment:+.2f} ± {entry.sem_coverage_enrichment:.2f}")
```

### SystemCoverageResult Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `aa_class_enrichment_dict()` | `dict[str, float]` | {aa_class: enrichment} mapping |
| `get_aa_class_enrichment(cls)` | `float` | Enrichment for specific AA class |
| `get_custom_group_enrichment(name)` | `float` | Enrichment for custom group |
| `custom_group_names()` | `list[str]` | List of custom group names |
| `user_partition_names()` | `list[str]` | List of user-defined partition names |
| `get_user_partition(name)` | `PartitionCoverageResult` | Get user-defined partition |
| `save(path)` / `load(path)` | | JSON serialization |

### PartitionCoverageEntry Fields

| Field | Type | Description |
|-------|------|-------------|
| `partition_element` | str | Group name within the partition |
| `total_contact_frames` | int | Sum of ALL polymer contacts to group |
| `coverage_share` | float | Fraction of all contacts to this group |
| `expected_share` | float | Expected share (surface availability) |
| `coverage_enrichment` | float | Zero-centered enrichment |
| `n_exposed_in_group` | int | Surface-exposed residues in group |
| `n_residues_in_group` | int | Total residues in group |
| `polymer_contributions` | dict[str, float] | Fraction from each polymer type |

### Why Overlapping Groups Cause Bugs

Understanding why partitions must be mutually exclusive helps avoid subtle errors:

```{admonition} Worked Example: The Overlap Bug
:class: warning

**Setup**: Consider two groups that overlap:
- `lid_helix_5`: residues {141, 142, 143, 144, 145}
- `lid_domain`: residues {140, 141, 142, 143, 144, 145, 146, 147}

**Expected shares** (if treated as a "partition"):
- `lid_helix_5`: 5 exposed residues → expected_share = 5/100 = 0.05
- `lid_domain`: 8 exposed residues → expected_share = 8/100 = 0.08
- **Sum: 0.13** (should be 1.0!)

The 5 overlapping residues are **double-counted** in expected_share.

**Coverage shares** (actual contacts):
- `lid_helix_5`: 1000 contact-frames → coverage_share = 1000/20000 = 0.05
- `lid_domain`: 1600 contact-frames → coverage_share = 1600/20000 = 0.08
- **Sum: 0.13** (contacts aren't double-counted if we sum over groups)

Wait—but contacts **are** double-counted! Those 5 overlapping residues contribute
to **both** groups. So coverage_share also sums to >1.0? No—the denominator
(total contacts) is shared, making coverage_shares additive only if groups
don't overlap.

**Result**: enrichment = (0.05 / 0.05) - 1 = 0, but the "expected" baseline
is wrong because the partition doesn't sum to 1.0. Comparisons between groups
become meaningless.

**Solution**: Use mutually exclusive partitions. PolyzyMD validates this at
config load time.
```

## Visualization

Both binding preference and system coverage generate publication-ready plots
when you run `polyzymd compare contacts`. These plots are automatically
enabled by default but can be controlled via `plot_settings`.

### Available Plots

| Plot Type | Description | Settings Field |
|-----------|-------------|----------------|
| **Binding Preference Heatmap** | Enrichment by (polymer type × protein group) across conditions | `generate_enrichment_heatmap` |
| **Binding Preference Bars** | Grouped bars per polymer type, comparing conditions | `generate_enrichment_bars` |
| **System Coverage Heatmap** | Aggregate coverage by protein group across conditions | `generate_system_coverage_heatmap` |
| **System Coverage Bars** | Grouped bars showing aggregate coverage across conditions | `generate_system_coverage_bars` |
| **User Partition Bars** | One grouped bar chart per user-defined partition | `generate_user_partition_bars` |

### Configuring Plot Settings

`````{tab-set}
````{tab-item} YAML
```yaml
# comparison.yaml
plot_settings:
  output_dir: "figures/"
  format: "png"
  dpi: 300
  
  contacts:
    # Binding preference plots
    generate_enrichment_heatmap: true
    generate_enrichment_bars: true
    figsize_enrichment_heatmap: [12, 8]  # Auto-calculated if null
    figsize_enrichment_bars: [10, 6]
    enrichment_colormap: "RdBu_r"  # Diverging colormap
    show_enrichment_error: true
    
    # System coverage plots (AA class partition)
    generate_system_coverage_heatmap: true
    generate_system_coverage_bars: true
    figsize_system_coverage_heatmap: [10, 6]  # Auto-calculated if null
    figsize_system_coverage_bars: [10, 6]
    show_system_coverage_error: true
    
    # User-defined partition plots
    generate_user_partition_bars: true
    figsize_user_partition_bars: [10, 6]
    show_user_partition_error: true
```
````

````{tab-item} Python
```python
from polyzymd.compare.config import PlotSettings, ContactsPlotSettings

plot_settings = PlotSettings(
    output_dir="figures/",
    format="png",
    dpi=300,
    contacts=ContactsPlotSettings(
        generate_enrichment_heatmap=True,
        generate_system_coverage_heatmap=True,
        enrichment_colormap="RdBu_r",
        show_enrichment_error=True,
        show_system_coverage_error=True,
        # User partition plots
        generate_user_partition_bars=True,
        show_user_partition_error=True,
    ),
)
```
````
`````

### Understanding the Heatmaps

**Binding Preference Heatmap:**
- One subplot per condition
- Rows: Protein groups (aromatic, polar, etc.)
- Columns: Polymer types (SBMA, EGMA, etc.)
- Color: Zero-centered diverging colormap
  - **Red (positive)**: Preferential binding
  - **White (zero)**: Neutral
  - **Blue (negative)**: Avoidance

**System Coverage Heatmap:**
- Single heatmap comparing conditions
- Rows: Protein groups
- Columns: Conditions (100% SBMA, 50/50 copolymer, etc.)
- Shows aggregate coverage across all polymer types

### Understanding the Bar Charts

**Binding Preference Bars:**
- One plot per polymer type
- Groups: Protein groups
- Bars within group: One per condition
- Error bars: SEM across replicates

**System Coverage Bars:**
- Single plot comparing all conditions
- Groups: AA class groups (aromatic, polar, etc.)
- Bars within group: One per condition
- Shows how different copolymer compositions collectively cover each protein group

**User Partition Bars:**
- One plot per user-defined partition
- Groups: Partition elements (e.g., lid_helix_5, lid_helix_10, rest_of_protein)
- Bars within group: One per condition
- Useful for comparing coverage of specific structural regions across conditions

### Disabling Specific Plots

To disable specific plot types, set the corresponding field to `false`:

```yaml
plot_settings:
  contacts:
    # Only generate bar charts, not heatmaps
    generate_enrichment_heatmap: false
    generate_system_coverage_heatmap: false
    generate_enrichment_bars: true
    generate_system_coverage_bars: true
    generate_user_partition_bars: true  # Enable user partition plots
```

## Troubleshooting

### "No exposed residues in group"

This warning appears when a protein group has no surface-exposed residues
(all are buried). Solutions:

1. Lower `surface_exposure_threshold` to include more residues
2. Verify the residue IDs in your custom `protein_groups` are correct
3. Check that your enzyme PDB has correct atom coordinates

### "Enzyme PDB not found"

The SASA calculation requires a PDB file. Ensure:

1. `enzyme_pdb_for_sasa` path is relative to `comparison.yaml` location
2. The file exists and has correct format
3. Alternatively, place a symlink in the `structures/` directory

### Enrichment values of -1.0

An enrichment of exactly -1.0 means **complete avoidance** (no contacts
observed for this polymer-group pair). This is expected for:
- Buried residue groups with no surface exposure
- Polymer types that physically cannot reach certain regions

Check the `n_exposed_in_group` and `n_residues_contacted` fields to diagnose.

### High SEM values

Large uncertainty (SEM) indicates high variability across replicates. This
could mean:

1. Insufficient sampling — consider longer simulations
2. True variability in binding — the polymer explores multiple binding modes
3. Rare events dominating one replicate — check per-replicate values

## See Also

- [Contacts Analysis Quick Start](analysis_contacts_quickstart.md) — basic contacts analysis
- [Contacts Analysis Cookbook](analysis_contacts_cookbook.md) — worked examples
- [Surface Exposure Calculation](analysis_contacts_cookbook.md#recipe-4-active-site-contacts-vs-surface-contacts) — SASA filtering
- [Statistics Best Practices](analysis_statistics_best_practices.md) — uncertainty quantification
- [Comparing Conditions](analysis_compare_conditions.md) — multi-condition workflows
