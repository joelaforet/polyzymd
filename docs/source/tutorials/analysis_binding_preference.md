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

The **expected share** can be normalized in two ways (see [Dual Normalization](#dual-normalization) below).

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

### Dual Normalization

Two normalization methods are provided to distinguish **chemical affinity**
from **geometric/steric effects**:

**Residue-based normalization** (default, `by="residue"`):

```{math}
\text{Expected Share} = \frac{\text{polymer residue count}}{\text{total polymer residues}}
```

This matches the **experimental viewpoint** where polymer concentrations are
specified in terms of monomer units (e.g., "0.2mM SBMA : 0.8mM EGMA = 1:4 ratio").

**Atom-based normalization** (`by="atoms"`):

```{math}
\text{Expected Share} = \frac{\text{polymer heavy atoms}}{\text{total polymer heavy atoms}}
```

This accounts for **monomer size differences**, since larger monomers have
more surface area and thus more contact opportunities.

### Interpreting Dual Metrics

| `enrichment_by_residue` | `enrichment_by_atoms` | Interpretation |
|-------------------------|----------------------|----------------|
| **Positive** | **Positive** | Strong evidence of **chemical preference** |
| **Negative** | **Negative** | Strong evidence of **avoidance** |
| **Positive** | **~0** | Enrichment explained by **larger monomer size** |
| **~0** | **Positive** | Smaller monomer "**punches above its weight**" |

```{tip}
Use residue-based enrichment for direct comparison with experimental expectations.
Use atom-based enrichment to reveal true chemical affinity vs. steric effects.
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

# Extract polymer composition from trajectory (required for dual normalization)
from polyzymd.analysis.contacts.binding_preference import extract_polymer_composition
polymer_composition = extract_polymer_composition(universe)  # universe from trajectory

# Compute binding preference
result = compute_binding_preference(
    contact_result=contacts,
    surface_exposure=surface_exposure,
    protein_groups=protein_groups,
    polymer_composition=polymer_composition,
)

# Access enrichment values (default: by residue)
print("Enrichment matrix (by residue):")
for polymer_type, groups in result.enrichment_matrix(by="residue").items():
    print(f"\n{polymer_type}:")
    for group, enrichment in groups.items():
        marker = "+" if enrichment > 0 else "-" if enrichment < 0 else "="
        print(f"  {group}: {enrichment:+.2f} {marker}")

# Compare with atom-based normalization
print("\nEnrichment matrix (by atoms):")
for polymer_type, groups in result.enrichment_matrix(by="atoms").items():
    print(f"\n{polymer_type}:")
    for group, enrichment in groups.items():
        marker = "+" if enrichment > 0 else "-" if enrichment < 0 else "="
        print(f"  {group}: {enrichment:+.2f} {marker}")
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
Enrichment normalized by: residue count

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
| `enrichment_normalization` | str | `"residue"` | Normalization method: `"residue"` or `"atoms"` |
| `polymer_type_selections` | dict | `null` | Custom polymer type definitions (see below) |

### Enrichment Normalization

The `enrichment_normalization` setting controls which normalization method is
used for display and plotting:

```yaml
analysis_settings:
  contacts:
    compute_binding_preference: true
    enrichment_normalization: "residue"  # Default: matches experimental ratios
    # enrichment_normalization: "atoms"  # Alternative: accounts for monomer size
```

Both normalizations are always computed and stored in the output. This setting
only controls which one is displayed in plots and console output.

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
      "polymer_residue_count": 50,
      "total_polymer_residues": 100,
      "expected_share_by_residue": 0.5,
      "enrichment_by_residue": 0.328,
      "polymer_heavy_atom_count": 750,
      "total_polymer_heavy_atoms": 1150,
      "expected_share_by_atoms": 0.652,
      "enrichment_by_atoms": -0.254
    }
  ],
  "polymer_composition": {
    "residue_counts": {"SBM": 50, "EGM": 50},
    "heavy_atom_counts": {"SBM": 750, "EGM": 400}
  },
  "n_frames": 10000,
  "total_exposed_residues": 81,
  "surface_exposure_threshold": 0.2,
  "schema_version": 2
}
```

```{note}
The schema version has been updated to `2` to reflect the new dual-normalization
fields. Both `enrichment_by_residue` and `enrichment_by_atoms` are always
computed and stored, regardless of the `enrichment_normalization` display setting.
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
    enrichment_normalization: "residue"  # Use residue-based for comparison
    
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
| `enrichment_matrix(by="residue")` | `dict[str, dict[str, float]]` | Enrichment values (zero-centered) |
| `enrichment_matrix(by="atoms")` | `dict[str, dict[str, float]]` | Atom-normalized enrichment |
| `contact_fraction_matrix()` | `dict[str, dict[str, float]]` | Mean contact fractions |
| `contact_share_matrix()` | `dict[str, dict[str, float]]` | Contact shares |
| `get_enrichment(polymer, group)` | `float` | Single enrichment (residue-based) |
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
| `polymer_residue_count` | int | Residues of this polymer type |
| `total_polymer_residues` | int | Total polymer residues (all types) |
| `expected_share_by_residue` | float | Expected share by residue count |
| `enrichment_by_residue` | float | Zero-centered enrichment (residue-based) |
| `polymer_heavy_atom_count` | int | Heavy atoms of this polymer type |
| `total_polymer_heavy_atoms` | int | Total polymer heavy atoms (all types) |
| `expected_share_by_atoms` | float | Expected share by atom count |
| `enrichment_by_atoms` | float | Zero-centered enrichment (atom-based) |

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
