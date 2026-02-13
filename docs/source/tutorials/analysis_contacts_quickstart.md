# Polymer-Protein Contacts Analysis: Quick Start

Analyze polymer-protein contact frequencies and coverage in your MD simulations.

```{note}
This command analyzes contacts between polymer chains (Chain C) and protein 
residues (Chain A), following PolyzyMD's chain convention. It computes per-residue 
contact fractions with proper statistical treatment.
```

## TL;DR

```bash
# Single replicate analysis
polyzymd analyze contacts -c config.yaml -r 1 --eq-time 10ns

# Custom cutoff distance (default: 4.0 Å)
polyzymd analyze contacts -c config.yaml -r 1 --eq-time 10ns --cutoff 5.0

# Force recompute (ignore cache)
polyzymd analyze contacts -c config.yaml -r 1 --eq-time 10ns --recompute
```

## Prerequisites

Before running contacts analysis, you need:

1. **Completed production simulation** - with polymer chains present
2. **A config.yaml file** - pointing to your simulation output
3. **Trajectory files** - in the scratch directory specified in your config
4. **solvated_system.pdb** - topology with correct chain assignments

## PolyzyMD Chain Convention

The contacts analyzer uses PolyzyMD's standard chain assignment:

| Chain | Contents |
|-------|----------|
| A | Protein/Enzyme |
| B | Substrate/Ligand |
| C | Polymers |
| D+ | Solvent (water, ions, co-solvents) |

By default, the command analyzes contacts between `segid C` (polymers) and 
`protein` (chain A).

## Basic Usage

### Single Replicate Analysis

```bash
cd /path/to/your/project
polyzymd analyze contacts -c config.yaml -r 1 --eq-time 10ns
```

**Output:**
```
Loading configuration from: config.yaml
Contact Analysis: MyEnzyme_Polymer_Study
  Replicates: 1
  Equilibration: 10ns
  Cutoff: 4.0 Å
  Polymer selection: segid C
  Protein selection: protein
  Grouping: residue
  Processing replicate 1... done (134/181 residues contacted, 74.0% coverage, 18.0% mean contact)

Contact Analysis Complete
  Contacted residues: 134/181
  Coverage: 74.0%
  Mean contact fraction: 18.0%
  Results saved: /path/to/project/analysis/contacts/contacts_rep1.json
```

### Key Metrics

- **Coverage**: Fraction of protein residues that had at least one contact with 
  polymer during the trajectory
- **Mean contact fraction**: Average fraction of frames where each residue was 
  in contact with polymer (averaged across all residues)
- **Contacted residues**: Count of residues with any polymer contact

## Using analysis.yaml

For reproducible, version-controlled analysis configuration, use `analysis.yaml` 
instead of CLI flags. Place this file alongside your `config.yaml`.

### Create analysis.yaml

```yaml
# analysis.yaml
replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

contacts:
  enabled: true
  polymer_selection: "chainID C"      # MDAnalysis selection for polymer
  protein_selection: "protein"        # MDAnalysis selection for protein
  cutoff: 4.5                         # Contact distance in Angstroms
  polymer_types: ["SBM", "EGM"]       # Optional: filter by monomer type
  grouping: "aa_class"                # aa_class | secondary_structure | none
  compute_residence_times: true       # Enable residence time statistics
```

### Run All Enabled Analyses

```bash
# Initialize a template (if starting fresh)
polyzymd analyze init

# Run all analyses defined in analysis.yaml
polyzymd analyze run

# Force recompute
polyzymd analyze run --recompute
```

### Configuration Options

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `false` | Whether to run contact analysis |
| `polymer_selection` | str | `"chainID C"` | MDAnalysis selection for polymer atoms |
| `protein_selection` | str | `"protein"` | MDAnalysis selection for protein atoms |
| `cutoff` | float | `4.5` | Distance cutoff in Angstroms |
| `polymer_types` | list | `null` | Filter by polymer residue names |
| `grouping` | str | `"aa_class"` | How to group protein residues |
| `compute_residence_times` | bool | `true` | Compute residence time statistics |

```{tip}
**When to use analysis.yaml vs CLI:** Use `analysis.yaml` for standard, 
reproducible workflows. Use CLI flags (`polyzymd analyze contacts ...`) for 
one-off analyses or when exploring different parameters.
```

## Residence Time Analysis

Residence time measures how long polymer segments remain in contact with protein 
residues. This is crucial for understanding binding kinetics and comparing 
different polymer types.

### Enable Residence Time Statistics

```bash
polyzymd analyze contacts -c config.yaml -r 1 --eq-time 10ns --residence-times
```

**Output:**
```
Contact Analysis Complete
  Contacted residues: 138/181
  Coverage: 76.2%
  Mean contact fraction: 18.3%
  Residence time by polymer type:
    EGM: mean=9.12 frames, max=1406 frames (6066 events)
    SBM: mean=8.93 frames, max=842 frames (5401 events)
```

### Multi-Replicate Aggregation

When analyzing multiple replicates, residence times are aggregated with proper 
statistical uncertainty:

```bash
polyzymd analyze contacts -c config.yaml -r 1-3 --eq-time 10ns --residence-times
```

**Output:**
```
Aggregated Contact Analysis Complete
  Contact fraction: 19.8% ± 1.9%
  Residence time by polymer type:
    EGM: 8.14 ± 0.56 frames
    SBM: 9.60 ± 0.53 frames
```

### Interpreting Residence Times

- **Mean residence time**: Average duration of individual contact events (in frames)
- **Max residence time**: Longest continuous contact observed
- **Total events**: Number of separate contact events for this polymer type

Longer residence times suggest stronger or more stable polymer-protein interactions. 
Comparing residence times across polymer types reveals differences in binding 
behavior - for example, zwitterionic polymers (SBMA) often show different residence 
times than PEG-like polymers (EGMA).

### Converting to Physical Time

Residence times are reported in frames by default. To convert to picoseconds, 
multiply by your trajectory's timestep:

```python
# If your trajectory saves every 100 ps:
residence_time_ps = mean_frames * 100  # ps
residence_time_ns = mean_frames * 0.1  # ns
```

The JSON output also includes `mean_ps` and `max_ps` fields computed using the 
timestep from your configuration.

## Command Options

| Option | Default | Description |
|--------|---------|-------------|
| `-c, --config` | Required | Path to config.yaml |
| `-r, --replicates` | `"1"` | Replicate specification: `'1-5'`, `'1,3,5'`, or `'1'` |
| `--eq-time` | `"0ns"` | Equilibration time to skip |
| `--cutoff` | `4.0` | Contact distance cutoff in Angstroms |
| `--polymer-selection` | `"segid C"` | MDAnalysis selection for polymers |
| `--protein-selection` | `"protein"` | MDAnalysis selection for protein |
| `--residence-times` | False | Compute and display residence time statistics by polymer type |
| `--recompute` | False | Force recompute even if cached |
| `-o, --output-dir` | Auto | Custom output directory |

## Output Files

Results are saved as JSON in the project's analysis directory:

```
project/
└── analysis/
    └── contacts/
        └── contacts_rep1.json
```

The JSON file contains per-residue contact data including:
- Contact events (start frame, duration)
- Contact fractions
- Amino acid classifications
- Polymer chain interactions

## Custom Selections

### Analyze specific polymer types

```bash
# Only SBMA monomers
polyzymd analyze contacts -c config.yaml -r 1 --eq-time 10ns \
    --polymer-selection "segid C and resname SBM"

# Only EGMA monomers  
polyzymd analyze contacts -c config.yaml -r 1 --eq-time 10ns \
    --polymer-selection "segid C and resname EGM"
```

### Analyze specific protein regions

```bash
# Only aromatic residues
polyzymd analyze contacts -c config.yaml -r 1 --eq-time 10ns \
    --protein-selection "protein and (resname TRP PHE TYR)"

# Active site region
polyzymd analyze contacts -c config.yaml -r 1 --eq-time 10ns \
    --protein-selection "protein and (resid 75-80 or resid 130-140)"
```

## Answering Scientific Questions

The contacts analysis module is designed to answer questions like:

| Question | Approach |
|----------|----------|
| Do zwitterionic polymers preferentially contact aromatic residues? | `interaction_matrix()` with custom polymer grouping |
| Which protein surface regions does my polymer bind? | `coverage_by_group()` |
| Does SBMA have longer residence times than EGMA? | `residence_time_summary()` |

### Quick Example: Interaction Matrix

The `interaction_matrix()` method computes contact metrics for each combination
of polymer type and protein amino acid class:

```python
from polyzymd.analysis.contacts.results import ContactResult

result = ContactResult.load("analysis/contacts/contacts_rep1.json")

# Get contact fraction by (polymer_type, protein_AA_class)
matrix = result.interaction_matrix(metric="contact_fraction")

# Compare polymer types contacting aromatic residues
print(f"SBMA-aromatic: {matrix['SBM']['aromatic']:.1%}")
print(f"EGMA-aromatic: {matrix['EGM']['aromatic']:.1%}")
```

**Example output:**
```
SBMA-aromatic: 45.4%
EGMA-aromatic: 37.0%
```

### Quick Example: Coverage by AA Group

See what fraction of each amino acid class your polymer contacts:

```python
coverage = result.coverage_by_group()
for group, frac in sorted(coverage.items(), key=lambda x: -x[1]):
    print(f"{group}: {frac:.1%}")
```

**Example output:**
```
charged_negative: 100.0%
aromatic: 100.0%
charged_positive: 100.0%
polar: 93.5%
nonpolar: 86.2%
```

```{tip}
For complete worked examples including custom polymer groupings, residence time
comparisons, and complex queries, see the [Contacts Analysis Cookbook](analysis_contacts_cookbook.md).
```

## Technical Details

### Contact Definition

A contact is defined when any heavy atom of a polymer residue is within the 
cutoff distance of any heavy atom of a protein residue. This uses MDAnalysis's 
`capped_distance` function with KDTree-based neighbor searching for O(N) 
performance.

### Contact Events

Contacts are stored as compressed events `(start_frame, duration)` rather than 
frame-by-frame booleans. This provides efficient storage and enables residence 
time analysis.

### Statistical Treatment

Per-residue contact fractions include proper uncertainty quantification:
- Statistical inefficiency (g) is computed per-residue following Chodera et al. (2007)
- Effective sample sizes account for temporal autocorrelation
- When aggregating across replicates, uncertainties are propagated correctly

## Troubleshooting

### "solvated_system.pdb not found"

The contacts analyzer requires `solvated_system.pdb` in the run directory 
(scratch). This file is created during `polyzymd build` and contains the 
correct chain assignments. Do not use `production_N_topology.pdb` as it may 
have lost chain information.

### "No polymer atoms selected"

Check that your polymer selection is correct. Use MDAnalysis syntax:
- `segid C` - all atoms in segment C (default polymer chain)
- `resname SBM EGM` - atoms with these residue names
- Verify chain assignment with: `polyzymd validate -c config.yaml`

### Slow performance

For very large systems or long trajectories:
- Use `--eq-time` to skip equilibration frames
- Results are cached - subsequent runs load from JSON
- Consider using frame striding in custom scripts

## See Also

- [Contacts Analysis Cookbook](analysis_contacts_cookbook.md) - worked examples for scientific questions
- [RMSF Analysis Quick Start](analysis_rmsf_quickstart.md) - complementary stability analysis
- [Catalytic Triad Analysis](analysis_triad_quickstart.md) - active site geometry
- [Comparing Conditions](analysis_compare_conditions.md) - statistical comparison across polymers
