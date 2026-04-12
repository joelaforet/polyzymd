# Polymer-Protein Contacts Analysis: Quick Start

Analyze polymer-protein contact frequencies and coverage in your MD simulations.

```{note}
This command analyzes contacts between polymer chains (Chain C) and protein 
residues (Chain A), following PolyzyMD's chain convention. It computes per-residue 
contact fractions with proper statistical treatment.
```

## TL;DR

```bash
# Initialize a comparison workspace
polyzymd compare init -n contacts_study
cd contacts_study

# Configure plugins.contacts in comparison.yaml, then run contacts
polyzymd compare run contacts

# Or run all enabled analyses
polyzymd compare run-all
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

By default, contacts analysis uses the selections configured in
`plugins.contacts` in `comparison.yaml`.

## Basic Usage

### Run Contacts Analysis

`````{tab-set}
````{tab-item} YAML (Recommended)
Create a `comparison.yaml` file and configure `plugins.contacts`:

```yaml
# comparison.yaml
name: "contacts_study"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]

  - label: "SBMA"
    config: "../sbma_100/config.yaml"
    replicates: [1, 2, 3]

plugins:
  contacts:
    polymer_selection: "chainID C"
    protein_selection: "protein"
    cutoff: 4.5
    compute_residence_times: true
```

Then run contacts analysis:

```bash
polyzymd compare run contacts
```
````

````{tab-item} CLI
```bash
cd /path/to/contacts_study
polyzymd compare run contacts -f comparison.yaml
```
````

````{tab-item} Python
```python
# Contacts analysis runs through the comparison pipeline.
# See the YAML tab for the recommended workflow.
#
# To load results programmatically after running:
import json
from pathlib import Path

result = json.loads(Path("analysis/<condition>/contacts/run_1/result.json").read_text())

print(f"Coverage: {result['coverage_fraction']:.1%}")
```
````
`````

```{note}
Contacts analysis requires at least 2 replicates for aggregation and
statistical comparison. Configure `replicates: [1, 2, 3]` (or more) in
each condition for robust results.
```

**Output:**
```
Loading configuration from: config.yaml
Contact Analysis: MyEnzyme_Polymer_Study
  Replicates: 3
  Equilibration: 10ns
  Cutoff: 4.5 Å
  Polymer selection: chainID C
  Protein selection: protein
  Grouping: aa_class
  Processing replicate 1... done (134/181 residues contacted, 74.0% coverage, 18.0% mean contact)
  Processing replicate 2... done (136/181 residues contacted, 75.1% coverage, 18.5% mean contact)
  Processing replicate 3... done (132/181 residues contacted, 72.9% coverage, 17.6% mean contact)

Contact Analysis Complete
  Contacted residues: 134/181 (mean)
  Coverage: 74.0%
  Mean contact fraction: 18.0%
  Results saved: analysis/<condition>/contacts/run_<N>/result.json
```

### Key Metrics

- **Coverage**: Fraction of protein residues that had at least one contact with 
  polymer during the trajectory
- **Mean contact fraction**: Average fraction of frames where each residue was 
  in contact with polymer (averaged across all residues)
- **Contacted residues**: Count of residues with any polymer contact

## Using comparison.yaml

For reproducible, version-controlled contacts analysis, configure
`plugins.contacts` in `comparison.yaml` and run through `polyzymd compare`.

### Multi-Replicate Analysis with Full Options

`````{tab-set}
````{tab-item} YAML (Recommended)
Create a `comparison.yaml` file with all contacts options:

```yaml
# comparison.yaml
name: "contacts_study"
control: "No Polymer"

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]

  - label: "SBMA"
    config: "../sbma_100/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

plugins:
  contacts:
    polymer_selection: "chainID C"      # MDAnalysis selection for polymer
    protein_selection: "protein"        # MDAnalysis selection for protein
    cutoff: 4.5                          # Contact distance in Angstroms
    polymer_types: ["SBM", "EGM"]      # Optional: filter by monomer type
    grouping: "aa_class"                # aa_class | secondary_structure | none
    compute_residence_times: true        # Enable residence time statistics
```

Then run contacts analysis:

```bash
# Run contacts only
polyzymd compare run contacts

# Or run all enabled analyses
polyzymd compare run-all
```
````

````{tab-item} CLI
```bash
# Run contacts using settings from comparison.yaml
polyzymd compare run contacts -f comparison.yaml

# Run all enabled analyses
polyzymd compare run-all -f comparison.yaml
```
````

````{tab-item} Python
```python
# Contacts analysis runs through the comparison pipeline.
# See the YAML tab for the recommended workflow.
#
# To load aggregated results programmatically after running:
import json
from pathlib import Path

aggregated = json.loads(
    Path("analysis/<condition>/contacts/aggregated/result.json").read_text()
)
print(
    f"Contact fraction: {aggregated['mean_contact_fraction']:.1%} "
    f"± {aggregated['mean_contact_fraction_sem']:.1%}"
)
```
````
`````

### Configuration Options

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `polymer_selection` | str | `"chainID C"` | MDAnalysis selection for polymer atoms |
| `protein_selection` | str | `"protein"` | MDAnalysis selection for protein atoms |
| `cutoff` | float | `4.5` | Distance cutoff in Angstroms |
| `polymer_types` | list | `null` | Filter by polymer residue names |
| `grouping` | str | `"aa_class"` | How to group protein residues |
| `compute_residence_times` | bool | `true` | Compute residence time statistics |

```{tip}
**Workflow recommendation:** Configure contacts settings in `comparison.yaml`
and run `polyzymd compare run contacts` for reproducible analyses.
```

## Residence Time Analysis

Residence time measures how long polymer segments remain in contact with protein 
residues. This is crucial for understanding binding kinetics and comparing 
different polymer types.

### Enable Residence Time Statistics

`````{tab-set}
````{tab-item} YAML (Recommended)
```yaml
# comparison.yaml (excerpt)
plugins:
  contacts:
    compute_residence_times: true  # Enable residence time statistics
```

```bash
polyzymd compare run contacts
```
````

````{tab-item} CLI
```bash
polyzymd compare run contacts -f comparison.yaml
```
````

````{tab-item} Python
```python
# Contacts analysis runs through the comparison pipeline.
# See the YAML tab for the recommended workflow.
#
# To load results programmatically after running:
import json
from pathlib import Path

result = json.loads(Path("analysis/<condition>/contacts/run_1/result.json").read_text())

# Access residence time data
for polymer_type, stats in result["residence_times"].items():
    print(f"{polymer_type}: mean={stats['mean']:.2f} frames, max={stats['max']} frames ({stats['n_events']} events)")
```
````
`````

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

`````{tab-set}
````{tab-item} YAML (Recommended)
```yaml
# comparison.yaml (excerpt)
conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]
  - label: "SBMA"
    config: "../sbma_100/config.yaml"
    replicates: [1, 2, 3]

plugins:
  contacts:
    compute_residence_times: true
```

```bash
polyzymd compare run contacts
```
````

````{tab-item} CLI
```bash
polyzymd compare run contacts -f comparison.yaml
```
````

````{tab-item} Python
```python
# Contacts analysis runs through the comparison pipeline.
# See the YAML tab for the recommended workflow.
#
# To load aggregated results programmatically after running:
import json
from pathlib import Path

aggregated = json.loads(
    Path("analysis/<condition>/contacts/aggregated/result.json").read_text()
)

print(f"Contact fraction: {aggregated['mean_contact_fraction']:.1%} ± {aggregated['mean_contact_fraction_sem']:.1%}")
for polymer_type, stats in aggregated["residence_time_by_polymer_type"].items():
    print(f"{polymer_type}: {stats['mean']:.2f} ± {stats['std']:.2f} frames")
```
````
`````

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

## Command Pattern

Contacts analysis is configured in `comparison.yaml` under `plugins.contacts`.
Run with:

```bash
polyzymd compare run contacts -f comparison.yaml
```

## Output Files

Results are saved as JSON in the project's analysis directory:

```
project/
└── analysis/
    └── <condition>/
        └── contacts/
            ├── run_1/
            │   └── result.json
            ├── run_2/
            │   └── result.json
            ├── run_3/
            │   └── result.json
            └── aggregated/
                └── result.json
```

The JSON file contains per-residue contact data including:
- Contact events (start frame, duration)
- Contact fractions
- Amino acid classifications
- Polymer chain interactions

## Custom Selections

### Analyze specific polymer types

`````{tab-set}
````{tab-item} YAML (Recommended)
```yaml
# comparison.yaml (excerpt) - Only SBMA monomers
plugins:
  contacts:
    polymer_selection: "segid C and resname SBM"  # Only SBMA
    protein_selection: "protein"
```

```bash
polyzymd compare run contacts
```

To analyze EGMA instead, change the selection:

```yaml
contacts:
  polymer_selection: "segid C and resname EGM"  # Only EGMA
```
````

````{tab-item} CLI
```bash
# Only SBMA monomers (configured in comparison.yaml)
polyzymd compare run contacts -f comparison.yaml

# Only EGMA monomers: update polymer_selection in comparison.yaml, then rerun
polyzymd compare run contacts -f comparison.yaml
```
````

````{tab-item} Python
```python
# Contacts analysis runs through the comparison pipeline.
# See the YAML tab for the recommended workflow.
#
# To load results programmatically after running:
import json
from pathlib import Path

result = json.loads(Path("analysis/<condition>/contacts/run_1/result.json").read_text())
print(f"Coverage: {result['coverage_fraction']:.1%}")
```
````
`````

### Analyze specific protein regions

`````{tab-set}
````{tab-item} YAML (Recommended)
```yaml
# comparison.yaml (excerpt) - Only aromatic residues
plugins:
  contacts:
    polymer_selection: "segid C"
    protein_selection: "protein and (resname TRP PHE TYR)"  # Aromatics only
```

```bash
polyzymd compare run contacts
```

For active site analysis:

```yaml
contacts:
  protein_selection: "protein and (resid 75-80 or resid 130-140)"  # Active site
```
````

````{tab-item} CLI
```bash
# Only aromatic residues (configured in comparison.yaml)
polyzymd compare run contacts -f comparison.yaml

# Active site region: update protein_selection in comparison.yaml, then rerun
polyzymd compare run contacts -f comparison.yaml
```
````

````{tab-item} Python
```python
# Contacts analysis runs through the comparison pipeline.
# See the YAML tab for the recommended workflow.
#
# To load results programmatically after running:
import json
from pathlib import Path

result = json.loads(Path("analysis/<condition>/contacts/run_1/result.json").read_text())
print(f"Mean contact fraction: {result['mean_contact_fraction']:.1%}")
```
````
`````

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
from polyzymd.analyses.contacts._results import ContactResult

result = ContactResult.load("analysis/<condition>/contacts/run_1/result.json")

# Get contact fraction by (polymer_type, protein_AA_class)
matrix = result.interaction_matrix(metric="contact_fraction")

# Compare polymer types contacting aromatic residues
print(f"SBMA-aromatic: {matrix['SBM']['aromatic']:.1%}")
print(f"EGMA-aromatic: {matrix['EGM']['aromatic']:.1%}")
```

```{note}
The `interaction_matrix()` and `coverage_by_group()` methods require loading
results via `ContactResult.load()` (from the private `_results` module)
rather than plain JSON. This API may change in a future release — prefer
the CLI workflow for standard comparisons.
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

```{seealso}
For a detailed explanation of autocorrelation, statistical inefficiency, and
the LiveCoMS methodology, see the
[Statistics Best Practices Guide](../explanation/analysis_statistics_best_practices.md).
```

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

- [Binding Preference Analysis](analysis_binding_preference.md) - enrichment-based AA class preferences
- [Contacts Analysis Cookbook](analysis_contacts_cookbook.md) - worked examples for scientific questions
- [RMSF Analysis Quick Start](analysis_rmsf_quickstart.md) - complementary stability analysis
- [Catalytic Triad Analysis](analysis_triad_quickstart.md) - active site geometry
- [Comparing Conditions](analysis_compare_conditions.md) - statistical comparison across polymers
