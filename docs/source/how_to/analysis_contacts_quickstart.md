# Polymer-Protein Contacts Analysis: Quick Start

Analyze polymer-protein contact frequencies and coverage for one or more
conditions using the `contacts` plugin.

## TL;DR

```bash
# Configure plugins.contacts in comparison.yaml, then run:
polyzymd compare run contacts -f comparison.yaml --eq-time 10ns

# Run all enabled analyses in the same workflow
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Force recompute
polyzymd compare run contacts -f comparison.yaml --eq-time 10ns --recompute
```

## Prerequisites

Before running contacts analysis, make sure you have:

1. Completed production trajectories for each replicate
2. A `comparison.yaml` with conditions and `plugins.contacts`
3. Topology with valid chain IDs and polymer atoms
4. At least 2 replicates per condition if you want robust comparison stats

## Chain convention used by contacts

| Chain | Contents |
|-------|----------|
| A | Protein/enzyme |
| B | Substrate/ligand |
| C | Polymer |
| D+ | Solvent and ions |

The default contacts setup expects polymer on chain C and protein on chain A.

## Basic usage

### 1) Configure `comparison.yaml`

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
    polymer_selection: "chainID C"
    protein_selection: "protein"
    cutoff: 4.5
    grouping: "aa_class"
    compute_residence_times: true
```

### 2) Run contacts

```bash
polyzymd compare run contacts -f comparison.yaml --eq-time 10ns
```

Expected output includes per-replicate progress and aggregated summary metrics
(coverage and mean contact fraction).

### 3) Run all enabled plugins (optional)

```bash
polyzymd compare run-all -f comparison.yaml --eq-time 10ns
```

## Key metrics to check first

- **Coverage**: fraction of protein residues contacted at least once
- **Mean contact fraction**: average per-residue fraction of frames in contact
- **Residence time (optional)**: average duration of individual contact events

## Common tasks

### Enable residence time statistics

Residence times are enabled by default, but it is fine to set this explicitly:

```yaml
plugins:
  contacts:
    compute_residence_times: true
```

Then run:

```bash
polyzymd compare run contacts -f comparison.yaml
```

### Analyze one polymer type only

```yaml
plugins:
  contacts:
    polymer_selection: "chainID C and resname SBM"
    protein_selection: "protein"
```

For EGMA-only analysis, switch to `resname EGM`.

### Restrict analysis to a protein region

```yaml
plugins:
  contacts:
    polymer_selection: "chainID C"
    protein_selection: "protein and (resname TRP PHE TYR)"
```

For an active-site slice, use a residue range selection such as:

```yaml
protein_selection: "protein and (resid 75-80 or resid 130-140)"
```

### Run with reproducible cache behavior

```bash
# Use cache if present
polyzymd compare run contacts -f comparison.yaml --eq-time 10ns

# Ignore cache and recompute
polyzymd compare run contacts -f comparison.yaml --eq-time 10ns --recompute
```

### Use a fuller contacts configuration

If you want one place to set the most common contacts options:

```yaml
plugins:
  contacts:
    polymer_selection: "chainID C"
    protein_selection: "protein"
    cutoff: 4.5
    polymer_types: ["SBM", "EGM"]
    grouping: "aa_class"
    compute_residence_times: true
    fdr_alpha: 0.05
    min_effect_size: 0.5
    top_residues: 10
```

This is usually enough for cross-condition comparison without extra tuning.

### Add user-defined protein groups and partitions

Use this when you want plots and summaries for specific regions:

```yaml
plugins:
  contacts:
    protein_groups:
      active_site: [77, 133, 156]
      binding_patch: [45, 46, 47, 82, 84]
      distal_surface: [12, 13, 14, 190, 191, 192]
    protein_partitions:
      functional_regions: [active_site, binding_patch, distal_surface]
```

After running, these partitions are used in partition-level contacts plots.

### Generate contacts plots after running

```bash
polyzymd compare plot-all -f comparison.yaml
```

You will get contact-fraction and residence-time profiles plus grouped bar
plots for AA classes and (if configured) user partitions.

For the full list of plot outputs and plot settings, see
{doc}`../reference/analysis_contacts_reference`.

### Run only contacts in a multi-plugin config

If your `comparison.yaml` enables several plugins, you can run only contacts:

```bash
polyzymd compare run contacts -f comparison.yaml --eq-time 10ns
```

Later, run all enabled plugins:

```bash
polyzymd compare run-all -f comparison.yaml --eq-time 10ns
```

## Quick output checks

After a run, verify these two things first:

1. **Coverage and mean contact fraction** in the aggregated result
2. **Replicate count used** in aggregation

Minimal check pattern:

```bash
ls analysis/<condition>/contacts/aggregated/
```

Then inspect key values programmatically:

```python
import json
from pathlib import Path

agg = json.loads(Path("analysis/<condition>/contacts/aggregated/result.json").read_text())
print(f"n_replicates={agg['n_replicates']}")
print(f"coverage={agg['coverage_mean']:.3f} ± {agg['coverage_sem']:.3f}")
print(
    "mean_contact_fraction="
    f"{agg['mean_contact_fraction']:.3f} ± {agg['mean_contact_fraction_sem']:.3f}"
)
```

If residence times were enabled, also check:

```python
for ptype, stats in agg.get("residence_time_by_polymer_type", {}).items():
    print(f"{ptype}: mean={stats[0]:.2f} frames, sem={stats[1]:.2f}")
```

## Programmatic post-processing (JSON)

After CLI execution, load result files directly:

```python
import json
from pathlib import Path

replicate_result = json.loads(
    Path("analysis/<condition>/contacts/run_1/result.json").read_text()
)
print(f"Coverage: {replicate_result['coverage_fraction']:.1%}")

aggregated_result = json.loads(
    Path("analysis/<condition>/contacts/aggregated/result.json").read_text()
)
print(
    "Mean contact fraction: "
    f"{aggregated_result['mean_contact_fraction']:.1%} "
    f"± {aggregated_result['mean_contact_fraction_sem']:.1%}"
)
```

For worked Python recipes (interaction matrices, group-level summaries, custom
queries), use {doc}`analysis_contacts_cookbook`.

## Compare conditions

Use the same comparison workflow as other stable analyses:

```bash
polyzymd compare run contacts -f comparison.yaml --eq-time 10ns
```

The plugin compares conditions with dual primary metrics:

- coverage
- mean contact fraction

For multi-plugin comparison workflow details, see
{doc}`analysis_compare_conditions`.

## Reference and troubleshooting

For complete lookup documentation, including:

- full configuration field tables
- output directory structure and JSON schemas
- full plot catalog and `plot_settings.contacts` options
- common CLI options
- troubleshooting fixes

see {doc}`../reference/analysis_contacts_reference`.

## Next steps

- {doc}`analysis_compare_conditions`
- {doc}`analysis_contacts_cookbook`
- {doc}`analysis_rmsf_quickstart`
- {doc}`analysis_triad_quickstart`
- {doc}`../reference/analysis_contacts_reference`
