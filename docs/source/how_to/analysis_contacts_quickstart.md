# Polymer-Protein Contacts Analysis: Quick Start

Analyze polymer-protein contact frequencies and coverage for one or more
conditions using the `contacts` plugin.

:::{admonition} Environment Setup
:class: tip

All analysis commands below assume you have activated the PolyzyMD analysis
pixi environment:

```bash
pixi shell -e analysis
```

Alternatively, prefix each command with `pixi run -e analysis`.
:::

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
    polymer_selection: "chainid C"
    protein_selection: "chainid A"
    cutoff: 4.5
```

### 2) Run contacts

```bash
polyzymd compare run contacts -f comparison.yaml --eq-time 10ns
```

Expected output includes per-replicate progress and one line per observable,
each naming its unit, its interval and the number of replicates behind it.

### 3) Run all enabled plugins (optional)

```bash
polyzymd compare run-all -f comparison.yaml --eq-time 10ns
```

## Key numbers to check first

- `contact_count`: residue pairs in contact per frame
- `coverage_per_frame`: share of the protein in contact per frame
- `coverage_any_frame`: share touched at any point in the window
- `contact_fraction`: which residues, as a profile over residue IDs
- `residence_time_distribution`: whether contacts are brief or long lived

## Common tasks

### Analyse one polymer type only

```yaml
plugins:
  contacts:
    protein_selection: "protein"
    polymer_selection: "resname SBM EGM"
    polymer_types: ["SBM"]
```

### Use the heavy-atom contact criterion

Hydrogens count toward the cutoff by default. To use the literature convention
instead, which lowers every contact count:

```yaml
plugins:
  contacts:
    heavy_atoms_only: true
```

### Change the residence-time bins

The default edges double from one frame of a 40 ps trajectory. Give your own
when your frames are spaced differently or your contacts are longer lived:

```yaml
plugins:
  contacts:
    residence_time_edges_ns: [0.0, 0.1, 0.5, 1.0, 5.0, 25.0]
```

Every condition of a comparison has to use the same edges, because the bins are
the index of a profile and the framework averages profiles element by element.

### Recompute after changing a setting

```bash
polyzymd compare run contacts -f comparison.yaml --eq-time 10ns --recompute
```

The framework reuses a replicate only when the polyzymd version, the plugin
source, the settings, the config and the input files all still match, so the
flag is rarely needed.

### Run only contacts in a multi-plugin config

```bash
polyzymd compare run contacts -f comparison.yaml --eq-time 10ns
```

Later, run all enabled plugins:

```bash
polyzymd compare run-all -f comparison.yaml --eq-time 10ns
```

## Read the results from Python

```python
import json
from pathlib import Path

aggregate = json.loads(
    Path("analysis/<condition>/contacts/aggregated/result.json").read_text()
)
for observable in aggregate["payload"]["observables"]:
    if observable["kind"] == "profile":
        continue
    print(
        f"{observable['name']}: {observable['mean']:.3g} {observable['unit']}"
        f" (n={observable['n_replicates']})"
    )
```

The per-residue profiles and the per-frame series are in
`run_<replicate>/observables.npz`, and the contact events are in
`run_<replicate>/sidecars/contact_events.npz` under the key `contact_events`.

## Compare conditions

```bash
polyzymd compare run contacts -f comparison.yaml --eq-time 10ns
```

Every scalar observable is tested against the control condition on
replicate-level values, and every test in the run shares one Benjamini-Hochberg
family. For multi-plugin comparison workflow details, see
{doc}`analysis_compare_conditions`.

## Reference and troubleshooting

For the settings tables, the observables, the output paths and the event table
layout, see {doc}`../reference/analysis_contacts_reference`.

## Next steps

- {doc}`analysis_compare_conditions`
- {doc}`analysis_rmsf_quickstart`
- {doc}`analysis_triad_quickstart`
- {doc}`../reference/analysis_contacts_reference`
