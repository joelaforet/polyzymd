# Distance Analysis: Quick Start

Compute inter-atomic distances with proper statistical handling in under
5 minutes.

```{note}
**Want to understand the statistics?** This guide focuses on getting results
quickly. For proper uncertainty quantification (autocorrelation correction,
SEM vs. SD), see
{doc}`../explanation/analysis_statistics_best_practices`.
```

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
# Configure distance pairs in comparison.yaml, then run:
polyzymd compare run distances -f comparison.yaml --eq-time 10ns

# Run all enabled analyses in the same workflow
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Force recompute and machine-readable output
polyzymd compare run distances -f comparison.yaml --eq-time 10ns --recompute --format json
```

## Prerequisites

Before running distance analysis, you need:

1. **Completed production simulation(s)** — at least one replicate
2. **Comparison config** — `comparison.yaml` with conditions and plugin settings
3. **Trajectory files** — in the scratch directory specified in config

Verify your setup:

```bash
# Check that trajectories exist
ls $(polyzymd info -c config.yaml --scratch-dir)/production_*/
```

## What distance analysis provides

For each configured pair, the plugin reports:

| Observable | Description |
|------------|-------------|
| Mean distance | Distance in angstrom, averaged over the production window of one replicate |
| Contact fraction | Fraction of frames in which the pair is below its threshold |

Both are reported per condition as a mean over replicates with its SEM and a
95 percent Student t interval, and are tested against the control condition on
replicate-level values.

```{tip}
**When to use distances vs. contacts vs. triad:**
- **Distances**: Specific atom pairs with continuous distance values
- **Contacts**: All residue-residue contacts at an interface (binary count)
- **Triad**: Pre-defined catalytic geometry with simultaneous contact analysis
```

## Basic Usage

`````{tab-set}

````{tab-item} YAML (Recommended)
Define distance pairs in `comparison.yaml`:

```yaml
# comparison.yaml
name: "distance_quickstart"
control: "no_polymer"

conditions:
  - label: "no_polymer"
    config: "configs/no_polymer.yaml"
    replicates: [1, 2, 3]
  - label: "with_polymer"
    config: "configs/with_polymer.yaml"
    replicates: [1, 2, 3]

plugins:
  distances:
    pairs:
      - label: "Ser77-His156"
        selection_a: "protein and resid 77 and name OG"
        selection_b: "protein and resid 156 and name NE2"
      - label: "His156-Asp133"
        selection_a: "protein and resid 156 and name ND1"
        selection_b: "midpoint(protein and resid 133 and name OD1 OD2)"
```

Run analysis:

```bash
# Run distances only
polyzymd compare run distances -f comparison.yaml --eq-time 10ns

# Run all enabled analyses
polyzymd compare run-all -f comparison.yaml --eq-time 10ns

# Force recompute
polyzymd compare run distances -f comparison.yaml --eq-time 10ns --recompute
```
````

````{tab-item} CLI
```bash
polyzymd compare run distances -f comparison.yaml --eq-time 10ns
```

Expected behavior:

```text
Loading comparison config from: comparison.yaml
Running plugin: distances
  Equilibration: 10ns
  Conditions: no_polymer, with_polymer
Distance comparison complete
```

Run all enabled analyses:

```bash
polyzymd compare run-all -f comparison.yaml --eq-time 10ns
```
````

`````

## Add Contact-Style Thresholds

Set a threshold to report the fraction of frames below a cutoff (useful for
hydrogen-bond-style geometry checks).

```yaml
plugins:
  distances:
    threshold: 3.5
    pairs:
      - label: "Ser77-His156"
        selection_a: "protein and resid 77 and name OG"
        selection_b: "protein and resid 156 and name NE2"
```

```bash
polyzymd compare run distances -f comparison.yaml --eq-time 10ns
```

## Write Robust Selections

PolyzyMD supports standard MDAnalysis selections plus helper syntax like
`midpoint(...)`, `com(...)`, and `pdbindex N`.

```{warning}
**Chain-aware selections are required**

Residue numbers restart by chain in PolyzyMD systems. A selection like
`resid 141-148` can match multiple chains.

For protein residues, include `protein and ...`:

```yaml
# Incorrect
selection_a: "com(resid 141-148)"

# Correct
selection_a: "com(protein and resid 141-148)"
```
```

Common patterns:

```yaml
# Midpoint of Asp carboxylate oxygens
selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"

# Center of mass of ligand
selection_b: "com(resname LIG)"

# Single atom
selection_a: "protein and resid 77 and name OG"

# Atom by PDB serial number
selection_a: "pdbindex 2740"
```

## Keep PBC on

Distance analysis applies the minimum image convention by default
(`use_pbc: true`), using the box stored in each frame. Leave it on unless you
know your trajectory is already unwrapped, because it is what keeps a pair from
being measured the long way around the box.

Distances are not aligned. `align_trajectory` and the `alignment_*` fields are
accepted for one release, ignored, and raise a `DeprecationWarning`. A distance does not
change when the whole system is rotated or translated, so alignment could only
cost time, and the in-memory alignment PolyzyMD used to run rotated the
coordinates without rotating the box, which broke the minimum image convention
for long pairs.

To override PBC, see {doc}`../reference/analysis_distances_reference`.

## Compare Distances Across Conditions

Use the same command after defining conditions and pairs in `comparison.yaml`:

```bash
polyzymd compare run distances -f comparison.yaml --eq-time 10ns
```

This provides:

- A per-pair mean distance and contact fraction for every condition, each with
  its SEM, its interval and its replicate count
- A test of every pair against the control, with Cohen's d and a
  Benjamini-Hochberg adjusted p-value over the whole run

For broader multi-plugin workflows, see
{doc}`analysis_compare_conditions`.

## Reference and Troubleshooting

For full field tables, output JSON schemas, plot types, CLI options, and
troubleshooting fixes, see
{doc}`../reference/analysis_distances_reference`.

## Next Steps

- Compare multiple analyses: {doc}`analysis_compare_conditions`
- Catalytic triad workflow: {doc}`analysis_triad_quickstart`
- Statistical interpretation: {doc}`../explanation/analysis_statistics_best_practices`
- Contacts workflow: {doc}`analysis_contacts_quickstart`
