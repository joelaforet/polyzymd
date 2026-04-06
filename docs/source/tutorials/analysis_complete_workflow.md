# Tutorial: Analyze a Study from Finished Simulations

:::{admonition} Modern workflow available
:class: note

This tutorial documents the manual, condition-by-condition analysis workflow.
For most users, the **comparison pipeline** is now the recommended approach —
it handles compute, aggregation, comparison, and plotting automatically:

```bash
polyzymd compare init -n my_study
# Edit comparison.yaml with your conditions
polyzymd compare run-all --plot
```

See {doc}`analysis_compare_conditions` for the full guide.
:::

This tutorial walks through one complete PolyzyMD analysis story:

- three simulation conditions already exist
- you create one `analysis.yaml`
- you run per-condition analyses
- you compare the conditions
- you finish with `polyzymd compare plot-all` as the smoke test

By the end, you will have a working comparison workspace with JSON results and
figures for a small three-condition study.

```{important}
This tutorial uses the stable `v1.2.0` comparison stack: RMSF, contacts,
distances, and catalytic triad. Experimental workflows are linked at the end,
but they are not part of the main tutorial path.
```

## Before You Start

You need:

- completed **OpenMM** production trajectories for each condition (DCD format
  in PolyzyMD's standard directory layout; GROMACS XTC support is planned for
  [v1.2.1](https://github.com/joelaforet/polyzymd/issues/47))
- one `config.yaml` per condition
- a topology such as `solvated_system.pdb` already produced during the build

## The Study We Will Analyze

We will assume a project laid out like this:

```text
my_enzyme_study/
├── noPoly_enzyme_DMSO/
│   ├── config.yaml
│   └── scratch/
├── SBMA_100_enzyme_DMSO/
│   ├── config.yaml
│   └── scratch/
└── EGMA_100_enzyme_DMSO/
    ├── config.yaml
    └── scratch/
```

The `scratch/` directories may be symlinks to large trajectory storage on your
cluster. PolyzyMD resolves those paths through each condition's `config.yaml`.

<!-- IMAGE OPPORTUNITY: Add a campaign directory-tree diagram showing the three
conditions plus the later comparison workspace. -->

## Step 1: Create `analysis.yaml` for the First Condition

Start in the control condition:

```bash
cd my_enzyme_study/noPoly_enzyme_DMSO
polyzymd analyze init
```

Edit the generated `analysis.yaml` so it enables a small stable analysis set:

```yaml
replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"

rmsf:
  enabled: true
  selection: "protein and name CA"
  reference_mode: "average"

catalytic_triad:
  enabled: true
  name: "Ser-His-Asp"
  threshold: 3.5
  pairs:
    - label: "Ser77-His156"
      selection_a: "protein and resid 77 and name OG"
      selection_b: "protein and resid 156 and name NE2"
    - label: "His156-Asp133"
      selection_a: "protein and resid 156 and name ND1"
      selection_b: "midpoint(protein and resid 133 and name OD1 OD2)"

distances:
  enabled: true
  pairs:
    - label: "Substrate-Ser77"
      selection_a: "resname SUB and name C1"
      selection_b: "protein and resid 77 and name OG"

contacts:
  enabled: true
  polymer_selection: "chainID C"
  protein_selection: "protein"
  cutoff: 4.5
  compute_residence_times: true
```

## Step 2: Reuse the Same `analysis.yaml` for the Other Conditions

From the study root:

```bash
cd ../
cp noPoly_enzyme_DMSO/analysis.yaml SBMA_100_enzyme_DMSO/
cp noPoly_enzyme_DMSO/analysis.yaml EGMA_100_enzyme_DMSO/
```

This works because the condition-specific trajectory paths come from each
condition's own `config.yaml`.

## Step 3: Run Per-Condition Analyses

Run the same command in each condition directory:

```bash
cd noPoly_enzyme_DMSO
polyzymd analyze run

cd ../SBMA_100_enzyme_DMSO
polyzymd analyze run

cd ../EGMA_100_enzyme_DMSO
polyzymd analyze run
```

After each run, expect an `analysis/` directory with subdirectories such as:

```text
analysis/
├── rmsf/
├── catalytic_triad/
├── distances/
└── contacts/
```

For the no-polymer control, contacts may be skipped automatically because there
are no polymer atoms to analyze.

## Step 4: Create the Comparison Workspace

Go back to the study root and initialize a comparison project:

```bash
cd ..
polyzymd compare init -n polymer_stabilization_study
cd polymer_stabilization_study
```

Now edit `comparison.yaml` to point at the three conditions:

```yaml
name: "polymer_stabilization_study"
description: "Effect of SBMA vs EGMA polymer conjugation on enzyme stability"
control: "No Polymer"

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

defaults:
  equilibration_time: "10ns"

plugins:
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

  distances:
    pairs:
      - label: "Substrate-Ser77"
        selection_a: "resname SUB and name C1"
        selection_b: "protein and resid 77 and name OG"

  contacts:
    polymer_selection: "chainID C"
    protein_selection: "protein"
    cutoff: 4.5
    compute_residence_times: true
```

## Step 5: Validate the Comparison Config

```bash
polyzymd compare validate
```

You should see a passing summary that lists the three conditions and the
enabled analyses.

## Step 6: Run the Cross-Condition Comparison

For the tutorial, use the batch runner:

```bash
polyzymd compare run-all
```

This runs every enabled comparison and writes canonical cache files into
`comparison/<analysis>/result.json`.

If you prefer to inspect one comparison first, a good sanity check is:

```bash
polyzymd compare run rmsf
```

## Step 7: Generate the Figures

Now run the plotting smoke test:

```bash
polyzymd compare plot-all --list-available
polyzymd compare plot-all
```

If those commands succeed, your comparison workspace is in good shape.

<!-- IMAGE OPPORTUNITY: Add one example comparison figure here so the tutorial
has a visual payoff immediately before the final success state. -->

## What Success Looks Like

At this point you should have:

```text
polymer_stabilization_study/
├── comparison.yaml
├── comparison/
│   ├── rmsf/
│   │   └── result.json
│   ├── contacts/
│   │   └── result.json
│   ├── distances/
│   │   └── result.json
│   └── catalytic_triad/
│       └── result.json
└── figures/
    ├── rmsf_comparison.png
    ├── rmsf_profile.png
    ├── triad_kde_panel.png
    └── ...
```

That is the tutorial success state: the canonical comparison caches exist, the
figures exist, and `polyzymd compare plot-all` completes without error.

## What to Do Next

- Use [How to Compare Simulation Conditions](analysis_compare_conditions.md) for
  a shorter operational version of this workflow
- Use [Comparison and Plotting Reference](../reference/analysis_comparison_reference.md)
  for CLI, config, and output lookup
- Explore metric-specific guides:
  - [Run RMSF Analysis](analysis_rmsf_quickstart.md)
  - [Run Contacts Analysis](analysis_contacts_quickstart.md)
  - [Run Distance Analysis](analysis_distances_quickstart.md)
  - [Run Catalytic Triad Analysis](analysis_triad_quickstart.md)

Experimental workflows remain available, but they are intentionally outside the
main tutorial path for this release.
