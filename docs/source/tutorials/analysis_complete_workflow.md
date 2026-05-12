# Tutorial: Analyze a Study from Finished Simulations

This tutorial walks through one complete PolyzyMD analysis story:

- three simulation conditions already exist
- you create one `comparison.yaml`
- you compare the conditions
- you finish with `polyzymd compare plot-all` as the smoke test

By the end, you will have a working comparison workspace with JSON results and
figures for a small three-condition study.

## What You Will Learn

- How to initialize a comparison workspace with `polyzymd compare init`
- How to write a `comparison.yaml` that defines conditions and analysis plugins
- How to run cross-condition comparisons and generate figures
- What the output directory structure looks like after a successful run

## Prerequisites

Before starting, make sure you have:

- Completed production trajectories for at least three conditions (DCD format
  in PolyzyMD's standard directory layout)
- One `config.yaml` per condition
- A topology such as `solvated_system.pdb` already produced during the build
- PolyzyMD installed in a pixi environment (see {doc}`installation`)

If you have not run a single-condition analysis yet, complete
{doc}`first_analysis` first.

```{important}
This tutorial uses the stable `v1.3.0` comparison stack: RMSD, Rg, RMSF,
contacts, distances, catalytic triad, secondary structure, SASA, and hydrogen
bonds. Experimental workflows are linked at the end, but they are not part of
the main tutorial path.
```

```{important}
**Resource requirements:** Workspace setup and validation commands are
lightweight. Commands that load trajectories, such as `polyzymd compare run`,
`run-all`, and plotting over large cached results, can require substantial RAM,
CPU/GPU time, and scratch I/O. On shared HPC systems, run them inside an
allocated job or interactive compute session, not on a login node. If a command
is killed or runs out of memory, request more resources or use
`polyzymd compare submit`.
```

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

## Step 1: Create the Comparison Workspace

From the study root, initialize a comparison project and move into it:

```bash
cd my_enzyme_study
pixi run -e build polyzymd compare init -n polymer_stabilization_study
cd polymer_stabilization_study
```

Now edit `comparison.yaml` to point at the three conditions and define the
analysis settings:

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
    polymer_selection: "chainid C"
    protein_selection: "chainid A"
    cutoff: 4.5
    compute_residence_times: true
```

## Step 2: Validate the Comparison Config

```bash
pixi run -e build polyzymd compare validate
```

You should see a passing summary that lists the three conditions and the
enabled analyses.

## Step 3: Run the Cross-Condition Comparison

For the tutorial, use the batch runner:

```bash
pixi run -e build polyzymd compare run-all
```

This runs every enabled comparison and writes canonical cache files into
`comparison/<analysis>/result.json`.

```{tip}
**On an HPC cluster?** For large studies, submit each analysis as a SLURM
job DAG instead of running interactively:

    pixi run -e build polyzymd compare submit sasa --partition <part> --mem 8G

This parallelizes across replicates and conditions. See
{doc}`../how_to/hpc_execution` for the complete HPC workflow.
```

If you prefer to inspect one comparison first, a good sanity check is:

```bash
pixi run -e build polyzymd compare run rmsf
```

## Step 4: Generate the Figures

Now run the plotting smoke test:

```bash
pixi run -e build polyzymd compare plot-all --list-available
pixi run -e build polyzymd compare plot-all
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

- Use [How to Compare Simulation Conditions](../how_to/analysis_compare_conditions.md) for
  a shorter operational version of this workflow
- Use [Comparison and Plotting Reference](../reference/analysis_comparison_reference.md)
  for CLI, config, and output lookup
- Explore metric-specific guides:
  - [Run RMSF Analysis](../how_to/analysis_rmsf_quickstart.md)
  - [Run Contacts Analysis](../how_to/analysis_contacts_quickstart.md)
  - [Run Distance Analysis](../how_to/analysis_distances_quickstart.md)
  - [Run Catalytic Triad Analysis](../how_to/analysis_triad_quickstart.md)
- For removed experimental analyses, see
  [Experimental analyses](../reference/experimental_analyses_archive.md); they
  are not active v1.3 workflows.
