# Tutorial: Analyze a Study from Finished Simulations

This tutorial walks through one complete PolyzyMD analysis story:

- three simulation conditions already exist
- you compare the protein-polymer hydrogen bonds with one `polyzymd analyze` command
- you read the result and find the stored values and figures
- you add RMSF and polymer-protein contacts for the same conditions

By the end, you will have validated comparisons, stored per-replicate results
and figures for a small three-condition study.

## What You Will Learn

- How to compare conditions with `polyzymd analyze NAME -c ... -c ...`
- How to read the per-condition values and the comparisons against the control
- How to pick another result of the same analysis with `--run`
- What the output directory structure looks like after a successful run

## Prerequisites

Before starting, make sure you have:

- Completed production trajectories for at least three conditions (DCD format
  in PolyzyMD's standard directory layout)
- One `config.yaml` per condition
- The topology written by the build: `system.prmtop`, or `solvated_system.pdb`
  for older runs (see {doc}`../reference/data_requirements`)
- PolyzyMD installed in a pixi environment (see {doc}`../get_started/installation`)

If you have not run a single-condition analysis yet, complete
{doc}`first_analysis` first.

```{important}
**Resource requirements:** `polyzymd analyze` loads trajectories and can
require substantial RAM, CPU time, and scratch I/O. On shared HPC systems, run
it inside an allocated job or interactive compute session, not on a login
node; {doc}`../how_to/hpc_execution` shows a batch script.
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
conditions plus the results folder. -->

## Step 1: Compare the Hydrogen Bonds

From the study root, make a folder for the results and run the comparison
there. Protein-polymer hydrogen bonds exist only in the two conditions with a
polymer, so compare those two; the first `-c` is the control, and `--label`
names the conditions in the same order:

```bash
cd my_enzyme_study
mkdir polymer_stabilization_study
cd polymer_stabilization_study
pixi run -e analysis polyzymd analyze hydrogen_bonds \
  -c ../SBMA_100_enzyme_DMSO/config.yaml \
  -c ../EGMA_100_enzyme_DMSO/config.yaml \
  --label "100% SBMA" --label "100% EGMA" \
  --replicates 1-3 --eq 10ns --set d_a_cutoff=3.0
```

The command discards the first 10 ns of every replicate and, on every
production frame, counts the hydrogen bonds between the protein (`chainid A`)
and the polymer (`chainid C`) with MDAnalysis `HydrogenBondAnalysis`, a donor
within 3.0 Å of the acceptor and a donor-hydrogen-acceptor angle of at least
150°. It prints one line per condition with the mean number of hydrogen bonds
per frame over the replicates, its 95% interval and every replicate value,
one line comparing 100% EGMA with 100% SBMA by Welch's t test, and a
`verdict:` line. Every `warning:` line is part of the result.

## Step 2: Pick Another Result

`all_runs` in the JSON report (`--format json`) lists every result of the
analysis. Pick the per-residue occupancy, the fraction of frames in which
each protein residue has a hydrogen bond to the polymer:

```bash
pixi run -e analysis polyzymd analyze hydrogen_bonds \
  -c ../SBMA_100_enzyme_DMSO/config.yaml \
  -c ../EGMA_100_enzyme_DMSO/config.yaml \
  --label "100% SBMA" --label "100% EGMA" \
  --replicates 1-3 --eq 10ns --set d_a_cutoff=3.0 \
  --run protein_polymer_residues
```

The comparison is made at every residue, with the Benjamini-Hochberg
correction over all of them; see {doc}`../how_to/hydrogen_bonds` for every
result and its figures.

## Step 3: Check the Outputs

At this point you should have:

```text
polymer_stabilization_study/
├── polyzymd_results/
│   ├── hydrogen_bonds_protein_polymer/
│   │   ├── 100_SBMA/
│   │   │   ├── replicate_1/
│   │   │   │   ├── record.json
│   │   │   │   └── values.npz
│   │   │   └── ...
│   │   └── 100_EGMA/
│   └── residue_hbond_occupancy_protein_polymer/
│       └── ...
└── figures/
    └── hydrogen_bonds/
        ├── hbonds_protein_polymer_mean_hbonds_comparison.png
        ├── hbonds_protein_polymer_residues_profile.png
        └── hbonds_protein_polymer_residues_difference.png
```

Each folder name is the result or condition label with every run of
characters other than letters, digits, `.`, `+` and `-` replaced by `_`.
`record.json` holds what was measured and on which inputs, and `values.npz`
the replicate's values. A later run with the same settings reads them back
instead of loading the trajectories again.

<!-- IMAGE OPPORTUNITY: Add one example comparison figure here so the tutorial
has a visual payoff immediately before the final success state. -->

## Step 4: Add RMSF and Contacts for the Same Study

RMSF takes all three conditions, with the no-polymer condition as the
control:

```bash
pixi run -e analysis polyzymd analyze rmsf \
  -c ../noPoly_enzyme_DMSO/config.yaml \
  -c ../SBMA_100_enzyme_DMSO/config.yaml \
  -c ../EGMA_100_enzyme_DMSO/config.yaml \
  --label "No Polymer" --label "100% SBMA" --label "100% EGMA" \
  --replicates 1-3 --eq 10ns
```

The report compares each condition's core RMSF with the control, and the
figures go to `figures/rmsf/`. See {doc}`../how_to/analysis_rmsf_quickstart`
for the reference, the core and the per-residue comparison.

Polymer-protein contacts run the same way, over the two conditions that have a
polymer. The command below compares the fraction of protein residues in contact
with the polymer on at least one frame, where a residue is in contact when the
polymer occludes its solvent-accessible surface; `--run mean_lifetime` reports
how long contacts last instead:

```bash
pixi run -e analysis polyzymd analyze contacts \
  -c ../SBMA_100_enzyme_DMSO/config.yaml \
  -c ../EGMA_100_enzyme_DMSO/config.yaml \
  --label "100% SBMA" --label "100% EGMA" \
  --replicates 1-3 --eq 10ns
```

See {doc}`../how_to/analysis_contacts_quickstart` for the other results and the
distance method.

## What to Do Next

- Use [How to Compare Simulation Conditions](../how_to/analysis_compare_conditions.md) for
  a shorter operational version of this workflow
- Use [Get a validated number with one command](../how_to/analysis_agent_protocol.md)
  to read every line of the report
- Explore metric-specific guides:
  - [Run RMSF Analysis](../how_to/analysis_rmsf_quickstart.md)
  - [Run Contacts Analysis](../how_to/analysis_contacts_quickstart.md)
  - [Run Distance Analysis](../how_to/analysis_distances_quickstart.md)
  - [Measure a Catalytic Triad on the Analysis API](../how_to/analysis_triad_quickstart.md)
- For removed experimental analyses, see
  [Experimental analyses](../reference/experimental_analyses_archive.md); they
  are not active v1.3 workflows.
