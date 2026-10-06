# Tutorial: Run Your First Analysis

This tutorial walks you from finished trajectory files to your first analysis
result. You will run the RMSF analysis on a single simulation condition with
`polyzymd analyze`, read the report, and see where the results and figures end
up on disk.

## What You Will Learn

- How to run the RMSF analysis with `polyzymd analyze rmsf`
- How to read the report
- Where the stored results and figures are, and how they are reused

## Prerequisites

Before starting, make sure you have:

- A completed production simulation with at least 1 replicate
- The `config.yaml` file from that simulation
- Trajectory files in the expected directory layout (see
  {doc}`../reference/data_requirements`)
- PolyzyMD installed in a pixi environment (see {doc}`../get_started/installation`)

If you have not run a simulation yet, complete
{doc}`../get_started/quickstart` first.

```{important}
**Resource requirements:** `polyzymd analyze` loads trajectories, which can
require substantial RAM, CPU time and scratch I/O. On shared HPC systems, run it
inside an allocated job or interactive compute session, not on a login node.
```

## Step 1: Run the RMSF Analysis

From the directory where you want the results, run:

```bash
pixi run -e analysis polyzymd analyze rmsf \
  -c /path/to/my_simulation/config.yaml --label "My Simulation" --eq 10ns
```

- **`-c`** points to the simulation's `config.yaml`. This is how PolyzyMD finds
  the topology and the production trajectory of every replicate on disk.
- **`--label`** names the condition in the report. Without it, the condition is
  named after the folder holding the config.
- **`--eq`** is the equilibration window: the time at the start of each
  replicate's production trajectory that is left out. Adjust it to your system.

By default the analysis measures the Cα atoms of the protein
(`protein and name CA`), superposes every production frame on the replicate's
most representative frame, and gives each residue's RMSF: how much it
fluctuates about its mean position.

## Step 2: Read the Report

You should see output similar to:

```text
# polyzymd analyze rmsf  metric core_rmsf  unit A  run core_rmsf  eq 10ns  conditions 1  replicates 1  protocol rmsf/2
My Simulation  n 1  mean 0.8214  sem na  ci95 na  values 0.8214
warning: condition My Simulation has one replicate, so it has no interval
verdict: My Simulation core_rmsf 0.8214 A (no interval, n 1)
```

- The header names the analysis, the reported result (`core_rmsf`), its unit,
  the equilibration window and the number of replicates.
- The condition line gives the number of replicates `n`, the mean, the standard
  error, the 95 percent interval and every replicate value.
- `core_rmsf` combines the residues into one number per replicate: the root of
  their mean square fluctuation.

With one replicate there is no standard error or interval, so `sem` and `ci95`
read `na`. This tutorial uses one replicate so you can complete the workflow
quickly, but uncertainty and comparisons need at least 2 replicates per
condition.

To see the value of every residue instead, ask for the profile:

```bash
pixi run -e analysis polyzymd analyze rmsf \
  -c /path/to/my_simulation/config.yaml --label "My Simulation" --eq 10ns --run rmsf
```

This prints one line per residue, starting with the residue ID. Add
`--format json` to any run for the full report, with every field documented in
{doc}`../reference/analysis_protocol_report`.

```{tip}
If you see an error about a missing working directory or trajectory, check the
`config` path and that your trajectory files exist on disk. See
{doc}`../how_to/troubleshooting` for common fixes.
```

## Step 3: Find Your Results

After the run, the directory holds:

```text
.
├── polyzymd_results/
│   └── rms_decomposition/
│       └── My_Simulation/
│           └── replicate_1/
│               ├── record.json    # what was measured, with which inputs and settings
│               └── values.npz     # the per-residue values of this replicate
└── figures/
    └── rmsf/
        ├── rmsf_profile.png
        ├── offset_profile.png
        ├── rmsd_per_residue_profile.png
        ├── rms_decomposition.png
        └── rmsf_comparison.png
```

- **`record.json`** holds the function and a hash of its source, its
  arguments, the config hash, the path, size and modification time of the
  topology and every trajectory file, the frames and times used, the residue
  labels and the software versions, so the value can be traced back to what
  produced it.
- **`values.npz`** holds the replicate's per-residue values.
- The **figures** show each residue's RMSF, its offset from the reference and
  its RMS deviation from the reference (`rmsd_per_residue`), and the three core values.

The second run above, with `--run rmsf`, did not read the trajectory again: a
stored result is reused when the inputs and settings match. Pass `--recompute`
to measure again anyway, and `--no-plots` to skip the figures.

## What's Next

Now that you have run one analysis on one condition, here are some natural next
steps:

- {doc}`../how_to/analysis_rmsf_quickstart` --- compare conditions, choose the
  reference, and define the core and regions
- {doc}`../how_to/analysis_compare_conditions` --- compare several conditions
- {doc}`../how_to/study_api` --- run your own function on every
  replicate from Python
- {doc}`../reference/data_requirements` --- directory layout reference and
  path resolution rules
