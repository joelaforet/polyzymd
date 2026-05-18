# Tutorial: Run Your First Analysis

This tutorial walks you from finished trajectory files to your first analysis
result. You will run RMSF analysis on a single simulation condition using the
comparison pipeline, and see where the results end up on disk.

## What You Will Learn

- How to create a comparison project for a single condition
- How to run RMSF analysis using `polyzymd compare run`
- How to read the output and find result files

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
**Resource requirements:** `polyzymd compare init`, `validate`, `status`, and
`--help` are lightweight. Commands that load trajectories, such as
`polyzymd compare run` and `run-all`, can require substantial RAM, CPU/GPU time,
and scratch I/O. On shared HPC systems, run them inside an allocated job or
interactive compute session, not on a login node. If an analysis is killed or
runs out of memory, request more resources or use `polyzymd compare submit`.
```

## Step 1: Create a Comparison Project

From the directory where you keep your simulation projects, run:

```bash
pixi run -e build polyzymd compare init -n my_first_analysis
cd my_first_analysis
```

This creates a small project scaffold:

```text
my_first_analysis/
├── comparison.yaml    # Analysis configuration (you will edit this)
├── comparison/        # Cross-condition comparison outputs
├── figures/           # Where plots are saved
└── structures/        # Optional shared structure files
```

The generated `comparison.yaml` is a template with placeholder values. You
will replace them in the next step. Per-replicate and per-condition analysis
artifacts are created under `analysis/` when you run an analysis.

## Step 2: Edit comparison.yaml

Open `comparison.yaml` in your editor and replace the contents with a minimal
single-condition configuration:

```yaml
name: "my_first_analysis"
description: "First analysis run"
control: null

conditions:
  - label: "My Simulation"
    config: "/path/to/my_simulation/config.yaml"
    replicates: [1]

defaults:
  equilibration_time: "10ns"

plugins:
  rmsf:
    selection: "protein and name CA"
```

Here is what each section does:

- **`name`** and **`description`** --- identify this comparison project.
- **`control`** --- the label of the control condition for statistical tests.
  Set to `null` when you only have one condition.
- **`conditions`** --- a list of simulation conditions to analyze. Each entry
  needs a `label`, a path to that simulation's `config.yaml`, and which
  `replicates` to include.
- **`defaults.equilibration_time`** --- how much time at the start of each
  trajectory to discard before analysis. Adjust to match your system's
  equilibration period.
- **`plugins.rmsf`** --- settings for the RMSF analysis plugin. The
  `selection` field is an MDAnalysis atom selection string.

```{important}
The `config` path must point to the simulation project's `config.yaml`. This
is how PolyzyMD locates your topology and trajectory files on disk. Relative
paths are resolved from the directory containing `comparison.yaml`.
```

For the full list of configuration fields, see
{doc}`../reference/analysis_comparison_reference`.

## Step 3: Run RMSF Analysis

Run the analysis with:

```bash
pixi run -e build polyzymd compare run rmsf -f comparison.yaml --eq-time 10ns
```

```{note}
The `--eq-time` flag overrides `defaults.equilibration_time` from your YAML
file. If you omit `--eq-time`, the value from `comparison.yaml` is used. This
is handy for quickly testing different equilibration cutoffs without editing the
YAML each time.
```

```{tip}
**On an HPC cluster?** Use `polyzymd compare submit` instead of `compare run`
to dispatch analysis as SLURM jobs. This is especially important for expensive
analyses (SASA, contacts, hydrogen bonds) on large studies. See
{doc}`../how_to/hpc_execution` for the full workflow.
```

You should see output similar to:

```text
Comparison: my_first_analysis
Plugin: rmsf
Conditions: 1
Equilibration: 10ns

[My Simulation] Computing replicate 1...
  Loading trajectory (skipping first 10 ns)...
  RMSF computed (142 residues, 490 frames)
[My Simulation] Aggregating 1 replicate...

RMSF Comparison Complete
  My Simulation: mean RMSF = 0.621 Å
  SEM: n/a (single replicate)
  Statistical comparisons: not testable until each condition has at least 2 replicates
```

```{tip}
If you see `RMSF Comparison Complete` with a mean value, the analysis succeeded.
If you see an error about a missing working directory or trajectory, check that
the `config` path in `comparison.yaml` is correct and that your trajectory
files exist on disk. See {doc}`../how_to/troubleshooting` for common fixes.
```

## Step 4: Find Your Results

After the run completes, your project directory looks like this:

```text
my_first_analysis/
├── comparison.yaml
├── analysis/
│   └── My_Simulation/
│       └── rmsf/
│           ├── run_1/
│           │   └── result.json                # ReplicateArtifact
│           └── aggregated/
│               └── result.json                # ConditionArtifact
├── comparison/
│   └── rmsf/
│       └── result.json                        # Comparison artifact/summary
├── figures/
└── structures/
```

The key files are:

- **`analysis/My_Simulation/rmsf/run_1/result.json`** --- a
  `ReplicateArtifact` for replicate 1, including per-replicate RMSF values for
  every residue in the selection after discarding the first 10 ns.
- **`analysis/My_Simulation/rmsf/aggregated/result.json`** --- a
  `ConditionArtifact` with statistics across replicates. With one replicate,
  the values closely mirror the replicate artifact.
- **`comparison/rmsf/result.json`** --- the comparison artifact/summary with mean
  RMSF and ranking information. Singleton SEM is unavailable and suppressed
  until at least 2 replicates contribute to a condition.

```{note}
PolyzyMD artifacts are more than bare result values. They also carry metadata,
provenance, warnings, and references to sidecar files when an analysis needs to
store larger tables or arrays outside the main JSON file.
```

```{note}
This tutorial intentionally uses one replicate so you can complete the workflow
quickly. RMSF supports this smoke-test mode, but a single replicate cannot
estimate between-replicate uncertainty. SEM is unavailable, and condition-level
statistical comparisons are not testable until each condition has at least 2
replicates.
```

## Step 5: Add Plotting (Optional)

To generate figures alongside the analysis, re-run with the `--plot` flag on
the `run-all` command:

```bash
pixi run -e build polyzymd compare run-all -f comparison.yaml --eq-time 10ns --plot
```

Or generate plots separately after the analysis has already been cached:

```bash
pixi run -e build polyzymd compare plot-all -f comparison.yaml
```

Figures are saved to the `figures/` directory:

```text
figures/
└── rmsf/
    ├── rmsf_comparison.png
    └── rmsf_profile.png
```

With a single condition, the comparison chart is a simple one-bar summary and
the profile plot shows per-residue RMSF values. Error bars and statistical
comparisons become more useful when you add additional conditions and
replicates.

## What's Next

Now that you have run one analysis on one condition, here are some natural next
steps:

- {doc}`../how_to/analysis_compare_conditions` --- Add a second condition
  and run a statistical comparison
- {doc}`analysis_complete_workflow` --- Full multi-condition workflow with
  multiple analysis types
- {doc}`../how_to/analysis_rmsf_quickstart` --- RMSF-specific options
  (reference modes, selections, troubleshooting)
- {doc}`../reference/data_requirements` --- Directory layout reference and
  path resolution rules
