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
- PolyzyMD installed in a pixi environment (see {doc}`installation`)

If you have not run a simulation yet, complete
{doc}`quickstart` first.

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
├── comparison/        # Where result JSON files are written
├── figures/           # Where plots are saved
└── structures/        # Optional shared structure files
```

The generated `comparison.yaml` is a template with placeholder values. You
will replace them in the next step.

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

```{tip}
If you are working with older PolyzyMD projects, you may see
`analysis_settings:` instead of `plugins:` in existing YAML files. Both keys
are accepted — `analysis_settings` is a legacy alias that still works but emits
a deprecation warning.
```

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
  My Simulation: mean RMSF = 0.621 ± 0.015 Å
```

```{tip}
If you see `RMSF Analysis Complete` with a mean value, the analysis succeeded.
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
│           │   └── rmsf_eq10ns.json           # Per-replicate result
│           └── aggregated/
│               └── result.json               # Combined result
├── comparison/
│   └── rmsf/
│       └── result.json                       # Comparison summary
├── figures/
└── structures/
```

The key files are:

- **`rmsf_eq10ns.json`** --- per-replicate RMSF values for every residue in
  the selection, computed after discarding the first 10 ns.
- **`aggregated/result.json`** --- aggregated statistics across replicates (with
  one replicate, this matches the per-replicate file).
- **`comparison/rmsf/result.json`** --- the comparison-level summary with mean
  RMSF, standard error, and ranking information.

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
    └── rmsf_profile.png
```

With a single condition, the profile plot shows per-residue RMSF values.
Comparison bar charts appear when you add a second condition.

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
