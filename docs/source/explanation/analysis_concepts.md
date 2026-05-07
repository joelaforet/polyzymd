# Analysis System Concepts

Your simulations are done. You have trajectories on disk and you want to
measure something — RMSF, contacts, distances, whatever. This page explains how
PolyzyMD's analysis system is put together so that when you run a command or
read an output file, you know what happened and where to look.

## The analysis pipeline

Every analysis in PolyzyMD follows the same four-stage pipeline:

```text
replicate stage  →  aggregate  →  compare  →  plot
```

Here is what each stage does:

| Stage | Scope | What it produces |
|-------|-------|------------------|
| **replicate stage** | One replicate of one condition | Raw per-replicate result (e.g., RMSF values for each residue) |
| **aggregate** | All replicates of one condition | Means, standard errors, and distributions across replicates |
| **compare** | All conditions together | Statistical tests, effect sizes, and rankings |
| **plot** | All conditions together | Figures saved as PNG files |

You don't call these stages yourself. When you run `polyzymd compare run`, the
CLI walks through the pipeline automatically. But knowing the stages helps when
you need to debug ("Which stage failed?") or interpret output ("Is this a
per-replicate file or an aggregated file?").

## `comparison.yaml` — the control file

The `comparison.yaml` file is the single input that defines an analysis run. It
tells PolyzyMD what simulations to analyze, what to measure, and how to compare
the results.

Here is a minimal example:

```yaml
name: "polymer_stability_study"

conditions:
  - label: "No Polymer"
    config: "../no_polymer/config.yaml"
    replicates: [1, 2, 3]
  - label: "100% SBMA"
    config: "../sbma_100/config.yaml"
    replicates: [1, 2, 3]

control: "No Polymer"

defaults:
  equilibration_time: "10ns"

plugins:
  rmsf:
    selection: "protein and name CA"
  contacts: {}
```

The key sections are:

### `conditions`

Each entry points to a simulation's `config.yaml` and lists which replicate
numbers to include. The `label` is a human-readable name that shows up in plots
and result files.

### `control`

Which condition to use as the baseline for statistical comparisons. Set this to
the label of your reference condition (typically an unmodified or no-polymer
system). If you only have one condition or don't want relative comparisons, set
it to `null` or leave it out.

### `defaults.equilibration_time`

How much trajectory to discard from the beginning of each run. Early frames
are typically not equilibrated, so the pipeline skips them. Specify as a string
with units: `"10ns"`, `"5000ps"`, etc. The default is `"10ns"`.

### `plugins`

Which analyses to run and their settings. Each key is a plugin name (like
`rmsf` or `contacts`), and the value is a settings block for that plugin. An
empty block `{}` means "run with defaults." Only plugins listed here are
executed — if you don't include `sasa`, SASA won't be computed.

For the complete schema with all fields, see
{doc}`../reference/comparison_yaml`.

## Conditions and replicates

These two terms come up everywhere in the analysis output, so it helps to be
precise about what they mean in PolyzyMD.

A **condition** is one simulation setup. Examples: "No Polymer", "SBMA-100",
"PEG-50". Each condition has its own `config.yaml` that defines the system
(which protein, which polymer, which force field, etc.).

A **replicate** is an independent run of the same condition. Replicates are
identified by number — 1, 2, 3, and so on. Each replicate uses the same
`config.yaml` but starts from a different random seed or initial configuration,
producing an independent trajectory.

The pipeline processes data in this order:

1. **Per-replicate**: the runner-backed replicate stage runs once for each
   replicate of each condition. If you have 2 conditions with 3 replicates
   each, that is 6 replicate-stage calls.
2. **Per-condition**: `aggregate` runs once per condition, combining the
   replicate results. That's 2 aggregate calls.
3. **Cross-condition**: `compare` and `plot` each run once, looking at all
   conditions together.

## Plugins — the analysis modules

PolyzyMD ships with 13 analysis plugins. Each plugin is a self-contained
module that knows how to compute one type of measurement, aggregate it, compare
across conditions, and generate plots.

The available plugins are:

| Plugin name | What it measures |
|-------------|-----------------|
| `rmsd` | Root-mean-square deviation over time |
| `rg` | Radius of gyration over time |
| `rmsf` | Root-mean-square fluctuation per residue |
| `contacts` | Intermolecular contacts between protein and other components |
| `distances` | Distances between specified atom groups |
| `catalytic_triad` | Catalytic triad geometry (active-site distances) |
| `secondary_structure` | Secondary structure content (helix, sheet, coil fractions) |
| `sasa` | Solvent-accessible surface area |
| `hydrogen_bonds` | Hydrogen bond occupancy and lifetimes |
| `polymer_bridging` | Polymer bridging topology between protein regions (experimental) |

Each plugin has a `Settings` model with configurable parameters. Most
parameters have sensible defaults, so you often just need `plugin_name: {}` in
your `comparison.yaml` to get started.

You configure plugins in the `plugins:` block. For example, to run RMSF with a
custom selection and contacts with defaults:

```yaml
plugins:
  rmsf:
    selection: "protein and name CA"
  contacts: {}
```

## Statistical comparison

When you have two or more conditions, the compare stage produces statistical
output so you can assess whether differences are meaningful. Here is what
PolyzyMD computes:

- **Pairwise t-tests** between each pair of conditions, with
  Benjamini–Hochberg FDR correction to account for multiple comparisons.
- **Effect sizes** (Cohen's d) for each pair, so you can see not just whether
  a difference is significant but how large it is.
- **ANOVA** when there are three or more conditions, as an omnibus test before
  the pairwise comparisons.
- **Rankings** of conditions from best to worst on each metric (where "best"
  depends on the metric — lower RMSF is better, higher helix fraction is
  better).

The comparison results are saved as JSON and also printed to the terminal when
you run `polyzymd compare run`. For details on interpreting these outputs, see
{doc}`../reference/analysis_comparison_reference`.

## Output structure

After running `polyzymd compare run`, your project directory will contain:

```text
comparison_project/
├── comparison.yaml
├── analysis/
│   └── <condition_label>/
│       └── <plugin_name>/
│           ├── run_<N>/        # Per-replicate results
│           └── aggregated/     # Cross-replicate aggregation
├── comparison/
│   └── <plugin_name>/
│       └── result.json         # Cross-condition comparison
└── figures/
    └── <plugin_name>/
        └── *.png               # Plots
```

The three output directories map directly to the pipeline stages:

- **`analysis/`** holds the compute and aggregate output. Each condition gets
  its own subdirectory, and within that, each plugin gets a directory with
  per-replicate (`run_1/`, `run_2/`, ...) and `aggregated/` results.
- **`comparison/`** holds the compare output. One `result.json` per plugin
  with the statistical tests and rankings.
- **`figures/`** holds the plot output. One subdirectory per plugin with PNG
  files.

Results are cached: if you rerun the pipeline without changing settings, the
compute stage skips replicates that already have results on disk.

## See also

- {doc}`../tutorials/first_analysis` — Hands-on tutorial for running your
  first analysis
- {doc}`../how_to/analysis_compare_conditions` — Practical guide to setting up
  a multi-condition comparison
- {doc}`../reference/comparison_yaml` — Full `comparison.yaml` schema reference
- {doc}`../reference/analysis_comparison_reference` — Plugin listing and
  statistical terms reference
