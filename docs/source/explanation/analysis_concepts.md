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
| **replicate stage** | One replicate of one condition | `ReplicateArtifact` at `analysis/<condition_label>/<plugin_name>/run_<N>/result.json` |
| **aggregate** | All replicates of one condition | `ConditionArtifact` at `analysis/<condition_label>/<plugin_name>/aggregated/result.json` |
| **compare** | All conditions together | `ComparisonArtifact` or active custom comparison result at `comparison/<plugin_name>/result.json` |
| **plot** | All conditions together | Figures saved in the configured format, with `png` as the default and `pdf` or `svg` also supported |

Each artifact stores a validated payload plus metadata, provenance, warnings,
and references to sidecar files when an analysis needs large tables or arrays
outside the main JSON document.

Trajectory-native plugins generally create `MDAAnalysisJob` objects for their
per-replicate computation. The corresponding collectors translate completed
jobs into `ReplicateArtifact` objects. PolyzyMD then owns the surrounding
workflow: condition aggregation, cross-condition comparison, artifact storage,
and plot orchestration.

Plots are deliberately downstream of this artifact layer. They read cached
artifacts and sidecars only; they do not reload trajectories or rerun the
analysis calculation.

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
and result files. When labels appear in directory names, PolyzyMD sanitizes them
for the filesystem; for example, `100% SBMA` may be written as `100_SBMA` in
paths while remaining `100% SBMA` in summaries and plots.

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

A **replicate** is a separate run of the same condition, intended to sample the
same setup independently. Replicates are identified by number — 1, 2, 3, and so
on. Each replicate uses the same
`config.yaml` but usually starts from a different random seed, initial velocity
assignment, or starting configuration. These choices help separate trajectories,
but they do not guarantee statistical independence by themselves. Interpretation
also depends on equilibration, stationarity, decorrelation, and whether the
simulated timescales are long enough for the process being measured.

The pipeline processes data in this order:

1. **Per-replicate**: the compute stage runs once for each replicate of each
   condition and writes a `ReplicateArtifact`. If you have 2 conditions with 3
   replicates each, that is 6 replicate artifacts.
2. **Per-condition**: `aggregate` runs once per condition, combining replicate
   artifacts into a `ConditionArtifact`. That's 2 aggregate calls.
3. **Cross-condition**: `compare` runs once, looking at all conditions together
   and writing a `ComparisonArtifact` or an active custom comparison result.
   `plot` then reads those cached outputs and any referenced sidecars.

## Plugins — the analysis modules

PolyzyMD ships with 9 analysis plugins. Each plugin is a self-contained
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

Each plugin has a `Settings` model with configurable parameters. Most
parameters have sensible defaults, so you often just need `plugin_name: {}` in
your `comparison.yaml` to get started.

For contributors, the plugin boundary is the supported extension point: a
plugin defines its settings, replicate computation, aggregation behavior,
comparison behavior, plotting behavior, and formatting behavior without changing
the core orchestration code. The conceptual boundary is important because
PolyzyMD owns artifact storage and orchestration, while plugins own the
domain-specific measurement and interpretation logic. For a contributor-focused
walkthrough, see {doc}`../contributor_guide/extending_analyses`.

You configure plugins in the `plugins:` block. For example, to run RMSF with a
custom selection and contacts with defaults:

```yaml
plugins:
  rmsf:
    selection: "protein and name CA"
  contacts: {}
```

### Why hydrogen bonds exclude carbon

MDAnalysis finds hydrogen bonds from geometry alone: it pairs each hydrogen with
a nearby heavy atom, then keeps the triplets whose donor-acceptor distance and
D-H...A angle pass the cutoffs. Which atoms are allowed to be donors and
acceptors is therefore a scientific choice, not a detail. If that choice is the
whole selection, every aliphatic and aromatic carbon that carries a hydrogen
becomes a donor and every atom becomes an acceptor, so C-H...O and N-H...C
contacts are counted alongside real hydrogen bonds.

The IUPAC definition requires the donor to be more electronegative than
hydrogen and the acceptor to carry a lone pair or a pi cloud, which carbon
generally does not (Arunan et al. 2011). PolyzyMD therefore restricts donors and
acceptors to nitrogen and oxygen by default. The practical reason matters as
much as the formal one: the share of short C-H...O geometries depends on polymer
chemistry, so counting them biases one condition relative to another instead of
shifting every condition by the same amount. Sulfur is a genuine but weaker
donor and acceptor, so it is available through `donor_acceptor_elements` rather
than on by default.

Arunan, E., et al. (2011). Definition of the hydrogen bond (IUPAC
Recommendations 2011). Pure and Applied Chemistry, 83(8), 1637-1641.
doi:10.1351/PAC-REC-10-01-02

## Statistical comparison

When you have two or more conditions, the compare stage produces statistical
output so you can assess whether differences are meaningful. There are two
comparison paths:

- **Default scalar/artifact comparison**: plugins that expose scalar metrics can
  use the framework's default comparison behavior. In that path, PolyzyMD can
  compute pairwise tests, effect sizes, optional omnibus statistics, and metric
  rankings from the condition artifacts.
- **Custom comparison**: plugins with richer result structures can implement
  their own comparison behavior. These plugins still write comparison output,
  but they may not produce the same tests, tables, or rankings as the default
  scalar path.

Every comparison plugin computes:

- **Pairwise t-tests** between each pair of conditions, with
  Benjamini–Hochberg FDR correction over one family per analysis run. The
  family holds every pairwise test that run produced, across all of its
  metrics and all of its condition pairs.
- **Effect sizes** (Cohen's d and Hedges' g) for each pair, so you can see not
  just whether a difference is significant but how large it is.
- **ANOVA** when there are three or more conditions. It is reported as an
  omnibus statement about whether any condition differs at all. It is not
  adjusted, and the pairwise tests run whether or not it reaches significance,
  so it gates nothing.
- **Rankings** of conditions according to each metric's directionality. These
  rankings are screening aids for follow-up interpretation, not biological truth
  by themselves.

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
│           ├── run_<N>/
│           │   └── result.json # ReplicateArtifact
│           └── aggregated/
│               └── result.json # ConditionArtifact
├── comparison/
│   └── <plugin_name>/
│       └── result.json         # ComparisonArtifact or active custom result
└── figures/
    └── <plugin_name>/
        └── *.<format>          # Plots; png by default, pdf/svg supported
```

The three output directories map directly to the pipeline stages:

- **`analysis/`** holds the compute and aggregate output. Each condition gets
  its own filesystem-sanitized subdirectory, and within that, each plugin gets
  a directory with `ReplicateArtifact` files in `run_1/`, `run_2/`, ... and a
  `ConditionArtifact` in `aggregated/result.json`.
- **`comparison/`** holds the compare output. One `result.json` per plugin
  stores a `ComparisonArtifact` or an active custom comparison result. Default
  scalar comparisons include framework-generated tests and rankings; custom
  comparison outputs may use plugin-specific summaries.
- **`figures/`** holds the plot output. One subdirectory per plugin with PNG
  files by default, or another configured format such as PDF or SVG. Plots are
  generated from cached artifacts and sidecars only.

Results are cached: if you rerun the pipeline without changing settings, the
compute stage skips replicates that already have canonical artifacts on disk.

## Why an unfinished segment is not read

A production segment that is still being written has a trajectory file that
ends wherever the last flush landed. Reading it gives a window that is short
for a reason nothing in the result records: the number that comes out is a real
average over fewer nanoseconds than the run directory appears to contain, and
nobody can tell afterwards which it was. The failure does not announce itself,
because a short window still produces a plausible number and the provenance
describes the run rather than the subset that was actually read.

The engine therefore consults the status it already writes for each segment and
refuses to read one marked as running or failed, recording which segments it
left out. Segments marked interrupted are a different case and are kept, since
the continuation chain resumes from an interrupted segment's saved state, so
its frames are part of the same time line. That choice has a cost: an
interrupted segment's frame count is not recorded, so the last one in a chain
may be partial. The alternative, dropping it, would put a hole in the middle of
the time line, which is worse.

## See also

- {doc}`../tutorials/first_analysis` — Hands-on tutorial for running your
  first analysis
- {doc}`../how_to/analysis_compare_conditions` — Practical guide to setting up
  a multi-condition comparison
- {doc}`../reference/comparison_yaml` — Full `comparison.yaml` schema reference
- {doc}`../reference/analysis_comparison_reference` — Plugin listing and
  statistical terms reference
