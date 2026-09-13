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
| `secondary_structure` | Secondary structure content (helix, strand, coil and unassigned fractions) |
| `sasa` | Solvent-accessible surface area |
| `hydrogen_bonds` | Hydrogen bond counts per partition and residue-pair occupancy |

Each plugin has a `Settings` model with configurable parameters. Most
parameters have sensible defaults, so you often just need `plugin_name: {}` in
your `comparison.yaml` to get started.

For contributors, the plugin boundary is the supported extension point: a
plugin defines its settings, replicate computation, aggregation behavior,
comparison behavior, plotting behavior, and formatting behavior without changing
the core orchestration code. The conceptual boundary is important because
PolyzyMD owns artifact storage and orchestration, while plugins own the
domain-specific measurement and interpretation logic. For a contributor-focused
walkthrough, see {doc}`../contributor_guide/analysis_plugins/index`.

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
## Periodic boundaries, whole molecules, and alignment

Molecular dynamics runs in a periodic box, and the coordinates written to a
trajectory are usually wrapped back into the primary cell. A molecule that
drifts across a face of the box is then stored in pieces, with some atoms at one
edge and the rest at the opposite edge. Nothing in the file says whether that
happened, so an observable computed from raw coordinates can be right on one
frame and badly wrong on the next.

Which observables care depends on what they measure. A radius of gyration is a
spread about a center of mass, so a molecule stored in two pieces reports a
radius of roughly half a box length rather than its real size. Alignment and
RMSD have the same problem, because fitting a structure that is split in two
fits a shape the molecule never adopted. A pair distance is different. The
minimum image convention gives the distance to the nearest periodic copy of the
second point, which is the physically meaningful separation whether or not the
molecules were wrapped, so distances need the box rather than whole molecules.
PolyzyMD therefore makes wholeness an explicit choice at load time through
`pbc_policy`, records the choice in provenance, and refuses to unwrap a topology
that has no bonds, because unwrapping walks the bond graph.

Alignment does not belong anywhere in a distance calculation, for two reasons.
The first is that a distance is invariant under rotation and translation, so
removing rigid-body motion cannot change a correct answer, and paying for a full
in-memory copy of the trajectory to do it buys nothing. The second is that it
was actively harmful here. MDAnalysis applies the minimum image convention by
building box vectors from the stored `(a, b, c, alpha, beta, gamma)` and
assuming those vectors describe the lattice of the coordinates it is handed.
Aligning in memory rotates every coordinate but copies the box unchanged, so the
two no longer agree. For a pair separated by less than half a box length in
every Cartesian component the wrapping is a no-op and the answer survives, which
is why catalytic-triad distances looked fine. For a long pair, such as a domain
center of mass against a polymer center of mass in a 90 Angstrom box, a
component can be folded against the wrong lattice vector and the reported
distance is wrong. Distances and the catalytic triad now measure the coordinates
as the trajectory stores them.

### References

- Michaud-Agrawal, N., Denning, E. J., Woolf, T. B., and Beckstein, O. (2011).
  MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
  *Journal of Computational Chemistry*, 32(10), 2319-2327.
  doi:10.1002/jcc.21787
- MDAnalysis `MDAnalysis.lib.distances` documentation, on the `box` argument and
  the minimum image convention:
  <https://docs.mdanalysis.org/stable/documentation_pages/lib/distances.html>
- MDAnalysis `MDAnalysis.transformations.wrap` documentation, on `unwrap` and
  its bond requirement:
  <https://docs.mdanalysis.org/stable/documentation_pages/transformations/wrap.html>

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

## Why a cached result is checked against its inputs

A replicate that gained three segments since its result was computed will keep
reporting the old window for as long as the cache is reused, so a figure
regenerated a week later still describes last week's data. The failure is quiet
in the same way an unfinished segment is quiet: the cached number is a real
average over the frames that were read, and the provenance describes the run
rather than the subset, so nothing in the output says the two have drifted
apart.

Reuse is therefore conditional rather than automatic. A cached result is reused
only when every input file it names still has the recorded size and
modification time, when the set of trajectory files has not changed, and when
the settings and equilibration window are the ones it was computed under. Size
and modification time are weaker than a content hash, but they are cheap on a
multi-gigabyte trajectory and they catch the case that actually happens, which
is a file that grew. A result that cannot prove any of this, including one
written before the framework recorded a cache key, is recomputed rather than
trusted; where the command that found it cannot recompute, it stops and says
which file changed.

## Why a solvent-accessible surface area depends on how it was batched

Some measurements are not a pure function of a frame. `mdtraj.shrake_rupley`,
which the `sasa` plugin calls, gives the same coordinates slightly different
areas depending on how many frames are in the array it is handed, so the
plugin's `chunk_size` shifts every total by about 0.1 percent. That is far below
the differences the analysis is used to detect, but it is an offset between two
runs rather than noise that averages away, so it can only be ignored when both
sides of a comparison were computed the same way. PolyzyMD does not hide this by
forcing one chunk size, because the setting exists to keep a large system inside
memory; it records the value in each observable's metadata instead, so a
comparison assembled from mismatched runs can be recognized rather than
believed.

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

## Reading two numbers that measure the same thing differently

A plugin often reports several scalars about one phenomenon, and they answer
different questions. Contacts is the clearest case. A higher `contact_count`
than the control means more protein-polymer residue pairs are in contact at any
moment; a higher `coverage_per_frame` means those contacts are spread over more
of the protein rather than concentrated. A formulation can raise the count and
leave coverage flat by binding one patch harder, so the pair has to be read
together. The `contact_fraction` profile then names the patch, and
`residence_time_distribution` says whether the contacts are many brief touches
or a few long ones, which no scalar can distinguish.

This is why the framework tests every observable a plugin declares rather than
one headline metric per plugin, and why a plugin marks an observable that is a
function of the others `tested=False`: reporting it is useful, testing it would
enlarge the correction family without adding information.

## See also

- {doc}`../tutorials/first_analysis` — Hands-on tutorial for running your
  first analysis
- {doc}`../how_to/analysis_compare_conditions` — Practical guide to setting up
  a multi-condition comparison
- {doc}`../reference/comparison_yaml` — Full `comparison.yaml` schema reference
- {doc}`../reference/analysis_comparison_reference` — Plugin listing and
  statistical terms reference
