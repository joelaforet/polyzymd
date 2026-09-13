# Comparison and Plotting Reference

Use this page when you need quick lookup information for `polyzymd compare`,
`comparison.yaml`, output paths, or plotting behavior.

## Comparison Project Layout

`polyzymd compare init -n my_study` creates a workspace like this:

```text
my_study/
├── comparison.yaml
├── comparison/
├── figures/
└── structures/
```

## Core `comparison.yaml` Fields

```yaml
name: "polymer_stability_study"
description: "Optional human-readable summary"
control: "No Polymer"  # optional

conditions:
  - label: "No Polymer"
    config: "no_polymer/config.yaml"
    replicates: [1, 2, 3]

defaults:
  equilibration_time: "10ns"
  ttest_method: "student"    # or "welch"
  posthoc_method: "ttest_bh" # or "tukey_hsd"
  fdr_alpha: 0.05

plugins:
  rmsf:
    selection: "protein and name CA"
```

## Hypothesis testing across plugins

`ttest_method`, `posthoc_method` and `fdr_alpha` from the `defaults:` block
reach every plugin, because every plugin is compared by the same code:
`polyzymd.analyses.contract.compare_observables`. There is no second comparison
path and no per-plugin comparison code.

One observable is tested at a time, each condition against the control, on the
replicate values. The rules are:

- **One run, one family.** Every pairwise test the run produced, across every
  observable and every condition pair, is corrected together with the
  Benjamini-Hochberg step-up procedure. Each comparison entry carries both
  `p_value` and `p_adjusted`, and `significant` is read from the adjusted
  value. With `posthoc_method: "tukey_hsd"` the family-wise Tukey adjustment
  replaces Benjamini-Hochberg and `correction` says which one ran.
- **No ANOVA.** Nothing reports an omnibus F test any more. It gated nothing,
  and an uncorrected omnibus p-value next to corrected pairwise ones invited
  the reader to treat it as a gate.
- **An untested observable stays out of the family.** A plugin declares
  `tested=False` for a quantity that is a function of others it already
  reports, such as the last of a set of fractions that sums to one. It is still
  aggregated and reported with its uncertainty, and it does not enter the
  pairwise tests or inflate the adjusted p-values of the rest.
- **A profile is not tested pairwise.** It has no single value. A plugin gives
  it a comparable scalar with `reduce`, reported as `<name>_mean` or
  `<name>_total`, and that scalar is tested.
- **Effect sizes.** Each comparison carries `cohens_d`. The agent report adds
  `hedges_g`, Cohen's d multiplied by `J = 1 - 3 / (4 * (n1 + n2) - 9)`.
- **Direction labels require significance.** A direction word is only assigned
  when the corrected test is significant. Otherwise the field reads
  `"no significant change"`.

Full field tables are in {doc}`posthoc_testing`.

## Stable plugin keys

- `rmsd`
- `rg`
- `rmsf`
- `contacts`
- `distances`
- `catalytic_triad`
- `secondary_structure`
- `sasa`
- `hydrogen_bonds`

## What each plugin reports

Every plugin reports observables, and the `kind` of each one decides how it is
reduced per replicate, whether it is tested, and which figure is drawn. The
names below are the keys that appear in the comparison artifact.

| Plugin | Observables | Kinds | Unit |
|--------|-------------|-------|------|
| `rmsd` | `rmsd_<run slug>_ref_<reference mode>`, one per configured run | `mean_of_timeseries` | `A` |
| `rg` | `rg_<slug>`; in fragment mode also `rg_<slug>_fragments` and `rg_<slug>_distribution` | `mean_of_timeseries`, `profile`, `profile` | `A`, `A`, `1/A` |
| `rmsf` | one per-residue profile per configured selection, reduced to `<name>_mean` | `profile` with `reduce="mean_over_index"` | `A` |
| `contacts` | `contact_count`, `coverage_per_frame`, `coverage_any_frame`, `contact_fraction`, `residence_time_distribution` | `mean_of_timeseries`, `fraction`, `fraction`, `profile`, `profile` | `count`, `fraction`, `fraction`, `fraction`, `ns` |
| `distances` | `<pair label>` and `<pair label> <state>` per configured pair | `mean_of_timeseries`, `fraction` | `A`, `fraction` |
| `catalytic_triad` | `<pair label>`, `<pair label> within <cutoff> A`, and `simultaneous_contact_fraction` | `mean_of_timeseries`, `fraction`, `fraction` | `A`, `fraction`, `fraction` |
| `secondary_structure` | `ss_helix`, `ss_strand`, `ss_coil`, `ss_unassigned`, plus a per-residue `<label>_occupancy` profile | `fraction`, `profile` | `fraction` |
| `sasa` | `sasa_<label>` and `relative_sasa_<label>` per measured context | `mean_of_timeseries`, `profile` | `A^2`, `fraction` |
| `hydrogen_bonds` | `hbonds_<summary>` per configured summary, plus `pair_occupancy_<summary>` | `mean_of_timeseries`, `profile` | `count`, `fraction` |

The per-frame series behind every observable is written to
`run_<replicate>/observables.npz`, keyed by observable name, so a value can be
inspected frame by frame without rerunning the analysis.

## Path Rules

- relative paths in `config:` are resolved relative to `comparison.yaml`
- absolute paths are used as-is
- `replicates` must be an explicit list such as `[1, 2, 3]`

## Replicate Counts

All stable shipped analyses support `replicates: [1]` for smoke tests and
protocol validation. One-replicate runs compute aggregates and plots, but a
standard error, a confidence interval and a hypothesis test all need at least
two independent replicates per condition. A singleton pairwise comparison is
reported with `testable: false` and a `note` saying why, rather than as not
significant.

## Commands

| Command | Purpose |
|---------|---------|
| `polyzymd compare init -n NAME` | Create a comparison workspace |
| `polyzymd compare validate` | Check `comparison.yaml` before running |
| `polyzymd compare run TYPE` | Run one analysis plugin |
| `polyzymd compare run --list` | List available comparison types |
| `polyzymd compare run-all` | Run every enabled plugin in one pass |
| `polyzymd compare plot-all` | Generate configured figures |
| `polyzymd compare plot-all --list-available` | List available plots and experimental labels |
| `polyzymd compare submit ANALYSIS` | Submit a SLURM DAG for one analysis plugin |
| `polyzymd compare status ANALYSIS` | Show status of a submitted SLURM DAG |
| `polyzymd compare finalize ANALYSIS` | Run comparison + plotting from on-disk aggregated results |

## Common Stable Commands

All commands below assume you are inside the pixi environment
(`pixi shell -e analysis`) or are prefixed with `pixi run -e analysis`.

```bash
polyzymd compare run rmsd
polyzymd compare run rg
polyzymd compare run rmsf
polyzymd compare run contacts
polyzymd compare run distances
polyzymd compare run catalytic_triad
polyzymd compare run sasa
polyzymd compare run hydrogen_bonds
polyzymd compare run-all
polyzymd compare plot-all
```

## Output Locations

- per-replicate cache files are written under
  `analysis/<condition>/<analysis>/run_<replicate>/`
- per-condition aggregate files are written under
  `analysis/<condition>/<analysis>/aggregated/`
- cross-condition comparison JSON files are written to
  `comparison/<analysis>/result.json`
- figures are written under the configured `plot_settings.output_dir`, usually
  `figures/<analysis>/`
- `polyzymd compare init` scaffolds `comparison/`, `figures/`, and
  `structures/` next to `comparison.yaml`; `analysis/` is created and
  populated during analysis runs

Typical comparison cache paths:

```text
comparison/rmsd/result.json
comparison/rg/result.json
comparison/rmsf/result.json
comparison/contacts/result.json
comparison/distances/result.json
comparison/catalytic_triad/result.json
comparison/sasa/result.json
comparison/hydrogen_bonds/result.json
```

## Cache reuse and input freshness

Every per-replicate, per-condition, and comparison artifact records the files
it was computed from. The record lives at
`provenance.universe_policy.provenance` and holds the path, size in bytes, and
modification time in nanoseconds of the topology and of each trajectory
segment. The framework also stamps a cache key on every replicate artifact it
writes, under `metadata`:

| `metadata` key | Meaning |
|---|---|
| `settings_fingerprint` | 8 hex characters of the canonical settings JSON |
| `equilibration` | the equilibration window the frame selection used |

The software versions are top-level fields on the artifact envelope, not
metadata keys, and the framework records them on every artifact it writes,
including condition and comparison artifacts:

| Envelope field | Meaning |
|---|---|
| `polyzymd_version` | the PolyzyMD version that wrote the artifact |
| `mdanalysis_version` | the MDAnalysis version that wrote the artifact |

`polyzymd compare run` reuses a cached `run_<replicate>/result.json` only when
all of the following hold. Anything else recomputes that replicate.

- Every recorded input file still has the recorded size and modification time.
- The set of trajectory files the engine resolves now is the same set the
  result recorded, so a segment that has appeared or completed since counts as
  a change even though no recorded file was touched.
- The stored `settings_fingerprint` and `equilibration` match the running
  command. A result that records neither, including any artifact written before
  this key existed, is treated as stale rather than reused. The first run after
  upgrading therefore recomputes replicates whose artifacts predate the key.

The other commands behave as follows:

- `--recompute` skips every check and recomputes.
- `polyzymd compare finalize`, `plot-all`, and the worker commands have no
  compute stage. When a replicate result they load was computed from files that
  have since changed, or records a different cache key, they raise
  `StaleCacheError`, name the file, and ask for `--recompute`. These commands
  read artifacts only and never open a trajectory, so they check the identity
  the artifact recorded rather than re-resolving the engine layout; the check
  for segments that appeared since a result was written runs in `compare run`,
  which resolves the layout anyway.
- An aggregate read from `aggregated/result.json` is refused when any
  `run_<replicate>/result.json` beside it was written later. The check applies
  only to an aggregate read from disk, not to one being written.
- A cache written by a different PolyzyMD version is used with a warning, not
  refused, since a version difference alone does not invalidate a number.

## Incomplete production segments

The OpenMM engine reads the status recorded for each production segment in
`progress.json` and leaves out any segment marked `running` or `failed`, since
its DCD is still being appended to or stopped at an arbitrary step. The
excluded indices appear in `excluded_segments` on the trajectory layout, in
`TrajectoryInfo`, and in the `provenance.universe_policy.provenance` block of
every artifact, next to `segment_status`. A warning names them.

`TrajectoryLoader.load_universe`, `TrajectoryLoader.get_trajectory_info`, and
`UniverseProvider` take `require_complete`, which defaults to `True`. Set it to
`False` to read a campaign that is still running; the segments are then
included and listed in `incomplete_segments` on the layout instead.

What this does and does not guarantee:

- Segments marked `interrupted` are kept, because the continuation chain
  resumes from an interrupted segment's saved state and its frames belong to
  the time line. An interrupted segment is written with `samples_written = 0`
  and the comment that samples may be partial
  (`simulation/continuation.py:828`), so its frame count is not recorded
  anywhere. Treat the last interrupted segment of a chain as partial.
- Excluding a segment that other segments continue from leaves a hole in the
  concatenated time line. The lineage check reports it and names the excluded
  segments as a possible cause.
- The contiguity check reads only the first and last time of each segment plus
  its frame interval. It proves that segment boundaries line up. It does not
  detect a dropped or duplicated frame inside a segment.
- Recorded input paths are absolute. A run directory copied or moved to another
  path does not match its recorded identity, so analyses recompute rather than
  reuse, and commands that only read caches raise instead.

The GROMACS engine records no per-file status, because its layout is a single
production XTC rather than a chain of segments. It accepts `require_complete`
for interface parity and ignores it, and its `segment_status` and
`excluded_segments` are always empty.

For what to do when a run is still in flight, see
{doc}`../how_to/analysis_compare_conditions`.

## Plotting Smoke Test

For a final smoke test after comparisons finish:

```bash
polyzymd compare plot-all --list-available
polyzymd compare plot-all
```

## Uncertainty fields

Every observable aggregate in a condition or comparison artifact carries these
fields, under `payload["observables"]`. The replicate is the sampling unit for
all of them.

| Field | Type | Description |
|-------|------|-------------|
| `mean` | `float` | Mean of the replicate values |
| `sem` | `float \| null` | Standard error of the mean across replicates, `s / sqrt(n)` with `ddof = 1`. `null` when a single replicate makes it inestimable |
| `std` | `float \| null` | Sample standard deviation across replicates. `null` for a single replicate |
| `n_replicates` | `int` | Number of replicates |
| `replicate_values` | `list[float]` | The per-replicate values the mean was taken over |
| `unit` | `str \| null` | Physical unit of `mean`, for example `"A"`, `"A^2"`, `"fraction"`, `"%"`, `"ns"`. `null` marks a dimensionless metric |
| `ci95_low` | `float \| null` | Lower limit of the 95 percent confidence interval. `null` for a single replicate |
| `ci95_high` | `float \| null` | Upper limit of the 95 percent confidence interval. `null` for a single replicate |
| `ci_method` | `str \| null` | Interval method, `"student_t"` when an interval exists |

The interval is `mean +/- t(0.975, n - 1) * sem`. The coverage factor is 4.303
at `n = 3` and 2.776 at `n = 5`. The interval is symmetric and is not clipped
to a physical range, so for a fraction, an occupancy or a coverage near 0 or 1
a limit can fall outside `[0, 1]`. Read such a limit as an indication of spread
rather than a bound.

Each condition artifact payload and each comparison artifact payload also
carries an `uncertainty` block:

```json
{
  "uncertainty": {
    "kind": "sem_across_replicates",
    "n": 3,
    "coverage": 0.95,
    "method": "student_t"
  }
}
```

`kind` names the standard uncertainty, `n` the replicate count behind the
narrowest interval in the artifact, `coverage` the probability the interval
covers, and `method` how the coverage factor was obtained.

## Statistical Terms

- `p-value`: significance of the observed difference under the null hypothesis
- `Cohen's d`: effect size magnitude, the mean difference over the pooled
  standard deviation
- `Hedges' g`: Cohen's d after the small-sample bias correction
  `J = 1 - 3 / (4 * (n1 + n2) - 9)`
- `SEM`: standard error of the mean across replicates
- `95% CI`: two-sided Student t confidence interval across replicates,
  `mean +/- t(0.975, n - 1) * SEM`
- `Benjamini-Hochberg (BH)`: step-up procedure for controlling the false
  discovery rate across multiple hypothesis tests
- `Adjusted p-value (p_adj)`: p-value corrected for multiple comparisons via
  the BH procedure (for `ttest_bh`) or family-wise Tukey adjustment (for
  `tukey_hsd`)
- `False Discovery Rate (FDR)`: expected proportion of false positives among
  rejected hypotheses
- `Effect size threshold`: minimum Cohen's d required for a pairwise difference
  to be considered practically significant

For interpretation guidance rather than lookup, see:

- [Statistical Best Practices for Analysis](../explanation/analysis_statistics_best_practices.md)
- [How to Compare Simulation Conditions](../how_to/analysis_compare_conditions.md)
- [Post-Hoc Testing Reference](posthoc_testing.md) — full post-hoc method details, output fields, and edge cases
