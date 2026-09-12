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

## Hypothesis Testing Across Plugins

`ttest_method`, `posthoc_method` and `fdr_alpha` from the `defaults:` block
reach every plugin's comparison step, including the plugins listed below as
"custom" in the plugin summary table. Asking for `ttest_method: "welch"` runs
Welch's unequal-variance t-test in `rmsd`, `rg`, `sasa`, `contacts`,
`distances` and the default scalar pipeline alike.

The multiple-comparison family is defined once for the whole package:

- **One run, one family.** Every pairwise test the run produced, across all of
  its metrics and all of its condition pairs, is corrected together with the
  Benjamini-Hochberg step-up procedure. Each pairwise result carries both
  `p_value` and `p_value_adjusted`, and `significant` is read from the
  adjusted value.
- **ANOVA is omnibus and uncorrected.** Its `p_value` is reported raw, its
  `p_value_adjusted` is always `null`, and no pairwise test is gated on it.
  Its `significant` flag compares the raw p-value with `fdr_alpha`.
- **Effect sizes carry `hedges_g` next to `cohens_d`.** Hedges' g is Cohen's d
  multiplied by `J = 1 - 3 / (4 * (n1 + n2) - 9)`. The Cohen adjective in
  `effect_size_interpretation` is `null` when `n1 + n2 < 10`.
- **Direction labels require significance.** A label such as `"stabilizing"`,
  `"increased"`, `"closer"` or `"exposure"` is only assigned when the
  corrected test is significant. Otherwise the field reads
  `"no significant change"`.

Full field tables are in {doc}`posthoc_testing`.

## Per-Plugin Statistical Settings

Some plugins support per-plugin statistical settings configured under the
`plugins:` block in `comparison.yaml`. These control false discovery rate
correction, effect-size filtering, and output truncation for cross-condition
comparisons.

### Canonical YAML Example

```yaml
plugins:
  contacts:
    cutoff: 4.5
    fdr_alpha: 0.05
    min_effect_size: 0.5
    top_residues: 10
```

### Settings Support Matrix

| Setting | contacts | Default |
|---------|----------|---------|
| `fdr_alpha` | ✓ | 0.05 |
| `min_effect_size` | ✓ | 0.5 |
| `top_residues` | ✓ | 10 |

### Setting Descriptions

- **`fdr_alpha`** — Significance threshold for pairwise comparisons. When
  `posthoc_method` is `"ttest_bh"`, this controls the Benjamini-Hochberg false
  discovery rate. When `posthoc_method` is `"tukey_hsd"`, this is the
  family-wise alpha threshold. Also used as the ANOVA significance threshold.
  Lower values are more conservative.
- **`min_effect_size`** — Minimum Cohen's d required for practical
  significance. Pairs that meet or exceed this threshold are highlighted with
  "†" in formatted output; all pairs are shown regardless.
- **`top_residues`** — Maximum number of contacted residues shown per
  condition, ranked by aggregated `contact_fraction_mean`. Affects both saved
  JSON and CLI output.

## Stable Plugin Keys

Stable analysis plugins:

- `rmsd`
- `rg`
- `rmsf`
- `contacts`
- `distances`
- `catalytic_triad`
- `secondary_structure`
- `sasa`
- `hydrogen_bonds`

## Plugin Summary Table

| Plugin | Default compare? | Primary metric | Key feature | Statistical method |
|--------|-----------------|----------------|-------------|-------------------|
| `rmsd` | No (custom) | `mean_rmsd` | Backbone stability over time | FDR-corrected per-run pairwise t-tests + omnibus ANOVA |
| `rg` | No (custom) | `mean_rg` | Protein compactness | FDR-corrected per-run pairwise t-tests + omnibus ANOVA |
| `rmsf` | Yes | `mean_rmsf` | Per-residue flexibility | FDR-corrected pairwise t-tests + omnibus ANOVA |
| `contacts` | No (custom) | Coverage + contact fraction | Per-residue contact mapping | FDR-corrected pairwise t-tests + omnibus ANOVA |
| `distances` | No (custom) | Multiple distance metrics | Named distance pairs | FDR-corrected per-pair t-tests + omnibus ANOVA |
| `catalytic_triad` | Yes | `simultaneous_contact_fraction` | Active-site geometry | FDR-corrected pairwise t-tests + omnibus ANOVA |
| `secondary_structure` | Yes | `ss_helix`, `ss_strand`, `ss_coil`, `ss_unassigned` | Secondary structure content | FDR-corrected pairwise t-tests over replicates |
| `sasa` | No (custom) | Per-run mean SASA | Multi-run target/context model | FDR-corrected per-run pairwise t-tests + omnibus ANOVA |
| `hydrogen_bonds` | Custom loader with default-style scalar statistics | `mean_hbonds_per_frame` per summary | Flexible named groups + summaries + composition analysis | FDR-corrected pairwise t-tests + ANOVA per configured summary |

## Path Rules

- relative paths in `config:` are resolved relative to `comparison.yaml`
- absolute paths are used as-is
- `replicates` must be an explicit list such as `[1, 2, 3]`

## Replicate Counts

All stable shipped analyses support `replicates: [1]` for smoke tests and
protocol validation. One-replicate runs compute aggregate metrics and plots, but
inferential statistics, FDR correction, and uncertainty bands require at least
two independent replicates per condition. Singleton pairwise tests and ANOVA are
reported as not testable rather than significant.

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

## Plugin-Specific Metadata Fields

Some plugins include additional metadata in their comparison output beyond the
standard ranking and statistical fields. These fields are **additive
diagnostics** — they do not affect rankings, p-values, or effect sizes.

### RMSD Convergence Output

The RMSD plugin includes per-run convergence diagnostics generated by the
sliding-window convergence heuristic in `analyses/shared/convergence.py`.
These fields appear in the per-condition summaries within
`comparison/rmsd/result.json`:

| Field | Type | Description |
|-------|------|-------------|
| `convergence_fraction` | `float` | Fraction of replicates that converged (0.0–1.0) |
| `n_converged_replicates` | `int` | Count of replicates where sustained convergence was detected |
| `mean_convergence_time_ns` | `float \| null` | Mean convergence time across converged replicates (ns) |
| `median_convergence_time_ns` | `float \| null` | Median convergence time across converged replicates (ns) |

:::{note}
Convergence metadata is purely informational. It does not influence the RMSD
ranking, pairwise t-tests, ANOVA, or effect-size calculations. Use it to
identify conditions where one or more replicates failed to reach a stable
plateau, which may warrant longer production runs or additional replicates.
:::

## Uncertainty fields

Every metric summary in a condition or comparison artifact carries these
fields. The replicate is the sampling unit for all of them.

| Field | Type | Description |
|-------|------|-------------|
| `mean` | `float` | Mean of the replicate values |
| `sem` | `float \| null` | Standard error of the mean across replicates, `s / sqrt(n)` with `ddof = 1`. `null` when a single replicate makes it inestimable |
| `std` | `float \| null` | Sample standard deviation across replicates. `null` for a single replicate |
| `n` | `int` | Number of replicates |
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

The default scalar comparison also writes `<metric>_unit`, `<metric>_ci95_low`,
`<metric>_ci95_high` and `<metric>_ci_method` into each entry of
`condition_summaries`, alongside the existing `<metric>_mean`, `<metric>_sem`
and `<metric>_replicate_values`.

## Statistical Terms

- `p-value`: significance of the observed difference under the null hypothesis
- `Cohen's d`: effect size magnitude, the mean difference over the pooled
  standard deviation
- `Hedges' g`: Cohen's d after the small-sample bias correction
  `J = 1 - 3 / (4 * (n1 + n2) - 9)`
- `ANOVA`: omnibus test across multiple conditions, reported uncorrected and
  gating nothing
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
