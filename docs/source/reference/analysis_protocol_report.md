# ProtocolReport schema

`polyzymd.analyses.protocols.ProtocolReport` is the return value of
`protocols.analyze` and the document printed by `polyzymd analyze --format
json`. It is a Pydantic model, so `ProtocolReport.model_validate_json(text)`
reads back exactly what `model_dump_json()` wrote.

## ProtocolReport

| Field | Type | Meaning |
|---|---|---|
| `analysis` | `str` | Canonical analysis name, for example `rg`. |
| `protocol_version` | `str` | Version of the report protocol, `"2"` for every report of the study API. With `analysis` it identifies the code that defined the metric; it changes when the meaning, unit or estimator of a reported metric changes. |
| `metric` | `str` | Name of the reported values, for example `mean_rg`. |
| `unit` | `str \| None` | Unit of `metric`, for example `A` or `%`. `None` marks a dimensionless metric, such as a fraction. |
| `run` | `str \| None` | Selected result when the analysis reports several: for sasa each context's total and per-residue SASA; for distances each pair's mean distance (`<label>`) and fraction below threshold; for hydrogen_bonds each summary's counts, lifetimes, per-residue and per-pair occupancies; for rmsf and rmsd_per_residue the core, region and plain-mean values and the per-residue profiles. `None` when the analysis reports one result. |
| `all_metrics` | `list[str]` | Other metric names of the result, `metric` first. Empty in the reports of the shipped analyses, which name each result in `all_runs` instead. |
| `all_runs` | `list[str]` | Every result the analysis measured, the reported one first. Empty when it measures one. Select another with `--run LABEL`. |
| `equilibration` | `str` | Equilibration window discarded from the start of every replicate, for example `10ns`. Applied uniformly to every replicate of every condition. |
| `stride` | `int` | Every `stride`-th production frame was measured, 1 by default. |
| `frames_per_replicate` | `dict[str, int \| list[int] \| None]` | Production frames each replicate of a condition contributed after the equilibration window and stride, keyed by condition label, one number per replicate in replicate order. |
| `conditions` | `list[ConditionReport]` | One entry per condition, in the order the configs were given; for a labelled result such as a per-residue profile, one entry per condition and label. |
| `pairwise` | `list[PairwiseReport]` | One entry per comparison of the primary metric, or per comparison and label for a labelled result. Empty for a single condition. |
| `warnings` | `list[str]` | Sampling and interval warnings: conditions with one replicate or the same value in every replicate, intervals that cross the bounds of a bounded quantity, untestable comparisons, and pymbar equilibration diagnostics. |
| `provenance` | `ProtocolProvenance` | Versions, config hashes and output paths. |
| `verdict` | `list[str]` | One sentence per pairwise comparison, or one sentence describing the single condition. |

`ProtocolReport.to_agent_text()` renders the report as plain text with no
borders and no blank lines, one line for each condition, comparison, warning
and verdict. No line is dropped, however many conditions the report holds.
For a labelled result the comparisons are summarised instead: one line per
compared condition gives the number of labels, the number tested, the family
size and how many labels are significantly lower and higher, followed by lines
listing each significant label. A final `note:` line says how many per-label
rows the JSON form holds; every one of them is kept there.

## ConditionReport

| Field | Type | Meaning |
|---|---|---|
| `label` | `str` | Condition label, from `--label` or the config's parent directory name. |
| `entry` | `str \| None` | Label of this row in a labelled result, such as a residue ID. `None` for a result with one value per replicate. |
| `n_replicates` | `int` | Number of replicates behind `mean`. This is the sample size for every test. |
| `mean` | `float` | Mean of the primary metric across replicates. |
| `sem` | `float \| None` | Standard error of that mean across replicates, `s / sqrt(n)` with `ddof = 1`. `None` for one replicate, where it does not exist. |
| `ci95` | `tuple[float, float] \| None` | Limits of the 95 percent Student t interval on the mean. `None` for one replicate. |
| `ci_method` | `str \| None` | `student_t` when the interval came from `replicate_values`; `not_estimable` when every replicate has the same value. `None` when no interval exists. |
| `replicate_values` | `list[float]` | The per-replicate values behind the mean, in replicate order. |
| `replicates` | `list[int]` | The replicate number of each entry of `replicate_values`. |
| `statistical_inefficiency` | `list[float]` | For a value reduced from a per-frame series, the pymbar statistical inefficiency of each replicate's series. Empty otherwise. |
| `n_effective` | `list[float]` | The effective sample size of each replicate's series, alongside `statistical_inefficiency`. |
| `eq_detected_frame` | `list[int]` | Start of the equilibrated region that pymbar `detect_equilibration` finds in each checked replicate's production series, as a production frame index from 0. A diagnostic: it changes no value. |
| `eq_detected_ns` | `list[float]` | The same start as simulation time in ns. |

## PairwiseReport

| Field | Type | Meaning |
|---|---|---|
| `a` | `str` | Control condition label. |
| `b` | `str` | Compared condition label. |
| `entry` | `str \| None` | Label compared in a labelled result, such as a residue ID. `None` for a result with one value per replicate. |
| `delta` | `float` | `mean(b) - mean(a)`, in the metric's unit. |
| `delta_ci95` | `tuple[float, float] \| None` | 95 percent Student t interval on `delta`, uncorrected for multiplicity. Pooled variance with `n_a + n_b - 2` degrees of freedom for `student_t`, separate variances with Welch-Satterthwaite degrees of freedom (passed to the quantile unrounded) for `welch_t`. `None` when a condition has fewer than two replicate values, or when both conditions have zero variance. |
| `p` | `float \| None` | Unadjusted p value of the two-sample test. |
| `p_adjusted` | `float \| None` | Benjamini-Hochberg adjusted p value. `None` for a row that was not tested. |
| `test` | `str` | `welch_t` (the default of `polyzymd analyze` and `compare()`) or `student_t`. |
| `correction` | `str` | `BH` (Benjamini-Hochberg). |
| `family_size` | `int \| None` | Number of tests in the Benjamini-Hochberg family this row was corrected in: the conditions compared with the control for this one outcome, or for a labelled result every tested label of every compared condition. `None` for a row that was not tested. Printed as `family <m>` on the comparison line. |
| `cohens_d` | `float \| None` | Standardised mean difference, the difference of the means over the pooled standard deviation, oriented like `delta`: positive means `b` is larger. |
| `hedges_g` | `float \| None` | `cohens_d` with the Hedges small-sample correction, oriented like `cohens_d`. |
| `direction` | `str` | `increased` or `decreased` for a significant row, otherwise `no significant change`. |
| `significant` | `bool` | Whether `p_adjusted` is at most 0.05. Always `False` when `testable` is `False`. |
| `testable` | `bool` | `False` when a condition has fewer than two replicates, or both conditions have the same value in every replicate, which makes the test undefined rather than non-significant. |

## ProtocolProvenance

| Field | Type | Meaning |
|---|---|---|
| `polyzymd_version` | `str` | Version of the package that ran the protocol. |
| `mdanalysis_version` | `str \| None` | Installed MDAnalysis version, `None` when it cannot be determined. |
| `config_hashes` | `dict[str, str]` | `polyzymd.analyses.identity.compute_config_hash` of each simulation config, keyed by condition label: the first 16 hex characters of the SHA-256 of the config fields that locate and describe its trajectories. |
| `settings_fingerprint` | `str \| None` | `None`: no shipped analysis sets it. |
| `settings` | `dict` | Settings a study-API analysis ran with. For rmsf and rmsd_per_residue: every setting, the resolved `reference_mode`, and under `residues` the residue IDs of the core and of each region. Empty for other analyses. |
| `study` | `dict \| None` | Set for a run from a study file: `path` and `sha256` of `study.yaml`, the `run`, and `git`, with the study folder's `commit`, its `uncommitted` files and `inputs_uncommitted`, those outside `results/` and `data.local.yaml`; `git` is `None` outside a repository. |
| `output_paths` | `dict[str, str]` | `results` is the `polyzymd_results/<name>/` folder holding every replicate's stored values and record; `figures` is the directory holding the generated plots, absent with `--no-plots`. |

## Verdict vocabulary

The first words of a verdict sentence come from a fixed set, so a caller can
branch on them without parsing the rest.

| Word | Meaning |
|---|---|
| `larger` | `b` differs from `a` after correction and `delta` is positive |
| `smaller` | `b` differs from `a` after correction and `delta` is negative |
| `no significant difference` | the test ran and `p_adjusted` did not clear alpha |
| `no test recorded` | the row has no multiplicity-corrected p value, so it describes a difference without deciding it |
| `changed` | the difference is significant but the two means are equal at the stored precision |
| `not testable` | a condition has fewer than two replicates, so no test exists |

A single-condition report has no comparison, and its one sentence states the
label, the metric, the mean with its unit, the interval and the replicate
count.

## See also

- {doc}`cli_reference` for `polyzymd analyze` and the `agent` format.
- {doc}`../how_to/analysis_agent_protocol` for the recipe.
- {doc}`../explanation/analysis_entry_points` for when to use the command
  line and when the Python study API.
