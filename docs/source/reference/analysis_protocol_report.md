# ProtocolReport schema

`polyzymd.analyses.protocols.ProtocolReport` is the return value of
`protocols.analyze` and the document printed by `polyzymd analyze --format
json`. It is a Pydantic model, so `ProtocolReport.model_validate_json(text)`
reads back exactly what `model_dump_json()` wrote.

## ProtocolReport

| Field | Type | Meaning |
|---|---|---|
| `analysis` | `str` | Canonical analysis name, for example `rg`. |
| `protocol_version` | `str` | The plugin's `Analysis.protocol_version`. With `analysis` it identifies the code that defined the metric. Every plugin starts at `"1"` and bumps it when the meaning, unit or estimator of a reported metric changes. |
| `metric` | `str` | Primary metric key: the first key the plugin's `extract_metrics()` returns. |
| `unit` | `str \| None` | Unit of `metric`, for example `A` or `%`. `None` marks a dimensionless metric, and also a plugin that declares no unit. |
| `run` | `str \| None` | Selected run or pair label, for a plugin that measures the same metric on several selections (rg on `Protein` and `Polymer Oligomers`, sasa on four contexts, distances on each atom pair). `None` when the plugin reports one run. |
| `all_metrics` | `list[str]` | Every metric key the plugin reported, `metric` first. Only `metric` is summarised in `conditions` and `pairwise`. |
| `all_runs` | `list[str]` | Every run or pair label the plugin reported, `run` first. Empty when the plugin reports one run. Select another with `--run LABEL`. |
| `equilibration` | `str` | Equilibration window discarded from the start of every replicate, for example `10ns`. Applied uniformly to every replicate of every condition. |
| `frames_per_replicate` | `dict[str, int \| None]` | Frames each replicate of a condition contributed, keyed by condition label, from the condition artifact's frame-selection provenance. `None` for a plugin that records no frame selection. |
| `conditions` | `list[ConditionReport]` | One entry per condition, in the order the configs were given. |
| `pairwise` | `list[PairwiseReport]` | One entry per comparison of the primary metric. Empty for a single condition. |
| `warnings` | `list[str]` | Sampling warnings first, then warnings carried by the comparison and condition artifacts. Deduplicated. |
| `provenance` | `ProtocolProvenance` | Versions, config hashes and output paths. |
| `verdict` | `list[str]` | One sentence per pairwise comparison, or one sentence describing the single condition. |

`ProtocolReport.to_agent_text(max_lines=25)` renders the report as at most 25
lines of plain text with no borders and no blank lines. Condition and
comparison lines are dropped first when a report does not fit, and the dropped
count is stated on the last line.

## ConditionReport

| Field | Type | Meaning |
|---|---|---|
| `label` | `str` | Condition label, from `--label` or the config's parent directory name. |
| `n_replicates` | `int` | Number of replicates behind `mean`. This is the sample size for every test. |
| `mean` | `float` | Mean of the primary metric across replicates. |
| `sem` | `float \| None` | Standard error of that mean across replicates, `s / sqrt(n)` with `ddof = 1`. `None` for one replicate, where it does not exist. |
| `ci95` | `tuple[float, float] \| None` | Limits of the 95 percent Student t interval on the mean. `None` for one replicate. |
| `ci_method` | `str \| None` | `student_t` when the interval came from `replicate_values`; `student_t_from_sem` when the plugin stored only a mean and a standard error and the interval was rebuilt as `mean` plus or minus `t(0.975, n - 1)` times `sem`. `None` when no interval exists. |
| `replicate_values` | `list[float]` | The per-replicate values behind the mean, in replicate order. Empty for a plugin that stores only summary statistics, such as contacts. |

## PairwiseReport

| Field | Type | Meaning |
|---|---|---|
| `a` | `str` | Control condition label. |
| `b` | `str` | Compared condition label. |
| `delta` | `float` | `mean(b) - mean(a)`, in the metric's unit. |
| `delta_ci95` | `tuple[float, float] \| None` | 95 percent Student t interval on `delta`, uncorrected for multiplicity. Pooled variance with `n_a + n_b - 2` degrees of freedom for `student_t`, separate variances with Welch-Satterthwaite degrees of freedom (passed to the quantile unrounded) for `welch_t`. `None` for `tukey_hsd`, whose interval is a studentised-range interval rather than a t interval, when a condition has fewer than two replicate values, or when the plugin stored no replicate values. |
| `p` | `float \| None` | Unadjusted p value of the two-sample test. |
| `p_adjusted` | `float \| None` | p value after the correction named by `correction`. `None` when the plugin stored no corrected value, which makes the comparison a description rather than a decision; the verdict then reads `no test recorded`. |
| `test` | `str` | `student_t`, `welch_t` or `tukey_hsd`. |
| `correction` | `str` | `BH` (Benjamini-Hochberg), `tukey_hsd`, or the configured post-hoc name. |
| `cohens_d` | `float \| None` | Standardised mean difference, oriented like `delta`: positive means `b` is larger. The framework's own `PairwiseResult.cohens_d` uses the opposite sign and is flipped here. |
| `hedges_g` | `float \| None` | Small-sample-corrected standardised mean difference, oriented like `cohens_d`, when the plugin reports one. Otherwise `None`. |
| `direction` | `str` | The plugin's own direction word, for example `increased`. |
| `significant` | `bool` | Whether `p_adjusted` cleared the configured alpha (0.05 by default). Always `False` when `testable` is `False`. |
| `testable` | `bool` | `False` when a condition has fewer than two replicates, which makes the test undefined rather than non-significant. |

## ProtocolProvenance

| Field | Type | Meaning |
|---|---|---|
| `polyzymd_version` | `str` | Version of the package that ran the protocol. |
| `mdanalysis_version` | `str \| None` | Installed MDAnalysis version, `None` when it cannot be determined. |
| `config_hashes` | `dict[str, str]` | SHA-256 of each simulation config file, keyed by condition label. |
| `settings_fingerprint` | `str \| None` | Fingerprint of the resolved plugin settings, the same one the aggregate cache is validated against. |
| `output_paths` | `dict[str, str]` | `comparison_result` is the cached comparison JSON; `figures` is the directory holding the generated plots. Either may be absent. |

## Verdict vocabulary

The first words of a verdict sentence come from a fixed set, so a caller can
branch on them without parsing the rest.

| Word | Meaning |
|---|---|
| `larger` | `b` differs from `a` after correction and `delta` is positive |
| `smaller` | `b` differs from `a` after correction and `delta` is negative |
| `no significant difference` | the test ran and `p_adjusted` did not clear alpha |
| `no test recorded` | the plugin stored no multiplicity-corrected p value, so the row describes a difference without deciding it |
| `changed` | the difference is significant but the two means are equal at the stored precision |
| `not testable` | a condition has fewer than two replicates, so no test exists |

A single-condition report has no comparison, and its one sentence states the
label, the metric, the mean with its unit, the interval and the replicate
count.

## See also

- {doc}`cli_reference` for `polyzymd analyze` and the `agent` format.
- {doc}`../how_to/analysis_agent_protocol` for the recipe.
- {doc}`../explanation/analysis_entry_points` for when to use this instead of
  `compare run`.
