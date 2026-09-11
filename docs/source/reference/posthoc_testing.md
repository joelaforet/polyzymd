# Post-Hoc Testing Reference

```{contents}
:local:
:depth: 2
```

## Overview

PolyzyMD performs post-hoc pairwise comparisons automatically during
`polyzymd compare run`. Two methods are available: BH-corrected t-tests
(default) and Tukey's HSD. Both methods compute effect sizes and
percent-change for every pair.

Every analysis plugin uses the same machinery. `ttest_method`,
`posthoc_method` and `fdr_alpha` from `comparison.yaml` reach every plugin's
`compare` step, and the correction family is defined in one place,
`polyzymd.analyses.shared.inferential_statistics.apply_family_correction`.

---

## Available Methods

| Method | `posthoc_method` value | When to use | Assumptions |
|--------|------------------------|-------------|-------------|
| BH-corrected t-tests | `"ttest_bh"` | Specific pairs of interest or heterogeneous sample sizes. **Default.** | Independence between pairs; equal variance (Student) or relaxed (Welch). |
| Tukey's HSD | `"tukey_hsd"` | All conditions are equally important; balanced design preferred. | Equal variance; equal (or similar) sample sizes across conditions. |

---

## Configuration

Set post-hoc options in the `defaults:` block of `comparison.yaml`:

```yaml
defaults:
  posthoc_method: "ttest_bh"   # or "tukey_hsd"
  ttest_method: "student"      # or "welch" (only used when posthoc_method is ttest_bh)
  fdr_alpha: 0.05              # significance threshold (used by both ttest_bh and tukey_hsd)
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `posthoc_method` | `"ttest_bh"` or `"tukey_hsd"` | `"ttest_bh"` | Selects which post-hoc procedure to use for pairwise comparisons. |
| `ttest_method` | `"student"` or `"welch"` | `"student"` | Controls the variance assumption for the two-sample t-test. Only used when `posthoc_method` is `"ttest_bh"`. |
| `fdr_alpha` | float (0, 1] | `0.05` | Significance threshold for pairwise comparisons and ANOVA. Used as the BH false-discovery-rate threshold when `posthoc_method` is `"ttest_bh"`, and as the family-wise alpha threshold when `posthoc_method` is `"tukey_hsd"`. |

---

## Method Details

### BH-Corrected t-Tests (`ttest_bh`)

- Runs independent two-sample t-tests for each pair of conditions.
- `ttest_method: "student"` assumes equal variances (`scipy.stats.ttest_ind` with `equal_var=True`).
- `ttest_method: "welch"` relaxes that assumption (`equal_var=False`).
- Raw p-values from all pairs across all metrics are collected into a single family, then adjusted via the Benjamini-Hochberg step-up procedure.
- A pair is significant when `p_adj <= fdr_alpha`.
- Cohen's d and Hedges' g are computed for every pair. The interpretation label is `"negligible"` (|d| < 0.2), `"small"` (0.2 <= |d| < 0.5), `"medium"` (0.5 <= |d| < 0.8) or `"large"` (|d| >= 0.8), and is `null` when `n1 + n2 < 10`.

#### The correction family

One analysis run is one family. The family holds every pairwise test the run
produced: every metric the plugin compares, and every condition pair, in one
Benjamini-Hochberg step-up procedure. Running RMSD with two named runs over
four conditions against a control therefore corrects 6 tests together, not
3 and 3.

The family never spans separate `polyzymd compare run` invocations, and it
never includes the ANOVA.

#### Effect size

| Field | Meaning |
|-------|---------|
| `cohens_d` | Mean difference divided by the pooled standard deviation. |
| `hedges_g` | `cohens_d` multiplied by the Hedges (1981) correction `J = 1 - 3 / (4 * (n1 + n2) - 9)`. With three replicates per condition, `J` is about 0.80. |
| `effect_size_interpretation` | The Cohen (1988) adjective, or `null` when `n1 + n2 < 10`. |

Quote `hedges_g` at the replicate counts molecular dynamics produces. Cohen's
d is biased upward at small `n`, and its standard error is of order one, so
the adjective is withheld below ten combined replicates rather than reporting
noise as "large".

```{note}
When a `control` label is set in `comparison.yaml`, only control-vs-treatment pairs are tested. Otherwise, all unique pairs are tested.
```

### Tukey's HSD (`tukey_hsd`)

- Tests all pairs simultaneously using `scipy.stats.tukey_hsd`.
- Controls the family-wise error rate (FWER) rather than FDR.
- A pair is significant when the Tukey-adjusted `p_value <= fdr_alpha`.
- `p_value_adjusted` mirrors `p_value` for Tukey results (Tukey p-values are already family-wise corrected).
- Best with balanced designs (equal replicates per condition).
- Cohen's d is still computed for each pair.
- The `t_statistic` field is set to `NaN` for Tukey results since the test does not produce a t-statistic.

---

## ANOVA

A one-way ANOVA (`scipy.stats.f_oneway`) is always run alongside post-hoc tests when there are 3 or more conditions. It tests whether at least one condition's mean differs from the others.

The ANOVA is an omnibus test and is reported uncorrected. It is not a member
of the pairwise family, its `p_value_adjusted` is always `null`, and no
pairwise test is gated on it. Its `significant` flag compares the raw
p-value with `fdr_alpha`. This rule is the same in every plugin.

| Field | Type | Description |
|-------|------|-------------|
| `f_statistic` | float | F-statistic from the one-way ANOVA. |
| `p_value` | float | Raw p-value for the omnibus test. Never adjusted. |
| `p_value_adjusted` | null | Always `null`. The ANOVA is outside the correction family. |
| `significant` | bool | Whether the raw `p_value <= fdr_alpha`. |

ANOVA does **not** determine which pairs differ -- that is the role of post-hoc tests. ANOVA is skipped when fewer than 3 conditions are present, and returns `NaN` statistics if any group has fewer than 2 observations.

---

## Output Fields

These fields appear in comparison JSON files and are used by the CLI formatter.

### Pairwise result fields

| Field | Type | Description |
|-------|------|-------------|
| `condition_a` | str | Label of first condition (typically control). |
| `condition_b` | str | Label of second condition (typically treatment). |
| `metric` | str | Name of the metric being compared. |
| `t_statistic` | float | T-test statistic. `NaN` for Tukey HSD results. |
| `p_value` | float | Raw p-value from the pairwise test. |
| `p_value_adjusted` | float or null | Adjusted p-value (BH-adjusted for `ttest_bh`; mirrors `p_value` for `tukey_hsd` since Tukey p-values are already family-wise corrected). `null` when not available. |
| `posthoc_method` | str | `"ttest_bh"` or `"tukey_hsd"`. |
| `cohens_d` | float | Effect size (positive = `condition_a` mean > `condition_b` mean). |
| `hedges_g` | float or null | `cohens_d` after the Hedges (1981) small-sample correction. |
| `effect_size_interpretation` | str or null | `"negligible"`, `"small"`, `"medium"` or `"large"`. `null` when `n1 + n2 < 10`. |
| `direction` | str | Direction of change, for example `"stabilizing"`, `"increased"` or `"closer"` depending on the plugin. Reads `"no significant change"` whenever `significant` is false. |
| `significant` | bool | Whether `p_adj <= alpha` (uses adjusted p-value when available, raw otherwise). |
| `percent_change` | float | Percent change from `condition_a` to `condition_b`. |

### Comparison-level fields

| Field | Type | Description |
|-------|------|-------------|
| `fdr_alpha` | float | The alpha threshold used for significance (BH FDR for `ttest_bh`; FWER for `tukey_hsd`). Also used as the ANOVA significance threshold. |
| `ttest_method` | str | `"student"` or `"welch"`. |
| `posthoc_method` | str | `"ttest_bh"` or `"tukey_hsd"`. |

---

## CLI Significance Markers

The CLI formatter annotates pairwise rows with significance markers:

| Marker | Meaning |
|--------|---------|
| `*` | `p_adj <= fdr_alpha` (default 0.05) |
| `**` | `p_adj <= 0.01` |
| `***` | `p_adj <= 0.001` |

```{note}
The contacts formatter uses multi-level markers (`**`, `***`). The default
scalar formatter uses a single `*` for any significant result.
```

Some plugins also use:

| Marker | Meaning |
|--------|---------|
| `†` (dagger) | Cohen's d meets the `min_effect_size` threshold (practical significance). Currently used by the contacts plugin. |

---

## Edge Cases

| Scenario | Behavior |
|----------|----------|
| Fewer than 2 replicates in a group | t-test returns `NaN` for t-statistic and p-value; Cohen's d returns `NaN`. |
| Equal values across all replicates in both groups | p-value = 1.0, Cohen's d = 0.0, interpretation `null` at `n1 + n2 < 10`. |
| Fewer than ten combined replicates | `effect_size_interpretation` is `null`; `cohens_d` and `hedges_g` are still reported. |
| Comparison not significant after correction | `direction` reads `"no significant change"`. |
| Single condition | No pairwise tests are generated; ANOVA is skipped. |
| Two conditions | Pairwise tests run normally; ANOVA is skipped (requires >= 3 conditions). |
| Tukey HSD with fewer than 2 groups or fewer than 2 observations per group | Returns empty results (no pairs generated). |
| Zero control mean | Percent change returns `inf` or `-inf`; `NaN` if both means are non-finite. |

---

## See Also

- {doc}`comparison_yaml` -- full `comparison.yaml` schema reference
- {doc}`analysis_comparison_reference` -- comparison CLI commands and plugin summary
- {doc}`../explanation/analysis_statistics_best_practices` -- autocorrelation, FDR concepts, and interpretation guidance

## References

- Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B*, 57(1), 289-300. doi:10.1111/j.2517-6161.1995.tb02031.x
- Cohen, J. (1988). *Statistical Power Analysis for the Behavioral Sciences*, 2nd edition. Lawrence Erlbaum Associates.
- Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect size and related estimators. *Journal of Educational Statistics*, 6(2), 107-128. doi:10.3102/10769986006002107
- Tukey, J. W. (1949). Comparing individual means in the analysis of variance. *Biometrics*, 5(2), 99-114. doi:10.2307/3001913
- Welch, B. L. (1947). The generalization of "Student's" problem when several different population variances are involved. *Biometrika*, 34(1-2), 28-35. doi:10.1093/biomet/34.1-2.28
