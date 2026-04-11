# Post-Hoc Testing Reference

```{contents}
:local:
:depth: 2
```

## Overview

PolyzyMD performs post-hoc pairwise comparisons automatically during
`polyzymd compare run`. Two methods are available: BH-corrected t-tests
(default) and Tukey's HSD. Both methods compute Cohen's d effect sizes and
percent-change for every pair.

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
  fdr_alpha: 0.05              # BH threshold (only used when posthoc_method is ttest_bh)
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `posthoc_method` | `"ttest_bh"` or `"tukey_hsd"` | `"ttest_bh"` | Selects which post-hoc procedure to use for pairwise comparisons. |
| `ttest_method` | `"student"` or `"welch"` | `"student"` | Controls the variance assumption for the two-sample t-test. Only used when `posthoc_method` is `"ttest_bh"`. |
| `fdr_alpha` | float (0, 1] | `0.05` | Significance threshold for Benjamini-Hochberg correction. Also used as the alpha threshold for Tukey's HSD when that method is selected. |

---

## Method Details

### BH-Corrected t-Tests (`ttest_bh`)

- Runs independent two-sample t-tests for each pair of conditions.
- `ttest_method: "student"` assumes equal variances (`scipy.stats.ttest_ind` with `equal_var=True`).
- `ttest_method: "welch"` relaxes that assumption (`equal_var=False`).
- Raw p-values from all pairs across all metrics are collected into a single family, then adjusted via the Benjamini-Hochberg step-up procedure.
- A pair is significant when `p_adj <= fdr_alpha`.
- Cohen's d is computed for effect size, with interpretation labels: `"negligible"` (|d| < 0.2), `"small"` (0.2 <= |d| < 0.5), `"medium"` (0.5 <= |d| < 0.8), `"large"` (|d| >= 0.8).

```{note}
When a `control` label is set in `comparison.yaml`, only control-vs-treatment pairs are tested. Otherwise, all unique pairs are tested.
```

### Tukey's HSD (`tukey_hsd`)

- Tests all pairs simultaneously using `scipy.stats.tukey_hsd`.
- Controls the family-wise error rate (FWER) rather than FDR.
- A pair is significant when the Tukey-adjusted `p_value <= 0.05` (uses `fdr_alpha` as the threshold).
- Best with balanced designs (equal replicates per condition).
- Cohen's d is still computed for each pair.
- The `t_statistic` field is set to `NaN` for Tukey results since the test does not produce a t-statistic.

---

## ANOVA

A one-way ANOVA (`scipy.stats.f_oneway`) is always run alongside post-hoc tests when there are 3 or more conditions. It tests whether at least one condition's mean differs from the others.

| Field | Type | Description |
|-------|------|-------------|
| `f_statistic` | float | F-statistic from the one-way ANOVA. |
| `p_value` | float | P-value for the omnibus test. |
| `significant` | bool | Whether `p_value < 0.05`. |

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
| `p_value_adjusted` | float or null | Adjusted p-value (BH for `ttest_bh`; Tukey-adjusted for `tukey_hsd`). `null` when not available. |
| `posthoc_method` | str | `"ttest_bh"` or `"tukey_hsd"`. |
| `cohens_d` | float | Effect size (positive = `condition_a` mean > `condition_b` mean). |
| `effect_size_interpretation` | str | `"negligible"`, `"small"`, `"medium"`, or `"large"`. |
| `direction` | str | Interpreted direction of change (e.g. `"higher"`, `"lower"`, or `"undetermined"`). |
| `significant` | bool | Whether `p_adj <= alpha` (uses adjusted p-value when available, raw otherwise). |
| `percent_change` | float | Percent change from `condition_a` to `condition_b`. |

### Comparison-level fields

| Field | Type | Description |
|-------|------|-------------|
| `fdr_alpha` | float | The alpha threshold used for BH correction. |
| `ttest_method` | str | `"student"` or `"welch"`. |
| `posthoc_method` | str | `"ttest_bh"` or `"tukey_hsd"`. |

---

## CLI Significance Markers

The CLI formatter annotates pairwise rows with significance markers:

| Marker | Meaning |
|--------|---------|
| `*` | `p_adj <= 0.05` |
| `**` | `p_adj <= 0.01` |
| `***` | `p_adj <= 0.001` |

```{note}
The multi-level markers (`**`, `***`) are used by the contacts, binding free energy, and polymer affinity formatters. The default scalar formatter uses a single `*` for any significant result.
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
| Equal values across all replicates in both groups | p-value = 1.0, Cohen's d = 0.0 (`"negligible"`). |
| Single condition | No pairwise tests are generated; ANOVA is skipped. |
| Two conditions | Pairwise tests run normally; ANOVA is skipped (requires >= 3 conditions). |
| Tukey HSD with fewer than 2 groups or fewer than 2 observations per group | Raises `ValueError`. |
| Zero control mean | Percent change returns `inf` or `-inf`; `NaN` if both means are non-finite. |

---

## See Also

- {doc}`comparison_yaml` -- full `comparison.yaml` schema reference
- {doc}`analysis_comparison_reference` -- comparison CLI commands and plugin summary
- {doc}`../explanation/analysis_statistics_best_practices` -- autocorrelation, FDR concepts, and interpretation guidance
