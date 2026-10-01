# Comparison Tests Reference

```{contents}
:local:
:depth: 2
```

## Overview

`polyzymd analyze` with two or more `-c`, and `ReplicateValues.compare()` on
the study API, compare every condition with the control, the first `-c` or
the `control=` argument. The replicate is the sampling unit: each condition's
sample is its one value per replicate. Each comparison row is one two-sample
t test of `mean(b) - mean(a)`, and the p values of one call are corrected
together with the Benjamini-Hochberg procedure. The rows are the `pairwise`
entries of the {doc}`ProtocolReport <analysis_protocol_report>`.

| Setting | `polyzymd analyze` | `ReplicateValues.compare()` |
|---|---|---|
| Control | first `-c` | `control=`, default the first condition of the study |
| Test | Welch's t test | `test="welch"` (default) or `test="student"` |
| Correction | Benjamini-Hochberg over every tested row of the report | Benjamini-Hochberg over every tested row of the call |
| Significance threshold | adjusted p at most 0.05 | adjusted p at most 0.05 |

No omnibus ANOVA is run and no Tukey HSD is offered: each condition is
compared with the control only, not with every other condition.

---

## The test

- `welch_t` runs `scipy.stats.ttest_ind(b, a, equal_var=False)`, which does
  not assume equal variances; `student_t` runs it with `equal_var=True`.
- `delta_ci95` is the 95 percent interval on `mean(b) - mean(a)` from the
  same `ttest_ind` call: a pooled variance for Student's t, separate
  variances with Welch-Satterthwaite degrees of freedom for Welch's t. It
  carries no multiplicity correction, so a row can be non-significant after
  the correction while its interval excludes zero.
- A row is **not testable** when a condition has fewer than two replicates,
  or when both conditions have the same value in every replicate. It then has
  no `p` and takes no part in the correction, and the verdict reads
  `not testable`.

---

## The correction family

One call is one family. For one number per replicate, the family is the
conditions compared with the control: four treated conditions against one
control correct four tests together. For labelled values, such as a
per-residue profile, the family is every tested label of every compared
condition, because scanning a profile for the labels that changed is a search
over that set. `family_size` on each row records the size of its family.
Labels passed as `untested=` are summarised but left out of the tests and the
family.

The family never spans separate commands or calls. For why the boundary is
drawn there, see {doc}`../explanation/analysis_statistics_best_practices`.

---

## Effect size and direction

| Field | Meaning |
|-------|---------|
| `cohens_d` | Mean difference divided by the pooled standard deviation, oriented like `delta`: positive when condition `b` is larger. |
| `hedges_g` | `cohens_d` multiplied by the Hedges (1981) correction `J = 1 - 3 / (4 * (n1 + n2) - 9)`. With three replicates per condition, `J` is about 0.80. |
| `direction` | `increased` or `decreased` when the row is significant, otherwise `no significant change`. |

For which number to quote, see
{doc}`../explanation/analysis_statistics_best_practices`.

---

## Edge cases

| Scenario | Behavior |
|----------|----------|
| One replicate in a condition | The row is not testable: no `p`, no `delta_ci95`; the report warns that the condition has no interval. |
| The same value in every replicate of both conditions | The row is not testable. |
| One condition | No comparison rows; the verdict summarises the condition. |
| A label in `untested` | Summarised in `conditions`, absent from `pairwise` and from the family. |

---

## See Also

- {doc}`analysis_protocol_report` -- every field of the report and its comparison rows
- {doc}`../explanation/analysis_api` -- `summary()` and `compare()` on the study API
- {doc}`../explanation/analysis_statistics_best_practices` -- autocorrelation, FDR concepts, and interpretation guidance

## References

- Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B*, 57(1), 289-300. doi:10.1111/j.2517-6161.1995.tb02031.x
- Benjamini, Y. (2010). Discovering the false discovery rate. *Journal of the Royal Statistical Society B*, 72(4), 405-416. doi:10.1111/j.1467-9868.2010.00746.x
- Hedges, L. V. (1981). Distribution theory for Glass's estimator of effect size and related estimators. *Journal of Educational Statistics*, 6(2), 107-128. doi:10.3102/10769986006002107
- Welch, B. L. (1947). The generalization of "Student's" problem when several different population variances are involved. *Biometrika*, 34(1-2), 28-35. doi:10.1093/biomet/34.1-2.28
