---
name: livecoms-check
description: Check analysis code against the LiveCoMS sampling and uncertainty best practices before committing, covering the replicate as sampling unit, required result fields, citations, and known-answer expectations.
---

# LiveCoMS check

Run this before committing any change under `src/polyzymd/analyses/`, and use it
as the rubric when reviewing someone else's analysis change. The standard is
Grossfield, Patrone, Roe, Schultz, Siderius and Zuckerman (2018). Best Practices
for Quantifying Sampling Quality and Uncertainty in Molecular Simulations.
Living Journal of Computational Molecular Science 1:5067.
doi:10.33011/livecoms.1.1.5067.

Read the diff, answer every question below, and fix what fails. Do not commit
with an unanswered item.

## Checklist items that apply to analysis code

1. Do not cherry-pick. Use all available data unless there is an objective
   reason not to, and apply any sampling metric uniformly to every simulation.
2. Remove an equilibration portion and analyse only the production part. The
   offset is one user-set value applied to every condition and replicate.
3. Run multiple simulations and compare results across them.
4. Estimate the number of statistically independent samples. Sequential frames
   are correlated and the correlation must be measured, not assumed.
5. Quantify uncertainty in each observable with a confidence interval, not a
   bare standard error.
6. Do the semiquantitative checks that can rule out sufficient sampling. Look at
   the time series, compare runs, split one run in half and treat the halves as
   two runs.
7. Report the whole uncertainty procedure. Say in the figure or table what the
   error bar means, for example "error bars are 95 percent confidence intervals
   over 3 replicates, production window t >= 10 ns".

## Rules this repository holds you to

- The replicate is the sampling unit for every cross-condition test and every
  condition-level uncertainty. Never run a hypothesis test on frames.
- Equilibration is applied uniformly. A diagnostic reports, it never selects
  data.
- Every reported uncertainty states what it is. Every metric carries a unit.
- Prefer the correct cited method over the convenient one.
- Write the failing test first, then fix. Known-answer tests live in
  `tests/analyses/scientific/`.
- Do not silently degrade. An invalid input raises a typed error from
  `polyzymd.analyses.exceptions`. It does not return zero.

## Fields every result must carry

A metric model without these is incomplete.

| Field | Meaning |
|-------|---------|
| `unit` | Physical unit of the value, for example `angstrom` or `dimensionless` |
| `ci95_low` | Lower bound of the 95 percent confidence interval |
| `ci95_high` | Upper bound of the 95 percent confidence interval |
| `ci_method` | How the interval was built, for example `t_interval` or `bootstrap_bca` |

Each artifact also carries an `uncertainty` block.

| Key | Meaning |
|-----|---------|
| `kind` | What the number is, for example `sem_across_replicates` or `ci95_t` |
| `n` | Number of sampling units the estimate rests on, replicates unless stated |
| `coverage` | Nominal coverage, for example `0.95` |

With n = 3 replicates the two-sided 95 percent coverage factor is 4.30, not
1.96. A half-width of one standard error covers about 61 percent, so reporting a
bare standard error as if it were a 95 percent interval is wrong by a factor of
four. State the factor you used.

## Citation convention

Every module that implements or adapts a published method carries a NumPy-style
`References` section in its module docstring, and the plugin `__init__.py` also
carries a one-line statement of the method it runs. Format each entry as
author, year, title, journal, volume and page, then the DOI.

```
References
----------
Chodera, J. D., Swope, W. C., Pitera, J. W., Seok, C., and Dill, K. A. (2007).
    Use of the weighted histogram analysis method for the analysis of simulated
    and parallel tempering simulations. Journal of Chemical Theory and
    Computation 3:26. doi:10.1021/ct0502864
```

The full list of what to cite is in section 10 of
`docs/planning/analyses_audit_2026-09-11.md`, and the rendered page is
`docs/source/explanation/references.md`. Do not cite a method the code does not
implement. A citation for block averaging or bootstrapping belongs in the code
only once that code exists.

## Known answers the tests must reproduce

If your change touches an estimator, check it against these before you commit.

| Input | Expected |
|-------|----------|
| White noise, iid normal | Statistical inefficiency g = 1, within 10 percent |
| AR(1) with coefficient phi | g = (1 + phi) / (1 - phi) |
| n iid normal draws with standard deviation sigma | Standard error of the mean = sigma / sqrt(n) |
| A point cloud and a rotated copy of it | RMSD after superposition below 1e-8 |
| One sphere of radius r with probe radius p | Solvent-accessible surface area = 4 pi (r + p) squared |

An estimator that reports g = 3 on white noise is broken. That was the headline
defect of the 2026-09-11 audit.

## Before you commit

- Ruff and black are clean.
- `tests/analyses/scientific` passes.
- Every new number has a unit and a stated uncertainty.
- Every new method has a reference, and every reference names a method the code
  actually runs.
- The commit message is a conventional commit with an imperative subject near 50
  characters, plain prose in the body, and no attribution trailer.
