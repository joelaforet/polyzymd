# Statistics Best Practices for MD Analysis

This page explains how to interpret statistical uncertainty in PolyzyMD
analysis results. It follows the spirit of the LiveCoMS best-practices
recommendations for molecular-simulation uncertainty quantification
(Grossfield et al., 2018), while staying cautious about what finite trajectories
can prove.

```{contents} On This Page
:local:
:depth: 2
```

## Why MD statistics need extra care

Molecular dynamics trajectories are correlated in time. Consecutive frames are
usually similar because they are snapshots from one continuous dynamical path,
not independent experimental repeats. This does not make trajectory averages
useless, but it does mean that the number of saved frames is not the number of
independent observations.

Autocorrelation mainly affects interpretation in three ways:

- It reduces the effective amount of independent information in a trajectory.
- It makes uncertainty estimates and convergence checks more difficult.
- It can make finite-trajectory fluctuation estimates unreliable, especially
  when dynamics are slow, rare events are missed, or the trajectory is not
  stationary.

```{warning}
Naive uncertainty estimates based on the raw frame count can be much too
confident. Treat statistical significance from short, correlated, or
non-stationary trajectories as provisional unless independent replicates and
convergence diagnostics support the conclusion.
```

PolyzyMD analyses may account for correlation directly where a plugin supports
it. A generalized metric-type system for automatic correlation handling is a
planned design direction, not a universal implemented contract.

## Key concepts

### Autocorrelation

An autocorrelation function measures how much an observable resembles itself at
later time lags:

$$C(\tau) = \frac{\langle (x(t) - \mu)(x(t+\tau) - \mu) \rangle}{\sigma^2}$$

At zero lag, $C(0) = 1$. As the observable decorrelates, $C(\tau)$ approaches
zero. Slow decay means that many saved frames contain overlapping information.

### Statistical inefficiency and effective sample size

For a correlated time series, the effective number of independent samples is
often described using the statistical inefficiency $g$:

$$N_{\mathrm{eff}} = \frac{N}{g}$$

For a simple exponentially decaying correlation, a useful approximation is:

$$g \approx 1 + \frac{2\tau}{\Delta t}$$

where $\tau$ is an integrated correlation time and $\Delta t$ is the spacing
between saved observations. The exact estimate depends on the observable, the
trajectory length, and how the correlation function is truncated or blocked.

Spacing frames by about $2\tau$ is therefore best understood as a heuristic for
selecting approximately independent frames, not as a universal formula for the
true number of independent samples.

```{tip}
A trajectory with 10,000 saved frames may contain far fewer effective samples
for a slow observable. The relevant timescale is the decorrelation time of the
observable being analyzed, not the output frequency alone.
```

### Stationarity and convergence

Correlation corrections assume that the sampled data represent one stationary
distribution after equilibration. If the trajectory drifts between metastable
states, remains trapped in one basin, or contains large relaxation transients,
then a corrected SEM can still be misleading.

For this reason, convergence is not just a numerical property of an error bar.
It is a scientific judgment about whether the simulation sampled the states
that matter for the question being asked.

## Means and fluctuations need different interpretation

### Averages

For observables such as mean distance, mean RMSD, contact fraction, or average
energy, all frames can contribute information about the mean if the trajectory
is equilibrated and representative. Autocorrelation mainly reduces precision:
the mean may be estimated from all frames, but the uncertainty should reflect
$N_{\mathrm{eff}}$ or an equivalent blocking/replicate-based estimate rather
than the raw frame count.

Conceptually:

```text
mean = average(all equilibrated frames)
uncertainty = spread adjusted for correlation or estimated across replicates
```

This is why a mean computed from a dense trajectory can be useful while its
naive frame-based SEM is not.

### Fluctuations

For observables that are themselves fluctuation measures, such as RMSF, energy
variance, compressibility-like quantities, or variance-derived order
parameters, finite correlated trajectories are harder to interpret. The issue is
not that every correlated variance estimate is universally biased downward in a
simple, predictable way. The problem is broader: finite, correlated, and
possibly non-stationary trajectories may not sample the full distribution of
motions that define the fluctuation.

Conservative strategies include:

- estimating fluctuations from approximately independent blocks or subsamples;
- comparing results across independent replicates;
- checking whether the estimate changes with trajectory length or block size;
- reporting caveats when slow motions are undersampled.

Subsampling frames at roughly $2\tau$ spacing can help avoid treating adjacent
frames as independent, but it does not rescue a trajectory that never visited
important conformational states.

## How to interpret PolyzyMD analysis results

PolyzyMD's analysis plugins are responsible for choosing statistically
appropriate handling for the quantities they report. Current plugins may use
plugin-specific logic for correlation-aware uncertainty, replicate aggregation,
or conservative subsampling. The details can differ because contact fractions,
distances, RMSF, secondary structure, and other observables answer different
scientific questions.

When reading output, focus on these questions:

Mean-like observables
: Does the uncertainty account for correlation or independent replicates rather
  than the raw number of frames?

Fluctuation-like observables
: Was the trajectory long enough to sample the motions that define the
  fluctuation? Do replicate or block estimates agree?

Condition comparisons
: Are the compared conditions based on comparable simulation protocols,
  equilibration choices, and replicate counts?

Significance markers
: Were multiple comparisons corrected, and are the effect sizes meaningful in
  addition to being statistically significant?

For implementation details of shared analysis utilities, see
{doc}`../api/analyses_shared`.

## Replicates and condition-level uncertainty

Independent replicates are often the most interpretable basis for
condition-level uncertainty. A replicate-level SEM estimates how much replicate
means vary across independent starts, seeds, or prepared systems.

This is preferred when the replicates are:

- genuinely independent;
- generated with comparable protocols;
- analyzed after appropriate equilibration removal;
- long enough to sample the relevant states for each observable.

Replicate SEM is not magically robust to bad sampling. If all replicates are
trapped in the same metastable basin, or if different replicates have not
reached comparable stationary behavior, the between-replicate spread may still
understate or misrepresent uncertainty.

## Practical interpretation principles

### Compare timescales before trusting uncertainty

Ask whether the process of interest decorrelates many times within the analyzed
window. If a loop rearrangement, contact transition, or polymer conformation
changes on a timescale comparable to the whole production trajectory, a small
error bar is not strong evidence of convergence.

### Prefer effect sizes with uncertainty over p-values alone

A statistically significant difference can be scientifically small, while a
large apparent effect can be uncertain if replicate counts are low. Interpret
PolyzyMD comparison output by considering the direction, magnitude, confidence
intervals or SEMs, and the physical plausibility of the change.

### Treat diagnostics as evidence, not proof

Autocorrelation estimates, block averages, and replicate SEMs are diagnostics.
They summarize available sampling; they cannot reveal unsampled states that the
trajectory never visited. Use them alongside structural inspection and domain
knowledge.

## Multiple comparison correction

### Why it matters

When comparing *N* conditions pairwise, the number of hypothesis tests grows as
*N*(*N*−1)/2. With 5 conditions this is already 10 tests; with 10 conditions,
45 tests. Without correction, the probability of at least one false positive
increases rapidly, even when every null hypothesis is true.

### Benjamini-Hochberg FDR control

When `posthoc_method` is `"ttest_bh"`, PolyzyMD uses the
Benjamini-Hochberg procedure to control the false discovery rate: the expected
proportion of false positives among rejected hypotheses. This is less
conservative than controlling the probability of any false positive, and is
often appropriate when many related comparisons are screened together.

Conceptually, the procedure ranks p-values, compares them to rank-dependent
thresholds, and marks discoveries only up to the largest rank that satisfies the
threshold. Adjusted p-values and CLI significance markers are based on this
correction rather than on raw p-values alone.

### Tukey HSD as an alternative

When `posthoc_method` is `"tukey_hsd"`, PolyzyMD uses Tukey's Honestly
Significant Difference test. Tukey HSD compares all pairs simultaneously and
controls the family-wise error rate. It is most appropriate for balanced designs
with similar replicate counts and approximately equal variance across
conditions.

Choose the post-hoc method based on the scientific question:

- Use BH-corrected t-tests when specific pairs are of interest or sample sizes
  are heterogeneous.
- Use Tukey HSD when all pairwise contrasts are part of one balanced comparison
  family.

### Within-trajectory uncertainty and multiple testing are separate

Autocorrelation handling and multiple-comparison correction address different
problems. Autocorrelation affects the uncertainty of an observable within or
across trajectories. Multiple-comparison correction controls how many false
positives are expected when many hypotheses are tested.

Both matter. A corrected p-value does not fix poor sampling, and a careful SEM
does not control the false-positive rate across many condition pairs.

For field-level details on post-hoc methods, configuration keys, and output
fields, see {doc}`../reference/posthoc_testing`.

## References

This explanation follows guidance from:

> **Grossfield, A., et al.** (2018). "Best Practices for Quantification of
> Uncertainty and Sampling Quality in Molecular Simulations."
> *Living Journal of Computational Molecular Science*, 1(1), 5067.
> [DOI: 10.33011/livecoms.1.1.5067](https://doi.org/10.33011/livecoms.1.1.5067)

Useful background includes:

> **Chodera, J. D., et al.** (2007). "Use of the Weighted Histogram Analysis
> Method for the Analysis of Simulated and Parallel Tempering Simulations."
> *Journal of Chemical Theory and Computation*, 3(1), 26-41.

> **Flyvbjerg, H., & Petersen, H. G.** (1989). "Error estimates on averages
> of correlated data." *Journal of Chemical Physics*, 91(1), 461-466.
