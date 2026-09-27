# Establishing Convergence in MD Simulations

Understanding when a molecular dynamics simulation has converged, and what
convergence means in practice, is essential for drawing reliable conclusions.

```{important}
Automated equilibration detection is a **diagnostic**, not proof that a
trajectory has converged, equilibrated, or sampled ergodically. Treat it as one
line of evidence alongside visual inspection, agreement among independent
replicates, uncertainty analysis, and scientific judgment about the system and
observable being studied. In PolyzyMD it never refuses a calculation, never
moves the equilibration window you set, and never changes a value or a
statistic.
```

```{versionadded} 1.3.0
The pymbar equilibration diagnostic was added with the study API.
```

## What Is Convergence in MD?

A simulation has *converged* when its observable of interest has stopped
drifting and is sampling from a stationary distribution. In the RMSD context,
this means the protein's deviation from a reference structure has settled into
a fluctuating plateau rather than continuing to increase or decrease.

Convergence is **not the same as equilibration**. Equilibration refers to the
initial transient period after simulation launch, during which the system
relaxes from its starting configuration. Convergence refers to the state of
the production region itself: whether the trajectory has sampled long enough
that running averages are stable and the statistical properties of the
observable are no longer evolving.

## Why It Matters

Conclusions drawn from non-converged simulations of an equilibrium property are
unreliable. If the RMSD is still drifting upward, the mean RMSD and its
uncertainty will change depending on how much data you include. Effect sizes
between conditions may appear significant or insignificant depending on where
you truncate the timeseries.

Grossfield et al. (2018) emphasize that quantifying uncertainty requires
sampling from a stationary distribution. If the distribution itself is still
evolving, as it is during a drift, standard error estimates understate the
true uncertainty.

```{note}
Not every study measures an equilibrium property. A thermal unfolding
simulation, for example, is not expected to reach a stationary state, and its
metrics are still worth computing and comparing between conditions as
replicate-level quantities. PolyzyMD computes every metric whatever the
diagnostic says; the diagnostic tells you how to read the result.
```

## Visual Indicators of Convergence

Before any automated tool, researchers assess convergence by inspecting
timeseries plots. Signs that a trajectory has converged include:

- **Plateau in the timeseries.** The observable fluctuates around a stable
  mean rather than trending in one direction.
- **Stable running averages.** A running mean computed over successively
  longer windows stops changing appreciably.
- **Decorrelation time stabilization.** The estimated autocorrelation time
  of the observable converges to a consistent value rather than growing.

These visual checks remain valuable even when automated diagnostics are
available. Automated methods can miss patterns, such as oscillations between
two metastable states, that are obvious to a trained eye.

## PolyzyMD's equilibration diagnostic

For every time series measured through the study API, PolyzyMD runs
`pymbar.timeseries.detect_equilibration` (Chodera 2016) on each replicate's
production frames, the frames left after the equilibration window you set.
The method tries a range of starting frames and picks the one that maximises
the number of effective samples in the remaining data. PolyzyMD tries about
100 evenly spaced starts, so the detected start is known to about 1 percent of
the series length.

The detected start is reported for every replicate, and the condition line of
the report ends with `eq_detected <ns>`, the latest detected start among the
condition's replicates. Two kinds of message can follow:

- **A replicate with at least 20 effective samples** whose detected start lies
  more than 10 percent of its production frames past the start of the window
  gets a warning such as
  `condition B replicate 2: pymbar detect_equilibration puts the start of the equilibrated region at 45.6 ns, production frame 457 of 2000, after the equilibration window; the window may be too short for it`.
  That replicate may still have been relaxing after your window.
- **Replicates with fewer than 20 effective samples** get one line per
  condition instead, such as
  `condition X: replicates 1, 3, 4 have fewer than 20 effective samples, so the start of an equilibrated region cannot be detected reliably; values and statistics are unaffected`.

The 20-sample threshold comes from how the method behaves on correlated data.
On stationary synthetic series with 100 or more effective samples, the
detected start passed 10 percent of the series in 1 to 7 percent of cases; with
about 10 effective samples it did so in 45 percent. On the stored RMSD series of
the LipA 363 K study, with 3 to 18 effective samples per replicate, the method
put the start in the last few percent of 28 of 30 replicates, on a short tail
whose statistical inefficiency is close to 1. On such data a late start says
the replicate holds too few independent samples to judge, not that the window
is too short.

Pass `--no-eq-check` to `polyzymd analyze`, or `detect_equilibration=False`
to `Timeseries.reduce`, to skip the diagnostic entirely. Values, intervals and
tests are identical with it on or off.

## Limitations

- **Not a proof of ergodic sampling.** A stationary-looking series does not
  guarantee that the simulation has explored all relevant conformational
  states. The system could be trapped in a metastable basin.

- **Needs many effective samples.** With few independent samples per
  replicate the detected start is unreliable, which is why PolyzyMD reports it
  without judging it below 20 effective samples.

- **Single-observable limitation.** Equilibration of RMSD does not imply
  equilibration of other observables (e.g., hydrogen bond occupancy, active
  site geometry). Different metrics may relax at different rates, which is why
  the diagnostic runs separately on every measured series.

- **Multiple independent replicates remain essential.** Every interval and
  test in PolyzyMD comes from the spread between independent replicates, and
  no single-replicate diagnostic replaces that.

## Relationship to Equilibration Time

The equilibration window and the diagnostic address related but distinct
concerns:

- **Equilibration** removes transient artifacts from the start of the
  simulation, the period during which the system relaxes from its initial
  configuration. The window is a fixed cutoff you set with `--eq` or
  `equilibration=`, applied to every replicate of every condition before
  analysis begins.

- **The diagnostic** asks whether the production-region observable, after the
  window, still looks like it is relaxing. It can help *inform* the choice of
  window, but PolyzyMD never changes the window itself.

In practice the two are complementary: set the window conservatively from
visual inspection of the time series, then read the diagnostic as an
independent consistency check.

## References

**Grossfield A, Patrone PN, Roe DR, Schultz AJ, Siderius DW, Zuckerman DM.**
(2018) "Best Practices for Quantification of Uncertainty and Sampling Quality
in Molecular Simulations." *Living Journal of Computational Molecular Science*
1(1):5067. [doi:10.33011/livecoms.1.1.5067](https://doi.org/10.33011/livecoms.1.1.5067)

The authoritative guide for uncertainty quantification in MD, including
discussion of convergence assessment, autocorrelation, and effective sample
sizes.

**Chodera JD.** (2016) "A Simple Method for Automated Equilibration Detection
in Molecular Simulations." *Journal of Chemical Theory and Computation*
12(4):1799-1805. [doi:10.1021/acs.jctc.5b00784](https://doi.org/10.1021/acs.jctc.5b00784)

The method behind `pymbar.timeseries.detect_equilibration`.

**Knapp B, Frantal S, Greshake B, Schwarz R, et al.** (2018) "Is an Intuitive
Convergence Definition of Molecular Dynamics Simulations Solely Based on the
Root Mean Square Deviation Possible?" *Journal of Computational Biology*
25:1069-1077.

Analysis of RMSD-based convergence criteria and their reliability.

## See Also

- {doc}`/how_to/analysis_rmsd_quickstart`: RMSD quick start
- {doc}`/explanation/analysis_rmsd_best_practices`: RMSD interpretation and best practices
- {doc}`/explanation/analysis_statistics_best_practices`: Statistical foundations for MD analysis
- [pymbar timeseries](https://pymbar.readthedocs.io/en/latest/timeseries.html):
  `detect_equilibration` and `statistical_inefficiency`
