# Establishing Convergence in MD Simulations

Understanding when a molecular dynamics simulation has converged — and what
convergence means in practice — is essential for drawing reliable conclusions.

```{versionadded} 1.3.0
Automated convergence detection was added alongside the RMSD analysis plugin.
```

## What Is Convergence in MD?

A simulation has *converged* when its observable of interest has stopped
drifting and is sampling from a stationary distribution. In the RMSD context,
this means the protein's deviation from a reference structure has settled into
a fluctuating plateau rather than continuing to increase or decrease.

Convergence is **not the same as equilibration**. Equilibration refers to the
initial transient period after simulation launch, during which the system
relaxes from its starting configuration. Convergence refers to the state of
the production region itself — whether the trajectory has sampled long enough
that running averages are stable and the statistical properties of the
observable are no longer evolving.

## Why It Matters

Conclusions drawn from non-converged simulations are unreliable. If the RMSD
is still drifting upward, the mean RMSD and its uncertainty will change
depending on how much data you include. Effect sizes between conditions may
appear significant or insignificant depending on where you truncate the
timeseries.

Grossfield et al. (2019) emphasize that quantifying uncertainty requires
sampling from a stationary distribution. If the distribution itself is still
evolving — as it is during a drift — standard error estimates understate the
true uncertainty.

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
available. Automated methods can miss patterns — such as oscillations between
two metastable states — that are obvious to a trained eye.

## PolyzyMD's Sliding-Window Approach

PolyzyMD implements a sliding-window slope heuristic for convergence detection.
The algorithm operates on any 1D timeseries (typically RMSD vs. time) and
proceeds as follows:

1. **Divide the timeseries into overlapping windows.** Each window spans a
   fixed duration (default: 15 ns) and successive windows are offset by a
   step size (default: 5 ns).

2. **Compute the mean observable in each window.** This smooths out
   frame-to-frame noise while preserving slow drift.

3. **Estimate the slope between successive window means.** The slope captures
   the rate of change in the smoothed signal.

4. **Check for sustained low slope.** If the absolute slope remains below a
   threshold for a sustained duration, the timeseries is declared converged.
   The convergence time is the start of the first sustained plateau.

This approach is robust to brief transient excursions — a single window with a
slightly elevated slope does not reset the clock unless it exceeds the
threshold. The requirement for *sustained* low slope avoids false positives
from momentary pauses in an otherwise drifting trajectory.

## Default Parameters and When to Tune Them

| Parameter | Default | Description |
|-----------|---------|-------------|
| `convergence_window_size_ns` | 15.0 | Width of each averaging window (ns) |
| `convergence_step_size_ns` | 5.0 | Step between successive window starts (ns) |
| `convergence_slope_threshold` | 0.0005 | Maximum absolute slope (Å/ns) to qualify as "flat" |
| `convergence_sustained_for_ns` | 15.0 | Required duration below threshold to declare convergence (ns) |

**Guidance for tuning:**

- **Very long simulations (> 500 ns):** Increase `convergence_window_size_ns`
  and `convergence_sustained_for_ns` proportionally. A 15 ns window in a
  1 μs trajectory may be too sensitive to short-timescale fluctuations.

- **High-precision comparisons:** Decrease `convergence_slope_threshold` to
  require a flatter plateau before declaring convergence.

- **Noisy observables:** Increase `convergence_slope_threshold` to tolerate
  larger fluctuations. Polymer RMSD, for example, tends to be noisier than
  protein backbone RMSD.

- **Short simulations (< 50 ns):** Decrease `convergence_window_size_ns` and
  `convergence_sustained_for_ns` so the algorithm has enough data to assess.
  Be aware that shorter windows reduce the reliability of the assessment.

## Limitations

The sliding-window heuristic is a practical diagnostic, not a theoretical
proof of convergence. Important limitations include:

- **Not a proof of ergodic sampling.** A flat RMSD timeseries does not
  guarantee that the simulation has explored all relevant conformational
  states. The system could be trapped in a metastable basin.

- **Insensitive to slow conformational drift.** If the drift rate is below the
  slope threshold, the algorithm will declare convergence even though the
  observable is still changing — just slowly.

- **Parameter dependent.** The four tunable parameters (window size, step size,
  slope threshold, sustained duration) introduce subjective choices. Different
  parameter values can yield different convergence conclusions for the same
  trajectory.

- **Single-observable limitation.** Convergence in RMSD does not imply
  convergence in other observables (e.g., hydrogen bond occupancy, active site
  geometry). Different metrics may converge at different rates.

- **Multiple independent replicates remain essential.** Even with a converged
  RMSD timeseries in a single replicate, independent replicates are needed to
  quantify system-level variability and confirm reproducibility.

## Relationship to Equilibration Time

Convergence detection and equilibration time (`--eq-time`) address related but
distinct concerns:

- **Equilibration** removes transient artifacts from the start of the
  simulation — the period during which the system relaxes from its initial
  configuration. The equilibration time is a fixed cutoff applied *before*
  analysis begins.

- **Convergence** confirms that the production region (after equilibration) is
  stationary. It answers whether the remaining data is reliable for computing
  time-averaged properties.

Convergence detection can help *inform* the choice of `--eq-time`. If the
convergence diagnostic reports a convergence time of 12 ns but you set
`--eq-time 5ns`, the first 7 ns of your "production" data may still contain
drift. Conversely, if convergence is detected at 5 ns and you set
`--eq-time 20ns`, you may be discarding usable data.

In practice, the two are complementary: set `--eq-time` conservatively based
on visual inspection of the RMSD timeseries, then use convergence detection as
an independent consistency check.

## Beyond RMSD

The convergence utility (`analyses/shared/convergence.py`) accepts any 1D
timeseries — it is not specific to RMSD. The `find_convergence_time()`
function takes arrays of time values and signal values, making it applicable to
radius of gyration, solvent-accessible surface area, or any other scalar
observable that should plateau when the system reaches equilibrium.

In future PolyzyMD versions, convergence detection may be integrated into
additional analysis plugins. The algorithmic foundation is already
general-purpose.

## References

**Grossfield A, Patrone PN, Roe DR, Schultz AJ, Siderius DW, Zuckerman DM.**
(2019) "Best Practices for Quantification of Uncertainty and Sampling Quality
in Molecular Simulations." *Living Journal of Computational Molecular Science*
1(1):5067. [doi:10.33011/livecoms.1.1.5067](https://doi.org/10.33011/livecoms.1.1.5067)

The authoritative guide for uncertainty quantification in MD, including
discussion of convergence assessment, autocorrelation, and effective sample
sizes.

**Knapp B, Frantal S, Greshake B, Schwarz R, et al.** (2018) "Is an Intuitive
Convergence Definition of Molecular Dynamics Simulations Solely Based on the
Root Mean Square Deviation Possible?" *Journal of Computational Biology*
25:1069-1077.

Analysis of RMSD-based convergence criteria and their reliability, motivating
the use of sliding-window and sustained-plateau approaches over simple visual
inspection alone.

## See Also

- {doc}`/how_to/analysis_rmsd_quickstart` — RMSD quick start with convergence configuration
- {doc}`/explanation/analysis_rmsd_best_practices` — RMSD interpretation and best practices
- {doc}`/explanation/analysis_statistics_best_practices` — Statistical foundations for MD analysis
