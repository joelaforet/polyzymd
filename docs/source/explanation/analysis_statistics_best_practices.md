# Statistics Best Practices for MD Analysis

This guide explains the statistical methods used in PolyzyMD's analysis module,
following the **LiveCoMS Best Practices** for uncertainty quantification in
molecular simulations (Grossfield et al., 2018).

```{contents} On This Page
:local:
:depth: 2
```

## Why Statistics Matter in MD

Molecular dynamics trajectories are **highly correlated in time**. Unlike
independent experimental measurements, consecutive MD frames show the system
in nearly identical configurations. This correlation has profound implications
for how we analyze and report results.

```{warning}
Ignoring correlation in MD analysis leads to:
- **Overconfident uncertainties** (SEM too small by factors of 10-100×)
- **Biased variance estimates** (RMSF, heat capacity systematically wrong)
- **Irreproducible conclusions** (apparent significance that doesn't replicate)
```

PolyzyMD automatically handles these issues using methods from the
`polyzymd.analyses.shared.autocorrelation` module.

---

## Key Concepts

### Autocorrelation Function (ACF)

The autocorrelation function measures how correlated a signal is with itself
at different time lags:

$$C(\tau) = \frac{\langle (x(t) - \mu)(x(t+\tau) - \mu) \rangle}{\sigma^2}$$

- $C(0) = 1$ (perfectly correlated with itself)
- $C(\tau) \to 0$ as $\tau \to \infty$ (decorrelation at long times)

```python
from polyzymd.analyses.shared.autocorrelation import compute_acf

# Compute ACF from a timeseries (e.g., RMSD, distance)
acf_result = compute_acf(
    timeseries,
    timestep=10.0,      # ps between frames
    timestep_unit="ps",
)

# acf_result.acf contains the autocorrelation values
# acf_result.lags contains the time lags
```

### Correlation Time (τ)

The **correlation time** τ is the characteristic timescale for decorrelation.
Frames separated by more than 2τ are approximately independent.

PolyzyMD estimates τ using numerical integration of the ACF (the most robust
method for noisy data):

$$\tau = \int_0^{\infty} C(t) \, dt$$

```python
from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time

tau_result = estimate_correlation_time(
    acf_result,
    method="integration",  # Most robust for MD data
)

print(f"Correlation time: {tau_result.tau:.1f} {tau_result.tau_unit}")
print(f"Independent samples: {tau_result.n_independent}")
```

### Statistical Inefficiency (g)

The **statistical inefficiency** quantifies how much correlation inflates
variance estimates:

$$g = 1 + 2\sum_{t=1}^{N} C(t) \left(1 - \frac{t}{N}\right)$$

The effective number of independent samples is:

$$N_{\text{eff}} = \frac{N}{g}$$

```{tip}
A typical 100 ns trajectory with 10,000 frames might have only 50-200
effective independent samples, depending on the observable and system dynamics.
```

---

## The Critical Distinction: Means vs. Variances

**This is the most important concept in MD statistics.**

Different quantities require different statistical treatments depending on
whether they involve first moments (means) or second moments (variances).

### First Moments: Use All Data, Correct Uncertainty

For quantities that are **averages** (first moments):
- Mean distance
- Mean RMSD
- Contact fraction
- Average energy

**Correlation affects precision but NOT accuracy.**

The sample mean $\bar{x}$ is an unbiased estimator regardless of correlation.
Correlated samples still provide valid information about the mean—they just
don't provide as much *independent* information.

**Correct approach**: Use all frames, but correct the standard error:

$$\text{SEM} = \frac{\sigma}{\sqrt{N_{\text{eff}}}} = \frac{\sigma}{\sqrt{N/g}}$$

This is exactly what PolyzyMD's distances plugin does internally in
`compute_replicate()`:

```python
# From analyses/distances/__init__.py (DistanceCalculator internals)
# SEM = std / sqrt(n_independent)
if n_independent_frames > 0:
    sem_distance = float(std_dist / np.sqrt(n_independent_frames))
```

### Second Moments: Subsample to Independent Frames

For quantities that measure **fluctuations** (second moments):
- RMSF (root mean square fluctuation)
- Heat capacity (energy variance)
- Compressibility (volume variance)
- Order parameters from variance

**Correlation introduces systematic BIAS, not just imprecision.**

Correlated samples systematically **underestimate** the true variance because
consecutive frames show the system in similar configurations.

#### Physical Intuition: The Pendulum Analogy

Imagine measuring the position variance of a swinging pendulum:

| Sampling Rate | What You See | Apparent Variance |
|--------------|--------------|-------------------|
| 1000 Hz | Smooth continuous motion | **Small** (consecutive samples similar) |
| 0.1 Hz | Independent snapshots at random phases | **Correct** (true swing amplitude) |

The fast sampling doesn't capture the full range of motion between samples—it
oversamples the slow transitions and underestimates fluctuations.

**Correct approach**: Use only independent frames (spaced by ≥2τ):

```python
# From analyses/rmsf plugin internals
frame_indices = get_independent_indices(
    n_frames=n_frames_total,
    correlation_time=correlation_time,
    timestep=timestep,
    start_frame=start_frame,
)
```

---

## How PolyzyMD Implements This

### Distance analysis

The distances plugin computes distances for **all frames** in
`compute_replicate()` and uses autocorrelation to correct the uncertainty:

```python
# Compute distances for every frame
for i, ts in enumerate(u.trajectory):
    dist = np.linalg.norm(pos2 - pos1)
    distances.append(dist)

# Compute autocorrelation
tau_result = estimate_correlation_time(distances_arr, ...)

# Correct SEM using effective sample size
sem_distance = std_dist / np.sqrt(n_independent_frames)
```

**Why this works**: The mean distance estimate benefits from all data points.
More frames = more precise mean estimate. The corrected SEM properly reflects
our actual uncertainty.

### RMSF analysis

The RMSF plugin **subsamples** to independent frames in `compute_replicate()`
before computing fluctuations:

```python
# Compute RMSD timeseries for ACF
rmsd_timeseries = self._compute_rmsd_timeseries(u, atoms, start_frame)

# Estimate correlation time
tau_result = estimate_correlation_time(acf_result, n_frames=n_frames_after_eq)

# Select only independent frames (spaced by 2τ)
frame_indices = get_independent_indices(
    n_frames=n_frames_total,
    correlation_time=correlation_time,
    timestep=timestep,
    start_frame=start_frame,
)

# Compute RMSF using ONLY independent frames
rmsf_values = self._compute_rmsf(u, atoms, frame_indices)
```

**Why this is necessary**: RMSF measures fluctuations (variance). Using
correlated frames would systematically underestimate protein flexibility.

---

## Summary Table

| Quantity | Statistical Type | Effect of Correlation | PolyzyMD Approach |
|----------|-----------------|----------------------|-------------------|
| Mean distance | 1st moment | Reduces precision (no bias) | All frames + corrected SEM |
| Contact fraction | 1st moment | Reduces precision (no bias) | All frames + corrected SEM |
| RMSF | 2nd moment | **Introduces bias** | Subsample to independent frames |
| Distance variance | 2nd moment | **Introduces bias** | Subsample to independent frames |

---

## Practical Recommendations

### 1. Check Your Correlation Time

Always examine the correlation time relative to your trajectory length:

```python
from polyzymd.analyses.shared.autocorrelation import (
    compute_acf,
    estimate_correlation_time,
    check_statistical_reliability,
)

acf = compute_acf(observable, timestep=10.0, timestep_unit="ps")
tau = estimate_correlation_time(acf, n_frames=len(observable))

print(f"Correlation time: {tau.tau:.1f} {tau.tau_unit}")
print(f"Independent samples: {tau.n_independent}")
print(f"Statistically reliable: {tau.is_reliable}")
```

```{warning}
If `n_independent < 10`, your results may not be statistically reliable.
Consider:
1. Running longer simulations
2. Using multiple independent replicates
3. Reporting results with appropriate caveats
```

### 2. Use Multiple Replicates

Independent replicates provide truly uncorrelated samples. PolyzyMD's
aggregation functions properly combine replicates:

```bash
# Aggregate across configured independent replicates
polyzymd compare run distances -f comparison.yaml --eq-time 100ns
```

Replicates are configured per condition in `comparison.yaml`.

The aggregated SEM is computed from the **variance across replicate means**,
which is statistically robust regardless of within-trajectory correlation.

### 3. Report Uncertainties Correctly

Always report uncertainties that account for correlation:

```{admonition} Good Practice
:class: tip

"The mean Ser77-His133 distance was 3.42 ± 0.15 Å (SEM, N_eff = 47 independent
samples from 5 replicates of 100 ns each)."
```

```{admonition} Poor Practice
:class: warning

"The mean distance was 3.42 ± 0.02 Å" (using naive SEM with N = 10,000 frames)
```

### 4. Understand What You're Measuring

Before running analysis, ask yourself:

1. **Am I computing a mean or a variance?**
2. **What timescale does this process occur on?**
3. **Is my trajectory long enough to sample this process multiple times?**

For slow processes (e.g., loop conformational changes with τ ~ 10 ns), even
a 100 ns trajectory may only contain ~10 independent samples.

---

## Multiple Comparison Correction

### Why It Matters

When comparing *N* conditions pairwise, the number of hypothesis tests grows as
*N*(*N*−1)/2. With 5 conditions this is already 10 tests; with 10 conditions,
45 tests. Without correction, the probability of at least one false positive
(the **familywise error rate**) increases rapidly — even when every null
hypothesis is true.

### Benjamini-Hochberg Procedure

PolyzyMD uses the **Benjamini-Hochberg (BH)** step-up procedure to control the
**false discovery rate** (FDR) — the expected proportion of false positives
among rejected hypotheses. The algorithm:

1. Rank all *m* p-values in ascending order: $p_{(1)} \le p_{(2)} \le \dots \le p_{(m)}$.
2. For each rank *k*, compare $p_{(k)}$ to the threshold $(k / m) \times \alpha$.
3. Find the largest *k* where $p_{(k)} \le (k / m) \times \alpha$.
4. Reject all hypotheses with rank ≤ *k*.

Adjusted p-values are computed with step-down monotonicity enforcement so that
$p_{\text{adj},(k)} \le p_{\text{adj},(k+1)}$ always holds.

### How PolyzyMD Applies BH Correction

Each plugin defines its own hypothesis family:

- **Contacts** — All pairwise t-test p-values (both coverage and
  contact_fraction metrics) form **one hypothesis family**. ANOVA p-values
  form a separate small family (typically 2 tests).
- **Binding free energy** — Same-temperature pairwise entries form **one
  family per temperature group**.
- **Polymer affinity** — Same pattern as binding free energy: one family per
  temperature group.

### The `fdr_alpha` Parameter

`fdr_alpha` is configured per plugin in `comparison.yaml` under the `plugins:`
block. The default is 0.05. Lowering the value (e.g., 0.01) makes the
correction more conservative; raising it allows more discoveries at the cost of
a higher expected false positive rate.

```yaml
plugins:
  contacts:
    fdr_alpha: 0.01  # More conservative than default
```

### Raw vs Adjusted p-Values

Both raw and BH-adjusted p-values are shown in formatted output. The
`significant` field in saved JSON and the significance markers (`*`, `**`,
`***`) in CLI tables are based on the **adjusted** value. This ensures that
reported significance accounts for the number of comparisons performed.

### Relationship to Within-Trajectory Statistics

FDR correction and autocorrelation handling address **different** sources of
statistical error:

| Problem | Scope | PolyzyMD Solution |
|---------|-------|-------------------|
| **Autocorrelation** | Within a single trajectory — consecutive frames are not independent | Correct SEM via $N_{\text{eff}}$ or subsample by 2τ (see earlier sections) |
| **Multiple comparisons** | Across conditions — many hypothesis tests inflate false positives | BH correction on pairwise p-values |

Both corrections are important. Autocorrelation handling ensures that the
per-condition uncertainty is not artificially small; FDR correction ensures
that the cross-condition significance claims are not inflated by the number of
tests performed.

---

## LiveCoMS References

This implementation follows the guidelines from:

> **Grossfield, A., et al.** (2018). "Best Practices for Quantification of
> Uncertainty and Sampling Quality in Molecular Simulations."
> *Living Journal of Computational Molecular Science*, 1(1), 5067.
> [DOI: 10.33011/livecoms.1.1.5067](https://doi.org/10.33011/livecoms.1.1.5067)

Key sections:
- **Section 3.1**: Uncertainty in averages (first moments)
- **Section 3.2**: Uncertainty in fluctuations (second moments)
- **Section 4**: Practical recommendations for MD analysis

Additional references:

> **Chodera, J. D., et al.** (2007). "Use of the Weighted Histogram Analysis
> Method for the Analysis of Simulated and Parallel Tempering Simulations."
> *Journal of Chemical Theory and Computation*, 3(1), 26-41.

> **Flyvbjerg, H., & Petersen, H. G.** (1989). "Error estimates on averages
> of correlated data." *Journal of Chemical Physics*, 91(1), 461-466.

---

## API Reference

The statistical functions are available in `polyzymd.analyses.shared.autocorrelation`:

| Function | Description |
|----------|-------------|
| `compute_acf()` | Compute autocorrelation function via FFT |
| `estimate_correlation_time()` | Estimate τ from ACF or timeseries |
| `get_independent_indices()` | Get frame indices spaced by 2τ |
| `statistical_inefficiency()` | Compute g directly from timeseries |
| `statistical_inefficiency_multiple()` | Compute g from multiple timeseries |
| `n_effective()` | Compute N_eff = N / g |
| `check_statistical_reliability()` | Warn if N_eff < 10 |

See the {doc}`/api/core` documentation for full API details.
