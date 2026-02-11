# RMSF Analysis: Statistical Best Practices

A comprehensive guide to statistically rigorous RMSF analysis, following
LiveCoMS recommendations for uncertainty quantification in molecular dynamics
simulations.

```{note}
**Just need quick results?** See the [Quick Start Guide](analysis_rmsf_quickstart.md)
for copy-paste commands and minimal setup.
```

## Introduction

Molecular dynamics trajectories are **temporally correlated**. This means
consecutive frames are not independent samples, and naive statistical analysis
will dramatically underestimate uncertainties. This guide explains:

- Why correlation matters for RMSF analysis
- How PolyzyMD quantifies correlation (autocorrelation analysis)
- What the warning messages mean and what to do about them
- How to properly compare conditions with statistical rigor

The recommendations in this guide follow the LiveCoMS best practices article
by Grossfield et al. (2018), which provides the definitive treatment of
uncertainty quantification in MD simulations.

## What is RMSF?

**Root Mean Square Fluctuation (RMSF)** measures how much each atom or residue
fluctuates around its average position during a simulation:

$$
\text{RMSF}_i = \sqrt{\frac{1}{T} \sum_{t=1}^{T} \left( \mathbf{r}_i(t) - \langle \mathbf{r}_i \rangle \right)^2}
$$

Where:
- $\mathbf{r}_i(t)$ is the position of atom/residue $i$ at time $t$
- $\langle \mathbf{r}_i \rangle$ is the time-averaged position
- $T$ is the number of frames

### What RMSF Measures

| High RMSF | Low RMSF |
|-----------|----------|
| Flexible regions | Rigid regions |
| Loops, termini | Core β-sheets |
| Solvent-exposed | Buried residues |
| Functionally mobile | Structurally constrained |

### Relationship to B-factors

RMSF is related to crystallographic B-factors (temperature factors):

$$
B_i = \frac{8\pi^2}{3} \langle \Delta r_i^2 \rangle = \frac{8\pi^2}{3} \text{RMSF}_i^2
$$

This allows comparison between simulation and experimental data, though
crystal packing effects mean the correspondence is imperfect.

## The Correlation Problem

### Why Consecutive Frames Aren't Independent

Consider a 200 ns simulation saved every 100 ps, giving 2000 frames. You might
expect 2000 independent measurements of protein flexibility. **This is wrong.**

The protein doesn't "forget" its conformation between frames. If frame 500
shows a loop in position A, frame 501 will almost certainly show the same loop
in nearly the same position. The frames are **correlated**.

### Visualizing Correlation

Imagine measuring RMSD over time. The autocorrelation function (ACF) shows
how correlated the RMSD is with itself at different time lags:

```
ACF
1.0 |*
    | *
    |  *
    |   **
    |     ***
    |        ****
    |            *****
0.0 |________________***********
    0         τ              Time lag
```

The **correlation time (τ)** is the characteristic decay time. Frames
separated by less than τ are substantially correlated; frames separated by
more than 2τ are approximately independent.

### Consequences of Ignoring Correlation

If you ignore correlation and treat all N frames as independent:

| Quantity | Naive Estimate | True Value |
|----------|----------------|------------|
| Sample size | N = 2000 | N_ind ≈ 5-10 |
| Standard error | σ/√2000 ≈ 0.02σ | σ/√5 ≈ 0.45σ |
| Uncertainty | Underestimated by ~20× | Correct |

**This leads to false confidence in small differences that may not be real.**

## Autocorrelation Analysis in PolyzyMD

PolyzyMD automatically performs autocorrelation analysis to quantify the
effective number of independent samples.

### The Calculation

1. **Compute RMSD timeseries** after alignment to reference structure
2. **Calculate autocorrelation function (ACF)** using FFT
3. **Integrate ACF** to get correlation time τ
4. **Compute statistical inefficiency**: $g = 1 + 2\tau/\Delta t$
5. **Estimate independent samples**: $N_{\text{ind}} = N / g$

### The LiveCoMS Formula

The number of independent samples is computed as:

$$
N_{\text{ind}} = \frac{N}{1 + 2\sum_{j=1}^{N_{\max}} C_j}
$$

Where $C_j$ is the normalized autocorrelation at lag $j$. This is equivalent
to:

$$
N_{\text{ind}} = \frac{N}{g} \quad \text{where} \quad g = 1 + 2\tau/\Delta t
$$

And τ is obtained by integrating the ACF:

$$
\tau = \int_0^{t_{\max}} C(t) \, dt \approx \sum_j C_j \cdot \Delta t
$$

### Integration Cutoff

The ACF integration stops at the first zero crossing (or when ACF < 0.05).
This prevents integrating noise at long lag times, which would artificially
inflate τ.

### Example Output

When you run RMSF analysis, PolyzyMD reports:

```
Correlation time: 15394 ps (15.4 ns)
Statistical inefficiency: 308.9
Independent samples: 6 (from 2000 frames)
```

This means:
- The RMSD decorrelates over ~15 ns timescales
- Each independent sample represents ~300 correlated frames
- You effectively have only 6 independent measurements

## Understanding the Warnings

### The "Low Statistical Reliability" Warning

When N_independent < 10, PolyzyMD issues this warning:

```
WARNING: Low statistical reliability: only 6 independent samples
(recommended >= 10). Correlation time τ = 15394 ps is comparable to
or longer than the trajectory sampling window. Consider:
(1) extending simulation time,
(2) using multiple independent trajectories, or
(3) interpreting results with caution.
See Grossfield et al. (2018) LiveCoMS 1:5067.
```

### What This Warning Means

The warning indicates that your trajectory may not have sampled enough
independent conformations to precisely estimate equilibrium properties. It
does **not** mean:

- Your simulation is broken
- Your data is useless
- You must start over

### Decision Tree: What To Do

```
Got "Low statistical reliability" warning?
│
├─ Are you using multiple independent replicates?
│  ├─ YES → You're probably fine. Replicate-based statistics are robust.
│  └─ NO  → Consider running 3-5 replicates instead of extending.
│
├─ What kind of conclusion do you need?
│  ├─ Qualitative (stable vs unstable) → Warning is informational only.
│  └─ Quantitative (precise mean ± SEM) → Need more sampling.
│
└─ Are you comparing conditions?
   ├─ YES → Focus on replicate-to-replicate variation for significance.
   └─ NO  → Report uncertainty based on replicate spread if available.
```

### When the Warning is Not a Problem

1. **You have multiple replicates**: Inter-replicate variation provides a
   robust uncertainty estimate independent of within-trajectory correlation.

2. **You're making qualitative conclusions**: "The protein is stable" doesn't
   require precise RMSF values.

3. **The effect size is large**: A 50% difference in RMSF is meaningful even
   with uncertain absolute values.

## Replicates vs Longer Simulations

### The LiveCoMS Recommendation

> "Multiple independent simulations are preferable to a single long simulation"
> — Grossfield et al. (2018)

### Why Replicates Are Better

| Multiple Replicates | Single Long Simulation |
|--------------------|------------------------|
| Truly independent samples | Samples remain correlated |
| Tests reproducibility | May be trapped in metastable state |
| Detects rare events/outliers | Rare events may dominate or be missed |
| Parallelizable | Sequential |
| Robust to system issues | Single point of failure |

### The Independence Argument

Even with τ = 20 ns, a 200 ns simulation yields only ~10 independent samples.
Running that same simulation 3 times gives you 3 **truly independent**
measurements of the mean RMSF:

```
Replicate 1: mean RMSF = 0.755 Å
Replicate 2: mean RMSF = 0.693 Å
Replicate 3: mean RMSF = 0.696 Å

Replicate mean: 0.715 Å
Replicate std:  0.035 Å
Replicate SEM:  0.020 Å (= std/√3)
```

This replicate-based uncertainty is **valid regardless of within-trajectory
correlation** because the replicates started from different initial velocities.

### How Many Replicates?

| Replicates | Power to Detect | Practical Guidance |
|------------|-----------------|-------------------|
| 3 | Large effects (d > 2) | Minimum for publication |
| 5 | Medium effects (d > 1.3) | Recommended standard |
| 10 | Small effects (d > 0.9) | High-precision studies |

For most enzyme-polymer studies, **3-5 replicates per condition** is a
reasonable balance between statistical power and computational cost.

## Handling Incomplete Data

During active research, simulations are often in progress. You may want to analyze
available data without waiting for all replicates to complete. PolyzyMD handles
this gracefully by skipping missing or failed replicates with informative warnings.

### Missing Replicates

When a requested replicate's data cannot be found, PolyzyMD logs a warning and
continues with the remaining replicates:

```
WARNING: Skipping replicate 2: trajectory data not found.
Working directory not found: /scratch/user/project/run_2
Has replicate 2 been simulated?
```

### Failed Replicates

If a replicate exists but analysis fails (e.g., corrupted trajectory, missing
atoms), PolyzyMD logs the error and continues:

```
WARNING: Skipping replicate 3: analysis failed with error: 
Could not read DCD file: unexpected end of file
```

### Aggregation Summary

When aggregating with incomplete data, PolyzyMD reports which replicates
were successfully analyzed:

```
WARNING: Aggregating 2 of 3 requested replicates. Skipped: [2]
```

### Minimum Requirements

**At least 2 successful replicates are required for aggregation.** This is
because SEM calculation requires n ≥ 2. If fewer than 2 replicates succeed,
PolyzyMD raises an error:

```
ValueError: Aggregation requires at least 2 successful replicates, but only
1 succeeded. Failed replicates: [2, 3]
```

### Output File Naming

Output filenames reflect the actual replicates used, not the requested range:

| Requested | Successful | Filename |
|-----------|------------|----------|
| 1-3 | 1, 2, 3 | `rmsf_reps1-3_eq100ns.json` |
| 1-3 | 1, 3 | `rmsf_reps1_3_eq100ns.json` |
| 1-5 | 1, 2, 4 | `rmsf_reps1_2_4_eq100ns.json` |

Note: Contiguous ranges use a hyphen (`1-3`), non-contiguous use underscores (`1_3`).

### Best Practices for Incomplete Data

1. **Investigate missing data**: If replicates consistently fail, check:
   - Did the simulation complete? Check SLURM logs for timeouts or errors.
   - Are trajectory files in the expected location? Verify paths in config.yaml.
   - Is the trajectory corrupted? Try loading it manually with MDAnalysis.

2. **Document which replicates were used**: When publishing results from
   incomplete data, clearly state which replicates contributed to aggregated
   statistics. The JSON output includes a `replicates` field for this purpose.

3. **Re-run with complete data**: Once all simulations finish, re-run
   analysis with `--recompute` to include all replicates:
   ```bash
   polyzymd analyze rmsf -c config.yaml -r 1-5 --recompute
   ```

4. **Consider statistical implications**: Results from 2 replicates have
   larger uncertainty than 3+. Be appropriately cautious when interpreting
   results with fewer replicates than planned.

### Example: Analyzing During Active Simulations

A common workflow when simulations are still running:

```bash
# Request all 5 planned replicates, but only 3 have completed
polyzymd analyze rmsf -c config.yaml -r 1-5 --eq-time 100ns

# Output shows:
# Skipping replicate 4: trajectory data not found...
# Skipping replicate 5: trajectory data not found...
# Aggregating 3 of 5 requested replicates. Skipped: [4, 5]
#
# Results saved to: rmsf_reps1-3_eq100ns.json
```

Later, when all simulations complete:

```bash
# Re-run to include all replicates
polyzymd analyze rmsf -c config.yaml -r 1-5 --eq-time 100ns --recompute

# Now all 5 are included:
# Results saved to: rmsf_reps1-5_eq100ns.json
```

## Comparing Conditions

### Experimental Design

To compare two conditions (e.g., polymer vs no-polymer):

1. Run **N replicates** of each condition (N ≥ 3)
2. Compute RMSF for each replicate
3. Use **replicate means** as your data points (not per-frame values)
4. Perform statistical test on replicate means

### The Statistical Test

Use a two-sample t-test on the per-replicate mean RMSF values:

```python
from scipy import stats
import numpy as np

# Per-replicate mean RMSF (not per-frame!)
condition_A = [0.755, 0.693, 0.696]  # No polymer
condition_B = [0.558, 0.738, 0.496]  # With polymer

# Two-sample t-test
t_stat, p_value = stats.ttest_ind(condition_A, condition_B)
print(f"t-statistic: {t_stat:.3f}")
print(f"p-value: {p_value:.4f}")
```

### Effect Size

The p-value tells you if the difference is **statistically significant**, but
not if it's **meaningful**. Effect size quantifies the magnitude:

```python
# Cohen's d
mean_A = np.mean(condition_A)
mean_B = np.mean(condition_B)
pooled_std = np.sqrt((np.var(condition_A, ddof=1) + np.var(condition_B, ddof=1)) / 2)
cohens_d = (mean_A - mean_B) / pooled_std
print(f"Cohen's d: {cohens_d:.2f}")
```

| Cohen's d | Interpretation |
|-----------|---------------|
| < 0.2 | Negligible effect |
| 0.2 - 0.5 | Small effect |
| 0.5 - 0.8 | Medium effect |
| > 0.8 | Large effect |

### Power Analysis

With small sample sizes, you can only detect large effects. The **minimum
detectable difference** depends on:

- Sample size per group (N)
- Within-group variability (σ)
- Desired significance level (α, typically 0.05)
- Desired power (typically 0.80)

```python
# Minimum detectable difference (approximate)
# For t-test with α=0.05, power=0.80
critical_t = {3: 2.78, 5: 2.31, 10: 2.10}  # df = 2*(N-1)
se_diff = pooled_std * np.sqrt(2/N)
min_detectable = critical_t[N] * se_diff
print(f"Minimum detectable difference: {min_detectable:.3f} Å")
```

### Worked Example

From a real analysis comparing polymer vs no-polymer:

```
RESULTS SUMMARY
===============
                    No Polymer    With Polymer
Replicate 1         0.755 Å       0.558 Å
Replicate 2         0.693 Å       0.738 Å
Replicate 3         0.696 Å       0.496 Å
-----------------------------------------------
Mean                0.715 Å       0.597 Å
Std                 0.035 Å       0.126 Å

STATISTICAL ANALYSIS
====================
Difference: 0.117 Å (16.4% reduction with polymer)
t-statistic: 1.56
p-value: 0.194
Cohen's d: 1.27 (large effect)

INTERPRETATION
==============
- Large effect size (d = 1.27) suggests a real difference
- p > 0.05 means we cannot reject null hypothesis
- High variability in polymer condition (replicate 2 is outlier)
- CONCLUSION: Suggestive but not conclusive; need more replicates
```

### Non-Significant p-value with Large Effect Size

This situation (p > 0.05 but large Cohen's d) is common with small samples.
It means:

1. The effect **might be real** but you lack power to detect it
2. More replicates would likely achieve significance
3. You should **not** conclude "no effect" — only "insufficient evidence"

The appropriate response is to run additional replicates, not to abandon
the hypothesis.

## Interpreting RMSF Values

### Typical Ranges

| RMSF (Å) | Protein Region | Interpretation |
|----------|---------------|----------------|
| 0.3 - 0.5 | Core β-sheets, buried helices | Very rigid |
| 0.5 - 1.0 | Surface helices, structured loops | Normal flexibility |
| 1.0 - 2.0 | Exposed loops, active site flaps | Flexible |
| 2.0 - 5.0 | Termini, disordered regions | Highly flexible |
| > 5.0 | Disordered tails, unfolded regions | May indicate instability |

### Factors Affecting RMSF

**Temperature**: Higher T → higher RMSF (more kinetic energy)

**Solvent**: 
- Water → normal flexibility
- Organic co-solvents → may increase or decrease depending on interactions
- Crowding agents → typically reduce RMSF

**Simulation parameters**:
- Force field → affects flexibility
- Timestep → too large can cause instability
- Thermostat → Langevin adds friction, affects dynamics

**System composition**:
- Polymers can restrict or enhance motion
- Ligand binding often rigidifies active site
- Oligomerization affects interface flexibility

### Active Site Flexibility

For enzyme studies, active site RMSF is particularly relevant:

- **Lower RMSF** may indicate:
  - Better substrate positioning
  - More stable catalytic geometry
  - Potentially higher activity (if not too rigid)

- **Higher RMSF** may indicate:
  - Conformational sampling for substrate binding
  - Induced fit dynamics
  - Potentially reduced activity if catalytic residues misalign

The "optimal" flexibility depends on the enzyme mechanism and should be
interpreted in context of experimental activity data.

## Common Pitfalls

### Pitfall 1: Using Unaligned Trajectories

**Symptom**: RMSF values > 10 Å, sometimes 50+ Å

**Cause**: Without alignment, RMSF includes global translation and rotation
of the entire protein through the simulation box.

**Solution**: PolyzyMD automatically aligns trajectories. If you see very
high values, check that:
- Your selection string matches atoms in the system
- The reference structure is valid
- Trajectory files are complete

### Pitfall 2: Insufficient Equilibration

**Symptom**: RMSF values drift over time; different equilibration times
give different results

**Cause**: Including non-equilibrated frames in the analysis

**Solution**: Use `--eq-time` to skip the equilibration period:

```bash
# Skip first 10% of a 200 ns simulation
polyzymd analyze rmsf -c config.yaml -r 1 --eq-time 20ns
```

Monitor RMSD vs time — equilibration is complete when RMSD plateaus.

### Pitfall 3: Over-Interpreting Small Differences

**Symptom**: Claiming significance for 0.05 Å differences

**Cause**: Not accounting for uncertainty

**Solution**: Always report uncertainty and perform statistical tests:

```python
# WRONG: "Condition A (0.715 Å) is more stable than B (0.720 Å)"
# RIGHT: "Condition A (0.715 ± 0.020 Å) and B (0.720 ± 0.025 Å) 
#         are not significantly different (p = 0.89)"
```

### Pitfall 4: Ignoring Replicate Variation

**Symptom**: Reporting within-trajectory SEM as the uncertainty

**Cause**: Treating correlated frames as independent

**Solution**: Use replicate-based statistics:

```python
# WRONG: SEM from 2000 correlated frames
sem_wrong = np.std(all_frames) / np.sqrt(2000)  # Way too small

# RIGHT: SEM from 3 independent replicates
replicate_means = [0.755, 0.693, 0.696]
sem_right = np.std(replicate_means, ddof=1) / np.sqrt(3)  # Correct
```

### Pitfall 5: Cherry-Picking Replicates

**Symptom**: Excluding "outlier" replicates without justification

**Cause**: Desire for cleaner results

**Solution**: Include all replicates unless there's a **technical** reason
to exclude (e.g., simulation crashed, obvious system error). Biological
variability is real and should be reported.

If one replicate differs substantially, investigate **why**:
- Did a conformational change occur?
- Did the ligand unbind?
- Is there a slow process being sampled?

## Advanced: Python API

### Accessing Autocorrelation Results

```python
from polyzymd.analysis import RMSFCalculator
import json

# Run analysis
calc = RMSFCalculator(config_path="config.yaml", equilibration="100ns")
result = calc.compute(replicate=1)

# Access autocorrelation data
print(f"Correlation time: {result.correlation_time:.1f} ps")
print(f"N independent: {result.n_independent_frames}")

# Load full result for detailed inspection
with open(result.source_file) as f:
    data = json.load(f)

# Per-residue data
for resid, rmsf in zip(data['residue_ids'], data['rmsf_per_residue']):
    print(f"Residue {resid}: {rmsf:.3f} Å")
```

### Custom Statistical Analysis

```python
import numpy as np
from scipy import stats

# Load multiple conditions
conditions = {
    'no_polymer': [0.755, 0.693, 0.696],
    'polymer_A': [0.558, 0.738, 0.496],
    'polymer_B': [0.612, 0.598, 0.621],
}

# Pairwise comparisons with Bonferroni correction
n_comparisons = 3
alpha = 0.05 / n_comparisons

for name_a, values_a in conditions.items():
    for name_b, values_b in conditions.items():
        if name_a >= name_b:
            continue
        t, p = stats.ttest_ind(values_a, values_b)
        sig = "**" if p < alpha else ""
        print(f"{name_a} vs {name_b}: p = {p:.4f} {sig}")
```

### Exporting for External Analysis

```python
import pandas as pd
import json

# Load aggregated result
with open("analysis/rmsf/aggregated/rmsf_reps1-3_eq10ns.json") as f:
    data = json.load(f)

# Create DataFrame for export
df = pd.DataFrame({
    'residue_id': data['residue_ids'],
    'residue_name': data['residue_names'],
    'mean_rmsf': data['mean_rmsf_per_residue'],
    'sem_rmsf': data['sem_rmsf_per_residue'],
})

# Export to CSV
df.to_csv("rmsf_results.csv", index=False)

# Export to Excel with multiple sheets
with pd.ExcelWriter("rmsf_analysis.xlsx") as writer:
    df.to_excel(writer, sheet_name="Per-Residue", index=False)
    
    summary = pd.DataFrame({
        'Metric': ['Mean RMSF', 'SEM', 'Min', 'Max', 'N_replicates'],
        'Value': [
            data['overall_mean_rmsf'],
            data['overall_sem_rmsf'],
            data['overall_min_rmsf'],
            data['overall_max_rmsf'],
            data['n_replicates'],
        ]
    })
    summary.to_excel(writer, sheet_name="Summary", index=False)
```

## References

### Primary Reference

**Grossfield A, Patrone PN, Roe DR, Schultz AJ, Siderius DW, Zuckerman DM.**
(2018) "Best Practices for Quantification of Uncertainty and Sampling Quality
in Molecular Simulations." *Living Journal of Computational Molecular Science*
1(1):5067. https://doi.org/10.33011/livecoms.1.1.5067

This is the definitive guide to uncertainty quantification in MD. The
methodology in PolyzyMD follows their recommendations.

### Additional References

**Flyvbjerg H, Petersen HG.** (1989) "Error estimates on averages of correlated
data." *Journal of Chemical Physics* 91:461-466.
https://doi.org/10.1063/1.457480

Classic paper on block averaging for correlated data.

**Chodera JD, Swope WC, Pitera JW, Seok C, Dill KA.** (2007) "Use of the
Weighted Histogram Analysis Method for the Analysis of Simulated and Parallel
Tempering Simulations." *Journal of Chemical Theory and Computation* 3:26-41.
https://doi.org/10.1021/ct0502864

Introduces statistical inefficiency for MD analysis.

**Knapp B, Frantal S, Greshake B, Schwarz R, et al.** (2018) "Is an Intuitive
Convergence Definition of Molecular Dynamics Simulations Solely Based on the
Root Mean Square Deviation Possible?" *Journal of Computational Biology*
25:1069-1077.

Discussion of RMSD-based convergence assessment.

## See Also

- [Quick Start Guide](analysis_rmsf_quickstart.md) — Get results fast
- [Reference Structure Selection](analysis_reference_selection.md) — Choose alignment reference
- [LiveCoMS Best Practices](https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5067) — Full methodology paper
