# RMSD Analysis: Best Practices

A guide to interpreting RMSD results, choosing reference structures,
understanding autocorrelation, and comparing conditions rigorously.

```{versionadded} 1.3.0
The RMSD analysis plugin was added in PolyzyMD 1.3.0.
```

```{note}
**Just need quick results?** See the [Quick Start Guide](../how_to/analysis_rmsd_quickstart.md)
for copy-paste commands and minimal setup.
```

```{seealso}
**For foundational statistical concepts** (autocorrelation, correlation time,
the difference between means vs. variances), see the
[Statistics Best Practices Guide](analysis_statistics_best_practices.md).

This page focuses on **RMSD-specific** guidance: what the values mean,
how to interpret timeseries behavior, and how to compare conditions.
```

## What is RMSD?

**Root Mean Square Deviation (RMSD)** measures the average distance between
atoms in a structure and a reference structure after optimal superposition:

$$
\text{RMSD}(t) = \sqrt{\frac{1}{N} \sum_{i=1}^{N} \left\| \mathbf{r}_i(t) - \mathbf{r}_i^{\text{ref}} \right\|^2}
$$

Where:
- $\mathbf{r}_i(t)$ is the position of atom $i$ at time $t$
- $\mathbf{r}_i^{\text{ref}}$ is the position of atom $i$ in the reference structure
- $N$ is the number of atoms in the selection

Unlike RMSF (which averages over time to give one value per residue), RMSD
gives **one value per frame** — a timeseries that tracks how the structure
evolves relative to its reference.

### What RMSD Measures

| RMSD Behavior | Structural Interpretation |
|---------------|--------------------------|
| Low and stable | Structure stays near reference; rigid or well-equilibrated |
| Gradually increasing | Conformational drift — may indicate slow unfolding |
| Plateau after rise | Equilibration followed by stable sampling |
| Sudden jump | Conformational transition or partial unfolding event |
| Oscillating | Sampling between multiple conformational states |

## Interpreting RMSD Values

### Typical Ranges for Enzyme Systems

| RMSD (Å) | Interpretation | Typical Origin |
|----------|----------------|----------------|
| 0.5 – 1.5 | Very stable | Well-folded globular protein, strong crystal contacts |
| 1.5 – 2.5 | Normal | Typical for soluble enzymes at 300 K |
| 2.5 – 3.5 | Moderate deviation | Flexible loops, domain motions, lid opening |
| 3.5 – 5.0 | Large deviation | Significant conformational change, partial unfolding |
| > 5.0 | Very large | Major structural rearrangement or unfolding |

```{note}
These ranges assume Cα RMSD for proteins of 100–400 residues. Larger
proteins tend to have higher RMSD due to more flexible regions. Always
compare like-with-like: same selection, same reference mode, same atoms.
```

### Selection Matters

The choice of atoms for RMSD calculation strongly affects the result:

| Selection | Typical RMSD | Best For |
|-----------|-------------|----------|
| `protein and name CA` | 1–3 Å | Overall backbone stability |
| `protein and backbone` | 1–3 Å | Backbone conformation (includes N, C, O) |
| `protein and name CA and resid 50:150` | 0.5–2 Å | Core region, excluding flexible termini |
| Active site residues | 0.5–2 Å | Catalytic geometry maintenance |
| `chainid C and not name H*` | Variable | Polymer conformation tracking |

## RMSD vs Time: What to Look For

### Plateau Behavior (Ideal)

```text
RMSD
 3 |          ___________
   |         /
 2 |        /
   |       /
 1 |      /
   |_____/
 0 +----------------------→ Time
   0    10   20   30   40 ns
```

The initial rise is **equilibration** — the system relaxing from the starting
structure. The plateau indicates **stable sampling** of an equilibrium
ensemble. Set `--eq-time` to skip the rise phase.

### Conformational Drift

```text
RMSD
 5 |                    /
   |                   /
 4 |                  /
   |                 /
 3 |                /
   |_______________/
 0 +----------------------→ Time
```

A continuously rising RMSD suggests the system has not equilibrated or is
slowly unfolding. Possible responses:
- Extend the simulation to see if it eventually plateaus
- Check temperature, pressure, and force field parameters
- The system may genuinely be unstable under these conditions

### Sudden Jumps

A sharp RMSD increase mid-trajectory typically indicates a conformational
transition. Check the structure at the jump to understand what happened:
- Loop flipping or lid opening
- Domain rearrangement
- Ligand unbinding
- Partial unfolding

```{tip}
When you observe a jump, load the trajectory in a molecular viewer (e.g.,
VMD, PyMOL) and examine frames around the transition. This often reveals
biologically meaningful events.
```

### Oscillations

Regular RMSD oscillations suggest the system samples between two or more
metastable states. This is often seen with:
- Hinge-bending motions
- Active site lid dynamics
- Allosteric transitions

For oscillating systems, report the **range and period** of oscillation rather
than just the mean RMSD.

## How PolyzyMD Handles Autocorrelation

RMSD timeseries are correlated — adjacent frames are similar because MD
evolves continuously. PolyzyMD automatically accounts for this:

1. **Computes RMSD timeseries** after alignment to the reference structure
2. **Estimates correlation time (τ)** via autocorrelation function integration
3. **Computes effective sample size** — `n_independent = n_frames / (2τ)`
4. **Reports autocorrelation-corrected SEM** — `SEM = σ / √n_independent`

### Example Autocorrelation Output

```text
Run: Protein Backbone
  Correlation time: 4521 ps (4.5 ns)
  Statistical inefficiency: 562.7
  Independent samples: 16 (from 9000 frames)
  SEM (corrected): 0.078 Å
```

This means:
- RMSD values decorrelate over ~4.5 ns timescales
- You effectively have 16 independent measurements from 9000 frames
- The reported SEM properly accounts for this correlation

```{seealso}
For the mathematical details of autocorrelation functions and the LiveCoMS
recommendations, see the
[Statistics Best Practices Guide](analysis_statistics_best_practices.md).
```

## Multi-Run Analysis: Why and When

### Why Multiple Runs?

Different RMSD selections answer different questions:

| Run Label | Selection | Question |
|-----------|-----------|----------|
| "Protein Backbone" | `protein and name CA` | Overall protein stability? |
| "Active Site" | Catalytic residues CA | Catalytic geometry maintenance? |
| "Polymer Core" | `chainid C and not name H*` | Polymer conformation stability? |
| "Crystal Deviation" | `protein and name CA` (external ref) | Distance from functional state? |

### When to Use Multi-Run

- **Always** include a global protein backbone run as a baseline
- **Add active site runs** when studying enzyme catalysis
- **Add polymer runs** when studying enzyme-polymer conjugates
- **Add external reference runs** when comparing to a known functional state

### Independent Ranking

Each run is ranked independently across conditions. This prevents averaging
RMSD from structurally different selections (which would be meaningless):

```text
Rankings:
  Protein Backbone: With Polymer < No Polymer (stabilizing)
  Active Site:      With Polymer < No Polymer (stabilizing)
  Polymer Core:     No Polymer — (single condition only)
```

## External Reference for Catalytic Competence

When studying enzyme catalysis across multiple conditions (polymer baths,
temperatures, mutants), the standard reference modes (`centroid`, `average`)
use a **condition-specific** reference — each condition's trajectory determines
its own reference structure.

The `external` reference mode uses a **condition-independent** reference,
typically a crystal structure representing the catalytically competent
geometry. RMSD then measures deviation from the functional state:

$$
\text{RMSD}^{\text{ext}}(t) = \sqrt{\frac{1}{N} \sum_{i=1}^{N} \left\| \mathbf{r}_i(t) - \mathbf{r}_i^{\text{crystal}} \right\|^2}
$$

**Interpretation changes with external reference:**

| Metric | Standard RMSD (centroid/average) | External Reference RMSD |
|--------|----------------------------------|------------------------|
| Low value | Structure stays near its own equilibrium | Structure stays near crystal geometry |
| High value | Structure deviates from its equilibrium | Structure deviates from functional state |
| Condition comparison | Which condition is more structurally stable? | Which condition best maintains catalytic geometry? |

```{tip}
**Which reference mode for enzymes?** Use `centroid` or `average` to answer
"how stable is the protein?" Use `external` to answer "how well does the
protein maintain its catalytically competent geometry?" These are complementary
questions — consider running both as separate runs in the same analysis.
```

## Replicates vs Longer Simulations

### The LiveCoMS Recommendation

> "Multiple independent simulations are preferable to a single long simulation"
> — Grossfield et al. (2018)

### Why Replicates Matter for RMSD

| Multiple Replicates | Single Long Simulation |
|--------------------|------------------------|
| Truly independent starting points | Frames remain correlated |
| Tests reproducibility of drift/plateau | May be trapped in metastable state |
| Robust SEM from replicate means | SEM requires autocorrelation correction |
| Parallelizable | Sequential |

### How Many Replicates?

| Replicates | Statistical Power | Practical Guidance |
|------------|-----------------|-------------------|
| 1 | Descriptive only | Exploratory — no inferential statistics |
| 3 | Large effects (d > 2) | Minimum for publication |
| 5 | Medium effects (d > 1.3) | Recommended standard |

```{note}
With only 1 replicate, PolyzyMD still computes RMSD and includes the condition
in descriptive summaries and rankings. Replicate SEM is unavailable and shown
as n/a because variability across independent simulations cannot be estimated
from a singleton. Pairwise inferential tests require at least 2 replicates per
condition.
```

## Comparing Conditions

### What PolyzyMD Computes

For each RMSD run, the comparison produces:

| Statistic | Description |
|-----------|-------------|
| **Ranking** | Conditions sorted by mean RMSD (lowest = most stable) |
| **Percent change** | Relative to control condition |
| **Direction** | `stabilizing` (< −1%), `destabilizing` (> +1%), or `unchanged` |
| **t-statistic** | Two-sample t-test on replicate means |
| **p-value** | Two-tailed significance |
| **Cohen's d** | Effect size magnitude |
| **ANOVA** | Omnibus F-test when 3+ conditions (per-run) |

### Direction Labels

PolyzyMD classifies the direction of change based on percent change in mean
RMSD relative to control:

| Percent Change | Direction | Meaning |
|---------------|-----------|---------|
| < −1% | `stabilizing` | Treatment reduces structural deviation |
| > +1% | `destabilizing` | Treatment increases structural deviation |
| −1% to +1% | `unchanged` | No meaningful difference |

### Interpreting the Comparison

```text
RMSD Comparison — Protein Backbone
===================================
Ranking (lower = more stable):
  1. 100% SBMA:   1.612 ± 0.028 Å
  2. No Polymer:  1.856 ± 0.034 Å
  3. 100% EGMA:   2.103 ± 0.051 Å

100% SBMA vs No Polymer:
  Change: -13.1% (stabilizing), p=0.0089*, d=2.41 (large)

100% EGMA vs No Polymer:
  Change: +13.3% (destabilizing), p=0.0142*, d=1.98 (large)

ANOVA: F=18.42, p=0.0031* (significant across all conditions)
```

**Reading this output:**
- SBMA polymer significantly stabilizes the enzyme (lower RMSD)
- EGMA polymer significantly destabilizes it (higher RMSD)
- The ANOVA confirms at least one condition differs from the others
- Large Cohen's d values mean these are substantial effects

## Common Pitfalls

### 1. Insufficient Equilibration

**Symptom:** RMSD mean and comparison results change with different
`--eq-time` values.

**Cause:** Including the initial equilibration rise biases the mean upward.

**Solution:** Plot the RMSD timeseries and visually identify when the plateau
begins. Set `--eq-time` to skip the rise phase:

```bash
polyzymd compare run rmsd -f comparison.yaml --eq-time 20ns
```

### 2. Comparing Different Selections

**Symptom:** RMSD values are not comparable across runs or publications.

**Cause:** Different atom selections yield different RMSD magnitudes.

**Solution:** Always report the exact selection string. Compare only runs with
identical selections.

### 3. Over-Interpreting Small Differences

**Symptom:** Claiming significance for 0.05 Å differences.

**Cause:** Not accounting for uncertainty.

**Solution:** Always report uncertainty and check statistical significance:

```text
# WRONG: "Condition A (1.856 Å) is less stable than B (1.861 Å)"
# RIGHT: "Condition A (1.856 ± 0.034 Å) and B (1.861 ± 0.028 Å)
#         are not significantly different (p = 0.91, unchanged)"
```

### 4. Ignoring Timeseries Shape

**Symptom:** Reporting only mean RMSD without inspecting the timeseries.

**Cause:** Two conditions can have the same mean RMSD but very different
dynamics (one plateaued, one drifting).

**Solution:** Always examine the RMSD timeseries plots. Use
`polyzymd compare plot-all` to generate them automatically.

### 5. Using All-Atom RMSD Without Justification

**Symptom:** Very high RMSD values even for stable proteins.

**Cause:** Side-chain motions dominate all-atom RMSD, obscuring backbone
changes.

**Solution:** Use Cα or backbone atoms unless you have a specific reason to
include side chains:

```yaml
selection: "protein and name CA"         # Recommended default
# NOT: selection: "protein"              # All-atom, noisy
```

### 6. Ignoring Replicate Variation

**Symptom:** Reporting within-trajectory SEM as the total uncertainty.

**Cause:** Treating autocorrelation-corrected SEM as sufficient.

**Solution:** Use replicate-based statistics when available. The replicate
SEM captures system-level variability that within-trajectory analysis cannot:

```text
Replicate 1: mean RMSD = 1.823 Å
Replicate 2: mean RMSD = 1.891 Å
Replicate 3: mean RMSD = 1.854 Å

Replicate mean: 1.856 Å
Replicate SEM:  0.034 Å  ← This is the gold standard uncertainty
```

### 7. Choosing Wrong Reference Mode

**Symptom:** Unexpected or hard-to-interpret comparison results.

**Cause:** Using `centroid` when `external` is more appropriate (or vice versa).

**Solution:** Match reference mode to your scientific question:
- Stability question → `centroid` or `average`
- Functional geometry question → `external` with crystal structure

### 8. Treating Automated Convergence as Ground Truth

**Symptom:** Trusting the `converged: true` flag without further inspection.

**Cause:** The sliding-window heuristic is parameter-dependent and can
miss slow drift, metastable trapping, or convergence issues in observables
other than RMSD.

**Solution:** Use convergence detection as one input among several. Always
inspect the RMSD timeseries visually, run multiple independent replicates,
and check convergence of other relevant observables (Rg, SASA, active site
distances). See {doc}`/explanation/convergence_detection` for a full
discussion of limitations.

## RMSD as an Equilibration Diagnostic

RMSD is the standard diagnostic for determining when equilibration is
complete. The principle is simple: when RMSD stops increasing and reaches
a stable plateau, the system is equilibrated.

As of v1.3.0, PolyzyMD also provides **automated convergence detection** — a
sliding-window slope heuristic that quantifies whether the RMSD timeseries has
plateaued. This runs automatically alongside every RMSD computation and reports
convergence status in the result JSON. See
{doc}`/explanation/convergence_detection` for the full conceptual discussion.

### Practical Workflow

1. Run RMSD with `--eq-time 0ns` (no equilibration skip) to see the full
   timeseries
2. Identify the plateau onset visually from the timeseries plot
3. Re-run with `--eq-time <plateau_onset>` for production analysis

```bash
# Step 1: Full timeseries
polyzymd compare run rmsd -f comparison.yaml --eq-time 0ns
polyzymd compare plot-all -f comparison.yaml

# Step 2: Inspect plots, identify plateau at ~10 ns

# Step 3: Production analysis
polyzymd compare run rmsd -f comparison.yaml --eq-time 10ns --recompute
```

```{tip}
If RMSD never plateaus within your simulation time, the system may need longer
equilibration, or the conditions may genuinely prevent equilibration (e.g.,
protein unfolding). Consider extending the simulation or investigating the
cause.
```

### Automated Convergence Detection

```{versionadded} 1.3.0
```

In addition to the manual workflow above, PolyzyMD automatically runs a
sliding-window convergence diagnostic on every RMSD timeseries. The algorithm
divides the timeseries into overlapping windows, computes the slope between
successive window means, and declares convergence when the slope remains below
a threshold for a sustained duration.

The convergence result is stored in the per-replicate JSON (`converged`,
`convergence_time_ns`, `convergence_message`) and aggregated across replicates
(`n_converged_replicates`, `convergence_fraction`,
`mean_convergence_time_ns`). Optional convergence diagnostic plots — showing
the RMSD timeseries alongside the sliding slope trace — can be enabled with
`show_convergence_plots: true` in `plot_settings.rmsd`.

**This is a diagnostic tool, not a definitive convergence proof.** The
heuristic can miss slow drift below the slope threshold, and convergence in
RMSD does not guarantee convergence of other observables. Always use multiple
replicates and visual inspection alongside automated diagnostics.

For a full conceptual treatment — including the algorithm, parameters, tuning
guidance, and limitations — see {doc}`/explanation/convergence_detection`.

## References

### Primary Reference

**Grossfield A, Patrone PN, Roe DR, Schultz AJ, Siderius DW, Zuckerman DM.**
(2018) "Best Practices for Quantification of Uncertainty and Sampling Quality
in Molecular Simulations." *Living Journal of Computational Molecular Science*
1(1):5067. https://doi.org/10.33011/livecoms.1.1.5067

### Additional References

**Knapp B, Frantal S, Greshake B, Schwarz R, et al.** (2018) "Is an Intuitive
Convergence Definition of Molecular Dynamics Simulations Solely Based on the
Root Mean Square Deviation Possible?" *Journal of Computational Biology*
25:1069-1077.

Discussion of RMSD-based convergence assessment and its limitations.

**Maiorov VN, Crippen GM.** (1994) "Significance of Root-Mean-Square Deviation
in Comparing Three-dimensional Structures of Globular Proteins." *Journal of
Molecular Biology* 235(2):625-634.
https://doi.org/10.1006/jmbi.1994.1017

Foundational work on RMSD as a structural similarity measure.

**Sargsyan K, Grauffel C, Bhagdev C.** (2017) "How Molecular Size Impacts RMSD
Applications in Molecular Dynamics Simulations." *Journal of Chemical Theory
and Computation* 13(4):1518-1524.
https://doi.org/10.1021/acs.jctc.7b00028

Analysis of how protein size affects expected RMSD values.

## See Also

- [Quick Start Guide](../how_to/analysis_rmsd_quickstart.md) — Get results fast
- [Convergence Detection](convergence_detection.md) — Conceptual guide to convergence: algorithm, parameters, and limitations
- [Statistics Best Practices](analysis_statistics_best_practices.md) — Foundational statistics for MD
- [RMSF Best Practices](analysis_rmsf_best_practices.md) — Per-residue fluctuation analysis
- [Reference Structure Selection](analysis_reference_selection.md) — Choose alignment reference
- [Compare Simulation Conditions](../how_to/analysis_compare_conditions.md) — Full comparison workflow
- [LiveCoMS Best Practices](https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5067) — Full methodology paper
