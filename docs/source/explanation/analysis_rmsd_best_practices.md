# RMSD Interpretation: Use, Limits, and Cautions

RMSD is a useful structural similarity diagnostic, but it is not proof of
thermodynamic equilibration, statistical convergence, or biological stability by
itself. This page explains what RMSD can and cannot support when interpreting
PolyzyMD trajectories, with emphasis on reference choice, atom selection,
autocorrelation, and cautious condition-level comparison.

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

This page focuses on **RMSD-specific** interpretation: what the values mean,
which assumptions they depend on, and where the conclusions can be ambiguous.
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

Unlike RMSF, which averages over time to give one value per residue, RMSD gives
one value per frame. The result is a timeseries describing how a selected set of
atoms moves relative to a chosen reference.

### What RMSD can and cannot tell you

RMSD measures distance from a chosen reference for a chosen atom selection. It
does not directly measure free energy, functional activity, or thermodynamic
stability. The same RMSD value can arise from different molecular motions, and a
single RMSD timeseries can hide local rearrangements that matter biologically.

| RMSD behavior | Cautious structural interpretation |
|---------------|------------------------------------|
| Low and stable | Selected atoms remain close to the chosen reference over the sampled interval. |
| Gradually increasing | Possible drift away from the reference; the cause is not unique. |
| Plateau after rise | Suggests structural stationarity for the selected atoms/reference, not full thermodynamic equilibration. |
| Sudden jump | May reflect a transition, alignment artifact, domain motion, ligand event, or unfolding; inspect structures. |
| Oscillating | May indicate repeated conformational motion or reference/alignment sensitivity; not necessarily two-state behavior. |

## Interpreting RMSD Values

### Rough Cα RMSD heuristics for folded globular proteins

The following values are rough heuristics for **Cα RMSD of small-to-medium
folded globular proteins**. They are not universal quality thresholds and should
not be used to label a system as stable or unstable without additional context.

| Cα RMSD (Å) | Possible interpretation | Common contributors |
|-------------|------------------------|---------------------|
| 0.5 – 1.5 | Close to reference for the selected atoms | Rigid core, short trajectory, restrained or crystal-like geometry |
| 1.5 – 2.5 | Modest deviation from reference | Typical backbone fluctuations for many compact proteins |
| 2.5 – 3.5 | Larger deviation from reference | Flexible loops, termini, lid opening, domain motion |
| 3.5 – 5.0 | Large reference-relative change | Domain rearrangement, alignment sensitivity, partial unfolding |
| > 5.0 | Very large reference-relative change | Major rearrangement, unfolding, or different conformational basin |

```{note}
These heuristics assume comparable atom selections, alignment choices,
reference modes, protein sizes, simulation lengths, and force-field contexts.
Always compare like-with-like: same selection, same reference mode, same atoms.
```

### Selection matters

The choice of atoms for RMSD calculation strongly affects the result:

| Selection | Typical use | Interpretation caution |
|-----------|-------------|------------------------|
| `protein and name CA` | Global backbone similarity | Flexible termini and domain motions can dominate. |
| `protein and backbone` | Backbone conformation | Includes more atoms than Cα and can change the scale. |
| `protein and name CA and resid 50:150` | Core-region similarity | Excludes regions that may be scientifically important. |
| Active site residues | Local catalytic-geometry proxy | Low RMSD does not prove catalytic competence. |
| `chainid C and not name H*` | Polymer conformation relative to reference | Polymer RMSD can be highly reference- and alignment-dependent. |

## RMSD vs Time: Interpreting Patterns Cautiously

### Plateau-like behavior

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

An initial rise followed by an apparent plateau can suggest that the selected
atoms have reached a reference-relative stationary regime. This is useful, but
limited: it does not prove thermodynamic equilibration, convergence of other
observables, or adequate sampling of all relevant conformations.

### Conformational drift

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

A continuously rising RMSD suggests ongoing movement away from the reference for
the selected atoms. Possible explanations include slow relaxation, domain
motion, unfolding, reference mismatch, alignment choices, or insufficient
sampling. RMSD alone usually cannot identify which explanation is correct.

### Sudden jumps

A sharp RMSD increase mid-trajectory indicates a rapid change in
reference-relative geometry, but the molecular cause is non-unique. It may
reflect loop flipping, lid opening, domain rearrangement, ligand motion,
alignment sensitivity, imaging artifacts, or partial unfolding.

```{tip}
When you observe a jump, load the trajectory in a molecular viewer and examine
frames around the transition. Visual inspection can distinguish chemically
meaningful events from alignment, imaging, or selection artifacts.
```

### Oscillations

Regular RMSD oscillations can be consistent with repeated motion between
reference-relative geometries, but they do not by themselves establish discrete
metastable states. Hinge bending, active-site lid dynamics, allosteric motion,
alignment choices, and periodic boundary artifacts can all produce oscillatory
patterns.

For oscillating systems, the range and timescale of the oscillation often convey
more than the mean RMSD alone.

## How PolyzyMD Handles Autocorrelation

RMSD timeseries are autocorrelated because adjacent MD frames are not
independent samples. A trajectory with many saved frames can still contain far
fewer statistically independent observations.

PolyzyMD reports uncertainty in terms of statistical inefficiency where
possible. Conceptually, the effective sample size is `N_eff = N / g`, where `g`
is the statistical inefficiency. For a simple integrated autocorrelation-time
estimate, `g ≈ 1 + 2τ/dt`, with `τ` the integrated autocorrelation time and `dt`
the frame spacing. Larger `g` means stronger correlation and fewer effective
samples.

This correction helps avoid treating adjacent frames as independent, but it
does not replace independent replicate simulations or guarantee convergence of
the underlying conformational ensemble.

```{seealso}
For the mathematical details of autocorrelation functions and the LiveCoMS
recommendations, see the
[Statistics Best Practices Guide](analysis_statistics_best_practices.md).
```

## Multi-Run Analysis: Why It Helps Interpretation

Different RMSD selections answer different questions:

| Run Label | Selection | Question |
|-----------|-----------|----------|
| "Protein Backbone" | `protein and name CA` | How close is the global backbone to this reference? |
| "Active Site" | Catalytic residues CA | How close is the local active-site geometry to this reference? |
| "Polymer Core" | `chainid C and not name H*` | How close is the polymer conformation to this reference? |
| "Crystal Deviation" | `protein and name CA` (external ref) | How close is the protein to an external structural state? |

Each run is ranked independently across conditions. This prevents averaging
RMSD from structurally different selections, which would be difficult to
interpret:

```text
Rankings:
  Protein Backbone: With Polymer < No Polymer (closer to reference)
  Active Site:      With Polymer < No Polymer (closer to reference)
  Polymer Core:     No Polymer — (single condition only)
```

## External Reference for Catalytic Competence

When studying enzyme catalysis across multiple conditions, the standard
reference modes (`centroid`, `average`) use a **condition-specific** reference:
each condition's trajectory determines its own reference structure.

The `external` reference mode uses a **condition-independent** reference,
typically a crystal structure representing a specific geometry of interest. RMSD
then measures deviation from that external structure:

$$
\text{RMSD}^{\text{ext}}(t) = \sqrt{\frac{1}{N} \sum_{i=1}^{N} \left\| \mathbf{r}_i(t) - \mathbf{r}_i^{\text{crystal}} \right\|^2}
$$

**Interpretation changes with external reference:**

| Metric | Standard RMSD (centroid/average) | External Reference RMSD |
|--------|----------------------------------|------------------------|
| Low value | Structure stays near its own trajectory-derived reference | Structure stays near the external geometry |
| High value | Structure deviates from its trajectory-derived reference | Structure deviates from the external geometry |
| Condition comparison | Which condition remains closer to its chosen internal reference? | Which condition remains closer to the external structure? |

```{tip}
**Which reference mode for enzymes?** Use `centroid` or `average` for
trajectory-internal reference-relative motion. Use `external` to ask whether a
trajectory remains close to a specific known structure. External-reference RMSD
does not make "closer" inherently better or more stable unless the external
structure is justified as the relevant state for the scientific question.
```

## Replicates vs Longer Simulations

### The LiveCoMS recommendation

> "Multiple independent simulations are preferable to a single long simulation"
> — Grossfield et al. (2018)

### Why replicates matter for RMSD

| Multiple Replicates | Single Long Simulation |
|--------------------|------------------------|
| Independent starting points | Frames remain correlated |
| Tests reproducibility of drift/plateau patterns | May remain trapped in one metastable state |
| Supports uncertainty from replicate means | Requires autocorrelation correction within trajectory |
| Parallelizable | Sequential |

```{note}
With only 1 replicate, PolyzyMD still computes RMSD and includes the condition
in descriptive summaries and rankings. Replicate SEM is unavailable because
variability across independent simulations cannot be estimated from a singleton.
Pairwise inferential tests require at least 2 replicates per condition.
```

## Comparing Conditions

### What PolyzyMD computes

For each RMSD run, the comparison produces:

| Statistic | Description |
|-----------|-------------|
| **Ranking** | Conditions sorted by mean RMSD (lowest = closest to the chosen reference) |
| **Percent change** | Relative to control condition |
| **Direction** | Plugin labels such as `stabilizing`, `destabilizing`, or `unchanged`; interpret as reference-relative unless separately justified |
| **t-statistic** | Two-sample t-test on replicate means |
| **p-value** | Two-tailed significance |
| **Cohen's d** | Effect size magnitude |
| **ANOVA** | Omnibus F-test when 3+ conditions (per-run) |

### Direction labels

PolyzyMD classifies the direction of change based on percent change in mean RMSD
relative to control. For RMSD, these labels are shorthand and should be read as
changes in closeness to the chosen reference, not proof of biological stability.

| Percent Change | Direction | Meaning |
|---------------|-----------|---------|
| < −1% | `stabilizing` | Treatment reduces reference-relative deviation |
| > +1% | `destabilizing` | Treatment increases reference-relative deviation |
| −1% to +1% | `unchanged` | No meaningful difference by this threshold |

### Interpreting the comparison

When one condition has lower mean RMSD than another, the most direct statement
is that it stayed closer to the chosen reference for the selected atoms over the
analyzed interval. Stronger claims, such as improved stability or functional
preservation, require supporting evidence from the scientific context and other
observables.

PolyzyMD writes canonical RMSD artifacts through the analysis lifecycle. The
stable locations are:

- `analysis/<sanitized_condition_label>/rmsd/run_<N>/result.json` for
  replicate-level artifacts
- `analysis/<sanitized_condition_label>/rmsd/aggregated/result.json` for
  condition-level artifacts
- `comparison/rmsd/result.json` for comparison artifacts

Treat artifact contents as structured payloads and provenance that may refer to
sidecars for larger data. Avoid depending on undocumented raw JSON field names
unless they are described in reference documentation.

## Common Pitfalls

### 1. Treating a plateau as proof of equilibration

**Symptom:** A plateau-like RMSD trace is described as complete equilibration.

**Caution:** A plateau suggests stationarity of the selected atoms relative to
the chosen reference. Other coordinates, slow modes, ligand states, solvent
structure, or functional observables may still be unequilibrated.

**Better interpretation:** "RMSD reached an apparent plateau for this selection
and reference after the initial relaxation period."

### 2. Comparing different selections

**Symptom:** RMSD values are not comparable across runs or publications.

**Cause:** Different atom selections yield different RMSD magnitudes.

**Better interpretation:** Always report the exact selection string. Compare
only runs with identical selections, references, and alignment conventions.

### 3. Over-interpreting small differences

**Symptom:** Claiming significance for 0.05 Å differences.

**Cause:** Not accounting for uncertainty.

**Better interpretation:** Report uncertainty and avoid implying meaningful
structural differences when confidence intervals overlap substantially or
replicate variation dominates:

```text
# WRONG: "Condition A (1.856 Å) is less stable than B (1.861 Å)"
# RIGHT: "Condition A (1.856 ± 0.034 Å) and B (1.861 ± 0.028 Å)
#         are not significantly different (p = 0.91, unchanged)"
```

### 4. Ignoring timeseries shape

**Symptom:** Reporting only mean RMSD without inspecting the timeseries.

**Cause:** Two conditions can have the same mean RMSD but very different
dynamics, such as one plateau-like trace and one drifting trace.

**Better interpretation:** Inspect the timeseries shape before reducing the
trajectory to a mean. Similar means can arise from stationary, drifting, or
multi-regime trajectories.

### 5. Using all-atom RMSD without justification

**Symptom:** Very high RMSD values even for compact proteins.

**Cause:** Side-chain motions can dominate all-atom RMSD, obscuring backbone
changes.

**Better interpretation:** Use Cα, backbone, all-atom, or local selections
according to the scientific question. Side-chain-rich selections are valid when
side-chain rearrangements are the intended observable, but their RMSD scale is
not interchangeable with Cα RMSD.

### 6. Ignoring replicate variation

**Symptom:** Reporting within-trajectory SEM as the total uncertainty.

**Cause:** Treating autocorrelation-corrected SEM as sufficient.

**Better interpretation:** Use independent replicate statistics when available.
Within-trajectory uncertainty can account for adjacent-frame correlation, but
replicate-to-replicate variability better reflects sensitivity to initial
conditions and sampling path.

### 7. Choosing the wrong reference mode

**Symptom:** Unexpected or hard-to-interpret comparison results.

**Cause:** Using `centroid` when `external` is more appropriate for the
scientific question, or interpreting external-reference RMSD as inherently
better when it is merely closer to the supplied structure.

**Better interpretation:** Match reference mode to your scientific question:

- Trajectory-internal reference-relative motion → `centroid` or `average`
- Closeness to a specified structural state → `external` with a justified
  reference structure

### 8. Treating automated convergence as ground truth

**Symptom:** Trusting an automated convergence diagnostic without further
inspection.

**Cause:** A sliding-window heuristic is parameter-dependent and can miss slow
drift, metastable trapping, or convergence issues in observables other than
RMSD.

**Better interpretation:** Use convergence diagnostics as one input among
several. Inspect the RMSD timeseries, run multiple independent replicates when
possible, and check other relevant observables such as Rg, SASA, contacts, or
active-site distances. See {doc}`/explanation/convergence_detection` for a full
discussion of limitations.

## RMSD as one equilibration diagnostic

RMSD is commonly used as an equilibration diagnostic because large structural
relaxations often appear as changes in reference-relative distance. Its role is
diagnostic, not definitive. A plateau can support the claim that the selected
atoms are no longer drifting relative to the reference on the observed
timescale, but it does not establish thermodynamic equilibration or convergence
of all relevant observables.

```{tip}
If RMSD never appears stationary within the simulation time, possible
explanations include slow relaxation, reference mismatch, large-amplitude domain
motion, unfolding, or simply insufficient sampling. Distinguish these by
inspecting structures and complementary observables.
```

### Automated convergence detection

```{versionadded} 1.3.0
```

PolyzyMD can run a sliding-window convergence diagnostic on RMSD timeseries. The
diagnostic evaluates whether reference-relative RMSD changes remain below a
configured threshold over a sustained interval. The resulting information is
stored as part of the canonical RMSD artifact payload and provenance, with
condition-level summaries represented in aggregated artifacts. Larger timeseries
or plot-ready data may be represented through sidecars referenced by the
artifact.

**This is a diagnostic tool, not a definitive convergence proof.** The
heuristic can miss slow drift below the slope threshold, and convergence in
RMSD does not guarantee convergence of other observables. Always use multiple
replicates and visual inspection alongside automated diagnostics.

For command-oriented usage, see the
[RMSD Quick Start Guide](../how_to/analysis_rmsd_quickstart.md). For a full
conceptual treatment of convergence diagnostics — including the algorithm,
parameters, tuning guidance, and limitations — see
{doc}`/explanation/convergence_detection`.

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
