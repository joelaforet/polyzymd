# RMSF Analysis: Statistical Best Practices

Root mean square fluctuation (RMSF) is useful for asking where a protein is
more rigid or flexible, but it is easy to over-interpret. This page explains
how PolyzyMD treats RMSF statistically, how to interpret warnings, and what
current RMSF artifacts mean for contributors.

```{note}
**Need commands rather than interpretation guidance?** See the
[RMSF quickstart](../how_to/analysis_rmsf_quickstart.md) for copy-paste CLI
examples and minimal setup.
```

```{seealso}
For foundational concepts such as autocorrelation, correlation time, and
statistical inefficiency, see
[Statistics Best Practices](analysis_statistics_best_practices.md). This page
focuses on RMSF-specific interpretation.
```

## What RMSF measures

RMSF measures how much each atom or residue fluctuates around a reference
position during a trajectory:

$$
\text{RMSF}_i = \sqrt{\frac{1}{T} \sum_{t=1}^{T} \left( \mathbf{r}_i(t) - \langle \mathbf{r}_i \rangle \right)^2}
$$

where $\mathbf{r}_i(t)$ is the position of atom or residue $i$ at time $t$,
$\langle \mathbf{r}_i \rangle$ is its time-averaged position, and $T$ is the
number of frames used in the estimate.

This is the standard PolyzyMD interpretation for non-external reference modes:
`centroid`, `average`, and `frame` affect alignment/reference generation, but
RMSF is computed as fluctuation around the mean positions of the aligned
analyzed trajectory. Those modes should not be read as direct deviation from a
centroid frame or a selected trajectory frame. See
[reference structure selection](analysis_reference_selection.md) for the
mode-by-mode interpretation.

Low RMSF often indicates a relatively rigid region, such as a buried core or
structured secondary element. High RMSF often indicates a flexible region, such
as a loop, terminus, or mobile binding-site element. These are interpretations,
not automatic conclusions: RMSF depends on alignment choice, reference choice,
atom selection, equilibration, force field, and sampling quality.

RMSF is related to crystallographic B-factors by

$$
B_i = \frac{8\pi^2}{3} \langle \Delta r_i^2 \rangle = \frac{8\pi^2}{3} \text{RMSF}_i^2
$$

but crystal packing, experimental model refinement, and simulation conditions
mean the correspondence is approximate.

## Why correlation matters for RMSF

MD trajectories contain correlated frames. Correlation affects how quickly RMSF
estimates converge and how much uncertainty should be assigned to them. It is
not enough to count every saved frame as an independent observation.

PolyzyMD's current RMSF strategy is conservative: it estimates a correlation
time and subsamples approximately independent frames before computing RMSF. This
reduces the risk of treating dense, correlated trajectory output as more
informative than it is. It also means short trajectories with long correlation
times may produce only a small number of effective samples.

Conceptually, PolyzyMD:

1. aligns the trajectory to the chosen reference;
2. estimates correlation from a trajectory-level timeseries;
3. selects frames spaced far enough apart to be treated as approximately
   independent for the current RMSF calculation;
4. computes RMSF from that reduced frame set.

Example diagnostic output may look like this:

```text
Correlation time: 15394 ps (15.4 ns)
Statistical inefficiency: 308.9
Independent samples: 6 (from 2000 frames)
```

This does not mean the trajectory is invalid. It means the RMSF estimate has
less independent information than the raw frame count suggests.

## Interpreting reliability warnings

PolyzyMD warns when the effective number of independent samples is small, for
example:

```text
WARNING: Low statistical reliability: only 6 independent samples
(recommended >= 10). Correlation time τ = 15394 ps is comparable to
or longer than the trajectory sampling window. Consider:
(1) extending simulation time,
(2) using multiple independent trajectories, or
(3) interpreting results with caution.
See Grossfield et al. (2018) LiveCoMS 1:5067.
```

Treat this as a sampling and uncertainty warning. It does not by itself prove
that a simulation is broken, but it should make you ask whether the conclusion
depends on poorly converged fluctuations.

Useful follow-up questions include:

- Do independent replicates show similar RMSF patterns?
- Are replicate means stable, or does one replicate dominate the conclusion?
- Is the trajectory stationary after the equilibration period, or do RMSD/RMSF
  summaries drift over time?
- Is the claimed effect large compared with replicate-to-replicate variation?
- Is the conclusion qualitative, or does it require a precise uncertainty
  estimate?

Replicates can make an RMSF result much more credible, but they are not a magic
fix. Different random seeds or initial velocities help explore independent
trajectory histories, yet they do not guarantee independent equilibrium
sampling if all simulations remain trapped in the same metastable basin or if
equilibration is incomplete.

## Replicates and trajectory length

LiveCoMS-style guidance generally favors multiple independent simulations over
placing all sampling effort into one long trajectory, especially when estimating
uncertainty. For RMSF, this is useful because replicate-to-replicate variation
shows whether the observed flexibility pattern is reproducible.

Multiple replicates help because they:

- test reproducibility across independently initialized simulations;
- reveal outlier trajectories or rare conformational events;
- provide condition-level uncertainty from replicate summaries;
- can be run in parallel.

A single longer trajectory can still be valuable, especially for slow processes
that are not reached in shorter runs. The right balance depends on the system,
expected timescales, and the scientific claim. For many enzyme-polymer studies,
3-5 replicates per condition is a practical starting point, not a universal
rule.

## Incomplete data and current artifacts

PolyzyMD can aggregate RMSF results when some requested replicates are missing
or fail analysis. One successful replicate may be useful for descriptive checks
or smoke tests, but between-replicate SEM and inferential comparisons require
multiple successful replicates. If only a subset of planned replicates is
available, interpret the result as provisional and document which replicates
contributed.

Current RMSF outputs use the analysis artifact lifecycle rather than legacy
top-level filenames. The stable entry points are:

- per-replicate result:
  `analysis/<sanitized_condition_label>/rmsf/run_<N>/result.json`
- per-condition aggregate:
  `analysis/<sanitized_condition_label>/rmsf/aggregated/result.json`
- cross-condition comparison:
  `comparison/rmsf/result.json`

Large arrays, per-residue tables, or other bulky data may live in sidecar files
referenced by these artifacts. Treat the artifact JSON files as the stable
entry points for consumers; do not assume every detailed array is embedded as a
top-level JSON field.

## Comparing conditions

RMSF comparisons should be based on replicate-level summaries, not on treating
all frames as independent observations. A typical interpretation workflow is:

1. compute per-replicate RMSF summaries for each condition;
2. aggregate those summaries within each condition;
3. compare condition-level distributions using replicate-level values;
4. interpret p-values together with effect size, uncertainty, and physical
   plausibility.

With small sample sizes, a large apparent effect can coexist with a
non-significant p-value. This means the current data are not sufficient to
reject the null hypothesis at the chosen threshold; it does not prove there is
no effect. Additional replicates may clarify reproducibility and uncertainty,
but they do not guarantee statistical significance.

When reviewing RMSF differences, prefer cautious language:

- "The polymer condition shows lower mean RMSF in these replicates" rather
  than "the polymer stabilizes the enzyme".
- "The effect is suggestive but uncertain" rather than "more replicates would
  make it significant".
- "Replicate 2 samples a different state" rather than "replicate 2 is bad",
  unless there is a documented technical failure.

## Interpreting RMSF magnitudes

The following ranges are rough heuristics for Cα RMSF in folded proteins under
typical simulation conditions. They are not universal thresholds and should not
be applied blindly to all atoms, intrinsically disordered regions, nucleic
acids, polymers, ligands, or externally referenced deviation metrics.

| Approximate Cα RMSF | Common interpretation |
|---------------------|-----------------------|
| 0.3-0.5 Å | Very rigid folded core or constrained secondary structure |
| 0.5-1.0 Å | Moderate flexibility in structured regions |
| 1.0-2.0 Å | Flexible loops, flaps, or exposed regions |
| 2.0-5.0 Å | Highly mobile termini or disordered segments |
| >5.0 Å | Possible disorder, unfolding, poor alignment, or reference mismatch |

Always interpret these magnitudes alongside structure, alignment selection,
temperature, solvent, force field, and replicate behavior.

For enzyme active sites, lower RMSF is not automatically better. A rigid active
site may preserve catalytic geometry, but some enzymes require conformational
breathing, induced fit, or loop motion. Active-site RMSF is most useful when
combined with geometry-specific analyses, substrate positioning, and
experimental activity data.

## External-reference RMSF-like deviations

PolyzyMD supports reference modes that change the scientific meaning of the
reported values. Standard non-external RMSF asks how much a residue fluctuates
around the aligned trajectory mean; `centroid`, `average`, and `frame` change how
the trajectory is aligned and how the alignment/reference structure is generated.
External mode is the special fixed-reference path: it uses mapped external
coordinates as the RMSF reference positions and measures an RMSF-like
per-residue deviation from a fixed external structure, such as a crystal model of
a catalytically competent state.

That external-reference quantity is useful, but it is not the same as standard
RMSF around the aligned trajectory mean. A low value means the residue remains
close to the chosen external structure; a high value means it departs from that
structure. This can be appropriate for questions about maintaining catalytic
geometry, while standard non-external RMSF is better for questions about
flexibility within each sampled ensemble.

Use both views when they answer different questions:

- standard non-external RMSF: "Which regions are flexible in this condition
  after the chosen alignment?"
- external-reference deviation: "Which condition remains closest to a chosen
  functional structure?"

For setup details, see the
[external PDB section](analysis_reference_selection.md#external-pdb-fixed-reference-rmsf-like-deviation) of the
reference selection guide.

## Common interpretation pitfalls

### Treating correlated frames as independent

The raw number of saved frames is not the number of independent samples. Use
replicate-level summaries and PolyzyMD's artifact outputs rather than computing
SEM from every frame as if each were independent.

### Ignoring stationarity

If RMSD or structural summaries drift after the equilibration cutoff, RMSF may
combine multiple regimes into one number. In that case, the main issue is not
only uncertainty but whether the analyzed window represents a stable ensemble.

### Over-interpreting small differences

Differences of a few hundredths of an Å can be smaller than uncertainty from
replicate variation, alignment choices, or reference selection. Report
uncertainty and avoid mechanistic conclusions from tiny differences alone.

### Cherry-picking replicates

Exclude a replicate only for a documented technical reason, such as a corrupted
trajectory or failed simulation. A conformational transition or ligand unbinding
event may be scientifically important rather than an error.

### Comparing incompatible reference definitions

Do not compare standard trajectory-mean RMSF and external-reference deviation as
if they were the same metric. They answer different questions and should be
labeled accordingly.

## Contributor notes

Contributors should keep RMSF documentation and downstream tooling aligned with
the artifact lifecycle. Programmatic consumers should start from canonical
artifact paths and follow payload or sidecar references rather than relying on
legacy filenames or stale top-level fields.

For default comparison results, condition summaries may be available under a
comparison artifact payload such as
`ComparisonArtifact.payload["condition_summaries"]`. The exact payload shape can
evolve with the analysis framework, so contributors should prefer documented
artifact contracts and public analysis APIs over private module imports.

## References

**Grossfield A, Patrone PN, Roe DR, Schultz AJ, Siderius DW, Zuckerman DM.**
(2018) "Best Practices for Quantification of Uncertainty and Sampling Quality
in Molecular Simulations." *Living Journal of Computational Molecular Science*
1(1):5067. https://doi.org/10.33011/livecoms.1.1.5067

**Flyvbjerg H, Petersen HG.** (1989) "Error estimates on averages of correlated
data." *Journal of Chemical Physics* 91:461-466.
https://doi.org/10.1063/1.457480

**Chodera JD, Swope WC, Pitera JW, Seok C, Dill KA.** (2007) "Use of the
Weighted Histogram Analysis Method for the Analysis of Simulated and Parallel
Tempering Simulations." *Journal of Chemical Theory and Computation* 3:26-41.
https://doi.org/10.1021/ct0502864

**Knapp B, Frantal S, Greshake B, Schwarz R, et al.** (2018) "Is an Intuitive
Convergence Definition of Molecular Dynamics Simulations Solely Based on the
Root Mean Square Deviation Possible?" *Journal of Computational Biology*
25:1069-1077.

## See Also

- [RMSF quickstart](../how_to/analysis_rmsf_quickstart.md) — commands and minimal setup
- [RMSF implementation verification](analysis_rmsf_verification.md) — benchmark against MDAnalysis and GROMACS
- [Statistics best practices](analysis_statistics_best_practices.md) — foundational statistics for MD
- [Reference structure selection](analysis_reference_selection.md) — choose alignment and reference modes
- [LiveCoMS best practices](https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5067) — methodology paper
