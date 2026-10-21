# Catalytic triad: interpretation and best practices

A catalytic triad measurement summarizes active-site geometry from MD
trajectories. It is useful for asking whether a serine protease, lipase, or
esterase active site tends to preserve the geometric arrangement associated
with catalysis.

It is not direct evidence of catalytic activity. The hydrogen-bond fractions,
distances and contact fractions are **geometric proxies**. Interpret them
together with the replicate uncertainty, the substrate pose, the hydrogen-bond
geometry, the protonation state and any experimental activity data.

```{note}
PolyzyMD has no separate triad analysis: the triad is a routine on the analysis
API. For the code, see {doc}`../how_to/analysis_triad_quickstart`. For what each
function measures, see {doc}`../reference/analysis_functions`.
```

## What the metric represents

A classical Ser-His-Asp/Glu catalytic triad depends on a hydrogen-bond network
that helps position histidine and activate the serine nucleophile. PolyzyMD does
not model reactivity or proton transfer. The routine of
{doc}`../how_to/analysis_triad_quickstart` measures the network in two ways.

**Hydrogen bonds.** `functions.hbond_count(group_a, group_b)` counts, on every
frame, the hydrogen bonds between two groups with MDAnalysis
`HydrogenBondAnalysis`: donor within 3.5 Å of the acceptor and a
donor-hydrogen-acceptor angle of at least 150° by default. One `study.timeseries`
call per triad bond, such as Ser OG-HG to His NE2 and His ND1-HD1 to Asp OD1 or
OD2, gives each bond's count per frame, and `Timeseries.transform` combines
them into a series that is 1 on the frames where every bond is formed.

**Distances.** `functions.pair_distance` measures a heavy-atom or
point-to-point distance per frame, and `polyzymd analyze distances --set
pairs=triad.yaml` measures every listed pair with the fraction of frames below
its threshold. `functions.all_below` combines stored distance series into a
series that is 1 on the frames where every pair is below its threshold.

Either way, the mean of the combined series over one replicate's production
frames is the **simultaneous fraction**:

$$
f_{\text{simultaneous}} = \frac{1}{N} \sum_{t=1}^{N} \prod_{i=1}^{M} \mathbb{1}[\text{pair } i \text{ formed at frame } t]
$$

where $N$ is the number of analyzed frames and $M$ the number of bonds or pairs.
A pair is formed when its hydrogen bond is counted, or when its distance
$d_i(t)$ is below the threshold $\theta_i$.

The simultaneous fraction is stricter than per-pair fractions. Two pairs can
each be formed 50% of the time but never in the same frames. In that case, the
simultaneous fraction is 0%, which is often more relevant to an intact triad
interpretation than either per-pair fraction alone.

## Running the routine

{doc}`../how_to/analysis_triad_quickstart` gives the code. Each replicate
contributes one value per result: the simultaneous fraction, from 0 to 1, and
each bond's fraction or each pair's mean distance and fraction below threshold.
`compare()` compares every condition with the control by Welch's t test on those
replicate values, with the Benjamini-Hochberg correction over the conditions
compared. The per-frame series of every replicate are stored with a record of
the selections, thresholds and input files, so combining them again, or
reporting another result, does not read the trajectory again.

```{tip}
Decide before looking at the results which outcome carries your conclusion.
For a claim about catalytic competence, make the simultaneous fraction the
primary outcome and read the individual pairs as supporting evidence. The
simultaneous fraction already combines the pairs, so it needs no correction
across them, and each outcome is corrected over its own condition comparisons.
Reporting whichever outcome happens to give the smallest p value is a form of
selective reporting that no correction repairs.
```

## Thresholds are heuristics, not activity cutoffs

The default threshold of 3.5 Å is a practical heavy-atom cutoff for hydrogen-bond
like contacts. It is not a universal boundary between active and inactive
enzyme states.

Useful ways to think about threshold choices:

- **Around 3.0 Å** is stricter and emphasizes close, well-formed contacts.
- **Around 3.5 Å** is a common heavy-atom proxy for hydrogen-bond-like contact.
- **Around 4.0 Å** is more permissive and may include weak or transient
  interactions.

Choose the thresholds before you compare conditions. Avoid tuning a
threshold after seeing the results just to make a preferred condition look
active or inactive. That kind of post-hoc threshold selection makes the metric
circular and can overstate the evidence.

When a threshold is uncertain, report how the result changes with it. For
example, state whether the same qualitative ordering appears at 3.0, 3.5, and
4.0 Å, rather than selecting only the cutoff that gives the clearest story.

## Heavy-atom distances are only hydrogen-bond proxies

Triad distances are usually measured between heavy atoms or user-defined
points such as `midpoint(...)`. This needs no hydrogens. It is only a proxy
for hydrogen bonding. `hbond_count` also checks the angle at the hydrogen.

Important cautions:

- A short N···O or O···O heavy-atom distance does not guarantee a productive
  hydrogen bond; angle and donor-hydrogen placement matter.
- A distance slightly above the threshold does not prove the active site is
  catalytically inactive; transient geometry, force-field behavior, and sampling
  limitations can all matter.
- The histidine tautomer and protonation state decide which nitrogen to
  measure, and how to read the Ser-His and Asp/Glu-His contacts.
- Asp/Glu atom choices matter. A midpoint of the carboxylate atoms can be useful
  for symmetric monitoring, but it is not the same as tracking a specific
  oxygen involved in a particular hydrogen bond.
- Substrate pose matters. A preserved Ser-His-Asp/Glu geometry is more
  convincing when the substrate is also positioned consistently with the
  proposed mechanism.

For mechanistic claims, combine the triad metric with direct inspection of
active-site snapshots, the hydrogen-bond counts of `hbond_count`, substrate
distance/orientation analyses, and experiment.

## Replicate uncertainty matters more than frame count

Frame-level contact states are temporally correlated. A trajectory with many
closely spaced frames does not provide the same evidence as many independent
samples. PolyzyMD summarizes replicate-level values and condition-level
uncertainty so comparisons are not based solely on frame counts.

Best interpretive practice:

- Prefer at least three independent replicates per condition for conclusions.
- Treat one-replicate output as descriptive or suitable for smoke tests, not as
  strong comparative evidence.
- Interpret large replicate-to-replicate variation as a signal that the active
  site may occupy multiple metastable states or that more sampling is needed.
- Compare conditions using replicate summaries and uncertainty, not raw frame
  counts.

See {doc}`analysis_statistics_best_practices` for the broader statistical
context.

## Reading common result patterns

### High simultaneous contact, low uncertainty

This is consistent with a stable triad geometry under the chosen selections and
threshold. It is strongest when per-pair distances are also reasonable, substrate
pose is compatible with catalysis, and replicates agree.

### High per-pair contact but low simultaneous contact

This suggests the contact network is not intact in the same frames. It may
indicate alternating conformational states or a flexible active site. Per-pair
plots and distance distributions are more informative than the scalar metric
alone.

### Low contact dominated by one pair

This often points to a specific disrupted interaction, incorrect atom choice, or
incorrect residue numbering. It is a diagnostic clue, not by itself proof of a
mechanistic cause.

### Large differences with broad uncertainty

Large apparent effects can be hypothesis-generating even when uncertainty is
high, but they should be described cautiously. Strong conclusions require
replicate support and ideally orthogonal evidence.

## Worked example interpretation: keep conclusions tentative

Suppose a LipA polymer study reports that a pure EGMA condition has a higher
simultaneous contact fraction than several mixed-polymer conditions, while the
mixed conditions show one or both pair distances shifted upward.

A cautious interpretation would be:

- The simulations suggest that EGMA may preserve the monitored triad geometry
  better than the mixed-polymer conditions under the chosen model and threshold.
- The mixed conditions may sample active-site geometries in which the monitored
  hydrogen-bond proxy distances are less often simultaneously close.
- If one pair, such as Asp-His, is especially shifted, that pair is a useful
  target for structural inspection.

Avoid stronger conclusions unless they are supported by the full evidence base.
For example, do not state that a polymer composition "preserves activity" or
"disrupts catalysis" from the contact fraction alone. Those claims need
replicate uncertainty, substrate pose consistency, hydrogen-bond geometry, and
experimental activity or other mechanistic validation.

## Figures

`polyzymd analyze distances` draws, for each pair, the distribution of its
per-frame distance in every condition, pooled over replicates with one thin
curve per replicate and the threshold marked (`distance_kde_<pair>`), a bar chart
of each fraction with every replicate value shown (`distance_fraction_<result>`),
`distance_threshold_bars`, every pair's fraction in one chart, and
`distance_kde_panel`, one distribution panel per pair. In Python,
`series.plot_distribution()` draws a stored series' distribution and
`series.reduce("mean").plot()` the condition means with every replicate value,
for the simultaneous fraction too. Use the distributions to see whether a
difference in a fraction comes from a broad shift, a small subpopulation, or one
limiting pair. The figures are interpretive aids; they do not replace
statistical uncertainty or structural validation.

## Common interpretation pitfalls

**Treating geometry as activity.**
: A preserved triad geometry is compatible with activity, but activity also
  depends on substrate binding, chemical step feasibility, solvent, protonation,
  and other factors.

**Using bare residue numbers without checking the topology.**
: Residue numbering and chain assignment can differ across prepared systems.
  Prefer chain-aware or protein-restricted selections and verify atom names.

**Choosing atom names without considering histidine chemistry.**
: `ND1` and `NE2` have different roles depending on tautomer/protonation and
  enzyme family. Confirm that the selected nitrogen matches the intended
  interaction.

**Overfitting the threshold.**
: A threshold chosen after inspecting condition rankings can make the result
  circular. Predefine thresholds or report a sensitivity analysis.

**Ignoring pair-level diagnostics.**
: The simultaneous fraction is compact but lossy. Always inspect which pair or
  distribution drives a change before making a mechanistic claim.

## References

**Hedstrom L.** (2002) "Serine Protease Mechanism and Specificity."
*Chemical Reviews* 102:4501-4524. https://doi.org/10.1021/cr000033x

**Blow DM.** (1976) "Structure and Mechanism of Chymotrypsin."
*Accounts of Chemical Research* 9:145-152.

**Grossfield A, Patrone PN, Roe DR, Schultz AJ, Siderius DW, Zuckerman DM.**
(2018) "Best Practices for Quantification of Uncertainty and Sampling Quality
in Molecular Simulations." *Living Journal of Computational Molecular Science*
1(1):5067. https://doi.org/10.33011/livecoms.1.1.5067

**Jeffrey GA, Saenger W.** (1991) *Hydrogen Bonding in Biological Structures.*
Springer-Verlag.

## See also

- {doc}`../how_to/analysis_triad_quickstart`
- {doc}`../reference/analysis_functions`
- {doc}`analysis_statistics_best_practices`
- {doc}`../how_to/analysis_compare_conditions`
