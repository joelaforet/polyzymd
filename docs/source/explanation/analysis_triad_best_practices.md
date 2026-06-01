# Catalytic triad analysis: interpretation and best practices

Catalytic triad analysis in PolyzyMD summarizes active-site geometry from MD
trajectories. It is useful for asking whether a serine protease, lipase, or
esterase active site tends to preserve the geometric arrangement associated
with catalysis.

It is not direct evidence of catalytic activity. The reported distances and
contact fractions are **geometric proxies**. They should be interpreted with
replicate uncertainty, substrate pose, hydrogen-bond geometry, protonation
state, and experimental activity data whenever those are available.

```{note}
For command examples and setup steps, see the
{doc}`../how_to/analysis_triad_quickstart`. For field-level lookup, see
{doc}`../reference/analysis_triad_reference`.
```

## What the metric represents

A classical Ser-His-Asp/Glu catalytic triad depends on a hydrogen-bond network
that helps position histidine and activate the serine nucleophile. PolyzyMD does
not model reactivity or proton transfer directly in this analysis. Instead, it
tracks user-defined heavy-atom or point-to-point distances such as:

- Asp/Glu carboxylate to His nitrogen
- His nitrogen to Ser hydroxyl oxygen
- other active-site distances chosen for a specific enzyme family

For each trajectory frame, the plugin asks whether every configured pair is
below the contact threshold. The main scalar metric is the **simultaneous
contact fraction**:

$$
f_{\text{contact}} = \frac{1}{N} \sum_{t=1}^{N} \prod_{i=1}^{M} \mathbb{1}[d_i(t) < \theta]
$$

where $N$ is the number of analyzed frames, $M$ is the number of configured
pairs, $d_i(t)$ is the distance for pair $i$ at frame $t$, and $\theta$ is the
distance threshold.

The simultaneous fraction is stricter than per-pair contact fractions. Two
pairs can each be in contact 50% of the time but never be in contact in the
same frames. In that case, the simultaneous contact fraction is 0%, which is
often more relevant to an intact triad interpretation than either per-pair
fraction alone.

## Configuration shape and plugin name

The plugin name is `catalytic_triad`. Use this canonical name in CLI commands
and `comparison.yaml` configuration.

Current `comparison.yaml` files define triad settings under `plugins:`:

```yaml
plugins:
  catalytic_triad:
    name: "LipA_catalytic_triad"
    description: "Ser-His-Asp catalytic triad of Lipase A"
    threshold: 3.5
    pairs:
      - label: "Asp133-His156"
        selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"
        selection_b: "protein and resid 156 and name ND1"
      - label: "His156-Ser77"
        selection_a: "protein and resid 156 and name NE2"
        selection_b: "protein and resid 77 and name OG"
```

Older examples that used a top-level `catalytic_triad:` block are stale. Keep
plugin settings under `plugins.catalytic_triad` so the comparison workflow can
discover and configure the plugin consistently.

## Artifact lifecycle and output interpretation

Catalytic triad analysis uses the current PolyzyMD analysis artifact lifecycle:

1. Per-replicate MDAnalysis jobs compute distance profiles and contact metrics.
2. Per-replicate artifacts are written for each condition and replicate.
3. Condition aggregation combines replicate artifacts without re-reading
   trajectories.
4. Cross-condition comparison reads condition artifacts and writes a comparison
   artifact.
5. Plotting reads cached artifacts and sidecars; it should not rerun the
   trajectory analysis.

Canonical artifact paths are:

- Per replicate:
  `analysis/<sanitized_condition_label>/catalytic_triad/run_<N>/result.json`
- Per condition:
  `analysis/<sanitized_condition_label>/catalytic_triad/aggregated/result.json`
- Cross-condition comparison:
  `comparison/catalytic_triad/result.json`

Condition labels from `comparison.yaml` are sanitized before they become
filesystem path components. For example, a label such as `75% SBMA / 25% EGMA`
is written under a filesystem-safe directory name rather than the literal label.
Use the label stored inside the artifact when you need the human-readable
condition name.

Large arrays, such as distance time series or aggregated distance profiles, may
be stored in artifact sidecars below `sidecars/`. Treat `result.json` as the
entry point and sidecars as validated data referenced by the artifact.

## Thresholds are heuristics, not activity cutoffs

The default threshold of 3.5 Å is a practical heavy-atom cutoff for hydrogen-bond
like contacts. It is not a universal boundary between active and inactive
enzyme states.

Useful ways to think about threshold choices:

- **Around 3.0 Å** is stricter and emphasizes close, well-formed contacts.
- **Around 3.5 Å** is a common heavy-atom proxy for hydrogen-bond-like contact.
- **Around 4.0 Å** is more permissive and may include weak or transient
  interactions.

Choose thresholds before comparing conditions whenever possible. Avoid tuning a
threshold after seeing the results just to make a preferred condition look
active or inactive. That kind of post-hoc threshold selection makes the metric
circular and can overstate the evidence.

When a threshold is uncertain, report sensitivity analyses honestly. For
example, note whether the same qualitative ordering appears at 3.0, 3.5, and
4.0 Å, rather than selecting only the cutoff that gives the clearest story.

## Heavy-atom distances are only hydrogen-bond proxies

PolyzyMD triad distances are usually measured between heavy atoms or
user-defined points such as `midpoint(...)`. This is robust and convenient, but
it is only a proxy for hydrogen bonding.

Important cautions:

- A short N···O or O···O heavy-atom distance does not guarantee a productive
  hydrogen bond; angle and donor-hydrogen placement matter.
- A distance slightly above the threshold does not prove the active site is
  catalytically inactive; transient geometry, force-field behavior, and sampling
  limitations can all matter.
- Histidine tautomer and protonation state affect which nitrogen should be used
  and how the Ser-His and Asp/Glu-His contacts should be interpreted.
- Asp/Glu atom choices matter. A midpoint of the carboxylate atoms can be useful
  for symmetric monitoring, but it is not the same as tracking a specific
  oxygen involved in a particular hydrogen bond.
- Substrate pose matters. A preserved Ser-His-Asp/Glu geometry is more
  convincing when the substrate is also positioned consistently with the
  proposed mechanism.

For mechanistic claims, combine the triad metric with direct inspection of
active-site snapshots, hydrogen-bond angle checks when available, substrate
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

## Plot behavior

The current plugin plot lifecycle reads existing artifacts and generates
high-level comparison figures. The primary outputs are:

- `triad_kde_panel.<format>` — per-pair distance distributions across conditions,
  with the configured threshold shown as a visual reference.
- `triad_threshold_bars.<format>` — grouped summaries of fractions below threshold,
  including the simultaneous contact metric and per-pair contact behavior.

The plot file format is configurable through PolyzyMD plot settings. Supported
formats include `png`, `pdf`, and `svg`; `png` is the default. For example, the
default filenames are `triad_kde_panel.png` and `triad_threshold_bars.png`, while
PDF output would use `triad_kde_panel.pdf` and `triad_threshold_bars.pdf`.

Use these plots to understand whether a scalar difference is driven by a broad
distributional shift, a small subpopulation, or one limiting pair. The plots are
interpretive aids; they do not replace statistical uncertainty or structural
validation.

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
- {doc}`../reference/analysis_triad_reference`
- {doc}`analysis_statistics_best_practices`
- {doc}`../how_to/analysis_compare_conditions`
