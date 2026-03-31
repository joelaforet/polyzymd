# Interpreting Polymer Bridging Results

```{warning}
Polymer bridging is **experimental**. This page documents how to read and
reason about the outputs. Definitions and interpretation guidance may change
as the methodology matures.
```

This page explains the concepts behind polymer bridging analysis — what the
numbers mean, what they do not mean, and how to think about them in the context
of enzyme-polymer conjugate design.

## What Problem Does This Solve?

Standard contact analysis tells you *how many* polymer-protein contacts exist,
but not *how those contacts are distributed among polymer chains*. Consider two
scenarios that produce identical total contact counts:

- **Scenario A:** 10 polymer chains each contact 1 protein residue → 10
  contacts, all single-site.
- **Scenario B:** 1 polymer chain contacts 10 protein residues → 10 contacts,
  one chain bridging across the surface.

These produce the same total contact count but represent fundamentally
different binding modes. Polymer bridging analysis distinguishes them by
analyzing contacts at the **per-fragment, per-frame** level.

## The Observation Model

The analysis produces a flat table of observations, conceptually:

| Frame | Fragment | Protein residues contacted | Eligible valency | Anchor | Peripheral |
|-------|----------|---------------------------|------------------|--------|------------|
| 100   | frag_0   | {Lys42, Phe88}           | 2                | Phe88  | Lys42      |
| 100   | frag_1   | {Glu15}                  | 1                | —      | —          |
| 101   | frag_0   | {Lys42, Phe88, Asp92}    | 3                | Phe88  | Lys42, Asp92 |
| 101   | frag_1   | {}                       | 0 (no contact)   | —      | —          |

Frame 101, frag_1 produces no observation (no contact). Frame 101, frag_0
produces one observation with eligible valency 3 (high-valency). All
statistics are computed over the set of all observations.

### Why Per-Frame, Not Per-Event?

A bridging "event" in physical terms is a sustained period where a polymer
fragment bridges two or more residues. The current implementation does not
track event persistence — each frame is independent. This means a 500-frame
bridging event counts 500 times, biasing statistics toward long-lived bridges.

This is a known limitation. It means:

- **Multisite fraction** is biased toward long-lived bridging (a feature if you
  care about cumulative adhesion, a bug if you care about event frequency).
- **Anchor/peripheral probabilities** over-represent stable configurations.
- Transient, single-frame bridging events are diluted rather than suppressed.

For most conjugate design questions — *"how much of the interaction is
multisite?"* — the frame-weighted view is informative. For event-counting
questions — *"how often does bridging initiate?"* — you would need
residence-time analysis (not yet implemented).

## Understanding the C-Alpha Distance Filter

### Why Filter at All?

Without filtering, a polymer contacting residues Lys42 and Arg43 (which are
sequential neighbors, C-alpha distance ~3.8 A) is classified as multisite.
This is technically correct but not physically interesting — the polymer is
merely spanning a single local patch.

Setting `min_ca_distance_angstrom` to, e.g., 8.0 A requires that at least one
pair of contacted residues be separated by >= 8 A in **that specific frame**.
This focuses the analysis on cases where the polymer genuinely bridges across
the protein surface.

### Why Dynamic (Per-Frame)?

The filter uses frame-wise C-alpha coordinates, not static crystal-structure
distances. This matters because:

1. **Loop motions** can bring distant residues closer or push nearby residues
   apart.
2. **Domain motions** can dramatically change inter-residue distances.
3. A static threshold based on the crystal structure would be wrong when the
   protein flexes.

The downside is that the effective stringency varies with protein dynamics. In
a highly flexible region, the same pair may alternate between qualifying and
not qualifying frame-to-frame, reducing the apparent multisite fraction
compared to a static filter.

### Choosing a Threshold

| Threshold | Effect | Use When |
|-----------|--------|----------|
| `0.0` | No filtering; any 2+ residues = multisite | Exploratory; counting all multi-residue contacts |
| `5.0–8.0 A` | Moderate; filters sequential neighbors | Default starting point for most proteins |
| `10.0–15.0 A` | Stringent; requires substantial spatial separation | Looking specifically for long-range bridging |
| `> 20.0 A` | Very stringent; may produce zero events | Only for large proteins with clear domain separation |

There is no single "correct" value. Report the threshold used and consider
running the analysis at 2–3 thresholds to assess sensitivity.

## Reading the Chemistry-Aware Outputs

### What the Probabilities Are

All chemistry outputs are **descriptive frequencies** over the multivalent
observation population. "Anchor protein class probability: aromatic = 0.45"
means that 45% of multivalent observations had an aromatic residue as the
anchor.

### What They Are Not

- **Not enrichment.** A probability of 0.45 for aromatics does not mean
  aromatics are preferentially involved in bridging. If 45% of the
  surface-exposed protein residues are aromatic, 0.45 is exactly the null
  expectation. Compare with surface composition (from contacts or binding
  preference analysis) before interpreting.

- **Not mechanism.** Observing that SBM monomers frequently anchor to aromatic
  residues does not prove a causal cation-pi or hydrophobic interaction. It
  may reflect spatial proximity in the initial placement, force-field bias,
  or sampling limitations.

- **Not normalized.** If condition A has 10 multivalent observations and
  condition B has 1000, both produce probability distributions that sum to 1.
  The statistical weight of the two conditions differs enormously, but this
  is not visible in the probabilities alone. Check `contacting_observations`
  and `multisite_observations` counts.

### Reading the Heatmaps

The anchor-to-peripheral matrix answers a conditional question: *"Given that
the anchor is class X, what classes are the peripheral contacts?"*

Each row is independently normalized to sum to 1.0. This means:

- You **can** compare cells within a row (e.g., "when the anchor is aromatic,
  peripheral contacts are 40% polar and 30% charged_negative").
- You **cannot** directly compare cells across rows unless the row totals
  (number of observations per anchor class) are similar.
- Empty rows indicate that no multivalent observation had that anchor class.

The same logic applies to the polymer-anchor-to-protein-anchor matrix.

## Limitations and Scientific Caveats

### What This Analysis Cannot Tell You

1. **Whether bridging is good or bad for conjugate function.** Multisite
   attachment could stabilize the enzyme (scaffold effect) or could lock
   the protein into an unfavorable conformation. The analysis describes
   the phenomenon; functional consequences require additional evidence
   (RMSF, triad geometry, activity assays).

2. **Whether condition A has "more" bridging than condition B in an absolute
   sense.** The multisite *fraction* can increase either because multisite
   events increase or because single-site events decrease. Both shift the
   fraction. Check the raw observation counts alongside the fractions.

3. **Equilibrium properties.** MD simulations of enzyme-polymer conjugates
   are rarely at equilibrium. Bridging statistics reflect the sampled
   trajectory, not a converged ensemble. Multiple independent replicates
   with different initial conditions help, but do not guarantee convergence.

### Known Methodological Limitations

- **No residence-time decomposition.** Cannot distinguish persistent bridging
  from rapid attachment/detachment.
- **Anchor selection is heuristic.** Based on minimum atom distance, which
  may not correspond to the energetically dominant contact.
- **No excluded-volume correction.** Larger polymer fragments have more
  opportunities for multisite contact purely due to chain length, independent
  of any interaction specificity.
- **Single-structure heatmaps.** The current plotting code shows
  cross-classification heatmaps for only the first condition. Programmatic
  loading is needed for cross-condition matrix comparison.

## Worked Example: Interpreting a Result

Suppose you have two conditions:

| Metric | 100% SBMA | SBMA-EGPMA 5% |
|--------|-----------|---------------|
| Multisite fraction | 0.31 +/- 0.02 | 0.49 +/- 0.03 |
| Mean contacts / oligomer | 1.41 +/- 0.03 | 1.72 +/- 0.05 |
| High-valency fraction | 0.05 +/- 0.01 | 0.13 +/- 0.01 |
| Anchor: aromatic | 0.38 | 0.44 |
| Anchor: charged_negative | 0.22 | 0.18 |

**What you can say:**

- The mixed copolymer shows significantly more multisite attachment
  (p < 0.01, large effect size). The difference is robust across replicates.
- The valency distribution shifts toward higher valency in the copolymer.
- Aromatic residues are the most common anchor class in both conditions,
  but this may simply reflect surface composition.

**What you should not say (without further analysis):**

- "The EGPMA comonomer *causes* more bridging." (Correlation; the copolymer
  differs in multiple ways.)
- "Aromatic residues *drive* bridging via pi-stacking." (No evidence of
  mechanism from frequencies alone.)
- "The copolymer provides better enzyme stabilization through bridging."
  (Bridging is not inherently stabilizing.)

## See Also

- {doc}`../tutorials/analysis_polymer_bridging` — configuration, CLI, and output reference
- {doc}`../tutorials/analysis_binding_preference` — surface-normalized enrichment for
  comparison
- {doc}`../tutorials/analysis_statistics_best_practices` — replicate planning and
  statistical interpretation
