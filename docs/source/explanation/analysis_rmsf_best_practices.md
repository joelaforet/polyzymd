# RMSF analysis: statistical best practices

Root mean square fluctuation (RMSF) is useful for asking where a protein is
more rigid or more flexible, but it is easy to over-interpret. This page
explains what the RMSF analysis measures, how to combine residues into one
number, how residues and conditions are compared, and which claims each
quantity supports.

```{note}
**Need commands rather than interpretation guidance?** See the
[RMSF quickstart](../how_to/analysis_rmsf_quickstart.md) for copy-paste CLI
and Python examples.
```

```{seealso}
For foundational concepts such as autocorrelation, statistical inefficiency and
the replicate as the sampling unit, see
[Statistics Best Practices](analysis_statistics_best_practices.md). This page
focuses on RMSF-specific interpretation.
```

## Fluctuation, offset and deviation

After every frame is superposed on a reference structure, each atom *i* has
three per-replicate values:

$$
\text{RMSF}_i^2 = \left\langle \lvert \mathbf{r}_i(t) - \langle \mathbf{r}_i \rangle \rvert^2 \right\rangle,
\qquad
\text{offset}_i = \lvert \langle \mathbf{r}_i \rangle - \mathbf{r}_i^{\mathrm{ref}} \rvert,
\qquad
\text{deviation}_i^2 = \left\langle \lvert \mathbf{r}_i(t) - \mathbf{r}_i^{\mathrm{ref}} \rvert^2 \right\rangle
= \text{RMSF}_i^2 + \text{offset}_i^2 .
$$

The averages run over the production frames. RMSF is the fluctuation within
the sampled ensemble, the quantity `gmx rmsf -o` gives. The offset is how far
the ensemble's mean structure has moved from the reference. The deviation, the
quantity `gmx rmsf -od` gives, combines the two. Two conditions can have the
same deviation from a crystal structure for opposite reasons: one fluctuates
more about a crystal-like mean, while the other fluctuates less about a mean
that has drifted. Reporting the three together separates those cases.

Low RMSF often indicates a relatively rigid region, such as a buried core or
structured secondary element. High RMSF often indicates a flexible region,
such as a loop, terminus or mobile binding-site element. These are
interpretations, not automatic conclusions: RMSF depends on the fitted atoms,
the equilibration window, the force field and the sampling.

RMSF is related to crystallographic B-factors by

$$
B_i = \frac{8\pi^2}{3} \, \text{RMSF}_i^2
$$

(Kuzmanic and Zagrovic 2010), but crystal packing, refinement models and
simulation conditions make the correspondence approximate.

## Lower RMSF is not higher stability

Rigidity and thermal stability are not reliably linked. Karshikoff, Nilsson
and Ladenstein (2015) review experimental and simulation evidence that
"thermal tolerance of a protein is not necessarily correlated with the
suppression of internal fluctuations and mobility". On *B. subtilis* lipase A,
a network-rigidity measure correlated only fairly with the thermostability of
16 variants, R² = 0.46 (Rathi, Jaeger and Gohlke 2015). A lower RMSF in one
condition supports "this region moves less". It does not by itself support
"the protein is more stable".

Match each quantity to the claim it can support:

| Claim | Quantity |
|---|---|
| "Region X moves less" | RMSF over a core or a named region, and the per-residue RMSF difference |
| "The structure stays closer to the crystal" | Offset and deviation from an `external` reference, over the core or a region such as the active site |
| "The protein resists unfolding at this temperature" | Measures of native structure over time, such as the fraction of native contacts and secondary-structure retention, rather than RMSF (Best, Hummer and Eaton 2013) |

At high temperature, a replicate that partly unfolds mixes states. Its mean
position is then an average over those states, and its RMSF grows with when
the unfolding happened rather than measuring an equilibrium fluctuation. Look
at the per-replicate values in every report and figure before comparing
condition means.

## Combining residues into one number

Each replicate gives one value per residue. The analysis turns them into one
number per replicate in two ways:

| Result | Formula over residues $i \in S$ | What it is |
|---|---|---|
| `core_rmsf`, `<region>_rmsf` | $F = \sqrt{\tfrac{1}{\lvert S \rvert}\sum_i \overline{\text{MSF}}_i}$ | Root of the mean square fluctuation of the set, which maps onto the mean B-factor, $\langle B \rangle = 8\pi^2 F^2/3$ |
| `mean_rmsf` | $\tfrac{1}{\lvert S \rvert}\sum_i \text{RMSF}_i$ | Plain mean of the residue values |

Here $\overline{\text{MSF}}_i$ is the mean over the residue's atoms of the
squared per-atom RMSF. The same forms give `core_offset`, `core_rms_deviation`
and the region values. Only the root-mean-square form keeps the decomposition
exact for the whole set: $F_{\text{deviation}}^2 = F_{\text{RMSF}}^2 +
F_{\text{offset}}^2$. The plain mean is always smaller than or equal to $F$,
so a report should name the aggregate it uses.

Squaring gives mobile residues more weight. A residue that fluctuates by 3 Å
contributes nine times as much to $F^2$ as one at 1 Å, so a few loop or
terminal residues can dominate $F$ over all residues. The remedy is choosing
the residue set, not the formula:

- Define the core before looking at results, for example the crystal
  structure's helix and strand residues with the termini left out, with
  `--set core=...`.
- Report the loops, the lid or the active site as named regions with
  `--set regions=...`, a handful fixed in advance, which keeps the number of
  outcomes small.
- Fit on the same atoms as the core with `alignment_selection`. The frames are
  superposed by the fitted atoms, so a mobile terminus in the fit adds apparent
  motion to the core.

## Comparing conditions

Every comparison uses one value per replicate: the replicate is the sampling
unit, so the number of frames never narrows an interval. For one-value results
(the core, region and plain-mean values), each condition is compared with the
control by Welch's t test, and the Benjamini-Hochberg correction runs over the
compared conditions.

For a per-residue profile, every residue of every compared condition is one
test, and all of them form one Benjamini-Hochberg family, since together they
answer one question: where does this profile differ from the control? The text
report gives, for each condition, how many residues are significantly lower
and higher and lists them, with the family size. Every per-residue row stays
in the JSON report.

With three to five replicates, the power at each residue is low. The count of
significant residues is a lower bound on how many differ. It is not evidence
that the others do not. A large apparent effect can coexist with a
non-significant p value: the data are then not sufficient to reject the null
hypothesis at the chosen threshold, which does not prove there is no effect.

When reviewing RMSF differences, prefer cautious language:

- "The polymer condition shows lower core RMSF in these replicates" rather
  than "the polymer stabilizes the enzyme".
- "The effect is suggestive but uncertain" rather than "more replicates would
  make it significant".
- "Replicate 2 samples a different state" rather than "replicate 2 is bad",
  unless there is a documented technical failure.

## Replicates and trajectory length

LiveCoMS guidance generally favours several independent simulations over one
long trajectory when estimating uncertainty (Grossfield et al. 2018). For
RMSF, replicate-to-replicate variation shows whether a flexibility pattern is
reproducible. Several replicates:

- test reproducibility across independently initialized simulations;
- reveal outlier trajectories or rare conformational events;
- give condition-level uncertainty from replicate values;
- can run in parallel.

A single longer trajectory can still be valuable for slow processes that
shorter runs do not reach. Different seeds do not guarantee independent
equilibrium sampling if every replicate stays in the same metastable basin.

PolyzyMD computes RMSF from every production frame. Correlated frames do not
bias an average of this kind, so discarding frames only adds noise. The
statistical inefficiency and the equilibration point that pymbar detects for
each replicate are reported as information and never change the calculation
(see {doc}`convergence_detection`).

## Interpreting RMSF magnitudes

The following ranges are rough heuristics for Cα RMSF in folded proteins under
typical simulation conditions. They are not universal thresholds, and they do
not apply to offsets, deviations from an external structure, intrinsically
disordered regions, polymers or ligands.

| Approximate Cα RMSF | Common interpretation |
|---------------------|-----------------------|
| 0.3-0.5 Å | Very rigid folded core or constrained secondary structure |
| 0.5-1.0 Å | Moderate flexibility in structured regions |
| 1.0-2.0 Å | Flexible loops, flaps, or exposed regions |
| 2.0-5.0 Å | Highly mobile termini or disordered segments |
| >5.0 Å | Possible disorder, unfolding, poor alignment, or reference mismatch |

For enzyme active sites, lower RMSF is not automatically better. A rigid
active site may preserve catalytic geometry, but some enzymes need
conformational breathing, induced fit or loop motion. Active-site RMSF is most
useful together with the offset of the active site from a competent
structure, catalytic distances and experimental activity data.

## Common interpretation pitfalls

### Treating correlated frames as independent

The raw number of saved frames is not the number of independent samples. Use
the replicate values that the reports give rather than a standard error over
every frame.

### Ignoring stationarity

If RMSD or structural summaries drift after the equilibration window, RMSF
combines several regimes into one number. The question is then not only
uncertainty but whether the window represents one ensemble.

### Over-interpreting small differences

Differences of a few hundredths of an Å can be smaller than the variation
between replicates, fit choices or reference selection. Report the interval
and avoid mechanistic conclusions from tiny differences alone.

### Cherry-picking replicates

Leave out a replicate only for a documented technical reason, such as a
corrupted trajectory or a failed simulation. A conformational transition may
be scientifically important rather than an error.

### Comparing values with different definitions

Do not compare RMSF with deviation from an external structure, or `core_*`
values with `mean_*` values, as if they were the same quantity. Name the
quantity, the reference, the fitted atoms and the residue set with every
number.

## References

**Best RB, Hummer G, Eaton WA.** (2013) "Native contacts determine protein
folding mechanisms in atomistic simulations." *PNAS* 110:17874-17879.
https://doi.org/10.1073/pnas.1311599110

**Grossfield A, Patrone PN, Roe DR, Schultz AJ, Siderius DW, Zuckerman DM.**
(2018) "Best Practices for Quantification of Uncertainty and Sampling Quality
in Molecular Simulations." *Living Journal of Computational Molecular Science*
1(1):5067. https://doi.org/10.33011/livecoms.1.1.5067

**Karshikoff A, Nilsson L, Ladenstein R.** (2015) "Rigidity versus
flexibility: the dilemma of understanding protein thermal stability." *FEBS
Journal* 282:3899-3917. https://doi.org/10.1111/febs.13343

**Kuzmanic A, Zagrovic B.** (2010) "Determination of ensemble-average pairwise
root mean-square deviation from experimental B-factors." *Biophysical
Journal* 98:861-871. https://doi.org/10.1016/j.bpj.2009.11.011

**Rathi PC, Jaeger KE, Gohlke H.** (2015) "Structural rigidity and protein
thermostability in variants of lipase A from *Bacillus subtilis*." *PLoS ONE*
10:e0130289. https://doi.org/10.1371/journal.pone.0130289

## See Also

- [RMSF quickstart](../how_to/analysis_rmsf_quickstart.md) — commands and minimal setup
- [Reference structure selection](analysis_reference_selection.md) — choose the reference and fitted atoms
- [RMSF implementation verification](analysis_rmsf_verification.md) — what the implementation is checked against
- [Statistics best practices](analysis_statistics_best_practices.md) — foundational statistics for MD
