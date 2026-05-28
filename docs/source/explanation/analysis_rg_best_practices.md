# Rg Analysis: Best Practices

An interpretation guide for Radius of Gyration (Rg): what it measures, what it
does not prove by itself, and how to compare conditions without over-reading a
single compactness metric.

```{versionadded} 1.3.0
The Rg analysis plugin was added in PolyzyMD 1.3.0.
```

```{note}
**Just need quick results?** See the [Quick Start Guide](../how_to/analysis_rg_quickstart.md)
for copy-paste commands and minimal setup.
```

```{seealso}
**For foundational statistical concepts** (autocorrelation, correlation time,
the difference between means vs. variances), see the
[Statistics Best Practices Guide](analysis_statistics_best_practices.md).

This page focuses on **Rg-specific** guidance: what the values mean,
how to interpret timeseries behavior, and how to compare conditions.
```

## What is Radius of Gyration?

**Radius of Gyration (Rg)** measures the mass-weighted root mean square
distance of atoms from the center of mass of a molecular selection:

$$
R_g = \sqrt{\frac{1}{M} \sum_{i=1}^{N} m_i \left\| \mathbf{r}_i - \mathbf{r}_{\text{cm}} \right\|^2}
$$

Where:
- $m_i$ is the mass of atom $i$
- $\mathbf{r}_i$ is the position of atom $i$
- $\mathbf{r}_{\text{cm}} = \frac{1}{M}\sum_i m_i \mathbf{r}_i$ is the center of mass
- $M = \sum_i m_i$ is the total mass
- $N$ is the number of atoms in the selection

Rg is a measure of **spatial spread around the center of mass**. Lower values
are consistent with a more compact selected structure; higher values are
consistent with expansion. These interpretations are non-unique: a similar Rg
can arise from different conformations, and an Rg change does not by itself
identify the structural mechanism. Treat Rg as one compactness observable to
interpret alongside trajectory visualization and complementary analyses such as
contacts, SASA, secondary structure, RMSD, and RMSF.

```{warning}
Rg assumes the selected atoms form the molecule or fragment you intend to
measure. If a selected molecule is split across periodic boundaries, the center
of mass and distances to it can be distorted. Make molecules/fragments whole
before interpreting Rg, especially for polymers, oligomers, or wrapped protein
complexes.
```

### What Rg Measures

| Rg Behavior | Cautious Interpretation |
|-------------|-------------------------|
| Low and stable | Consistent compactness for the selected atoms |
| Gradually increasing | Possible expansion, swelling, unfolding, or domain separation |
| Gradually decreasing | Possible compaction, collapse, or reorganization |
| Plateau after change | Possible equilibration to a new size regime |
| Sudden jump up | Possible conformational transition, unwrapping, or PBC artifact |
| Sudden jump down | Possible collapse, wrapping, aggregation, or PBC artifact |
| Oscillating | Possible exchange between compact and extended states |

## Interpreting Rg Values

```{important}
Rg is highly system-specific — it depends on protein size, shape, fold topology,
and which atoms are included in the selection. There are no universal "good" or
"bad" Rg values. Always compare Rg across conditions using the **same selection
and atom types** rather than comparing against generic reference ranges.
```

### Selection Matters

The choice of atoms for Rg calculation affects the result:

| Selection | Best For |
|-----------|----------|
| `protein` | Overall protein compactness (all atoms) |
| `protein and name CA` | Backbone compactness (less noise from sidechains) |
| `protein and name CA and resid 20:250` | Core region, excluding flexible termini |
| `chainid C` | Polymer compactness/extension |
| `protein or chainid C` | Combined enzyme-polymer system size |

Under the PolyzyMD chain convention, `chainid C` is the standard polymer
selection for a single conjugated polymer chain. Use `resname`-based selections
when you intentionally want to include a multi-fragment polymer population, such
as several polymer chains or monomer chemistries, rather than because residue
names are generally more reliable than the chain convention.

### Rg Scales with Protein Size

Unlike RMSD (which is relatively size-independent for similar fold types), Rg
scales roughly as:

$$
R_g \propto N^{\nu}
$$

for proteins, where $N$ is the number of residues and $\nu$ is the Flory
exponent. For compact globular proteins, the theoretical expectation is
$\nu = 1/3$ (solid sphere packing), though empirical fits to PDB structures
give $\nu \approx 0.38$–$0.40$ due to imperfect packing, voids, and surface
roughness ([Dima & Thirumalai, 2004](https://doi.org/10.1021/jp037128y)).
This means larger proteins have inherently larger Rg values. When comparing
Rg across different proteins, normalize by the expected Rg for the protein
size.

## Rg vs Time: What to Look For

### Stable Plateau (Ideal)

```text
Rg (Å)
  20 |  ___________________________
     | /
  19 |/
     |
  18 |
     |
  17 +-----------------------------→ Time
     0    10   20   30   40   50 ns
```

A stable Rg plateau suggests the selected atoms maintain consistent compactness
throughout the analyzed interval. The initial transient (if present) may be
equilibration.

### Expansion (Possible Unfolding)

```text
Rg (Å)
  25 |                          /
     |                         /
  23 |                        /
     |                       /
  21 |                ______/
     |_______________/
  19 +-----------------------------→ Time
```

A rising Rg trend suggests the selected atoms are expanding. This may reflect
unfolding, swelling, domain separation, ligand/polymer unwrapping, or a PBC
artifact if the selection is not whole. Possible responses:

- Inspect frames in a molecular viewer before assigning a mechanism.
- Compare with RMSD/RMSF, secondary structure, contacts, and SASA when relevant.
- Verify force-field, topology, and trajectory imaging choices.

### Compaction

```text
Rg (Å)
  20 |_______________
     |               \
  19 |                \
     |                 \_______
  18 |
     |
  17 +-----------------------------→ Time
```

A decreasing Rg trend indicates that the selected atoms are becoming more
compact. This can be consistent with tighter packing, polymer wrapping,
hydrophobic collapse, or loss of extended structural elements. It is not, by
itself, evidence of stabilization; check whether native contacts, secondary
structure, RMSD/RMSF, or active-site geometry support that interpretation.

### Sudden Jumps

A sharp Rg change mid-trajectory typically indicates a conformational
transition. Check the structure at the jump to understand what happened:

- Domain rearrangement or hinge motion
- Partial unfolding or refolding
- Ligand unbinding leading to structural change
- Polymer wrapping or unwrapping
- Molecule or fragment wrapping across periodic boundaries

```{tip}
When you observe a jump, load the trajectory in a molecular viewer (e.g.,
VMD, PyMOL) and examine frames around the transition. Compare with the RMSD
timeseries — a jump in Rg should correlate with RMSD changes if the same
region is affected.
```

### Oscillations

Regular Rg oscillations suggest the system samples between compact and
extended conformational states. This is often seen with:

- Breathing motions in multi-domain enzymes
- Allosteric transitions
- Polymer wrapping/unwrapping cycles

For oscillating systems, report the **range and period** of oscillation rather
than just the mean Rg.

## How PolyzyMD Handles Autocorrelation

Rg timeseries are correlated — adjacent frames are similar because MD
evolves continuously. PolyzyMD automatically accounts for this:

1. **Computes Rg timeseries** using MDAnalysis `AtomGroup.radius_of_gyration()`
2. **Estimates correlation time (τ)** via autocorrelation function integration
3. **Computes effective sample size** — `n_independent = n_frames / (2τ)`
4. **Reports autocorrelation-corrected SEM** — `SEM = σ / √n_independent`

### Example Autocorrelation Output

```text
Run: Whole Protein
  Correlation time: 3821 ps (3.8 ns)
  Statistical inefficiency: 473.7
  Independent samples: 19 (from 9000 frames)
  SEM (corrected): 0.098 Å
```

This means:
- Rg values decorrelate over ~3.8 ns timescales
- You effectively have 19 independent measurements from 9000 frames
- The reported SEM properly accounts for this correlation

```{seealso}
For the mathematical details of autocorrelation functions and the LiveCoMS
recommendations, see the
[Statistics Best Practices Guide](analysis_statistics_best_practices.md).
```

## Multi-Run Analysis: Why and When

### Why Multiple Runs?

Different Rg selections answer different questions:

| Run Label | Selection | Question |
|-----------|-----------|----------|
| "Whole Protein" | `protein` | Overall protein compactness? |
| "Protein Backbone" | `protein and name CA` | Backbone compactness (less side-chain noise)? |
| "Core Region" | Core residues only | Is the structured core stable? |
| "Polymer" | `chainid C` | Is the polymer extended or collapsed? |
| "Enzyme+Polymer" | `protein or chainid C` | Overall conjugate compactness? |

### When to Use Multi-Run

- **Always** include at least one whole-protein or backbone Rg run as a baseline
- **Add core-region runs** when flexible termini or loops dominate the signal
- **Add polymer runs** when studying enzyme-polymer conjugate behavior
- **Add combined runs** when the relative sizes of enzyme and polymer matter

### Independent Ranking

Each run is ranked independently across conditions. This prevents averaging
Rg from structurally different selections (which would be meaningless):

```text
Rankings:
  Whole Protein:     With Polymer < No Polymer (compaction)
  Protein Backbone:  With Polymer < No Polymer (compaction)
  Polymer:           100% SBMA < 100% EGMA (more compact polymer)
```

## Why No Alignment or Reference?

Rg is **intrinsically translation and rotation invariant**. The quantity being
measured — the mass-weighted spread of atoms around their center of mass —
does not depend on the absolute position or orientation of the molecule in
the simulation box.

Mathematically, this is because:
- The center of mass $\mathbf{r}_{\text{cm}}$ moves with the molecule
- The distances $\|\mathbf{r}_i - \mathbf{r}_{\text{cm}}\|$ are internal
  coordinates

This gives Rg several practical advantages over RMSD:
- **No alignment artifacts** — RMSD can be affected by imperfect alignment
- **No reference structure needed** — no need to choose centroid, average, or external
- **Simpler configuration** — only `label` and `selection` are required
- **Complementary information** — Rg and RMSD together give a more complete picture

```{note}
This does not mean Rg is "better" than RMSD — they measure different things.
RMSD tells you *how much* the structure has changed from a specific reference.
Rg tells you *how compact* the structure is, regardless of what it looked
like before. Use both for comprehensive structural analysis.
```

## Replicates vs Longer Simulations

### The LiveCoMS Recommendation

> "Multiple independent simulations are preferable to a single long simulation"
> — Grossfield et al. (2018)

### Why Replicates Matter for Rg

| Multiple Replicates | Single Long Simulation |
|--------------------|------------------------|
| Truly independent starting points | Frames remain correlated |
| Tests reproducibility of compactness | May be trapped in metastable state |
| Robust SEM from replicate means | SEM requires autocorrelation correction |
| Parallelizable | Sequential |

### How Many Replicates?

| Replicates | Statistical Power | Practical Guidance |
|------------|-----------------|-------------------|
| 1 | Descriptive only | Exploratory — no inferential statistics |
| 3 | Large effects (d > 2) | Minimum for publication |
| 5 | Medium effects (d > 1.3) | Recommended standard |

```{note}
With only 1 replicate, PolyzyMD still computes Rg and includes the condition in
descriptive summaries and rankings. Replicate SEM is unavailable and shown as
n/a because variability across independent simulations cannot be estimated from
a singleton. Pairwise inferential tests require at least 2 replicates per
condition.
```

## Comparing Conditions

### What PolyzyMD Computes

For each Rg run, the comparison produces:

| Statistic | Description |
|-----------|-------------|
| **Ranking** | Conditions sorted by mean Rg (lowest = most compact) |
| **Percent change** | Relative to control condition |
| **Direction** | `compaction` (< −1%), `expansion` (> +1%), or `unchanged` |
| **t-statistic** | Two-sample t-test on replicate means |
| **p-value** | Two-tailed significance |
| **Cohen's d** | Effect size magnitude |
| **ANOVA** | Omnibus F-test when 3+ conditions (per-run) |

### Direction Labels

PolyzyMD classifies the direction of change based on percent change in mean
Rg relative to control:

| Percent Change | Direction | Meaning |
|---------------|-----------|---------|
| < −1% | `compaction` | Mean Rg is lower than control for this selection |
| > +1% | `expansion` | Mean Rg is higher than control for this selection |
| −1% to +1% | `unchanged` | Mean Rg change is within the direction threshold |

### Interpreting the Comparison

```text
Rg Comparison — Whole Protein
================================
Ranking (lower = more compact):
  1. 100% SBMA:   17.812 ± 0.038 Å
  2. No Polymer:  18.256 ± 0.044 Å
  3. 100% EGMA:   18.891 ± 0.061 Å

100% SBMA vs No Polymer:
  Change: -2.4% (compaction), p=0.0123*, d=1.87 (large)

100% EGMA vs No Polymer:
  Change: +3.5% (expansion), p=0.0078*, d=2.14 (large)

ANOVA: F=22.31, p=0.0018* (significant across all conditions)
```

**Reading this output:**
- SBMA is associated with lower whole-protein Rg for this selection.
- EGMA is associated with higher whole-protein Rg for this selection.
- The ANOVA confirms at least one condition differs from the others
- Large Cohen's d values mean the replicate-level Rg differences are large
  relative to replicate variation.

Lower or higher Rg alone does not establish stabilization or destabilization of
the native fold. Use the comparison as evidence for a compactness difference,
then check structural mechanisms with visualization and complementary analyses.

## Common Pitfalls

### 1. Insufficient Equilibration

**Symptom:** Rg mean and comparison results change with different
`--eq-time` values.

**Cause:** Including the initial equilibration phase biases the mean.

**Solution:** Plot the Rg timeseries (or RMSD timeseries) and visually
identify when the plateau begins. Set `--eq-time` to skip the transient:

```bash
polyzymd compare run rg -f comparison.yaml --eq-time 20ns
```

### 2. Comparing Different Selections

**Symptom:** Rg values are not comparable across runs or publications.

**Cause:** Different atom selections yield different Rg magnitudes.

**Solution:** Always report the exact selection string. Compare only runs with
identical selections.

### 3. Over-Interpreting Small Differences

**Symptom:** Claiming significance for 0.05 Å Rg differences.

**Cause:** Not accounting for uncertainty.

**Solution:** Always report uncertainty and check statistical significance:

```text
# WRONG: "Condition A (18.256 Å) is less compact than B (18.291 Å)"
# RIGHT: "Condition A (18.256 ± 0.044 Å) and B (18.291 ± 0.038 Å)
#         are not significantly different (p = 0.62, unchanged)"
```

### 4. Ignoring Timeseries Shape

**Symptom:** Reporting only mean Rg without inspecting the timeseries.

**Cause:** Two conditions can have the same mean Rg but very different
dynamics (one stable, one drifting upward then returning).

**Solution:** Always examine the Rg timeseries plots. Use
`polyzymd compare plot-all` to generate them automatically.

### 5. Confusing Rg with RMSD

**Symptom:** Expecting RMSD-like values (1–3 Å) from Rg analysis.

**Cause:** Rg and RMSD measure fundamentally different quantities. Rg values
are typically much larger (12–25 Å for whole proteins).

**Solution:** Understand that Rg is an absolute size measure, while RMSD is a
relative deviation measure. A 1 Å change in Rg is typically a smaller
*relative* change than a 1 Å change in RMSD.

### 6. Ignoring Replicate Variation

**Symptom:** Reporting within-trajectory SEM as the total uncertainty.

**Cause:** Treating autocorrelation-corrected SEM as sufficient.

**Solution:** Use replicate-based statistics when available. The replicate
SEM captures system-level variability that within-trajectory analysis cannot:

```text
Replicate 1: mean Rg = 18.234 Å
Replicate 2: mean Rg = 18.291 Å
Replicate 3: mean Rg = 18.244 Å

Replicate mean: 18.256 Å
Replicate SEM:  0.044 Å  ← This is the gold standard uncertainty
```

### 7. Using Inappropriate Selections

**Symptom:** Rg timeseries is noisy or dominated by flexible regions.

**Cause:** Including highly flexible termini, disordered loops, or solvent
atoms in the selection.

**Solution:** Match your selection to your scientific question:

- Whole-protein Rg → `"protein"` or `"protein and name CA"`
- Core stability → exclude flexible termini with specific residue ranges
- Polymer behavior → `"chainid C"`

## Fragment Mode Best Practices

```{versionadded} 1.3.0
```

When your selection contains multiple disconnected molecules (e.g., many
polymer chains in solution), use `calculation_mode: "fragments"` to compute
per-fragment Rg and reduce to a meaningful per-frame average. Without
fragment mode, the whole-group Rg is dominated by the spatial separation
between molecules rather than individual chain conformations.

### Selection Strategy

For a standard PolyzyMD enzyme-polymer conjugate, `chainid C` remains the
canonical polymer-chain selection. Use `resname`-based selections when your
scientific question is explicitly about a multi-fragment polymer population,
for example several independent chains or multiple polymer residue names in the
same condition.

```{code-block} yaml
:caption: Fragment-mode run nested under plugins.rg.runs

plugins:
  rg:
    runs:
      - label: polymer_population_rg
        selection: "resname SBM or resname EGM or resname EGP"
        calculation_mode: fragments
```

### Verify Fragment Count

Before running large production analyses, verify that MDAnalysis detects the
expected number of fragments with a quick test:

```python
import MDAnalysis as mda

u = mda.Universe("topology.pdb", "trajectory.dcd")
ag = u.select_atoms("resname SBM or resname EGM or resname EGP")
print(f"Atoms: {len(ag)}, Fragments: {len(ag.fragments)}")
```

If the fragment count does not match the expected number of independent
polymer chains, check your topology for unexpected bonds bridging chains.

### When to Use Each Mode

| Scenario | Recommended mode |
|----------|-----------------|
| Single protein chain | `selection` (default) |
| Single polymer chain | `selection` |
| Many polymer chains in solution | `fragments` |
| Oligomer populations | `fragments` |
| Protein + single polymer combined | `selection` |

### Fragment Weighting

- **`equal`** (default): Arithmetic mean — all fragments contribute equally
  regardless of size. Best when fragments are similar in length and you want
  each chain to contribute equally to the per-frame reduced value.
- **`mass`**: Mass-weighted mean — heavier fragments contribute more. Best
  when fragment sizes vary significantly and you want the average to reflect
  the total material, not just the chain count.

### Statistical Comparison with Fragment Mode

The **reduced Rg timeseries** (per-frame mean across fragments) is the
primary metric used for cross-condition statistical comparison (t-tests,
ANOVA, ranking). This is stored in `rg_values` in the NPZ sidecar and
drives the mean, SEM, and correlation time reported in JSON results.

The **fragment Rg distribution** is supplementary and descriptive. Fragment
values pooled across frames and trajectories are correlated and should not be
treated as independent evidence for statistical significance. Use fragment
distributions to understand *why* replicate-level comparisons may differ, not
as a substitute for reduced or replicate-level inference.

## Interpreting Distribution Plots

Distribution plots provide a deeper view of Rg behavior beyond mean and SEM.

### Reduced Rg Distribution

The reduced distribution shows the spread of per-frame Rg values (one value
per frame). Because each frame's value is already an average over multiple
fragments (in fragment mode), this distribution is relatively **narrow** —
a consequence of the central limit theorem.

Use reduced distributions to:

- Compare overall conformational states across conditions
- Identify bimodal behavior (two distinct conformational states)
- Assess whether conditions produce overlapping or distinct Rg ranges

### Fragment Rg Distribution

The fragment distribution pools ALL individual fragment Rg values across
all frames and all replicates. It captures the **full range of sizes**
that individual chains adopt, including rare extended or collapsed
conformations that average out in the reduced series.

Use fragment distributions to:

- Detect conformational heterogeneity within a population
- Identify subpopulations of chains with distinct sizes
- Understand the physical origin of differences seen in reduced distributions

Because these values are correlated within trajectories and across fragments in
the same simulation box, interpret this distribution as descriptive structure,
not as a collection of independent samples for p-values.

### Comparing Reduced and Fragment Distributions

| Observation | Possible Interpretation |
|-------------|-------------------------|
| Reduced distributions differ, fragment distributions also differ | Many chains may shift conformational state in the same direction |
| Reduced distributions differ, fragment distributions overlap | The reduced mean may be influenced by a subset of chains or frames |
| Reduced distributions overlap, fragment distributions differ | Individual chains may sample different states that average out |
| Both distributions overlap | No clear descriptive distribution difference |

```{tip}
If reduced distributions overlap but fragment distributions differ, this
suggests individual chains are sampling different conformational states that
cancel out in the average. This is a sign of **conformational heterogeneity**
that merits visual inspection of trajectories.
```

## Rg as a Folding Clue

Rg is a classical compactness measure used in folding studies. For idealized
polymer ensembles, the relationship between Rg and chain length follows rough
scaling expectations
(Flory, *Principles of Polymer Chemistry*, Cornell University Press, 1953;
de Gennes, *Scaling Concepts in Polymer Physics*, Cornell University Press, 1979;
[Kohn et al., 2004](https://doi.org/10.1073/pnas.0403643101)):

- **Folded globular expectation:** $R_g \propto N^{1/3}$ for ideal compact
  packing. Empirical fits to protein structures are often closer to
  $\nu \approx 0.38$–$0.40$ because real proteins have imperfect packing,
  cavities, and rough surfaces ([Dima & Thirumalai, 2004](https://doi.org/10.1021/jp037128y)).
- **Self-avoiding coil expectation:** $R_g \propto N^{0.588}$, often
  approximated as a Flory exponent of $3/5$ for unfolded chains in good solvent
  ([Kohn et al., 2004](https://doi.org/10.1073/pnas.0403643101)).
- **Fully extended geometric limit:** $R_g \propto N^{1.0}$ for a stretched
  chain-like object.

These are polymer-physics expectations, not universal diagnostics for finite
enzyme or enzyme-polymer MD trajectories. Finite size, domain architecture,
glycosylation or conjugation, solvent quality, force-field behavior, and the
chosen atom selection can all blur these categories.

Monitoring Rg during simulation can provide clues about:

- **Possible unfolding or swelling:** Rg increases away from a compact baseline.
- **Possible refolding or compaction:** Rg decreases from an extended baseline.
- **Possible molten-globule-like behavior:** Rg is moderately elevated and
  fluctuating, but secondary-structure and contact analyses are needed before
  assigning this state.

```{tip}
For enzyme-polymer conjugate studies, comparing protein Rg across conditions
can suggest whether a polymer is associated with compaction or expansion. It
does not by itself prove stabilization or destabilization of the native fold;
combine it with native contacts, secondary structure, RMSD/RMSF, SASA, and
visual inspection.
```

## Complementary Use with RMSD

Rg and RMSD provide complementary structural information. Using both together
gives a more complete picture:

| Rg Trend | RMSD Trend | Possible Explanation |
|----------|-----------|----------------------|
| Stable | Stable | Similar size and similar reference-relative structure |
| Increasing | Increasing | Possible unfolding, domain movement, or major conformational change |
| Stable | Increasing | Local rearrangement or domain rotation without overall size change |
| Increasing | Stable | Expansion that preserves the aligned reference-relative features |
| Decreasing | Increasing | Compaction coupled to structural reorganization |
| Decreasing | Stable | Mild compaction without large reference-relative deviation |

```{important}
When Rg and RMSD disagree, investigate further. For example, a stable Rg
with increasing RMSD could mean a domain rotation that changes local structure
without changing overall compactness. A molecular viewer is essential for
interpreting such cases.
```

Rg is also complementary to contacts, SASA, secondary structure, and RMSF. For
example, a lower protein Rg with preserved native contacts and secondary
structure supports a different interpretation than a lower Rg accompanied by
lost helices or collapsed non-native contacts.

## References

### Primary References

**Flory PJ.** (1969) *Statistical Mechanics of Chain Molecules.* Wiley
Interscience, New York.

Foundational work establishing the theoretical framework for polymer chain
dimensions, including Rg scaling laws.

**Lobanov MY, Bogatyreva NS, Galzitskaya OV.** (2008) "Radius of Gyration
as an Indicator of Protein Structure Compactness." *Molecular Biology*
42(4):623-628. https://doi.org/10.1134/S0026893308040195

Systematic analysis of Rg as a compactness metric for protein structures,
including empirical scaling relationships.

### Additional References

**Grossfield A, Patrone PN, Roe DR, Schultz AJ, Siderius DW, Zuckerman DM.**
(2018) "Best Practices for Quantification of Uncertainty and Sampling Quality
in Molecular Simulations." *Living Journal of Computational Molecular Science*
1(1):5067. https://doi.org/10.33011/livecoms.1.1.5067

Best practices for handling autocorrelation and uncertainty in MD observables
including Rg timeseries.

**Vitalis A, Pappu RV.** (2009) "Methods for Monte Carlo Simulations of
Biomacromolecules." *Annual Reports in Computational Chemistry* 5:49-76.

Discussion of Rg as an order parameter for conformational sampling quality.

## See Also

- [Quick Start Guide](../how_to/analysis_rg_quickstart.md) — Get results fast
- [Statistics Best Practices](analysis_statistics_best_practices.md) — Foundational statistics for MD
- [RMSD Best Practices](analysis_rmsd_best_practices.md) — Complementary structural deviation analysis
- [RMSF Best Practices](analysis_rmsf_best_practices.md) — Per-residue fluctuation analysis
- [Compare Simulation Conditions](../how_to/analysis_compare_conditions.md) — Full comparison workflow
- [LiveCoMS Best Practices](https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5067) — Full methodology paper
