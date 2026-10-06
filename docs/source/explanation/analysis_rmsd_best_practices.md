# RMSD: what it measures and how to read it

The root mean square deviation (RMSD) measures how far a set of atoms is from
a reference structure at each frame. It gives one value per frame: a time
series. RMSD is a distance from one chosen structure. It does not measure
stability, free energy or activity.

For the commands, see {doc}`../how_to/analysis_rmsd_quickstart`. For
replicates, correlation and intervals, see
{doc}`analysis_statistics_best_practices`.

## Definition

$$
\text{RMSD}(t) = \sqrt{\frac{1}{N} \sum_{i=1}^{N} \left\| \mathbf{r}_i(t) - \mathbf{r}_i^{\text{ref}} \right\|^2}
$$

- $\mathbf{r}_i(t)$ is the position of atom $i$ at frame $t$, after
  superposition.
- $\mathbf{r}_i^{\text{ref}}$ is the position of atom $i$ in the reference.
- $N$ is the number of atoms in the selection.

RMSF averages over time and gives one value per residue. RMSD averages over
atoms and gives one value per frame. For the relation between RMSD, RMSF,
offset and the per-residue deviation, see
{ref}`Fluctuation, offset and deviation <rmsf-fluctuation-offset-deviation>`.

## What the shipped analysis does

The `rmsd` analysis does these steps for each replicate:

1. It builds the reference for the `selection` atoms. `reference_mode`
   selects how. See the next section.
2. On each production frame, it calls MDAnalysis `rms.rmsd` with
   `center=True` and `superposition=True`. This moves both sets of atoms to
   their centers of geometry and rotates the frame's `selection` atoms onto
   the reference.
3. It returns the RMSD in Å, with each atom given equal weight.

Each replicate's value is the mean RMSD over its frames after the
{term}`equilibration window`. The default `selection` and
`alignment_selection` are both `protein and name CA`.

## The reference decides the question

`reference_mode` takes one of four values. If you give `reference_file` and
no mode, the mode is `external`. Otherwise the default is `centroid`.

| Mode | The reference is | The question it answers |
|---|---|---|
| `centroid` | The production frame closest to the replicate's average structure | How far does the structure move from its own typical state? |
| `average` | The mean positions of the replicate's frames | How far does each frame sit from the replicate's mean structure? |
| `frame` | Production frame `reference_frame`, counted from 1 after the equilibration window | How far does the structure move from this time point? |
| `external` | The `selection` atoms of `reference_file` | How far is the structure from a known state, such as a crystal structure? |

In `centroid` and `average` mode, PolyzyMD superposes the
`alignment_selection` atoms to find the reference. Each replicate has its own
reference. A low RMSD then means that a replicate stays near its own typical
structure, which can differ between conditions.

In `external` mode, every replicate of every condition uses the same
structure. A low RMSD means that the structure stays near that file. "Closer
to the crystal" is better only if the crystal is the state your question
needs. {doc}`analysis_reference_selection` gives each mode in detail.

## The selection decides the value

Compare RMSD values only when they use the same selection and the same
reference mode. Report both with every value.

| Selection | What it measures |
|---|---|
| `protein and name CA` | Distance of the backbone trace from the reference |
| `protein and backbone` | The same, with more atoms per residue, so a different scale |
| `protein and name CA and resid 50:150` | Distance of a core, without loops or termini |
| Active-site residues | Distance of the local geometry; a low value does not show catalytic competence |
| `chainid C and not name H*` | Distance of the polymer from its reference; it depends strongly on the reference |

Side chains move more than the backbone. An all-atom RMSD is therefore larger
than a Cα RMSD of the same frames. Use it only when side-chain rearrangement
is the quantity you want.

## RMSD depends on protein size

A larger protein gives a larger RMSD for the same kind of motion
([Sargsyan et al., 2017](https://doi.org/10.1021/acs.jctc.7b00028)). The
ranges below are rough values for the Cα RMSD of small and medium folded
proteins. They are not thresholds for "stable" or "unstable".

| Cα RMSD (Å) | Common contributors |
|---|---|
| 0.5 to 1.5 | A rigid core, a short simulation, or restraints |
| 1.5 to 2.5 | Normal backbone motion of many compact proteins |
| 2.5 to 3.5 | Flexible loops, termini, lid opening or domain motion |
| 3.5 to 5.0 | Domain rearrangement or partial unfolding |
| above 5.0 | Unfolding, or a different conformational state |

## Reading the time series

The `rmsd` analysis draws the time series of each condition. Look at it before
you compare means. Two conditions can have the same mean RMSD, with one at a
plateau and one still drifting.

Each pattern below has more than one possible cause. Look at the frames
around a change in a molecular viewer before you name a mechanism.

| Pattern | Causes to check |
|---|---|
| Rise, then plateau | Relaxation away from the starting structure, then a stable state for this selection |
| Continuous rise | Slow relaxation, domain motion, unfolding, or a reference that does not match |
| Sudden jump | A loop flip, lid opening, domain motion, or a molecule split across the box |
| Oscillation | Repeated motion such as hinge bending, or a reference between two states |

A plateau shows that this selection no longer drifts from this reference on
the time scale you simulated. It does not show that other coordinates, slow
modes or the solvent have reached equilibrium.

## RMSD as an equilibration check

PolyzyMD runs `pymbar.timeseries.detect_equilibration` on each replicate's
RMSD series. It reports where the equilibrated region starts. It warns when a
replicate with at least 20 effective samples still relaxes after your
equilibration window. It reports a replicate with fewer effective samples as
too correlated to judge. The check never changes the window or any value.
See {doc}`convergence_detection`.

The frames of a time series are correlated. PolyzyMD reports each replicate's
{term}`statistical inefficiency` g and {term}`n_eff`. Intervals and tests use
one value per replicate, so g does not change them.

## Reporting a difference

PolyzyMD reports each condition's mean RMSD with a
{term}`95 % confidence interval` and n, and a Welch t test against the control
with a {term}`Benjamini-Hochberg` adjusted p value. Report these numbers, not
the mean ± one standard error.

```text
WRONG: Condition A (1.856 Å) is less stable than condition B (1.861 Å).
RIGHT: A 1.86 Å (95 % CI 1.71–2.01, n = 3) and B 1.86 Å
       (95 % CI 1.74–1.98, n = 3); Welch p_adj = 0.91.
```

A lower RMSD means "closer to the reference". It does not mean "more stable".

## References

**Grossfield A, Patrone PN, Roe DR, Schultz AJ, Siderius DW, Zuckerman DM.**
(2018) "Best practices for quantification of uncertainty and sampling quality
in molecular simulations." *Living Journal of Computational Molecular Science*
1(1):5067. https://doi.org/10.33011/livecoms.1.1.5067

**Knapp B, Frantal S, Greshake B, Schwarz R, et al.** (2018) "Is an intuitive
convergence definition of molecular dynamics simulations solely based on the
root mean square deviation possible?" *Journal of Computational Biology*
25:1069-1077.

**Maiorov VN, Crippen GM.** (1994) "Significance of root-mean-square deviation
in comparing three-dimensional structures of globular proteins." *Journal of
Molecular Biology* 235(2):625-634. https://doi.org/10.1006/jmbi.1994.1017

**Sargsyan K, Grauffel C, Lim C.** (2017) "How molecular size impacts RMSD
applications in molecular dynamics simulations." *Journal of Chemical Theory
and Computation* 13(4):1518-1524. https://doi.org/10.1021/acs.jctc.7b00028

## See also

- {doc}`../how_to/analysis_rmsd_quickstart` — run the analysis
- {doc}`analysis_reference_selection` — choose the reference mode
- {doc}`analysis_rmsf_best_practices` — fluctuation, offset and deviation
- {doc}`convergence_detection` — the equilibration check
- {doc}`analysis_statistics_best_practices` — replicates, intervals and tests
