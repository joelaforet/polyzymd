# Rg: what it measures and how to read it

The radius of gyration (Rg) measures the size of a group of atoms. It is the
mass-weighted root mean square distance of the atoms from their center of
mass. A lower Rg means a more compact selection. Rg alone does not tell you
which motion changed the size.

For the commands, see {doc}`../how_to/analysis_rg_quickstart`. For replicates,
correlation and intervals, see {doc}`analysis_statistics_best_practices`.

## Definition

$$
R_g = \sqrt{\frac{1}{M} \sum_{i=1}^{N} m_i \left\| \mathbf{r}_i - \mathbf{r}_{\text{cm}} \right\|^2}
$$

- $m_i$ is the mass of atom $i$, and $M = \sum_i m_i$ is the total mass.
- $\mathbf{r}_i$ is the position of atom $i$.
- $\mathbf{r}_{\text{cm}} = \frac{1}{M}\sum_i m_i \mathbf{r}_i$ is the center
  of mass.
- $N$ is the number of atoms in the selection.

The shipped `rg` analysis calls MDAnalysis
`AtomGroup.radius_of_gyration()` on every production frame. It gives one value
in Å per frame. Each replicate's value is the mean over its frames after the
{term}`equilibration window`. The default selection is `protein`.

## Rg needs no reference and no alignment

Rg depends only on the distances of the atoms from their own center of mass.
These distances do not change when the molecule moves or rotates. Therefore
the `rg` analysis takes no reference structure and does not superpose frames.
Only the atom selection decides the result.

## Periodic boundaries

The `rg` analysis uses the coordinates as the trajectory stores them. It does
not make molecules whole. If the box boundary splits a selected molecule, the
center of mass is wrong and Rg is too large. The result is a sudden jump in
the time series.

OpenMM's DCD reporter keeps each molecule whole when it wraps molecules into
the box. A selection of several molecules, such as several polymer chains,
can still be split between periodic images. To measure such a selection, make
the molecules whole first, for example with MDAnalysis
`transformations.unwrap` in your own function (see
{doc}`../how_to/study_api`).

## The selection decides the value

Compare Rg values only when they use the same selection. Report the selection
string with every value.

| Selection | What it measures |
|---|---|
| `protein` | The size of the whole protein, side chains included |
| `protein and name CA` | The size of the backbone trace, with less side-chain noise |
| `protein and name CA and resid 20:250` | The size of a core, without flexible termini |
| `chainid C` | The size of the polymer, in the PolyzyMD chain convention |
| `protein or chainid C` | The size of the protein and polymer together |

In the PolyzyMD chain convention, chain C holds every polymer chain. See
{doc}`residue_assignment`. Use a `resname` selection to measure only some
monomer types.

## Rg depends on protein size

Rg grows with the number of residues $N$:

$$
R_g \propto N^{\nu}
$$

A compact sphere gives $\nu = 1/3$. Fits to folded proteins in the PDB give
$\nu \approx 0.38$ to $0.40$, because real proteins have cavities and rough
surfaces ([Dima and Thirumalai, 2004](https://doi.org/10.1021/jp037128y)).
Unfolded chains in good solvent give $\nu \approx 0.588$
([Kohn et al., 2004](https://doi.org/10.1073/pnas.0403643101)).

Therefore a larger protein has a larger Rg. Do not compare Rg values of two
different proteins directly. Compare each protein with its own control.

These exponents describe ideal polymers. A finite enzyme, a conjugated
polymer or a multi-domain protein does not follow them exactly. Use them to
see whether a change is large, not to classify a state.

## Reading the time series

The `rg` analysis draws the time series of each condition. Look at it before
you compare means. Two conditions can have the same mean Rg, with one stable
and one drifting.

Each pattern below has more than one possible cause. Look at the frames
around a change in a molecular viewer before you name a mechanism.

| Pattern | Causes to check |
|---|---|
| Stable | The selection keeps the same size |
| Slow rise | Unfolding, swelling, domain separation, or a polymer that moves away |
| Slow fall | Compaction, polymer wrapping, or collapse |
| Sudden jump | A conformational change, or a molecule split across the box |
| Oscillation | Exchange between compact and extended states, such as hinge motion |

For an oscillating series, report the range and the period, not only the
mean.

## Rg together with RMSD

Rg measures size. RMSD measures distance from a reference structure. A
domain rotation can change the RMSD and leave Rg unchanged. Read the two
together:

| Rg | RMSD | Causes to check |
|---|---|---|
| Stable | Stable | Same size and same structure as the reference |
| Rises | Rises | Unfolding or a large domain movement |
| Stable | Rises | A rearrangement that keeps the size, such as a domain rotation |
| Rises | Stable | Expansion of parts outside the RMSD selection |
| Falls | Rises | Compaction with a change of structure |
| Falls | Stable | Small compaction without a change of structure |

A lower protein Rg does not show that a polymer stabilizes the fold. Check
the native contacts, the secondary structure and the RMSF as well. A lower Rg
with native contacts kept means something different from a lower Rg with
helices lost.

## Reporting a difference

PolyzyMD compares conditions with one mean Rg per replicate. It reports each
condition's mean with a {term}`95 % confidence interval` and n, and a Welch t
test against the control with a {term}`Benjamini-Hochberg` adjusted p value.
Report these numbers, not the mean ± one standard error. With three
replicates, one standard error covers far less than 95 % of the likely values.

```text
WRONG: Condition A (18.26 Å) is less compact than condition B (18.29 Å).
RIGHT: A 18.26 Å (95 % CI 18.07–18.45, n = 3) and B 18.29 Å
       (95 % CI 18.13–18.45, n = 3); Welch p_adj = 0.62.
```

Rg differences of a few hundredths of an ångström are smaller than the
spread between replicates of most proteins.

## References

**Dima RI, Thirumalai D.** (2004) "Asymmetry in the shapes of folded and
denatured states of proteins." *Journal of Physical Chemistry B*
108:6564-6570. https://doi.org/10.1021/jp037128y

**Flory PJ.** (1969) *Statistical Mechanics of Chain Molecules.* Wiley
Interscience, New York.

**Kohn JE, Millett IS, Jacob J, et al.** (2004) "Random-coil behavior and the
dimensions of chemically unfolded proteins." *PNAS* 101:12491-12496.
https://doi.org/10.1073/pnas.0403643101

**Lobanov MY, Bogatyreva NS, Galzitskaya OV.** (2008) "Radius of gyration as
an indicator of protein structure compactness." *Molecular Biology*
42(4):623-628. https://doi.org/10.1134/S0026893308040195

## See also

- {doc}`../how_to/analysis_rg_quickstart` — run the analysis
- {doc}`analysis_statistics_best_practices` — replicates, intervals and tests
- {doc}`analysis_rmsd_best_practices` — distance from a reference
- {doc}`convergence_detection` — the equilibration check
