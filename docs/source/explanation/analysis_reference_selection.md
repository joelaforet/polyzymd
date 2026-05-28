# How RMSF reference selection changes interpretation

Root mean square fluctuation (RMSF) is not an absolute property of a protein.
It is a measurement of motion **relative to a chosen reference** after removing
overall translation and rotation by structural alignment. Changing the
reference changes the scientific question the number answers.

For new contributors, this is the most important point: two RMSF analyses can
use the same trajectory and the same atoms but produce different values because
they define a different baseline for "fluctuation." Those differences are not
implementation noise; they are part of the interpretation.

## The reference defines the question

RMSF asks how far an atom or residue is, on average, from a reference position:

$$
\text{RMSF}_i = \sqrt{\frac{1}{T} \sum_{t=1}^{T}
\left(\mathbf{r}_i(t) - \mathbf{r}_i^{\text{ref}}\right)^2}
$$

The meaning of $\mathbf{r}_i^{\text{ref}}$ depends on the reference choice. A
trajectory-derived reference asks about motion around states sampled by the
simulation. A specific-frame reference asks about deviation from one chosen
conformation. An external reference asks how the simulated ensemble departs from
an independently supplied structure.

PolyzyMD exposes these choices through settings such as `reference_mode`,
`reference_frame`, `reference_file`, `alignment_selection`, and
`centroid_selection`. Most users set them in `comparison.yaml` and run RMSF
through the CLI; this page explains why the choice matters rather than how to
configure every option.

## Centroid reference: fluctuation around a sampled central structure

The default centroid reference chooses a **real sampled frame** that is closest
to an aligned mean or cluster center. In PolyzyMD's current RMSF implementation,
centroid mode uses k-means clustering with `k=1` over the centroid selection and
then selects the sampled frame nearest that center.

This is useful because the reference is an actual conformation from the
trajectory, not a synthetic average structure. It often gives an intuitive
baseline for a stable, single-basin simulation: residues with larger RMSF values
are the residues that depart more from a representative sampled structure.

The caveat is that `k=1` should not be interpreted as "the most populated
conformational state." In a multimodal trajectory, a single cluster center can
fall between basins or be biased by transitions. The selected frame is closest
to the global center under the chosen alignment and atom selection; it may not
represent the dominant basin.

## Average reference: fluctuation around the trajectory mean

An average reference uses the mean position of each selected atom across the
aligned trajectory. This asks how much each atom fluctuates around the
trajectory mean.

For a stationary trajectory sampling one conformational basin, this can
approximate a thermal-like fluctuation measure. That interpretation becomes
weaker when the trajectory samples multiple long-lived states. In that case,
the RMSF includes both within-state motion and between-state conformational
heterogeneity. The average structure may also be geometrically unphysical
because averaged coordinates are not required to preserve realistic bond
lengths, angles, or side-chain conformations.

Average-based RMSF is therefore best understood as fluctuation around the
sampled mean, not as a guaranteed measurement of pure thermal fluctuations.

## Frame reference: deviation from a chosen simulated conformation

A frame reference uses one specified frame from the trajectory as the baseline.
This changes the interpretation from "how flexible is this region around a
central sampled structure?" to "how much does this region depart from this
particular conformation?"

That can be scientifically useful when the selected frame has independent
meaning, such as a catalytically competent active-site geometry, a ligand-bound
pose, or a conformation immediately before a transition. The resulting RMSF
values describe persistence or loss of that chosen geometry over the analyzed
ensemble.

The weakness is that a single frame may contain transient noise. A frame chosen
because it is visually appealing or because it occurs at a convenient time point
can overstate biological meaning. Contributors should describe why a selected
`reference_frame` is scientifically meaningful whenever they use this mode.

## External PDB

An external reference uses coordinates from an external structure file, commonly
a prepared crystal or model structure, as the aligned baseline. This asks how
the simulated ensemble deviates from that external conformation.

This is not the same interpretation as conventional trajectory-mean RMSF. The
value combines two effects:

1. fluctuations within the simulated ensemble; and
2. systematic offset between the simulation ensemble and the external
   structure after alignment.

That combination can be exactly what you want when the scientific question is
whether different conditions preserve a known functional structure. It is less
appropriate if the question is only local flexibility within each simulated
condition.

External-reference RMSF should also not be equated directly with experimental
B-factors. Qualitative comparisons may be informative, but they require careful
atom selection, structure preparation, comparable methodology, and explicit
attention to crystal-packing effects, refinement models, temperature,
occupancy, unresolved regions, and the limits of interpreting crystallographic
disorder as simulation fluctuation.

## Alignment selection is part of the scientific definition

Alignment removes whole-protein translation and rotation before RMSF is
computed. The atoms used for this superposition define what counts as internal
motion.

For example, aligning on all C-alpha atoms emphasizes motion relative to the
global backbone. Aligning on a stable domain can make motion in another domain
appear larger because the analysis treats the stable domain as the reference
body. This is not wrong, but it changes the question from whole-protein
flexibility to domain-relative displacement.

The same principle applies to `centroid_selection`: the atoms used to choose the
centroid frame influence which sampled structure is considered central. A
centroid chosen from all protein atoms can be influenced by side-chain or loop
motions, while a backbone-only centroid emphasizes the folded core.

## External references require structural equivalence

For an external reference to be meaningful, the selected atoms in the external
structure must correspond to the selected atoms in the trajectory. Atom count
alone is insufficient. Contributors should consider at least the following
sources of mismatch:

- atom ordering and atom-selection equivalence;
- residue mapping and residue numbering;
- chain IDs;
- missing residues or unresolved loops;
- alternate locations in experimental structures;
- protonation and tautomer states;
- residue and atom naming conventions;
- terminal patches, caps, or other end-state differences.

If these details differ, a numerically successful alignment can still compare
the wrong atoms or embed a systematic structural artifact in every RMSF value.

## Choosing a reference by scientific question

The best reference is the one that matches the claim you want the RMSF plot to
support.

| Scientific question | Reference interpretation |
|---------------------|--------------------------|
| Which residues move most around a representative sampled structure? | Use a centroid-style trajectory reference, with caution for multimodal trajectories. |
| How much does each residue fluctuate around the sampled mean? | Use an average reference, and state whether the trajectory appears stationary and single-basin. |
| How strongly does the ensemble depart from a selected functional conformation? | Use a justified frame reference. |
| Which conditions preserve or depart from an independently known structure? | Use an external reference, while treating the result as deviation from that structure rather than conventional RMSF. |

When comparing conditions, keep the reference logic consistent across the
comparison. A condition-independent external reference can make offsets from the
same structure visible. A condition-specific trajectory reference can emphasize
within-condition flexibility but may hide differences in mean structure between
conditions. Neither choice is universally better; they answer different
scientific questions.
