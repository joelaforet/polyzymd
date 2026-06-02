# How RMSF reference selection changes interpretation

Root mean square fluctuation (RMSF) is not an absolute property of a protein.
In PolyzyMD's standard non-external modes, it measures motion around the
**trajectory mean positions after alignment**. The selected `reference_mode`
controls how the trajectory is aligned and how the alignment/reference structure
is generated; it does not make centroid, average, or frame modes compute direct
deviations from those fixed coordinates.

For new contributors, this is the most important point: two RMSF analyses can
use the same trajectory and the same atoms but produce different values because
alignment choices change the coordinates whose mean and fluctuations are
estimated. Only `external` mode supplies fixed RMSF reference positions and
should be interpreted as an RMSF-like deviation from an external structure.

## The reference defines the question

In standard PolyzyMD RMSF modes (`centroid`, `average`, and `frame`), RMSF asks
how far an atom or residue is, on average, from its mean position in the aligned
analyzed trajectory:

$$
\text{RMSF}_i = \sqrt{\frac{1}{T} \sum_{t=1}^{T}
\left(\mathbf{r}_i^{\text{aligned}}(t) - \left\langle \mathbf{r}_i^{\text{aligned}} \right\rangle\right)^2}
$$

Here $\left\langle \mathbf{r}_i^{\text{aligned}} \right\rangle$ is computed from
the aligned frames that enter the RMSF calculation. For non-external modes,
`reference_mode` affects how frames are superposed before this mean is computed.
It does not replace the mean with centroid-frame coordinates or a selected
trajectory frame.

External mode is different. It can use mapped external coordinates as fixed RMSF
reference positions:

$$
\text{RMSF-like deviation}_i = \sqrt{\frac{1}{T} \sum_{t=1}^{T}
\left(\mathbf{r}_i^{\text{aligned}}(t) - \mathbf{r}_i^{\text{external}}\right)^2}
$$

That external quantity includes both fluctuations within the simulated ensemble
and systematic displacement from the external structure.

PolyzyMD exposes these choices through settings such as `reference_mode`,
`reference_frame`, `reference_file`, `alignment_selection`, and
`centroid_selection`. Most users set them in `comparison.yaml` and run RMSF
through the CLI; this page explains why the choice matters rather than how to
configure every option.

## Centroid mode: sampled-frame alignment, trajectory-mean RMSF

The default centroid mode chooses a **real sampled frame** that is closest to an
aligned mean or cluster center for alignment/reference generation. In PolyzyMD's
current RMSF implementation, centroid mode uses k-means clustering with `k=1`
over the centroid selection and then selects the sampled frame nearest that
center.

This is useful because the alignment reference is an actual conformation from
the trajectory, not a synthetic average structure. After alignment, however,
standard PolyzyMD RMSF is still computed as fluctuation around the aligned
trajectory mean positions. Residues with larger RMSF values are residues with
larger fluctuations around that aligned mean, not necessarily residues with the
largest direct deviation from the centroid frame.

The caveat is that `k=1` should not be interpreted as "the most populated
conformational state." In a multimodal trajectory, a single cluster center can
fall between basins or be biased by transitions. The selected frame is closest
to the global center under the chosen alignment and atom selection; it may not
represent the dominant basin.

## Average mode: average-structure alignment, trajectory-mean RMSF

Average mode uses a trajectory-derived average structure for alignment/reference
generation. Because standard RMSF is also computed around the aligned trajectory
mean positions, this mode has the most direct interpretation as fluctuation
around the sampled mean.

For a stationary trajectory sampling one conformational basin, this can
approximate a thermal-like fluctuation measure. That interpretation becomes
weaker when the trajectory samples multiple long-lived states. In that case,
the RMSF includes both within-state motion and between-state conformational
heterogeneity. The average structure may also be geometrically unphysical
because averaged coordinates are not required to preserve realistic bond
lengths, angles, or side-chain conformations.

Average-based RMSF is therefore best understood as fluctuation around the
sampled mean, not as a guaranteed measurement of pure thermal fluctuations.

## Frame mode: selected-frame alignment, trajectory-mean RMSF

Frame mode uses one specified frame from the trajectory for alignment/reference
generation. In the current non-external RMSF path, it does **not** compute RMSF
as direct deviation from that frame's coordinates. After alignment, PolyzyMD
computes standard RMSF around the mean positions of the aligned analyzed
trajectory.

That can still be scientifically useful when the selected frame has independent
meaning, such as a catalytically competent active-site geometry, a ligand-bound
pose, or a conformation immediately before a transition, because the chosen
frame defines the alignment basis. The resulting RMSF values describe
fluctuation around the aligned trajectory mean under that alignment choice, not
persistence or loss of the selected frame geometry as a fixed coordinate target.

The weakness is that a single frame may contain transient noise. A frame chosen
because it is visually appealing or because it occurs at a convenient time point
can overstate biological meaning. Contributors should describe why a selected
`reference_frame` is scientifically meaningful whenever they use this mode.

## External PDB: fixed-reference RMSF-like deviation

An external reference uses coordinates from an external structure file, commonly
a prepared crystal or model structure, as mapped fixed RMSF reference positions.
This asks how the simulated ensemble deviates from that external conformation
after alignment.

This is not the same interpretation as standard trajectory-mean RMSF. The value
combines two effects:

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
centroid frame influence which sampled structure is used for alignment/reference
generation. A centroid chosen from all protein atoms can be influenced by
side-chain or loop motions, while a backbone-only centroid emphasizes the folded
core.

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
support. Start from the question, then choose the mode whose interpretation fits
that question.

### To ask which residues fluctuate most after representative-frame alignment

Choose `centroid`.

This is appropriate when you want alignment/reference generation based on a real
sampled conformation that represents the trajectory under the chosen centroid
selection. Interpret the RMSF values as fluctuations around the aligned
trajectory mean, not as direct deviations from the centroid frame.

Use extra caution for multimodal trajectories: the selected centroid frame may
sit near a global center rather than represent the most populated state.

### To ask how residues fluctuate around the sampled mean

Choose `average`.

This gives the most direct trajectory-mean interpretation among the
non-external modes because both the alignment/reference generation and the RMSF
calculation are tied to trajectory-derived mean structure. It is most
straightforward for stationary, single-basin trajectories.

If the trajectory samples multiple long-lived conformations, the result mixes
within-state fluctuation with between-state heterogeneity.

### To ask how fluctuations look under a meaningful trajectory-frame alignment

Choose `frame`.

This is useful when a particular trajectory frame has independent scientific
meaning, such as a catalytically competent geometry, a ligand-bound pose, or a
pre-transition conformation. The selected frame defines the alignment basis.

Do not interpret `frame` mode as direct deviation from that frame. In the
standard non-external RMSF path, PolyzyMD still reports fluctuations around the
aligned trajectory mean.

### To ask whether conditions preserve an independently known structure

Choose `external`.

This is the fixed-reference case. It is appropriate when the scientific claim is
about preservation of, or departure from, a prepared crystal structure or other
external model.

Interpret the result as an RMSF-like deviation from fixed external coordinates,
not as conventional trajectory-mean RMSF. The value combines ensemble
fluctuation with systematic offset from the external structure.

### When comparing conditions

Keep the reference logic consistent with the comparison claim.

A condition-independent external reference makes offsets from the same structure
visible across conditions. A condition-specific non-external reference emphasizes
within-condition fluctuation after alignment, but can hide differences in mean
structure between conditions.

Neither choice is universally better. They answer different scientific
questions.
