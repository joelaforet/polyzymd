# How the reference changes RMSF, offset and deviation

The RMSF analysis superposes every production frame on a reference structure
and then measures three things for each atom: its RMSF about its mean position,
the offset of that mean position from the reference, and its root mean square
deviation from the reference. The `reference_mode` decides two things:

1. what the frames are superposed on, by the `alignment_selection` atoms;
2. what the offset and the deviation are measured from.

RMSF is always the fluctuation about the mean position of the superposed
frames. The reference changes it only through the superposition, so RMSF from
two references that fit the same atoms is usually close. The offset and the
deviation depend on the reference directly.

## The four references

| Mode | Reference positions | Offset and deviation answer |
|---|---|---|
| `centroid` | The production frame whose `alignment_selection` atoms have the smallest RMSD to their iterative average structure (MDAnalysis `align.iterative_average`) | How far is the ensemble from its most representative sampled frame? |
| `average` | The mean positions after the frames are superposed on the first production frame, averaged, and superposed again on that average | How far is each frame from the ensemble mean? The offset is then close to zero |
| `frame` | Production frame `reference_frame`, counted from 1 after the equilibration window | How far is the ensemble from one chosen frame? |
| `external` | The atoms of `selection` and `alignment_selection` read from `reference_file` | How far is the ensemble from an independently known structure, such as a crystal? |

The `average` and `centroid` references are built separately for every
replicate from its own production frames. An `external` reference is the same
for every replicate and every condition, and its file's SHA-256 hash is
recorded with the result.

## Centroid mode

The centroid is a real sampled conformation, not a synthetic average, so its
geometry is physical. In a trajectory that samples several long-lived states,
the frame closest to the overall average can sit between them rather than in
the most populated one. Its offset and deviation then measure distance from
that in-between frame.

## Average mode

The average structure is the ensemble mean, so the offset of each atom is close
to zero and the deviation is close to the RMSF. The average structure may be
geometrically unphysical, since averaged coordinates need not keep realistic
bond lengths or side-chain conformations, but that does not affect the
fluctuation about it. For a stationary trajectory in one basin, RMSF about the
average approximates a thermal fluctuation. Across several long-lived states it
mixes within-state motion with the differences between states.

## Frame mode

Use a chosen frame when it has independent meaning, such as a catalytically
competent geometry or a ligand-bound pose, and say why it was chosen. One frame
also carries that frame's thermal noise, so its offset includes noise that a
crystal or average reference does not.

## External mode

An external structure asks whether each condition keeps a known structure. The
deviation from it combines the fluctuation within the ensemble with the
systematic offset of the ensemble's mean from that structure, and the
decomposition reports them separately. This is the mode for comparing
conditions against one common structure, since the reference is the same for
all of them.

Deviation from an external structure is not an experimental B-factor.
Qualitative comparison can be informative, but crystal packing, refinement
models, temperature, occupancy and unresolved regions all differ between a
crystal and a simulation.

## The fitted atoms are part of the definition

Superposition removes translation and rotation before anything is measured, and
the `alignment_selection` atoms define what counts as internal motion. Fitting
on all Cα atoms measures motion relative to the whole backbone. Fitting on a
stable domain makes motion in another domain look larger, because the stable
domain becomes the frame of reference. That is not wrong, but it changes the
question from whole-protein flexibility to motion relative to that domain.

When a result combines a set of residues, such as the `core` or a region, fit
on the same set to measure motion within it. A mobile terminus in the fit adds
apparent motion to every other residue.

## External references need matching atoms

For an external reference to mean anything, the selected atoms in the file
must be the same atoms as in the trajectory. PolyzyMD checks that the atom
counts match and raises an error if they do not. A matching count does not
guarantee matching atoms, so also check:

- atom ordering and atom-selection equivalence;
- residue mapping and residue numbering;
- chain IDs;
- missing residues or unresolved loops;
- alternate locations in experimental structures;
- protonation and tautomer states;
- residue and atom naming conventions;
- terminal patches, caps, or other end-state differences.

If these differ, a numerically successful superposition can compare the wrong
atoms and add a systematic offset to every value.

## Choosing a reference by question

| Question | Mode |
|---|---|
| Which residues are flexible in this condition? | Any mode; read `rmsf`. `average` or `centroid` keep the reference inside the sampled ensemble |
| How far does each condition's structure drift from a known structure? | `external`; read `offset` and `rms_deviation` |
| Does a condition keep a specific sampled geometry? | `frame`, with the frame's meaning stated |

When comparing conditions, a common `external` reference makes differences in
mean structure between conditions visible as offsets. A per-replicate
reference (`average`, `centroid` or `frame`) measures each ensemble against
itself, which shows fluctuation but hides where the mean structure went.
