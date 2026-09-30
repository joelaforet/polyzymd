# Contacts implementation verification

This page records how `polyzymd analyze contacts` counts a contact under each
method, and what its values were checked against.

## What the code does

`polyzymd.analyses.functions.residue_occlusion`, behind `method=occlusion`,
does the following for each frame:

1. Compute each protein residue's SASA with the protein alone, as
   `residue_sasa` does: MDTraj `shrake_rupley(mode="atom")` on the protein
   atoms, one frame per call, summed per residue.
2. Move each polymer molecule, a bonded fragment of the polymer selection, by
   the box vector that brings its centroid to the periodic image nearest the
   protein centroid (MDAnalysis `lib.distances.minimize_vectors`).
3. Compute each protein residue's SASA again with the protein and the moved
   polymer atoms, which cover the protein without being counted.
4. Divide both by the residue's maximum ASA from Tien et al. (2013). The
   residue is exposed when the first ratio is at least `threshold`, and in
   contact when it is exposed and the second ratio is below `threshold`. The
   occluded area is `max(0, alone - with)`.

It then averages over the production frames. `residue_contacts`, behind
`method=distance`, finds atom pairs within the cutoff with MDAnalysis
`lib.distances.capped_distance` on each frame.

## Agreement with an independent occlusion script

The analysis scripts of a study of polymer-coated lipases compute the same
occlusion with their own code, which calls MDTraj directly and moves each
intact oligomer to its nearest image the same way. In September 2026 the
per-frame, per-residue SASA values those scripts stored for one 363 K
*B. subtilis* lipase A replicate with a 50:50 SBMA-EGMA polymer (179 residues,
2,697 protein atoms, 7,722 polymer atoms) were compared with
`residue_occlusion` on the same frames:

| Check | Result |
|---|---|
| The script's first frame of a 10-frame MDTraj call, protein alone and with polymer | Every residue agrees within 3 × 10⁻⁵ Å² |
| Its later frames of that call, against `residue_occlusion`, one frame per call | Every residue about 0.1 percent larger in the script |
| Its later frames, against the same 10 frames given to MDTraj in one call on one thread | Every residue of all 10 frames agrees within 3 × 10⁻⁵ Å² |
| Contact on each residue and frame, 10 frames | 1,790 of 1,790 residue-frames agree |

The script passed 10 frames to each MDTraj call, and MDTraj 1.11.1 returns
about 0.1 percent too much area for the later frames of a call (see
{doc}`analysis_sasa_verification`). Reproducing its batching gives its values
exactly, so the two implementations compute the same quantity, and
`residue_occlusion` gives the values without that bias.

## Agreement with the distance count used before this version

Before this version contacts were counted by distance, on all atoms within
4.5 Å. On 716 production frames of the same replicate, `--set method=distance
--set cutoff=4.5 --set heavy_atoms=false` gave a coverage of 0.9217877094972067
and a mean contact fraction of 0.5482974938360227, and every residue's contact
fraction and per-monomer-type fraction came out exactly equal to that count.

## What the tests check

`tests/analyses/test_residue_contacts.py` checks, on small synthetic systems,
the contact fractions against counts made by hand for both methods, the exact
cutoff, the minimum image, the nearest-image move of whole polymer molecules,
the exposure condition of occlusion, residues without a maximum ASA, the
per-type rows, and the results, figures, settings and refusals of
`polyzymd analyze contacts`.

## Scope

These checks cover the tested selections and trajectories. The occlusion
values depend on MDTraj's radius table and on Tien et al.'s maximum ASA
values, which exist only for the 20 standard amino acids and their
protonation states. The nearest-image move assumes each polymer molecule is
whole in the trajectory; the protein is used as loaded.

## See Also

- [Contacts Quick Start](../how_to/analysis_contacts_quickstart.md) — commands and settings
- [SASA implementation verification](analysis_sasa_verification.md) — the MDTraj batching defect
- [Shipped analysis functions](../reference/analysis_functions.md) — what each function measures
