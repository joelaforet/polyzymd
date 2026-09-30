# RMSF implementation verification

This page records how the RMSF analysis computes its values, what they are
checked against, and what those checks do and do not establish.

## What the code does

`polyzymd.analyses.functions.rms_decomposition` and the single-profile
functions `rmsf` and `residue_rmsd` do the following for each replicate:

1. Superpose every production frame on the reference by the fitted atoms,
   with the same rotation as MDAnalysis `align.AlignTraj`, on a copy of the
   coordinates. The trajectory itself is never modified.
2. Load the superposed positions of the measured atoms into an in-memory
   universe and compute each atom's RMSF with MDAnalysis
   `analysis.rms.RMSF`.
3. From the same superposed positions, compute each atom's root mean square
   deviation from the reference position, as `gmx rmsf -od` does, and the
   distance of its mean position from the reference position (the offset).
4. Average each per-atom value over the residue's atoms, and also the squares
   of the per-atom values, which give the core and region values.

The production frames are passed as explicit indices, so the frames after the
equilibration window, including any boundary frames the loader leaves out of a
restart chain, are exactly the frames used.

## What the tests check

`tests/analyses/test_rmsf.py` checks, on synthetic trajectories of noisy,
rotated and translated copies of a small structure:

- `rmsf` equals MDAnalysis `rms.RMSF` after `AlignTraj` onto the same
  reference, and `residue_rmsd` and the offset equal the definitions computed
  by hand from the same aligned coordinates, to 1e-5 Å;
- the squared deviation equals the squared RMSF plus the squared offset for
  every atom, and the per-residue mean squares satisfy
  `ms_deviation = msf + ms_offset`;
- fitting on different atoms gives different values, and the trajectory's
  coordinates are unchanged afterwards;
- the core and region values are the root of the mean of `msf` over their
  residues, and `core_residue_rmsd² = core_rmsf² + core_offset²` for every
  replicate;
- a reference with the wrong number of atoms, a core that selects no residues
  and a reserved region name are refused;
- per-residue comparisons line replicates up by residue ID, correct over every
  residue of every compared condition, and keep every row in the JSON report.

On two local replicates of *B. subtilis* lipase A (Cα of residues 4 to 174,
against the 1ISP crystal), the identity held to within 3e-7 Å² per atom in all
four reference modes, and to within 3e-8 Å² for the core and region values.

## Agreement with GROMACS

`scripts/benchmarks/rmsf_gromacs_parity.py` compares the per-residue values
with GROMACS `gmx rmsf`. For each reference mode, PolyzyMD builds the reference
with `pz.reference` and computes `rmsf`, `residue_rmsd` and `offset`. The same
production frames of the same atoms, not superposed, are written to a TRR file
and the reference positions to a GROMOS96 file, and `gmx rmsf -res` fits every
frame to that reference and writes the fluctuation (`-o`) and the deviation
from the reference (`-od`). The GROMACS offset is
$\sqrt{\text{deviation}^2 - \text{RMSF}^2}$.

Run in September 2026 with GROMACS 2026.0 and MDAnalysis 2.10.0 on two
*B. subtilis* lipase A replicates at 363 K (Cα of residues 4 to 174, 171
residues, against the 1ISP crystal for `external`, production frame 50 for
`frame`):

| Replicate | Mode | Frames | `rmsf` vs `-o`: r, mean \|Δ\|, max \|Δ\| | `residue_rmsd` vs `-od`: r, mean \|Δ\|, max \|Δ\| |
|---|---|---|---|---|
| No polymer | `external` | 6,418 | 0.99999996, 2.7e-4 Å, 5.0e-4 Å | 0.99999999, 2.4e-4 Å, 5.2e-4 Å |
| No polymer | `average` | 6,418 | 0.99999996, 2.6e-4 Å, 5.1e-4 Å | 0.99999996, 2.6e-4 Å, 5.1e-4 Å |
| No polymer | `centroid` | 6,418 | 0.99999996, 2.5e-4 Å, 5.0e-4 Å | 0.99999997, 2.5e-4 Å, 5.4e-4 Å |
| No polymer | `frame` | 6,418 | 0.99999997, 2.3e-4 Å, 5.0e-4 Å | 0.99999999, 2.6e-4 Å, 5.1e-4 Å |
| SBMA-EGMA 50:50 | `external` | 716 | 0.99999952, 2.6e-4 Å, 5.0e-4 Å | 0.99999994, 2.5e-4 Å, 5.3e-4 Å |
| SBMA-EGMA 50:50 | `average` | 716 | 0.99999959, 2.4e-4 Å, 5.0e-4 Å | 0.99999959, 2.4e-4 Å, 5.0e-4 Å |
| SBMA-EGMA 50:50 | `centroid` | 716 | 0.99999951, 2.6e-4 Å, 5.0e-4 Å | 0.99999970, 2.4e-4 Å, 5.1e-4 Å |
| SBMA-EGMA 50:50 | `frame` | 716 | 0.99999956, 2.4e-4 Å, 5.0e-4 Å | 0.99999984, 2.6e-4 Å, 5.1e-4 Å |

`gmx rmsf` writes its values to 4 decimals in nm, so it rounds each value by
up to 0.00005 nm, 5e-4 Å. Every maximum difference is at that rounding, so the
values agree to the precision GROMACS writes. The offset computed from the
rounded GROMACS values agrees to within 5e-3 Å in `external`, `centroid` and
`frame` modes, the larger difference coming from taking the root of a
difference of two rounded squares. In `average` mode both programs give an
offset of almost zero (below 1.2e-4 Å), since the reference is the mean
structure itself.

With one atom per residue, GROMACS fits with equal masses, as PolyzyMD does.
For selections with several atoms of different masses per residue, GROMACS
weights its fit by mass and PolyzyMD does not, so the two are not expected to
match.

## Scope

These checks cover the tested selections, reference modes and trajectory
formats. They are not a guarantee for every topology, periodic-boundary
treatment or external-reference atom mapping. In particular, an external
reference is checked only for its atom count; whether its atoms correspond to
the trajectory's is up to the structure preparation (see
{doc}`analysis_reference_selection`).

## See Also

- [RMSF Quick Start](../how_to/analysis_rmsf_quickstart.md) — commands and minimal setup
- [RMSF Best Practices](analysis_rmsf_best_practices.md) — interpretation and pitfalls
- [Reference Structure Selection](analysis_reference_selection.md) — what each reference measures against
