# RMSF implementation verification

This page records how the RMSF analysis computes its values, what they are
checked against, and what those checks do and do not establish.

## What the code does

`polyzymd.analyses.functions.rms_decomposition` and the single-profile
functions `rmsf` and `rms_deviation` do the following for each replicate:

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
  reference, and `rms_deviation` and the offset equal the definitions computed
  by hand from the same aligned coordinates, to 1e-5 Å;
- the squared deviation equals the squared RMSF plus the squared offset for
  every atom, and the per-residue mean squares satisfy
  `ms_deviation = msf + ms_offset`;
- fitting on different atoms gives different values, and the trajectory's
  coordinates are unchanged afterwards;
- the core and region values are the root of the mean of `msf` over their
  residues, and `core_rms_deviation² = core_rmsf² + core_offset²` for every
  replicate;
- a reference with the wrong number of atoms, a core that selects no residues
  and a reserved region name are refused;
- per-residue comparisons line replicates up by residue ID, correct over every
  residue of every compared condition, and keep every row in the JSON report.

On two local replicates of *B. subtilis* lipase A (Cα of residues 4 to 174,
against the 1ISP crystal), the identity held to within 3e-7 Å² per atom in all
four reference modes, and to within 3e-8 Å² for the core and region values.

## Agreement of MDAnalysis RMSF with GROMACS

Because the fluctuation comes from MDAnalysis `rms.RMSF`, its agreement with
other programs is that of MDAnalysis. An earlier benchmark on lipase A (171 Cα
atoms, 2,500 frames of a 1,000 ns NPT run), with every method given the same
aligned coordinates, found:

| Comparison | Pearson *r* | Mean \|delta\| (Å) | Max \|delta\| (Å) |
|---|---|---|---|
| MDAnalysis `RMSF` vs GROMACS `gmx rmsf` (v2026.0), average-structure alignment | 0.99999998 | 0.000250 | 0.000579 |
| MDAnalysis `RMSF` vs GROMACS `gmx rmsf` (v2026.0), centroid-frame alignment | 0.99999998 | 0.000249 | 0.000616 |

The GROMACS difference is explained by the GRO coordinate format, which stores
positions to 0.001 nm. The benchmark script and its inputs are not tracked in
the repository.

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
