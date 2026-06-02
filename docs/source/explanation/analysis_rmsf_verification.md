# RMSF Implementation Verification

PolyzyMD uses a custom NumPy-based RMSF implementation rather than relying
directly on MDAnalysis's built-in `RMSF` analysis class. This page summarizes
the design rationale and benchmark evidence from representative LipA cases. It
records what those checks support, and just as importantly, what they do not
prove about every possible topology, atom selection, trajectory format, or
reference-mapping workflow.

## Why a Custom Implementation?

MDAnalysis provides a built-in
[`RMSF`](https://docs.mdanalysis.org/stable/documentation_pages/analysis/rms.html#MDAnalysis.analysis.rms.RMSF)
class, but it cannot support two features PolyzyMD requires:

1. **Autocorrelation-based frame subsampling.** PolyzyMD selects approximately
   decorrelated frames using an autocorrelation-based heuristic
   (following [Grossfield et al., 2018](https://doi.org/10.33011/livecoms.1.1.5067)),
   producing arbitrary non-uniform frame indices. MDAnalysis's `RMSF`
   class only accepts uniform `start/stop/step` slicing.

2. **External reference positions.** In `external` reference mode,
   RMSF is computed as deviations from a crystal structure's atomic
   positions rather than the trajectory's own average. MDAnalysis's
   `RMSF` class always uses the trajectory average internally.

The core computation is mathematically identical across all
implementations:

$$
\text{RMSF}_i = \sqrt{\frac{1}{N}\sum_{t=1}^{N} \left|\mathbf{r}_i(t) - \mathbf{r}_i^{\text{ref}}\right|^2}
$$

## Implementation boundary

The supported contributor interface is the public `RMSFAnalysis` plugin and the
general PolyzyMD analysis lifecycle documented in the contributor guide. The RMSF
plugin package contains private implementation files that are useful provenance
for maintainers, but they are not stable extension points and should not be
imported by contributor plugins.

Conceptually, the RMSF workflow is part of the MDAnalysis job/artifact lifecycle:
it constructs trajectory-native jobs, computes an internal RMSD time series for
autocorrelation analysis, applies the selected frame policy, and writes
canonical replicate and condition artifacts. The numerical RMSF kernel remains a
plain NumPy calculation inside that lifecycle.

## Benchmark Methodology

The benchmark evidence includes 3-way comparisons for trajectory-average
reference modes and an independent reference-kernel check for external-reference
mode:

| Method | Implementation | Notes |
|--------|---------------|-------|
| **PolyzyMD** | Custom NumPy (two-pass) | The code under test |
| **MDAnalysis** | `MDAnalysis.analysis.rms.RMSF` class | Industry-standard Python MD library |
| **GROMACS** | `gmx rmsf` (v2026.0) | Industry-standard C simulation engine |

### Test System

- **Protein:** Lipase A (LipA) from *Bacillus subtilis*, 181 residues
- **Trajectory:** 2,500 frames from a 1,000 ns NPT production run
- **Selection:** 171 C-alpha atoms (residues 5--175, excluding flexible termini)
- **Topology:** `solvated_system.pdb` (11,636 atoms)

### Modes Tested

Three RMSF reference modes were checked in the benchmarked workflow:

| Mode | Alignment Reference | RMSF Reference | Methods |
|------|-------------------|---------------|---------|
| **Average** | Trajectory average structure | Trajectory average positions | 3-way |
| **Centroid** | K-Means centroid frame (k=1) | Trajectory average positions | 3-way |
| **External** | External crystal PDB | External PDB positions | independent reference-kernel check |

GROMACS is excluded from external mode because `gmx rmsf` always computes RMSF
relative to the trajectory average — the `-s` file is used only for fitting, not
as the RMSF reference. There is no option to override this behavior. Native
MDAnalysis `RMSF` also does not validate PolyzyMD's external-reference mode,
because the MDAnalysis class uses trajectory-average positions internally. The
external-reference path was therefore checked against an independent
NumPy/reference-kernel calculation using the same aligned coordinate arrays and
external reference mapping.

### Procedure

For average and centroid modes:

1. Load the trajectory and align all frames to the reference using
   MDAnalysis `AlignTraj` (in-memory).
2. Compute per-residue C-alpha RMSF using each method on the **same
   aligned coordinates**.
3. For GROMACS: export the aligned trajectory to GRO/XTC, run
   `gmx rmsf -nofit -res`, convert output from nm to Angstroms.
4. Compare per-residue values: Pearson correlation, mean absolute
   difference, maximum absolute difference.

For external-reference mode, the same aligned coordinate arrays were compared
against an independent NumPy/reference-kernel calculation rather than against
native MDAnalysis `RMSF` or GROMACS.

## Results

### Average Mode (3-way)

| Comparison | Pearson *r* | Mean \|delta\| (Angstrom) | Max \|delta\| (Angstrom) |
|---|---|---|---|
| PolyzyMD vs MDAnalysis | 1.00000000 | 0.000000 | 0.000000 |
| PolyzyMD vs GROMACS | 0.99999998 | 0.000250 | 0.000579 |
| MDAnalysis vs GROMACS | 0.99999998 | 0.000250 | 0.000579 |

### Centroid Mode (3-way)

| Comparison | Pearson *r* | Mean \|delta\| (Angstrom) | Max \|delta\| (Angstrom) |
|---|---|---|---|
| PolyzyMD vs MDAnalysis | 1.00000000 | 0.000000 | 0.000000 |
| PolyzyMD vs GROMACS | 0.99999998 | 0.000249 | 0.000616 |
| MDAnalysis vs GROMACS | 0.99999998 | 0.000249 | 0.000616 |

### External Mode (independent check)

| Comparison | Pearson *r* | Mean \|delta\| (Angstrom) | Max \|delta\| (Angstrom) |
|---|---|---|---|
| PolyzyMD vs independent reference kernel | 1.00000000 | 0.000000 | 0.000000 |

## Interpretation

**PolyzyMD and MDAnalysis produce bit-identical results** for the benchmarked
average and centroid modes. This is expected — both run in the same Python
process on the same NumPy arrays, and the RMSF formula is deterministic once the
alignment and frame set are fixed.

**GROMACS differs by < 0.0007 Angstrom** in the average and centroid
benchmarks. This small discrepancy is fully explained by the GRO coordinate
format, which truncates positions to 3 decimal places in nanometers (0.001 nm =
0.01 Angstrom precision). The residual pattern is consistent with rounding
noise, not a systematic bias.

**Conclusion:** The custom RMSF implementation is validated for the benchmarked
LipA C-alpha selection, aligned workflow, and tested reference modes within the
tolerances shown above. These results support the design choice to keep a custom
kernel while adding features PolyzyMD requires, especially non-uniform frame
subsampling and external-reference handling. They are not a universal guarantee
for every topology, atom selection, periodic-boundary treatment, reference
mapping, or trajectory format.

## Benchmark availability

The original benchmark script and input/output assets are not currently tracked
in the repository. This page records the validation result and its scope until a
maintained, reproducible benchmark is added to the project.

## See Also

- [RMSF Quick Start](../how_to/analysis_rmsf_quickstart.md) — Get RMSF results in 5 minutes
- [RMSF Best Practices](analysis_rmsf_best_practices.md) — Statistical rigor, interpretation, pitfalls
- [Reference Structure Selection](analysis_reference_selection.md) — Average, centroid, and external reference modes
