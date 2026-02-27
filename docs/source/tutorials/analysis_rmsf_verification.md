# RMSF Implementation Verification

PolyzyMD uses a custom NumPy-based RMSF implementation rather than
MDAnalysis's built-in `RMSF` analysis class. This page documents the
3-way benchmark that validates the custom implementation against two
industry-standard tools.

## Why a Custom Implementation?

MDAnalysis provides a built-in
[`RMSF`](https://docs.mdanalysis.org/stable/documentation_pages/analysis/rms.html#MDAnalysis.analysis.rms.RMSF)
class, but it cannot support two features PolyzyMD requires:

1. **Autocorrelation-based frame subsampling.** PolyzyMD identifies
   statistically independent frames via autocorrelation analysis
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

## Source Code

The RMSF and RMSD calculations live in a single file:

- **RMSF:** [`_compute_rmsf()`](https://github.com/joelaforet/polyzymd/blob/feature/analysis-module/src/polyzymd/analysis/rmsf/calculator.py#L633-L683)
  in `src/polyzymd/analysis/rmsf/calculator.py` (lines 633--683).
  Two-pass algorithm: (1) compute average positions from selected frames
  (or use external reference), (2) accumulate squared deviations.

- **RMSD:** [`_compute_rmsd_timeseries()`](https://github.com/joelaforet/polyzymd/blob/feature/analysis-module/src/polyzymd/analysis/rmsf/calculator.py#L613-L631)
  in the same file (lines 613--631). Used internally for autocorrelation
  analysis, not exposed to users.

Both methods are plain NumPy with no MDAnalysis analysis class
dependencies. The `MDAnalysis.analysis.rms.RMSF` and
`MDAnalysis.analysis.rms.RMSD` imports that existed previously were
confirmed dead code and removed in commit `19a247e`.

## Benchmark Methodology

A 3-way benchmark was performed comparing:

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

Three RMSF reference modes were benchmarked:

| Mode | Alignment Reference | RMSF Reference | Methods |
|------|-------------------|---------------|---------|
| **Average** | Trajectory average structure | Trajectory average positions | 3-way |
| **Centroid** | K-Means centroid frame (k=1) | Trajectory average positions | 3-way |
| **External** | External crystal PDB | External PDB positions | 2-way |

GROMACS is excluded from external mode because `gmx rmsf` always
computes RMSF relative to the trajectory average — the `-s` file is
used only for fitting, not as the RMSF reference. There is no option
to override this behavior.

### Procedure

For each mode:

1. Load the trajectory and align all frames to the reference using
   MDAnalysis `AlignTraj` (in-memory).
2. Compute per-residue C-alpha RMSF using each method on the **same
   aligned coordinates**.
3. For GROMACS: export the aligned trajectory to GRO/XTC, run
   `gmx rmsf -nofit -res`, convert output from nm to Angstroms.
4. Compare per-residue values: Pearson correlation, mean absolute
   difference, maximum absolute difference.

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

### External Mode (2-way)

| Comparison | Pearson *r* | Mean \|delta\| (Angstrom) | Max \|delta\| (Angstrom) |
|---|---|---|---|
| PolyzyMD vs MDAnalysis | 1.00000000 | 0.000000 | 0.000000 |

## Interpretation

**PolyzyMD and MDAnalysis produce bit-identical results** in all three
modes. This is expected — both run in the same Python process on the
same NumPy arrays, and the RMSF formula is deterministic.

**GROMACS differs by < 0.0007 Angstrom** (max absolute deviation across all
residues and modes). This small discrepancy is fully explained by the
GRO coordinate format, which truncates positions to 3 decimal places
in nanometers (0.001 nm = 0.01 Angstrom precision). The residual pattern
is consistent with rounding noise, not a systematic bias.

**Conclusion:** PolyzyMD's custom RMSF implementation is verified correct
against both MDAnalysis (bit-identical) and GROMACS (limited only by
coordinate format precision). The custom implementation is functionally
equivalent to industry-standard tools while supporting the additional
features (non-uniform frame subsampling, external reference positions)
that PolyzyMD requires.

## Reproducing the Benchmark

The benchmark script is located at `scripts/benchmark_rmsf_3way.py` in the
repository root. To reproduce:

```bash
mamba run -n polyzymd-env python scripts/benchmark_rmsf_3way.py
```

Requirements:
- `polyzymd-env` conda environment (MDAnalysis, NumPy, scikit-learn, matplotlib)
- GROMACS >= 2020 (`gmx rmsf`) — configure the `GMX` variable at the top of the script
- A simulation directory with topology (PDB) and trajectory (DCD) files

The script outputs per-mode CSV data and a publication-quality benchmark
figure to `scripts/benchmark_rmsf_output/`.

## See Also

- [RMSF Quick Start](analysis_rmsf_quickstart.md) — Get RMSF results in 5 minutes
- [RMSF Best Practices](analysis_rmsf_best_practices.md) — Statistical rigor, interpretation, pitfalls
- [Reference Structure Selection](analysis_reference_selection.md) — Average, centroid, and external reference modes
