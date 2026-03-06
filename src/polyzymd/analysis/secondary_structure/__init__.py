"""Secondary structure analysis module.

This module provides DSSP-based secondary structure analysis for
MD trajectories using mdtraj as the backend.

Classes
-------
SecondaryStructureCalculator
    Main class for DSSP computation and aggregation across replicates.

Notes
-----
Secondary structure is assigned using ``mdtraj.compute_dssp(simplified=True)``,
which returns three states per residue per frame:

- ``H`` — helix (alpha, 3-10, and pi helices)
- ``E`` — extended strand (beta sheet)
- ``C`` — coil (everything else, including turns and bends)

Per-frame assignments are stored as integer-encoded NPZ sidecars
(0=C, 1=H, 2=E) to keep JSON results small. Persistence fractions
(fraction of frames each residue spends in each state) are stored in
the JSON result for cross-replicate aggregation.
"""

from polyzymd.analysis.secondary_structure.calculator import SecondaryStructureCalculator

__all__ = [
    "SecondaryStructureCalculator",
]
