# Fast mixed conjugation E2E fixture

This directory contains committed, small inputs for the opt-in fast mixed
conjugation checkpoint. The fixture builds cleaned 1UBQ with:

- An attached ABC 3-mer at `A:LYS63:NZ`
- An attached N-glycan-like GlcNAc branch at `A:ASN60:ND2`
- One free SBMA 3-mer
- TIP3P water with neutralizing ions and ff14SB plus OpenFF Sage 2.3.0

The tests write all generated artifacts to pytest temporary directories. No
runtime output, trajectories, checkpoints, or exported systems should be
committed from this fixture.

The helper rewrites runtime cache and output paths into the pytest temporary
directory before loading the YAML config so this fixture directory remains
read-only during opt-in runs.

`1ubq.pdb` was created from the RCSB 1UBQ PDB entry by retaining only `ATOM`,
`TER`, and `END` records so the test has a small, self-contained protein input.
