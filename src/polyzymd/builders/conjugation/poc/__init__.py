"""Proof-of-concept scripts for polymer-protein conjugation.

These files validated the end-to-end conjugation workflow:

- conjugation_poc.py : Marimo notebook that builds the conjugate
  (protein loading, polymer generation, RDKit graph surgery,
  mixed-FF parameterization, OpenMM minimization)

- generate_polymer.py : Subprocess helper for polymer generation
  via polyzymd/polymerist in a Python 3.11 mbuild environment

- run_conjugate_simulation.py : Full simulation workflow
  (solvation, staged equilibration, NPT production, validation)

See GitHub issue #51 for results and the v1.4 integration plan.
"""
