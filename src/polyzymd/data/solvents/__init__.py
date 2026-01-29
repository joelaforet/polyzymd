"""Pre-parameterized solvent molecules with cached partial charges.

This subpackage contains SDF files for common solvents with pre-computed
AM1BCC partial charges. Using pre-computed charges ensures:

1. All copies of a solvent molecule have identical parameters
2. Fast parameterization (no repeated AM1BCC calculations)
3. Reproducible simulations across different runs

The solvent_molecules module provides the main interface for loading
these molecules.
"""
