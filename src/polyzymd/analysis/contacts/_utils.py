"""Shared utility functions for contact analysis calculators.

This module contains helper functions extracted from ContactAnalyzer
and ParallelContactAnalyzer to avoid code duplication.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup, Residue


def identify_polymer_chains(
    query_residues: AtomGroup,
) -> list[list[tuple[int, str, Residue]]]:
    """Identify polymer chains and their segments from an atom group.

    Groups residues by connected fragments (molecular graph components)
    and returns structured chain information.

    Parameters
    ----------
    query_residues : AtomGroup
        MDAnalysis AtomGroup containing the polymer residues to classify.

    Returns
    -------
    list[list[tuple[int, str, Residue]]]
        List of chains, where each chain is a list of
        (chain_idx, resname, residue) tuples.
    """
    all_atoms = query_residues.atoms
    fragments = all_atoms.fragments if all_atoms.fragments else [all_atoms]

    chains = []
    for chain_idx, frag in enumerate(fragments):
        chain_residues = []
        for res in frag.residues:
            chain_residues.append((chain_idx, res.resname, res))
        if chain_residues:
            chains.append(chain_residues)

    return chains
