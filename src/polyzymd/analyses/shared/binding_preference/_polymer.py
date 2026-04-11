"""Polymer composition extraction from topology selections."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from ._models import PolymerComposition

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

logger = logging.getLogger(__name__)


def extract_polymer_composition(
    universe: "Universe",
    polymer_type_selections: dict[str, str] | None = None,
    polymer_chain: str = "C",
) -> PolymerComposition:
    """Extract polymer composition (residue and heavy atom counts) from trajectory.

    This function counts the number of residues and heavy atoms for each
    polymer type in the system. The composition data is used for dual
    normalization of binding preference enrichment ratios.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe with loaded topology
    polymer_type_selections : dict[str, str], optional
        Mapping of type name to MDAnalysis selection string.
        If None, auto-detects from ``chainID <polymer_chain>`` using unique
        resnames.
    polymer_chain : str
        Chain ID for auto-detection when *polymer_type_selections* is None.
        Defaults to ``"C"`` (PolyzyMD chain convention).

    Returns
    -------
    PolymerComposition
        Composition data with residue_counts and heavy_atom_counts per type.
        Heavy atoms are defined as all atoms with element != 'H'.

    Notes
    -----
    For auto-detection (polymer_type_selections=None), the function:
    1. Selects all atoms in the specified polymer chain
    2. Groups residues by resname
    3. Counts residues and heavy atoms per resname

    For explicit selections, each selection string defines which atoms
    belong to that polymer type. The function counts residues (unique
    resids) and heavy atoms within each selection.

    Examples
    --------
    >>> # Auto-detect from chain C
    >>> composition = extract_polymer_composition(universe)
    >>> print(composition.residue_counts)
    {'SBM': 50, 'EGM': 50}
    >>> print(composition.heavy_atom_counts)
    {'SBM': 750, 'EGM': 400}

    >>> # Explicit selections
    >>> selections = {"SBMA": "chainID C and resname SBM"}
    >>> composition = extract_polymer_composition(universe, selections)
    """
    residue_counts: dict[str, int] = {}
    heavy_atom_counts: dict[str, int] = {}

    if polymer_type_selections is not None:
        # Use explicit selections
        for type_name, selection in polymer_type_selections.items():
            try:
                atoms = universe.select_atoms(selection)
                if len(atoms) == 0:
                    logger.warning(
                        f"Polymer type '{type_name}' selection matched no atoms: {selection}"
                    )
                    continue

                # Count unique residues
                n_residues = len(atoms.residues)
                residue_counts[type_name] = n_residues

                # Count heavy atoms (non-hydrogen)
                # MDAnalysis stores element info - filter out hydrogens
                try:
                    heavy_atoms = atoms.select_atoms("not element H")
                    n_heavy = len(heavy_atoms)
                except (AttributeError, ValueError, TypeError):
                    # Fallback: count atoms with mass > 1.1 (heavier than H)
                    n_heavy = sum(1 for a in atoms if a.mass > 1.1)

                heavy_atom_counts[type_name] = n_heavy

                logger.debug(
                    f"Polymer type '{type_name}': {n_residues} residues, {n_heavy} heavy atoms"
                )

            except (ValueError, AttributeError, TypeError) as e:
                logger.warning(f"Failed to extract composition for '{type_name}': {e}")

    else:
        # Auto-detect from polymer chain
        try:
            polymer_atoms = universe.select_atoms(f"chainID {polymer_chain}")
            if len(polymer_atoms) == 0:
                logger.warning(f"No atoms found in chainID {polymer_chain} for polymer composition")
                return PolymerComposition()

            # Group by resname
            resnames = set(polymer_atoms.residues.resnames)

            for resname in resnames:
                # Select atoms of this resname within polymer chain
                type_atoms = universe.select_atoms(f"chainID {polymer_chain} and resname {resname}")

                # Count residues
                n_residues = len(type_atoms.residues)
                residue_counts[resname] = n_residues

                # Count heavy atoms
                try:
                    heavy_atoms = type_atoms.select_atoms("not element H")
                    n_heavy = len(heavy_atoms)
                except (AttributeError, ValueError, TypeError):
                    # Fallback: count atoms with mass > 1.1
                    n_heavy = sum(1 for a in type_atoms if a.mass > 1.1)

                heavy_atom_counts[resname] = n_heavy

                logger.debug(
                    f"Auto-detected polymer '{resname}': {n_residues} residues, "
                    f"{n_heavy} heavy atoms"
                )

        except (ValueError, AttributeError, TypeError) as e:
            logger.warning(f"Failed to extract polymer composition: {e}")
            return PolymerComposition()

    composition = PolymerComposition(
        residue_counts=residue_counts,
        heavy_atom_counts=heavy_atom_counts,
    )

    logger.info(
        f"Polymer composition extracted: {composition.total_residues} total residues, "
        f"{composition.total_heavy_atoms} total heavy atoms across {len(residue_counts)} types"
    )

    return composition


