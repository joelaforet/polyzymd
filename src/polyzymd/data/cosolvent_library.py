"""Library of common co-solvents with physical properties.

This module provides a lookup table of common co-solvents used in
molecular dynamics simulations, including their SMILES strings,
densities, and molar masses.

The library enables users to specify co-solvents by name without
needing to provide SMILES or density values manually.

Example:
    >>> from polyzymd.data.cosolvent_library import get_cosolvent
    >>> dmso = get_cosolvent("dmso")
    >>> print(f"DMSO density: {dmso.density} g/mL")
    DMSO density: 1.1 g/mL

To extend the library, add new entries to COSOLVENT_LIBRARY dict.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple


@dataclass(frozen=True)
class CoSolventData:
    """Physical properties for a co-solvent.

    Attributes:
        name: Display name for the co-solvent.
        smiles: SMILES string for molecule generation.
        density: Liquid density in g/mL at ~25°C.
        molar_mass: Molar mass in g/mol.
        common_names: Alternative names for lookup (lowercase).
    """

    name: str
    smiles: str
    density: float  # g/mL
    molar_mass: float  # g/mol
    common_names: Tuple[str, ...] = ()


# =============================================================================
# Co-solvent Library
# =============================================================================
# Add new co-solvents here. Keys should be lowercase with no spaces.
# Densities are at approximately 25°C unless otherwise noted.
#
# To add a new co-solvent:
#   1. Find the SMILES string (e.g., from PubChem)
#   2. Find the liquid density at ~25°C
#   3. Calculate or look up the molar mass
#   4. Add common alternative names for easier lookup
# =============================================================================

COSOLVENT_LIBRARY: Dict[str, CoSolventData] = {
    # Polar aprotic solvents
    "dmso": CoSolventData(
        name="DMSO",
        smiles="CS(=O)C",
        density=1.10,
        molar_mass=78.13,
        common_names=("dimethylsulfoxide", "dimethyl sulfoxide"),
    ),
    "dmf": CoSolventData(
        name="DMF",
        smiles="CN(C)C=O",
        density=0.95,
        molar_mass=73.09,
        common_names=("dimethylformamide", "n,n-dimethylformamide"),
    ),
    "acetonitrile": CoSolventData(
        name="Acetonitrile",
        smiles="CC#N",
        density=0.786,
        molar_mass=41.05,
        common_names=("mecn", "acn"),
    ),
    # Chaotropes / Denaturants
    "urea": CoSolventData(
        name="Urea",
        smiles="C(=O)(N)N",
        density=1.32,
        molar_mass=60.06,
        common_names=("carbamide",),
    ),
    # Alcohols
    "ethanol": CoSolventData(
        name="Ethanol",
        smiles="CCO",
        density=0.789,
        molar_mass=46.07,
        common_names=("etoh", "ethyl alcohol"),
    ),
    "methanol": CoSolventData(
        name="Methanol",
        smiles="CO",
        density=0.792,
        molar_mass=32.04,
        common_names=("meoh", "methyl alcohol"),
    ),
    "isopropanol": CoSolventData(
        name="Isopropanol",
        smiles="CC(C)O",
        density=0.786,
        molar_mass=60.10,
        common_names=("ipa", "2-propanol", "isopropyl alcohol"),
    ),
    # Polyols
    "glycerol": CoSolventData(
        name="Glycerol",
        smiles="C(C(CO)O)O",
        density=1.261,
        molar_mass=92.09,
        common_names=("glycerin", "glycerine"),
    ),
    "ethylene_glycol": CoSolventData(
        name="Ethylene glycol",
        smiles="C(CO)O",
        density=1.114,
        molar_mass=62.07,
        common_names=("eg", "ethanediol"),
    ),
    # Other common co-solvents
    "acetone": CoSolventData(
        name="Acetone",
        smiles="CC(=O)C",
        density=0.784,
        molar_mass=58.08,
        common_names=("propanone",),
    ),
    "thf": CoSolventData(
        name="THF",
        smiles="C1CCOC1",
        density=0.883,
        molar_mass=72.11,
        common_names=("tetrahydrofuran",),
    ),
    "dioxane": CoSolventData(
        name="1,4-Dioxane",
        smiles="C1COCCO1",
        density=1.033,
        molar_mass=88.11,
        common_names=("dioxan",),
    ),
}


def get_cosolvent(name: str) -> Optional[CoSolventData]:
    """Look up co-solvent by name (case-insensitive).

    Searches the library by primary key first, then by common names.

    Args:
        name: Co-solvent name to look up.

    Returns:
        CoSolventData if found, None otherwise.

    Example:
        >>> get_cosolvent("DMSO")
        CoSolventData(name='DMSO', smiles='CS(=O)C', density=1.1, ...)
        >>> get_cosolvent("dimethyl sulfoxide")  # Also works
        CoSolventData(name='DMSO', smiles='CS(=O)C', density=1.1, ...)
        >>> get_cosolvent("unknown_solvent")
        None
    """
    # Normalize the input
    key = name.lower().replace(" ", "").replace("-", "").replace("_", "")

    # Direct lookup
    if key in COSOLVENT_LIBRARY:
        return COSOLVENT_LIBRARY[key]

    # Search common names
    for data in COSOLVENT_LIBRARY.values():
        normalized_common = [
            n.lower().replace(" ", "").replace("-", "").replace("_", "") for n in data.common_names
        ]
        if key in normalized_common:
            return data

    return None


def list_cosolvents() -> Dict[str, CoSolventData]:
    """Return a copy of the full co-solvent library.

    Returns:
        Dictionary of all available co-solvents.
    """
    return dict(COSOLVENT_LIBRARY)
