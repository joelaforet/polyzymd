"""Pre-parameterized solvent molecules with cached partial charges.

This module provides solvent molecules with pre-computed partial charges,
ensuring that all copies of a solvent molecule in a simulation have
identical force field parameters.

Why This Matters
----------------
When using AM1BCC or other charge assignment methods, each call can produce
slightly different charges (within numerical noise). If we generate charges
independently for each solvent molecule, this leads to:

1. Non-identical molecules that should be chemically identical
2. Slow parameterization (AM1BCC is expensive)
3. Potential issues with OpenFF Interchange

By pre-computing charges once and caching them in SDF files, we ensure
reproducible, fast, and physically correct simulations.

Usage
-----
>>> from polyzymd.data.solvent_molecules import get_solvent_molecule
>>> dmso = get_solvent_molecule("dmso")
>>> # dmso already has partial charges assigned
>>> print(dmso.partial_charges)

For library solvents (dmso, ethanol, etc.), pre-computed charges are loaded
instantly from bundled SDF files.

For custom solvents, AM1BCC charges are computed once and cached in
~/.polyzymd/solvent_cache/ for future use.

Adding New Solvents
-------------------
See docs/source/tutorials/contributing.md for instructions on adding new
solvents to the library.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

from openff.toolkit import Molecule
from openff.units import unit as offunit

LOGGER = logging.getLogger(__name__)

# Path to bundled SDF files (shipped with package)
_SOLVENTS_DIR = Path(__file__).parent / "solvents"

# User cache directory for custom solvents
_USER_CACHE_DIR = Path.home() / ".polyzymd" / "solvent_cache"

# In-memory cache for loaded molecules (avoid repeated file I/O)
_loaded_molecules: Dict[str, Molecule] = {}

# =============================================================================
# Water Models - Hardcoded charges from literature
# =============================================================================

# TIP3P water charges (Jorgensen et al., 1983)
# These are well-established literature values, not computed
TIP3P_CHARGES = {
    "O": -0.834,
    "H": 0.417,
}

# SPC/E water charges (Berendsen et al., 1987)
SPCE_CHARGES = {
    "O": -0.8476,
    "H": 0.4238,
}

# TIP4P-Ew charges (Horn et al., 2004)
# Note: TIP4P models have a virtual site, but for basic use we can
# approximate with 3-site model or let the force field handle it
TIP4PEW_CHARGES = {
    "O": 0.0,  # Charge on virtual site
    "H": 0.52422,
    "M": -1.04844,  # Virtual site
}


def get_solvent_molecule(
    name: str,
    smiles: Optional[str] = None,
    residue_name: Optional[str] = None,
    cache: bool = True,
) -> Molecule:
    """Get a solvent molecule with pre-computed partial charges.

    This function provides a unified interface for obtaining solvent molecules
    with consistent partial charges. It follows this lookup order:

    1. In-memory cache (already loaded this session)
    2. Bundled library SDFs (src/polyzymd/data/solvents/)
    3. User cache (~/.polyzymd/solvent_cache/)
    4. Generate from SMILES, compute AM1BCC charges, cache for future use

    Args:
        name: Solvent identifier (e.g., "dmso", "ethanol", "tip3p").
              For water models, use "tip3p", "spce", etc.
        smiles: SMILES string. Required if solvent is not in the library.
        residue_name: 3-letter residue name for topology. If not provided,
                     uses first 3 characters of name (uppercase).
        cache: If True, cache newly generated molecules to disk.

    Returns:
        OpenFF Molecule with partial charges assigned.

    Raises:
        ValueError: If solvent not found and no SMILES provided.

    Example:
        >>> # Library solvent - instant load
        >>> dmso = get_solvent_molecule("dmso")

        >>> # Custom solvent - computed and cached on first use
        >>> custom = get_solvent_molecule(
        ...     name="my_solvent",
        ...     smiles="CCOC(=O)C",
        ...     residue_name="MYS",
        ... )
    """
    # Normalize name
    name_key = name.lower().strip()

    # Check in-memory cache first
    if name_key in _loaded_molecules:
        LOGGER.debug(f"Returning {name_key} from in-memory cache")
        return _loaded_molecules[name_key]

    # Handle water models specially (hardcoded literature charges)
    if name_key in ("tip3p", "water_tip3p", "water"):
        mol = _create_tip3p_water()
        _loaded_molecules[name_key] = mol
        return mol

    if name_key in ("spce", "spc_e", "spc/e"):
        mol = _create_spce_water()
        _loaded_molecules[name_key] = mol
        return mol

    # Set default residue name
    if residue_name is None:
        residue_name = name[:3].upper()

    # Try bundled library
    library_path = _SOLVENTS_DIR / f"{name_key}.sdf"
    if library_path.exists():
        LOGGER.info(f"Loading {name_key} from library: {library_path}")
        mol = _load_molecule_from_sdf(library_path)
        _loaded_molecules[name_key] = mol
        return mol

    # Try user cache
    _USER_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = _USER_CACHE_DIR / f"{name_key}.sdf"
    if cache_path.exists():
        LOGGER.info(f"Loading {name_key} from user cache: {cache_path}")
        mol = _load_molecule_from_sdf(cache_path)
        _loaded_molecules[name_key] = mol
        return mol

    # Need to generate - requires SMILES
    if smiles is None:
        # Try to get from cosolvent library
        from polyzymd.data.cosolvent_library import get_cosolvent

        cosolvent_data = get_cosolvent(name)
        if cosolvent_data is not None:
            smiles = cosolvent_data.smiles
            LOGGER.info(f"Found SMILES for {name_key} in cosolvent library")
        else:
            raise ValueError(
                f"Solvent '{name}' not found in library and no SMILES provided. "
                f"Please provide a SMILES string to generate the molecule."
            )

    # Generate molecule with charges
    LOGGER.info(f"Generating {name_key} from SMILES and computing AM1BCC charges...")
    mol = _generate_charged_molecule(smiles, residue_name)

    # Cache to disk if requested
    if cache:
        LOGGER.info(f"Caching {name_key} to: {cache_path}")
        _save_molecule_to_sdf(mol, cache_path)

    _loaded_molecules[name_key] = mol
    return mol


def _create_tip3p_water() -> Molecule:
    """Create a TIP3P water molecule with literature charges.

    TIP3P charges are well-established values from:
    Jorgensen, W.L., et al. J. Chem. Phys. 79, 926 (1983)

    Returns:
        OpenFF Molecule for TIP3P water.
    """
    water = Molecule.from_smiles("O")
    water.name = "water_TIP3P"

    # Assign literature charges
    charges = []
    for atom in water.atoms:
        charges.append(TIP3P_CHARGES[atom.symbol])
    water.partial_charges = charges * offunit.elementary_charge

    # Set residue metadata
    for atom in water.atoms:
        atom.metadata["residue_name"] = "HOH"

    LOGGER.debug("Created TIP3P water with literature charges")
    return water


def _create_spce_water() -> Molecule:
    """Create an SPC/E water molecule with literature charges.

    SPC/E charges are from:
    Berendsen, H.J.C., et al. J. Phys. Chem. 91, 6269 (1987)

    Returns:
        OpenFF Molecule for SPC/E water.
    """
    water = Molecule.from_smiles("O")
    water.name = "water_SPCE"

    # Assign literature charges
    charges = []
    for atom in water.atoms:
        charges.append(SPCE_CHARGES[atom.symbol])
    water.partial_charges = charges * offunit.elementary_charge

    # Set residue metadata
    for atom in water.atoms:
        atom.metadata["residue_name"] = "HOH"

    LOGGER.debug("Created SPC/E water with literature charges")
    return water


def _generate_charged_molecule(smiles: str, residue_name: str) -> Molecule:
    """Generate a molecule from SMILES and compute AM1BCC charges.

    This is used for custom solvents not in the library. The charges are
    computed once and can be cached for future use.

    Args:
        smiles: SMILES string for the molecule.
        residue_name: 3-letter residue name.

    Returns:
        OpenFF Molecule with AM1BCC partial charges.
    """
    # Create molecule from SMILES
    mol = Molecule.from_smiles(smiles)

    # Generate a 3D conformer (required for AM1BCC)
    mol.generate_conformers(n_conformers=1)

    # Compute AM1BCC charges
    # This uses the OpenFF toolkit's charge assignment
    mol.assign_partial_charges(partial_charge_method="am1bcc")

    # Set residue metadata
    for atom in mol.atoms:
        atom.metadata["residue_name"] = residue_name

    LOGGER.info(
        f"Generated molecule with {mol.n_atoms} atoms, "
        f"total charge: {sum(mol.partial_charges.magnitude):.4f}"
    )

    return mol


def _load_molecule_from_sdf(path: Path) -> Molecule:
    """Load a molecule from an SDF file with embedded charges.

    The SDF file should contain partial charges in the atom.dprop.PartialCharge
    property and residue metadata in the metadata property.

    Args:
        path: Path to the SDF file.

    Returns:
        OpenFF Molecule with charges and metadata restored.
    """
    mol = Molecule.from_file(str(path), file_format="sdf")

    # Restore metadata from SDF properties
    if "metadata" in mol.properties:
        try:
            metadata_dict = json.loads(mol.properties["metadata"])
            for idx_str, atom_meta in metadata_dict.items():
                idx = int(idx_str)
                if idx < mol.n_atoms:
                    for key, value in atom_meta.items():
                        mol.atom(idx).metadata[key] = value
        except (json.JSONDecodeError, KeyError) as e:
            LOGGER.warning(f"Could not restore metadata from {path}: {e}")

    return mol


def _save_molecule_to_sdf(mol: Molecule, path: Path) -> None:
    """Save a molecule to an SDF file with embedded charges and metadata.

    The charges are automatically saved by OpenFF Toolkit. We also save
    atom metadata (like residue_name) as a JSON property.

    Args:
        mol: Molecule to save.
        path: Output path.
    """
    # Ensure parent directory exists
    path.parent.mkdir(parents=True, exist_ok=True)

    # Package metadata into a JSON property
    metadata_dict = {}
    for i, atom in enumerate(mol.atoms):
        if atom.metadata:
            metadata_dict[str(i)] = dict(atom.metadata)

    if metadata_dict:
        mol.properties["metadata"] = json.dumps(metadata_dict)

    # Save to SDF - OpenFF Toolkit automatically includes partial charges
    mol.to_file(str(path), file_format="sdf")

    LOGGER.debug(f"Saved molecule to {path}")


def list_available_solvents() -> Dict[str, str]:
    """List all available pre-parameterized solvents.

    Returns:
        Dictionary mapping solvent names to their source (library or cache).
    """
    available = {}

    # Water models (always available)
    available["tip3p"] = "built-in"
    available["spce"] = "built-in"

    # Library solvents
    if _SOLVENTS_DIR.exists():
        for sdf_file in _SOLVENTS_DIR.glob("*.sdf"):
            available[sdf_file.stem] = "library"

    # User cache
    if _USER_CACHE_DIR.exists():
        for sdf_file in _USER_CACHE_DIR.glob("*.sdf"):
            if sdf_file.stem not in available:
                available[sdf_file.stem] = "user_cache"

    return available


def clear_cache(name: Optional[str] = None) -> None:
    """Clear cached solvent molecules.

    Args:
        name: Specific solvent to clear. If None, clears all cached solvents.
    """
    if name is not None:
        # Clear specific solvent
        name_key = name.lower().strip()
        if name_key in _loaded_molecules:
            del _loaded_molecules[name_key]

        cache_path = _USER_CACHE_DIR / f"{name_key}.sdf"
        if cache_path.exists():
            cache_path.unlink()
            LOGGER.info(f"Removed {name_key} from cache")
    else:
        # Clear all
        _loaded_molecules.clear()
        if _USER_CACHE_DIR.exists():
            for sdf_file in _USER_CACHE_DIR.glob("*.sdf"):
                sdf_file.unlink()
            LOGGER.info("Cleared all cached solvents")
