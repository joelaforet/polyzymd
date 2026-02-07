"""
Default reaction templates for polymer generation.

This module provides access to bundled ATRP (Atom-Transfer Radical Polymerization)
reaction templates for methacrylate-based monomers.

The reaction templates define:
- Initiation: Chlorination of the vinyl group
- Polymerization: Chain extension via radical coupling
- Termination: Restoration of terminal alkene

Example:
    >>> from polyzymd.data.reactions import get_atrp_reactions
    >>> reactions = get_atrp_reactions()
    >>> print(reactions["initiation"])
    PosixPath('/path/to/atrp_initiation.rxn')

Made by PolyzyMD, by Joseph R. Laforet Jr.
"""

from pathlib import Path
from typing import Dict

# Directory containing reaction files
_REACTIONS_DIR = Path(__file__).parent


def get_atrp_reaction_paths() -> Dict[str, Path]:
    """Get paths to the default ATRP reaction template files.

    Returns:
        Dictionary with keys 'initiation', 'polymerization', 'termination'
        mapping to the corresponding .rxn file paths.

    Raises:
        FileNotFoundError: If any reaction file is missing.
    """
    paths = {
        "initiation": _REACTIONS_DIR / "atrp_initiation.rxn",
        "polymerization": _REACTIONS_DIR / "atrp_polymerization.rxn",
        "termination": _REACTIONS_DIR / "atrp_termination.rxn",
    }

    # Validate all files exist
    for name, path in paths.items():
        if not path.exists():
            raise FileNotFoundError(f"ATRP {name} reaction file not found: {path}")

    return paths


def get_atrp_initiation_path() -> Path:
    """Get path to ATRP initiation reaction template."""
    return _REACTIONS_DIR / "atrp_initiation.rxn"


def get_atrp_polymerization_path() -> Path:
    """Get path to ATRP polymerization reaction template."""
    return _REACTIONS_DIR / "atrp_polymerization.rxn"


def get_atrp_termination_path() -> Path:
    """Get path to ATRP termination reaction template."""
    return _REACTIONS_DIR / "atrp_termination.rxn"
