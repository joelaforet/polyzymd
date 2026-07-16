"""
Default reaction templates for polymer generation.

This module provides access to bundled ATRP (Atom-Transfer Radical Polymerization)
reaction templates for methacrylate-based monomers.

The reaction templates define:
- Initiation: Chlorination of the vinyl group
- Polymerization: Chain extension via radical coupling
- Termination: Restoration of terminal alkene

Example:
    >>> from polyzymd.data.reactions import get_atrp_reaction_paths
    >>> reactions = get_atrp_reaction_paths()
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


def is_default_atrp_reaction_set(
    initiation: Path | str,
    polymerization: Path | str,
    termination: Path | str,
) -> bool:
    """Return whether paths are the bundled default ATRP reaction set.

    Parameters
    ----------
    initiation : pathlib.Path or str
        Candidate initiation reaction path.
    polymerization : pathlib.Path or str
        Candidate polymerization reaction path.
    termination : pathlib.Path or str
        Candidate termination reaction path.

    Returns
    -------
    bool
        ``True`` when all three paths resolve to the bundled defaults.
    """
    defaults = get_atrp_reaction_paths()
    candidates = {
        "initiation": _normalize_reaction_candidate(
            initiation, role="initiation", defaults=defaults
        ),
        "polymerization": _normalize_reaction_candidate(
            polymerization,
            role="polymerization",
            defaults=defaults,
        ),
        "termination": _normalize_reaction_candidate(
            termination, role="termination", defaults=defaults
        ),
    }
    return all(candidates[key].resolve() == defaults[key].resolve() for key in defaults)


def _normalize_reaction_candidate(
    candidate: Path | str,
    *,
    role: str,
    defaults: Dict[str, Path],
) -> Path:
    """Normalize a reaction candidate before default-set comparison.

    Parameters
    ----------
    candidate : pathlib.Path or str
        Candidate reaction path or literal default token.
    role : str
        Reaction role used to select the packaged default path.
    defaults : dict[str, pathlib.Path]
        Packaged default reaction paths by role.

    Returns
    -------
    pathlib.Path
        Packaged default path for literal ``default`` candidates, otherwise the
        candidate converted to a path.
    """
    text = str(candidate).strip()
    if text.lower() == "default":
        return defaults[role]
    return Path(candidate)
