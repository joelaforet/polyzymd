"""Amino acid classification and SASA reference data.

This module provides centralized reference data for amino acid properties:
- Maximum accessible surface area (maxASA) from Tien et al. 2013
- Standard amino acid classification by physicochemical properties
- Default MDAnalysis selection strings for each AA class

These constants are used by:
- Surface exposure filtering (binding preference analysis)
- Protein grouping in contact analysis
- Template generation for analysis configs

References
----------
Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO.
Maximum allowed solvent accessibilities of residues in proteins.
PLoS One. 2013 Nov 21;8(11):e80635.
doi: 10.1371/journal.pone.0080635. PMID: 24278298; PMCID: PMC3836772.
"""

from __future__ import annotations

from enum import Enum
from typing import Final

# Canonical ordering of amino acid classes for consistent display across plots.
# Used by binding preference, contacts, and binding free energy plotters.
CANONICAL_AA_CLASS_ORDER: Final[list[str]] = [
    "aromatic",
    "polar",
    "nonpolar",
    "charged_positive",
    "charged_negative",
]


class AAClass(str, Enum):
    """Standard amino acid classifications."""

    AROMATIC = "aromatic"
    POLAR = "polar"
    NONPOLAR = "nonpolar"
    CHARGED_POSITIVE = "charged_positive"
    CHARGED_NEGATIVE = "charged_negative"
    UNKNOWN = "unknown"


# =============================================================================
# Maximum Accessible Surface Area (maxASA)
# =============================================================================

# Tien et al. 2013 empirical maxASA values (Angstrom^2)
# These represent the maximum solvent-accessible surface area
# for each amino acid type in extended tripeptide conformations.
MAX_ASA_TABLE: Final[dict[str, float]] = {
    "ALA": 121.0,
    "ARG": 265.0,
    "ASN": 187.0,
    "ASP": 187.0,
    "CYS": 148.0,
    "GLU": 214.0,
    "GLN": 214.0,
    "GLY": 97.0,
    "HIS": 216.0,
    "ILE": 195.0,
    "LEU": 191.0,
    "LYS": 230.0,
    "MET": 203.0,
    "PHE": 228.0,
    "PRO": 154.0,
    "SER": 143.0,
    "THR": 163.0,
    "TRP": 264.0,
    "TYR": 255.0,
    "VAL": 165.0,
}


# =============================================================================
# Amino Acid Classification
# =============================================================================

# Standard classification by physicochemical properties
# Note: Histidine is classified as aromatic (contains imidazole ring)
AA_CLASSIFICATION_TABLE: Final[dict[str, str]] = {
    # Aromatic (contain aromatic rings)
    "PHE": "aromatic",
    "TRP": "aromatic",
    "TYR": "aromatic",
    "HIS": "aromatic",  # Imidazole ring
    # Polar (uncharged, can form H-bonds)
    "SER": "polar",
    "THR": "polar",
    "ASN": "polar",
    "GLN": "polar",
    "CYS": "polar",
    # Nonpolar (hydrophobic)
    "ALA": "nonpolar",
    "VAL": "nonpolar",
    "ILE": "nonpolar",
    "LEU": "nonpolar",
    "MET": "nonpolar",
    "GLY": "nonpolar",
    "PRO": "nonpolar",
    # Charged positive (basic)
    "ARG": "charged_positive",
    "LYS": "charged_positive",
    # Charged negative (acidic)
    "ASP": "charged_negative",
    "GLU": "charged_negative",
}

# All standard 3-letter amino acid codes
STANDARD_AA_CODES: Final[list[str]] = list(AA_CLASSIFICATION_TABLE.keys())


# =============================================================================
# MDAnalysis Selection Strings
# =============================================================================

# Default MDAnalysis selections for each AA class
# These are used in analysis.yaml templates and binding preference analysis
DEFAULT_AA_CLASS_SELECTIONS: Final[dict[str, str]] = {
    "aromatic": "protein and resname PHE TRP TYR HIS",
    "polar": "protein and resname SER THR ASN GLN CYS",
    "nonpolar": "protein and resname ALA VAL ILE LEU MET GLY PRO",
    "charged_positive": "protein and resname ARG LYS",
    "charged_negative": "protein and resname ASP GLU",
}

# Residue names for each class (for programmatic access)
AA_CLASS_RESIDUES: Final[dict[str, list[str]]] = {
    "aromatic": ["PHE", "TRP", "TYR", "HIS"],
    "polar": ["SER", "THR", "ASN", "GLN", "CYS"],
    "nonpolar": ["ALA", "VAL", "ILE", "LEU", "MET", "GLY", "PRO"],
    "charged_positive": ["ARG", "LYS"],
    "charged_negative": ["ASP", "GLU"],
}


# =============================================================================
# Helper Functions
# =============================================================================


def get_aa_class(resname: str) -> str:
    """Get amino acid classification for a residue name.

    Parameters
    ----------
    resname : str
        3-letter amino acid code (case-insensitive)

    Returns
    -------
    str
        Classification: 'aromatic', 'polar', 'nonpolar',
        'charged_positive', 'charged_negative', or 'unknown'

    Examples
    --------
    >>> get_aa_class("PHE")
    'aromatic'
    >>> get_aa_class("lys")
    'charged_positive'
    >>> get_aa_class("UNK")
    'unknown'
    """
    # Handle common protonation state variants
    variants = {
        "HIE": "HIS",
        "HID": "HIS",
        "HIP": "HIS",
        "HSE": "HIS",
        "HSD": "HIS",
        "HSP": "HIS",
        "CYSH": "CYS",
        "CYX": "CYS",
        "ASH": "ASP",
        "GLH": "GLU",
        "LYN": "LYS",
    }
    normalized = resname.upper().strip()
    normalized = variants.get(normalized, normalized)

    return AA_CLASSIFICATION_TABLE.get(normalized, "unknown")


def get_max_asa(resname: str) -> float | None:
    """Get maximum accessible surface area for a residue name.

    Parameters
    ----------
    resname : str
        3-letter amino acid code (case-insensitive)

    Returns
    -------
    float or None
        Maximum ASA in Angstrom^2, or None if residue not in table

    Examples
    --------
    >>> get_max_asa("ALA")
    121.0
    >>> get_max_asa("TRP")
    264.0
    >>> get_max_asa("UNK")  # Returns None for unknown residues
    """
    normalized = resname.upper().strip()
    return MAX_ASA_TABLE.get(normalized)


def get_residues_for_class(aa_class: str) -> list[str]:
    """Get all residue names belonging to an amino acid class.

    Parameters
    ----------
    aa_class : str
        One of: 'aromatic', 'polar', 'nonpolar',
        'charged_positive', 'charged_negative'

    Returns
    -------
    list[str]
        List of 3-letter amino acid codes in this class

    Raises
    ------
    ValueError
        If aa_class is not a valid classification

    Examples
    --------
    >>> get_residues_for_class("aromatic")
    ['PHE', 'TRP', 'TYR', 'HIS']
    """
    if aa_class not in AA_CLASS_RESIDUES:
        valid = list(AA_CLASS_RESIDUES.keys())
        raise ValueError(f"Unknown AA class '{aa_class}'. Valid: {valid}")
    return AA_CLASS_RESIDUES[aa_class].copy()


def get_selection_for_class(aa_class: str) -> str:
    """Get MDAnalysis selection string for an amino acid class.

    Parameters
    ----------
    aa_class : str
        One of: 'aromatic', 'polar', 'nonpolar',
        'charged_positive', 'charged_negative'

    Returns
    -------
    str
        MDAnalysis selection string

    Raises
    ------
    ValueError
        If aa_class is not a valid classification

    Examples
    --------
    >>> get_selection_for_class("aromatic")
    'protein and resname PHE TRP TYR HIS'
    """
    if aa_class not in DEFAULT_AA_CLASS_SELECTIONS:
        valid = list(DEFAULT_AA_CLASS_SELECTIONS.keys())
        raise ValueError(f"Unknown AA class '{aa_class}'. Valid: {valid}")
    return DEFAULT_AA_CLASS_SELECTIONS[aa_class]
