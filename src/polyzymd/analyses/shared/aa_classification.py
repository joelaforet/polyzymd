"""Maximum accessible surface area of each amino acid, from Tien et al. 2013.

:func:`get_max_asa` supplies the maximum ASA that
:func:`polyzymd.analyses.functions.residue_occlusion` compares each residue's
SASA against.

References
----------
Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO.
Maximum allowed solvent accessibilities of residues in proteins.
PLoS One. 2013 Nov 21;8(11):e80635.
doi: 10.1371/journal.pone.0080635. PMID: 24278298; PMCID: PMC3836772.
"""

from __future__ import annotations

from typing import Final

#: Residue names of protonation states and disulfide-bonded cysteine, mapped to
#: the standard residue name.
PROTONATION_VARIANTS: Final[dict[str, str]] = {
    "HIE": "HIS",
    "HID": "HIS",
    "HIP": "HIS",
    "HSE": "HIS",
    "HSD": "HIS",
    "HSP": "HIS",
    "CYSH": "CYS",
    "CYX": "CYS",
    "CYM": "CYS",
    "ASH": "ASP",
    "GLH": "GLU",
    "LYN": "LYS",
}


# =============================================================================
# Maximum Accessible Surface Area (maxASA)
# =============================================================================

# Tien et al. 2013 theoretical maxASA values (Angstrom^2), Table 1 column 1:
# the largest ASA of the residue X in a Gly-X-Gly tripeptide over the allowed
# backbone conformations. The authors recommend these for normalizing ASA.
THEORETICAL_MAX_ASA_TABLE: Final[dict[str, float]] = {
    "ALA": 129.0,
    "ARG": 274.0,
    "ASN": 195.0,
    "ASP": 193.0,
    "CYS": 167.0,
    "GLU": 223.0,
    "GLN": 225.0,
    "GLY": 104.0,
    "HIS": 224.0,
    "ILE": 197.0,
    "LEU": 201.0,
    "LYS": 236.0,
    "MET": 224.0,
    "PHE": 240.0,
    "PRO": 159.0,
    "SER": 155.0,
    "THR": 172.0,
    "TRP": 285.0,
    "TYR": 263.0,
    "VAL": 174.0,
}

# Tien et al. 2013 empirical maxASA values (Angstrom^2), Table 1 column 2:
# the largest ASA of each residue type observed in a set of protein structures.
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
# Helper Functions
# =============================================================================


def get_max_asa(resname: str, table: str = "theoretical") -> float | None:
    """Return the maximum accessible surface area of a residue name, in Å².

    Parameters
    ----------
    resname : str
        3-letter amino acid code, case-insensitive. Protonation states such as
        ``HID`` or ``ASH`` and ``CYX`` take the value of the standard residue.
    table : {"theoretical", "empirical"}
        Column of Table 1 of Tien et al. 2013: ``theoretical`` (default, which
        the authors recommend for normalizing ASA) or ``empirical``.

    Returns
    -------
    float or None
        Maximum ASA in Å², or None if the residue is not in the table.

    Examples
    --------
    >>> get_max_asa("ALA")
    129.0
    >>> get_max_asa("HID", table="empirical")
    216.0
    >>> get_max_asa("UNK")  # Returns None for unknown residues
    """
    tables = {"theoretical": THEORETICAL_MAX_ASA_TABLE, "empirical": MAX_ASA_TABLE}
    if table not in tables:
        raise ValueError(f"table must be 'theoretical' or 'empirical', got {table!r}")
    normalized = resname.upper().strip()
    return tables[table].get(PROTONATION_VARIANTS.get(normalized, normalized))
