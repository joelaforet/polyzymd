"""Diagnostic utilities for analysis troubleshooting.

This module provides helpers for generating actionable error messages
when analysis operations fail, helping users quickly identify and fix
issues with their inputs.

Functions
---------
get_selection_diagnostics
    Generate diagnostic info for failed atom selections
get_residue_info
    Get information about a specific residue
get_protein_residue_range
    Get the range of residue IDs in the protein
format_diagnostic_message
    Format diagnostic lines into a readable message
warn_if_multi_chain_selection
    Warn if a selection matched atoms from multiple chains
"""

from __future__ import annotations

import logging
import re
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from MDAnalysis import Universe
    from MDAnalysis.core.groups import AtomGroup

LOGGER = logging.getLogger(__name__)


def get_protein_residue_range(universe: "Universe") -> Optional[tuple[int, int]]:
    """Get the min and max residue IDs for protein atoms.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe

    Returns
    -------
    tuple[int, int] or None
        (min_resid, max_resid) or None if no protein atoms found
    """
    protein = universe.select_atoms("protein")
    if len(protein) == 0:
        return None
    resids = protein.residues.resids
    return int(resids.min()), int(resids.max())


def get_residue_info(
    universe: "Universe",
    resid: int,
) -> Optional[dict]:
    """Get information about a specific residue.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe
    resid : int
        Residue ID to query

    Returns
    -------
    dict or None
        Dictionary with 'resname', 'resid', 'atom_names' keys,
        or None if residue doesn't exist
    """
    residue_sel = universe.select_atoms(f"resid {resid}")
    if len(residue_sel) == 0:
        return None

    return {
        "resid": resid,
        "resname": residue_sel.residues[0].resname,
        "atom_names": sorted(set(residue_sel.atoms.names)),
    }


def format_diagnostic_message(lines: list[str], header: str = "Diagnostic info:") -> str:
    """Format diagnostic lines into a readable message.

    Parameters
    ----------
    lines : list[str]
        Individual diagnostic lines
    header : str
        Header text for the diagnostic block

    Returns
    -------
    str
        Formatted diagnostic message
    """
    if not lines:
        return ""
    return f"{header}\n  - " + "\n  - ".join(lines)


def get_selection_diagnostics(
    universe: "Universe",
    selection: str,
    context: Optional[str] = None,
) -> str:
    """Generate diagnostic info for a failed atom selection.

    Analyzes the selection string and universe to provide actionable
    guidance on why a selection matched no atoms.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe
    selection : str
        The selection string that failed
    context : str, optional
        Additional context (e.g., "for pair 'His156-Ser77'")

    Returns
    -------
    str
        Formatted diagnostic message with suggestions

    Examples
    --------
    >>> diag = get_selection_diagnostics(u, "resid 999 and name OG")
    >>> print(diag)
    Diagnostic info:
      - Protein residue IDs range from 1 to 267
      - Residue 999 does not exist in this structure

    >>> diag = get_selection_diagnostics(u, "resid 77 and name XYZ")
    >>> print(diag)
    Diagnostic info:
      - Residue 77 is SER
      - Available atoms: C, CA, CB, N, O, OG
      - Atom name 'XYZ' not found in residue 77
    """
    lines = []

    # Add context if provided
    if context:
        lines.append(context)

    # Get protein residue range
    residue_range = get_protein_residue_range(universe)
    if residue_range:
        min_res, max_res = residue_range
        lines.append(f"Protein residue IDs range from {min_res} to {max_res}")

    # Check if selection contains a resid
    resid_match = re.search(r"resid\s+(\d+)", selection)
    if resid_match:
        target_resid = int(resid_match.group(1))
        res_info = get_residue_info(universe, target_resid)

        if res_info is None:
            lines.append(f"Residue {target_resid} does not exist in this structure")
        else:
            # Residue exists - provide atom info
            lines.append(f"Residue {target_resid} is {res_info['resname']}")
            lines.append(f"Available atoms: {', '.join(res_info['atom_names'])}")

            # Check if they specified an atom name that doesn't exist
            # Handle both single names and multiple names (e.g., "name OD1 OD2")
            name_matches = re.findall(r"name\s+([\w\s]+?)(?:\s+and|\s*$)", selection)
            if name_matches:
                # Parse atom names from the match
                specified_names = []
                for match in name_matches:
                    specified_names.extend(match.split())

                missing_names = [
                    name for name in specified_names if name not in res_info["atom_names"]
                ]
                if missing_names:
                    for name in missing_names:
                        lines.append(f"Atom name '{name}' not found in residue {target_resid}")

    # Check for common selection keywords that might be misspelled
    common_keywords = ["protein", "backbone", "name", "resid", "resname", "segid"]
    selection_lower = selection.lower()
    for keyword in common_keywords:
        # Look for close misspellings
        if keyword not in selection_lower:
            # Check for common typos
            typos = _get_common_typos(keyword)
            for typo in typos:
                if typo in selection_lower:
                    lines.append(f"Possible typo: '{typo}' - did you mean '{keyword}'?")
                    break

    # If we couldn't provide specific diagnostics, give generic advice
    if len(lines) == 0 or (len(lines) == 1 and context):
        lines.append("Check your selection syntax and residue/atom naming")
        lines.append("Common selections: 'protein', 'backbone', 'protein and name CA'")

    return format_diagnostic_message(lines)


def _get_common_typos(keyword: str) -> list[str]:
    """Get common typos for MDAnalysis selection keywords.

    Parameters
    ----------
    keyword : str
        Correct keyword

    Returns
    -------
    list[str]
        List of common misspellings
    """
    typo_map = {
        "protein": ["protien", "protine", "prtein", "protin"],
        "backbone": ["backbon", "backboen", "bakcbone"],
        "resid": ["resdi", "redid", "ressid"],
        "resname": ["resnam", "resnme", "resanme"],
        "name": ["nmae", "naem", "nam"],
        "segid": ["segdi", "segi"],
    }
    return typo_map.get(keyword, [])


def validate_equilibration_time(
    equilibration_ns: float,
    trajectory_ns: float,
    warn_threshold: float = 0.5,
) -> tuple[bool, Optional[str]]:
    """Validate equilibration time against trajectory length.

    Parameters
    ----------
    equilibration_ns : float
        Equilibration time in nanoseconds
    trajectory_ns : float
        Total trajectory length in nanoseconds
    warn_threshold : float
        Fraction of trajectory that triggers a warning (default 0.5 = 50%)

    Returns
    -------
    tuple[bool, str or None]
        (is_valid, error_or_warning_message)
        - If equilibration >= trajectory: (False, error_message)
        - If equilibration > threshold: (True, warning_message)
        - Otherwise: (True, None)

    Examples
    --------
    >>> valid, msg = validate_equilibration_time(500.0, 100.0)
    >>> valid
    False
    >>> "exceeds" in msg
    True

    >>> valid, msg = validate_equilibration_time(80.0, 100.0)
    >>> valid
    True
    >>> "80.0%" in msg
    True
    """
    if equilibration_ns >= trajectory_ns:
        return (
            False,
            f"Equilibration time ({equilibration_ns:.1f} ns) exceeds or equals "
            f"trajectory length ({trajectory_ns:.1f} ns). "
            f"Reduce --eq-time or check your trajectory.",
        )

    fraction = equilibration_ns / trajectory_ns
    if fraction > warn_threshold:
        percentage = fraction * 100
        return (
            True,
            f"Warning: Skipping {percentage:.1f}% of trajectory for equilibration "
            f"({equilibration_ns:.1f} ns of {trajectory_ns:.1f} ns). "
            f"This leaves only {trajectory_ns - equilibration_ns:.1f} ns for analysis.",
        )

    return (True, None)


def warn_if_multi_chain_selection(
    atoms: "AtomGroup",
    selection: str,
    context: str = "",
) -> bool:
    """Warn if a selection matched atoms from multiple chains.

    Residue numbers restart at 1 for each chain in PolyzyMD systems.
    A selection like ``resid 141-148`` without chain restriction will
    match residues from protein (chain A), polymers (chain C), and
    potentially water (chain D+), leading to incorrect calculations.

    This function checks if the selected atoms span multiple chains
    and logs a warning if so.

    Parameters
    ----------
    atoms : AtomGroup
        The atoms selected by MDAnalysis
    selection : str
        The original selection string (for error message)
    context : str, optional
        Additional context for the warning message (e.g., "for pair 'Lid Domain'")

    Returns
    -------
    bool
        True if atoms span multiple chains (warning was issued),
        False if atoms are from a single chain (no warning)

    Examples
    --------
    >>> atoms = u.select_atoms("resid 141-148")  # May match multiple chains!
    >>> warn_if_multi_chain_selection(atoms, "resid 141-148", "for Lid Domain distance")
    True  # Warning logged

    >>> atoms = u.select_atoms("protein and resid 141-148")  # Single chain
    >>> warn_if_multi_chain_selection(atoms, "protein and resid 141-148")
    False  # No warning

    Notes
    -----
    PolyzyMD uses a chain convention:
    - Chain A: Protein (enzyme)
    - Chain B: Substrate (small molecule)
    - Chain C: Polymer
    - Chain D+: Solvent, ions

    When analyzing protein residues, always use ``protein and resid X``
    or ``chainid A and resid X`` to avoid accidentally including atoms
    from polymers or water that happen to have the same residue numbers.
    """
    if len(atoms) == 0:
        return False

    # Get unique chain IDs
    try:
        chain_ids = set(atom.chainID for atom in atoms)
    except AttributeError:
        # Fallback to segment IDs if chainID not available
        chain_ids = set(atoms.segments.segids)

    if len(chain_ids) > 1:
        # Get residue names to help user understand what was selected
        resnames = set(atoms.resnames)
        resnames_str = ", ".join(sorted(resnames)[:5])
        if len(resnames) > 5:
            resnames_str += f", ... ({len(resnames)} total)"

        context_str = f" {context}" if context else ""

        LOGGER.warning(
            f"Selection '{selection}' matched atoms from multiple chains: {sorted(chain_ids)}.{context_str}\n"
            f"  Residue types selected: {resnames_str}\n"
            f"  This may cause incorrect calculations (e.g., COM includes polymer/water atoms).\n"
            f"  Consider adding 'protein and' or 'chainid A and' to restrict to a single chain.\n"
            f"  Example: 'com(protein and resid 141-148)' instead of 'com(resid 141-148)'"
        )
        return True

    return False
