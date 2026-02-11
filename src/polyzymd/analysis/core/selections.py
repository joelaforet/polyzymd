"""Special selection syntax for distance calculations.

This module provides parsing of extended selection syntax for defining
atom positions in distance calculations. It supports:

1. Standard MDAnalysis selections: "resid 133 and name OD1"
2. Midpoint of multiple atoms: "midpoint(resid 133 and name OD1 OD2)"
3. Center of mass of a group: "com(resid 50-75)"

Examples
--------
>>> from polyzymd.analysis.core.selections import parse_selection, get_position
>>>
>>> # Standard selection - single atom
>>> ag = parse_selection(universe, "resid 77 and name OG")
>>> pos = get_position(ag)  # Returns position of single atom
>>>
>>> # Midpoint of Asp carboxyl oxygens
>>> ag = parse_selection(universe, "midpoint(resid 133 and name OD1 OD2)")
>>> pos = get_position(ag)  # Returns midpoint of OD1 and OD2
>>>
>>> # Center of mass of lid domain
>>> ag = parse_selection(universe, "com(resid 50-75)")
>>> pos = get_position(ag)  # Returns COM of residues 50-75

Notes
-----
The `midpoint()` syntax is particularly useful for catalytic residues where
the functional position is between two atoms (e.g., Asp carboxyl oxygens,
Glu carboxyl oxygens).

The `com()` syntax is useful for domain motions where you want to track
the center of mass of a group of residues (e.g., lid opening in lipases).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup
    from MDAnalysis.core.universe import Universe


class SelectionMode(str, Enum):
    """Mode for position calculation from atom selection."""

    SINGLE = "single"  # Single atom position
    CENTROID = "centroid"  # Center of geometry (default for multiple atoms)
    MIDPOINT = "midpoint"  # Explicit midpoint (same as centroid but explicit)
    COM = "com"  # Center of mass


@dataclass
class ParsedSelection:
    """Result of parsing a selection string.

    Attributes
    ----------
    selection : str
        The MDAnalysis selection string (without wrapper function)
    mode : SelectionMode
        How to compute the position
    original : str
        The original input string
    """

    selection: str
    mode: SelectionMode
    original: str

    def __str__(self) -> str:
        return self.original


# Regex patterns for special syntax
_MIDPOINT_PATTERN = re.compile(r"^midpoint\s*\(\s*(.+)\s*\)$", re.IGNORECASE)
_COM_PATTERN = re.compile(r"^com\s*\(\s*(.+)\s*\)$", re.IGNORECASE)


def parse_selection_string(selection: str) -> ParsedSelection:
    """Parse a selection string to extract mode and MDAnalysis selection.

    Parameters
    ----------
    selection : str
        Selection string, possibly with special syntax:
        - "resid 77 and name OG" - standard MDAnalysis
        - "midpoint(resid 133 and name OD1 OD2)" - midpoint mode
        - "com(resid 50-75)" - center of mass mode

    Returns
    -------
    ParsedSelection
        Parsed selection with mode and clean selection string

    Examples
    --------
    >>> parsed = parse_selection_string("midpoint(resid 133 and name OD1 OD2)")
    >>> parsed.mode
    <SelectionMode.MIDPOINT: 'midpoint'>
    >>> parsed.selection
    "resid 133 and name OD1 OD2"
    """
    selection = selection.strip()

    # Check for midpoint() syntax
    midpoint_match = _MIDPOINT_PATTERN.match(selection)
    if midpoint_match:
        inner = midpoint_match.group(1).strip()
        return ParsedSelection(
            selection=inner,
            mode=SelectionMode.MIDPOINT,
            original=selection,
        )

    # Check for com() syntax
    com_match = _COM_PATTERN.match(selection)
    if com_match:
        inner = com_match.group(1).strip()
        return ParsedSelection(
            selection=inner,
            mode=SelectionMode.COM,
            original=selection,
        )

    # Standard MDAnalysis selection
    return ParsedSelection(
        selection=selection,
        mode=SelectionMode.SINGLE,
        original=selection,
    )


def select_atoms(universe: "Universe", selection: str) -> "AtomGroup":
    """Select atoms from universe using potentially special syntax.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe
    selection : str
        Selection string (standard or special syntax)

    Returns
    -------
    AtomGroup
        Selected atoms

    Raises
    ------
    ValueError
        If selection matches no atoms, with diagnostic info
    """
    from polyzymd.analysis.core.diagnostics import get_selection_diagnostics

    parsed = parse_selection_string(selection)
    atoms = universe.select_atoms(parsed.selection)

    if len(atoms) == 0:
        diag = get_selection_diagnostics(universe, selection)
        raise ValueError(f"Selection '{selection}' matched no atoms.\n\n{diag}")

    return atoms


def get_position(
    atoms: "AtomGroup",
    mode: SelectionMode = SelectionMode.SINGLE,
) -> NDArray[np.float64]:
    """Get position from atom group based on mode.

    Parameters
    ----------
    atoms : AtomGroup
        MDAnalysis AtomGroup
    mode : SelectionMode
        How to compute position:
        - SINGLE: Position of single atom (error if multiple)
        - CENTROID/MIDPOINT: Center of geometry
        - COM: Center of mass

    Returns
    -------
    NDArray[np.float64]
        3D position vector [x, y, z]

    Raises
    ------
    ValueError
        If mode is SINGLE but multiple atoms selected
    """
    if mode == SelectionMode.SINGLE:
        if len(atoms) == 1:
            return atoms.positions[0].astype(np.float64)
        else:
            # Default to centroid for multiple atoms
            return atoms.center_of_geometry().astype(np.float64)

    elif mode in (SelectionMode.CENTROID, SelectionMode.MIDPOINT):
        return atoms.center_of_geometry().astype(np.float64)

    elif mode == SelectionMode.COM:
        return atoms.center_of_mass().astype(np.float64)

    else:
        raise ValueError(f"Unknown selection mode: {mode}")


def get_position_from_selection(
    universe: "Universe",
    selection: str,
) -> NDArray[np.float64]:
    """Get position from selection string in one step.

    This is a convenience function that combines parsing, selection,
    and position calculation.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe
    selection : str
        Selection string (standard or special syntax)

    Returns
    -------
    NDArray[np.float64]
        3D position vector [x, y, z]

    Examples
    --------
    >>> # Single atom
    >>> pos = get_position_from_selection(u, "resid 77 and name OG")
    >>>
    >>> # Midpoint of Asp carboxyl
    >>> pos = get_position_from_selection(u, "midpoint(resid 133 and name OD1 OD2)")
    """
    from polyzymd.analysis.core.diagnostics import get_selection_diagnostics

    parsed = parse_selection_string(selection)
    atoms = universe.select_atoms(parsed.selection)

    if len(atoms) == 0:
        diag = get_selection_diagnostics(universe, selection)
        raise ValueError(f"Selection '{selection}' matched no atoms.\n\n{diag}")

    return get_position(atoms, parsed.mode)


def validate_selection(universe: "Universe", selection: str) -> dict:
    """Validate a selection string and return diagnostic info.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe
    selection : str
        Selection string to validate

    Returns
    -------
    dict
        Diagnostic information:
        - valid: bool
        - n_atoms: int
        - mode: str
        - atoms: list of atom info dicts
        - error: str (if invalid)
        - diagnostics: str (detailed diagnostics if invalid)
    """
    from polyzymd.analysis.core.diagnostics import get_selection_diagnostics

    try:
        parsed = parse_selection_string(selection)
        atoms = universe.select_atoms(parsed.selection)

        if len(atoms) == 0:
            diag = get_selection_diagnostics(universe, selection)
            return {
                "valid": False,
                "error": f"Selection matched no atoms: {parsed.selection}",
                "diagnostics": diag,
                "mode": parsed.mode.value,
                "n_atoms": 0,
            }

        atom_info = []
        for atom in atoms[:10]:  # Limit to first 10
            atom_info.append(
                {
                    "name": atom.name,
                    "resname": atom.resname,
                    "resid": atom.resid,
                    "index": atom.index,
                }
            )

        return {
            "valid": True,
            "n_atoms": len(atoms),
            "mode": parsed.mode.value,
            "atoms": atom_info,
            "truncated": len(atoms) > 10,
        }

    except Exception as e:
        return {
            "valid": False,
            "error": str(e),
            "n_atoms": 0,
        }


def format_selection_for_label(selection: str) -> str:
    """Convert selection string to a short label for filenames/display.

    Parameters
    ----------
    selection : str
        Selection string (standard or special syntax)

    Returns
    -------
    str
        Short label (e.g., "Asp133_mid" or "Ser77_OG")

    Examples
    --------
    >>> format_selection_for_label("midpoint(resid 133 and name OD1 OD2)")
    "res133_mid"
    >>> format_selection_for_label("resid 77 and name OG")
    "res77_OG"
    """
    parsed = parse_selection_string(selection)

    # Extract resid
    resid_match = re.search(r"resid\s+(\d+)", parsed.selection, re.IGNORECASE)
    resid_str = f"res{resid_match.group(1)}" if resid_match else ""

    # Extract atom name(s)
    name_match = re.search(r"name\s+(\w+(?:\s+\w+)*)", parsed.selection, re.IGNORECASE)
    if name_match:
        names = name_match.group(1).split()
        if len(names) == 1:
            name_str = names[0].upper()
        else:
            name_str = "mid" if parsed.mode == SelectionMode.MIDPOINT else "cog"
    else:
        name_str = "com" if parsed.mode == SelectionMode.COM else "grp"

    if resid_str and name_str:
        return f"{resid_str}_{name_str}"
    elif resid_str:
        return resid_str
    else:
        return name_str or "sel"
