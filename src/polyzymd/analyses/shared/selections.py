"""Special selection syntax for distance calculations.

This module provides parsing of extended selection syntax for defining
atom positions in distance calculations. It supports:

1. Standard MDAnalysis selections: "resid 133 and name OD1"
2. Midpoint of multiple atoms: "midpoint(resid 133 and name OD1 OD2)"
3. Center of mass of a group: "com(resid 50-75)"

Examples
--------
>>> from polyzymd.analyses.shared.selections import parse_selection_string, get_position
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

import logging
import re
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup
    from MDAnalysis.core.universe import Universe

LOGGER = logging.getLogger(__name__)


# =============================================================================
# Selection Translation (PolyzyMD → MDAnalysis)
# =============================================================================


def translate_selection(selection: str) -> str:
    """Translate PolyzyMD selection keywords to MDAnalysis equivalents.

    This allows users to use the same selection syntax in analysis as they
    use in config.yaml for restraints and other atom selections.

    Translations
    ------------
    - ``pdbindex N`` → ``bynum N``: the N-th atom of the system, counted
      from 1

    ``pdbindex`` means the same as in a config's restraints
    (:mod:`polyzymd.core.restraints`): the atom at position N - 1 of the
    built system, which is the PDB ATOM serial PolyzyMD writes and PyMOL
    shows. MDAnalysis's ``bynum`` selects by that position. ``id`` would
    select by the serial read from the topology file, which differs from the
    position when serials wrap above 99,999 atoms.

    Parameters
    ----------
    selection : str
        Selection string with possible PolyzyMD-specific keywords

    Returns
    -------
    str
        Selection string with MDAnalysis-compatible keywords

    Examples
    --------
    >>> translate_selection("pdbindex 100 and name CA")
    "bynum 100 and name CA"

    >>> translate_selection("midpoint(pdbindex 100 and name OD1 OD2)")
    "midpoint(bynum 100 and name OD1 OD2)"
    """
    translated = re.sub(r"\bpdbindex\b", "bynum", selection, flags=re.IGNORECASE)

    if translated != selection:
        LOGGER.debug("Translated selection keyword: 'pdbindex' → 'bynum' (position from 1)")

    return translated


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

    Also translates PolyzyMD-specific keywords (like ``pdbindex``) to their
    MDAnalysis equivalents (like ``id``).

    Parameters
    ----------
    selection : str
        Selection string, possibly with special syntax:
        - "resid 77 and name OG" - standard MDAnalysis
        - "midpoint(resid 133 and name OD1 OD2)" - midpoint mode
        - "com(resid 50-75)" - center of mass mode
        - "pdbindex 100 and name CA" - PolyzyMD pdbindex (translated to bynum)

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

    >>> parsed = parse_selection_string("pdbindex 100 and name CA")
    >>> parsed.selection
    "id 100 and name CA"
    """
    selection = selection.strip()

    # Check for midpoint() syntax
    midpoint_match = _MIDPOINT_PATTERN.match(selection)
    if midpoint_match:
        inner = midpoint_match.group(1).strip()
        return ParsedSelection(
            selection=translate_selection(inner),
            mode=SelectionMode.MIDPOINT,
            original=selection,
        )

    # Check for com() syntax
    com_match = _COM_PATTERN.match(selection)
    if com_match:
        inner = com_match.group(1).strip()
        return ParsedSelection(
            selection=translate_selection(inner),
            mode=SelectionMode.COM,
            original=selection,
        )

    # Standard MDAnalysis selection (also translate)
    return ParsedSelection(
        selection=translate_selection(selection),
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
    from polyzymd.analyses.shared.diagnostics import get_selection_diagnostics

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
        if len(atoms) > 1:
            raise ValueError(
                "SelectionMode.SINGLE requires exactly one atom, "
                f"but selection matched {len(atoms)} atoms. "
                "Use SelectionMode.MIDPOINT or SelectionMode.COM for multi-atom selections"
            )
        raise ValueError("Cannot compute position for empty atom selection")

    elif mode in (SelectionMode.CENTROID, SelectionMode.MIDPOINT):
        return atoms.center_of_geometry().astype(np.float64)

    elif mode == SelectionMode.COM:
        return atoms.center_of_mass().astype(np.float64)

    else:
        raise ValueError(f"Unknown selection mode: {mode}")
