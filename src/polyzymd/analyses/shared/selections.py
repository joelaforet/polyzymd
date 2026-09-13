"""Extended selection syntax for pair-distance endpoints.

A pair endpoint is one point, and this module says which point a selection
string means. Three forms are understood:

1. A standard MDAnalysis selection of one atom, ``"resid 77 and name OG"``.
2. The midpoint of several atoms, ``"midpoint(resid 133 and name OD1 OD2)"``,
   which is how a carboxyl group is usually treated.
3. The centre of mass of a group, ``"com(resid 50-75)"``, which is how a
   domain such as a lipase lid is usually treated.

:func:`parse_selection_string` splits the wrapper from the selection and
:func:`get_position` reduces the selected atoms to one position.
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
    - ``pdbindex N`` → ``id N`` (PDB ATOM serial number)

    The ``pdbindex`` keyword refers to the 1-indexed atom serial number
    from the PDB ATOM record (column 7-11), which is what PyMOL displays
    as "id". In MDAnalysis, this is accessed via the ``id`` selection keyword.

    Note: MDAnalysis also has ``bynum`` which is 1-indexed *positional*
    (i.e., bynum 1 = first atom, bynum 2 = second atom), but this does NOT
    correspond to PDB serial numbers when there are gaps in numbering.
    We use ``id`` because it matches actual PDB serial numbers.

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
    "id 100 and name CA"

    >>> translate_selection("midpoint(pdbindex 100 and name OD1 OD2)")
    "midpoint(id 100 and name OD1 OD2)"
    """
    # pdbindex N → id N (PDB ATOM serial number)
    translated = re.sub(r"\bpdbindex\b", "id", selection, flags=re.IGNORECASE)

    if translated != selection:
        LOGGER.debug("Translated selection keyword: 'pdbindex' → 'id' (PDB ATOM serial number)")

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
        - "pdbindex 100 and name CA" - PolyzyMD pdbindex (translated to id)

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
