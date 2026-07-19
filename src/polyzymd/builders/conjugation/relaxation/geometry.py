"""Pure geometry helpers for conjugate relaxation diagnostics."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np

LOGGER = logging.getLogger(__name__)
_POSITION_CONVERSION_ERRORS = (AttributeError, TypeError, ValueError)


def positions_to_numpy(positions: Any, unit_module: Any | None = None) -> np.ndarray:
    """Return positions as a floating point nanometer array.

    Parameters
    ----------
    positions : Any
        Coordinate container or OpenMM-like quantity.
    unit_module : Any or None, optional
        Unit module exposing ``nanometer`` for unit-aware extraction, by default ``None``.

    Returns
    -------
    numpy.ndarray
        Coordinate array in nanometers when units are available.
    """
    conversion_error: Exception | None = None
    if unit_module is not None and hasattr(positions, "value_in_unit"):
        try:
            return np.asarray(positions.value_in_unit(unit_module.nanometer), dtype=float)
        except _POSITION_CONVERSION_ERRORS as exc:
            conversion_error = exc
    if hasattr(positions, "m_as"):
        try:
            return np.asarray(positions.m_as("nanometer"), dtype=float)
        except _POSITION_CONVERSION_ERRORS as exc:
            conversion_error = exc
    if conversion_error is not None:
        LOGGER.warning(
            "Falling back to raw np.asarray() for positions of type %s after unit-aware "
            "coordinate conversion failed: %s",
            type(positions).__name__,
            conversion_error,
        )
    return np.asarray(positions, dtype=float)


def coordinate_span_nm(positions: Any, unit_module: Any | None = None) -> float:
    """Return the maximum per-axis coordinate span in nanometers."""
    coords = positions_to_numpy(positions, unit_module)
    if coords.size == 0:
        return 0.0
    return float(np.max(np.ptp(coords, axis=0)))


def require_finite_coordinates(
    positions: Any,
    unit_module: Any | None = None,
    *,
    label: str = "coordinates",
) -> np.ndarray:
    """Return coordinates after verifying shape and finite values."""
    coords = positions_to_numpy(positions, unit_module)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise RuntimeError(f"{label} have invalid shape: {coords.shape}")
    if not np.all(np.isfinite(coords)):
        raise RuntimeError(f"{label} contain non-finite values")
    return coords


def replace_pdb_coordinates(
    template_pdb_path: Path | str,
    coordinates_angstrom: Any,
    output_path: Path | str,
) -> Path:
    """Write a PDB by replacing only ATOM/HETATM coordinate columns.

    Parameters
    ----------
    template_pdb_path : pathlib.Path or str
        Authoritative PDB whose atom, residue, chain, TER, LINK, and CONECT
        records must be preserved.
    coordinates_angstrom : Any
        Replacement coordinate array with one finite XYZ triplet per atom.
    output_path : pathlib.Path or str
        Destination for the coordinate-updated PDB.

    Returns
    -------
    pathlib.Path
        Path to the written PDB.

    Raises
    ------
    ValueError
        If the coordinate count, shape, finiteness, or fixed-width coordinate
        range is invalid for the template atom table.
    """
    template_path = Path(template_pdb_path)
    output = Path(output_path)
    coords = np.asarray(coordinates_angstrom, dtype=float)
    atom_lines = tuple(
        line
        for line in template_path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith(("ATOM", "HETATM"))
    )
    if coords.shape != (len(atom_lines), 3):
        raise ValueError(
            "Coordinate replacement requires one XYZ triplet per PDB atom: "
            f"expected {(len(atom_lines), 3)}, got {coords.shape}"
        )
    if not np.all(np.isfinite(coords)):
        raise ValueError("Replacement coordinates contain non-finite values")

    output.parent.mkdir(parents=True, exist_ok=True)
    coord_index = 0
    out_lines: list[str] = []
    with template_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(("ATOM", "HETATM")):
                x_coord, y_coord, z_coord = (
                    _format_pdb_coordinate(value) for value in coords[coord_index]
                )
                padded = f"{line:<80}"
                line = f"{padded[:30]}{x_coord}{y_coord}{z_coord}{padded[54:]}"
                coord_index += 1
            out_lines.append(line)
    output.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return output


def _format_pdb_coordinate(value: float) -> str:
    """Return one fixed-width PDB coordinate field.

    Parameters
    ----------
    value : float
        Coordinate value in angstroms.

    Returns
    -------
    str
        Coordinate formatted as an exactly 8-character PDB field.

    Raises
    ------
    ValueError
        If the coordinate is non-finite or cannot fit after three-decimal
        rounding.
    """
    if not np.isfinite(value):
        raise ValueError("Replacement coordinates contain non-finite values")
    field = f"{value:8.3f}"
    if len(field) != 8:
        raise ValueError(
            "Replacement coordinate cannot fit in an 8-character PDB field: "
            f"value={value!r}, formatted={field!r}"
        )
    return field
