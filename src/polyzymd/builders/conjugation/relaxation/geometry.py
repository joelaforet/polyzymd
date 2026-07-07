"""Pure geometry helpers for conjugate relaxation diagnostics."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


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
    if unit_module is not None and hasattr(positions, "value_in_unit"):
        return np.asarray(positions.value_in_unit(unit_module.nanometer), dtype=float)
    if hasattr(positions, "m_as"):
        return np.asarray(positions.m_as("nanometer"), dtype=float)
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
    """Write a PDB by replacing only ATOM/HETATM coordinate columns."""
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
                x_coord, y_coord, z_coord = coords[coord_index]
                padded = f"{line:<80}"
                line = f"{padded[:30]}{x_coord:8.3f}{y_coord:8.3f}{z_coord:8.3f}{padded[54:]}"
                coord_index += 1
            out_lines.append(line)
    output.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return output
