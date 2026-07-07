"""Compact diagnostics for restrained OpenMM product checks."""

from __future__ import annotations

import json
import logging
import math
from pathlib import Path
from typing import Any

import numpy as np

from polyzymd.builders.conjugation.relaxation._linkages import _crosslink_bond_diagnostics
from polyzymd.builders.conjugation.relaxation.models import (
    CrosslinkBondDiagnostic,
    GeometryPairDiagnostic,
    OpenMMValidationPhaseDiagnostics,
    ProductGeometryDiagnostics,
)

LOGGER = logging.getLogger(__name__)
_POSITION_CONVERSION_ERRORS = (AttributeError, TypeError, ValueError)


def _is_heavy_atom(atom: Any) -> bool:
    """Return whether an OpenMM-like atom is non-hydrogen."""
    element = getattr(atom, "element", None)
    symbol = getattr(element, "symbol", None)
    if symbol is not None:
        return str(symbol).upper() != "H"
    return not str(getattr(atom, "name", "")).upper().startswith("H")


def validate_finite_energy(value: float, *, label: str = "energy") -> None:
    """Validate that an energy-like scalar is finite."""
    if not math.isfinite(float(value)):
        raise RuntimeError(f"OpenMM validation produced non-finite {label}: {value}")


def validate_finite_positions(
    positions: Any,
    unit_module: Any | None = None,
    *,
    label: str,
    max_span_nm: float | None = 50.0,
) -> float:
    """Validate that coordinates are finite and contained in a realistic span.

    Parameters
    ----------
    positions : Any
        Coordinate array or OpenMM-like quantity in nanometers.
    unit_module : Any or None, optional
        OpenMM unit module used to extract nanometer values, by default ``None``.
    label : str
        Human-readable label included in validation errors.
    max_span_nm : float or None, optional
        Maximum allowed per-axis coordinate span in nanometers. Set to ``None``
        to skip the span check, by default 50.0.

    Returns
    -------
    float
        Maximum per-axis coordinate span in nanometers.
    """
    array = _positions_to_numpy(positions, unit_module)
    if array.size == 0:
        raise RuntimeError(f"OpenMM validation produced empty {label}")
    if not np.all(np.isfinite(array)):
        raise RuntimeError(f"OpenMM validation produced non-finite {label}")
    if array.ndim != 2 or array.shape[1] != 3:
        raise RuntimeError(f"OpenMM validation produced invalid {label} shape: {array.shape}")
    span_nm = float(np.max(np.ptp(array, axis=0)))
    if max_span_nm is not None and span_nm > max_span_nm:
        raise RuntimeError(
            f"OpenMM validation produced unrealistic coordinate span for {label}: "
            f"{span_nm:.3f} nm exceeds {max_span_nm:.3f} nm. "
            "The restrained OpenMM relaxation likely became unstable; refusing to solvate "
            "expanded post-relaxation coordinates."
        )
    return span_nm


def analyze_product_geometry(
    topology: Any,
    positions: Any,
    unit_module: Any | None = None,
    *,
    crosslinked_pdb_path: Path | str | None = None,
    attachment_specs: tuple[Any, ...] = (),
    heavy_heavy_close_nm: float = 0.12,
    h_heavy_close_nm: float = 0.08,
    max_pairs: int = 50,
) -> ProductGeometryDiagnostics:
    """Measure geometry diagnostics before restrained OpenMM validation.

    Parameters
    ----------
    topology : Any
        OpenMM topology-like object with atoms and bonds.
    positions : Any
        Coordinates in nanometers or an OpenMM-compatible quantity.
    unit_module : Any or None, optional
        OpenMM unit module used for coordinate conversion, by default ``None``.
    crosslinked_pdb_path : pathlib.Path, str, or None, optional
        Product PDB used to measure crosslink lengths, by default ``None``.
    attachment_specs : tuple of Any, optional
        Attachment specs carrying resolved plans, by default ``()``.
    heavy_heavy_close_nm : float, optional
        Heavy-heavy close-contact threshold, by default 0.12.
    h_heavy_close_nm : float, optional
        Hydrogen-heavy close-contact threshold, by default 0.08.
    max_pairs : int, optional
        Maximum number of pair diagnostics to retain per category, by default 50.

    Returns
    -------
    ProductGeometryDiagnostics
        Structured geometry report.
    """
    coords = _positions_to_numpy(positions, unit_module)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise RuntimeError(f"Product geometry coordinates have invalid shape: {coords.shape}")
    if not np.all(np.isfinite(coords)):
        raise RuntimeError("Product geometry coordinates contain non-finite values")
    atom_records = tuple(topology.atoms())
    identities = tuple(_topology_atom_identity(atom) for atom in atom_records)
    heavy = tuple(_is_heavy_atom(atom) for atom in atom_records)
    bonded_pairs = _topology_bond_index_pairs(topology)
    bonded_set = {tuple(sorted(pair)) for pair in bonded_pairs}
    min_heavy_heavy, min_h_heavy, contacts = _close_contact_diagnostics(
        coords,
        heavy,
        identities,
        bonded_set=bonded_set,
        heavy_heavy_close_nm=heavy_heavy_close_nm,
        h_heavy_close_nm=h_heavy_close_nm,
        max_pairs=max_pairs,
    )
    outliers = _bonded_distance_outliers(coords, identities, bonded_pairs, max_pairs=max_pairs)
    crosslinks = _crosslink_bond_diagnostics(crosslinked_pdb_path, attachment_specs)
    return ProductGeometryDiagnostics(
        atom_count=int(coords.shape[0]),
        coordinate_span_nm=float(np.max(np.ptp(coords, axis=0))) if coords.size else 0.0,
        min_heavy_heavy_distance_nm=min_heavy_heavy,
        min_h_heavy_distance_nm=min_h_heavy,
        close_contacts=tuple(contacts),
        bonded_distance_outliers=tuple(outliers),
        crosslink_bonds=tuple(crosslinks),
    )


def _validation_phase_diagnostics(
    phase: str,
    topology: Any,
    positions: Any,
    unit_module: Any | None,
    *,
    potential_energy_kj_mol: float | None,
    max_force_kj_mol_nm: float | None = None,
    attachment_specs: tuple[Any, ...] = (),
) -> OpenMMValidationPhaseDiagnostics:
    """Collect finite, span, contact, bond, and linkage diagnostics for a phase."""
    coords = _positions_to_numpy(positions, unit_module)
    has_nan = bool(np.any(np.isnan(coords))) if coords.size else False
    has_inf = bool(np.any(np.isinf(coords))) if coords.size else False
    finite = bool(np.all(np.isfinite(coords))) if coords.size else False
    span_nm = None
    min_heavy_heavy = None
    min_h_heavy = None
    contacts: list[GeometryPairDiagnostic] = []
    outliers: list[GeometryPairDiagnostic] = []
    crosslinks: list[CrosslinkBondDiagnostic] = []
    if finite and coords.ndim == 2 and coords.shape[1] == 3:
        atom_records = tuple(topology.atoms())
        identities = tuple(_topology_atom_identity(atom) for atom in atom_records)
        heavy = tuple(_is_heavy_atom(atom) for atom in atom_records)
        bonded_pairs = _topology_bond_index_pairs(topology)
        bonded_set = {tuple(sorted(pair)) for pair in bonded_pairs}
        span_nm = float(np.max(np.ptp(coords, axis=0))) if coords.size else 0.0
        min_heavy_heavy, min_h_heavy, contacts = _close_contact_diagnostics(
            coords,
            heavy,
            identities,
            bonded_set=bonded_set,
            heavy_heavy_close_nm=0.12,
            h_heavy_close_nm=0.08,
            max_pairs=50,
        )
        outliers = _bonded_distance_outliers(coords, identities, bonded_pairs, max_pairs=50)
        crosslinks = _crosslink_bond_diagnostics_from_topology(topology, coords, attachment_specs)
    return OpenMMValidationPhaseDiagnostics(
        phase=phase,
        coordinate_span_nm=span_nm,
        coordinates_are_finite=finite,
        has_nan=has_nan,
        has_inf=has_inf,
        potential_energy_kj_mol=potential_energy_kj_mol,
        max_force_kj_mol_nm=max_force_kj_mol_nm,
        min_heavy_heavy_distance_nm=min_heavy_heavy,
        min_h_heavy_distance_nm=min_h_heavy,
        close_contacts=tuple(contacts),
        bonded_distance_outliers=tuple(outliers),
        crosslink_bonds=tuple(crosslinks),
    )


def _invalid_phase_reason(
    phase: OpenMMValidationPhaseDiagnostics,
    *,
    max_span_nm: float,
) -> str | None:
    """Return the phase name when coordinates or energies first become invalid."""
    if not phase.coordinates_are_finite or phase.has_nan or phase.has_inf:
        return phase.phase
    if phase.coordinate_span_nm is not None and phase.coordinate_span_nm > max_span_nm:
        return phase.phase
    energy = phase.potential_energy_kj_mol
    if energy is not None and not math.isfinite(float(energy)):
        return phase.phase
    return None


def _first_invalid_phase(
    phases: tuple[OpenMMValidationPhaseDiagnostics, ...],
    *,
    max_span_nm: float,
) -> str | None:
    """Return the first invalid phase name from ordered diagnostics."""
    for phase in phases:
        reason = _invalid_phase_reason(phase, max_span_nm=max_span_nm)
        if reason is not None:
            return reason
    return None


def _close_contact_diagnostics(
    coords: np.ndarray,
    heavy: tuple[bool, ...],
    identities: tuple[str, ...],
    *,
    bonded_set: set[tuple[int, int]],
    heavy_heavy_close_nm: float,
    h_heavy_close_nm: float,
    max_pairs: int,
) -> tuple[float | None, float | None, list[GeometryPairDiagnostic]]:
    """Return minimum nonbonded distances and close-contact pairs."""
    min_heavy_heavy: float | None = None
    min_h_heavy: float | None = None
    contacts: list[GeometryPairDiagnostic] = []
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            if (i, j) in bonded_set:
                continue
            distance_nm = float(np.linalg.norm(coords[i] - coords[j]))
            both_heavy = heavy[i] and heavy[j]
            one_h_one_heavy = heavy[i] != heavy[j]
            if both_heavy:
                min_heavy_heavy = (
                    distance_nm if min_heavy_heavy is None else min(min_heavy_heavy, distance_nm)
                )
            if one_h_one_heavy:
                min_h_heavy = distance_nm if min_h_heavy is None else min(min_h_heavy, distance_nm)
            category = None
            if both_heavy and distance_nm < heavy_heavy_close_nm:
                category = "heavy-heavy-close-contact"
            elif one_h_one_heavy and distance_nm < h_heavy_close_nm:
                category = "h-heavy-close-contact"
            if category is not None and len(contacts) < max_pairs:
                contacts.append(_pair_diagnostic(i, j, distance_nm, identities, category))
    contacts.sort(key=lambda item: item.distance_nm)
    return min_heavy_heavy, min_h_heavy, contacts


def _bonded_distance_outliers(
    coords: np.ndarray,
    identities: tuple[str, ...],
    bonded_pairs: tuple[tuple[int, int], ...],
    *,
    max_pairs: int,
) -> list[GeometryPairDiagnostic]:
    """Return topology bond distances outside broad covalent bounds."""
    outliers: list[GeometryPairDiagnostic] = []
    for atom_i, atom_j in bonded_pairs:
        if atom_i >= len(coords) or atom_j >= len(coords):
            continue
        distance_nm = float(np.linalg.norm(coords[atom_i] - coords[atom_j]))
        if distance_nm < 0.06 or distance_nm > 0.25:
            outliers.append(
                _pair_diagnostic(atom_i, atom_j, distance_nm, identities, "bonded-distance-outlier")
            )
    outliers.sort(key=lambda item: abs(item.distance_nm - 0.15), reverse=True)
    return outliers[:max_pairs]


def _positions_to_numpy(positions: Any, unit_module: Any | None) -> np.ndarray:
    """Convert coordinate containers to a NumPy array for finite checks."""
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


def _topology_bond_index_pairs(topology: Any) -> tuple[tuple[int, int], ...]:
    """Return atom-index pairs from an OpenMM topology-like object."""
    pairs: list[tuple[int, int]] = []
    bonds = getattr(topology, "bonds", None)
    if bonds is None:
        return ()
    for bond in bonds():
        try:
            atom_i, atom_j = bond
            pairs.append((int(atom_i.index), int(atom_j.index)))
        except (AttributeError, TypeError, ValueError):
            continue
    return tuple(pairs)


def _topology_atom_identity(atom: Any) -> str:
    """Format an OpenMM atom-like object for diagnostics."""
    residue = getattr(atom, "residue", None)
    chain = getattr(residue, "chain", None)
    chain_id = str(getattr(chain, "id", "") or "").strip()
    residue_name = str(getattr(residue, "name", "") or "").strip()
    residue_id = str(getattr(residue, "id", "") or "").strip()
    atom_name = str(getattr(atom, "name", "") or "").strip()
    atom_index = getattr(atom, "index", None)
    return f"{atom_index}:{chain_id}:{residue_name}{residue_id}:{atom_name}"


def _pair_diagnostic(
    atom_i: int,
    atom_j: int,
    distance_nm: float,
    identities: tuple[str, ...],
    category: str,
) -> GeometryPairDiagnostic:
    """Build a pair diagnostic with atom identity text."""
    return GeometryPairDiagnostic(
        atom_i=atom_i,
        atom_j=atom_j,
        distance_nm=distance_nm,
        distance_angstrom=distance_nm * 10.0,
        atom_i_identity=identities[atom_i] if atom_i < len(identities) else None,
        atom_j_identity=identities[atom_j] if atom_j < len(identities) else None,
        category=category,
    )


def _state_max_force_kj_mol_nm(state: Any, openmm_unit: Any) -> float | None:
    """Return the maximum force norm from an OpenMM state when available."""
    try:
        forces = state.getForces(asNumpy=True)
    except (AttributeError, TypeError, ValueError, RuntimeError) as exc:
        LOGGER.debug("Could not collect max force diagnostic: %s", exc)
        return None
    force_array = np.asarray(
        forces.value_in_unit(openmm_unit.kilojoule_per_mole / openmm_unit.nanometer),
        dtype=float,
    )
    if force_array.size == 0 or not np.all(np.isfinite(force_array)):
        return None
    return float(np.max(np.linalg.norm(force_array, axis=1)))
