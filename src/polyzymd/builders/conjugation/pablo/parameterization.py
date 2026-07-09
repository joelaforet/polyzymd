"""OpenFF Interchange parameterization helpers for conjugation POC workflows."""

from __future__ import annotations

import logging
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation._linkage import parse_pdb_atom_records

LOGGER = logging.getLogger(__name__)

_OPENFF_NONBONDED_LOGGER = "openff.interchange.smirnoff._nonbonded"

_POSITION_CONVERSION_ERRORS = (AttributeError, TypeError, ValueError)

DEFAULT_CONJUGATION_FORCE_FIELD_NAMES: tuple[str, ...] = (
    "ff14sb_off_impropers_0.0.4.offxml",
    "openff-2.0.0.offxml",
)


class InterchangeParameterizationSettings(BaseModel):
    """Settings for OpenFF ForceField-to-Interchange creation."""

    force_field_names: tuple[str, ...] = DEFAULT_CONJUGATION_FORCE_FIELD_NAMES
    allow_partial_bond_orders_from_molecules: bool = True


class InterchangeParameterizationResult(BaseModel):
    """Parameterization result with the heavy Interchange object excluded."""

    model_config = {"arbitrary_types_allowed": True}

    success: bool
    interchange: Any | None = Field(default=None, exclude=True)
    force_field_names: tuple[str, ...]
    topology_type: str
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)


def create_interchange_from_pablo_topology(
    topology: Any,
    *,
    settings: InterchangeParameterizationSettings | None = None,
    force_field: Any | None = None,
    charge_from_molecules: Sequence[Any] | None = None,
    require_charge_templates: bool = False,
) -> InterchangeParameterizationResult:
    """Create an OpenFF Interchange from a Pablo topology.

    Parameters
    ----------
    topology : Any
        OpenFF topology returned by Pablo.
    settings : InterchangeParameterizationSettings or None, optional
        Force field settings, by default ``None``.
    force_field : Any or None, optional
        Optional injected force field object for tests, by default ``None``.
    charge_from_molecules : sequence of Any or None, optional
        Charged OpenFF molecule templates for ordinary non-conjugate components
        forwarded to Interchange, by default ``None``.
    require_charge_templates : bool, optional
        When ``True``, reject empty charge templates before OpenFF can assign
        implicit charges. Pablo residue-template parameterization for a
        covalent conjugate can set this to ``False`` because product-state
        residue charges come from Pablo residue definitions, by default
        ``False``.

    Returns
    -------
    InterchangeParameterizationResult
        Parameterization result containing the Interchange object on success.

    Raises
    ------
    RuntimeError
        If OpenFF imports or parameterization fail.
    """
    return create_interchange_from_openff_topology(
        topology,
        settings=settings,
        force_field=force_field,
        charge_from_molecules=charge_from_molecules,
        require_charge_templates=require_charge_templates,
        success_diagnostic="OpenFF Interchange was created from the Pablo topology",
        failure_subject="Pablo topology",
    )


def load_combined_smirnoff_force_field(
    settings: InterchangeParameterizationSettings | None = None,
    *,
    force_field_cls: Any | None = None,
) -> Any:
    """Load the combined SMIRNOFF force field stack for conjugates.

    Parameters
    ----------
    settings : InterchangeParameterizationSettings or None, optional
        Force field settings, by default ``None``.
    force_field_cls : Any or None, optional
        Optional ForceField class injection for tests, by default ``None``.

    Returns
    -------
    Any
        OpenFF ``ForceField`` instance loaded with ff14SB and Sage defaults.

    Raises
    ------
    RuntimeError
        If OpenFF Toolkit is unavailable.
    """
    parameterization_settings = settings or InterchangeParameterizationSettings()
    if force_field_cls is None:
        try:
            from openff.toolkit import ForceField
        except ImportError as exc:
            raise RuntimeError(
                "OpenFF Toolkit is required for conjugation parameterization. Run this helper in "
                "the build or a simulation pixi environment."
            ) from exc
        force_field_cls = ForceField

    return force_field_cls(*parameterization_settings.force_field_names)


def create_interchange_from_openff_topology(
    topology: Any,
    *,
    settings: InterchangeParameterizationSettings | None = None,
    force_field: Any | None = None,
    charge_from_molecules: Sequence[Any] | None = None,
    require_charge_templates: bool = True,
    success_diagnostic: str = "OpenFF Interchange was created from the OpenFF topology",
    failure_subject: str = "OpenFF topology",
) -> InterchangeParameterizationResult:
    """Create an Interchange from an OpenFF topology with explicit charges.

    Parameters
    ----------
    topology : Any
        OpenFF topology to parameterize.
    settings : InterchangeParameterizationSettings or None, optional
        Force field settings, by default ``None``.
    force_field : Any or None, optional
        Optional injected force field object for tests, by default ``None``.
    charge_from_molecules : sequence of Any or None, optional
        Charged OpenFF molecule templates forwarded via
        ``charge_from_molecules``, by default ``None``.
    require_charge_templates : bool, optional
        When ``True``, reject empty charge templates to avoid accidental
        full-protein AM1-BCC charge assignment, by default ``True``.
    success_diagnostic : str, optional
        Diagnostic text for successful parameterization.
    failure_subject : str, optional
        Human-readable topology label for error messages.

    Returns
    -------
    InterchangeParameterizationResult
        Parameterization result containing the Interchange object on success.

    Raises
    ------
    RuntimeError
        If Interchange parameterization fails.
    ValueError
        If required charge templates are missing or incomplete.
    """
    parameterization_settings = settings or InterchangeParameterizationSettings()
    if force_field is None:
        force_field = load_combined_smirnoff_force_field(parameterization_settings)

    charge_templates = deduplicate_charge_templates(charge_from_molecules or ())
    if require_charge_templates and not charge_templates:
        raise ValueError(
            "Conjugation Interchange creation requires explicit charged templates for the "
            "linked conjugate. Refusing to let OpenFF fall back to expensive full-protein "
            "charge assignment."
        )

    try:
        kwargs: dict[str, Any] = {}
        if charge_templates:
            kwargs["charge_from_molecules"] = list(charge_templates)
        interchange = _create_interchange_quietly(force_field, topology, kwargs)
    except Exception as exc:  # noqa: BLE001 - third-party parameterization errors need context
        raise RuntimeError(
            f"OpenFF Interchange parameterization failed for the {failure_subject}. "
            "Noncanonical product residue support may require user-supplied custom residue "
            f"definitions, charges, or force field coverage. Original error: {exc}"
        ) from exc

    if getattr(topology, "box_vectors", None) is not None:
        interchange.box = topology.box_vectors

    return InterchangeParameterizationResult(
        success=True,
        interchange=interchange,
        force_field_names=parameterization_settings.force_field_names,
        topology_type=type(topology).__name__,
        diagnostics=(success_diagnostic,),
    )


def deduplicate_charge_templates(charge_templates: Sequence[Any]) -> tuple[Any, ...]:
    """Validate and deduplicate charged template molecules.

    Parameters
    ----------
    charge_templates : sequence of Any
        OpenFF molecule-like objects with ``partial_charges`` populated.

    Returns
    -------
    tuple of Any
        Deduplicated templates in first-seen order.

    Raises
    ------
    ValueError
        If any template lacks partial charges.
    """
    deduplicated: list[Any] = []
    seen_keys: set[str] = set()
    for template in charge_templates:
        if template is None:
            continue
        if getattr(template, "partial_charges", None) is None:
            raise ValueError(
                "Charge templates passed to conjugation Interchange creation must already have "
                "partial_charges; refusing to trigger implicit charge assignment."
            )
        key = _charge_template_key(template)
        if key in seen_keys:
            continue
        seen_keys.add(key)
        deduplicated.append(template)
    return tuple(deduplicated)


def set_topology_positions_from_pdb(topology: Any, pdb_path: Path | str) -> Any:
    """Set OpenFF topology positions from a same-order PDB coordinate file.

    Parameters
    ----------
    topology : Any
        OpenFF topology whose atom order matches the PDB file.
    pdb_path : pathlib.Path or str
        PDB file containing relaxed coordinates.

    Returns
    -------
    Any
        The same topology object with updated positions.

    Raises
    ------
    ValueError
        If the PDB atom count differs from the topology atom count or any
        coordinate is non-finite.
    """
    import numpy as np
    from openff.units import Quantity

    atoms = parse_pdb_atom_records(Path(pdb_path))
    topology_atom_count = int(topology.n_atoms)
    if len(atoms) != topology_atom_count:
        raise ValueError(
            "Relaxed PDB atom count does not match OpenFF topology: "
            f"PDB has {len(atoms)} atoms but topology has {topology_atom_count} atoms"
        )
    positions_angstrom = np.asarray([[atom.x, atom.y, atom.z] for atom in atoms], dtype=float)
    if not np.all(np.isfinite(positions_angstrom)):
        raise ValueError(f"Relaxed PDB contains non-finite coordinates: {pdb_path}")

    scale = _infer_pdb_coordinate_scale(topology, positions_angstrom)
    if scale != 1.0:
        LOGGER.warning(
            "Relaxed PDB coordinates in %s appear scaled relative to the OpenFF topology; "
            "applying %.6g before storing positions",
            pdb_path,
            scale,
        )
        positions_angstrom = positions_angstrom * scale

    # Store explicit nanometer quantities because OpenFF topology positions are
    # later consumed by box builders and OpenMM conversion paths in nanometers
    topology.set_positions(Quantity(positions_angstrom / 10.0, "nanometer"))
    return topology


def _infer_pdb_coordinate_scale(topology: Any, positions_angstrom: Any) -> float:
    """Infer whether PDB coordinate magnitudes need unit-scale correction.

    Parameters
    ----------
    topology : Any
        OpenFF topology before the relaxed coordinates are transferred.
    positions_angstrom : Any
        Parsed PDB coordinates, nominally in Angstrom.

    Returns
    -------
    float
        Multiplicative scale that converts the parsed coordinates to Angstrom.
    """
    import numpy as np

    reference_positions = _topology_positions_as_angstrom(topology)
    if reference_positions is None:
        return 1.0

    reference_span = np.ptp(reference_positions, axis=0)
    pdb_span = np.ptp(np.asarray(positions_angstrom, dtype=float), axis=0)
    usable_axes = (reference_span > 1.0e-6) & (pdb_span > 1.0e-6)
    if int(np.count_nonzero(usable_axes)) < 2:
        return 1.0

    span_ratios = pdb_span[usable_axes] / reference_span[usable_axes]
    median_ratio = float(np.median(span_ratios))
    if median_ratio <= 5.0:
        return 1.0

    candidates = (1.0, 0.1, 0.01, 0.001)

    def _score(scale: float) -> float:
        scaled_ratios = span_ratios * scale
        return float(np.median(np.abs(np.log10(scaled_ratios))))

    unscaled_score = _score(1.0)
    best_scale = min(candidates, key=_score)
    if best_scale == 1.0 or _score(best_scale) >= unscaled_score * 0.5:
        return 1.0

    scaled_ratios = span_ratios * best_scale
    if not bool(np.all((scaled_ratios > 0.2) & (scaled_ratios < 5.0))):
        return 1.0
    return best_scale


def _topology_positions_as_angstrom(topology: Any) -> Any | None:
    """Return current topology positions as an Angstrom array when available.

    Parameters
    ----------
    topology : Any
        OpenFF topology-like object with optional positions.

    Returns
    -------
    Any or None
        Coordinate array in Angstrom, or ``None`` if positions are unavailable.
    """
    import numpy as np

    get_positions = getattr(topology, "get_positions", None)
    if not callable(get_positions):
        return None
    try:
        positions = get_positions()
    except _POSITION_CONVERSION_ERRORS as exc:
        LOGGER.warning(
            "Could not read topology positions from %s; coordinate-scale inference will be "
            "skipped: %s",
            type(topology).__name__,
            exc,
        )
        return None
    if positions is None:
        return None
    if hasattr(positions, "m_as"):
        try:
            return np.asarray(positions.m_as("angstrom"), dtype=float)
        except _POSITION_CONVERSION_ERRORS as exc:
            LOGGER.warning(
                "Could not convert topology positions of type %s to Angstrom; "
                "coordinate-scale inference will be skipped: %s",
                type(positions).__name__,
                exc,
            )
            return None
    try:
        return np.asarray(positions, dtype=float)
    except _POSITION_CONVERSION_ERRORS as exc:
        LOGGER.warning(
            "Could not coerce topology positions of type %s to an Angstrom array; "
            "coordinate-scale inference will be skipped: %s",
            type(positions).__name__,
            exc,
        )
        return None


def _create_interchange_quietly(force_field: Any, topology: Any, kwargs: dict[str, Any]) -> Any:
    """Call ``create_interchange`` while suppressing noisy charge logs."""
    with suppress_openff_library_charge_info():
        return force_field.create_interchange(topology, **kwargs)


@contextmanager
def suppress_openff_library_charge_info() -> Iterator[None]:
    """Temporarily suppress per-atom OpenFF LibraryCharges INFO records.

    OpenFF Interchange emits one INFO record per atom from its SMIRNOFF
    nonbonded handler while applying library charges. The scoped override keeps
    OpenFF warnings and errors visible, restores the exact previous logger
    level, and does not alter PolyzyMD loggers.

    Yields
    ------
    None
        Control while the OpenFF nonbonded logger is raised to ``WARNING``.
    """
    nonbonded_logger = logging.getLogger(_OPENFF_NONBONDED_LOGGER)
    previous_level = nonbonded_logger.level
    nonbonded_logger.setLevel(logging.WARNING)
    try:
        yield
    finally:
        nonbonded_logger.setLevel(previous_level)


def _charge_template_key(template: Any) -> str:
    """Return a stable deduplication key for a charge template."""
    to_smiles = getattr(template, "to_smiles", None)
    if callable(to_smiles):
        try:
            return str(to_smiles(mapped=True))
        except TypeError:
            return str(to_smiles())
    return f"{type(template).__name__}:{id(template)}"


def _formal_charge_value(formal_charge: Any) -> float:
    """Return a formal charge value in elementary charge units."""
    conversion_error: Exception | None = None
    if hasattr(formal_charge, "m_as"):
        try:
            return float(formal_charge.m_as("elementary_charge"))
        except _formal_charge_conversion_errors() as exc:
            conversion_error = exc
    if hasattr(formal_charge, "value_in_unit"):
        try:
            from openff.units.openmm import to_openmm

            openmm_charge = to_openmm(formal_charge)
            from openmm import unit as openmm_unit

            return float(openmm_charge.value_in_unit(openmm_unit.elementary_charge))
        except _formal_charge_conversion_errors() as exc:
            conversion_error = exc
    try:
        scalar_charge = float(formal_charge)
    except (TypeError, ValueError) as scalar_exc:
        if conversion_error is not None:
            raise ValueError(
                "Could not convert formal charge to elementary-charge units via unit-aware "
                "conversion or scalar fallback. Unit conversion error: "
                f"{conversion_error}. Scalar fallback error: {scalar_exc}"
            ) from scalar_exc
        raise ValueError(
            "Could not convert formal charge to elementary-charge units via scalar fallback: "
            f"{scalar_exc}"
        ) from scalar_exc
    if conversion_error is not None:
        LOGGER.warning(
            "Falling back to scalar formal charge conversion for type %s after unit-aware "
            "conversion failed: %s",
            type(formal_charge).__name__,
            conversion_error,
        )
    return scalar_charge


def _formal_charge_conversion_errors() -> tuple[type[BaseException], ...]:
    """Return expected formal-charge conversion exception classes."""
    errors: tuple[type[BaseException], ...] = (AttributeError, TypeError, ValueError, ImportError)
    try:
        from openff.units.errors import UnitConversionError
    except ImportError:
        pass
    else:
        errors += (UnitConversionError,)
    try:
        from openmm.unit import IncompatibleUnitError
    except (ImportError, AttributeError):
        pass
    else:
        errors += (IncompatibleUnitError,)
    return errors
