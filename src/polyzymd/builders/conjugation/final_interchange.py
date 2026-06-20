"""Final Interchange creation for solvated conjugated systems."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from typing import Any

from polyzymd.builders.conjugation.pablo.parameterization import (
    InterchangeParameterizationSettings,
    create_interchange_from_pablo_topology,
    deduplicate_charge_templates,
)

LOGGER = logging.getLogger(__name__)

_REFUSAL_MESSAGE = (
    "PolyzyMD refuses whole-conjugate AM1-BCC charge assignment for a polymer-protein "
    "conjugate. Final conjugated Interchange creation requires product-state Pablo "
    "residue/template provenance from the successful conjugate construction path."
)


def create_final_conjugated_interchange(
    builder: Any,
    *,
    product_state_pablo_library: Any,
    settings: InterchangeParameterizationSettings | None = None,
    parameterizer: Any | None = None,
) -> Any:
    """Create final solvated Interchange for a conjugated system.

    Parameters
    ----------
    builder : Any
        Prepared ``SystemBuilder`` with a solvated topology.
    product_state_pablo_library : Any
        Product-state Pablo library created during conjugate construction.
    settings : InterchangeParameterizationSettings or None, optional
        Conjugation parameterization settings, by default ``None``.
    parameterizer : Any or None, optional
        Optional parameterization function for tests, by default ``None``.

    Returns
    -------
    Any
        OpenFF Interchange object.

    Raises
    ------
    RuntimeError
        If product-state provenance is missing or parameterization fails.
    """
    topology = getattr(builder, "_solvated_topology", None)
    if topology is None:
        raise RuntimeError("System must be solvated before creating final conjugated Interchange")

    _validate_product_state_provenance(product_state_pablo_library)
    product_templates = _product_state_charge_templates(product_state_pablo_library)
    standard_templates = _standard_charge_templates(builder)
    charge_templates = deduplicate_charge_templates((*product_templates, *standard_templates))

    LOGGER.info(
        "Creating final conjugated Interchange with %d product-state and %d standard "
        "charge template(s)",
        len(product_templates),
        len(standard_templates),
    )
    parameterize = parameterizer or create_interchange_from_pablo_topology
    result = parameterize(
        topology,
        settings=settings,
        charge_from_molecules=charge_templates,
    )
    if not getattr(result, "success", False) or getattr(result, "interchange", None) is None:
        raise RuntimeError(
            "Final conjugated Interchange parameterization did not produce an interchange"
        )
    builder._interchange = result.interchange
    return result.interchange


def _validate_product_state_provenance(product_state_pablo_library: Any) -> None:
    """Validate that product-state Pablo residue provenance is available."""
    if product_state_pablo_library is None:
        raise RuntimeError(_REFUSAL_MESSAGE)

    definitions = tuple(getattr(product_state_pablo_library, "definitions", ()) or ())
    if not definitions:
        raise RuntimeError(_REFUSAL_MESSAGE)

    residue_names = _product_residue_names(product_state_pablo_library, definitions)
    if not residue_names:
        raise RuntimeError(_REFUSAL_MESSAGE)


def _product_state_charge_templates(product_state_pablo_library: Any) -> tuple[Any, ...]:
    """Return explicit charged product templates when the library provides them."""
    templates = _first_present_sequence(
        product_state_pablo_library,
        (
            "charge_from_molecules",
            "charge_templates",
            "charged_templates",
            "template_molecules",
            "molecules",
        ),
    )
    if not templates:
        return ()
    try:
        return deduplicate_charge_templates(templates)
    except ValueError as exc:
        raise RuntimeError(f"{_REFUSAL_MESSAGE} Uncharged conjugate template: {exc}") from exc


def _standard_charge_templates(builder: Any) -> tuple[Any, ...]:
    """Collect standard non-conjugate charge templates from a builder."""
    collector = getattr(builder, "collect_standard_charge_templates", None)
    if callable(collector):
        return tuple(collector())
    return ()


def _first_present_sequence(source: Any, names: Sequence[str]) -> tuple[Any, ...]:
    """Return the first non-empty sequence attribute from a source object."""
    for name in names:
        value = getattr(source, name, None)
        if value:
            return tuple(value)
    return ()


def _product_residue_names(
    product_state_pablo_library: Any, definitions: tuple[Any, ...]
) -> tuple[str, ...]:
    """Resolve product residue names from summaries or definitions."""
    names = tuple(
        str(name) for name in getattr(product_state_pablo_library, "residue_names", ()) or ()
    )
    if names:
        return names
    definition_names = []
    for definition in definitions:
        name = str(getattr(definition, "residue_name", "")).strip()
        if name:
            definition_names.append(name)
    return tuple(definition_names)
