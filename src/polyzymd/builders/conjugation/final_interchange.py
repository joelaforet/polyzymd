"""Final Interchange creation for solvated conjugated systems."""

from __future__ import annotations

import logging
from typing import Any

from polyzymd.builders.conjugation.pablo.charge_templates import (
    build_conjugate_charge_templates,
)
from polyzymd.builders.conjugation.pablo.parameterization import (
    InterchangeParameterizationSettings,
    create_interchange_from_pablo_topology,
    deduplicate_charge_templates,
)
from polyzymd.builders.conjugation.pablo.product_state import product_residue_names

LOGGER = logging.getLogger(__name__)

_REFUSAL_MESSAGE = (
    "PolyzyMD refuses whole-conjugate AM1-BCC charge assignment for a polymer-protein "
    "conjugate. Final conjugated Interchange creation requires product-state Pablo "
    "residue/template provenance from the successful conjugate construction path. "
    "Formal-charge templates are not validated/explicit partial charges."
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
        Product-state Pablo library created during conjugate construction. The
        library must provide residue definitions and validated/explicit
        partial-charge provenance for the final conjugate molecule template.
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

    residue_names = _validate_product_state_provenance(product_state_pablo_library)
    conjugate_templates = build_conjugate_charge_templates(topology, product_state_pablo_library)
    standard_templates = _standard_charge_templates(builder)
    charge_templates = deduplicate_charge_templates((*conjugate_templates, *standard_templates))

    LOGGER.info(
        "Creating final conjugated Interchange with product-state Pablo provenance for "
        "%d residue name(s), %d conjugate charge template(s), %d conjugate atom(s), "
        "and %d standard charge template(s)",
        len(residue_names),
        len(conjugate_templates),
        sum(
            int(getattr(template, "n_atoms", len(tuple(getattr(template, "atoms", ()) or ()))))
            for template in conjugate_templates
        ),
        len(standard_templates),
    )
    parameterize = parameterizer or create_interchange_from_pablo_topology
    result = parameterize(
        topology,
        settings=settings,
        charge_from_molecules=charge_templates,
        require_charge_templates=True,
    )
    if not getattr(result, "success", False) or getattr(result, "interchange", None) is None:
        raise RuntimeError(
            "Final conjugated Interchange parameterization did not produce an interchange"
        )
    builder._interchange = result.interchange
    return result.interchange


def _validate_product_state_provenance(product_state_pablo_library: Any) -> tuple[str, ...]:
    """Validate that product-state Pablo residue provenance is available.

    Returns
    -------
    tuple of str
        Product-state residue names covered by the Pablo library.
    """
    if product_state_pablo_library is None:
        raise RuntimeError(_REFUSAL_MESSAGE)

    if getattr(product_state_pablo_library, "residue_library", None) is None:
        raise RuntimeError(f"{_REFUSAL_MESSAGE} Product-state Pablo residue library is missing.")

    definitions = tuple(getattr(product_state_pablo_library, "definitions", ()) or ())
    if not definitions:
        raise RuntimeError(
            f"{_REFUSAL_MESSAGE} Product-state Pablo residue definitions are missing."
        )

    residue_names = _product_residue_names(product_state_pablo_library, definitions)
    if not residue_names:
        raise RuntimeError(f"{_REFUSAL_MESSAGE} Product-state residue names are missing.")
    return residue_names


def _standard_charge_templates(builder: Any) -> tuple[Any, ...]:
    """Collect standard non-conjugate charge templates from a builder."""
    collector = getattr(builder, "collect_standard_charge_templates", None)
    if callable(collector):
        return deduplicate_charge_templates(tuple(collector()))
    return ()


def _product_residue_names(
    product_state_pablo_library: Any, definitions: tuple[Any, ...]
) -> tuple[str, ...]:
    """Resolve product residue names from summaries or definitions."""
    return product_residue_names(product_state_pablo_library, definitions)
