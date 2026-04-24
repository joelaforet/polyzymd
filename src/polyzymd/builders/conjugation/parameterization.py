"""OpenFF Interchange parameterization helpers for conjugation POC workflows."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class InterchangeParameterizationSettings(BaseModel):
    """Settings for OpenFF ForceField-to-Interchange creation."""

    force_field_names: tuple[str, ...] = (
        "ff14sb_off_impropers_0.0.4.offxml",
        "openff-2.0.0.offxml",
    )
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

    Returns
    -------
    InterchangeParameterizationResult
        Parameterization result containing the Interchange object on success.

    Raises
    ------
    RuntimeError
        If OpenFF imports or parameterization fail.
    """
    parameterization_settings = settings or InterchangeParameterizationSettings()
    if force_field is None:
        try:
            from openff.toolkit import ForceField
        except ImportError as exc:
            raise RuntimeError(
                "OpenFF Toolkit is required for conjugation parameterization. Run this helper in "
                "the conjugation-py312 or cuda-12-4 pixi environment."
            ) from exc
        force_field = ForceField(*parameterization_settings.force_field_names)

    try:
        interchange = force_field.create_interchange(topology)
    except Exception as exc:  # noqa: BLE001 - third-party parameterization errors need context
        raise RuntimeError(
            "OpenFF Interchange parameterization failed for the Pablo topology. Noncanonical "
            "product residue support may require user-supplied custom residue definitions, "
            f"charges, or force field coverage. Original error: {exc}"
        ) from exc

    return InterchangeParameterizationResult(
        success=True,
        interchange=interchange,
        force_field_names=parameterization_settings.force_field_names,
        topology_type=type(topology).__name__,
        diagnostics=("OpenFF Interchange was created from the Pablo topology",),
    )
