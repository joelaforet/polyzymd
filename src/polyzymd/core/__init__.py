"""Core module containing fundamental data structures and utilities."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS: dict[str, tuple[str, str]] = {
    "PREDEFINED_GROUPS": ("polyzymd.core.atom_groups", "PREDEFINED_GROUPS"),
    "AtomGroupResolver": ("polyzymd.core.atom_groups", "AtomGroupResolver"),
    "SystemComponentInfo": ("polyzymd.core.atom_groups", "SystemComponentInfo"),
    "BarostatParameters": ("polyzymd.core.parameters", "BarostatParameters"),
    "IntegratorParameters": ("polyzymd.core.parameters", "IntegratorParameters"),
    "ReporterParameters": ("polyzymd.core.parameters", "ReporterParameters"),
    "SimulationParameters": ("polyzymd.core.parameters", "SimulationParameters"),
    "SimulationPhase": ("polyzymd.core.parameters", "SimulationPhase"),
    "ThermoParameters": ("polyzymd.core.parameters", "ThermoParameters"),
    "ThermostatParameters": ("polyzymd.core.parameters", "ThermostatParameters"),
    "KCAL_MOL_ANGSTROM2_TO_KJ_MOL_NM2": (
        "polyzymd.core.position_restraints",
        "KCAL_MOL_ANGSTROM2_TO_KJ_MOL_NM2",
    ),
    "PositionalRestraintForce": (
        "polyzymd.core.position_restraints",
        "PositionalRestraintForce",
    ),
    "add_position_restraints_to_system": (
        "polyzymd.core.position_restraints",
        "add_position_restraints_to_system",
    ),
    "create_position_restraints": (
        "polyzymd.core.position_restraints",
        "create_position_restraints",
    ),
    "remove_position_restraints_from_system": (
        "polyzymd.core.position_restraints",
        "remove_position_restraints_from_system",
    ),
    "AtomSelection": ("polyzymd.core.restraints", "AtomSelection"),
    "RestraintDefinition": ("polyzymd.core.restraints", "RestraintDefinition"),
    "RestraintFactory": ("polyzymd.core.restraints", "RestraintFactory"),
    "RestraintType": ("polyzymd.core.restraints", "RestraintType"),
}

__all__ = [
    "ThermostatParameters",
    "BarostatParameters",
    "ThermoParameters",
    "IntegratorParameters",
    "ReporterParameters",
    "SimulationParameters",
    "SimulationPhase",
    "RestraintType",
    "AtomSelection",
    "RestraintDefinition",
    "RestraintFactory",
    "PositionalRestraintForce",
    "create_position_restraints",
    "add_position_restraints_to_system",
    "remove_position_restraints_from_system",
    "KCAL_MOL_ANGSTROM2_TO_KJ_MOL_NM2",
    "PREDEFINED_GROUPS",
    "SystemComponentInfo",
    "AtomGroupResolver",
]


def __getattr__(name: str) -> Any:
    """Lazily import OpenMM-dependent core helpers."""
    try:
        module_name, attr_name = _EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc

    module = import_module(module_name)
    return getattr(module, attr_name)


def __dir__() -> list[str]:
    """Return module attributes for tab completion and introspection."""
    return sorted(set(globals()) | set(__all__))
