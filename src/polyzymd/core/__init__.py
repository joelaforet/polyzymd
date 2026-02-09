"""Core module containing fundamental data structures and utilities."""

from polyzymd.core.atom_groups import (
    PREDEFINED_GROUPS,
    AtomGroupResolver,
    SystemComponentInfo,
)
from polyzymd.core.parameters import (
    BarostatParameters,
    IntegratorParameters,
    ReporterParameters,
    SimulationParameters,
    SimulationPhase,
    ThermoParameters,
    ThermostatParameters,
)
from polyzymd.core.position_restraints import (
    KCAL_MOL_ANGSTROM2_TO_KJ_MOL_NM2,
    PositionalRestraintForce,
    add_position_restraints_to_system,
    create_position_restraints,
    remove_position_restraints_from_system,
)
from polyzymd.core.restraints import (
    AtomSelection,
    RestraintDefinition,
    RestraintFactory,
    RestraintType,
)

__all__ = [
    # Parameters
    "ThermostatParameters",
    "BarostatParameters",
    "ThermoParameters",
    "IntegratorParameters",
    "ReporterParameters",
    "SimulationParameters",
    "SimulationPhase",
    # Distance restraints
    "RestraintType",
    "AtomSelection",
    "RestraintDefinition",
    "RestraintFactory",
    # Position restraints
    "PositionalRestraintForce",
    "create_position_restraints",
    "add_position_restraints_to_system",
    "remove_position_restraints_from_system",
    "KCAL_MOL_ANGSTROM2_TO_KJ_MOL_NM2",
    # Atom groups
    "PREDEFINED_GROUPS",
    "SystemComponentInfo",
    "AtomGroupResolver",
]
