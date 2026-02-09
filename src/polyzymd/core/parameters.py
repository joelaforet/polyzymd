"""
Simulation parameter dataclasses for OpenMM integration.

This module provides dataclasses that bridge the configuration system
with OpenMM's simulation parameters, supporting serialization to JSON
for checkpoint/restart capabilities.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

from openmm import unit as openmm_unit
from openmm.unit import Quantity

# =============================================================================
# Unit Serialization Helpers
# =============================================================================


def quantity_to_dict(q: Quantity) -> Dict[str, Any]:
    """Convert an OpenMM Quantity to a serializable dictionary.

    Args:
        q: OpenMM Quantity with units

    Returns:
        Dictionary with value and unit information
    """
    # Get the unit as a string
    unit_str = str(q.unit)

    # Handle inverse units (e.g., /picosecond)
    if hasattr(q, "_value"):
        value = q._value
    else:
        value = q.value_in_unit(q.unit)

    return {"__class__": "Quantity", "__values__": {"value": value, "unit": unit_str}}


def quantity_from_dict(d: Dict[str, Any]) -> Quantity:
    """Restore an OpenMM Quantity from a serialized dictionary.

    Args:
        d: Dictionary with value and unit information

    Returns:
        OpenMM Quantity
    """
    values = d.get("__values__", d)
    value = values["value"]
    unit_str = values["unit"]

    # Map common unit strings to OpenMM units
    unit_mapping = {
        "kelvin": openmm_unit.kelvin,
        "atmosphere": openmm_unit.atmosphere,
        "atmospheres": openmm_unit.atmosphere,
        "femtosecond": openmm_unit.femtosecond,
        "femtoseconds": openmm_unit.femtoseconds,
        "picosecond": openmm_unit.picosecond,
        "picoseconds": openmm_unit.picoseconds,
        "nanosecond": openmm_unit.nanosecond,
        "nanoseconds": openmm_unit.nanoseconds,
        "nanometer": openmm_unit.nanometer,
        "nanometers": openmm_unit.nanometers,
        "angstrom": openmm_unit.angstrom,
        "angstroms": openmm_unit.angstroms,
    }

    # Handle inverse units
    if unit_str.startswith("/"):
        base_unit_str = unit_str[1:]
        if base_unit_str in unit_mapping:
            return value / unit_mapping[base_unit_str]

    if unit_str in unit_mapping:
        return value * unit_mapping[unit_str]

    # Try to get from openmm.unit module directly
    if hasattr(openmm_unit, unit_str):
        return value * getattr(openmm_unit, unit_str)

    raise ValueError(f"Unknown unit: {unit_str}")


# =============================================================================
# Parameter Dataclasses
# =============================================================================


@dataclass
class ThermostatParameters:
    """Parameters for temperature control.

    Attributes:
        temperature: Target temperature
        timescale: Coupling timescale (friction coefficient)
        thermostat: Type of thermostat to use
    """

    temperature: Quantity = field(default_factory=lambda: 300.0 * openmm_unit.kelvin)
    timescale: Quantity = field(default_factory=lambda: 1.0 / openmm_unit.picosecond)
    thermostat: str = "LangevinMiddle"

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "__class__": "ThermostatParameters",
            "__values__": {
                "temperature": quantity_to_dict(self.temperature),
                "timescale": quantity_to_dict(self.timescale),
                "thermostat": self.thermostat,
            },
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ThermostatParameters":
        """Deserialize from dictionary."""
        values = d.get("__values__", d)
        return cls(
            temperature=quantity_from_dict(values["temperature"]),
            timescale=quantity_from_dict(values["timescale"]),
            thermostat=values["thermostat"],
        )


@dataclass
class BarostatParameters:
    """Parameters for pressure control.

    Attributes:
        pressure: Target pressure
        temperature: Temperature for barostat (should match thermostat)
        update_frequency: Steps between barostat updates
        barostat: Type of barostat to use
    """

    pressure: Quantity = field(default_factory=lambda: 1.0 * openmm_unit.atmosphere)
    temperature: Quantity = field(default_factory=lambda: 300.0 * openmm_unit.kelvin)
    update_frequency: int = 25
    barostat: str = "MC"

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "__class__": "BarostatParameters",
            "__values__": {
                "pressure": quantity_to_dict(self.pressure),
                "temperature": quantity_to_dict(self.temperature),
                "update_frequency": self.update_frequency,
                "barostat": self.barostat,
            },
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "BarostatParameters":
        """Deserialize from dictionary."""
        values = d.get("__values__", d)
        return cls(
            pressure=quantity_from_dict(values["pressure"]),
            temperature=quantity_from_dict(values["temperature"]),
            update_frequency=values["update_frequency"],
            barostat=values["barostat"],
        )


@dataclass
class ThermoParameters:
    """Combined thermodynamic parameters.

    Attributes:
        thermostat_params: Temperature control parameters
        barostat_params: Pressure control parameters (None for NVT/NVE)
    """

    thermostat_params: ThermostatParameters = field(default_factory=ThermostatParameters)
    barostat_params: Optional[BarostatParameters] = None

    @property
    def ensemble(self) -> str:
        """Infer ensemble from parameters."""
        if self.barostat_params is not None:
            return "NPT"
        return "NVT"

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        result = {
            "__class__": "ThermoParameters",
            "__values__": {"thermostat_params": self.thermostat_params.to_dict()},
        }
        if self.barostat_params is not None:
            result["__values__"]["barostat_params"] = self.barostat_params.to_dict()
        return result

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ThermoParameters":
        """Deserialize from dictionary."""
        values = d.get("__values__", d)
        barostat = None
        if "barostat_params" in values:
            barostat = BarostatParameters.from_dict(values["barostat_params"])
        return cls(
            thermostat_params=ThermostatParameters.from_dict(values["thermostat_params"]),
            barostat_params=barostat,
        )


@dataclass
class IntegratorParameters:
    """Parameters for the MD integrator.

    Attributes:
        time_step: Integration time step
        total_time: Total simulation time
        num_samples: Number of frames to save
    """

    time_step: Quantity = field(default_factory=lambda: 2.0 * openmm_unit.femtosecond)
    total_time: Quantity = field(default_factory=lambda: 1.0 * openmm_unit.nanosecond)
    num_samples: int = 100

    @property
    def total_steps(self) -> int:
        """Calculate total number of integration steps."""
        return int(self.total_time / self.time_step)

    @property
    def reporting_interval(self) -> int:
        """Calculate steps between trajectory frames."""
        return max(1, self.total_steps // self.num_samples)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "__class__": "IntegratorParameters",
            "__values__": {
                "time_step": quantity_to_dict(self.time_step),
                "total_time": quantity_to_dict(self.total_time),
                "num_samples": self.num_samples,
            },
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "IntegratorParameters":
        """Deserialize from dictionary."""
        values = d.get("__values__", d)
        return cls(
            time_step=quantity_from_dict(values["time_step"]),
            total_time=quantity_from_dict(values["total_time"]),
            num_samples=values["num_samples"],
        )


# Default state data properties to report
DEFAULT_STATE_DATA_PROPS = {
    "step": True,
    "time": True,
    "potentialEnergy": True,
    "kineticEnergy": True,
    "totalEnergy": True,
    "temperature": True,
    "volume": True,
    "density": True,
    "speed": True,
    "progress": False,
    "remainingTime": False,
    "elapsedTime": False,
}


@dataclass
class ReporterParameters:
    """Parameters for simulation output/reporting.

    Attributes:
        traj_ext: Trajectory file extension
        state_data: Dictionary of state data properties to report
        report_checkpoint: Whether to save checkpoint files
        report_state: Whether to save state XML files
        report_trajectory: Whether to save trajectory
        report_state_data: Whether to save CSV state data
    """

    traj_ext: str = "dcd"
    state_data: Dict[str, bool] = field(default_factory=lambda: DEFAULT_STATE_DATA_PROPS.copy())
    report_checkpoint: bool = True
    report_state: bool = True
    report_trajectory: bool = True
    report_state_data: bool = True

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "__class__": "ReporterParameters",
            "__values__": {
                "traj_ext": self.traj_ext,
                "state_data": self.state_data,
                "report_checkpoint": self.report_checkpoint,
                "report_state": self.report_state,
                "report_trajectory": self.report_trajectory,
                "report_state_data": self.report_state_data,
            },
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ReporterParameters":
        """Deserialize from dictionary."""
        values = d.get("__values__", d)
        return cls(
            traj_ext=values.get("traj_ext", "dcd"),
            state_data=values.get("state_data", DEFAULT_STATE_DATA_PROPS.copy()),
            report_checkpoint=values.get("report_checkpoint", True),
            report_state=values.get("report_state", True),
            report_trajectory=values.get("report_trajectory", True),
            report_state_data=values.get("report_state_data", True),
        )


@dataclass
class SimulationParameters:
    """Complete simulation parameters combining all components.

    Attributes:
        thermo_params: Thermodynamic parameters
        integ_params: Integrator parameters
        reporter_params: Reporter parameters
    """

    thermo_params: ThermoParameters = field(default_factory=ThermoParameters)
    integ_params: IntegratorParameters = field(default_factory=IntegratorParameters)
    reporter_params: ReporterParameters = field(default_factory=ReporterParameters)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "__class__": "SimulationParameters",
            "__values__": {
                "thermo_params": self.thermo_params.to_dict(),
                "integ_params": self.integ_params.to_dict(),
                "reporter_params": self.reporter_params.to_dict(),
            },
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "SimulationParameters":
        """Deserialize from dictionary."""
        values = d.get("__values__", d)
        return cls(
            thermo_params=ThermoParameters.from_dict(values["thermo_params"]),
            integ_params=IntegratorParameters.from_dict(values["integ_params"]),
            reporter_params=ReporterParameters.from_dict(values["reporter_params"]),
        )

    def to_json(self, path: Union[str, Path]) -> None:
        """Save parameters to JSON file."""
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_json(cls, path: Union[str, Path]) -> "SimulationParameters":
        """Load parameters from JSON file."""
        with open(path, "r") as f:
            data = json.load(f)
        return cls.from_dict(data)


@dataclass
class SimulationPhase:
    """Represents a single simulation phase (e.g., equilibration, production).

    Attributes:
        name: Phase identifier
        parameters: Simulation parameters for this phase
    """

    name: str
    parameters: SimulationParameters

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {"name": self.name, "parameters": self.parameters.to_dict()}

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "SimulationPhase":
        """Deserialize from dictionary."""
        return cls(name=d["name"], parameters=SimulationParameters.from_dict(d["parameters"]))


def create_nvt_parameters(
    temperature: float,
    duration: float,
    samples: int,
    time_step: float = 2.0,
    thermostat_timescale: float = 1.0,
) -> SimulationParameters:
    """Factory function to create NVT simulation parameters.

    Args:
        temperature: Temperature in Kelvin
        duration: Simulation duration in nanoseconds
        samples: Number of trajectory frames to save
        time_step: Integration time step in femtoseconds
        thermostat_timescale: Thermostat timescale in picoseconds

    Returns:
        Configured SimulationParameters for NVT
    """
    return SimulationParameters(
        thermo_params=ThermoParameters(
            thermostat_params=ThermostatParameters(
                temperature=temperature * openmm_unit.kelvin,
                timescale=1.0 / (thermostat_timescale * openmm_unit.picosecond),
                thermostat="LangevinMiddle",
            )
        ),
        integ_params=IntegratorParameters(
            time_step=time_step * openmm_unit.femtosecond,
            total_time=duration * openmm_unit.nanosecond,
            num_samples=samples,
        ),
        reporter_params=ReporterParameters(),
    )


def create_npt_parameters(
    temperature: float,
    pressure: float,
    duration: float,
    samples: int,
    time_step: float = 2.0,
    thermostat_timescale: float = 1.0,
    barostat_frequency: int = 25,
) -> SimulationParameters:
    """Factory function to create NPT simulation parameters.

    Args:
        temperature: Temperature in Kelvin
        pressure: Pressure in atmospheres
        duration: Simulation duration in nanoseconds
        samples: Number of trajectory frames to save
        time_step: Integration time step in femtoseconds
        thermostat_timescale: Thermostat timescale in picoseconds
        barostat_frequency: Steps between barostat updates

    Returns:
        Configured SimulationParameters for NPT
    """
    return SimulationParameters(
        thermo_params=ThermoParameters(
            thermostat_params=ThermostatParameters(
                temperature=temperature * openmm_unit.kelvin,
                timescale=1.0 / (thermostat_timescale * openmm_unit.picosecond),
                thermostat="LangevinMiddle",
            ),
            barostat_params=BarostatParameters(
                pressure=pressure * openmm_unit.atmosphere,
                temperature=temperature * openmm_unit.kelvin,
                update_frequency=barostat_frequency,
                barostat="MC",
            ),
        ),
        integ_params=IntegratorParameters(
            time_step=time_step * openmm_unit.femtosecond,
            total_time=duration * openmm_unit.nanosecond,
            num_samples=samples,
        ),
        reporter_params=ReporterParameters(),
    )
