"""
Simulation runner for executing MD simulations with OpenMM.

This module handles running equilibration and production phases
with configurable parameters, reporters, and checkpoint management.
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Literal, Optional, Tuple, Union

import openmm
from openmm import unit as omm_unit
from openmm import XmlSerializer
from openmm.app import DCDReporter, PDBFile, StateDataReporter, Simulation

if TYPE_CHECKING:
    from polyzymd.config.schema import SimulationConfig, SimulationPhaseConfig
    from polyzymd.core.parameters import SimulationParameters

LOGGER = logging.getLogger(__name__)

# Phase types
PhaseType = Literal["equilibration", "production"]


class SimulationRunner:
    """Runner for executing OpenMM molecular dynamics simulations.

    This class manages:
    - Equilibration and production phase execution
    - Checkpoint saving and state data reporting
    - Energy minimization
    - Unique force group assignment for energy decomposition

    Example:
        >>> runner = SimulationRunner(
        ...     topology=omm_topology,
        ...     system=omm_system,
        ...     positions=omm_positions,
        ...     working_dir="output/",
        ... )
        >>> runner.minimize()
        >>> runner.run_equilibration(temperature=300, duration_ns=0.5)
        >>> runner.run_production(temperature=300, duration_ns=100)
    """

    def __init__(
        self,
        topology: Any,
        system: openmm.System,
        positions: Any,
        working_dir: Union[str, Path],
        platform: str = "CUDA",
    ) -> None:
        """Initialize the SimulationRunner.

        Args:
            topology: OpenMM Topology.
            system: OpenMM System.
            positions: Initial positions with units.
            working_dir: Working directory for output files.
            platform: Compute platform (CUDA, OpenCL, CPU).
        """
        self._topology = topology
        self._system = system
        self._positions = positions
        self._working_dir = Path(working_dir)
        self._platform_name = platform

        self._simulation: Optional[Simulation] = None
        self._current_positions = positions
        self._history: Dict[str, Any] = {}

        # Ensure working directory exists
        self._working_dir.mkdir(parents=True, exist_ok=True)

        # Apply unique force groups for energy decomposition
        self._impose_unique_force_groups()

    @property
    def simulation(self) -> Optional[Simulation]:
        """Get the current OpenMM Simulation object."""
        return self._simulation

    @property
    def working_dir(self) -> Path:
        """Get the working directory path."""
        return self._working_dir

    @property
    def history(self) -> Dict[str, Any]:
        """Get the simulation history."""
        return self._history

    def _impose_unique_force_groups(self) -> None:
        """Assign unique force groups to each force for energy decomposition."""
        from polymerist.mdtools.openmmtools.forcegroups import impose_unique_force_groups

        impose_unique_force_groups(self._system)
        LOGGER.debug("Assigned unique force groups")

    def _get_platform(self) -> openmm.Platform:
        """Get the compute platform.

        Returns:
            OpenMM Platform object.
        """
        try:
            platform = openmm.Platform.getPlatformByName(self._platform_name)
            LOGGER.info(f"Using {self._platform_name} platform")
        except Exception:
            LOGGER.warning(f"Platform {self._platform_name} not available, falling back to CPU")
            platform = openmm.Platform.getPlatformByName("CPU")
        return platform

    def _create_integrator(
        self,
        temperature: float,
        friction: float = 1.0,
        timestep: float = 2.0,
        thermostat: str = "LangevinMiddle",
    ) -> openmm.Integrator:
        """Create an integrator for the simulation.

        Args:
            temperature: Temperature in Kelvin.
            friction: Friction coefficient in 1/ps.
            timestep: Time step in femtoseconds.
            thermostat: Thermostat type.

        Returns:
            OpenMM Integrator.
        """
        temp = temperature * omm_unit.kelvin
        fric = friction / omm_unit.picosecond
        dt = timestep * omm_unit.femtosecond

        if thermostat == "LangevinMiddle":
            return openmm.LangevinMiddleIntegrator(temp, fric, dt)
        elif thermostat == "Langevin":
            return openmm.LangevinIntegrator(temp, fric, dt)
        else:
            LOGGER.warning(f"Unknown thermostat {thermostat}, using LangevinMiddle")
            return openmm.LangevinMiddleIntegrator(temp, fric, dt)

    def _add_barostat(
        self,
        pressure: float = 1.0,
        temperature: float = 300.0,
        frequency: int = 25,
    ) -> None:
        """Add a Monte Carlo barostat to the system.

        Args:
            pressure: Pressure in atmospheres.
            temperature: Temperature in Kelvin.
            frequency: Update frequency in steps.
        """
        barostat = openmm.MonteCarloBarostat(
            pressure * omm_unit.atmosphere,
            temperature * omm_unit.kelvin,
            frequency,
        )
        self._system.addForce(barostat)
        LOGGER.info(f"Added MC barostat: {pressure} atm, {temperature} K")

    def _remove_barostat(self) -> None:
        """Remove any barostat from the system."""
        forces_to_remove = []
        for i in range(self._system.getNumForces()):
            force = self._system.getForce(i)
            if isinstance(force, openmm.MonteCarloBarostat):
                forces_to_remove.append(i)

        for i in reversed(forces_to_remove):
            self._system.removeForce(i)
            LOGGER.debug("Removed barostat")

    def minimize(
        self,
        max_iterations: int = 1000,
        tolerance: float = 10.0,
    ) -> float:
        """Run energy minimization.

        Args:
            max_iterations: Maximum iterations (0 = until convergence).
            tolerance: Energy tolerance in kJ/mol/nm.

        Returns:
            Final potential energy in kJ/mol.
        """
        LOGGER.info("Running energy minimization")

        # Create temporary simulation for minimization
        integrator = openmm.VerletIntegrator(1.0 * omm_unit.femtosecond)
        platform = self._get_platform()

        simulation = Simulation(self._topology, self._system, integrator, platform)
        simulation.context.setPositions(self._current_positions)

        # Minimize
        simulation.minimizeEnergy(
            tolerance=tolerance * omm_unit.kilojoule_per_mole / omm_unit.nanometer,
            maxIterations=max_iterations,
        )

        # Get final state
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        energy = state.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
        self._current_positions = state.getPositions()

        LOGGER.info(f"Minimization complete: E = {energy:.2f} kJ/mol")

        return energy

    def run_equilibration(
        self,
        temperature: float,
        duration_ns: float,
        num_samples: int = 10,
        timestep_fs: float = 2.0,
        friction: float = 1.0,
        output_prefix: str = "equilibration",
    ) -> Dict[str, Any]:
        """Run NVT equilibration.

        Args:
            temperature: Temperature in Kelvin.
            duration_ns: Duration in nanoseconds.
            num_samples: Number of trajectory frames to save.
            timestep_fs: Time step in femtoseconds.
            friction: Friction coefficient in 1/ps.
            output_prefix: Prefix for output files.

        Returns:
            Dictionary with phase results.
        """
        LOGGER.info(f"Starting equilibration: {duration_ns} ns at {temperature} K (NVT)")

        # Remove any barostat for NVT
        self._remove_barostat()

        # Create output directory
        phase_dir = self._working_dir / output_prefix
        phase_dir.mkdir(exist_ok=True)

        # Calculate steps
        total_steps = int(duration_ns * 1e6 / timestep_fs)
        report_interval = max(1, total_steps // num_samples)

        # Create integrator and simulation
        integrator = self._create_integrator(
            temperature=temperature,
            friction=friction,
            timestep=timestep_fs,
        )
        platform = self._get_platform()

        self._simulation = Simulation(self._topology, self._system, integrator, platform)
        self._simulation.context.setPositions(self._current_positions)
        self._simulation.context.setVelocitiesToTemperature(temperature * omm_unit.kelvin)

        # Add reporters
        traj_path = phase_dir / f"{output_prefix}_trajectory.dcd"
        state_path = phase_dir / f"{output_prefix}_state_data.csv"
        pdb_path = phase_dir / f"{output_prefix}_topology.pdb"

        self._simulation.reporters.append(DCDReporter(str(traj_path), report_interval))
        self._simulation.reporters.append(
            StateDataReporter(
                str(state_path),
                report_interval,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                temperature=True,
                volume=True,
                density=True,
                speed=True,
            )
        )

        # Save topology
        with open(pdb_path, "w") as f:
            PDBFile.writeFile(
                self._topology,
                self._current_positions,
                f,
            )

        # Run simulation
        LOGGER.info(f"Running {total_steps} steps...")
        self._simulation.step(total_steps)

        # Get final state
        state = self._simulation.context.getState(
            getPositions=True, getVelocities=True, getEnergy=True
        )
        self._current_positions = state.getPositions()

        # Save checkpoint
        checkpoint_path = phase_dir / f"{output_prefix}_checkpoint.chk"
        self._simulation.saveCheckpoint(str(checkpoint_path))

        results = {
            "phase": "equilibration",
            "ensemble": "NVT",
            "temperature_K": temperature,
            "duration_ns": duration_ns,
            "total_steps": total_steps,
            "final_energy_kJ_mol": state.getPotentialEnergy().value_in_unit(
                omm_unit.kilojoule_per_mole
            ),
            "trajectory_path": str(traj_path),
            "checkpoint_path": str(checkpoint_path),
        }

        self._history["equilibration"] = results
        LOGGER.info("Equilibration complete")

        return results

    def run_production(
        self,
        temperature: float,
        duration_ns: float,
        num_samples: int = 2500,
        timestep_fs: float = 2.0,
        friction: float = 1.0,
        pressure: float = 1.0,
        barostat_frequency: int = 25,
        output_prefix: str = "production",
        segment_index: int = 0,
    ) -> Dict[str, Any]:
        """Run NPT production simulation.

        Args:
            temperature: Temperature in Kelvin.
            duration_ns: Duration in nanoseconds.
            num_samples: Number of trajectory frames to save.
            timestep_fs: Time step in femtoseconds.
            friction: Friction coefficient in 1/ps.
            pressure: Pressure in atmospheres.
            barostat_frequency: Barostat update frequency.
            output_prefix: Prefix for output files.
            segment_index: Segment index for daisy-chaining.

        Returns:
            Dictionary with phase results.
        """
        LOGGER.info(
            f"Starting production: {duration_ns} ns at {temperature} K, {pressure} atm (NPT)"
        )

        # Add barostat for NPT
        self._remove_barostat()
        self._add_barostat(
            pressure=pressure,
            temperature=temperature,
            frequency=barostat_frequency,
        )

        # Create output directory
        phase_name = f"{output_prefix}_{segment_index}"
        phase_dir = self._working_dir / phase_name
        phase_dir.mkdir(exist_ok=True)

        # Calculate steps
        total_steps = int(duration_ns * 1e6 / timestep_fs)
        report_interval = max(1, total_steps // num_samples)

        # Create integrator and simulation
        integrator = self._create_integrator(
            temperature=temperature,
            friction=friction,
            timestep=timestep_fs,
        )
        platform = self._get_platform()

        self._simulation = Simulation(self._topology, self._system, integrator, platform)
        self._simulation.context.setPositions(self._current_positions)

        # Set velocities if first segment, otherwise they'll be loaded from checkpoint
        if segment_index == 0 and "equilibration" not in self._history:
            self._simulation.context.setVelocitiesToTemperature(temperature * omm_unit.kelvin)

        # Add reporters
        traj_path = phase_dir / f"{phase_name}_trajectory.dcd"
        state_path = phase_dir / f"{phase_name}_state_data.csv"
        pdb_path = phase_dir / f"{phase_name}_topology.pdb"

        self._simulation.reporters.append(DCDReporter(str(traj_path), report_interval))
        self._simulation.reporters.append(
            StateDataReporter(
                str(state_path),
                report_interval,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=True,
                temperature=True,
                volume=True,
                density=True,
                speed=True,
            )
        )

        # Save topology
        with open(pdb_path, "w") as f:
            PDBFile.writeFile(
                self._topology,
                self._current_positions,
                f,
            )

        # Run simulation
        LOGGER.info(f"Running {total_steps} steps...")
        self._simulation.step(total_steps)

        # Get final state
        state = self._simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            getEnergy=True,
            getForces=True,
            getParameters=True,
            enforcePeriodicBox=True,
        )
        self._current_positions = state.getPositions()

        # Save checkpoint
        checkpoint_path = phase_dir / f"{phase_name}_checkpoint.chk"
        self._simulation.saveCheckpoint(str(checkpoint_path))

        # Save state XML (needed for continuation/daisy-chain)
        state_xml_path = phase_dir / f"{phase_name}_state.xml"
        with open(state_xml_path, "w") as f:
            f.write(XmlSerializer.serialize(state))
        LOGGER.info(f"Saved state to {state_xml_path}")

        # Save system XML (needed for continuation/daisy-chain)
        system_xml_path = phase_dir / f"{phase_name}_system.xml"
        with open(system_xml_path, "w") as f:
            f.write(XmlSerializer.serialize(self._system))
        LOGGER.info(f"Saved system to {system_xml_path}")

        # Save parameters JSON (needed for continuation/daisy-chain)
        params_dict = {
            "__class__": "SimulationParameters",
            "__values__": {
                "thermo_params": {
                    "__class__": "ThermoParameters",
                    "__values__": {
                        "temperature": {
                            "__class__": "Quantity",
                            "__values__": {"value": temperature, "unit": "kelvin"},
                        },
                        "thermostat_params": {
                            "__class__": "ThermostatParameters",
                            "__values__": {
                                "temperature": {
                                    "__class__": "Quantity",
                                    "__values__": {"value": temperature, "unit": "kelvin"},
                                },
                                "timescale": {
                                    "__class__": "Quantity",
                                    "__values__": {"value": friction, "unit": "1/picosecond"},
                                },
                            },
                        },
                        "barostat_params": {
                            "__class__": "BarostatParameters",
                            "__values__": {
                                "pressure": {
                                    "__class__": "Quantity",
                                    "__values__": {"value": pressure, "unit": "atmosphere"},
                                },
                                "temperature": {
                                    "__class__": "Quantity",
                                    "__values__": {"value": temperature, "unit": "kelvin"},
                                },
                                "frequency": barostat_frequency,
                            },
                        },
                    },
                },
                "integ_params": {
                    "__class__": "IntegratorParameters",
                    "__values__": {
                        "time_step": {
                            "__class__": "Quantity",
                            "__values__": {"value": timestep_fs, "unit": "femtosecond"},
                        },
                        "total_time": {
                            "__class__": "Quantity",
                            "__values__": {"value": duration_ns, "unit": "nanosecond"},
                        },
                        "num_samples": num_samples,
                    },
                },
                "reporter_params": {
                    "__class__": "ReporterParameters",
                    "__values__": {
                        "report_interval": report_interval,
                        "report_trajectory": True,
                        "report_state_data": True,
                    },
                },
            },
        }
        params_path = phase_dir / f"{phase_name}_parameters.json"
        with open(params_path, "w") as f:
            json.dump(params_dict, f, indent=2)
        LOGGER.info(f"Saved parameters to {params_path}")

        results = {
            "phase": "production",
            "segment": segment_index,
            "ensemble": "NPT",
            "temperature_K": temperature,
            "pressure_atm": pressure,
            "duration_ns": duration_ns,
            "total_steps": total_steps,
            "final_energy_kJ_mol": state.getPotentialEnergy().value_in_unit(
                omm_unit.kilojoule_per_mole
            ),
            "trajectory_path": str(traj_path),
            "checkpoint_path": str(checkpoint_path),
        }

        self._history[phase_name] = results
        LOGGER.info(f"Production segment {segment_index} complete")

        return results

    def run_from_config(
        self,
        config: "SimulationConfig",
        segment_index: int = 0,
    ) -> Dict[str, Any]:
        """Run simulation phases from a configuration.

        Args:
            config: SimulationConfig with phase settings.
            segment_index: Segment index for daisy-chaining.

        Returns:
            Combined results dictionary.
        """
        results = {}

        # Only run equilibration on first segment
        if segment_index == 0:
            eq_config = config.simulation_phases.equilibration
            eq_result = self.run_equilibration(
                temperature=config.thermodynamics.temperature,
                duration_ns=eq_config.duration,
                num_samples=eq_config.samples,
                timestep_fs=eq_config.time_step,
            )
            results["equilibration"] = eq_result

        # Run production
        prod_config = config.simulation_phases.production
        prod_result = self.run_production(
            temperature=config.thermodynamics.temperature,
            duration_ns=prod_config.duration,
            num_samples=prod_config.samples,
            timestep_fs=prod_config.time_step,
            pressure=config.thermodynamics.pressure,
            barostat_frequency=prod_config.barostat_frequency,
            segment_index=segment_index,
        )
        results["production"] = prod_result

        return results

    def save_history(self, path: Optional[Union[str, Path]] = None) -> None:
        """Save simulation history to JSON.

        Args:
            path: Output path (defaults to working_dir/simulation_history.json).
        """
        if path is None:
            path = self._working_dir / "simulation_history.json"
        else:
            path = Path(path)

        with open(path, "w") as f:
            json.dump(self._history, f, indent=2)

        LOGGER.info(f"Saved simulation history to {path}")

    def load_checkpoint(self, checkpoint_path: Union[str, Path]) -> None:
        """Load state from a checkpoint file.

        Args:
            checkpoint_path: Path to checkpoint file.
        """
        checkpoint_path = Path(checkpoint_path)

        if not checkpoint_path.exists():
            raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")

        if self._simulation is None:
            raise RuntimeError(
                "No active simulation. Create a simulation first with "
                "run_equilibration or run_production."
            )

        self._simulation.loadCheckpoint(str(checkpoint_path))

        # Update current positions
        state = self._simulation.context.getState(getPositions=True)
        self._current_positions = state.getPositions()

        LOGGER.info(f"Loaded checkpoint from {checkpoint_path}")
