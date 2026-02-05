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
    from polyzymd.config.schema import (
        EquilibrationStageConfig,
        SimulationConfig,
        SimulationPhaseConfig,
        SimulationPhasesConfig,
    )
    from polyzymd.core.atom_groups import AtomGroupResolver
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
        self._current_box_vectors = None  # Updated during NPT stages
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
        duration_ns: Optional[float] = None,
        num_samples: int = 10,
        timestep_fs: float = 2.0,
        friction: float = 1.0,
        output_prefix: str = "equilibration",
        config: Optional["SimulationPhasesConfig"] = None,
    ) -> Dict[str, Any]:
        """Run equilibration phase.

        Supports two modes:
        1. Parameter-based (legacy): Pass duration_ns, num_samples, etc.
        2. Config-based: Pass config for automatic mode selection

        When config is provided and uses staged equilibration, position
        restraints and temperature ramping are handled automatically.
        Component information is derived from the topology's chain IDs.

        Args:
            temperature: Temperature in Kelvin
            duration_ns: Duration in nanoseconds (required for legacy mode)
            num_samples: Number of trajectory frames to save
            timestep_fs: Time step in femtoseconds
            friction: Friction coefficient in 1/ps
            output_prefix: Prefix for output files
            config: SimulationPhasesConfig for config-based dispatch

        Returns:
            Dictionary with equilibration results

        Raises:
            ValueError: If neither config nor duration_ns is provided
        """
        # Config-based dispatch (preferred path)
        if config is not None:
            if config.uses_staged_equilibration:
                # Multi-stage equilibration with position restraints
                from polyzymd.core.atom_groups import AtomGroupResolver, SystemComponentInfo

                # Log all stages upfront for reproducibility
                LOGGER.info(
                    f"Starting multi-stage equilibration with "
                    f"{len(config.equilibration_stages)} stages:"
                )
                for i, stage in enumerate(config.equilibration_stages):
                    restraint_info = (
                        ", ".join(
                            f"{r.group}@{r.force_constant:.0f}" for r in stage.position_restraints
                        )
                        or "none"
                    )
                    if stage.is_temperature_ramping:
                        temp_info = f"{stage.temperature_start}K -> {stage.temperature_end}K"
                    else:
                        temp_info = f"{stage.temperature}K"
                    LOGGER.info(
                        f"  Stage {i}: {stage.name} - {stage.duration} ns, "
                        f"{stage.ensemble.value}, {temp_info}, restraints: [{restraint_info}]"
                    )

                component_info = SystemComponentInfo.from_topology(self._topology)
                resolver = AtomGroupResolver(self._topology, component_info)

                return self.run_staged_equilibration(
                    stages=config.equilibration_stages,
                    atom_group_resolver=resolver,
                    target_temperature=temperature,
                )
            else:
                # Simple equilibration via config
                eq_config = config.equilibration
                return self._run_simple_equilibration(
                    temperature=temperature,
                    duration_ns=eq_config.duration,
                    num_samples=eq_config.samples,
                    timestep_fs=eq_config.time_step or timestep_fs,
                    friction=friction,
                    output_prefix=output_prefix,
                )

        # Legacy parameter-based mode
        if duration_ns is None:
            raise ValueError(
                "duration_ns is required when config is not provided. "
                "Either pass duration_ns or pass a SimulationPhasesConfig."
            )

        return self._run_simple_equilibration(
            temperature=temperature,
            duration_ns=duration_ns,
            num_samples=num_samples,
            timestep_fs=timestep_fs,
            friction=friction,
            output_prefix=output_prefix,
        )

    def _run_simple_equilibration(
        self,
        temperature: float,
        duration_ns: float,
        num_samples: int = 10,
        timestep_fs: float = 2.0,
        friction: float = 1.0,
        output_prefix: str = "equilibration",
    ) -> Dict[str, Any]:
        """Run simple NVT equilibration (internal implementation).

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

    def run_equilibration_stage(
        self,
        stage: "EquilibrationStageConfig",
        reference_positions: Any,
        atom_group_resolver: "AtomGroupResolver",
        stage_index: int,
        default_timestep: float = 2.0,
        default_friction: float = 1.0,
    ) -> Dict[str, Any]:
        """Run a single equilibration stage with optional position restraints.

        This method runs one stage of a multi-stage equilibration protocol.
        It supports:
        - Position restraints on predefined atom groups
        - Temperature ramping (simulated annealing)
        - NVT or NPT ensembles

        Args:
            stage: EquilibrationStageConfig with stage settings
            reference_positions: Positions to restrain atoms to (typically post-minimization)
            atom_group_resolver: Resolver for predefined atom group names
            stage_index: Index of this stage (for output naming)
            default_timestep: Default time step in fs if not specified in stage
            default_friction: Default friction coefficient in 1/ps

        Returns:
            Dictionary with stage results
        """
        from polyzymd.core.position_restraints import (
            add_position_restraints_to_system,
            remove_position_restraints_from_system,
        )
        from polyzymd.config.schema import Ensemble

        stage_name = f"equilibration_{stage_index}_{stage.name}"
        LOGGER.info(f"Starting equilibration stage: {stage.name} ({stage.duration} ns)")

        # Create output directory for this stage
        phase_dir = self._working_dir / stage_name
        phase_dir.mkdir(exist_ok=True)

        # Get stage parameters (with defaults)
        timestep_fs = stage.time_step if stage.time_step is not None else default_timestep
        friction = default_friction
        thermostat_timescale = stage.thermostat_timescale if stage.thermostat_timescale else 1.0

        # Calculate steps and reporting interval
        total_steps = int(stage.duration * 1e6 / timestep_fs)
        report_interval = max(1, total_steps // stage.samples)

        # Handle ensemble - add/remove barostat
        if stage.ensemble == Ensemble.NPT:
            self._remove_barostat()
            pressure = 1.0  # Default pressure for NPT stages
            barostat_freq = stage.barostat_frequency if stage.barostat_frequency else 25
            start_temp = stage.get_start_temperature()
            self._add_barostat(
                pressure=pressure,
                temperature=start_temp,
                frequency=barostat_freq,
            )
        else:
            # NVT - ensure no barostat
            self._remove_barostat()

        # Add position restraints
        restraint_force_indices = []
        for restraint_config in stage.position_restraints:
            atom_indices = atom_group_resolver.resolve(restraint_config.group)
            if atom_indices:
                force_idx = add_position_restraints_to_system(
                    system=self._system,
                    atom_indices=atom_indices,
                    positions=reference_positions,
                    force_constant=restraint_config.force_constant,
                )
                if force_idx >= 0:
                    restraint_force_indices.append(force_idx)
                LOGGER.info(
                    f"Added position restraints to {len(atom_indices)} atoms "
                    f"in group '{restraint_config.group}' "
                    f"(k={restraint_config.force_constant:.1f} kJ/mol/nm^2)"
                )
            else:
                LOGGER.warning(
                    f"No atoms found for group '{restraint_config.group}' - skipping restraint"
                )

        # Create integrator with starting temperature
        start_temp = stage.get_start_temperature()
        integrator = self._create_integrator(
            temperature=start_temp,
            friction=friction,
            timestep=timestep_fs,
        )
        platform = self._get_platform()

        # Create simulation
        self._simulation = Simulation(self._topology, self._system, integrator, platform)

        # Set box vectors BEFORE positions - critical for NPT stage transitions
        # where box dimensions may have changed from previous stage
        if self._current_box_vectors is not None:
            self._simulation.context.setPeriodicBoxVectors(*self._current_box_vectors)

        self._simulation.context.setPositions(self._current_positions)
        self._simulation.context.setVelocitiesToTemperature(start_temp * omm_unit.kelvin)

        # Log initial energy for diagnostics
        _state = self._simulation.context.getState(getEnergy=True)
        _energy = _state.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
        LOGGER.info(f"Stage {stage_index} ({stage_name}): initial PE = {_energy:.2f} kJ/mol")

        # Set up reporters
        traj_path = phase_dir / f"{stage_name}_trajectory.dcd"
        state_path = phase_dir / f"{stage_name}_state_data.csv"
        pdb_path = phase_dir / f"{stage_name}_topology.pdb"

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

        # Save initial topology
        with open(pdb_path, "w") as f:
            PDBFile.writeFile(
                self._topology,
                self._current_positions,
                f,
            )

        # Run simulation with temperature ramping if needed
        if stage.is_temperature_ramping:
            LOGGER.info(
                f"Temperature ramping: {stage.temperature_start} K -> {stage.temperature_end} K "
                f"(increment={stage.temperature_increment} K every {stage.temperature_interval} fs)"
            )
            current_temp = stage.temperature_start
            steps_per_update = int(stage.temperature_interval / timestep_fs)

            # Calculate total temperature updates needed
            temp_range = stage.temperature_end - stage.temperature_start
            num_updates = int(temp_range / stage.temperature_increment)
            steps_for_ramping = num_updates * steps_per_update
            remaining_steps = total_steps - steps_for_ramping

            # Temperature ramping phase
            while current_temp < stage.temperature_end:
                integrator.setTemperature(current_temp * omm_unit.kelvin)
                self._simulation.step(steps_per_update)
                current_temp += stage.temperature_increment

            # Final temperature - run remaining steps
            integrator.setTemperature(stage.temperature_end * omm_unit.kelvin)
            if remaining_steps > 0:
                LOGGER.info(
                    f"Running {remaining_steps} steps at final temperature {stage.temperature_end} K"
                )
                self._simulation.step(remaining_steps)
        else:
            # Constant temperature - just run all steps
            LOGGER.info(f"Running {total_steps} steps at {stage.temperature} K")
            self._simulation.step(total_steps)

        # Get final state (including box vectors for NPT stages)
        state = self._simulation.context.getState(
            getPositions=True, getVelocities=True, getEnergy=True
        )
        self._current_positions = state.getPositions()
        self._current_box_vectors = state.getPeriodicBoxVectors()

        # Log final energy for diagnostics
        final_energy = state.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
        LOGGER.info(f"Stage {stage_index} ({stage_name}): final PE = {final_energy:.2f} kJ/mol")

        # Save checkpoint
        checkpoint_path = phase_dir / f"{stage_name}_checkpoint.chk"
        self._simulation.saveCheckpoint(str(checkpoint_path))

        # Remove position restraints from system for next stage
        if restraint_force_indices:
            LOGGER.info(
                f"Removing {len(restraint_force_indices)} position restraint force(s) for next stage"
            )
            remove_position_restraints_from_system(self._system, restraint_force_indices)

        # Build results
        final_temp = stage.get_final_temperature()
        results = {
            "stage_index": stage_index,
            "stage_name": stage.name,
            "ensemble": stage.ensemble.value,
            "duration_ns": stage.duration,
            "total_steps": total_steps,
            "temperature_start_K": stage.get_start_temperature(),
            "temperature_end_K": final_temp,
            "is_temperature_ramping": stage.is_temperature_ramping,
            "position_restraints": [
                {"group": r.group, "force_constant": r.force_constant}
                for r in stage.position_restraints
            ],
            "final_energy_kJ_mol": state.getPotentialEnergy().value_in_unit(
                omm_unit.kilojoule_per_mole
            ),
            "trajectory_path": str(traj_path),
            "checkpoint_path": str(checkpoint_path),
        }

        LOGGER.info(f"Equilibration stage '{stage.name}' complete")
        return results

    def run_staged_equilibration(
        self,
        stages: List["EquilibrationStageConfig"],
        atom_group_resolver: "AtomGroupResolver",
        target_temperature: float,
    ) -> Dict[str, Any]:
        """Run complete multi-stage equilibration protocol.

        This method executes a sequence of equilibration stages, each with
        potentially different:
        - Temperature (constant or ramping)
        - Position restraints on different atom groups
        - Thermodynamic ensemble (NVT/NPT)

        Positions carry over between stages, and restraint forces are
        added/removed as needed.

        Args:
            stages: List of EquilibrationStageConfig objects
            atom_group_resolver: Resolver for predefined atom group names
            target_temperature: Final target temperature (for logging)

        Returns:
            Dictionary with all stage results and summary
        """
        LOGGER.info(f"Starting multi-stage equilibration with {len(stages)} stages")

        # Store reference positions for restraints (post-minimization)
        reference_positions = self._current_positions

        results = {
            "type": "staged_equilibration",
            "num_stages": len(stages),
            "stages": [],
            "total_duration_ns": 0.0,
        }

        for i, stage in enumerate(stages):
            stage_result = self.run_equilibration_stage(
                stage=stage,
                reference_positions=reference_positions,
                atom_group_resolver=atom_group_resolver,
                stage_index=i,
            )
            results["stages"].append(stage_result)
            results["total_duration_ns"] += stage.duration

        # Get final energy
        if results["stages"]:
            results["final_energy_kJ_mol"] = results["stages"][-1]["final_energy_kJ_mol"]
            results["final_temperature_K"] = results["stages"][-1]["temperature_end_K"]

        self._history["equilibration"] = results
        LOGGER.info(
            f"Multi-stage equilibration complete: {len(stages)} stages, "
            f"{results['total_duration_ns']:.3f} ns total"
        )

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

        # Capture velocities and box vectors from equilibration before creating new Simulation
        # (creating new Simulation destroys the old context)
        # Box vectors are critical for NPT stages where box dimensions change
        equilibration_velocities = None
        equilibration_box_vectors = None
        if self._simulation is not None and segment_index == 0:
            state = self._simulation.context.getState(getVelocities=True)
            equilibration_velocities = state.getVelocities()
            equilibration_box_vectors = state.getPeriodicBoxVectors()

        self._simulation = Simulation(self._topology, self._system, integrator, platform)

        # Set box vectors BEFORE positions - critical for correct periodic boundary handling
        # Use captured box vectors from equilibration, or fall back to stored ones from staged equilibration
        if equilibration_box_vectors is not None:
            self._simulation.context.setPeriodicBoxVectors(*equilibration_box_vectors)
        elif self._current_box_vectors is not None:
            self._simulation.context.setPeriodicBoxVectors(*self._current_box_vectors)

        self._simulation.context.setPositions(self._current_positions)

        # Log initial energy for diagnostics
        _state = self._simulation.context.getState(getEnergy=True)
        _energy = _state.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
        LOGGER.info(f"Production segment {segment_index}: initial PE = {_energy:.2f} kJ/mol")

        # Set velocities for production
        # - If we have velocities from equilibration, use them (physical continuity)
        # - Otherwise generate new velocities at target temperature
        # Note: For daisy-chain continuation (segment > 0), ContinuationManager uses
        # loadState() which restores both positions and velocities from the XML state file
        if segment_index == 0:
            if equilibration_velocities is not None:
                self._simulation.context.setVelocities(equilibration_velocities)
                LOGGER.info("Using velocities preserved from equilibration")
            else:
                self._simulation.context.setVelocitiesToTemperature(temperature * omm_unit.kelvin)
                LOGGER.info("Initialized velocities from Maxwell-Boltzmann distribution")

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

        # Get final state (no enforcePeriodicBox to preserve molecular continuity)
        state = self._simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            getEnergy=True,
            getForces=True,
            getParameters=True,
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
                                    "__values__": {"value": friction, "unit": "/picosecond"},
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
