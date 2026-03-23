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
from openmm import XmlSerializer
from openmm import unit as omm_unit
from openmm.app import CheckpointReporter, DCDReporter, PDBFile, Simulation, StateDataReporter

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
        self._current_velocities = None  # Carried between equilibration stages
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
        from polyzymd.utils import impose_unique_force_groups

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

        # Get final state (including box vectors for proper handoff to equilibration)
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        energy = state.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
        self._current_positions = state.getPositions()
        self._current_box_vectors = state.getPeriodicBoxVectors()

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

        # Get final state (including box vectors for potential NPT follow-up)
        state = self._simulation.context.getState(
            getPositions=True, getVelocities=True, getEnergy=True
        )
        self._current_positions = state.getPositions()
        self._current_velocities = state.getVelocities()
        self._current_box_vectors = state.getPeriodicBoxVectors()

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
        resume_from_step: int = 0,
        resume_temperature: float | None = None,
    ) -> Dict[str, Any]:
        """Run a single equilibration stage with optional position restraints.

        This method runs one stage of a multi-stage equilibration protocol.
        It supports:
        - Position restraints on predefined atom groups
        - Temperature ramping (simulated annealing)
        - NVT or NPT ensembles
        - Graceful interruption and mid-stage resume

        Args:
            stage: EquilibrationStageConfig with stage settings
            reference_positions: Positions to restrain atoms to (typically post-minimization)
            atom_group_resolver: Resolver for predefined atom group names
            stage_index: Index of this stage (for output naming)
            default_timestep: Default time step in fs if not specified in stage
            default_friction: Default friction coefficient in 1/ps
            resume_from_step: Step to resume from within this stage (0 = start fresh).
                When resuming, the simulation must already be loaded from checkpoint.
            resume_temperature: Current temperature when resuming a temperature-ramping
                stage. Only used when ``resume_from_step > 0`` and the stage uses
                temperature ramping.

        Returns:
            Dictionary with stage results
        """
        from polyzymd.config.schema import Ensemble
        from polyzymd.core.position_restraints import (
            add_position_restraints_to_system,
            remove_position_restraints_from_system,
        )

        stage_name = f"equilibration_{stage_index}_{stage.name}"
        LOGGER.info(f"Starting equilibration stage: {stage.name} ({stage.duration} ns)")

        # Create output directory for this stage
        phase_dir = self._working_dir / stage_name
        phase_dir.mkdir(exist_ok=True)

        # Get stage parameters (with defaults)
        timestep_fs = stage.time_step if stage.time_step is not None else default_timestep
        thermostat_timescale = stage.thermostat_timescale if stage.thermostat_timescale else 1.0
        friction = 1.0 / thermostat_timescale  # friction (1/ps) = 1 / timescale (ps)

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

        # Velocity initialization: only generate fresh Maxwell-Boltzmann
        # velocities for the first stage. Subsequent stages inherit velocities
        # from the previous stage for physical continuity — matching the
        # GROMACS convention (gen_vel=yes only for stage 0).
        if stage_index == 0 or self._current_velocities is None:
            self._simulation.context.setVelocitiesToTemperature(start_temp * omm_unit.kelvin)
            LOGGER.info(f"Stage {stage_index}: initialized velocities at {start_temp} K")
        else:
            self._simulation.context.setVelocities(self._current_velocities)
            LOGGER.info(f"Stage {stage_index}: inherited velocities from previous stage")

        # Log initial energy
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

        # Periodic checkpoint reporter so state survives hard kills
        eq_chk_path = phase_dir / f"{stage_name}_checkpoint.chk"
        self._simulation.reporters.append(CheckpointReporter(str(eq_chk_path), report_interval))

        # Install signal handlers for graceful shutdown
        from polyzymd.simulation.signals import (
            GracefulExit,
            get_interrupt_signal,
            install_handlers,
            is_interrupted,
        )

        install_handlers()

        # Helper: save EQ_INTERRUPTED marker for mid-stage resume
        def _save_eq_interrupted(steps_done: int, current_temp: float) -> None:
            marker_path = phase_dir / "EQ_INTERRUPTED"
            marker_path.write_text(
                f"stage_index={stage_index}\n"
                f"stage_name={stage.name}\n"
                f"steps_completed={steps_done}\n"
                f"total_steps={total_steps}\n"
                f"current_temperature={current_temp}\n"
                f"is_temperature_ramping={stage.is_temperature_ramping}\n"
            )
            LOGGER.info(
                f"Saved EQ_INTERRUPTED marker: stage {stage_index}, "
                f"step {steps_done}/{total_steps}, temp {current_temp} K"
            )

        # Track steps completed (starting from resume point)
        steps_done = resume_from_step

        # Run simulation with temperature ramping if needed
        if stage.is_temperature_ramping:
            LOGGER.info(
                f"Temperature ramping: {stage.temperature_start} K -> {stage.temperature_end} K "
                f"(increment={stage.temperature_increment} K every {stage.temperature_interval} fs)"
            )
            steps_per_update = int(stage.temperature_interval / timestep_fs)

            # Calculate total temperature updates needed
            temp_range = stage.temperature_end - stage.temperature_start
            num_updates = int(temp_range / stage.temperature_increment)
            steps_for_ramping = num_updates * steps_per_update
            remaining_steps_at_final = total_steps - steps_for_ramping

            # Determine starting temperature — always begin from
            # temperature_start and let the fast-forward loop advance
            # current_temp by skipping already-completed chunks.  On
            # resume, resume_temperature is used only for logging.
            current_temp = stage.temperature_start
            if resume_from_step > 0 and resume_temperature is not None:
                LOGGER.info(
                    f"Resuming temperature ramp from step {resume_from_step}, "
                    f"saved temp {resume_temperature} K (fast-forwarding from {current_temp} K)"
                )

            # Temperature ramping phase — each update is a chunk we can
            # interrupt between.  Skip chunks already completed on resume.
            ramp_step_count = 0
            while current_temp < stage.temperature_end:
                chunk_end = ramp_step_count + steps_per_update
                if chunk_end <= resume_from_step:
                    # Already completed this chunk in a previous run
                    ramp_step_count = chunk_end
                    current_temp += stage.temperature_increment
                    continue

                integrator.setTemperature(current_temp * omm_unit.kelvin)
                if stage.ensemble == Ensemble.NPT:
                    self._simulation.context.setParameter(
                        openmm.MonteCarloBarostat.Temperature(),
                        current_temp * omm_unit.kelvin,
                    )

                # If resuming mid-chunk, only run the remainder
                steps_already = max(0, resume_from_step - ramp_step_count)
                steps_this_chunk = steps_per_update - steps_already
                self._simulation.step(steps_this_chunk)
                steps_done += steps_this_chunk
                ramp_step_count = chunk_end

                if is_interrupted():
                    LOGGER.warning(
                        f"Interrupt during equilibration stage {stage_index} "
                        f"at step {steps_done}/{total_steps} (ramping, T={current_temp:.1f} K)"
                    )
                    _save_eq_interrupted(steps_done, current_temp)
                    raise GracefulExit(
                        signal_number=get_interrupt_signal(), steps_completed=steps_done
                    )

                current_temp += stage.temperature_increment

            # Final temperature - run remaining steps in chunks
            integrator.setTemperature(stage.temperature_end * omm_unit.kelvin)
            if stage.ensemble == Ensemble.NPT:
                self._simulation.context.setParameter(
                    openmm.MonteCarloBarostat.Temperature(),
                    stage.temperature_end * omm_unit.kelvin,
                )
            current_temp = stage.temperature_end
            steps_at_final_done = steps_done - steps_for_ramping
            steps_at_final_remaining = max(0, remaining_steps_at_final - steps_at_final_done)

            if steps_at_final_remaining > 0:
                LOGGER.info(
                    f"Running {steps_at_final_remaining} steps at final "
                    f"temperature {stage.temperature_end} K"
                )
                chunk_size = min(report_interval, steps_at_final_remaining)
                while steps_at_final_remaining > 0:
                    this_chunk = min(chunk_size, steps_at_final_remaining)
                    self._simulation.step(this_chunk)
                    steps_done += this_chunk
                    steps_at_final_remaining -= this_chunk

                    if is_interrupted():
                        LOGGER.warning(
                            f"Interrupt during equilibration stage {stage_index} "
                            f"at step {steps_done}/{total_steps} (final temp)"
                        )
                        _save_eq_interrupted(steps_done, current_temp)
                        raise GracefulExit(
                            signal_number=get_interrupt_signal(),
                            steps_completed=steps_done,
                        )
        else:
            # Constant temperature - run in chunks with signal checking
            current_temp = stage.temperature
            steps_remaining = total_steps - resume_from_step
            if resume_from_step > 0:
                LOGGER.info(
                    f"Resuming constant-temp stage from step {resume_from_step}, "
                    f"{steps_remaining} steps remaining"
                )
            else:
                LOGGER.info(f"Running {total_steps} steps at {stage.temperature} K")

            chunk_size = min(report_interval, steps_remaining)
            while steps_remaining > 0:
                this_chunk = min(chunk_size, steps_remaining)
                self._simulation.step(this_chunk)
                steps_done += this_chunk
                steps_remaining -= this_chunk

                if is_interrupted():
                    LOGGER.warning(
                        f"Interrupt during equilibration stage {stage_index} "
                        f"at step {steps_done}/{total_steps}"
                    )
                    _save_eq_interrupted(steps_done, current_temp)
                    raise GracefulExit(
                        signal_number=get_interrupt_signal(), steps_completed=steps_done
                    )

        # Get final state (including box vectors for NPT stages)
        state = self._simulation.context.getState(
            getPositions=True, getVelocities=True, getEnergy=True
        )
        self._current_positions = state.getPositions()
        self._current_velocities = state.getVelocities()
        self._current_box_vectors = state.getPeriodicBoxVectors()

        # Log final energy for diagnostics
        final_energy = state.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
        LOGGER.info(f"Stage {stage_index} ({stage_name}): final PE = {final_energy:.2f} kJ/mol")

        # Save checkpoint
        checkpoint_path = phase_dir / f"{stage_name}_checkpoint.chk"
        self._simulation.saveCheckpoint(str(checkpoint_path))

        # Remove EQ_INTERRUPTED marker if present (stage completed successfully)
        eq_interrupted_marker = phase_dir / "EQ_INTERRUPTED"
        if eq_interrupted_marker.exists():
            eq_interrupted_marker.unlink()
            LOGGER.info("Removed EQ_INTERRUPTED marker — stage completed successfully")

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

    def _find_completed_eq_stages(
        self,
        stages: List["EquilibrationStageConfig"],
    ) -> List[int]:
        """Find equilibration stages whose checkpoints exist on disk.

        Scans ``equilibration_N_name/`` directories for checkpoint files.
        Stops at the first gap — stages must be contiguous from index 0.
        An interrupted stage (has ``EQ_INTERRUPTED`` marker) is NOT
        considered completed and stops the scan.

        Parameters
        ----------
        stages : list of EquilibrationStageConfig
            The full list of stages from the config.

        Returns
        -------
        list of int
            Indices of completed stages (contiguous from 0).
        """
        completed: List[int] = []
        for i, stage in enumerate(stages):
            stage_name = f"equilibration_{i}_{stage.name}"
            stage_dir = self._working_dir / stage_name
            chk = stage_dir / f"{stage_name}_checkpoint.chk"
            eq_marker = stage_dir / "EQ_INTERRUPTED"
            if chk.exists() and not eq_marker.exists():
                completed.append(i)
            else:
                break  # Stop at first gap — can't skip stages
        return completed

    def _find_interrupted_eq_stage(
        self,
        stages: List["EquilibrationStageConfig"],
        completed_indices: List[int],
    ) -> Dict[str, Any] | None:
        """Check if the next unfinished equilibration stage was interrupted mid-run.

        Looks for an ``EQ_INTERRUPTED`` marker in the stage directory that
        follows the last completed stage. If found, parses resume metadata.

        Parameters
        ----------
        stages : list of EquilibrationStageConfig
            The full list of stages from the config.
        completed_indices : list of int
            Indices of fully completed stages (from ``_find_completed_eq_stages``).

        Returns
        -------
        dict or None
            Dictionary with ``stage_index``, ``steps_completed``, ``total_steps``,
            and ``current_temperature`` if an interrupted stage is found;
            ``None`` otherwise.
        """
        next_idx = len(completed_indices)
        if next_idx >= len(stages):
            return None

        stage = stages[next_idx]
        stage_name = f"equilibration_{next_idx}_{stage.name}"
        stage_dir = self._working_dir / stage_name
        marker_path = stage_dir / "EQ_INTERRUPTED"

        if not marker_path.exists():
            return None

        # Parse the marker
        info: Dict[str, Any] = {"stage_index": next_idx}
        try:
            text = marker_path.read_text()
            for line in text.strip().splitlines():
                key, _, value = line.partition("=")
                key = key.strip()
                value = value.strip()
                if key == "steps_completed":
                    info["steps_completed"] = int(value)
                elif key == "total_steps":
                    info["total_steps"] = int(value)
                elif key == "current_temperature":
                    info["current_temperature"] = float(value)
                elif key == "is_temperature_ramping":
                    info["is_temperature_ramping"] = value.lower() == "true"
        except (ValueError, OSError) as exc:
            LOGGER.warning(f"Could not parse EQ_INTERRUPTED marker {marker_path}: {exc}")
            return None

        # Verify checkpoint exists (needed for resume)
        chk = stage_dir / f"{stage_name}_checkpoint.chk"
        if not chk.exists():
            LOGGER.warning(
                f"EQ_INTERRUPTED marker found for stage {next_idx} but no checkpoint — "
                f"will restart stage from beginning"
            )
            return None

        LOGGER.info(
            f"Found interrupted equilibration stage {next_idx} ({stage.name}): "
            f"{info.get('steps_completed', 0)}/{info.get('total_steps', '?')} steps, "
            f"T={info.get('current_temperature', '?')} K"
        )
        return info

    def _load_eq_stage_state(
        self,
        stage_index: int,
        stage_name: str,
    ) -> None:
        """Load positions, velocities, and box vectors from a completed equilibration checkpoint.

        Creates a temporary ``Simulation`` to deserialise the binary
        checkpoint, then stores the extracted state on ``self`` so
        subsequent stages or production can pick up seamlessly.
        The temporary simulation is discarded after extraction.

        Parameters
        ----------
        stage_index : int
            Index of the completed stage.
        stage_name : str
            Name of the completed stage (used in directory/file naming).
        """
        dir_name = f"equilibration_{stage_index}_{stage_name}"
        chk_path = self._working_dir / dir_name / f"{dir_name}_checkpoint.chk"

        if not chk_path.exists():
            raise FileNotFoundError(f"Equilibration checkpoint not found: {chk_path}")

        # Temporary simulation context to deserialise the binary checkpoint.
        # We use a dummy VerletIntegrator because we only need to extract
        # geometric state (positions, velocities, box vectors) — the next
        # run_equilibration_stage() call creates a proper LangevinMiddleIntegrator.
        integrator = openmm.VerletIntegrator(1.0 * omm_unit.femtosecond)
        platform = self._get_platform()
        temp_sim = Simulation(self._topology, self._system, integrator, platform)
        temp_sim.loadCheckpoint(str(chk_path))

        state = temp_sim.context.getState(getPositions=True, getVelocities=True)
        self._current_positions = state.getPositions()
        self._current_velocities = state.getVelocities()
        self._current_box_vectors = state.getPeriodicBoxVectors()
        # Discard the temporary simulation — it has a dummy integrator
        del temp_sim

        LOGGER.info(
            f"Loaded state from equilibration stage {stage_index} ({stage_name}) checkpoint"
        )

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
        added/removed as needed. If a previous run was interrupted, completed
        stages are detected on disk via their checkpoint files and skipped
        automatically.

        Args:
            stages: List of EquilibrationStageConfig objects
            atom_group_resolver: Resolver for predefined atom group names
            target_temperature: Final target temperature (for logging)

        Returns:
            Dictionary with all stage results and summary
        """
        import shutil

        LOGGER.info(f"Starting multi-stage equilibration with {len(stages)} stages")

        # Store reference positions for restraints BEFORE loading any checkpoint.
        # These are the post-minimization positions that restraint forces target.
        # Must be captured before _load_eq_stage_state() overwrites _current_positions.
        reference_positions = self._current_positions

        # Detect stages that already completed (checkpoint exists on disk)
        completed_indices = self._find_completed_eq_stages(stages)
        if completed_indices:
            last_idx = completed_indices[-1]
            last_stage = stages[last_idx]
            LOGGER.info(
                f"Found {len(completed_indices)} completed equilibration stage(s) "
                f"on disk — resuming after stage {last_idx} ({last_stage.name})"
            )
            self._load_eq_stage_state(last_idx, last_stage.name)

        # Check if the next stage was interrupted mid-run (has EQ_INTERRUPTED marker)
        interrupted_info = self._find_interrupted_eq_stage(stages, completed_indices)

        results: Dict[str, Any] = {
            "type": "staged_equilibration",
            "num_stages": len(stages),
            "stages": [],
            "total_duration_ns": 0.0,
        }

        for i, stage in enumerate(stages):
            if i in completed_indices:
                LOGGER.info(f"Skipping completed equilibration stage {i}: {stage.name}")
                results["stages"].append(
                    {
                        "stage_index": i,
                        "stage_name": stage.name,
                        "skipped": True,
                        "duration_ns": stage.duration,
                    }
                )
                results["total_duration_ns"] += stage.duration
                continue

            stage_dir_name = f"equilibration_{i}_{stage.name}"
            partial_dir = self._working_dir / stage_dir_name

            # Check if this stage was interrupted and can be resumed
            if interrupted_info is not None and interrupted_info["stage_index"] == i:
                # Resume mid-stage from checkpoint
                resume_step = interrupted_info.get("steps_completed", 0)
                resume_temp = interrupted_info.get("current_temperature")
                LOGGER.info(
                    f"Resuming interrupted equilibration stage {i} ({stage.name}) "
                    f"from step {resume_step}"
                )
                self._load_eq_stage_state(i, stage.name)

                stage_result = self.run_equilibration_stage(
                    stage=stage,
                    reference_positions=reference_positions,
                    atom_group_resolver=atom_group_resolver,
                    stage_index=i,
                    resume_from_step=resume_step,
                    resume_temperature=resume_temp,
                )
                results["stages"].append(stage_result)
                results["total_duration_ns"] += stage.duration
                continue

            # Clean up partial stage directory (exists but no checkpoint and
            # no EQ_INTERRUPTED marker — truly failed, restart from scratch)
            if partial_dir.exists():
                LOGGER.warning(
                    f"Stage {i} ({stage.name}) directory exists without checkpoint "
                    f"— removing partial output and re-running"
                )
                shutil.rmtree(partial_dir)

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
            # Find last non-skipped result
            last_result = None
            for r in reversed(results["stages"]):
                if not r.get("skipped"):
                    last_result = r
                    break
            if last_result:
                results["final_energy_kJ_mol"] = last_result["final_energy_kJ_mol"]
                results["final_temperature_K"] = last_result["temperature_end_K"]

        self._history["equilibration"] = results
        skipped = len(completed_indices)
        ran = len(stages) - skipped
        LOGGER.info(
            f"Multi-stage equilibration complete: {len(stages)} stages "
            f"({skipped} skipped, {ran} ran), "
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
        report_interval: int | None = None,
        checkpoint_interval_s: float = 60.0,
    ) -> Dict[str, Any]:
        """Run NPT production simulation.

        Args:
            temperature: Temperature in Kelvin.
            duration_ns: Duration in nanoseconds.
            num_samples: Number of trajectory frames to save. Ignored when
                ``report_interval`` is provided.
            timestep_fs: Time step in femtoseconds.
            friction: Friction coefficient in 1/ps.
            pressure: Pressure in atmospheres.
            barostat_frequency: Barostat update frequency.
            output_prefix: Prefix for output files.
            segment_index: Segment index for multi-segment production.
            report_interval: Fixed reporter interval in steps. When provided,
                this overrides the per-segment ``total_steps // num_samples``
                calculation to keep frame spacing uniform across segments.
            checkpoint_interval_s: Wall-time interval in seconds between
                portable restart checkpoints.  Also controls how frequently
                the loop checks for SLURM preemption signals.  Set to 0 to
                disable wall-time checkpoints (reverts to legacy behaviour).

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

        # Determine report interval: prefer the fixed global value if given,
        # otherwise fall back to per-segment calculation (legacy callers).
        if report_interval is None:
            report_interval = max(1, total_steps // num_samples)

        # Create integrator and simulation
        integrator = self._create_integrator(
            temperature=temperature,
            friction=friction,
            timestep=timestep_fs,
        )
        platform = self._get_platform()

        self._simulation = Simulation(self._topology, self._system, integrator, platform)

        # Set box vectors BEFORE positions - critical for correct periodic boundary handling
        if self._current_box_vectors is not None:
            self._simulation.context.setPeriodicBoxVectors(*self._current_box_vectors)

        self._simulation.context.setPositions(self._current_positions)

        # Log initial energy
        _state = self._simulation.context.getState(getEnergy=True)
        _energy = _state.getPotentialEnergy().value_in_unit(omm_unit.kilojoule_per_mole)
        LOGGER.info(f"Production segment {segment_index}: initial PE = {_energy:.2f} kJ/mol")

        # Save system XML early so it exists on disk even if the segment is
        # hard-killed (SIGKILL / OOM).  This is required for checkpoint-based
        # recovery: loadCheckpoint() needs a matching System object.  The file
        # is overwritten at segment completion with the final state.
        system_xml_path = phase_dir / f"{phase_name}_system.xml"
        with open(system_xml_path, "w") as f:
            f.write(XmlSerializer.serialize(self._system))
        LOGGER.info(f"Saved initial system to {system_xml_path}")

        # Set velocities for production
        # - If we have velocities from equilibration, use them (physical continuity)
        # - Otherwise generate new velocities at target temperature
        # Note: For continuation segments (segment > 0), ContinuationManager uses
        # loadState() which restores both positions and velocities from the XML state file
        if segment_index == 0:
            if self._current_velocities is not None:
                self._simulation.context.setVelocities(self._current_velocities)
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

        # Periodic checkpoint reporter — ensures a .chk file is written every
        # report_interval steps.  Without this, segment 0 has NO checkpoint on
        # disk until the very end, so a hard kill (SIGKILL / OOM / node failure)
        # would lose all progress.  Matches the behaviour already present in
        # continuation.py for segments >= 1.
        prod_chk_path = phase_dir / f"{phase_name}_checkpoint.chk"
        self._simulation.reporters.append(CheckpointReporter(str(prod_chk_path), report_interval))

        # Save topology
        with open(pdb_path, "w") as f:
            PDBFile.writeFile(
                self._topology,
                self._current_positions,
                f,
            )

        # Save parameters JSON before the simulation loop so it exists even
        # if the segment is interrupted (continuation.py requires this file)
        # Save parameters JSON (needed for continuation across segments)
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

        # Install signal handlers for graceful shutdown (SIGUSR1 / SIGTERM)
        from polyzymd.simulation.signals import (
            GracefulExit,
            get_interrupt_signal,
            install_handlers,
            is_interrupted,
            save_interrupted_state,
            save_restart_checkpoint,
        )

        install_handlers()

        # Mark this segment as RUNNING in progress.json so that
        # check-progress can distinguish actively running simulations
        # from interrupted ones.
        self._write_segment_started(segment_index, total_steps)

        # Run simulation with adaptive sub-chunks for interrupt responsiveness
        # and periodic wall-time restart checkpoints for preemption resilience.
        LOGGER.info(f"Running {total_steps} steps...")
        steps_done = 0
        import time as _time

        from polyzymd.simulation.progress import PROGRESS_UPDATE_INTERVAL_SECONDS

        _last_progress_write = _time.monotonic()
        _last_checkpoint_write = _time.monotonic()
        _loop_start = _time.monotonic()

        # Adaptive sub-chunk sizing: start with report_interval (the original
        # chunk_size).  After the first checkpoint interval elapses, measure
        # actual steps/second and adapt sub_chunk to target
        # checkpoint_interval / 4 seconds (~15s worth of steps).  This
        # ensures ~4 interrupt checks per checkpoint interval regardless of
        # system size or hardware speed.
        sub_chunk = min(report_interval, total_steps)
        _adapted = False

        try:
            while steps_done < total_steps:
                remaining = total_steps - steps_done
                this_chunk = min(sub_chunk, remaining)
                self._simulation.step(this_chunk)
                steps_done += this_chunk

                _now = _time.monotonic()

                # Adaptive sub-chunk calibration (once, after first interval)
                if (
                    not _adapted
                    and checkpoint_interval_s > 0
                    and (_now - _loop_start) >= checkpoint_interval_s
                ):
                    elapsed = _now - _loop_start
                    steps_per_sec = steps_done / elapsed if elapsed > 0 else 1.0
                    # Target sub-chunk duration = checkpoint_interval / 4
                    target_seconds = checkpoint_interval_s / 4.0
                    new_sub_chunk = max(10, int(steps_per_sec * target_seconds))
                    # Sub-chunk must be a divisor-friendly size relative to
                    # report_interval to avoid misaligned reporter writes.
                    # Round down to the nearest multiple of report_interval,
                    # or use report_interval itself if it's already smaller.
                    if new_sub_chunk >= report_interval:
                        new_sub_chunk = report_interval
                    else:
                        # Ensure sub_chunk divides evenly into report_interval
                        # so reporters fire at exact multiples.
                        # Find largest divisor of report_interval <= new_sub_chunk
                        best = new_sub_chunk
                        for candidate in [
                            report_interval // k
                            for k in range(1, report_interval // max(1, new_sub_chunk) + 2)
                        ]:
                            if candidate <= new_sub_chunk and report_interval % candidate == 0:
                                best = candidate
                                break
                        new_sub_chunk = max(10, best)

                    if new_sub_chunk != sub_chunk:
                        LOGGER.info(
                            f"Adaptive sub-chunk: {sub_chunk} -> {new_sub_chunk} steps "
                            f"(~{new_sub_chunk / steps_per_sec:.1f}s at "
                            f"{steps_per_sec:.0f} steps/s)"
                        )
                        sub_chunk = new_sub_chunk
                    _adapted = True

                # Periodically update RUNNING record so check-progress
                # can display real-time remaining nanoseconds.
                if _now - _last_progress_write >= PROGRESS_UPDATE_INTERVAL_SECONDS:
                    self._update_progress_running(
                        segment_index=segment_index,
                        steps_done=steps_done,
                        timestep_fs=timestep_fs,
                    )
                    _last_progress_write = _now

                # Wall-time restart checkpoint
                if (
                    checkpoint_interval_s > 0
                    and (_now - _last_checkpoint_write) >= checkpoint_interval_s
                    and steps_done < total_steps  # skip if we're about to finish
                ):
                    save_restart_checkpoint(
                        simulation=self._simulation,
                        output_dir=phase_dir,
                    )
                    _last_checkpoint_write = _now

                if is_interrupted():
                    LOGGER.warning(f"Interrupt detected at step {steps_done}/{total_steps}")
                    save_interrupted_state(
                        simulation=self._simulation,
                        output_dir=phase_dir,
                        segment_index=segment_index,
                        steps_completed=steps_done,
                        total_steps=total_steps,
                    )
                    # Update progress tracker (interrupted)
                    self._update_progress_interrupted(
                        segment_index=segment_index,
                        steps_done=steps_done,
                        total_steps=total_steps,
                        duration_ns=duration_ns,
                        timestep_fs=timestep_fs,
                    )
                    raise GracefulExit(
                        signal_number=get_interrupt_signal(), steps_completed=steps_done
                    )
        except GracefulExit:
            raise  # Re-raise so caller can handle exit code
        except Exception:
            # On unexpected crash, still try to save interrupted state
            try:
                save_interrupted_state(
                    simulation=self._simulation,
                    output_dir=phase_dir,
                    segment_index=segment_index,
                    steps_completed=steps_done,
                    total_steps=total_steps,
                )
            except Exception:
                LOGGER.error("Failed to save interrupted state after crash")
            raise

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

        # Save state XML (needed for continuation across segments)
        state_xml_path = phase_dir / f"{phase_name}_state.xml"
        with open(state_xml_path, "w") as f:
            f.write(XmlSerializer.serialize(state))
        LOGGER.info(f"Saved state to {state_xml_path}")

        # Save system XML (needed for continuation across segments)
        system_xml_path = phase_dir / f"{phase_name}_system.xml"
        with open(system_xml_path, "w") as f:
            f.write(XmlSerializer.serialize(self._system))
        LOGGER.info(f"Saved system to {system_xml_path}")

        # Update progress tracker (successful completion)
        self._update_progress_completed(
            segment_index=segment_index,
            total_steps=total_steps,
            num_samples=num_samples,
            duration_ns=duration_ns,
            timestep_fs=timestep_fs,
        )

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

    def _write_segment_started(
        self,
        segment_index: int,
        total_steps: int,
    ) -> None:
        """Write a RUNNING segment record to progress.json at segment start.

        This marks the segment as actively executing so that
        ``check-progress`` can distinguish a running simulation from
        one that was interrupted. The record is later updated to
        COMPLETED or INTERRUPTED by the corresponding handler.

        Parameters
        ----------
        segment_index : int
            Production segment index (0 for initial production).
        total_steps : int
            Total steps planned for this segment.
        """
        from polyzymd.simulation.progress import (
            PROGRESS_UPDATE_INTERVAL_SECONDS,
            SegmentRecord,
            SegmentStatus,
            SimulationStatus,
            _update_or_append_segment,
            load_progress,
            save_progress,
        )

        progress = load_progress(self._working_dir)
        if progress is None:
            LOGGER.warning("No progress file found — skipping segment-started write")
            return

        record = SegmentRecord(
            index=segment_index,
            steps_completed=0,
            steps_requested=total_steps,
            samples_written=0,
            status=SegmentStatus.RUNNING,
        )

        _update_or_append_segment(progress, record)
        progress.status = SimulationStatus.RUNNING

        save_progress(self._working_dir, progress)
        LOGGER.info(f"Marked segment {segment_index} as RUNNING in progress file")

    def _update_progress_running(
        self,
        segment_index: int,
        steps_done: int,
        timestep_fs: float,
    ) -> None:
        """Periodically update the RUNNING segment's step count.

        Called during the simulation loop to keep ``progress.json``
        up to date so that ``check-progress`` can show real-time
        remaining nanoseconds.

        Parameters
        ----------
        segment_index : int
            Production segment index.
        steps_done : int
            Steps completed so far in this segment.
        timestep_fs : float
            Integration timestep in femtoseconds.
        """
        from polyzymd.simulation.progress import (
            SegmentStatus,
            SimulationStatus,
            _update_or_append_segment,
            load_progress,
            save_progress,
        )

        progress = load_progress(self._working_dir)
        if progress is None:
            return

        # Find and update the existing RUNNING record
        for seg in progress.segments:
            if seg.index == segment_index and seg.status == SegmentStatus.RUNNING:
                seg.steps_completed = steps_done
                seg.duration_ns = (steps_done * timestep_fs) / 1e6
                break
        else:
            # No RUNNING record found — shouldn't happen, but be safe
            return

        progress.status = SimulationStatus.RUNNING
        save_progress(self._working_dir, progress)
        LOGGER.debug(f"Updated running progress: segment {segment_index}, {steps_done} steps")

    def _update_progress_completed(
        self,
        segment_index: int,
        total_steps: int,
        num_samples: int,
        duration_ns: float,
        timestep_fs: float,
    ) -> None:
        """Update progress file after successful segment completion.

        Parameters
        ----------
        segment_index : int
            Production segment index (0 for initial production).
        total_steps : int
            Steps completed in this segment.
        num_samples : int
            Samples written in this segment.
        duration_ns : float
            Simulation time of this segment in nanoseconds.
        timestep_fs : float
            Integration timestep in femtoseconds.
        """
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            SimulationStatus,
            _now_iso,
            _update_or_append_segment,
            load_progress,
            save_progress,
        )

        progress = load_progress(self._working_dir)
        if progress is None:
            LOGGER.warning("No progress file found — skipping progress update")
            return

        record = SegmentRecord(
            index=segment_index,
            steps_completed=total_steps,
            steps_requested=total_steps,
            samples_written=num_samples,
            status=SegmentStatus.COMPLETED,
            duration_ns=duration_ns,
        )
        record.finished_at = _now_iso()

        _update_or_append_segment(progress, record)

        # Update overall status
        if progress.is_complete:
            progress.status = SimulationStatus.COMPLETED
        else:
            progress.status = SimulationStatus.RUNNING

        save_progress(self._working_dir, progress)
        LOGGER.info(
            f"Progress updated: {progress.total_steps_completed}/"
            f"{progress.total_steps_requested} steps "
            f"({progress.fraction_complete():.1%})"
        )

    def _update_progress_interrupted(
        self,
        segment_index: int,
        steps_done: int,
        total_steps: int,
        duration_ns: float,
        timestep_fs: float,
    ) -> None:
        """Update progress file after segment interruption.

        Parameters
        ----------
        segment_index : int
            Production segment index.
        steps_done : int
            Steps completed before interruption.
        total_steps : int
            Steps that were planned for this segment.
        duration_ns : float
            Planned simulation time of this segment in nanoseconds.
        timestep_fs : float
            Integration timestep in femtoseconds.
        """
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            SimulationStatus,
            _update_or_append_segment,
            load_progress,
            save_progress,
        )

        progress = load_progress(self._working_dir)
        if progress is None:
            LOGGER.warning("No progress file found — skipping progress update")
            return

        # Calculate actual duration completed
        actual_duration_ns = (steps_done * timestep_fs) / 1e6

        record = SegmentRecord(
            index=segment_index,
            steps_completed=steps_done,
            steps_requested=total_steps,
            samples_written=0,  # Interrupted — samples may be partial
            status=SegmentStatus.INTERRUPTED,
            duration_ns=actual_duration_ns,
        )

        _update_or_append_segment(progress, record)
        progress.status = SimulationStatus.INTERRUPTED

        save_progress(self._working_dir, progress)
        LOGGER.info(
            f"Progress updated (interrupted): {steps_done}/{total_steps} steps "
            f"in segment {segment_index}"
        )

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

        # Update current positions, velocities, and box vectors
        state = self._simulation.context.getState(
            getPositions=True, getVelocities=True
        )
        self._current_positions = state.getPositions()
        self._current_velocities = state.getVelocities()
        self._current_box_vectors = state.getPeriodicBoxVectors()

        LOGGER.info(f"Loaded checkpoint from {checkpoint_path}")
