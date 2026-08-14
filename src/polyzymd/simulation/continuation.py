"""
Continuation manager for resuming MD simulations from checkpoints.

This module handles loading simulation state from previous segments
and continuing the simulation for self-resubmitting HPC workflows.
Each segment runs until completion or interruption (wall-time / preemption),
updates the progress tracker, and the SLURM script handles resubmission.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Optional, Union

if TYPE_CHECKING:
    import openmm
    from openmm.app import Simulation
    from openmm.unit import Quantity

LOGGER = logging.getLogger(__name__)


def _get_openmm_module() -> Any:
    """Import and return the OpenMM top-level module lazily."""

    import openmm

    return openmm


def _get_openmm_app_classes() -> tuple[Any, Any, Any, Any, Any]:
    """Import and return OpenMM app classes lazily."""

    from openmm.app import (
        CheckpointReporter,
        DCDReporter,
        PDBFile,
        Simulation,
        StateDataReporter,
    )

    return CheckpointReporter, DCDReporter, PDBFile, Simulation, StateDataReporter


def _get_xml_serializer() -> Any:
    """Import and return the OpenMM XML serializer lazily."""

    from openmm import XmlSerializer

    return XmlSerializer


def _get_openmm_unit_module() -> Any:
    """Import and return the OpenMM unit module lazily."""

    from openmm import unit as u

    return u


def quantity_from_dict(qdict: Dict[str, Any]) -> Quantity:
    """Convert serialized quantity dictionary back to OpenMM Quantity.

    Args:
        qdict: Dictionary with __values__ containing value and unit.

    Returns:
        OpenMM Quantity with appropriate units.
    """
    u = _get_openmm_unit_module()

    value = qdict["__values__"]["value"]
    unit_str = qdict["__values__"]["unit"]

    # Handle inverse units (e.g., "/picosecond")
    if unit_str.startswith("/"):
        base_unit = getattr(u, unit_str[1:])
        return value / base_unit

    # Map common unit variations
    unit_mapping = {
        "atmosphere": u.atmospheres,
        "atmospheres": u.atmospheres,
        "kelvin": u.kelvin,
        "femtosecond": u.femtoseconds,
        "femtoseconds": u.femtoseconds,
        "nanosecond": u.nanoseconds,
        "nanoseconds": u.nanoseconds,
        "picosecond": u.picoseconds,
        "picoseconds": u.picoseconds,
    }

    if unit_str in unit_mapping:
        return value * unit_mapping[unit_str]
    else:
        return value * getattr(u, unit_str)


class ContinuationManager:
    """Manager for continuing MD simulations from previous segments.

    This class handles loading state from previous production segments
    and continuing the simulation. It integrates with the progress
    tracking system to enable self-resubmitting idempotent jobs.

    Example
    -------
    >>> manager = ContinuationManager(
    ...     working_dir="simulation_output/",
    ...     segment_index=2,  # Continuing to segment 2
    ... )
    >>> manager.load_previous_state()
    >>> manager.run_segment(duration_ns=20.0, num_samples=250)
    """

    def __init__(
        self,
        working_dir: Union[str, Path],
        segment_index: int,
        platform: str = "CUDA",
        precision: str = "mixed",
        device_index: str | None = None,
    ) -> None:
        """Initialize the ContinuationManager.

        Parameters
        ----------
        working_dir : str or Path
            Working directory containing simulation outputs.
        segment_index : int
            Current segment index (0-based for first continuation after
            initial production, incrementing from there).
        """
        self._working_dir = Path(working_dir)
        self._segment_index = segment_index
        self._prev_segment = segment_index - 1
        self._platform_name = platform
        self._platform_precision = precision
        self._platform_device_index = device_index

        # State
        self._system: Optional[openmm.System] = None
        self._topology: Optional[Any] = None
        self._simulation: Optional[Simulation] = None
        self._param_dict: Optional[Dict[str, Any]] = None
        self._use_checkpoint_recovery: bool = False

    @property
    def working_dir(self) -> Path:
        """Get the working directory."""
        return self._working_dir

    @property
    def segment_index(self) -> int:
        """Get the current segment index."""
        return self._segment_index

    @property
    def simulation(self) -> Optional[Simulation]:
        """Get the OpenMM Simulation object."""
        return self._simulation

    def _find_solvated_pdb(self) -> Path:
        """Find the solvated PDB file in the working directory.

        Returns
        -------
        Path
            Path to the solvated PDB file.

        Raises
        ------
        FileNotFoundError
            If no suitable PDB file is found.
        """
        allowed_paths = (
            self._working_dir / "solvated_system.pdb",
            self._working_dir / "production_0" / "production_0_topology.pdb",
            self._working_dir / "production" / "production_topology.pdb",
        )
        for pdb_path in allowed_paths:
            if pdb_path.exists():
                return pdb_path

        # Arbitrary recursive PDB discovery is disallowed to avoid selecting decoys or inputs

        raise FileNotFoundError(f"Could not find solvated PDB file in {self._working_dir}")

    def _get_previous_paths(self) -> Dict[str, Path]:
        """Get paths to files from the previous segment.

        Returns
        -------
        dict
            Dictionary with paths to state, system, and parameter files.

            Recovery priority (portable state XML preferred over binary .chk):

            1. **Normal completion** — ``production_N_state.xml`` and
               ``production_N_system.xml`` exist.
            2. **Graceful interruption with interrupted state** —
               ``interrupted_state.xml`` saved by signal handler.  Uses
               ``loadState()`` (portable).
            3. **Graceful interruption with restart checkpoint** —
               ``restart_state.xml`` saved by wall-time checkpoint loop.
               Uses ``loadState()`` (portable).
            4. **Emergency interruption (.chk only)** —
               ``interrupted_checkpoint.chk`` exists but no state XML.
               Falls back to ``loadCheckpoint()`` (non-portable) only after
               XML state recovery is unavailable.
            5. **Hard kill/OOM/node failure** — periodic ``checkpoint.chk``
               from CheckpointReporter exists but no XML state files.
               Falls back to ``loadCheckpoint()`` (non-portable).
        """
        prev_dir = self._working_dir / f"production_{self._prev_segment}"

        state_path = prev_dir / f"production_{self._prev_segment}_state.xml"
        system_path = prev_dir / f"production_{self._prev_segment}_system.xml"
        checkpoint_path = prev_dir / f"production_{self._prev_segment}_checkpoint.chk"
        params_path = prev_dir / f"production_{self._prev_segment}_parameters.json"

        # Portable state XMLs from interruption handlers
        interrupted_state = prev_dir / "interrupted_state.xml"
        interrupted_system = prev_dir / "interrupted_system.xml"
        restart_state = prev_dir / "restart_state.xml"
        restart_system = prev_dir / "restart_system.xml"
        interrupted_chk = prev_dir / "interrupted_checkpoint.chk"

        use_checkpoint = False

        if state_path.exists():
            # Case 1: Normal completion — state.xml exists
            pass

        elif interrupted_state.exists():
            # Case 2: Graceful interruption — portable interrupted_state.xml
            LOGGER.info(
                f"Previous segment {self._prev_segment} was interrupted — "
                f"recovering from interrupted_state.xml (portable)"
            )
            state_path = interrupted_state
            if interrupted_system.exists():
                system_path = interrupted_system

        elif restart_state.exists():
            # Case 3: Interrupted between checkpoints — wall-time restart
            LOGGER.info(
                f"Previous segment {self._prev_segment} was interrupted — "
                f"recovering from restart_state.xml (portable wall-time checkpoint)"
            )
            state_path = restart_state
            if restart_system.exists():
                system_path = restart_system

        elif interrupted_chk.exists() and interrupted_system.exists():
            # Case 4: Emergency recovery when only binary .chk survived
            LOGGER.warning(
                f"Previous segment {self._prev_segment} was interrupted — "
                f"no portable state XML found, using emergency "
                f"interrupted_checkpoint.chk recovery (non-portable)"
            )
            system_path = interrupted_system
            checkpoint_path = interrupted_chk
            use_checkpoint = True

        elif checkpoint_path.exists() and system_path.exists():
            # Case 5: Emergency hard-kill/OOM/node-failure .chk recovery
            LOGGER.warning(
                f"Previous segment {self._prev_segment} appears hard-killed — "
                f"no state XML or interrupted files, recovering from "
                f"periodic checkpoint + early-saved system XML (non-portable)"
            )
            use_checkpoint = True

        elif checkpoint_path.exists():
            # Case 5b: Hard kill but no system.xml — cannot recover
            raise FileNotFoundError(
                f"Previous segment {self._prev_segment} has a checkpoint "
                f"({checkpoint_path}) but no system.xml ({system_path}) — "
                f"cannot recover.  The segment must be re-run from scratch."
            )

        return {
            "state": state_path,
            "system": system_path,
            "params": params_path,
            "checkpoint": checkpoint_path,
            "use_checkpoint": use_checkpoint,  # type: ignore[dict-item]
        }

    def load_previous_state(self) -> None:
        """Load state from the previous segment.

        This loads the system, topology, and parameters from the previous
        production segment.  Recovery prefers portable state XML files
        (``production_N_state.xml``, ``interrupted_state.xml``, or
        ``restart_state.xml``) over binary ``.chk`` checkpoints.  Only
        falls back to ``loadCheckpoint()`` when no portable state XML
        is available (legacy interrupted segments, hard-killed segments).

        Raises
        ------
        FileNotFoundError
            If required files are missing.
        """
        LOGGER.info(f"Loading state from segment {self._prev_segment}")

        paths = self._get_previous_paths()
        use_checkpoint = bool(paths.pop("use_checkpoint", False))

        # Check that required files exist
        for name, path in paths.items():
            if name == "state" and use_checkpoint:
                # State XML doesn't exist for interrupted/hard-killed segments;
                # we'll use the checkpoint instead in run_segment()
                continue
            if name == "checkpoint" and not use_checkpoint:
                # Checkpoint only required when recovering from interruption
                continue
            if not path.exists():
                raise FileNotFoundError(f"Required file not found: {path}")

        # Load system (either normal or interrupted system XML)
        LOGGER.info(f"Loading system from {paths['system']}")
        XmlSerializer = _get_xml_serializer()
        _, _, PDBFile, _, _ = _get_openmm_app_classes()
        with open(paths["system"], "r") as f:
            self._system = XmlSerializer.deserialize(f.read())

        # Load topology
        pdb_path = self._find_solvated_pdb()
        LOGGER.info(f"Loading topology from {pdb_path}")
        self._topology = PDBFile(str(pdb_path)).topology

        # Load parameters
        LOGGER.info(f"Loading parameters from {paths['params']}")
        with open(paths["params"], "r") as f:
            self._param_dict = json.load(f)

        # Store whether we need checkpoint recovery for run_segment()
        self._use_checkpoint_recovery = use_checkpoint

        LOGGER.info("Previous state loaded successfully")

    def _create_integrator(self) -> openmm.Integrator:
        """Create an integrator from the parameter dictionary.

        Returns
        -------
        openmm.Integrator
            OpenMM LangevinMiddleIntegrator.
        """
        if self._param_dict is None:
            raise RuntimeError("Parameters not loaded. Call load_previous_state first.")

        openmm = _get_openmm_module()

        integ_raw = self._param_dict["__values__"]["integ_params"]["__values__"]
        time_step = quantity_from_dict(integ_raw["time_step"])

        thermo_raw = self._param_dict["__values__"]["thermo_params"]["__values__"]
        thermostat_raw = thermo_raw["thermostat_params"]["__values__"]

        temperature = quantity_from_dict(thermostat_raw["temperature"])
        friction_coeff = quantity_from_dict(thermostat_raw["timescale"])

        return openmm.LangevinMiddleIntegrator(temperature, friction_coeff, time_step)

    def _add_barostat_if_needed(self) -> None:
        """Add barostat to the system if parameters specify NPT."""
        if self._system is None or self._param_dict is None:
            raise RuntimeError("System/parameters not loaded")

        openmm = _get_openmm_module()

        thermo_raw = self._param_dict["__values__"]["thermo_params"]["__values__"]

        if "barostat_params" not in thermo_raw:
            return

        # Check if barostat already exists
        has_barostat = any(
            isinstance(self._system.getForce(i), openmm.MonteCarloBarostat)
            for i in range(self._system.getNumForces())
        )

        if has_barostat:
            LOGGER.debug("Barostat already present")
            return

        barostat_raw = thermo_raw["barostat_params"]["__values__"]
        temperature = quantity_from_dict(barostat_raw["temperature"])
        pressure = quantity_from_dict(barostat_raw["pressure"])
        frequency = barostat_raw.get("update_frequency", 25)

        barostat = openmm.MonteCarloBarostat(pressure, temperature, frequency)
        self._system.addForce(barostat)
        LOGGER.info(f"Added barostat: {pressure} at {temperature}")

    def _setup_reporters(
        self,
        report_interval: int,
        output_dir: Path,
    ) -> None:
        """Setup reporters for the simulation.

        Parameters
        ----------
        report_interval : int
            Step interval between reporter outputs. Kept constant
            across all segments to ensure uniform frame spacing.
        output_dir : Path
            Output directory for this segment.
        """
        if self._simulation is None:
            raise RuntimeError("Simulation not created")

        CheckpointReporter, DCDReporter, _, _, StateDataReporter = _get_openmm_app_classes()

        # Trajectory reporter
        traj_path = output_dir / f"production_{self._segment_index}_trajectory.dcd"
        self._simulation.reporters.append(DCDReporter(str(traj_path), report_interval))

        # State data reporter
        state_path = output_dir / f"production_{self._segment_index}_state_data.csv"
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

        # Checkpoint reporter
        checkpoint_path = output_dir / f"production_{self._segment_index}_checkpoint.chk"
        self._simulation.reporters.append(CheckpointReporter(str(checkpoint_path), report_interval))

        LOGGER.info(f"Setup reporters with interval {report_interval}")

    def _save_final_state(self, output_dir: Path) -> None:
        """Save the final state and system after simulation.

        Parameters
        ----------
        output_dir : Path
            Output directory for this segment.
        """
        if self._simulation is None:
            raise RuntimeError("Simulation not available")

        XmlSerializer = _get_xml_serializer()

        # Save state (no enforcePeriodicBox to preserve molecular continuity)
        state_path = output_dir / f"production_{self._segment_index}_state.xml"
        state = self._simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            getForces=True,
            getEnergy=True,
            getParameters=True,
        )

        with open(state_path, "w") as f:
            f.write(XmlSerializer.serialize(state))

        # Save system
        system_path = output_dir / f"production_{self._segment_index}_system.xml"
        with open(system_path, "w") as f:
            f.write(XmlSerializer.serialize(self._simulation.system))

        LOGGER.info(f"Saved final state to {state_path}")
        LOGGER.info(f"Saved system to {system_path}")

    def _write_segment_started(self, total_steps: int) -> None:
        """Write a RUNNING segment record to progress.json at segment start.

        This marks the segment as actively executing so that
        ``check-progress`` can distinguish a running simulation from
        one that was interrupted. The record is later updated to
        COMPLETED or INTERRUPTED by the corresponding handler.

        Parameters
        ----------
        total_steps : int
            Total steps planned for this segment.
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
            LOGGER.warning("No progress file found — skipping segment-started write")
            return

        record = SegmentRecord(
            index=self._segment_index,
            steps_completed=0,
            steps_requested=total_steps,
            samples_written=0,
            status=SegmentStatus.RUNNING,
        )

        _update_or_append_segment(progress, record)
        progress.status = SimulationStatus.RUNNING

        save_progress(self._working_dir, progress)
        LOGGER.info(f"Marked segment {self._segment_index} as RUNNING in progress file")

    def _update_progress_running(
        self,
        steps_done: int,
        timestep_fs: float,
    ) -> None:
        """Periodically update the RUNNING segment's step count.

        Called during the simulation loop to keep ``progress.json``
        up to date so that ``check-progress`` can show real-time
        remaining nanoseconds.

        Parameters
        ----------
        steps_done : int
            Steps completed so far in this segment.
        timestep_fs : float
            Integration timestep in femtoseconds.
        """
        from polyzymd.simulation.progress import (
            SegmentStatus,
            SimulationStatus,
            load_progress,
            save_progress,
        )

        progress = load_progress(self._working_dir)
        if progress is None:
            return

        # Find and update the existing RUNNING record
        for seg in progress.segments:
            if seg.index == self._segment_index and seg.status == SegmentStatus.RUNNING:
                seg.steps_completed = steps_done
                seg.duration_ns = (steps_done * timestep_fs) / 1e6
                break
        else:
            # No RUNNING record found — shouldn't happen, but be safe
            return

        progress.status = SimulationStatus.RUNNING
        save_progress(self._working_dir, progress)
        LOGGER.debug(f"Updated running progress: segment {self._segment_index}, {steps_done} steps")

    def _update_progress_completed(
        self,
        total_steps: int,
        num_samples: int,
        duration_ns: float,
        timestep_fs: float,
    ) -> None:
        """Update progress file after successful segment completion.

        Parameters
        ----------
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
            _update_or_append_segment,
            load_progress,
            save_progress,
        )

        progress = load_progress(self._working_dir)
        if progress is None:
            LOGGER.warning("No progress file found — skipping progress update")
            return

        record = SegmentRecord(
            index=self._segment_index,
            steps_completed=total_steps,
            steps_requested=total_steps,
            samples_written=num_samples,
            status=SegmentStatus.COMPLETED,
            duration_ns=duration_ns,
        )
        from polyzymd.simulation.progress import _now_iso

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
        steps_done: int,
        total_steps: int,
        duration_ns: float,
        timestep_fs: float,
    ) -> None:
        """Update progress file after segment interruption.

        Parameters
        ----------
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
            index=self._segment_index,
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
            f"in segment {self._segment_index}"
        )

    def run_segment(
        self,
        duration_ns: float,
        num_samples: int = 250,
        timestep_fs: float = 2.0,
        *,
        report_interval: int,
        checkpoint_interval_s: float,
    ) -> Dict[str, Any]:
        """Run the continuation segment.

        Runs the simulation for the specified duration, saving trajectory
        frames at regular intervals. On completion, updates the progress
        tracker. On interruption (SIGUSR1/SIGTERM), saves interrupted state,
        updates the progress tracker, and raises ``GracefulExit``.

        Parameters
        ----------
        duration_ns : float
            Duration of this segment in nanoseconds.
        num_samples : int
            Number of trajectory frames to save.
        timestep_fs : float
            Time step in femtoseconds.
        report_interval : int
            Explicit reporter interval in steps.
        checkpoint_interval_s : float
            Wall-time interval in seconds between portable restart
            checkpoints. Also controls how frequently the loop checks
            for SLURM preemption signals.

        Returns
        -------
        dict
            Dictionary with segment results.
        """
        if self._system is None or self._topology is None:
            raise RuntimeError("State not loaded. Call load_previous_state first.")

        LOGGER.info(
            f"Starting segment {self._segment_index}: {duration_ns} ns, {num_samples} frames"
        )

        # Update parameters for this segment
        if self._param_dict:
            integ_values = self._param_dict["__values__"]["integ_params"]["__values__"]
            integ_values["total_time"] = {
                "__class__": "Quantity",
                "__values__": {"value": duration_ns, "unit": "nanosecond"},
            }
            integ_values["num_samples"] = num_samples

        # Add barostat if needed
        self._add_barostat_if_needed()

        # Create integrator and simulation
        integrator = self._create_integrator()
        _, _, _, Simulation, _ = _get_openmm_app_classes()
        from polyzymd.simulation.platform import resolve_platform

        selection = resolve_platform(
            self._platform_name,
            precision=self._platform_precision,
            device_index=self._platform_device_index,
        )
        self._simulation = Simulation(
            self._topology,
            self._system,
            integrator,
            selection.platform,
            selection.properties,
        )

        # Load state from previous segment
        paths = self._get_previous_paths()
        if self._use_checkpoint_recovery:
            # Cases 4/5: Only binary checkpoint available (non-portable).
            chk_path = paths["checkpoint"]
            LOGGER.info(f"Recovering from interrupted segment via checkpoint: {chk_path}")
            self._simulation.loadCheckpoint(str(chk_path))
        else:
            # Cases 1/2/3: Portable state XML (normal, interrupted, or restart).
            LOGGER.info(f"Loading state from {paths['state']}")
            self._simulation.loadState(str(paths["state"]))

        # Create output directory
        output_dir = self._working_dir / f"production_{self._segment_index}"
        output_dir.mkdir(exist_ok=True)

        # Save system XML early so it exists on disk even if the segment is
        # hard-killed (SIGKILL / OOM).  This is required for checkpoint-based
        # recovery: loadCheckpoint() needs a matching System object.  The file
        # is overwritten at segment completion by _save_final_state().
        system_xml_path = output_dir / f"production_{self._segment_index}_system.xml"
        XmlSerializer = _get_xml_serializer()
        with open(system_xml_path, "w") as f:
            f.write(XmlSerializer.serialize(self._system))
        LOGGER.info(f"Saved initial system to {system_xml_path}")

        # Calculate total steps
        total_steps = int(duration_ns * 1e6 / timestep_fs)

        if report_interval <= 0:
            raise ValueError("report_interval must be a positive integer")
        if checkpoint_interval_s <= 0:
            raise ValueError("checkpoint_interval_s must be positive")
        seg_report_interval = report_interval

        # Setup reporters
        self._setup_reporters(seg_report_interval, output_dir)

        # Save parameters for this segment
        if self._param_dict:
            param_path = output_dir / f"production_{self._segment_index}_parameters.json"
            with open(param_path, "w") as f:
                json.dump(self._param_dict, f, indent=2)

        # Install signal handlers for graceful shutdown (SIGUSR1 / SIGTERM)
        from polyzymd.simulation.signals import (
            GracefulExit,
            get_interrupt_signal,
            install_handlers,
            interrupted_state_save_exceptions,
            is_interrupted,
            save_interrupted_state,
            save_restart_checkpoint,
        )

        install_handlers()

        # Mark this segment as RUNNING in progress.json so that
        # check-progress can distinguish actively running simulations
        # from interrupted ones.
        self._write_segment_started(total_steps)

        # Run simulation with adaptive sub-chunks for interrupt responsiveness
        # and periodic wall-time restart checkpoints for preemption resilience.
        LOGGER.info(f"Running {total_steps} steps...")
        steps_done = 0
        import time as _time

        from polyzymd.simulation.progress import PROGRESS_UPDATE_INTERVAL_SECONDS

        _last_progress_write = _time.monotonic()
        _last_checkpoint_write = _time.monotonic()
        _loop_start = _time.monotonic()

        # Adaptive sub-chunk sizing: start with seg_report_interval (the
        # original chunk_size).  After the first checkpoint interval elapses,
        # measure actual steps/second and adapt sub_chunk to target
        # checkpoint_interval / 4 seconds (~15s worth of steps).  This
        # ensures ~4 interrupt checks per checkpoint interval regardless of
        # system size or hardware speed.
        sub_chunk = min(seg_report_interval, total_steps)
        _adapted = False

        try:
            while steps_done < total_steps:
                remaining = total_steps - steps_done
                this_chunk = min(sub_chunk, remaining)
                self._simulation.step(this_chunk)
                steps_done += this_chunk

                _now = _time.monotonic()

                # Adaptive sub-chunk calibration (once, after first interval)
                if not _adapted and (_now - _loop_start) >= checkpoint_interval_s:
                    elapsed = _now - _loop_start
                    steps_per_sec = steps_done / elapsed if elapsed > 0 else 1.0
                    # Target sub-chunk duration = checkpoint_interval / 4
                    target_seconds = checkpoint_interval_s / 4.0
                    new_sub_chunk = max(10, int(steps_per_sec * target_seconds))
                    # Sub-chunk must be a divisor-friendly size relative to
                    # seg_report_interval to avoid misaligned reporter writes.
                    if new_sub_chunk >= seg_report_interval:
                        new_sub_chunk = seg_report_interval
                    else:
                        # Find largest divisor of seg_report_interval <= new_sub_chunk
                        best = new_sub_chunk
                        for candidate in [
                            seg_report_interval // k
                            for k in range(1, seg_report_interval // max(1, new_sub_chunk) + 2)
                        ]:
                            if candidate <= new_sub_chunk and seg_report_interval % candidate == 0:
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
                        steps_done=steps_done,
                        timestep_fs=timestep_fs,
                    )
                    _last_progress_write = _now

                # Wall-time restart checkpoint
                if (
                    (_now - _last_checkpoint_write) >= checkpoint_interval_s
                    and steps_done < total_steps  # skip if we're about to finish
                ):
                    save_restart_checkpoint(
                        simulation=self._simulation,
                        output_dir=output_dir,
                    )
                    _last_checkpoint_write = _now

                if is_interrupted():
                    LOGGER.warning(f"Interrupt detected at step {steps_done}/{total_steps}")
                    save_interrupted_state(
                        simulation=self._simulation,
                        output_dir=output_dir,
                        segment_index=self._segment_index,
                        steps_completed=steps_done,
                        total_steps=total_steps,
                    )
                    # Update progress before raising
                    self._update_progress_interrupted(
                        steps_done=steps_done,
                        total_steps=total_steps,
                        duration_ns=duration_ns,
                        timestep_fs=timestep_fs,
                    )
                    raise GracefulExit(
                        signal_number=get_interrupt_signal(), steps_completed=steps_done
                    )
        except GracefulExit:
            raise  # Re-raise so caller can set exit code
        except interrupted_state_save_exceptions():
            # On unexpected crash, still try to save interrupted state
            try:
                save_interrupted_state(
                    simulation=self._simulation,
                    output_dir=output_dir,
                    segment_index=self._segment_index,
                    steps_completed=steps_done,
                    total_steps=total_steps,
                )
            except interrupted_state_save_exceptions() as save_exc:
                LOGGER.exception(
                    "Failed to save interrupted state after crash: %s",
                    save_exc,
                )
            raise

        # Save final state
        self._save_final_state(output_dir)

        # Update progress tracker
        self._update_progress_completed(
            total_steps=total_steps,
            num_samples=num_samples,
            duration_ns=duration_ns,
            timestep_fs=timestep_fs,
        )

        results = {
            "segment_index": self._segment_index,
            "duration_ns": duration_ns,
            "total_steps": total_steps,
            "num_samples": num_samples,
            "output_dir": str(output_dir),
        }

        LOGGER.info(f"Segment {self._segment_index} completed successfully")
        return results
