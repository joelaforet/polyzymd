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
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Union

import openmm
from openmm import XmlSerializer
from openmm import unit as u
from openmm.app import (
    CheckpointReporter,
    DCDReporter,
    PDBFile,
    Simulation,
    StateDataReporter,
)
from openmm.unit import Quantity

LOGGER = logging.getLogger(__name__)


def quantity_from_dict(qdict: Dict[str, Any]) -> Quantity:
    """Convert serialized quantity dictionary back to OpenMM Quantity.

    Args:
        qdict: Dictionary with __values__ containing value and unit.

    Returns:
        OpenMM Quantity with appropriate units.
    """
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
        patterns = [
            "*solvated*.pdb",
            "*_solvated.pdb",
            "solvated_*.pdb",
            "equilibration/*_topology.pdb",
            "production_0/*_topology.pdb",
        ]

        for pattern in patterns:
            pdb_files = list(self._working_dir.glob(pattern))
            if pdb_files:
                return pdb_files[0]

        # Fallback to any PDB
        pdb_files = list(self._working_dir.glob("**/*.pdb"))
        if pdb_files:
            return pdb_files[0]

        raise FileNotFoundError(f"Could not find solvated PDB file in {self._working_dir}")

    def _get_previous_paths(self) -> Dict[str, Path]:
        """Get paths to files from the previous segment.

        Returns
        -------
        dict
            Dictionary with paths to state, system, and parameter files.
            If the previous segment completed normally, ``state`` and ``system``
            point to ``production_N_state.xml`` / ``production_N_system.xml``.
            If the previous segment was interrupted (those files don't exist),
            ``system`` falls back to ``interrupted_system.xml`` and
            ``use_checkpoint`` is set to ``True`` to signal that
            ``load_previous_state`` should use ``loadCheckpoint()`` instead
            of ``loadState()``.
        """
        prev_dir = self._working_dir / f"production_{self._prev_segment}"

        state_path = prev_dir / f"production_{self._prev_segment}_state.xml"
        system_path = prev_dir / f"production_{self._prev_segment}_system.xml"
        checkpoint_path = prev_dir / f"production_{self._prev_segment}_checkpoint.chk"
        params_path = prev_dir / f"production_{self._prev_segment}_parameters.json"

        # Check if previous segment was interrupted (normal state files missing)
        interrupted_system = prev_dir / "interrupted_system.xml"
        use_checkpoint = False

        if not state_path.exists() and interrupted_system.exists():
            LOGGER.info(
                f"Previous segment {self._prev_segment} was interrupted — "
                f"recovering from checkpoint + interrupted system XML"
            )
            system_path = interrupted_system
            use_checkpoint = True

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
        production segment. If the previous segment was interrupted (no
        ``production_N_state.xml``), it loads the system from
        ``interrupted_system.xml`` and marks the checkpoint for recovery
        via ``loadCheckpoint()`` in ``run_segment()``.

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
            if name == "checkpoint":
                continue
            if name == "state" and use_checkpoint:
                # State XML doesn't exist for interrupted segments; we'll
                # use the checkpoint instead in run_segment()
                continue
            if not path.exists():
                raise FileNotFoundError(f"Required file not found: {path}")

        # Load system (either normal or interrupted system XML)
        LOGGER.info(f"Loading system from {paths['system']}")
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

        progress.segments.append(record)

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

        progress.segments.append(record)
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
        report_interval: int | None = None,
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
            Number of trajectory frames to save. Ignored when
            ``report_interval`` is provided.
        timestep_fs : float
            Time step in femtoseconds.
        report_interval : int or None
            Fixed reporter interval in steps. When provided, this
            overrides the per-segment ``total_steps // num_samples``
            calculation to keep frame spacing uniform across segments.

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
            self._param_dict["__values__"]["integ_params"]["__values__"]["total_time"] = {
                "__class__": "Quantity",
                "__values__": {"value": duration_ns, "unit": "nanosecond"},
            }
            self._param_dict["__values__"]["integ_params"]["__values__"][
                "num_samples"
            ] = num_samples

        # Add barostat if needed
        self._add_barostat_if_needed()

        # Create integrator and simulation
        integrator = self._create_integrator()
        self._simulation = Simulation(self._topology, self._system, integrator)

        # Load state from previous segment
        paths = self._get_previous_paths()
        if self._use_checkpoint_recovery:
            # Previous segment was interrupted — no state XML exists.
            # Use loadCheckpoint() with the reporter-interval checkpoint,
            # which is consistent with the partial trajectory DCD.
            chk_path = paths["checkpoint"]
            LOGGER.info(f"Recovering from interrupted segment via checkpoint: {chk_path}")
            self._simulation.loadCheckpoint(str(chk_path))
        else:
            LOGGER.info(f"Loading state from {paths['state']}")
            self._simulation.loadState(str(paths["state"]))

        # Create output directory
        output_dir = self._working_dir / f"production_{self._segment_index}"
        output_dir.mkdir(exist_ok=True)

        # Calculate total steps
        total_steps = int(duration_ns * 1e6 / timestep_fs)

        # Determine report interval: prefer the fixed global value if given,
        # otherwise fall back to per-segment calculation (legacy callers).
        if report_interval is not None:
            seg_report_interval = report_interval
        else:
            seg_report_interval = max(1, total_steps // num_samples)

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
            is_interrupted,
            save_interrupted_state,
        )

        install_handlers()

        # Run simulation in chunks so we can check for interrupt signals
        LOGGER.info(f"Running {total_steps} steps...")
        chunk_size = min(seg_report_interval, total_steps)
        steps_done = 0
        try:
            while steps_done < total_steps:
                remaining = total_steps - steps_done
                this_chunk = min(chunk_size, remaining)
                self._simulation.step(this_chunk)
                steps_done += this_chunk
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
        except Exception:
            # On unexpected crash, still try to save interrupted state
            try:
                save_interrupted_state(
                    simulation=self._simulation,
                    output_dir=output_dir,
                    segment_index=self._segment_index,
                    steps_completed=steps_done,
                    total_steps=total_steps,
                )
            except Exception:
                LOGGER.error("Failed to save interrupted state after crash")
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


def main() -> int:
    """Legacy entry point for continuation script.

    .. deprecated::
        Use ``polyzymd run-segment`` CLI command instead, which provides
        unified initial/continuation logic with progress tracking.

    Returns
    -------
    int
        Exit code (0 for success, 1 for failure, 99 for graceful interrupt).
    """
    import argparse

    from polyzymd.simulation.signals import EXIT_CODE_INTERRUPTED, GracefulExit

    parser = argparse.ArgumentParser(description="Continue MD simulation from previous segment")
    parser.add_argument(
        "-s",
        "--segment_index",
        type=int,
        required=True,
        help="Current segment index (1-based)",
    )
    parser.add_argument(
        "-w",
        "--working_dir",
        type=str,
        required=True,
        help="Working directory path",
    )
    parser.add_argument(
        "-t",
        "--segment_time",
        type=float,
        required=True,
        help="Time for this segment in nanoseconds",
    )
    parser.add_argument(
        "-n",
        "--num_samples",
        type=int,
        default=250,
        help="Number of frames to save for this segment",
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(
                Path(args.working_dir) / f"simulation_status_segment_{args.segment_index}.log"
            ),
            logging.StreamHandler(sys.stdout),
        ],
    )

    try:
        manager = ContinuationManager(
            working_dir=args.working_dir,
            segment_index=args.segment_index,
        )
        manager.load_previous_state()
        manager.run_segment(
            duration_ns=args.segment_time,
            num_samples=args.num_samples,
        )
        return 0

    except GracefulExit as e:
        LOGGER.warning(f"Graceful exit: {e}")
        return EXIT_CODE_INTERRUPTED

    except Exception as e:
        LOGGER.error(f"Error during simulation: {e}")
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
