"""
Continuation manager for resuming MD simulations from checkpoints.

This module handles loading simulation state from previous segments
and continuing production simulations for daisy-chain workflows.
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Union

import openmm
from openmm import XmlSerializer
from openmm.app import (
    CheckpointReporter,
    DCDReporter,
    PDBFile,
    Simulation,
    StateDataReporter,
)
from openmm.unit import Quantity
from openmm import unit as u

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
    and continuing the simulation, primarily for daisy-chain workflows
    on HPC clusters.

    Example:
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

        Args:
            working_dir: Working directory containing simulation outputs.
            segment_index: Current segment index (1-based).
        """
        self._working_dir = Path(working_dir)
        self._segment_index = segment_index
        self._prev_segment = segment_index - 1

        # State
        self._system: Optional[openmm.System] = None
        self._topology: Optional[Any] = None
        self._simulation: Optional[Simulation] = None
        self._param_dict: Optional[Dict[str, Any]] = None

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

        Returns:
            Path to the solvated PDB file.

        Raises:
            FileNotFoundError: If no suitable PDB file is found.
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

        Returns:
            Dictionary with paths to state, system, and parameter files.
        """
        prev_dir = self._working_dir / f"production_{self._prev_segment}"

        return {
            "state": prev_dir / f"production_{self._prev_segment}_state.xml",
            "system": prev_dir / f"production_{self._prev_segment}_system.xml",
            "params": prev_dir / f"production_{self._prev_segment}_parameters.json",
            "checkpoint": prev_dir / f"production_{self._prev_segment}_checkpoint.chk",
        }

    def load_previous_state(self) -> None:
        """Load state from the previous segment.

        This loads the system, topology, and parameters from the previous
        production segment.

        Raises:
            FileNotFoundError: If required files are missing.
        """
        LOGGER.info(f"Loading state from segment {self._prev_segment}")

        paths = self._get_previous_paths()

        # Check that required files exist
        for name, path in paths.items():
            if name != "checkpoint" and not path.exists():
                raise FileNotFoundError(f"Required file not found: {path}")

        # Load system
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

        LOGGER.info("Previous state loaded successfully")

    def _create_integrator(self) -> openmm.Integrator:
        """Create an integrator from the parameter dictionary.

        Returns:
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
        total_steps: int,
        num_samples: int,
        output_dir: Path,
    ) -> None:
        """Setup reporters for the simulation.

        Args:
            total_steps: Total steps for this segment.
            num_samples: Number of trajectory frames to save.
            output_dir: Output directory for this segment.
        """
        if self._simulation is None:
            raise RuntimeError("Simulation not created")

        report_interval = max(1, total_steps // num_samples)

        # Trajectory reporter - write unwrapped coordinates for proper post-processing
        # (enforcePeriodicBox=False prevents OpenMM from wrapping molecules,
        # which would break them across periodic boundaries)
        traj_path = output_dir / f"production_{self._segment_index}_trajectory.dcd"
        self._simulation.reporters.append(
            DCDReporter(str(traj_path), report_interval, enforcePeriodicBox=False)
        )

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

        Args:
            output_dir: Output directory for this segment.
        """
        if self._simulation is None:
            raise RuntimeError("Simulation not available")

        # Save state
        state_path = output_dir / f"production_{self._segment_index}_state.xml"
        state = self._simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            getForces=True,
            getEnergy=True,
            getParameters=True,
            enforcePeriodicBox=True,
        )

        with open(state_path, "w") as f:
            f.write(XmlSerializer.serialize(state))

        # Save system
        system_path = output_dir / f"production_{self._segment_index}_system.xml"
        with open(system_path, "w") as f:
            f.write(XmlSerializer.serialize(self._simulation.system))

        LOGGER.info(f"Saved final state to {state_path}")
        LOGGER.info(f"Saved system to {system_path}")

    def run_segment(
        self,
        duration_ns: float,
        num_samples: int = 250,
        timestep_fs: float = 2.0,
    ) -> Dict[str, Any]:
        """Run the continuation segment.

        Args:
            duration_ns: Duration of this segment in nanoseconds.
            num_samples: Number of trajectory frames to save.
            timestep_fs: Time step in femtoseconds.

        Returns:
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
            self._param_dict["__values__"]["integ_params"]["__values__"]["num_samples"] = (
                num_samples
            )

        # Add barostat if needed
        self._add_barostat_if_needed()

        # Create integrator and simulation
        integrator = self._create_integrator()
        self._simulation = Simulation(self._topology, self._system, integrator)

        # Load state from previous segment
        paths = self._get_previous_paths()
        LOGGER.info(f"Loading state from {paths['state']}")
        self._simulation.loadState(str(paths["state"]))

        # Create output directory
        output_dir = self._working_dir / f"production_{self._segment_index}"
        output_dir.mkdir(exist_ok=True)

        # Calculate total steps
        total_steps = int(duration_ns * 1e6 / timestep_fs)

        # Setup reporters
        self._setup_reporters(total_steps, num_samples, output_dir)

        # Save parameters for this segment
        if self._param_dict:
            param_path = output_dir / f"production_{self._segment_index}_parameters.json"
            with open(param_path, "w") as f:
                json.dump(self._param_dict, f, indent=2)

        # Run simulation
        LOGGER.info(f"Running {total_steps} steps...")
        self._simulation.step(total_steps)

        # Save final state
        self._save_final_state(output_dir)

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
    """Main entry point for continuation script.

    Returns:
        Exit code (0 for success, 1 for failure).
    """
    import argparse

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

    except Exception as e:
        LOGGER.error(f"Error during simulation: {e}")
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
