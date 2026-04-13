"""OpenMM engine adapter over existing PolyzyMD simulation workflow."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, ClassVar

from polyzymd.engines.base import EngineSubmitRequest, SimulationEngine, TrajectoryLayout
from polyzymd.simulation.progress import SimulationProgress


class OpenMMEngine(SimulationEngine):
    """Thin adapter for the OpenMM execution path."""

    name: ClassVar[str] = "openmm"

    def __init__(self, config: object):
        """Initialize the engine adapter.

        Parameters
        ----------
        config : object
            Simulation configuration object.
        """
        self._config = config

    @classmethod
    def from_config(cls, config: object) -> OpenMMEngine:
        """Create an OpenMM engine instance from configuration.

        Parameters
        ----------
        config : object
            Simulation configuration object.

        Returns
        -------
        OpenMMEngine
            Configured OpenMM engine.
        """
        return cls(config=config)

    def run_local(self, replicate: int, working_dir: Path, skip_build: bool = False) -> None:
        """Run the OpenMM simulation locally.

        Parameters
        ----------
        replicate : int
            Replicate index.
        working_dir : Path
            Working directory for simulation output.
        skip_build : bool, optional
            Whether to skip build and reuse existing artifacts.
        """
        from polyzymd.cli.main import _run_initial_segment

        prod = self._config.simulation_phases.production
        _run_initial_segment(
            sim_config=self._config,
            working_dir=working_dir,
            replicate=replicate,
            skip_build=skip_build,
            duration_ns=prod.duration,
            num_samples=prod.samples,
            timestep_fs=prod.time_step,
            report_interval=None,
            checkpoint_interval_s=prod.checkpoint_interval,
        )

    def prepare_submission(self, request: EngineSubmitRequest) -> Path:
        """Prepare a SLURM submission script for OpenMM.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.

        Returns
        -------
        Path
            Path to the generated SLURM script.
        """
        from polyzymd.workflow.slurm import SlurmScriptGenerator

        if request.slurm_config is None:
            raise ValueError("OpenMM submission requires slurm_config")

        script_dir = request.working_dir / "daisy_chain_scripts"
        script_dir.mkdir(parents=True, exist_ok=True)
        script_path = script_dir / f"run_rep{request.replicate}.sh"

        generator = SlurmScriptGenerator(config=request.slurm_config)
        script = generator.generate_job_script(
            config_path=str(request.config_path),
            replicate=request.replicate,
            working_dir=str(request.working_dir),
            job_name=request.job_name,
        )
        generator.save_script(script, script_path)
        return script_path

    def submit(self, request: EngineSubmitRequest) -> Any:
        """Submit an OpenMM replicate through the daisy-chain workflow.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.

        Returns
        -------
        Any
            Submission result object from the existing workflow.
        """
        from polyzymd.workflow.daisy_chain import DaisyChainConfig, DaisyChainSubmitter

        if request.slurm_config is None:
            raise ValueError("OpenMM submission requires slurm_config")

        dc_config = DaisyChainConfig.from_simulation_config(
            sim_config=self._config,
            slurm_config=request.slurm_config,
            replicates=[request.replicate],
            output_script_dir=request.working_dir / "daisy_chain_scripts",
            config_path=str(request.config_path),
        )
        submitter = DaisyChainSubmitter(self._config, dc_config)
        return submitter.submit_replicate(request.replicate)

    def load_or_scan_progress(self, working_dir: Path, replicate: int) -> SimulationProgress:
        """Load OpenMM progress from progress.json or filesystem scan.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.

        Returns
        -------
        SimulationProgress
            Current progress state for the replicate.
        """
        from polyzymd.simulation.progress import load_or_scan_progress

        prod = self._config.simulation_phases.production
        total_steps = int(prod.duration * 1e6 / prod.time_step)

        return load_or_scan_progress(
            working_dir=working_dir,
            config_path="",
            total_steps=total_steps,
            total_samples=prod.samples,
            timestep_fs=prod.time_step,
            replicate=replicate,
        )

    def resolve_trajectory_layout(self, working_dir: Path, replicate: int) -> TrajectoryLayout:
        """Resolve OpenMM trajectory and topology paths.

        Replicates the search order used by TrajectoryLoader so that the
        engine can serve as the authoritative file resolver.

        Topology search order:

        1. ``solvated_system.pdb`` in working_dir root
        2. ``production_0/production_0_topology.pdb``
        3. ``production/production_topology.pdb`` (legacy)
        4. Glob ``production_*/*_topology.pdb``
        5. Glob ``production/*_topology.pdb``
        6. Glob ``*.pdb`` in root

        Trajectory search order:

        1. ``production_N/production_N_trajectory.dcd`` (ordered by N)
        2. ``production/production_trajectory.dcd`` (legacy)
        3. Glob ``**/production*trajectory.dcd``

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index (unused, kept for interface parity).

        Returns
        -------
        TrajectoryLayout
            Resolved DCD/PDB layout for downstream analyses.
        """
        _ = replicate

        topology_path = self._find_openmm_topology(working_dir)
        trajectory_paths = self._find_openmm_trajectories(working_dir)

        return TrajectoryLayout(
            topology_path=topology_path,
            trajectory_paths=trajectory_paths,
            trajectory_format="dcd",
            topology_format="pdb",
        )

    @staticmethod
    def _find_openmm_topology(working_dir: Path) -> Path | None:
        """Find topology PDB using the canonical OpenMM search order.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.

        Returns
        -------
        Path or None
            Path to topology file, or None if not found.
        """
        # Primary topology path used by simulation workflow
        candidate = working_dir / "solvated_system.pdb"
        if candidate.exists():
            return candidate

        # Daisy-chain output from the first production segment
        candidate = working_dir / "production_0" / "production_0_topology.pdb"
        if candidate.exists():
            return candidate

        # Legacy single production directory layout
        candidate = working_dir / "production" / "production_topology.pdb"
        if candidate.exists():
            return candidate

        # Fallback glob patterns retained for backward compatibility
        for pattern in [
            "production_*/*_topology.pdb",
            "production/*_topology.pdb",
            "*.pdb",
        ]:
            hits = sorted(working_dir.glob(pattern))
            if hits:
                return hits[0]

        return None

    @staticmethod
    def _find_openmm_trajectories(working_dir: Path) -> list[Path]:
        """Find trajectory DCD files using the canonical OpenMM search order.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.

        Returns
        -------
        list[Path]
            Ordered list of trajectory files.
        """
        # Daisy-chain segmented trajectories ordered by production index
        prod_re = re.compile(r"production_(\d+)[/\\]production_\d+_trajectory\.dcd$")
        segments: dict[int, Path] = {}
        for file_path in working_dir.glob("production_*/production_*_trajectory.dcd"):
            match = prod_re.search(str(file_path))
            if match:
                segments[int(match.group(1))] = file_path

        if segments:
            return [segments[index] for index in sorted(segments.keys())]

        # Legacy single-file production trajectory path
        legacy = working_dir / "production" / "production_trajectory.dcd"
        if legacy.exists():
            return [legacy]

        # Last-resort recursive pattern for historical layouts
        return sorted(working_dir.glob("**/production*trajectory.dcd"))
