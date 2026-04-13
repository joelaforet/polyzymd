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

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.

        Returns
        -------
        TrajectoryLayout
            Resolved DCD/PDB layout for downstream analyses.
        """
        _ = replicate

        prod_dir_re = re.compile(r"^production_(\d+)$")
        production_dirs: list[tuple[int, Path]] = []
        for child in working_dir.iterdir():
            if child.is_dir():
                match = prod_dir_re.match(child.name)
                if match:
                    production_dirs.append((int(match.group(1)), child))

        production_dirs.sort(key=lambda item: item[0])

        trajectory_paths: list[Path] = []
        topology_path: Path | None = None
        for _, prod_dir in production_dirs:
            trajectory_paths.extend(sorted(prod_dir.glob("*_trajectory.dcd")))
            if topology_path is None:
                topo_candidates = sorted(prod_dir.glob("*_topology.pdb"))
                if topo_candidates:
                    topology_path = topo_candidates[0]

        if topology_path is None:
            fallback = working_dir / "solvated_system.pdb"
            if fallback.exists():
                topology_path = fallback

        return TrajectoryLayout(
            topology_path=topology_path,
            trajectory_paths=trajectory_paths,
            trajectory_format="dcd",
            topology_format="pdb",
        )
