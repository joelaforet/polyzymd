"""GROMACS engine implementation for PolyzyMD."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any, ClassVar

from polyzymd.engines.base import EngineSubmitRequest, SimulationEngine, TrajectoryLayout
from polyzymd.simulation.progress import SimulationProgress

from .binary import resolve_gromacs_binary
from .progress import load_or_scan_gromacs_progress
from .slurm import GromacsSlurmScriptGenerator


class GromacsEngine(SimulationEngine):
    """GROMACS execution adapter for local and scheduler workflows."""

    name: ClassVar[str] = "gromacs"

    def __init__(self, config: object, gmx_binary: str = "gmx"):
        """Initialize a GROMACS engine adapter.

        Parameters
        ----------
        config : object
            Simulation configuration object.
        gmx_binary : str, optional
            Resolved GROMACS executable.
        """
        self._config = config
        self._gmx_binary = gmx_binary

    @classmethod
    def from_config(cls, config: object) -> GromacsEngine:
        """Create a GROMACS engine from simulation config.

        Parameters
        ----------
        config : object
            Simulation configuration object.

        Returns
        -------
        GromacsEngine
            Configured GROMACS engine instance.
        """
        gmx_binary = resolve_gromacs_binary(config=config)
        return cls(config=config, gmx_binary=gmx_binary)

    def run_local(self, replicate: int, working_dir: Path, skip_build: bool = False) -> None:
        """Run local GROMACS workflow with exported files.

        Parameters
        ----------
        replicate : int
            Replicate index.
        working_dir : Path
            Working directory with GROMACS input files.
        skip_build : bool, optional
            Build-skip placeholder for interface parity.
        """
        from polyzymd.exporters.gromacs import GromacsRunner

        _ = replicate
        _ = skip_build

        prefix = self._config.generate_system_name()
        eq_mdps = sorted(path.name for path in working_dir.glob("eq_*.mdp"))

        runner = GromacsRunner(
            working_dir=working_dir,
            prefix=prefix,
            equilibration_mdps=eq_mdps,
            gmx_command=self._gmx_binary,
        )
        runner.run_full_workflow()

    def prepare_submission(self, request: EngineSubmitRequest) -> Path:
        """Prepare scheduler artifacts for a GROMACS job.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.

        Returns
        -------
        Path
            Path to the generated SLURM script.
        """
        if request.slurm_config is None:
            raise ValueError("GROMACS submission requires slurm_config")

        pixi_env = str(request.extra.get("pixi_env", "cuda-12-4"))
        skip_build = bool(request.extra.get("skip_build", False))

        request.working_dir.mkdir(parents=True, exist_ok=True)

        # Build and export GROMACS inputs when not already present
        prefix = self._config.generate_system_name()
        top_path = request.working_dir / f"{prefix}.top"
        gro_path = request.working_dir / f"{prefix}.gro"
        em_path = request.working_dir / "em.mdp"
        prod_path = request.working_dir / "prod.mdp"

        inputs_exist = (
            top_path.exists() and gro_path.exists() and em_path.exists() and prod_path.exists()
        )

        if not inputs_exist and skip_build:
            raise FileNotFoundError(
                "skip_build=True but required GROMACS inputs are missing "
                f"in {request.working_dir}. Expected: {top_path.name}, {gro_path.name}, "
                f"{em_path.name}, {prod_path.name}."
            )

        if not inputs_exist:
            from polyzymd.builders.system_builder import SystemBuilder
            from polyzymd.exporters.gromacs import GromacsExporter

            builder = SystemBuilder.from_config(self._config)
            interchange = builder.build_from_config(
                config=self._config,
                working_dir=request.working_dir,
                polymer_seed=request.replicate,
            )
            component_info = builder.get_component_info()
            exporter = GromacsExporter(interchange, self._config, component_info=component_info)
            exporter.export(
                output_dir=request.working_dir,
                prefix=prefix,
                gmx_command=self._gmx_binary,
            )

        eq_mdps = sorted(path.name for path in request.working_dir.glob("eq_*.mdp"))

        script_dir = request.working_dir / "daisy_chain_scripts"
        script_dir.mkdir(parents=True, exist_ok=True)
        script_path = script_dir / f"run_rep{request.replicate}.sh"

        generator = GromacsSlurmScriptGenerator(
            slurm_config=request.slurm_config,
            pixi_env=pixi_env,
            gmx_binary=self._gmx_binary,
            grompp_flags=self._config.gromacs.grompp_flags,
            mdrun_flags=self._config.gromacs.mdrun_flags,
            module_load=self._config.gromacs.module_load,
        )
        script = generator.generate_job_script(
            config_path=str(request.config_path),
            replicate=request.replicate,
            working_dir=str(request.working_dir),
            system_prefix=prefix,
            equilibration_mdps=eq_mdps,
            job_name=request.job_name,
        )
        generator.save_script(script, script_path)
        return script_path

    def submit(self, request: EngineSubmitRequest) -> Any:
        """Submit GROMACS jobs to scheduler.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.

        Returns
        -------
        Any
            Submission metadata with script path and optional SLURM job id.
        """
        script_path = self.prepare_submission(request)

        if shutil.which("sbatch") is None:
            return {
                "submitted": False,
                "script_path": script_path,
                "reason": "sbatch_not_available",
            }

        command = ["sbatch", str(script_path)]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"sbatch submission failed: {result.stderr.strip()}")

        return {
            "submitted": True,
            "script_path": script_path,
            "stdout": result.stdout.strip(),
            "stderr": result.stderr.strip(),
            "returncode": result.returncode,
        }

    def load_or_scan_progress(self, working_dir: Path, replicate: int) -> SimulationProgress:
        """Load or reconstruct GROMACS progress state.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.

        Returns
        -------
        SimulationProgress
            Current progress model for the replicate.
        """
        prod = self._config.simulation_phases.production
        total_steps = int(prod.duration * 1e6 / prod.time_step)

        return load_or_scan_gromacs_progress(
            working_dir=working_dir,
            config_path="",
            replicate=replicate,
            total_steps=total_steps,
            total_samples=prod.samples,
            timestep_fs=prod.time_step,
        )

    def resolve_trajectory_layout(self, working_dir: Path, replicate: int) -> TrajectoryLayout:
        """Resolve GROMACS trajectory layout for downstream analyses.

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index.

        Returns
        -------
        TrajectoryLayout
            Resolved XTC layout with PDB topology fallback to GRO.
        """
        _ = replicate

        trajectory_paths = sorted(working_dir.glob("*.xtc"))

        pdb_candidates = sorted(working_dir.glob("*.pdb"))
        gro_candidates = sorted(working_dir.glob("*.gro"))

        topology_path: Path | None = None
        topology_format = "pdb"
        if pdb_candidates:
            topology_path = pdb_candidates[0]
        elif gro_candidates:
            topology_path = gro_candidates[0]
            topology_format = "gro"

        return TrajectoryLayout(
            topology_path=topology_path,
            trajectory_paths=trajectory_paths,
            trajectory_format="xtc",
            topology_format=topology_format,
        )
