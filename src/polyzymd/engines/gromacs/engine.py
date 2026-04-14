"""GROMACS engine implementation for PolyzyMD."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any, ClassVar

from polyzymd.engines.base import EngineSubmitRequest, SimulationEngine, TrajectoryLayout
from polyzymd.simulation.progress import SimulationProgress
from polyzymd.workflow.slurm import SlurmConfig

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

    def _generate_prefix(self) -> str:
        """Generate filename prefix from simulation config.

        Replicates the prefix logic from ``GromacsExporter._generate_prefix``
        so that file paths are consistent between build and submit flows.

        Returns
        -------
        str
            System prefix (e.g. ``"CALB_SBMA-EGMA"``).
        """
        parts: list[str] = []

        enzyme = getattr(self._config, "enzyme", None)
        enzyme_name = getattr(enzyme, "name", None)
        if isinstance(enzyme_name, str) and enzyme_name:
            parts.append(enzyme_name)

        polymers = getattr(self._config, "polymers", None)
        polymers_enabled = getattr(polymers, "enabled", False)
        polymer_prefix = getattr(polymers, "type_prefix", None)
        if (
            isinstance(polymers_enabled, bool)
            and polymers_enabled
            and isinstance(polymer_prefix, str)
        ):
            parts.append(polymer_prefix)

        return "_".join(parts) if parts else "system"

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

        prefix = self._generate_prefix()
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
        prefix = self._generate_prefix()
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

        effective_slurm = self._resolve_slurm_config(request.slurm_config)
        effective_mdrun_flags = self._resolve_mdrun_flags(effective_slurm)

        generator = GromacsSlurmScriptGenerator(
            slurm_config=effective_slurm,
            pixi_env=pixi_env,
            gmx_binary=self._gmx_binary,
            grompp_flags=self._config.gromacs.grompp_flags,
            mdrun_flags=effective_mdrun_flags,
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

    def _resolve_slurm_config(self, base: SlurmConfig) -> SlurmConfig:
        """Override base SLURM config with GROMACS-specific hardware settings.

        Parameters
        ----------
        base : SlurmConfig
            SLURM config from preset or CLI.

        Returns
        -------
        SlurmConfig
            Config with GROMACS hardware overrides applied.
        """
        from dataclasses import replace

        gromacs_cfg = self._config.gromacs
        overrides: dict[str, Any] = {
            "ntasks": gromacs_cfg.ntmpi,
            "cpus_per_task": gromacs_cfg.ntomp,
            "memory": gromacs_cfg.memory,
        }
        if gromacs_cfg.gpu:
            overrides["gpus"] = max(base.gpus, 1)
        else:
            overrides["gpus"] = 0
        return replace(base, **overrides)

    def _resolve_mdrun_flags(self, effective_slurm: SlurmConfig) -> str:
        """Compose final mdrun flags from config + hardware settings.

        Appends ``-ntmpi`` and ``-ntomp`` only if the user has not already
        specified them in ``mdrun_flags``.

        The ``-ntmpi`` flag is only appended for thread-MPI builds (binary
        name ``gmx``). Real-MPI builds (``gmx_mpi``, ``gmx_mpi_d``, etc.)
        do **not** support ``-ntmpi`` — MPI rank count is controlled by the
        MPI launcher (``mpirun``/``srun``).

        Parameters
        ----------
        effective_slurm : SlurmConfig
            Resolved SLURM config with hardware overrides.

        Returns
        -------
        str
            Complete mdrun flags string.
        """
        import shlex
        from pathlib import Path

        raw = self._config.gromacs.mdrun_flags
        tokens = shlex.split(raw) if raw else []
        token_set = set(tokens)

        # Detect whether the binary is an MPI build (gmx_mpi, gmx_mpi_d, ...)
        # Thread-MPI builds use plain "gmx" or "gmx_d"; only they accept -ntmpi.
        binary_name = Path(self._gmx_binary).name
        is_mpi_build = "_mpi" in binary_name

        extras: list[str] = []
        if not is_mpi_build and "-ntmpi" not in token_set:
            extras.append(f"-ntmpi {effective_slurm.ntasks}")
        if "-ntomp" not in token_set:
            extras.append(f"-ntomp {effective_slurm.cpus_per_task}")

        parts = [raw] + extras
        return " ".join(part for part in parts if part).strip()

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

        Topology search order:

        1. ``<prefix>.pdb`` (from system name)
        2. Any ``*.pdb`` (sorted, first match)
        3. ``<prefix>.gro``
        4. Any ``*.gro`` (sorted, first match)

        Trajectory search order:

        1. ``prod.xtc`` (standard GROMACS production output)
        2. Any ``*.xtc`` (sorted)

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index (unused, kept for interface parity).

        Returns
        -------
        TrajectoryLayout
            Resolved XTC layout with PDB topology preferred over GRO.
        """
        _ = replicate

        prod_xtc = working_dir / "prod.xtc"
        if prod_xtc.exists():
            trajectory_paths = [prod_xtc]
        else:
            trajectory_paths = sorted(working_dir.glob("*.xtc"))

        topology_path: Path | None = None
        topology_format = "pdb"

        prefix = self._generate_prefix()
        if prefix:
            named_pdb = working_dir / f"{prefix}.pdb"
            if named_pdb.exists():
                topology_path = named_pdb

        if topology_path is None:
            pdb_candidates = sorted(working_dir.glob("*.pdb"))
            if pdb_candidates:
                topology_path = pdb_candidates[0]

        if topology_path is None and prefix:
            named_gro = working_dir / f"{prefix}.gro"
            if named_gro.exists():
                topology_path = named_gro
                topology_format = "gro"

        if topology_path is None:
            gro_candidates = sorted(working_dir.glob("*.gro"))
            if gro_candidates:
                topology_path = gro_candidates[0]
                topology_format = "gro"

        return TrajectoryLayout(
            topology_path=topology_path,
            trajectory_paths=trajectory_paths,
            trajectory_format="xtc",
            topology_format=topology_format,
        )
