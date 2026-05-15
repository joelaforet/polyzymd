"""GROMACS engine implementation for PolyzyMD."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Any, ClassVar

from polyzymd.engines.base import EngineSubmitRequest, SimulationEngine, TrajectoryLayout
from polyzymd.simulation.progress import SimulationProgress
from polyzymd.workflow.slurm import SlurmConfig

from .binary import is_mpi_binary, resolve_gromacs_binary
from .progress import load_or_scan_gromacs_progress
from .slurm import GromacsSlurmScriptGenerator


class GromacsEngine(SimulationEngine):
    """GROMACS execution adapter for local and scheduler workflows."""

    name: ClassVar[str] = "gromacs"
    engine_subdir: ClassVar[str] = "gromacs"

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
    def from_config(cls, config: object, defer_binary: bool = False) -> GromacsEngine:
        """Create a GROMACS engine from simulation config.

        Parameters
        ----------
        config : object
            Simulation configuration object.
        defer_binary : bool, optional
            When True, skip local PATH probing and use the configured
            binary name directly (or ``"gmx"`` fallback). This is used
            by scheduler submission paths where module loading happens
            on compute nodes.

        Returns
        -------
        GromacsEngine
            Configured GROMACS engine instance.

        Raises
        ------
        ValueError
            If GPU mode is enabled but the resolved binary is a real-MPI
            build (which typically lacks CUDA support).
        """
        gromacs_cfg = getattr(config, "gromacs", None)
        gpu = getattr(gromacs_cfg, "gpu", False) if gromacs_cfg else False
        if defer_binary:
            configured = getattr(gromacs_cfg, "gmx_binary", None) if gromacs_cfg else None
            gmx_binary = configured or "gmx"
            if gpu and is_mpi_binary(gmx_binary):
                logging.getLogger(__name__).warning(
                    "GPU mode with a real-MPI binary (%s) requires manual -ntmpi/-ntomp "
                    "flag management via mdrun_flags. Consider using thread-MPI (gmx) "
                    "for single-node GPU workflows.",
                    gmx_binary,
                )
        else:
            gmx_binary = resolve_gromacs_binary(config=config, gpu=gpu)
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
        if inputs_exist and not eq_mdps:
            logging.getLogger(__name__).warning(
                "Core GROMACS inputs found in %s but no equilibration MDPs (eq_*.mdp). "
                "The generated script will skip equilibration and run production from EM output.",
                request.working_dir,
            )

        script_dir = request.working_dir / "daisy_chain_scripts"
        script_dir.mkdir(parents=True, exist_ok=True)
        script_path = script_dir / f"run_rep{request.replicate}.sh"

        effective_slurm = self._resolve_slurm_config(request.slurm_config)
        effective_mdrun_flags = self._resolve_mdrun_flags(effective_slurm)
        mdrun_flags_eq, mdrun_flags_prod = self._resolve_stage_mdrun_flags(effective_slurm)

        generator = GromacsSlurmScriptGenerator(
            slurm_config=effective_slurm,
            pixi_env=pixi_env,
            gmx_binary=self._gmx_binary,
            grompp_flags=self._config.gromacs.grompp_flags,
            mdrun_flags=effective_mdrun_flags,
            mdrun_flags_eq=mdrun_flags_eq,
            mdrun_flags_prod=mdrun_flags_prod,
            command_prefix=self._config.gromacs.command_prefix,
            mpi_launcher_flags=self._config.gromacs.mpi_launcher_flags,
            module_load=self._config.gromacs.module_load,
            env_exports=self._config.gromacs.env_exports,
            setup_commands=self._config.gromacs.setup_commands,
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
            "ntasks": (
                gromacs_cfg.slurm_ntasks
                if gromacs_cfg.slurm_ntasks is not None
                else gromacs_cfg.ntmpi
            ),
            "cpus_per_task": gromacs_cfg.ntomp,
            "memory": gromacs_cfg.memory,
        }
        if gromacs_cfg.gpu:
            overrides["gpus"] = gromacs_cfg.gpus
        else:
            overrides["gpus"] = 0
        return replace(base, **overrides)

    def _resolve_mdrun_flags(
        self, effective_slurm: SlurmConfig, *, base_flags: str | None = None
    ) -> str:
        """Compose final mdrun flags from config + hardware settings.

        Appends ``-ntmpi`` and ``-ntomp`` only if the user has not already
        specified them in ``mdrun_flags``.

        The ``-ntmpi`` flag is only appended for thread-MPI builds (binary
        name ``gmx``). Real-MPI builds (``gmx_mpi``, ``gmx_mpi_d``, etc.)
        do **not** support ``-ntmpi`` — MPI rank count is controlled by the
        MPI launcher (``mpirun``/``srun``).

        When ``gpu`` is enabled in the GROMACS config, the offload flags
        ``-nb gpu``, ``-pme gpu``, ``-bonded gpu``, and ``-update gpu`` are
        appended automatically. Each flag is skipped if the user already
        specified that flag key (e.g., ``-nb cpu`` prevents auto-adding
        ``-nb gpu``).

        Parameters
        ----------
        effective_slurm : SlurmConfig
            Resolved SLURM config with hardware overrides.
        base_flags : str | None, optional
            Override for base flags instead of
            ``self._config.gromacs.mdrun_flags``.

        Returns
        -------
        str
            Complete mdrun flags string.
        """
        import shlex

        raw = base_flags if base_flags is not None else self._config.gromacs.mdrun_flags
        tokens = shlex.split(raw) if raw else []
        token_set = set(tokens)

        mpi_build = is_mpi_binary(self._gmx_binary)

        extras: list[str] = []
        if not mpi_build and "-ntmpi" not in token_set:
            extras.append(f"-ntmpi {self._config.gromacs.ntmpi}")
        if "-ntomp" not in token_set:
            extras.append(f"-ntomp {effective_slurm.cpus_per_task}")

        # GPU offload flags — auto-add when gpu:true, skip if user already
        # specified the flag key (e.g., user has "-nb cpu" → don't add "-nb gpu")
        if self._config.gromacs.gpu:
            _GPU_OFFLOAD_FLAGS = [
                ("-nb", "gpu"),
                ("-pme", "gpu"),
                ("-bonded", "gpu"),
                ("-update", "gpu"),
            ]
            for flag_key, flag_val in _GPU_OFFLOAD_FLAGS:
                if flag_key not in token_set:
                    extras.append(f"{flag_key} {flag_val}")

        parts = [raw] + extras
        return " ".join(part for part in parts if part).strip()

    def _resolve_stage_mdrun_flags(
        self, effective_slurm: SlurmConfig
    ) -> tuple[str | None, str | None]:
        """Resolve equilibration and production flag strings with fallback.

        Parameters
        ----------
        effective_slurm : SlurmConfig
            Resolved SLURM config with hardware overrides.

        Returns
        -------
        tuple[str | None, str | None]
            Stage-specific mdrun flags for equilibration and production.
            ``None`` indicates fallback to global ``mdrun_flags`` in script.
        """
        eq_raw = getattr(self._config.gromacs, "mdrun_flags_equilibration", None)
        prod_raw = getattr(self._config.gromacs, "mdrun_flags_production", None)

        eq_flags = self._resolve_mdrun_flags_for_raw(eq_raw, effective_slurm) if eq_raw else None
        prod_flags = (
            self._resolve_mdrun_flags_for_raw(prod_raw, effective_slurm) if prod_raw else None
        )
        return eq_flags, prod_flags

    def _resolve_mdrun_flags_for_raw(self, raw_flags: str, effective_slurm: SlurmConfig) -> str:
        """Resolve one mdrun flag string against current SLURM settings.

        Parameters
        ----------
        raw_flags : str
            Raw mdrun flag string from configuration.
        effective_slurm : SlurmConfig
            Resolved SLURM config with hardware overrides.

        Returns
        -------
        str
            Fully resolved mdrun flags for one stage.
        """
        return self._resolve_mdrun_flags(effective_slurm, base_flags=raw_flags)

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
        module_load = self._config.gromacs.module_load

        if not module_load and shutil.which("sbatch") is None:
            return {
                "submitted": False,
                "script_path": script_path,
                "reason": "sbatch_not_available",
            }

        from polyzymd.workflow.slurm_submit import run_sbatch

        result = run_sbatch(script_path, module_load=module_load)
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

        1. ``solvated_system.pdb`` (preferred system topology with chain IDs)
        2. ``<prefix>.pdb`` (from system name)
        3. ``<prefix>.gro``
        4. Any ``*.gro`` (sorted, first match)

        Trajectory search order:

        1. ``prod_centered.xtc`` (whole molecules, centered protein)
        2. ``prod_nojump.xtc`` (nojump only)
        3. ``prod.xtc`` (raw production output)
        4. Any ``*.xtc`` (sorted)

        Parameters
        ----------
        working_dir : Path
            Replicate working directory.
        replicate : int
            Replicate index (unused, kept for interface parity).

        Returns
        -------
        TrajectoryLayout
            Resolved XTC layout with preferred post-processed trajectories
            and PDB topology preferred over GRO.
        """
        _ = replicate
        logger = logging.getLogger(__name__)

        trajectory_paths: list[Path] = []
        preferred_trajectories = ["prod_centered.xtc", "prod_nojump.xtc", "prod.xtc"]
        for filename in preferred_trajectories:
            candidate = working_dir / filename
            if candidate.exists():
                trajectory_paths = [candidate]
                logger.info("Resolved trajectory for analysis: %s", candidate.name)
                break

        if not trajectory_paths:
            trajectory_paths = sorted(working_dir.glob("*.xtc"))
            if trajectory_paths:
                logger.info("Resolved trajectory for analysis: %s", trajectory_paths[0].name)

        topology_path: Path | None = None
        topology_format = "pdb"

        prefix = self._generate_prefix()
        solvated_system_pdb = working_dir / "solvated_system.pdb"
        if solvated_system_pdb.exists():
            topology_path = solvated_system_pdb

        if topology_path is None and prefix:
            named_pdb = working_dir / f"{prefix}.pdb"
            if named_pdb.exists():
                topology_path = named_pdb

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

        if topology_path is not None:
            logger.info("Resolved topology for analysis: %s", topology_path.name)

        return TrajectoryLayout(
            topology_path=topology_path,
            trajectory_paths=trajectory_paths,
            trajectory_format="xtc",
            topology_format=topology_format,
        )
