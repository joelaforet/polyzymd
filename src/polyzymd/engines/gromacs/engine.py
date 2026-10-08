"""GROMACS engine implementation for PolyzyMD."""

from __future__ import annotations

import logging
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

        from polyzymd.analyses.shared.gromacs import system_prefix

        prefix = system_prefix(self._config)
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

        The script is written to ``extra["script_path"]`` when given, and to
        ``daisy_chain_scripts/run_rep<N>.sh`` otherwise.

        Parameters
        ----------
        request : EngineSubmitRequest
            Submission request details.

        Returns
        -------
        Path
            Path to the generated SLURM script.

        Raises
        ------
        FileNotFoundError
            If the replicate has no GROMACS inputs from ``polyzymd build``.
        """
        if request.slurm_config is None:
            raise ValueError("GROMACS submission requires slurm_config")

        pixi_env = str(request.extra.get("pixi_env", "build"))
        self.check_build(request)
        from polyzymd.analyses.shared.gromacs import system_prefix

        prefix = system_prefix(self._config)

        eq_mdps = sorted(path.name for path in request.working_dir.glob("eq_*.mdp"))
        if not eq_mdps:
            logging.getLogger(__name__).warning(
                "Core GROMACS inputs found in %s but no equilibration MDPs (eq_*.mdp). "
                "The generated script will skip equilibration and run production from EM output.",
                request.working_dir,
            )

        script_path = Path(
            request.extra.get(
                "script_path",
                request.working_dir / "daisy_chain_scripts" / f"run_rep{request.replicate}.sh",
            )
        )

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
            output_file=str(
                self._config.output.get_slurm_logs_directory() / f"{request.job_name}.%j.out"
            ),
        )
        generator.save_script(script, script_path)
        return script_path

    def check_build(self, request: EngineSubmitRequest) -> None:
        """Refuse a replicate without the GROMACS inputs from ``polyzymd build``.

        Submission never builds: building runs OpenFF and Packmol, which do not
        belong on a login node.

        Raises
        ------
        FileNotFoundError
            With the ``polyzymd build`` command to run.
        """
        from polyzymd.analyses.shared.gromacs import system_prefix

        prefix = system_prefix(self._config)
        required = [f"{prefix}.top", f"{prefix}.gro", "em.mdp", "prod.mdp"]
        missing = [name for name in required if not (request.working_dir / name).exists()]
        if missing:
            raise FileNotFoundError(
                f"No GROMACS build for replicate {request.replicate} in {request.working_dir} "
                f"(missing {', '.join(missing)}). Build it first, in a compute job: "
                f"polyzymd build -c {request.config_path} -r {request.replicate} "
                "--format gromacs. Then submit again."
            )

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
        ``-nb gpu``, ``-pme gpu`` and ``-bonded gpu`` are appended
        automatically. Each flag is skipped if the user already specified that
        flag key (e.g., ``-nb cpu`` prevents auto-adding ``-nb gpu``).
        ``-update gpu`` is never added: GROMACS updates on the GPU only with
        ``integrator = md``, and the Langevin thermostats export ``sd``.

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

        from polyzymd.workflow.slurm_submit import run_sbatch

        result = run_sbatch(script_path)
        if result.returncode != 0:
            raise RuntimeError(f"sbatch submission failed: {result.stderr.strip()}")

        return {
            "submitted": True,
            "script_path": script_path,
            "stdout": result.stdout.strip(),
            "stderr": result.stderr.strip(),
            "returncode": result.returncode,
        }

    #: Production trajectories a GROMACS run leaves: the raw ``mdrun`` output and the
    #: post-processed whole-molecule and centred copies that analyses read.
    TRAJECTORY_NAMES: ClassVar[tuple[str, ...]] = (
        "prod.xtc",
        "prod_nojump.xtc",
        "prod_centered.xtc",
    )

    def trajectory_files(
        self, working_dir: Path, progress: SimulationProgress | None
    ) -> list[Path]:
        """Return the production XTC files that exist, once the run has completed.

        ``mdrun`` appends every restart to ``prod.xtc``, so the run's files are
        finished only when it is: a run whose ``progress.json`` records it as
        anything but completed gives none. A run without one gives its files.

        Parameters
        ----------
        working_dir : Path
            GROMACS working directory of one replicate.
        progress : SimulationProgress or None
            The run's progress, or None when it has no ``progress.json``.

        Returns
        -------
        list of Path
            The existing files of :attr:`TRAJECTORY_NAMES`, in that order.
        """
        from polyzymd.simulation.progress import SimulationStatus

        if progress is not None and progress.status != SimulationStatus.COMPLETED:
            return []
        return [
            Path(working_dir) / name
            for name in self.TRAJECTORY_NAMES
            if (Path(working_dir) / name).is_file()
        ]

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

    def resolve_trajectory_layout(
        self,
        working_dir: Path,
        replicate: int,
        *,
        require_complete: bool = True,
    ) -> TrajectoryLayout:
        """Resolve GROMACS trajectory layout for downstream analyses.

        Topology search order:

        1. ``prod.tpr``, the compiled run input, which carries every atom,
           bond, mass and charge with no atom limit and is what analyses
           should read
        2. ``solvated_system.pdb`` (viewer topology with chain IDs; its
           CONECT records are unreadable above 99,999 atoms)
        3. ``<prefix>.pdb`` (from system name)
        4. ``<prefix>.gro``
        5. Any ``*.gro`` (sorted, first match); a GRO file carries no bonds

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
        require_complete : bool, optional
            Accepted for interface parity and ignored. The GROMACS layout is a
            single production XTC rather than a chain of segments, and
            ``progress.json`` records no per-file status for it, so this engine
            has nothing to exclude. A production XTC that is still being
            written is therefore still read; check the job state before
            analyzing a live GROMACS run.

        Returns
        -------
        TrajectoryLayout
            Resolved XTC layout with preferred post-processed trajectories
            and PDB topology preferred over GRO. ``segment_status`` and
            ``excluded_segments`` are always empty.
        """
        _ = replicate
        _ = require_complete
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

        from polyzymd.analyses.shared.gromacs import system_prefix

        prefix = system_prefix(self._config)
        production_tpr = working_dir / "prod.tpr"
        if production_tpr.exists():
            topology_path = production_tpr
            topology_format = "tpr"

        solvated_system_pdb = working_dir / "solvated_system.pdb"
        if topology_path is None and solvated_system_pdb.exists():
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

        from polyzymd.analyses.shared.gromacs import topology_name

        return TrajectoryLayout(
            topology_path=topology_path,
            gromacs_topology_path=working_dir / topology_name(self._config),
            trajectory_paths=trajectory_paths,
            trajectory_format="xtc",
            topology_format=topology_format,
        )
