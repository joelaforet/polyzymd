"""
Job submission for HPC SLURM scheduler.

This module provides utilities for submitting self-resubmitting MD
simulation jobs to SLURM.  Each replicate gets a single job script
that calls ``polyzymd run-segment``, checks progress, and resubmits
itself until the simulation is complete.

.. versionchanged:: 1.1.0
    Replaced the legacy daisy-chain (dependency-chain) model with
    self-resubmitting jobs.  The public API (``submit_daisy_chain``,
    ``DaisyChainConfig``, ``DaisyChainSubmitter``) is preserved for
    backward compatibility but the internal behaviour is simplified.
"""

from __future__ import annotations

import logging
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

from polyzymd.config.schema import SimulationConfig
from polyzymd.utils.replicates import parse_replicate_range, validate_replicate_range
from polyzymd.workflow.slurm import (
    JobContext,
    SlurmConfig,
    SlurmScriptGenerator,
)

LOGGER = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# squeue-based duplicate detection
# ---------------------------------------------------------------------------


def check_existing_slurm_jobs(job_name: str) -> List[str]:
    """Query SLURM for RUNNING or PENDING jobs that match *job_name*.

    This is a best-effort check: if ``squeue`` is unavailable (e.g. in a
    non-SLURM environment or CI), a warning is logged and an empty list is
    returned so that submission proceeds unimpeded.

    Parameters
    ----------
    job_name : str
        The SLURM ``--job-name`` to search for (exact match).

    Returns
    -------
    list of str
        SLURM job IDs that are RUNNING or PENDING with the given name.
        Empty if ``squeue`` is unavailable or returns no matches.
    """
    try:
        result = subprocess.run(
            [
                "squeue",
                "--noheader",
                "--name",
                job_name,
                "--states",
                "RUNNING,PENDING",
                "--format",
                "%i",
            ],
            capture_output=True,
            text=True,
            timeout=15,
        )
    except FileNotFoundError:
        LOGGER.warning(
            "squeue not found — skipping duplicate-job check "
            "(this is expected outside of SLURM environments)"
        )
        return []
    except subprocess.TimeoutExpired:
        LOGGER.warning("squeue timed out — skipping duplicate-job check")
        return []
    except OSError as exc:
        LOGGER.warning(f"squeue failed ({exc}) — skipping duplicate-job check")
        return []

    if result.returncode != 0:
        LOGGER.warning(
            f"squeue returned exit code {result.returncode} — skipping duplicate-job check"
        )
        return []

    job_ids = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    return job_ids


def create_job_name(sim_config: SimulationConfig, replicate: int) -> str:
    """Create a descriptive SLURM job name for a replicate.

    Produces names like ``r1_310K_Fibronectin_SBMA-OEGMA_A75_B25``
    matching the directory naming convention.

    Parameters
    ----------
    sim_config : SimulationConfig
        Validated simulation configuration.
    replicate : int
        Replicate number.

    Returns
    -------
    str
        Formatted job name.
    """
    enzyme = sim_config.enzyme.name
    temp = int(sim_config.thermodynamics.temperature)

    polymer_info = ""
    if sim_config.polymers and sim_config.polymers.enabled:
        prefix = sim_config.polymers.type_prefix
        probs = {m.label: m.probability for m in sim_config.polymers.monomers}
        composition = "_".join(f"{lbl}{int(probs[lbl] * 100)}" for lbl in sorted(probs))
        polymer_info = f"_{prefix}_{composition}"

    return f"r{replicate}_{temp}K_{enzyme}{polymer_info}"


@dataclass
class DaisyChainConfig:
    """Configuration for job submission.

    Despite the legacy name, this now configures single self-resubmitting
    jobs (one per replicate) rather than dependency chains.

    Attributes
    ----------
    slurm_config : SlurmConfig
        SLURM job configuration.
    total_production_time_ns : float
        Total production time in nanoseconds.
    total_samples : int
        Total trajectory frames across the entire production run.
    equilibration_time_ns : float
        Equilibration time (informational only).
    replicates : list of int
        Replicate numbers to run.
    dry_run : bool
        If True, preview only. No scripts are written and no jobs are submitted.
    generate_only : bool
        If True, create scripts but don't submit.
    force : bool
        If True, skip the squeue duplicate-job check and submit even if
        a RUNNING/PENDING job already exists for the same replicate.
    output_script_dir : Path
        Directory for generated job scripts.
    config_path : str
        Path to the YAML configuration file.
    """

    slurm_config: SlurmConfig
    total_production_time_ns: float
    total_samples: int = 2500
    equilibration_time_ns: float = 0.5
    replicates: List[int] = field(default_factory=lambda: [1])
    dry_run: bool = False
    generate_only: bool = False
    force: bool = False
    output_script_dir: Path = Path("daisy_chain_scripts")
    config_path: str = "config.yaml"

    @classmethod
    def from_simulation_config(
        cls,
        sim_config: SimulationConfig,
        slurm_config: SlurmConfig,
        replicates: Union[str, List[int]] = "1",
        dry_run: bool = False,
        generate_only: bool = False,
        force: bool = False,
        output_script_dir: Union[str, Path] = "daisy_chain_scripts",
        config_path: str = "config.yaml",
    ) -> "DaisyChainConfig":
        """Create DaisyChainConfig from a SimulationConfig.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        slurm_config : SlurmConfig
            SLURM configuration.
        replicates : str or list of int
            Replicate range string (e.g. ``"1-5"``) or list of ints.
        dry_run : bool
            If True, preview only and write no files.
        generate_only : bool
            If True, create scripts without submitting.
        force : bool
            If True, skip duplicate-job check.
        output_script_dir : str or Path
            Directory for job scripts.
        config_path : str
            Path to the YAML configuration file.

        Returns
        -------
        DaisyChainConfig
            Configured instance.
        """
        # Parse replicates if string
        if isinstance(replicates, str):
            validate_replicate_range(replicates)
            replicate_list = parse_replicate_range(replicates)
        else:
            replicate_list = replicates

        return cls(
            slurm_config=slurm_config,
            total_production_time_ns=sim_config.simulation_phases.production.duration,
            total_samples=sim_config.simulation_phases.production.samples,
            equilibration_time_ns=sim_config.simulation_phases.total_equilibration_duration,
            replicates=replicate_list,
            dry_run=dry_run,
            generate_only=generate_only,
            force=force,
            output_script_dir=Path(output_script_dir),
            config_path=config_path,
        )


@dataclass
class SubmissionResult:
    """Result of job submission.

    Attributes
    ----------
    job_id : str
        SLURM job ID (or dummy ID for dry run).
    script_path : Path
        Path to the generated script.
    segment_index : int
        Always 0 in the self-resubmitting model (kept for compatibility).
    replicate : int
        Replicate number.
    is_dry_run : bool
        Whether this was a dry run.
    """

    job_id: str
    script_path: Path
    segment_index: int
    replicate: int
    is_dry_run: bool = False


class DaisyChainSubmitter:
    """Handles job submission for MD simulations.

    In the self-resubmitting model, each replicate gets a single job
    script.  The script calls ``polyzymd run-segment``, checks progress,
    and resubmits itself until the simulation is complete.

    Example
    -------
    >>> sim_config = SimulationConfig.from_yaml("config.yaml")
    >>> slurm_config = SlurmConfig.from_preset("aa100", email="user@example.com")
    >>> dc_config = DaisyChainConfig.from_simulation_config(
    ...     sim_config, slurm_config, replicates="1-3"
    ... )
    >>> submitter = DaisyChainSubmitter(sim_config, dc_config)
    >>> results = submitter.submit_all()
    """

    def __init__(
        self,
        sim_config: SimulationConfig,
        dc_config: DaisyChainConfig,
        pixi_env: str = "cuda-12-4",
        openff_logs: bool = False,
        skip_build: bool = False,
    ) -> None:
        """Initialize the submitter.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        dc_config : DaisyChainConfig
            Submission configuration.
        pixi_env : str
            Pixi environment name (e.g. ``"cuda-12-4"``, ``"cuda-12-6"``).
        openff_logs : bool
            Enable verbose OpenFF logs in generated scripts.
        skip_build : bool
            Skip system building in generated scripts.
        """
        self._sim_config = sim_config
        self._dc_config = dc_config
        self._openff_logs = openff_logs
        self._skip_build = skip_build
        self._generator = SlurmScriptGenerator(
            dc_config.slurm_config, pixi_env, openff_logs=openff_logs, skip_build=skip_build
        )

        # Track submitted jobs per replicate
        self._job_chains: Dict[int, List[SubmissionResult]] = {}

    @property
    def sim_config(self) -> SimulationConfig:
        """Get the simulation configuration."""
        return self._sim_config

    @property
    def dc_config(self) -> DaisyChainConfig:
        """Get the submission configuration."""
        return self._dc_config

    @property
    def job_chains(self) -> Dict[int, List[SubmissionResult]]:
        """Get the submission results for all replicates."""
        return self._job_chains

    def _create_job_name(self, replicate: int) -> str:
        """Create a descriptive job name for a replicate.

        Delegates to the module-level :func:`create_job_name` function.

        Parameters
        ----------
        replicate : int
            Replicate number.

        Returns
        -------
        str
            Formatted job name.
        """
        return create_job_name(self._sim_config, replicate)

    def _get_scratch_dir(self, replicate: int) -> str:
        """Get the scratch directory path for a replicate.

        Parameters
        ----------
        replicate : int
            Replicate number.

        Returns
        -------
        str
            Absolute scratch directory path.
        """
        scratch_dir = self._sim_config.get_working_directory(replicate)
        return str(scratch_dir.resolve())

    def generate_job_script(self, replicate: int) -> str:
        """Generate a self-resubmitting job script for a replicate.

        Parameters
        ----------
        replicate : int
            Replicate number.

        Returns
        -------
        str
            Complete SLURM batch script content.
        """
        job_name = self._create_job_name(replicate)
        logs_subdir = self._sim_config.output.slurm_logs_subdir
        output_file = f"{logs_subdir}/{job_name}.%j.out"

        return self._generator.generate_job_script(
            config_path=self._dc_config.config_path,
            replicate=replicate,
            working_dir=self._get_scratch_dir(replicate),
            job_name=job_name,
            output_file=output_file,
        )

    def _save_script(self, content: str, filename: str) -> Path:
        """Save a script to the output directory.

        Parameters
        ----------
        content : str
            Script content.
        filename : str
            Script filename.

        Returns
        -------
        Path
            Path to saved script.
        """
        output_dir = self._dc_config.output_script_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        script_path = output_dir / filename
        with open(script_path, "w") as f:
            f.write(content)

        os.chmod(script_path, 0o755)
        return script_path

    def _submit_job(
        self,
        script_path: Path,
        replicate: int,
    ) -> SubmissionResult:
        """Submit a job to SLURM.

        Parameters
        ----------
        script_path : Path
            Path to the job script.
        replicate : int
            Replicate number.

        Returns
        -------
        SubmissionResult
            Submission result with job information.
        """
        if self._dc_config.generate_only:
            job_id = f"GENERATED_{replicate}"
            LOGGER.info(f"[GENERATE ONLY] Script generated: {script_path}")
            return SubmissionResult(
                job_id=job_id,
                script_path=script_path,
                segment_index=0,
                replicate=replicate,
                is_dry_run=True,
            )

        if self._dc_config.dry_run:
            job_id = f"DRY_RUN_{replicate}"
            LOGGER.info(f"[DRY RUN] Would generate and submit {script_path}")
            return SubmissionResult(
                job_id=job_id,
                script_path=script_path,
                segment_index=0,
                replicate=replicate,
                is_dry_run=True,
            )

        # Use --export=NONE to start with clean environment, letting the
        # script's pixi shell-hook initialization work properly regardless
        # of submission context
        cmd = ["sbatch", "--export=NONE"]

        if self._dc_config.slurm_config.exclude:
            cmd.extend(["--exclude", self._dc_config.slurm_config.exclude])

        cmd.append(str(script_path))

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            job_id = result.stdout.strip().split()[-1]
            LOGGER.info(f"Submitted job {job_id} for replicate {replicate}")

            return SubmissionResult(
                job_id=job_id,
                script_path=script_path,
                segment_index=0,
                replicate=replicate,
                is_dry_run=False,
            )

        except subprocess.CalledProcessError as e:
            LOGGER.error(f"Error submitting job: {e}")
            LOGGER.error(f"STDOUT: {e.stdout}")
            LOGGER.error(f"STDERR: {e.stderr}")
            raise RuntimeError(f"Failed to submit job: {e.stderr}") from e

    def submit_replicate(self, replicate: int) -> SubmissionResult:
        """Generate and submit the job for a single replicate.

        Before submitting, checks ``squeue`` for existing RUNNING/PENDING
        jobs with the same job name.  If duplicates are found and
        ``force`` is not set, raises ``RuntimeError``.

        Parameters
        ----------
        replicate : int
            Replicate number.

        Returns
        -------
        SubmissionResult
            Submission result.

        Raises
        ------
        RuntimeError
            If a SLURM job is already RUNNING or PENDING for this
            replicate and ``force`` is False.
        """
        job_name = self._create_job_name(replicate)

        # Best-effort duplicate guard (only when actually submitting)
        if (
            not self._dc_config.dry_run
            and not self._dc_config.generate_only
            and not self._dc_config.force
        ):
            existing = check_existing_slurm_jobs(job_name)
            if existing:
                ids = ", ".join(existing)
                raise RuntimeError(
                    f"Replicate {replicate} already has RUNNING/PENDING SLURM "
                    f"job(s): {ids} (job name '{job_name}'). "
                    "Use --force to submit anyway."
                )

        LOGGER.info(f"Submitting self-resubmitting job for replicate {replicate}")

        if self._dc_config.dry_run:
            script_path = self._dc_config.output_script_dir / f"run_rep{replicate}.sh"
            result = SubmissionResult(
                job_id=f"DRY_RUN_{replicate}",
                script_path=script_path,
                segment_index=0,
                replicate=replicate,
                is_dry_run=True,
            )
            self._job_chains[replicate] = [result]
            LOGGER.info(f"[DRY RUN] Would generate script: {script_path}")
            LOGGER.info(f"[DRY RUN] Would submit replicate {replicate} with sbatch")
            return result

        script_content = self.generate_job_script(replicate)
        filename = f"run_rep{replicate}.sh"
        script_path = self._save_script(script_content, filename)

        result = self._submit_job(script_path=script_path, replicate=replicate)

        # Store as a single-element list for backward compatibility
        self._job_chains[replicate] = [result]
        return result

    def submit_all(self) -> Dict[int, List[SubmissionResult]]:
        """Submit jobs for all replicates.

        Returns
        -------
        dict
            Mapping of replicate numbers to their submission results
            (each value is a single-element list for compatibility).
        """
        self._print_submission_summary()

        for replicate in self._dc_config.replicates:
            self.submit_replicate(replicate)

        self._print_completion_summary()
        return self._job_chains

    def _print_submission_summary(self) -> None:
        """Print a summary before submission."""
        config = self._dc_config
        num_replicates = len(config.replicates)

        LOGGER.info("\nPreparing self-resubmitting simulation jobs")
        LOGGER.info(f"  Enzyme: {self._sim_config.enzyme.name}")

        if self._sim_config.polymers and self._sim_config.polymers.enabled:
            LOGGER.info(f"  Polymer: {self._sim_config.polymers.type_prefix}")
            LOGGER.info(f"  Polymer count: {self._sim_config.polymers.count}")

        LOGGER.info(f"  Temperature: {self._sim_config.thermodynamics.temperature} K")
        LOGGER.info(f"  Total production time: {config.total_production_time_ns} ns")
        LOGGER.info(f"  Total samples: {config.total_samples}")
        LOGGER.info(f"  Replicates: {config.replicates} ({num_replicates} total)")
        LOGGER.info(f"  Jobs to submit: {num_replicates} (one per replicate, self-resubmitting)")
        LOGGER.info("")
        LOGGER.info("SLURM Configuration:")
        LOGGER.info(f"  Partition: {config.slurm_config.partition}")
        if config.slurm_config.qos:
            LOGGER.info(f"  QoS: {config.slurm_config.qos}")
        if config.slurm_config.account:
            LOGGER.info(f"  Account: {config.slurm_config.account}")
        LOGGER.info(f"  Time limit: {config.slurm_config.time_limit}")
        LOGGER.info("")

        if config.dry_run:
            LOGGER.info("*** DRY RUN MODE - Preview only, no files will be written ***")
            LOGGER.info("")
        elif config.generate_only:
            LOGGER.info("*** GENERATE-ONLY MODE - Scripts will be created but not submitted ***")
            LOGGER.info("")

    def _print_completion_summary(self) -> None:
        """Print a summary after submission."""
        config = self._dc_config
        total_jobs = sum(len(chain) for chain in self._job_chains.values())

        if config.dry_run:
            LOGGER.info(f"\nDry run completed. {total_jobs} replicate(s) previewed.")
            LOGGER.info("No files were written and no jobs were submitted.")
        elif config.generate_only:
            LOGGER.info(f"\nScript generation complete. {total_jobs} job script(s) created.")
            LOGGER.info(f"Scripts saved to: {config.output_script_dir}")
            LOGGER.info("Review the scripts and run without --generate-only to submit them.")
        else:
            LOGGER.info(f"\nAll {total_jobs} job(s) submitted successfully!")
            LOGGER.info("\nSubmitted jobs:")

            for replicate, results in sorted(self._job_chains.items()):
                job_id = results[0].job_id
                LOGGER.info(f"  Replicate {replicate}: job {job_id} (self-resubmitting)")

            LOGGER.info("\nEach job will automatically resubmit until the simulation completes.")
            LOGGER.info("Monitor progress with: squeue -u $USER")
            LOGGER.info(
                "Check simulation status with: polyzymd check-progress -c <config> -r <rep>"
            )


def submit_daisy_chain(
    config_path: Union[str, Path],
    slurm_preset: str = "aa100",
    replicates: str = "1",
    email: str = "",
    dry_run: bool = False,
    generate_only: bool = False,
    force: bool = False,
    pixi_env: str = "cuda-12-4",
    output_dir: Optional[Union[str, Path]] = None,
    scratch_dir: Optional[Union[str, Path]] = None,
    projects_dir: Optional[Union[str, Path]] = None,
    time_limit: Optional[str] = None,
    memory: Optional[str] = None,
    account: Optional[str] = None,
    gpu_type: Optional[str] = None,
    openff_logs: bool = False,
    skip_build: bool = False,
) -> Dict[int, List[SubmissionResult]]:
    """Submit self-resubmitting simulation jobs from a YAML config.

    This is the main entry point called by ``polyzymd submit``.  Despite
    the legacy function name, it now submits one self-resubmitting job
    per replicate rather than a chain of dependent jobs.

    Parameters
    ----------
    config_path : str or Path
        Path to simulation YAML config.
    slurm_preset : str
        SLURM preset name (aa100, al40, blanca-shirts, bridges2, testing).
    replicates : str
        Replicate range string (e.g. ``"1-5"``, ``"1,3,5"``).
    email : str
        Email for job notifications.
    dry_run : bool
        If True, preview only and write no files.
    generate_only : bool
        If True, create scripts without submitting.
    force : bool
        If True, skip the squeue duplicate-job check.
    pixi_env : str
        Pixi environment name (e.g. ``"cuda-12-4"``, ``"cuda-12-6"``).
    output_dir : str or Path or None
        Directory for job scripts.
    scratch_dir : str or Path or None
        Override scratch directory for simulation output.
    projects_dir : str or Path or None
        Override projects directory for scripts/logs.
    time_limit : str or None
        Override SLURM time limit (format: ``HH:MM:SS``).
    memory : str or None
        Override SLURM memory allocation (e.g. ``"4G"``).
    account : str or None
        Override SLURM account / allocation ID.
    gpu_type : str or None
        Override GPU type for presets that use ``--gpus`` directive.
    openff_logs : bool
        Enable verbose OpenFF logs in generated scripts.
    skip_build : bool
        Skip system building in generated scripts.

    Returns
    -------
    dict
        Mapping of replicate numbers to submission results.

    Raises
    ------
    ValueError
        If the SLURM account is empty on a preset that requires one
        and neither ``dry_run`` nor ``generate_only`` is set.
    """
    # Load simulation config
    sim_config = SimulationConfig.from_yaml(config_path)

    # Apply CLI overrides for directories
    if scratch_dir:
        sim_config.output.scratch_directory = Path(scratch_dir)
    if projects_dir:
        sim_config.output.projects_directory = Path(projects_dir)

    # Determine output script directory
    if output_dir:
        script_output_dir = Path(output_dir)
    else:
        script_output_dir = sim_config.output.get_job_scripts_directory()

    # Create SLURM config from preset
    slurm_config = SlurmConfig.from_preset(slurm_preset, email=email)  # type: ignore[arg-type]

    # Record whether the preset itself ships with an empty account.
    preset_account_is_empty = not slurm_config.account

    # Apply CLI overrides
    if time_limit:
        slurm_config.time_limit = time_limit
    if memory:
        slurm_config.memory = memory
    if account:
        slurm_config.account = account
    if gpu_type:
        slurm_config.gpu_type = gpu_type

    # Guard: an empty account on presets that require one (e.g. Alpine) will
    # produce an invalid SBATCH script.  Skip the guard when the preset itself
    # ships with account="" (e.g. bridges2).
    if not slurm_config.account and not preset_account_is_empty:
        msg = (
            f"SLURM account is required but was not set for preset '{slurm_preset}'. "
            "Pass your allocation ID with --account <id>."
        )
        if dry_run or generate_only:
            LOGGER.warning(msg)
        else:
            raise ValueError(msg)

    # Create submission config
    dc_config = DaisyChainConfig.from_simulation_config(
        sim_config=sim_config,
        slurm_config=slurm_config,
        replicates=replicates,
        dry_run=dry_run,
        generate_only=generate_only,
        force=force,
        output_script_dir=script_output_dir,
        config_path=str(Path(config_path).resolve()),
    )

    # Create submitter and submit
    submitter = DaisyChainSubmitter(
        sim_config, dc_config, pixi_env=pixi_env, openff_logs=openff_logs, skip_build=skip_build
    )
    return submitter.submit_all()
