"""
Job submission for HPC SLURM scheduler.

This module provides utilities for submitting daisy-chain MD simulation
jobs to SLURM. In PolyzyMD, daisy-chain is the canonical term for serial
MD segments on preempted hardware: each replicate gets a job script that
calls ``polyzymd run-segment``, checks progress, and resubmits itself
until the simulation is complete.

.. versionchanged:: 1.1.0
    Standardized daisy-chain execution on self-resubmitting jobs, where
    each submission advances one serial MD segment before scheduling the
    next segment as needed.
"""

from __future__ import annotations

import getpass
import logging
import os
import re
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
from polyzymd.workflow.slurm_submit import make_log_folder

LOGGER = logging.getLogger(__name__)

_SLURM_JOB_NAME_UNSAFE = re.compile(r"[^A-Za-z0-9._-]+")
_SLURM_JOB_NAME_UNDERSCORES = re.compile(r"_+")


# ---------------------------------------------------------------------------
# squeue-based duplicate detection
# ---------------------------------------------------------------------------


def check_existing_slurm_jobs(
    run_dir: Union[str, Path], job_name: Optional[str] = None
) -> List[str]:
    """Query SLURM for RUNNING or PENDING jobs of the run in *run_dir*.

    A job matches as described in :func:`job_belongs_to_run`: by its working
    directory (``squeue`` field ``%Z``), or, for a chain submitted by an
    older PolyzyMD, by its name *job_name*.

    This is a best-effort check: if ``squeue`` is unavailable (e.g. in a
    non-SLURM environment or CI), a warning is logged and an empty list is
    returned so that submission proceeds unimpeded.

    Parameters
    ----------
    run_dir : str or Path
        The replicate's run directory.
    job_name : str, optional
        The replicate's job name, from :func:`create_job_name`.

    Returns
    -------
    list of str
        IDs of the user's RUNNING or PENDING jobs of the run.
        Empty if ``squeue`` is unavailable or returns no matches.
    """
    user = os.environ.get("USER") or getpass.getuser()
    try:
        result = subprocess.run(
            [
                "squeue",
                "--noheader",
                "--user",
                user,
                "--states",
                "RUNNING,PENDING",
                "--format",
                "%i|%j|%Z",
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

    job_ids = []
    for line in result.stdout.splitlines():
        job_id, _, rest = line.strip().partition("|")
        name, _, work_dir = rest.partition("|")
        if job_id and job_belongs_to_run(name, work_dir, run_dir, job_name):
            job_ids.append(job_id)
    return job_ids


def job_belongs_to_run(
    name: str,
    work_dir: Union[str, Path, None],
    run_dir: Union[str, Path],
    job_name: Optional[str],
) -> bool:
    """Return True when the job *name* working in *work_dir* runs *run_dir*.

    New jobs work in the run directory (GROMACS jobs in ``<run_dir>/gromacs``).
    Chains submitted by older PolyzyMD versions work in the folder they were
    submitted from, so a job named *job_name* also matches, unless it works in
    a folder with the run directory's name: that is the same-named run of
    another condition.
    """
    if is_within_run_dir(work_dir, run_dir):
        return True
    if not job_name or name != job_name or not work_dir:
        return False
    folder = Path(work_dir)
    if folder.name == "gromacs":
        folder = folder.parent
    return folder.name != Path(run_dir).name


def is_within_run_dir(path: Union[str, Path, None], run_dir: Union[str, Path]) -> bool:
    """Return True when *path* is *run_dir* or a folder inside it."""
    if not path:
        return False
    target = Path(run_dir).resolve()
    candidate = Path(path).resolve()
    return candidate == target or target in candidate.parents


def cancel_slurm_jobs(job_ids: List[str]) -> List[str]:
    """Cancel SLURM jobs by ID, best effort.

    Mirrors :func:`check_existing_slurm_jobs`: outside a SLURM environment
    (CI, a workstation) the absence of ``scancel`` is logged and reported as
    "nothing cancelled" rather than raised, so that the STOP-file half of
    ``polyzymd cancel`` still works.

    Parameters
    ----------
    job_ids : list of str
        SLURM job IDs to cancel.

    Returns
    -------
    list of str
        The job IDs that ``scancel`` accepted.
    """
    if not job_ids:
        return []
    try:
        result = subprocess.run(
            ["scancel", *job_ids],
            capture_output=True,
            text=True,
            timeout=30,
        )
    except FileNotFoundError:
        LOGGER.warning(
            "scancel not found — no SLURM jobs were cancelled "
            "(this is expected outside of SLURM environments)"
        )
        return []
    except subprocess.TimeoutExpired:
        LOGGER.warning("scancel timed out — jobs may still be queued")
        return []
    except OSError as exc:
        LOGGER.warning(f"scancel failed ({exc}) — jobs may still be queued")
        return []

    if result.returncode != 0:
        LOGGER.warning(f"scancel returned exit code {result.returncode}: {result.stderr.strip()}")
        return []
    return list(job_ids)


def _sanitize_slurm_job_name(job_name: object, replicate: int) -> str:
    """Sanitize a job name for SLURM and log-file use.

    Parameters
    ----------
    job_name : object
        Candidate job name.
    replicate : int
        Replicate number used for the fallback name.

    Returns
    -------
    str
        Sanitized job name.
    """
    sanitized = _SLURM_JOB_NAME_UNSAFE.sub("_", str(job_name).strip())
    sanitized = _SLURM_JOB_NAME_UNDERSCORES.sub("_", sanitized).strip("_")
    if not sanitized:
        sanitized = f"pzmd_r{replicate}"
    return sanitized


def create_job_name(sim_config: SimulationConfig, replicate: int) -> str:
    """Create a sanitized SLURM job name for a replicate.

    The name is :meth:`SimulationConfig.format_run_directory_name`, so SLURM
    job names match run directory names.

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
    return _sanitize_slurm_job_name(sim_config.format_run_directory_name(replicate), replicate)


@dataclass
class DaisyChainConfig:
    """Configuration for daisy-chain job submission.

    Daisy-chain submission is PolyzyMD's canonical SLURM workflow for
    serial MD segments on preempted hardware. Each replicate is managed by
    one self-resubmitting job script that advances the trajectory segment
    by segment until production is complete.

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
        Initial segment index for the self-resubmitting daisy-chain job.
    replicate : int
        Replicate number.
    is_dry_run : bool
        Whether this was a dry run.
    is_generated_only : bool
        Whether this was a generate-only script output.
    """

    job_id: str
    script_path: Path
    segment_index: int
    replicate: int
    is_dry_run: bool = False
    is_generated_only: bool = False


class DaisyChainSubmitter:
    """Handle daisy-chain job submission for MD simulations.

    In PolyzyMD's daisy-chain model, each replicate gets a single
    self-resubmitting job script. The script calls ``polyzymd run-segment``,
    checks progress, and resubmits itself to run serial MD segments until
    the simulation is complete.

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
        pixi_env: str = "sim-cuda-12-4",
        openff_logs: bool = False,
    ) -> None:
        """Initialize the submitter.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        dc_config : DaisyChainConfig
            Submission configuration.
        pixi_env : str
            Pixi environment name (e.g. ``"sim-cuda-12-4"``, ``"sim-cuda-12-6"``).
        openff_logs : bool
            Enable verbose OpenFF logs in generated scripts.
        """
        self._sim_config = sim_config
        self._dc_config = dc_config
        self._openff_logs = openff_logs
        # The simulation environments have no OpenFF, so jobs load the build.
        self._generator = SlurmScriptGenerator(
            dc_config.slurm_config, pixi_env, openff_logs=openff_logs, skip_build=True
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

    def _generate_job_script(self, replicate: int, job_name: str) -> str:
        """Generate a self-resubmitting job script with a precomputed job name.

        Parameters
        ----------
        replicate : int
            Replicate number.
        job_name : str
            Sanitized SLURM job name shared with duplicate detection and log paths.

        Returns
        -------
        str
            Complete SLURM batch script content.
        """
        logs_subdir = self._sim_config.output.slurm_logs_subdir
        output_file = f"{logs_subdir}/{job_name}.%j.out"

        return self._generator.generate_job_script(
            config_path=self._dc_config.config_path,
            replicate=replicate,
            working_dir=self._get_scratch_dir(replicate),
            job_name=job_name,
            output_file=output_file,
        )

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
        return self._generate_job_script(replicate, self._create_job_name(replicate))

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
                is_dry_run=False,
                is_generated_only=True,
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
                is_generated_only=False,
            )

        # The job starts in its run directory (``#SBATCH --chdir``), and sbatch
        # does not create the folder of its log.
        Path(self._get_scratch_dir(replicate)).mkdir(parents=True, exist_ok=True)
        make_log_folder(script_path)

        # Use --export=NONE to start with clean environment, letting the
        # script's pixi shell-hook initialization work properly regardless
        # of submission context
        cmd = ["sbatch", "--export=NONE"]

        if self._dc_config.slurm_config.exclude:
            cmd.extend(["--exclude", self._dc_config.slurm_config.exclude])

        cmd.append(str(script_path))

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            match = re.search(r"\b(\d+)\b", result.stdout)
            if not match:
                raise RuntimeError(f"Could not parse job ID from sbatch output: {result.stdout!r}")
            job_id = match.group(1)
            LOGGER.info(f"Submitted job {job_id} for replicate {replicate}")

            return SubmissionResult(
                job_id=job_id,
                script_path=script_path,
                segment_index=0,
                replicate=replicate,
                is_dry_run=False,
                is_generated_only=False,
            )

        except subprocess.CalledProcessError as e:
            LOGGER.error(f"Error submitting job: {e}")
            LOGGER.error(f"STDOUT: {e.stdout}")
            LOGGER.error(f"STDERR: {e.stderr}")
            raise RuntimeError(f"Failed to submit job: {e.stderr}") from e

    def submit_replicate(self, replicate: int) -> SubmissionResult:
        """Generate and submit the job for a single replicate.

        Before submitting, checks ``squeue`` for existing RUNNING/PENDING
        jobs of this replicate (see :func:`job_belongs_to_run`).  If duplicates are
        found and ``force`` is not set, raises ``RuntimeError``.

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
            run_dir = self._get_scratch_dir(replicate)
            existing = check_existing_slurm_jobs(run_dir, job_name)
            if existing:
                ids = ", ".join(existing)
                raise RuntimeError(
                    f"Replicate {replicate} already has RUNNING/PENDING SLURM "
                    f"job(s): {ids} (run directory {run_dir}). "
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
                is_generated_only=False,
            )
            self._job_chains[replicate] = [result]
            LOGGER.info(f"[DRY RUN] Would generate script: {script_path}")
            LOGGER.info(f"[DRY RUN] Would submit replicate {replicate} with sbatch")
            return result

        script_content = self._generate_job_script(replicate, job_name)
        filename = f"run_rep{replicate}.sh"
        script_path = self._save_script(script_content, filename)

        result = self._submit_job(script_path=script_path, replicate=replicate)

        # One active daisy-chain job owns serial segment progression
        self._job_chains[replicate] = [result]
        return result

    def submit_all(self) -> Dict[int, List[SubmissionResult]]:
        """Submit jobs for all replicates.

        Returns
        -------
        dict
            Mapping of replicate numbers to daisy-chain submission results.
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
    pixi_env: str = "sim-cuda-12-4",
    output_dir: str | Path | None = None,
    scratch_dir: str | Path | None = None,
    projects_dir: str | Path | None = None,
    time_limit: str | None = None,
    memory: str | None = None,
    account: str | None = None,
    partition: str | None = None,
    qos: str | None = None,
    gpu_type: str | None = None,
    constraint: str | None = None,
    nodelist: str | None = None,
    exclude: str | None = None,
    openff_logs: bool = False,
) -> Dict[int, List[SubmissionResult]]:
    """Submit daisy-chain simulation jobs from a YAML config.

    This is the main entry point called by ``polyzymd submit``. Daisy-chain
    is PolyzyMD's canonical term for serial MD segments on preempted hardware;
    this function submits one self-resubmitting job per replicate to advance
    those segments until completion.

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
        Pixi environment name (e.g. ``"sim-cuda-12-4"``, ``"sim-cuda-12-6"``).
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
    partition : str or None
        Override SLURM partition.
    qos : str or None
        Override SLURM QoS value.
    gpu_type : str or None
        Override GPU type for presets that use ``--gpus`` directive.
    constraint : str or None
        SLURM ``--constraint`` expression (e.g. ``"A40|A100"``).
    nodelist : str or None
        Optional SLURM ``--nodelist`` override.
    exclude : str or None
        Optional SLURM ``--exclude`` override.  Replaces (never appends to)
        the preset's excluded-node list.
    openff_logs : bool
        Enable verbose OpenFF logs in generated scripts.

    Returns
    -------
    dict
        Mapping of replicate numbers to submission results.

    Raises
    ------
    ValueError
        If the SLURM account is empty on a preset that requires one
        and neither ``dry_run`` nor ``generate_only`` is set.
    FileNotFoundError
        If a replicate has no valid build from ``polyzymd build`` (not
        checked for ``dry_run``).
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
    if partition:
        slurm_config.partition = partition
    if qos:
        slurm_config.qos = qos
    if gpu_type:
        slurm_config.gpu_type = gpu_type
    if constraint:
        slurm_config.constraint = constraint
    if nodelist is not None:
        slurm_config.nodelist = nodelist
    if exclude is not None:
        slurm_config.exclude = exclude or None

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

    # Jobs run in a simulation environment without OpenFF, so they cannot
    # build. Check every replicate before any job is written or submitted.
    if not dry_run:
        from polyzymd.simulation.artifact_integrity import validate_build_bundle

        for replicate in dc_config.replicates:
            working_dir = sim_config.get_working_directory(replicate)
            try:
                validate_build_bundle(working_dir, sim_config)
            except (OSError, RuntimeError, ValueError) as exc:
                raise FileNotFoundError(
                    f"Replicate {replicate} has no usable build in {working_dir}: {exc}\n"
                    f"Build it first, in a compute job: polyzymd build -c {config_path} "
                    f"-r {replicate}. Then submit again."
                ) from exc

    # Create submitter and submit
    submitter = DaisyChainSubmitter(
        sim_config, dc_config, pixi_env=pixi_env, openff_logs=openff_logs
    )
    return submitter.submit_all()
