"""
Daisy-chain job submission for HPC SLURM scheduler.

This module provides utilities for breaking long MD simulations into
smaller dependent jobs that are automatically chained together using
SLURM job dependencies.
"""

from __future__ import annotations

import logging
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

from polyzymd.config.schema import SimulationConfig
from polyzymd.workflow.slurm import (
    JobContext,
    SlurmConfig,
    SlurmScriptGenerator,
    parse_replicate_range,
    validate_replicate_range,
)

LOGGER = logging.getLogger(__name__)


@dataclass
class SegmentInfo:
    """Information about a single simulation segment.

    Attributes:
        index: Segment index (0-based for initial, 1+ for continuations)
        duration_ns: Duration of this segment in nanoseconds
        samples: Number of trajectory frames to save
        is_initial: Whether this is the initial (build + equilibration + first prod) segment
        cumulative_time_ns: Total simulated time up to and including this segment
    """

    index: int
    duration_ns: float
    samples: int
    is_initial: bool
    cumulative_time_ns: float


@dataclass
class DaisyChainConfig:
    """Configuration for daisy-chain submission.

    Attributes:
        slurm_config: SLURM job configuration
        total_production_time_ns: Total production time in nanoseconds
        total_segments: Number of segments to split production into
        total_samples: Total trajectory frames across all segments
        equilibration_time_ns: Equilibration time (only for initial segment)
        replicates: List of replicate numbers to run
        dry_run: If True, create scripts but don't submit
        output_script_dir: Directory for generated job scripts
        main_script: Path to the main simulation script (for initial job)
        main_script_args: Additional arguments for main script
    """

    slurm_config: SlurmConfig
    total_production_time_ns: float
    total_segments: int = 10
    total_samples: int = 2500
    equilibration_time_ns: float = 0.5
    replicates: List[int] = field(default_factory=lambda: [1])
    dry_run: bool = False
    output_script_dir: Path = Path("daisy_chain_scripts")
    main_script: Optional[str] = None
    main_script_args: Dict[str, Any] = field(default_factory=dict)

    @property
    def segment_duration_ns(self) -> float:
        """Get the duration of each segment in nanoseconds."""
        return self.total_production_time_ns / self.total_segments

    @property
    def samples_per_segment(self) -> int:
        """Get the number of frames per segment."""
        return self.total_samples // self.total_segments

    def get_segments(self) -> List[SegmentInfo]:
        """Generate segment information for all segments.

        Returns:
            List of SegmentInfo objects for each segment.
        """
        segments = []
        cumulative_time = 0.0

        for i in range(self.total_segments):
            duration = self.segment_duration_ns
            cumulative_time += duration

            segments.append(
                SegmentInfo(
                    index=i,
                    duration_ns=duration,
                    samples=self.samples_per_segment,
                    is_initial=(i == 0),
                    cumulative_time_ns=cumulative_time,
                )
            )

        return segments

    @classmethod
    def from_simulation_config(
        cls,
        sim_config: SimulationConfig,
        slurm_config: SlurmConfig,
        replicates: Union[str, List[int]] = "1",
        dry_run: bool = False,
        output_script_dir: Union[str, Path] = "daisy_chain_scripts",
        main_script: Optional[str] = None,
        main_script_args: Optional[Dict[str, Any]] = None,
    ) -> "DaisyChainConfig":
        """Create DaisyChainConfig from a SimulationConfig.

        Args:
            sim_config: Simulation configuration
            slurm_config: SLURM configuration
            replicates: Replicate range string (e.g., "1-5") or list of ints
            dry_run: If True, don't submit jobs
            output_script_dir: Directory for job scripts
            main_script: Main simulation script path
            main_script_args: Additional arguments for main script

        Returns:
            Configured DaisyChainConfig
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
            total_segments=sim_config.simulation_phases.segments,
            total_samples=sim_config.simulation_phases.production.samples,
            equilibration_time_ns=sim_config.simulation_phases.equilibration.duration,
            replicates=replicate_list,
            dry_run=dry_run,
            output_script_dir=Path(output_script_dir),
            main_script=main_script,
            main_script_args=main_script_args or {},
        )


@dataclass
class SubmissionResult:
    """Result of job submission.

    Attributes:
        job_id: SLURM job ID (or dummy ID for dry run)
        script_path: Path to the generated script
        segment_index: Segment index for this job
        replicate: Replicate number
        is_dry_run: Whether this was a dry run
    """

    job_id: str
    script_path: Path
    segment_index: int
    replicate: int
    is_dry_run: bool = False


class DaisyChainSubmitter:
    """Handles daisy-chain job submission for MD simulations.

    This class generates SLURM job scripts and submits them with proper
    dependencies so that continuation jobs run after their prerequisites.

    Example:
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
        conda_env: str = "polymerist-env",
    ) -> None:
        """Initialize the DaisyChainSubmitter.

        Args:
            sim_config: Simulation configuration
            dc_config: Daisy-chain configuration
            conda_env: Conda environment name
        """
        self._sim_config = sim_config
        self._dc_config = dc_config
        self._generator = SlurmScriptGenerator(dc_config.slurm_config, conda_env)

        # Track submitted jobs per replicate
        self._job_chains: Dict[int, List[SubmissionResult]] = {}

    @property
    def sim_config(self) -> SimulationConfig:
        """Get the simulation configuration."""
        return self._sim_config

    @property
    def dc_config(self) -> DaisyChainConfig:
        """Get the daisy-chain configuration."""
        return self._dc_config

    @property
    def job_chains(self) -> Dict[int, List[SubmissionResult]]:
        """Get the job chains for all replicates."""
        return self._job_chains

    def _create_job_name(self, segment_index: int, replicate: int) -> str:
        """Create a descriptive job name.

        Args:
            segment_index: Segment index
            replicate: Replicate number

        Returns:
            Formatted job name
        """
        enzyme = self._sim_config.enzyme.name
        temp = int(self._sim_config.thermodynamics.temperature)

        polymer_info = ""
        if self._sim_config.polymers and self._sim_config.polymers.enabled:
            prefix = self._sim_config.polymers.type_prefix
            # Get minority percentage
            probs = [m.probability for m in self._sim_config.polymers.monomers]
            minority_pct = int(min(probs) * 100)
            polymer_info = f"_{prefix}-{minority_pct}%"

        return f"s{segment_index}_r{replicate}_{temp}K_{enzyme}{polymer_info}"

    def _create_output_file_pattern(self, segment_index: int, replicate: int) -> str:
        """Create output file pattern for SLURM logs.

        SLURM logs go to the slurm_logs subdirectory within projects.

        Args:
            segment_index: Segment index
            replicate: Replicate number

        Returns:
            Output file pattern (relative to projects_dir)
        """
        job_name = self._create_job_name(segment_index, replicate)
        logs_subdir = self._sim_config.output.slurm_logs_subdir
        return f"{logs_subdir}/{job_name}.%A_%a.out"

    def _get_scratch_dir(self, replicate: int) -> str:
        """Get the scratch directory path for a replicate.

        This is where simulation output (trajectories, checkpoints) goes.

        Args:
            replicate: Replicate number

        Returns:
            Scratch directory path (absolute)
        """
        scratch_dir = self._sim_config.get_working_directory(replicate)
        return str(scratch_dir.resolve())

    def _get_projects_dir(self) -> str:
        """Get the projects directory path.

        This is where scripts, configs, and logs live.

        Returns:
            Projects directory path (absolute)
        """
        projects_dir = self._sim_config.get_projects_directory()
        return str(projects_dir.resolve())

    def _build_initial_script_args(self, replicate: int) -> Dict[str, Any]:
        """Build arguments for the initial simulation script.

        Args:
            replicate: Replicate number

        Returns:
            Dictionary of arguments
        """
        args = {
            "replicate": replicate,
            "config": str(self._dc_config.main_script_args.get("config", "config.yaml")),
            "segment_time": self._dc_config.segment_duration_ns,
            "segment_frames": self._dc_config.samples_per_segment,
            "num_segments": self._dc_config.total_segments,
        }

        # Add any extra args from config
        for key, value in self._dc_config.main_script_args.items():
            if key not in args:
                args[key] = value

        return args

    def generate_initial_script(self, replicate: int) -> str:
        """Generate the initial job script content.

        Args:
            replicate: Replicate number

        Returns:
            Script content string
        """
        context = JobContext(
            job_name=self._create_job_name(0, replicate),
            output_file=self._create_output_file_pattern(0, replicate),
            scratch_dir=self._get_scratch_dir(replicate),
            projects_dir=self._get_projects_dir(),
            segment_index=0,
            replicate_num=replicate,
        )

        python_script = self._dc_config.main_script or "polyzymd_run.py"
        python_args = self._build_initial_script_args(replicate)

        return self._generator.generate_initial_job(
            context=context,
            python_script=python_script,
            python_args=python_args,
        )

    def generate_continuation_script(self, segment_index: int, replicate: int) -> str:
        """Generate a continuation job script content.

        Args:
            segment_index: Segment index (1 or higher)
            replicate: Replicate number

        Returns:
            Script content string
        """
        context = JobContext(
            job_name=self._create_job_name(segment_index, replicate),
            output_file=self._create_output_file_pattern(segment_index, replicate),
            scratch_dir=self._get_scratch_dir(replicate),
            projects_dir=self._get_projects_dir(),
            segment_index=segment_index,
            replicate_num=replicate,
        )

        return self._generator.generate_continuation_job(
            context=context,
            segment_time=self._dc_config.segment_duration_ns,
            num_samples=self._dc_config.samples_per_segment,
        )

    def _save_script(self, content: str, filename: str) -> Path:
        """Save a script to the output directory.

        Args:
            content: Script content
            filename: Script filename

        Returns:
            Path to saved script
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
        segment_index: int,
        replicate: int,
        dependency_job_id: Optional[str] = None,
    ) -> SubmissionResult:
        """Submit a job to SLURM.

        Args:
            script_path: Path to the job script
            segment_index: Segment index
            replicate: Replicate number
            dependency_job_id: Job ID to depend on (for continuation jobs)

        Returns:
            SubmissionResult with job information
        """
        if self._dc_config.dry_run:
            job_id = f"DRY_RUN_{replicate}_{segment_index}"
            LOGGER.info(f"[DRY RUN] Would submit {script_path}")
            return SubmissionResult(
                job_id=job_id,
                script_path=script_path,
                segment_index=segment_index,
                replicate=replicate,
                is_dry_run=True,
            )

        # Build sbatch command
        cmd = ["sbatch"]

        if dependency_job_id:
            cmd.extend(["--dependency", f"afterok:{dependency_job_id}"])

        # Add exclude if configured
        if self._dc_config.slurm_config.exclude:
            cmd.extend(["--exclude", self._dc_config.slurm_config.exclude])

        cmd.append(str(script_path))

        # Submit
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            job_id = result.stdout.strip().split()[-1]
            LOGGER.info(f"Submitted job {job_id} from {script_path}")

            return SubmissionResult(
                job_id=job_id,
                script_path=script_path,
                segment_index=segment_index,
                replicate=replicate,
                is_dry_run=False,
            )

        except subprocess.CalledProcessError as e:
            LOGGER.error(f"Error submitting job: {e}")
            LOGGER.error(f"STDOUT: {e.stdout}")
            LOGGER.error(f"STDERR: {e.stderr}")
            raise RuntimeError(f"Failed to submit job: {e.stderr}") from e

    def submit_replicate_chain(self, replicate: int) -> List[SubmissionResult]:
        """Submit all jobs for a single replicate.

        Args:
            replicate: Replicate number

        Returns:
            List of SubmissionResults for all segments
        """
        LOGGER.info(f"Submitting job chain for replicate {replicate}")

        results: List[SubmissionResult] = []
        segments = self._dc_config.get_segments()

        for segment in segments:
            if segment.is_initial:
                # Initial job
                script_content = self.generate_initial_script(replicate)
                filename = f"initial_seg{segment.index}_rep{replicate}.sh"
                script_path = self._save_script(script_content, filename)

                result = self._submit_job(
                    script_path=script_path,
                    segment_index=segment.index,
                    replicate=replicate,
                    dependency_job_id=None,
                )

            else:
                # Continuation job
                script_content = self.generate_continuation_script(segment.index, replicate)
                filename = f"continue_seg{segment.index}_rep{replicate}.sh"
                script_path = self._save_script(script_content, filename)

                # Depend on previous segment
                prev_job_id = results[-1].job_id

                result = self._submit_job(
                    script_path=script_path,
                    segment_index=segment.index,
                    replicate=replicate,
                    dependency_job_id=prev_job_id,
                )

            results.append(result)

        self._job_chains[replicate] = results
        return results

    def submit_all(self) -> Dict[int, List[SubmissionResult]]:
        """Submit jobs for all replicates.

        Returns:
            Dictionary mapping replicate numbers to their job chains
        """
        self._print_submission_summary()

        for replicate in self._dc_config.replicates:
            self.submit_replicate_chain(replicate)

        self._print_completion_summary()
        return self._job_chains

    def _print_submission_summary(self) -> None:
        """Print a summary before submission."""
        config = self._dc_config
        num_replicates = len(config.replicates)
        total_jobs = num_replicates * config.total_segments

        print(f"\nPreparing {config.total_segments}-segment simulation jobs")
        print(f"  Enzyme: {self._sim_config.enzyme.name}")

        if self._sim_config.polymers and self._sim_config.polymers.enabled:
            print(f"  Polymer: {self._sim_config.polymers.type_prefix}")
            print(f"  Polymer count: {self._sim_config.polymers.count}")

        print(f"  Temperature: {self._sim_config.thermodynamics.temperature} K")
        print(f"  Total production time: {config.total_production_time_ns} ns")
        print(f"  Time per segment: {config.segment_duration_ns} ns")
        print(f"  Samples per segment: {config.samples_per_segment}")
        print(f"  Replicates: {config.replicates} ({num_replicates} total)")
        print(f"  Total jobs to submit: {total_jobs}")
        print(f"  Dependency chains: {num_replicates} independent chains")
        print()
        print(f"SLURM Configuration:")
        print(f"  Partition: {config.slurm_config.partition}")
        print(f"  QoS: {config.slurm_config.qos}")
        print(f"  Account: {config.slurm_config.account}")
        print(f"  Time limit: {config.slurm_config.time_limit}")
        print()

        if config.dry_run:
            print("*** DRY RUN MODE - Scripts will be created but not submitted ***")
            print()

    def _print_completion_summary(self) -> None:
        """Print a summary after submission."""
        config = self._dc_config
        total_jobs = sum(len(chain) for chain in self._job_chains.values())

        if config.dry_run:
            print(f"\nDry run completed. {total_jobs} job scripts created.")
            print(f"Scripts saved to: {config.output_script_dir}")
            print("Review the scripts and run without --dry-run to submit them.")
        else:
            print(f"\nAll {total_jobs} jobs submitted successfully!")
            print("\nDependency chains:")

            for replicate, results in sorted(self._job_chains.items()):
                job_ids = [r.job_id for r in results]
                print(f"  Replicate {replicate}: {' -> '.join(job_ids)}")

            print(f"\nMonitor progress with: squeue -u $USER")
            print(f"Check job details with: scontrol show job <job_id>")


def submit_daisy_chain(
    config_path: Union[str, Path],
    slurm_preset: str = "aa100",
    replicates: str = "1",
    email: str = "",
    dry_run: bool = False,
    main_script: Optional[str] = None,
    conda_env: str = "polymerist-env",
    output_dir: Optional[Union[str, Path]] = None,
    scratch_dir: Optional[Union[str, Path]] = None,
    projects_dir: Optional[Union[str, Path]] = None,
) -> Dict[int, List[SubmissionResult]]:
    """Convenience function to submit daisy-chain jobs from a YAML config.

    Args:
        config_path: Path to simulation YAML config
        slurm_preset: SLURM preset name (aa100, al40, blanca-shirts, testing)
        replicates: Replicate range string (e.g., "1-5", "1,3,5")
        email: Email for job notifications
        dry_run: If True, don't submit jobs
        main_script: Path to main simulation script
        conda_env: Conda environment name
        output_dir: Directory for job scripts (default: from config or "job_scripts")
        scratch_dir: Override scratch directory for simulation output
        projects_dir: Override projects directory for scripts/logs

    Returns:
        Dictionary mapping replicate numbers to submission results

    Example:
        >>> results = submit_daisy_chain(
        ...     config_path="simulation.yaml",
        ...     slurm_preset="aa100",
        ...     replicates="1-5",
        ...     email="user@example.com",
        ...     dry_run=True,
        ... )
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
    # Cast to PresetType for type checker (validated by argparse choices)
    from polyzymd.workflow.slurm import PresetType

    slurm_config = SlurmConfig.from_preset(slurm_preset, email=email)  # type: ignore[arg-type]

    # Create daisy-chain config
    dc_config = DaisyChainConfig.from_simulation_config(
        sim_config=sim_config,
        slurm_config=slurm_config,
        replicates=replicates,
        dry_run=dry_run,
        output_script_dir=script_output_dir,
        main_script=main_script,
        main_script_args={"config": str(config_path)},
    )

    # Create submitter and submit
    submitter = DaisyChainSubmitter(sim_config, dc_config, conda_env=conda_env)
    return submitter.submit_all()


def main() -> int:
    """Main entry point for daisy-chain submission CLI.

    Returns:
        Exit code (0 for success, 1 for failure).
    """
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Submit daisy-chained MD simulation jobs to SLURM")

    parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=True,
        help="Path to simulation YAML configuration file",
    )
    parser.add_argument(
        "-r",
        "--replicates",
        type=str,
        default="1",
        help="Replicate range (e.g., '1-5', '1,3,5'). Default: 1",
    )
    parser.add_argument(
        "--preset",
        type=str,
        choices=["aa100", "al40", "blanca-shirts", "testing"],
        default="aa100",
        help="SLURM partition preset. Default: aa100",
    )
    parser.add_argument(
        "--email",
        type=str,
        default="",
        help="Email for job notifications",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Generate scripts but don't submit them",
    )
    parser.add_argument(
        "--main-script",
        type=str,
        default=None,
        help="Path to main simulation script for initial job",
    )
    parser.add_argument(
        "--conda-env",
        type=str,
        default="polymerist-env",
        help="Conda environment name. Default: polymerist-env",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="daisy_chain_scripts",
        help="Output directory for job scripts. Default: daisy_chain_scripts",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args()

    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    try:
        submit_daisy_chain(
            config_path=args.config,
            slurm_preset=args.preset,
            replicates=args.replicates,
            email=args.email,
            dry_run=args.dry_run,
            main_script=args.main_script,
            conda_env=args.conda_env,
            output_dir=args.output_dir,
        )
        return 0

    except FileNotFoundError as e:
        LOGGER.error(f"Configuration file not found: {e}")
        return 1

    except ValueError as e:
        LOGGER.error(f"Invalid configuration: {e}")
        return 1

    except Exception as e:
        LOGGER.error(f"Error during submission: {e}")
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    import sys

    sys.exit(main())
