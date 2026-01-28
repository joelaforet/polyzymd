"""
SLURM job script generation for HPC cluster submission.

This module provides templates and utilities for generating SLURM
batch scripts for MD simulations.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Literal, Optional, Union

LOGGER = logging.getLogger(__name__)

# Preset types
PresetType = Literal["aa100", "al40", "blanca-shirts", "testing"]


@dataclass
class SlurmConfig:
    """Configuration for SLURM job submission.

    Attributes:
        partition: SLURM partition(s) to use.
        qos: Quality of service.
        account: Account for resource allocation.
        time_limit: Wall time limit (HH:MM:SS).
        email: Email for notifications.
        nodes: Number of nodes.
        ntasks: Number of tasks.
        memory: Memory allocation (e.g., "3G").
        gpus: Number of GPUs.
        exclude: Nodes to exclude.
    """

    partition: str = "aa100"
    qos: str = "normal"
    account: str = "ucb625_asc1"
    time_limit: str = "23:59:59"
    email: str = ""
    nodes: int = 1
    ntasks: int = 1
    memory: str = "3G"
    gpus: int = 1
    exclude: Optional[str] = None

    @classmethod
    def from_preset(cls, preset: PresetType, email: str = "") -> "SlurmConfig":
        """Create a SlurmConfig from a preset.

        Args:
            preset: Preset name.
            email: Email for notifications.

        Returns:
            SlurmConfig with preset values.
        """
        presets: Dict[PresetType, Dict] = {
            "aa100": {
                "partition": "aa100",
                "qos": "normal",
                "account": "ucb625_asc1",
                "time_limit": "23:59:59",
            },
            "al40": {
                "partition": "al40",
                "qos": "normal",
                "account": "ucb625_asc1",
                "time_limit": "23:59:59",
            },
            "blanca-shirts": {
                "partition": "blanca,blanca-shirts",
                "qos": "preemptable",
                "account": "blanca-shirts",
                "time_limit": "23:59:59",
                "exclude": "bgpu-bortz1",
            },
            "testing": {
                "partition": "atesting_a100",
                "qos": "testing",
                "account": "ucb625_asc1",
                "time_limit": "0:05:59",
            },
        }

        config_dict = presets.get(preset, presets["aa100"])
        return cls(email=email, **config_dict)


@dataclass
class JobContext:
    """Context for job script template rendering.

    Attributes:
        job_name: SLURM job name.
        output_file: Output file pattern (for SLURM logs).
        scratch_dir: Directory for simulation output (trajectories, checkpoints).
        projects_dir: Directory for scripts and logs.
        segment_index: Current segment index.
        replicate_num: Replicate number.
        extra_vars: Additional template variables.
    """

    job_name: str
    output_file: str
    scratch_dir: str  # Where simulation data goes (trajectories, checkpoints)
    projects_dir: str = "."  # Where scripts and logs live
    segment_index: int = 0
    replicate_num: int = 1
    extra_vars: Dict = field(default_factory=dict)

    # Legacy alias for backwards compatibility
    @property
    def working_dir(self) -> str:
        """Alias for scratch_dir for backwards compatibility."""
        return self.scratch_dir


class SlurmScriptGenerator:
    """Generator for SLURM batch scripts.

    Supports separate directories for:
    - projects_dir: Where scripts live and jobs are submitted from
    - scratch_dir: Where simulation output goes (trajectories, checkpoints)

    Example:
        >>> config = SlurmConfig.from_preset("aa100", email="user@example.com")
        >>> generator = SlurmScriptGenerator(config)
        >>> script = generator.generate_initial_job(
        ...     context=JobContext(
        ...         job_name="my_sim",
        ...         output_file="logs/output.log",
        ...         scratch_dir="/scratch/user/sim_output",
        ...         projects_dir="/projects/user/polyzymd",
        ...     ),
        ...     python_script="run_simulation.py",
        ...     python_args={"temperature": 300},
        ... )
    """

    # Template for initial simulation jobs
    # - Job is submitted from projects_dir
    # - SLURM logs go to projects_dir/slurm_logs/
    # - Simulation output goes to scratch_dir
    INITIAL_JOB_TEMPLATE = """#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --job-name=i_{job_name}
#SBATCH --output={output_file}
#SBATCH --qos={qos}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --mem={memory}
#SBATCH --time={time_limit}
#SBATCH --gres=gpu:{gpus}
#SBATCH --mail-type=FAIL
#SBATCH --mail-user={email}
#SBATCH --account={account}
{exclude_line}

# Exit immediately if any command fails
set -e

# =============================================================================
# PolyzyMD Initial Simulation Job
# Segment: {segment_index}
# =============================================================================

module purge
module load miniforge

# Initialize conda/mamba for non-interactive shell
eval "$(conda shell.bash hook)"
mamba activate {conda_env}

# Projects directory (scripts, configs, logs)
PROJECTS_DIR="{projects_dir}"

# Scratch directory (simulation output)
SCRATCH_DIR="{scratch_dir}"

# Ensure scratch directory exists
mkdir -p "$SCRATCH_DIR"

# Change to projects directory where config and scripts live
cd "$PROJECTS_DIR"

echo "Starting initial simulation segment {segment_index}"
echo "Projects dir: $PROJECTS_DIR"
echo "Scratch dir: $SCRATCH_DIR"
echo "Config: {config_path}"
echo "Replicate: {replicate}"
echo "Timestamp: $(date)"

# Run the initial simulation using polyzymd CLI
# This builds the system, runs equilibration, and runs the first production segment
polyzymd run -c "{config_path}" \\
    --replicate {replicate} \\
    --scratch-dir "$SCRATCH_DIR" \\
    --segment-time {segment_time} \\
    --segment-frames {segment_frames}

echo "Segment {segment_index} completed successfully at $(date)"
"""

    # Template for continuation jobs
    CONTINUATION_JOB_TEMPLATE = """#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --job-name=c_{job_name}
#SBATCH --output={output_file}
#SBATCH --qos={qos}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --mem={memory}
#SBATCH --time={time_limit}
#SBATCH --gres=gpu:{gpus}
#SBATCH --mail-type=FAIL
#SBATCH --mail-user={email}
#SBATCH --account={account}
{exclude_line}

# Exit immediately if any command fails
set -e

# =============================================================================
# PolyzyMD Continuation Job
# Segment: {segment_index}
# =============================================================================

module purge
module load miniforge

# Initialize conda/mamba for non-interactive shell
eval "$(conda shell.bash hook)"
mamba activate {conda_env}

# Projects directory (scripts, configs, logs)
PROJECTS_DIR="{projects_dir}"

# Scratch directory (simulation output - where previous segment data lives)
SCRATCH_DIR="{scratch_dir}"

# Change to projects directory
cd "$PROJECTS_DIR"

echo "Starting continuation segment {segment_index}"
echo "Projects dir: $PROJECTS_DIR"
echo "Scratch dir: $SCRATCH_DIR"
echo "Timestamp: $(date)"

# Continue simulation from previous segment using polyzymd CLI
# Reads checkpoint from previous segment in SCRATCH_DIR
# Writes new trajectory and checkpoint to SCRATCH_DIR
polyzymd continue \\
    -w "$SCRATCH_DIR" \\
    -s {segment_index} \\
    -t {segment_time} \\
    -n {num_samples}

echo "Segment {segment_index} completed successfully at $(date)"
"""

    def __init__(
        self,
        config: SlurmConfig,
        conda_env: str = "polymerist-env",
    ) -> None:
        """Initialize the generator.

        Args:
            config: SLURM configuration.
            conda_env: Conda environment name.
        """
        self._config = config
        self._conda_env = conda_env

    @property
    def config(self) -> SlurmConfig:
        """Get the SLURM configuration."""
        return self._config

    def generate_initial_job(
        self,
        context: JobContext,
        config_path: str,
        replicate: int,
        segment_time: float,
        segment_frames: int,
    ) -> str:
        """Generate an initial simulation job script.

        Args:
            context: Job context information.
            config_path: Path to the YAML configuration file.
            replicate: Replicate number.
            segment_time: Duration of first segment in nanoseconds.
            segment_frames: Number of frames to save in first segment.

        Returns:
            SLURM batch script content.
        """
        # Format exclude line
        exclude_line = ""
        if self._config.exclude:
            exclude_line = f"#SBATCH --exclude={self._config.exclude}"

        # Use context.projects_dir
        projects_dir = context.projects_dir if context.projects_dir != "." else "."

        return self.INITIAL_JOB_TEMPLATE.format(
            partition=self._config.partition,
            job_name=context.job_name,
            output_file=context.output_file,
            qos=self._config.qos,
            nodes=self._config.nodes,
            ntasks=self._config.ntasks,
            memory=self._config.memory,
            time_limit=self._config.time_limit,
            gpus=self._config.gpus,
            email=self._config.email,
            account=self._config.account,
            exclude_line=exclude_line,
            conda_env=self._conda_env,
            projects_dir=projects_dir,
            scratch_dir=context.scratch_dir,
            config_path=config_path,
            replicate=replicate,
            segment_time=segment_time,
            segment_frames=segment_frames,
            segment_index=context.segment_index,
        )

    def generate_continuation_job(
        self,
        context: JobContext,
        segment_time: float,
        num_samples: int,
    ) -> str:
        """Generate a continuation job script.

        Args:
            context: Job context information.
            segment_time: Duration of this segment in nanoseconds.
            num_samples: Number of frames to save.

        Returns:
            SLURM batch script content.
        """
        exclude_line = ""
        if self._config.exclude:
            exclude_line = f"#SBATCH --exclude={self._config.exclude}"

        # Use context.projects_dir
        projects_dir = context.projects_dir if context.projects_dir != "." else "."

        return self.CONTINUATION_JOB_TEMPLATE.format(
            partition=self._config.partition,
            job_name=context.job_name,
            output_file=context.output_file,
            qos=self._config.qos,
            nodes=self._config.nodes,
            ntasks=self._config.ntasks,
            memory=self._config.memory,
            time_limit=self._config.time_limit,
            gpus=self._config.gpus,
            email=self._config.email,
            account=self._config.account,
            exclude_line=exclude_line,
            conda_env=self._conda_env,
            projects_dir=projects_dir,
            scratch_dir=context.scratch_dir,
            segment_index=context.segment_index,
            segment_time=segment_time,
            num_samples=num_samples,
        )

    def save_script(
        self,
        script_content: str,
        output_path: Union[str, Path],
        make_executable: bool = True,
    ) -> Path:
        """Save a script to a file.

        Args:
            script_content: Script content.
            output_path: Output file path.
            make_executable: Whether to make the script executable.

        Returns:
            Path to the saved script.
        """
        import os

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w") as f:
            f.write(script_content)

        if make_executable:
            os.chmod(output_path, 0o755)

        LOGGER.info(f"Saved script to {output_path}")
        return output_path


def parse_replicate_range(replicate_range: str) -> List[int]:
    """Parse a SLURM array range into a list of replicate numbers.

    Args:
        replicate_range: SLURM array format (e.g., "1-5", "1,3,5", "1-10:2").

    Returns:
        List of replicate numbers.

    Example:
        >>> parse_replicate_range("1-5")
        [1, 2, 3, 4, 5]
        >>> parse_replicate_range("1,3,5")
        [1, 3, 5]
        >>> parse_replicate_range("1-10:2")
        [1, 3, 5, 7, 9]
    """
    replicates = []

    parts = replicate_range.split(",")

    for part in parts:
        part = part.strip()
        if "-" in part:
            if ":" in part:
                range_part, step = part.split(":")
                step = int(step)
            else:
                range_part = part
                step = 1

            start, end = map(int, range_part.split("-"))
            replicates.extend(range(start, end + 1, step))
        else:
            replicates.append(int(part))

    return sorted(list(set(replicates)))


def validate_replicate_range(replicate_range: str) -> bool:
    """Validate that a replicate range is in proper SLURM array format.

    Args:
        replicate_range: Range string to validate.

    Returns:
        True if valid.

    Raises:
        ValueError: If the format is invalid.
    """
    import re

    pattern = r"^(\d+(-\d+(:\d+)?)?)(,\d+(-\d+(:\d+)?)?)*$"
    if not re.match(pattern, replicate_range):
        raise ValueError(f"Invalid replicate range format: {replicate_range}")
    return True
