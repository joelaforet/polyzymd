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
PresetType = Literal["aa100", "al40", "blanca-shirts", "bridges2", "testing"]

# Valid GPU types for Bridges2 (PSC).  Adding a new type is a one-line change here.
BRIDGES2_GPU_TYPES: List[str] = ["v100-16", "v100-32", "l40s-48", "h100-80"]


@dataclass
class SlurmConfig:
    """Configuration for SLURM job submission.

    Attributes:
        partition: SLURM partition(s) to use.
        qos: Quality of service. Set to ``""`` to omit the ``--qos`` directive
            entirely (required for clusters such as Bridges2 that do not use QoS).
        account: Account / allocation ID for resource allocation.  Set to ``""``
            to omit the ``--account`` directive entirely (e.g. Bridges2, which
            infers the allocation from the submitting user's login).
        time_limit: Wall time limit (HH:MM:SS).
        email: Email address for SLURM failure notifications.  Set to ``""`` to
            omit both ``--mail-type`` and ``--mail-user`` directives.
        nodes: Number of nodes.
        ntasks: Number of tasks.  Ignored when ``gpu_directive_style == "gpus"``
            (Bridges2-style); those scripts emit ``#SBATCH -N {nodes}`` only.
        memory: Memory allocation (e.g. ``"3G"``).  Set to ``None`` to omit the
            ``--mem`` directive entirely (some clusters allocate memory per GPU
            and reject an explicit ``--mem`` request).
        gpus: Number of GPUs.
        exclude: Nodes to exclude (omitted when ``None``).
        gpu_type: Optional GPU type string used with the ``--gpus`` directive
            (e.g. ``"v100-32"`` for Bridges2).  When ``None`` the classic
            ``--gres=gpu:<N>`` directive is emitted instead.
        gpu_directive_style: ``"gres"`` (default, Alpine-style) or ``"gpus"``
            (Bridges2-style).  Controls which SBATCH GPU directive is written.
            Also governs which nodes/ntasks format is emitted.
        module_load: Module name passed to ``ml`` in the job script
            (e.g. ``"miniforge"`` for Alpine, ``"anaconda3/2024.10-1"`` for
            Bridges2).
        conda_command: Conda frontend used to activate the environment
            (``"mamba"`` for Alpine, ``"conda"`` for Bridges2).
    """

    partition: str = "aa100"
    qos: str = "normal"
    account: str = "ucb625_asc1"
    time_limit: str = "23:59:59"
    email: str = ""
    nodes: int = 1
    ntasks: int = 1
    memory: Optional[str] = "3G"
    gpus: int = 1
    exclude: Optional[str] = None
    # --- GPU directive fields (new in v1.0.1) ---
    gpu_type: Optional[str] = None
    gpu_directive_style: str = "gres"
    # --- Module / conda fields (new in v1.0.1) ---
    module_load: str = "miniforge"
    conda_command: str = "mamba"

    @classmethod
    def from_preset(cls, preset: PresetType, email: str = "") -> "SlurmConfig":
        """Create a SlurmConfig from a named preset.

        Args:
            preset: Preset name.
            email: Email for notifications.

        Returns:
            SlurmConfig with preset values.
        """
        presets: Dict[str, Dict] = {
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
            "bridges2": {
                "partition": "GPU-shared",
                # Bridges2 does not use QoS — omit the directive entirely.
                "qos": "",
                # Bridges2 infers allocation from the submitting user's login;
                # omit the --account directive entirely.
                "account": "",
                "time_limit": "24:00:00",
                # GPU-shared allocates resources per GPU; explicit --mem is
                # not required and may be rejected.  Set to None to omit.
                "memory": None,
                # Use the newer --gpus=<type>:<n> SBATCH syntax (also selects
                # -N 1 nodes format instead of --nodes + --ntasks).
                "gpu_type": "v100-32",
                "gpu_directive_style": "gpus",
                # Bridges2 uses anaconda3 + conda rather than miniforge + mamba.
                "module_load": "anaconda3/2024.10-1",
                "conda_command": "conda",
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
    #
    # {qos_line}, {mem_line}, {gpu_line}, {nodes_line}, {account_line},
    # {mail_line} are computed conditionally by helper methods so that clusters
    # which omit specific directives (e.g. Bridges2 omits --qos, --mem,
    # --account, and --ntasks) produce valid scripts without empty lines.
    INITIAL_JOB_TEMPLATE = """#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --job-name=i_{job_name}
#SBATCH --output={output_file}
{qos_line}
{nodes_line}
{mem_line}
#SBATCH --time={time_limit}
{gpu_line}
{mail_line}
{account_line}
{exclude_line}

# =============================================================================
# PolyzyMD Initial Simulation Job
# Segment: {segment_index}
# =============================================================================

# Load conda environment (ignore module warnings on some HPC systems)
module purge 2>/dev/null || true
ml {module_load} 2>/dev/null || true

# Initialize conda/mamba for non-interactive shell
eval "$(conda shell.bash hook)"
{conda_command} activate {conda_env}

# Enable strict error handling after environment setup
set -e

# Required for OpenFF Interchange.combine() functionality
export INTERCHANGE_EXPERIMENTAL=1

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
polyzymd{openff_logs_flag} run -c "{config_path}" \\
    --replicate {replicate} \\
    --scratch-dir "$SCRATCH_DIR" \\
    --segment-time {segment_time} \\
    --segment-frames {segment_frames}{skip_build_flag}

echo "Segment {segment_index} completed successfully at $(date)"
"""

    # Template for continuation jobs
    CONTINUATION_JOB_TEMPLATE = """#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --job-name=c_{job_name}
#SBATCH --output={output_file}
{qos_line}
{nodes_line}
{mem_line}
#SBATCH --time={time_limit}
{gpu_line}
{mail_line}
{account_line}
{exclude_line}

# =============================================================================
# PolyzyMD Continuation Job
# Segment: {segment_index}
# =============================================================================

# Load conda environment (ignore module warnings on some HPC systems)
module purge 2>/dev/null || true
ml {module_load} 2>/dev/null || true

# Initialize conda/mamba for non-interactive shell
eval "$(conda shell.bash hook)"
{conda_command} activate {conda_env}

# Enable strict error handling after environment setup
set -e

# Required for OpenFF Interchange.combine() functionality
export INTERCHANGE_EXPERIMENTAL=1

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
polyzymd{openff_logs_flag} continue \\
    -w "$SCRATCH_DIR" \\
    -s {segment_index} \\
    -t {segment_time} \\
    -n {num_samples}

echo "Segment {segment_index} completed successfully at $(date)"
"""

    def __init__(
        self,
        config: SlurmConfig,
        conda_env: str = "polyzymd-env",
        openff_logs: bool = False,
        skip_build: bool = False,
    ) -> None:
        """Initialize the generator.

        Args:
            config: SLURM configuration.
            conda_env: Conda environment name.
            openff_logs: Enable verbose OpenFF logs in generated scripts.
            skip_build: Skip system building in generated scripts (use pre-built system).
        """
        self._config = config
        self._conda_env = conda_env
        self._openff_logs = openff_logs
        self._skip_build = skip_build

    @property
    def config(self) -> SlurmConfig:
        """Get the SLURM configuration."""
        return self._config

    # ------------------------------------------------------------------
    # Internal helpers — compute optional SBATCH directive lines
    # ------------------------------------------------------------------

    def _gpu_line(self) -> str:
        """Return the appropriate GPU SBATCH directive for this config.

        Returns ``#SBATCH --gpus=<type>:<n>`` for clusters that use the newer
        ``--gpus`` syntax (e.g. Bridges2), or ``#SBATCH --gres=gpu:<n>`` for
        clusters that use the classic Generic RESources syntax (Alpine).
        """
        if self._config.gpu_directive_style == "gpus" and self._config.gpu_type:
            return f"#SBATCH --gpus={self._config.gpu_type}:{self._config.gpus}"
        return f"#SBATCH --gres=gpu:{self._config.gpus}"

    def _nodes_line(self) -> str:
        """Return the nodes/tasks SBATCH directive(s) appropriate for this config.

        Alpine-style (``gpu_directive_style == "gres"``) emits two lines::

            #SBATCH --nodes=N
            #SBATCH --ntasks=N

        Bridges2-style (``gpu_directive_style == "gpus"``) emits a single
        short-flag line::

            #SBATCH -N N
        """
        if self._config.gpu_directive_style == "gpus":
            return f"#SBATCH -N {self._config.nodes}"
        return f"#SBATCH --nodes={self._config.nodes}\n#SBATCH --ntasks={self._config.ntasks}"

    def _qos_line(self) -> str:
        """Return the QoS SBATCH directive, or an empty string to omit it."""
        return f"#SBATCH --qos={self._config.qos}" if self._config.qos else ""

    def _mem_line(self) -> str:
        """Return the memory SBATCH directive, or an empty string to omit it."""
        return f"#SBATCH --mem={self._config.memory}" if self._config.memory else ""

    def _account_line(self) -> str:
        """Return the account SBATCH directive, or an empty string to omit it.

        An empty account string means the cluster infers the allocation from the
        submitting user's login (e.g. Bridges2).
        """
        return f"#SBATCH --account={self._config.account}" if self._config.account else ""

    def _mail_line(self) -> str:
        """Return the mail-type + mail-user SBATCH directives, or empty string.

        Both ``--mail-type`` and ``--mail-user`` are omitted together when no
        email address is configured, keeping the script clean.
        """
        if self._config.email:
            return f"#SBATCH --mail-type=FAIL\n#SBATCH --mail-user={self._config.email}"
        return ""

    def _exclude_line(self) -> str:
        """Return the exclude SBATCH directive, or an empty string to omit it."""
        return f"#SBATCH --exclude={self._config.exclude}" if self._config.exclude else ""

    # ------------------------------------------------------------------

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
        # Use context.projects_dir
        projects_dir = context.projects_dir if context.projects_dir != "." else "."

        # Format openff_logs flag
        openff_logs_flag = " --openff-logs" if self._openff_logs else ""

        # Format skip_build flag (only for initial job - continuation jobs don't build)
        skip_build_flag = " \\\n    --skip-build" if self._skip_build else ""

        return self.INITIAL_JOB_TEMPLATE.format(
            partition=self._config.partition,
            job_name=context.job_name,
            output_file=context.output_file,
            qos_line=self._qos_line(),
            nodes_line=self._nodes_line(),
            mem_line=self._mem_line(),
            time_limit=self._config.time_limit,
            gpu_line=self._gpu_line(),
            mail_line=self._mail_line(),
            account_line=self._account_line(),
            exclude_line=self._exclude_line(),
            module_load=self._config.module_load,
            conda_command=self._config.conda_command,
            conda_env=self._conda_env,
            projects_dir=projects_dir,
            scratch_dir=context.scratch_dir,
            config_path=config_path,
            replicate=replicate,
            segment_time=segment_time,
            segment_frames=segment_frames,
            segment_index=context.segment_index,
            openff_logs_flag=openff_logs_flag,
            skip_build_flag=skip_build_flag,
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
        # Use context.projects_dir
        projects_dir = context.projects_dir if context.projects_dir != "." else "."

        # Format openff_logs flag
        openff_logs_flag = " --openff-logs" if self._openff_logs else ""

        return self.CONTINUATION_JOB_TEMPLATE.format(
            partition=self._config.partition,
            job_name=context.job_name,
            output_file=context.output_file,
            qos_line=self._qos_line(),
            nodes_line=self._nodes_line(),
            mem_line=self._mem_line(),
            time_limit=self._config.time_limit,
            gpu_line=self._gpu_line(),
            mail_line=self._mail_line(),
            account_line=self._account_line(),
            exclude_line=self._exclude_line(),
            module_load=self._config.module_load,
            conda_command=self._config.conda_command,
            conda_env=self._conda_env,
            projects_dir=projects_dir,
            scratch_dir=context.scratch_dir,
            segment_index=context.segment_index,
            segment_time=segment_time,
            num_samples=num_samples,
            openff_logs_flag=openff_logs_flag,
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
