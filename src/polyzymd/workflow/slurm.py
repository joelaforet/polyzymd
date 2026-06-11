"""
SLURM job script generation for HPC cluster submission.

This module provides templates and utilities for generating SLURM
batch scripts for self-resubmitting MD simulation jobs.

.. versionchanged:: 1.1.0
    Replaced conda/module-load environment activation with pixi.
    The ``module_load`` and ``conda_command`` fields on ``SlurmConfig``
    have been removed.  Environment activation is now handled by
    ``pixi shell-hook`` using the ``pixi_env`` parameter on
    ``SlurmScriptGenerator``.
"""

from __future__ import annotations

import logging
import re as _re
import shutil
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Dict, List, Literal, Optional, Union

from polyzymd.core.branding import FULL_CREDIT_LINE
from polyzymd.utils.templates import render_package_template

LOGGER = logging.getLogger(__name__)

# Preset types
PresetType = Literal["aa100", "al40", "blanca-shirts", "blanca-chbe-rdi", "bridges2", "testing"]

# Valid GPU types for Bridges2 (PSC).  Adding a new type is a one-line change here.
BRIDGES2_GPU_TYPES: List[str] = ["v100-16", "v100-32", "l40s-48", "h100-80"]

# Mapping from SLURM preset → default pixi environment name.
# Used by the CLI to pick the right environment automatically when
# the user doesn't pass ``--pixi-env`` explicitly.
PRESET_DEFAULT_PIXI_ENV: Dict[str, str] = {
    "aa100": "cuda-12-4",
    "al40": "cuda-12-4",
    "blanca-shirts": "cuda-12-4",
    "blanca-chbe-rdi": "cuda-12-4",
    "bridges2": "cuda-12-6",
    "testing": "cuda-12-4",
}

# Pattern allowing alphanumerics, common path chars, and SLURM-safe punctuation
# Intentionally excludes shell metacharacters, including bare $ and backslash
_SAFE_SCRIPT_VALUE = _re.compile(r"^[A-Za-z0-9._/,:\-@%=+ ]+$")
_SAFE_CONSTRAINT_VALUE = _re.compile(r"^[A-Za-z0-9._\-|&]+$")
_SAFE_NODELIST_VALUE = _re.compile(r"^[A-Za-z0-9._,\-\[\]]+$")
_SAFE_GPU_TYPE_VALUE = _re.compile(r"^[A-Za-z0-9._\-]+$")
_WORKFLOW_TEMPLATE_PACKAGE = "polyzymd.workflow"
_OPENMM_SELF_RESUBMITTING_TEMPLATE = "openmm_self_resubmitting.sh.jinja"


def _load_workflow_template_source(template_name: str) -> str:
    """Load a workflow template resource as text.

    Parameters
    ----------
    template_name : str
        Template filename within the workflow template resource directory.

    Returns
    -------
    str
        Raw template source text.
    """
    return (
        resources.files(_WORKFLOW_TEMPLATE_PACKAGE).joinpath("templates", template_name).read_text()
    )


def _validate_script_value(value: str, field_name: str) -> str:
    """Reject values containing shell metacharacters before bash interpolation.

    Parameters
    ----------
    value : str
        The value to validate.
    field_name : str
        Name of the field for error messages.

    Returns
    -------
    str
        The validated value unchanged.

    Raises
    ------
    ValueError
        If the value contains unsafe characters.
    """
    if value and not _SAFE_SCRIPT_VALUE.match(value):
        raise ValueError(
            f"SLURM script field '{field_name}' contains unsafe characters: {value!r}. "
            "Only alphanumerics, spaces, and -_./:,@%=+ are allowed."
        )
    return value


def _validate_constraint_value(value: str, field_name: str) -> str:
    """Validate a SLURM --constraint value, allowing ``|`` (OR) and ``&`` (AND).

    SLURM constraint expressions use ``|`` and ``&`` as boolean operators
    (e.g. ``"A40|A100"``), which are deliberately forbidden by
    ``_validate_script_value``.  This helper uses a separate regex that
    permits those two characters while still rejecting dangerous shell
    metacharacters (``; $ ` ( ) { } < > ' " \\ !``).

    Parameters
    ----------
    value : str
        The constraint value to validate.
    field_name : str
        Name of the field for error messages.

    Returns
    -------
    str
        The validated value unchanged.

    Raises
    ------
    ValueError
        If the value contains unsafe characters.
    """
    if value and not _SAFE_CONSTRAINT_VALUE.match(value):
        raise ValueError(
            f"SLURM constraint field '{field_name}' contains unsafe characters: {value!r}. "
            "Only alphanumerics, hyphens, dots, underscores, | (OR), and & (AND) are allowed."
        )
    return value


def _validate_nodelist_value(value: str, field_name: str = "nodelist") -> str:
    """Validate a SLURM nodelist value, allowing bracket hostlist syntax.

    Parameters
    ----------
    value : str
        The nodelist string to validate.
    field_name : str, optional
        Field name for error messages, by default "nodelist".

    Returns
    -------
    str
        The validated value.

    Raises
    ------
    ValueError
        If the value contains unsafe characters.
    """
    if value and not _SAFE_NODELIST_VALUE.match(value):
        raise ValueError(
            f"{field_name} contains unsafe characters: {value!r}. "
            "Only alphanumeric, '.', '_', '-', ',', '[', ']' are allowed."
        )
    return value


def _validate_gpu_type_value(value: str, field_name: str = "gpu_type") -> str:
    """Validate a GPU type string for SBATCH GRES rendering.

    Parameters
    ----------
    value : str
        GPU type value to validate.
    field_name : str, optional
        Field name used in validation error messages.

    Returns
    -------
    str
        The validated value.

    Raises
    ------
    ValueError
        If the value contains unsafe characters.
    """
    if value and not _SAFE_GPU_TYPE_VALUE.match(value):
        raise ValueError(
            f"{field_name} contains unsafe characters: {value!r}. "
            "Only alphanumeric, '.', '_', '-' are allowed."
        )
    return value


def _validate_positive_replicate(value: object) -> int:
    """Validate the replicate number used in generated shell scripts.

    Parameters
    ----------
    value : object
        Replicate value supplied by the caller.

    Returns
    -------
    int
        Positive integer replicate value.

    Raises
    ------
    ValueError
        If the replicate is not a positive integer.
    """
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"replicate must be a positive integer, got {value!r}")
    return value


def _discover_manifest_path() -> str:
    """Auto-detect the pixi workspace manifest (``pixi.toml``).

    Strategy: find the ``polyzymd`` executable on ``$PATH`` (installed by
    pixi into ``.pixi/envs/<name>/bin/``), walk up the directory tree from
    the resolved binary to locate ``pixi.toml``.

    Returns
    -------
    str
        Absolute path to ``pixi.toml``.

    Raises
    ------
    RuntimeError
        If the manifest cannot be found.
    """
    exe = shutil.which("polyzymd")
    if exe is None:
        raise RuntimeError(
            "Cannot auto-detect pixi manifest: 'polyzymd' is not on PATH. "
            "Are you running inside a pixi environment?"
        )
    # Resolve symlinks → .pixi/envs/<env>/bin/polyzymd
    exe_path = Path(exe).resolve()
    # Walk up looking for pixi.toml
    candidate = exe_path.parent
    for _ in range(10):  # prevent infinite traversal
        manifest = candidate / "pixi.toml"
        if manifest.is_file():
            return str(manifest)
        if candidate.parent == candidate:
            break
        candidate = candidate.parent

    raise RuntimeError(
        f"Cannot auto-detect pixi manifest: walked up from {exe_path} "
        "but no pixi.toml found. Ensure you cloned the polyzymd repo and "
        "installed via 'pixi install'."
    )


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
        cpus_per_task: Number of CPUs allocated per task.
        memory: Memory allocation (e.g. ``"3G"``).  Set to ``None`` to omit the
            ``--mem`` directive entirely (some clusters allocate memory per GPU
            and reject an explicit ``--mem`` request).
        gpus: Number of GPUs.
        exclude: Nodes to exclude (omitted when ``None``).
        nodelist: Optional SLURM ``--nodelist`` value.
        gpu_type: Optional GPU type string used with the ``--gpus`` directive
            (e.g. ``"v100-32"`` for Bridges2).  When ``None`` the classic
            ``--gres=gpu:<N>`` directive is emitted instead.
        gpu_directive_style: ``"gres"`` (default, Alpine-style) or ``"gpus"``
            (Bridges2-style).  Controls which SBATCH GPU directive is written.
            Also governs which nodes/ntasks format is emitted.
        constraint: Optional SLURM ``--constraint`` expression. Supports
            boolean expressions with ``|`` (OR) and ``&`` (AND), such as
            ``"A40|A100"``.
    """

    partition: str = "aa100"
    qos: str = "normal"
    account: str = "ucb625_asc1"
    time_limit: str = "23:59:59"
    email: str = ""
    nodes: int = 1
    ntasks: int = 1
    cpus_per_task: int = 1
    memory: Optional[str] = "3G"
    gpus: int = 1
    exclude: Optional[str] = None
    nodelist: Optional[str] = None
    # --- GPU directive fields ---
    gpu_type: Optional[str] = None
    gpu_directive_style: str = "gres"
    constraint: Optional[str] = None

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
            "blanca-chbe-rdi": {
                "partition": "blanca,blanca-chbe-rdi",
                "qos": "preemptable",
                "account": "blanca-chbe-rdi",
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
            },
            "testing": {
                "partition": "atesting_a100",
                "qos": "testing",
                "account": "ucb625_asc1",
                "time_limit": "0:05:59",
            },
        }

        if preset not in presets:
            raise ValueError(
                f"Unknown SLURM preset {preset!r}. Valid presets: {', '.join(sorted(presets))}"
            )
        config_dict = presets[preset]
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


class SlurmScriptGenerator:
    """Generator for SLURM batch scripts.

    Supports separate directories for:
    - projects_dir: Where scripts live and jobs are submitted from
    - scratch_dir: Where simulation output goes (trajectories, checkpoints)

    Example:
        >>> config = SlurmConfig.from_preset("aa100", email="user@example.com")
        >>> generator = SlurmScriptGenerator(config)
        >>> script = generator.generate_job_script(
        ...     config_path="/projects/user/config.yaml",
        ...     replicate=1,
        ...     working_dir="/scratch/user/sim_output",
        ... )
    """

    # =====================================================================
    # Self-resubmitting job template
    # =====================================================================
    # Every SLURM job runs the same script: it calls `polyzymd run-segment`
    # which inspects progress state, runs the next segment, and exits.
    # The bash wrapper then checks whether more work remains and resubmits
    # itself via `sbatch "$SLURM_JOB_SCRIPT"`.
    #
    # Exit codes from `polyzymd run-segment`:
    #   0  — segment completed normally (may or may not be final)
    #   2  — concurrent execution detected (another job already running);
    #        this duplicate chain should terminate without resubmitting
    #   99 — graceful interruption (wall-time signal); should resubmit
    #   other — unexpected failure; do NOT resubmit
    # =====================================================================
    JOB_TEMPLATE = _load_workflow_template_source(_OPENMM_SELF_RESUBMITTING_TEMPLATE)

    def __init__(
        self,
        config: SlurmConfig,
        pixi_env: str = "cuda-12-4",
        openff_logs: bool = False,
        skip_build: bool = False,
    ) -> None:
        """Initialize the generator.

        Args:
            config: SLURM configuration.
            pixi_env: Pixi environment name (e.g. ``"cuda-12-4"``, ``"cuda-12-6"``).
            openff_logs: Enable verbose OpenFF logs in generated scripts.
            skip_build: Skip system building in generated scripts (use pre-built system).
        """
        self._config = config
        self._pixi_env = pixi_env
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
        if self._config.gpus == 0:
            return ""
        if self._config.gpu_directive_style == "gpus" and self._config.gpu_type:
            return f"#SBATCH --gpus={self._config.gpu_type}:{self._config.gpus}"
        if self._config.gpu_type:
            return f"#SBATCH --gres=gpu:{self._config.gpu_type}:{self._config.gpus}"
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
        if self._config.gpu_directive_style == "gpus" and self._config.gpus > 0:
            return f"#SBATCH -N {self._config.nodes}"
        return f"#SBATCH --nodes={self._config.nodes}\n#SBATCH --ntasks={self._config.ntasks}"

    def _cpus_line(self) -> str:
        """Return the CPUs-per-task directive when requested."""
        if self._config.cpus_per_task > 1:
            return f"#SBATCH --cpus-per-task={self._config.cpus_per_task}"
        return ""

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

    def _nodelist_line(self) -> str:
        """Return the nodelist SBATCH directive, or empty string when unset."""
        return f"#SBATCH --nodelist={self._config.nodelist}" if self._config.nodelist else ""

    def _constraint_line(self) -> str:
        """Return the constraint SBATCH directive, or an empty string to omit it."""
        return f"#SBATCH --constraint={self._config.constraint}" if self._config.constraint else ""

    # ------------------------------------------------------------------
    # Self-resubmitting job generation
    # ------------------------------------------------------------------

    def generate_job_script(
        self,
        config_path: str,
        replicate: int,
        working_dir: str,
        job_name: str | None = None,
        output_file: str | None = None,
    ) -> str:
        """Generate a self-resubmitting SLURM job script.

        This produces a single script that handles the entire simulation
        lifecycle.  Each invocation calls ``polyzymd run-segment`` which
        determines what work remains, runs the next segment, and exits.
        The bash wrapper then checks progress and resubmits itself if
        more work is needed.

        Parameters
        ----------
        config_path : str
            Absolute path to the YAML configuration file.
        replicate : int
            Replicate number.
        working_dir : str
            Directory for simulation output (trajectories, checkpoints).
        job_name : str or None, optional
            SLURM job name.  Callers should use
            :func:`~polyzymd.workflow.daisy_chain.create_job_name` to
            produce descriptive names (e.g. ``r1_310K_Fibronectin_...``).
            Falls back to ``pzmd_r{replicate}`` if not provided.
        output_file : str or None, optional
            SLURM log file pattern.  Falls back to
            ``slurm_logs/{job_name}.%j.out`` relative to the directory
            where ``sbatch`` is invoked.

        Returns
        -------
        str
            Complete SLURM batch script content.
        """
        replicate = _validate_positive_replicate(replicate)

        if job_name is None:
            job_name = f"pzmd_r{replicate}"
        if output_file is None:
            output_file = f"slurm_logs/{job_name}.%j.out"

        openff_logs_flag = " --openff-logs" if self._openff_logs else ""
        skip_build_flag = " \\\n    --skip-build" if self._skip_build else ""

        # Auto-detect the pixi manifest path from the current installation.
        manifest_path = _discover_manifest_path()

        # Validate all interpolated string values against shell injection
        _validate_script_value(self._config.partition, "partition")
        _validate_script_value(job_name, "job_name")
        _validate_script_value(output_file, "output_file")
        _validate_script_value(self._config.time_limit, "time_limit")
        _validate_script_value(self._pixi_env, "pixi_env")
        _validate_script_value(str(manifest_path), "manifest_path")
        _validate_script_value(str(config_path), "config_path")
        _validate_script_value(str(working_dir), "working_dir")

        if self._config.qos:
            _validate_script_value(self._config.qos, "qos")
        if self._config.memory:
            _validate_script_value(self._config.memory, "memory")
        if self._config.account:
            _validate_script_value(self._config.account, "account")
        if self._config.email:
            _validate_script_value(self._config.email, "email")
        if self._config.exclude:
            _validate_script_value(self._config.exclude, "exclude")
        if self._config.nodelist:
            _validate_nodelist_value(self._config.nodelist, "nodelist")
        if self._config.constraint:
            _validate_constraint_value(self._config.constraint, "constraint")
        if self._config.gpu_type:
            _validate_gpu_type_value(self._config.gpu_type, "gpu_type")

        return render_package_template(
            _WORKFLOW_TEMPLATE_PACKAGE,
            _OPENMM_SELF_RESUBMITTING_TEMPLATE,
            {
                "partition": self._config.partition,
                "job_name": job_name,
                "output_file": output_file,
                "qos_line": self._qos_line(),
                "nodes_line": self._nodes_line(),
                "cpus_line": self._cpus_line(),
                "mem_line": self._mem_line(),
                "time_limit": self._config.time_limit,
                "gpu_line": self._gpu_line(),
                "mail_line": self._mail_line(),
                "account_line": self._account_line(),
                "exclude_line": self._exclude_line(),
                "nodelist_line": self._nodelist_line(),
                "constraint_line": self._constraint_line(),
                "pixi_env": self._pixi_env,
                "manifest_path": manifest_path,
                "config_path": config_path,
                "replicate": replicate,
                "working_dir": working_dir,
                "openff_logs_flag": openff_logs_flag,
                "skip_build_flag": skip_build_flag,
                "FULL_CREDIT_LINE": FULL_CREDIT_LINE,
            },
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
