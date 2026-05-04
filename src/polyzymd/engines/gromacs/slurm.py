"""GROMACS-specific SLURM job script generation."""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path

from polyzymd.core.branding import FULL_CREDIT_LINE
from polyzymd.utils.templates import render_package_template
from polyzymd.workflow.slurm import (
    SlurmConfig,
    _discover_manifest_path,
    _validate_constraint_value,
    _validate_gpu_type_value,
    _validate_nodelist_value,
    _validate_script_value,
)

from .binary import is_mpi_binary

LOGGER = logging.getLogger(__name__)
_TEMPLATE_PACKAGE = "polyzymd.engines.gromacs"
_SLURM_TEMPLATE = "slurm_self_resubmitting.sh.jinja"


class GromacsSlurmScriptGenerator:
    """Generate self-resubmitting SLURM scripts for GROMACS simulations.

    The generated script handles the complete lifecycle:

    - Activates the requested pixi environment
    - Optionally loads an HPC module (for example ``module load gromacs/2024``)
    - Resolves the GROMACS binary and runs EM plus staged equilibration
    - Runs production with ``-maxh`` and checkpoint restart via ``-cpi``
    - Updates ``progress.json`` with ``polyzymd _update-gromacs-progress``
    - Resubmits itself with ``sbatch`` until production is complete
    """

    MAXH_SAFETY_FACTOR = 0.90
    _COMMAND_PREFIX_PATTERN = re.compile(r"^[A-Za-z0-9._/,:\-@%=+ ]+$")
    _ENV_VAR_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

    _SETUP_COMMAND_TOKEN_PATTERN = re.compile(r"^[A-Za-z0-9._/,:\-@%=+]+$")
    _MODULE_ACTIONS_WITH_ARGS = {"load", "unload", "use", "unuse", "restore"}
    _MODULE_ACTIONS_WITHOUT_ARGS = {"purge", "reset"}

    def __init__(
        self,
        slurm_config: SlurmConfig,
        pixi_env: str = "cuda-12-4",
        gmx_binary: str = "gmx",
        grompp_flags: str = "-maxwarn 1",
        mdrun_flags: str = "",
        mdrun_flags_eq: str | None = None,
        mdrun_flags_prod: str | None = None,
        command_prefix: str | None = None,
        mpi_launcher_flags: str = "",
        module_load: str | None = None,
        env_exports: dict[str, str] | None = None,
        setup_commands: list[str] | None = None,
    ) -> None:
        """Initialize the GROMACS SLURM script generator.

        Parameters
        ----------
        slurm_config : SlurmConfig
            Scheduler resource configuration.
        pixi_env : str, optional
            Pixi environment name used for job execution.
        gmx_binary : str, optional
            GROMACS executable path or command name.
        grompp_flags : str, optional
            Extra flags appended to ``gmx grompp`` invocations.
        mdrun_flags : str, optional
            Extra flags appended to ``gmx mdrun`` invocations.
        mdrun_flags_eq : str | None, optional
            Equilibration-specific mdrun flags.
        mdrun_flags_prod : str | None, optional
            Production-specific mdrun flags.
        command_prefix : str | None, optional
            Prefix prepended to all GROMACS commands.
        mpi_launcher_flags : str, optional
            Extra flags passed to mpirun for real-MPI builds.
        module_load : str | None, optional
            Optional module command executed before simulation commands.
        env_exports : dict[str, str] | None, optional
            Environment variables exported in the script prior to GROMACS commands.
        setup_commands : list[str] | None, optional
            Simple environment setup commands executed after module loading and env exports.
        """
        self._config = slurm_config
        self._pixi_env = pixi_env
        self._gmx_binary = gmx_binary
        self._is_mpi_binary = is_mpi_binary(self._gmx_binary)
        self._grompp_flags = grompp_flags
        self._mdrun_flags = mdrun_flags
        self._mdrun_flags_eq = mdrun_flags_eq
        self._mdrun_flags_prod = mdrun_flags_prod
        self._command_prefix = command_prefix
        self._mpi_launcher_flags = mpi_launcher_flags
        self._module_load = module_load
        self._env_exports = env_exports or {}
        self._setup_commands = setup_commands or []

    def _gpu_line(self) -> str:
        """Return the GPU SBATCH directive for the configured cluster style."""
        if self._config.gpus == 0:
            return ""
        if self._config.gpu_directive_style == "gpus" and self._config.gpu_type:
            return f"#SBATCH --gpus={self._config.gpu_type}:{self._config.gpus}"
        if self._config.gpu_type:
            return f"#SBATCH --gres=gpu:{self._config.gpu_type}:{self._config.gpus}"
        return f"#SBATCH --gres=gpu:{self._config.gpus}"

    def _nodes_line(self) -> str:
        """Return nodes and tasks directives for the configured cluster style."""
        if self._config.gpu_directive_style == "gpus" and self._config.gpus > 0:
            return f"#SBATCH -N {self._config.nodes}"
        return f"#SBATCH --nodes={self._config.nodes}\n#SBATCH --ntasks={self._config.ntasks}"

    def _cpus_line(self) -> str:
        """Return CPUs-per-task SBATCH directive when greater than one."""
        if self._config.cpus_per_task > 1:
            return f"#SBATCH --cpus-per-task={self._config.cpus_per_task}"
        return ""

    def _qos_line(self) -> str:
        """Return the optional QoS SBATCH directive."""
        return f"#SBATCH --qos={self._config.qos}" if self._config.qos else ""

    def _mem_line(self) -> str:
        """Return the optional memory SBATCH directive."""
        return f"#SBATCH --mem={self._config.memory}" if self._config.memory else ""

    def _account_line(self) -> str:
        """Return the optional account SBATCH directive."""
        return f"#SBATCH --account={self._config.account}" if self._config.account else ""

    def _mail_line(self) -> str:
        """Return optional SLURM email directives."""
        if self._config.email:
            return f"#SBATCH --mail-type=FAIL,END\n#SBATCH --mail-user={self._config.email}"
        return ""

    def _exclude_line(self) -> str:
        """Return optional node exclusion SBATCH directive."""
        return f"#SBATCH --exclude={self._config.exclude}" if self._config.exclude else ""

    def _constraint_line(self) -> str:
        """Return the constraint SBATCH directive, or an empty string to omit it."""
        return f"#SBATCH --constraint={self._config.constraint}" if self._config.constraint else ""

    def _nodelist_line(self) -> str:
        """Return optional nodelist SBATCH directive."""
        return f"#SBATCH --nodelist={self._config.nodelist}" if self._config.nodelist else ""

    def _validate_setup_command(self, command: str) -> None:
        """Validate setup command content before script interpolation.

        Setup commands are rendered as shell lines, so only a small allowlist of
        common environment setup forms is accepted. This preserves simple module,
        source, conda activation, and export commands while rejecting shell
        metacharacters such as command substitution, redirection, and chaining.

        Parameters
        ----------
        command : str
            Setup command to validate.

        Raises
        ------
        ValueError
            If the command is not an allowed simple setup command.
        """
        if "\n" in command or "\r" in command:
            raise ValueError(
                "setup_commands entries must be single-line commands without embedded newlines"
            )
        if "\0" in command:
            raise ValueError("setup_commands entries must not contain NUL bytes")

        stripped = command.strip()
        if not stripped:
            raise ValueError("setup_commands entries must not be empty")

        tokens = stripped.split()
        if any(not self._SETUP_COMMAND_TOKEN_PATTERN.match(token) for token in tokens):
            raise ValueError(
                "setup_commands entries may only contain alphanumerics, spaces, and "
                "-_./:,@%=+. Shell metacharacters are not allowed."
            )

        command_name = tokens[0]
        if command_name == "module":
            self._validate_module_tokens(tokens, "setup_commands")
            return

        if command_name == "source":
            if len(tokens) < 2:
                raise ValueError("setup_commands source entries require a file path")
            return

        if command_name == "conda":
            if len(tokens) != 3 or tokens[1] != "activate":
                raise ValueError("setup_commands conda entries must be 'conda activate <env>'")
            return

        if command_name == "export":
            if len(tokens) < 2:
                raise ValueError("setup_commands export entries require at least one variable")
            for assignment in tokens[1:]:
                key = assignment.split("=", maxsplit=1)[0]
                if not self._ENV_VAR_NAME.match(key):
                    raise ValueError(
                        "setup_commands export entries must use valid shell variable names"
                    )
            return

        raise ValueError(
            "setup_commands entries must start with one of: module, source, conda activate, "
            "or export"
        )

    def _validate_module_load(self, command: str) -> None:
        """Validate the dedicated module load command.

        The value is rendered as a standalone shell line, so it must be a simple
        environment module command rather than any safe-character shell command.

        Parameters
        ----------
        command : str
            Module command to validate.

        Raises
        ------
        ValueError
            If the command is not an allowed module command.
        """
        if "\n" in command or "\r" in command:
            raise ValueError("module_load must be a single-line module command")
        if "\0" in command:
            raise ValueError("module_load must not contain NUL bytes")

        stripped = command.strip()
        if not stripped:
            raise ValueError("module_load must not be empty")

        tokens = stripped.split()
        if any(not self._SETUP_COMMAND_TOKEN_PATTERN.match(token) for token in tokens):
            raise ValueError(
                "module_load may only contain alphanumerics, spaces, and -_./:,@%=+. "
                "Shell metacharacters are not allowed."
            )
        if tokens[0] != "module":
            raise ValueError("module_load must start with 'module'")

        self._validate_module_tokens(
            tokens,
            "module_load",
            require_swap_pair=True,
            reject_option_like_args=True,
        )

    def _validate_module_tokens(
        self,
        tokens: list[str],
        field_name: str,
        require_swap_pair: bool = False,
        reject_option_like_args: bool = False,
    ) -> None:
        """Validate tokenized environment module command structure.

        Parameters
        ----------
        tokens : list[str]
            Whitespace-split module command tokens.
        field_name : str
            Field name used in error messages.
        require_swap_pair : bool, optional
            Require ``module swap`` to include exactly old and new module names.
        reject_option_like_args : bool, optional
            Reject module arguments that look like command-line options.

        Raises
        ------
        ValueError
            If the module command does not match the supported action forms.
        """
        allowed_actions = (
            self._MODULE_ACTIONS_WITH_ARGS | self._MODULE_ACTIONS_WITHOUT_ARGS | {"swap"}
        )
        if len(tokens) < 2 or tokens[1] not in allowed_actions:
            raise ValueError(
                f"{field_name} module entries must use load, unload, use, unuse, purge, "
                "reset, restore, or swap"
            )

        action = tokens[1]
        if action in self._MODULE_ACTIONS_WITHOUT_ARGS:
            if len(tokens) != 2:
                raise ValueError(f"{field_name} module {action} entries must not include args")
            return

        if action == "swap":
            if require_swap_pair and len(tokens) != 4:
                raise ValueError(f"{field_name} module swap entries require two arguments")
            if not require_swap_pair and len(tokens) < 3:
                raise ValueError(f"{field_name} module swap entries require an argument")
            self._validate_module_args(tokens[2:], field_name, reject_option_like_args)
            return

        if len(tokens) < 3:
            raise ValueError(f"{field_name} module {action} entries require an argument")
        self._validate_module_args(tokens[2:], field_name, reject_option_like_args)

    def _validate_module_args(
        self, args: list[str], field_name: str, reject_option_like_args: bool
    ) -> None:
        """Validate module command arguments.

        Parameters
        ----------
        args : list[str]
            Module arguments following the action token.
        field_name : str
            Field name used in error messages.
        reject_option_like_args : bool
            Reject arguments beginning with ``-`` when strict validation is requested.

        Raises
        ------
        ValueError
            If strict module argument validation rejects an argument.
        """
        if reject_option_like_args and any(arg.startswith("-") for arg in args):
            raise ValueError(f"{field_name} module arguments must not start with '-'")

    def _validate_command_prefix(self, command_prefix: str) -> None:
        """Validate command prefix content for script interpolation.

        Parameters
        ----------
        command_prefix : str
            Prefix command to validate.

        Raises
        ------
        ValueError
            If the command prefix contains unsafe characters.
        """
        if not self._COMMAND_PREFIX_PATTERN.match(command_prefix):
            raise ValueError(
                "SLURM script field 'command_prefix' contains unsafe characters: "
                f"{command_prefix!r}."
            )

    def generate_job_script(
        self,
        config_path: str,
        replicate: int,
        working_dir: str,
        system_prefix: str,
        equilibration_mdps: list[str],
        job_name: str | None = None,
        output_file: str | None = None,
    ) -> str:
        """Generate a self-resubmitting GROMACS SLURM script.

        Parameters
        ----------
        config_path : str
            Path to the simulation YAML config.
        replicate : int
            Replicate index.
        working_dir : str
            Directory containing exported GROMACS files and outputs.
        system_prefix : str
            Prefix used by exported ``.gro`` and ``.top`` files.
        equilibration_mdps : list[str]
            Ordered equilibration stage MDP filenames.
        job_name : str | None, optional
            Custom scheduler job name.
        output_file : str | None, optional
            Custom SLURM log path pattern.

        Returns
        -------
        str
            Full script content.
        """
        if job_name is None:
            job_name = f"pzmd_gmx_r{replicate}"
        if output_file is None:
            output_file = f"slurm_logs/{job_name}.%j.out"

        wall_hours = self._parse_wall_time_hours(self._config.time_limit)
        maxh_hours = wall_hours * self.MAXH_SAFETY_FACTOR
        energy_minimization_block = self._generate_energy_minimization_block()
        equilibration_block = self._generate_equilibration_block(equilibration_mdps)
        manifest_path = _discover_manifest_path()
        module_load_line = self._module_load.strip() if self._module_load else ""
        gmx_base = self._gmx_binary
        if self._command_prefix:
            gmx_base = f"{self._command_prefix} {self._gmx_binary}"

        if self._is_mpi_binary and self._command_prefix and self._mpi_launcher_flags:
            LOGGER.warning(
                "command_prefix is set alongside mpi_launcher_flags for a real-MPI "
                "binary (%s). The command_prefix will handle process launching, so "
                "mpi_launcher_flags (%s) will be ignored. Remove mpi_launcher_flags "
                "or fold them into command_prefix.",
                self._gmx_binary,
                self._mpi_launcher_flags,
            )

        if self._is_mpi_binary and not self._command_prefix:
            mpi_prefix = "mpirun"
            if self._mpi_launcher_flags:
                mpi_prefix = f"mpirun {self._mpi_launcher_flags}"
            mdrun_command = f"{mpi_prefix} $GMX mdrun"
        else:
            mdrun_command = "$GMX mdrun"

        _validate_script_value(self._config.partition, "partition")
        _validate_script_value(job_name, "job_name")
        _validate_script_value(output_file, "output_file")
        _validate_script_value(self._config.time_limit, "time_limit")
        _validate_script_value(self._pixi_env, "pixi_env")
        _validate_script_value(str(manifest_path), "manifest_path")
        _validate_script_value(str(config_path), "config_path")
        _validate_script_value(str(working_dir), "working_dir")
        _validate_script_value(self._gmx_binary, "gmx_binary")
        if self._command_prefix:
            self._validate_command_prefix(self._command_prefix)
        if self._mpi_launcher_flags:
            _validate_script_value(self._mpi_launcher_flags, "mpi_launcher_flags")
        _validate_script_value(system_prefix, "system_prefix")
        _validate_script_value(self._grompp_flags, "grompp_flags")
        _validate_script_value(self._mdrun_flags, "mdrun_flags")
        if self._mdrun_flags_eq is not None:
            _validate_script_value(self._mdrun_flags_eq, "mdrun_flags_eq")
        if self._mdrun_flags_prod is not None:
            _validate_script_value(self._mdrun_flags_prod, "mdrun_flags_prod")
        if module_load_line:
            self._validate_module_load(module_load_line)
        for key, value in self._env_exports.items():
            if not self._ENV_VAR_NAME.match(key):
                raise ValueError(
                    f"env_exports key '{key}' is not a valid shell variable name. "
                    "Must match [A-Za-z_][A-Za-z0-9_]*."
                )
            _validate_script_value(value, f"env_exports[{key}]")
        for setup_command in self._setup_commands:
            self._validate_setup_command(setup_command)
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
        for mdp_name in equilibration_mdps:
            _validate_script_value(mdp_name, "equilibration_mdp")

        env_exports_block = "\n".join(
            f'export {key}="{value}"' for key, value in self._env_exports.items()
        )
        setup_commands_block = "\n".join(self._setup_commands)

        return render_package_template(
            _TEMPLATE_PACKAGE,
            _SLURM_TEMPLATE,
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
                "gmx_binary": gmx_base,
                "system_prefix": system_prefix,
                "maxh_hours": f"{maxh_hours:.2f}",
                "grompp_flags": self._grompp_flags,
                "mdrun_flags": self._mdrun_flags,
                "mdrun_flags_eq": self._mdrun_flags_eq or "",
                "mdrun_flags_prod": self._mdrun_flags_prod or "",
                "mdrun_command": mdrun_command,
                "energy_minimization_block": energy_minimization_block,
                "equilibration_block": equilibration_block,
                "module_load_line": module_load_line,
                "env_exports_block": env_exports_block,
                "setup_commands_block": setup_commands_block,
                "FULL_CREDIT_LINE": FULL_CREDIT_LINE,
            },
        )

    def _generate_energy_minimization_block(self) -> str:
        """Generate single-stage energy minimization bash block.

        Returns
        -------
        str
            Rendered bash block.
        """
        lines = [
            "if [ ! -f em.gro ]; then",
            '    echo "=== Energy Minimization ==="',
            "    rm -f em.tpr",
            "        run_foreground $GMX grompp -f em.mdp -c ${{PREFIX}}.gro -r ${{PREFIX}}.gro -p ${{PREFIX}}.top -o em.tpr {grompp_flags}",
            "    if [ -f em.cpt ]; then",
            '        echo "Resuming energy minimization from checkpoint: em.cpt"',
            '        run_mdrun_stage "energy minimization" $MDRUN -deffnm em -cpi em.cpt -cpo em.cpt -append $MDRUN_FLAGS_EM -v',
            "    else",
            '        run_mdrun_stage "energy minimization" $MDRUN -deffnm em -cpo em.cpt $MDRUN_FLAGS_EM -v',
            "    fi",
            "    if [ ! -f em.gro ]; then",
            '        echo "FATAL: Energy minimization failed — em.gro not produced"',
            "        exit 1",
            "    fi",
            "",
            "    # Verify EM health (detect infinite forces from bad initial geometry)",
            '    if grep -qi "force.*not finite\\|inf.*atom" em.log 2>/dev/null; then',
            '        echo "FATAL: Energy minimization failed — infinite forces detected in em.log"',
            '        echo "This usually indicates severe atomic overlaps in the initial structure."',
            '        echo "Try increasing packing padding or box size, or check Packmol logs."',
            "        exit 1",
            "    fi",
            '    echo "Energy minimization complete."',
            "else",
            '    echo "Skipping energy minimization (em.gro exists)."',
            "fi",
        ]
        return "\n".join(lines).format(grompp_flags=self._grompp_flags)

    def _parse_wall_time_hours(self, time_limit: str) -> float:
        """Parse SLURM wall time syntax and return hours.

        Parameters
        ----------
        time_limit : str
            SLURM time string in ``HH:MM:SS`` or ``D-HH:MM:SS`` format.

        Returns
        -------
        float
            Total wall time in hours.

        Raises
        ------
        ValueError
            If the input does not match a supported format.
        """
        days = 0
        hms_part = time_limit
        if "-" in time_limit:
            day_part, hms_part = time_limit.split("-", maxsplit=1)
            days = int(day_part)

        parts = hms_part.split(":")
        if len(parts) != 3:
            raise ValueError("Invalid SLURM time limit format. Expected HH:MM:SS or D-HH:MM:SS")

        hours, minutes, seconds = (int(part) for part in parts)
        total_hours = days * 24 + hours + (minutes / 60.0) + (seconds / 3600.0)
        return total_hours

    def _generate_equilibration_block(self, equilibration_mdps: list[str]) -> str:
        """Generate staged equilibration bash with skip-if-done behavior.

        Parameters
        ----------
        equilibration_mdps : list[str]
            Ordered list of equilibration MDP filenames.

        Returns
        -------
        str
            Rendered bash block.
        """
        lines: list[str] = ['LAST_EQ="em"']

        if not equilibration_mdps:
            lines.append(
                'echo "No equilibration stages configured; using em output for production."'
            )
            return "\n".join(lines)

        for idx, mdp_name in enumerate(equilibration_mdps, start=1):
            stage = f"eq_{idx:02d}"
            lines.append("")
            lines.append(f"if [ ! -f {stage}.gro ]; then")
            lines.append(f'    echo "=== Equilibration {idx}: {mdp_name} ==="')
            lines.append("        if [ -f ${LAST_EQ}.cpt ]; then")
            lines.append(f"            rm -f {stage}.tpr")
            lines.append(
                f"            run_foreground $GMX grompp -f {mdp_name} -c ${{LAST_EQ}}.gro -r em.gro "
                f"-t ${{LAST_EQ}}.cpt "
                f"-p ${{PREFIX}}.top -o {stage}.tpr {self._grompp_flags}"
            )
            lines.append("        else")
            lines.append(f"            rm -f {stage}.tpr")
            lines.append(
                f"            run_foreground $GMX grompp -f {mdp_name} -c ${{LAST_EQ}}.gro -r em.gro "
                f"-p ${{PREFIX}}.top -o {stage}.tpr {self._grompp_flags}"
            )
            lines.append("        fi")
            lines.append(f"    if [ -f {stage}.cpt ]; then")
            lines.append(
                f'        echo "Resuming equilibration stage {idx} from checkpoint: {stage}.cpt"'
            )
            lines.append(
                f'        run_mdrun_stage "equilibration stage {idx}" '
                f"$MDRUN -deffnm {stage} -cpi {stage}.cpt -cpo {stage}.cpt -append $MDRUN_FLAGS_EQ -v"
            )
            lines.append("    else")
            lines.append(
                f'        run_mdrun_stage "equilibration stage {idx}" '
                f"$MDRUN -deffnm {stage} -cpo {stage}.cpt $MDRUN_FLAGS_EQ -v"
            )
            lines.append("    fi")
            lines.append(f"    if [ ! -f {stage}.gro ]; then")
            lines.append(
                f'        echo "FATAL: Equilibration stage {idx} failed — {stage}.gro not produced"'
            )
            lines.append("        exit 1")
            lines.append("    fi")
            lines.append(f'    echo "Equilibration stage {idx} complete."')
            lines.append("else")
            lines.append(f'    echo "Skipping equilibration stage {idx} ({stage}.gro exists)."')
            lines.append("fi")
            lines.append(f'LAST_EQ="{stage}"')

        return "\n".join(lines)

    def save_script(self, script_content: str, output_path: Path) -> Path:
        """Save script content to disk and mark it executable.

        Parameters
        ----------
        script_content : str
            Rendered script text.
        output_path : Path
            Destination path.

        Returns
        -------
        Path
            Saved path.
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(script_content)
        os.chmod(output_path, 0o755)
        LOGGER.info(f"Saved GROMACS SLURM script to {output_path}")
        return output_path
