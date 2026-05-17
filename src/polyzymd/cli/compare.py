"""CLI commands for the compare command group.

This module lives under ``polyzymd.cli`` and provides the ``polyzymd compare``
command group with subcommands for initializing comparison projects and running
comparisons.
"""

from __future__ import annotations

import json
import logging
import shutil
import sys
from pathlib import Path
from typing import Any

import click
import yaml
from pydantic import ValidationError

from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.cli._compare_utils import (
    common_compare_options,
    load_comparison_config,
    validate_and_report,
)
from polyzymd.cli.colors import echo_logo
from polyzymd.config.comparison import (
    ComparisonConfig,
    generate_comparison_template,
)
from polyzymd.core.archived_features import (
    format_archived_analysis_message,
    get_archived_analysis_plugin,
)
from polyzymd.core.branding import prepend_file_header
from polyzymd.core.experimental import (
    echo_experimental_warning,
    experimental_features_for_comparison_type,
    format_experimental_suffix,
    normalize_experimental_features,
)

LOGGER = logging.getLogger("polyzymd.cli.compare")


def _display_path(path: Path) -> str:
    """Return a human-friendly display path, relative to CWD if possible."""
    try:
        return str(path.absolute().relative_to(Path.cwd().absolute()))
    except ValueError:
        return str(path.absolute())


def _resolve_hpc_dir(config: ComparisonConfig, analysis_name: str) -> Path:
    """Return the analysis HPC artifact directory."""
    source = config.source_path
    if source is None:
        return Path("comparison") / analysis_name / "_hpc"
    return source.parent / "comparison" / analysis_name / "_hpc"


def _contract_error_analysis_name(exc: PluginContractError, fallback: str) -> str:
    """Infer the analysis name from a plugin contract error message.

    Parameters
    ----------
    exc : PluginContractError
        Contract exception raised by the analysis framework.
    fallback : str
        Analysis name to use when the message cannot be parsed.

    Returns
    -------
    str
        Best-effort analysis name for user-facing diagnostics.
    """
    message = str(exc)
    if "." in message:
        candidate = message.split(".", 1)[0].strip()
        if candidate:
            return candidate
    return fallback


def _plugin_contract_click_exception(
    analysis_name: str,
    exc: PluginContractError,
) -> click.ClickException:
    """Build a user-facing ClickException for plugin contract failures.

    Parameters
    ----------
    analysis_name : str
        Name of the analysis plugin that violated the framework contract.
    exc : PluginContractError
        Original contract exception.

    Returns
    -------
    click.ClickException
        Click exception with an actionable diagnostic.
    """
    return click.ClickException(
        f"Plugin contract error in analysis '{analysis_name}': {exc}\n"
        "This is likely a PolyzyMD/plugin bug, not missing trajectory data. "
        "Please report it with the command and configuration used."
    )


def _archived_analysis_click_exception(analysis_name: str, *, context: str) -> click.ClickException:
    """Build a Click exception for archived analysis requests.

    Parameters
    ----------
    analysis_name : str
        User-provided analysis name.
    context : str
        CLI context where the analysis was requested.

    Returns
    -------
    click.ClickException
        User-facing exception with archive recovery information.
    """
    return click.ClickException(format_archived_analysis_message(analysis_name, context=context))


def _resolve_submit_resources_with_hints(
    *,
    plugin,
    pixi_path: str,
    partition: str | None,
    qos: str | None,
    account: str | None,
    ntasks: int,
    cpus_per_task: int,
    mem: str,
    time_limit: str,
    max_retries: int,
    mail_user: str | None,
    cpus_from_cli: bool,
    mem_from_cli: bool,
    time_from_cli: bool,
):
    """Resolve submit resources with CLI and plugin hint precedence.

    Precedence for CPU, memory, and time is:

    - explicit CLI flag
    - plugin ``slurm_resource_hint``
    - global system default
    """
    from polyzymd.workflow.analysis_slurm import AnalysisSlurmResources

    defaults = AnalysisSlurmResources()
    hint = getattr(type(plugin), "slurm_resource_hint", None)

    resolved_mem = mem
    if not mem_from_cli:
        resolved_mem = hint.mem if hint is not None and hint.mem is not None else defaults.mem

    resolved_time = time_limit
    if not time_from_cli:
        resolved_time = hint.time if hint is not None and hint.time is not None else defaults.time

    resolved_cpus = cpus_per_task
    if not cpus_from_cli:
        resolved_cpus = (
            hint.cpus_per_task
            if hint is not None and hint.cpus_per_task is not None
            else defaults.cpus_per_task
        )

    return AnalysisSlurmResources(
        pixi_path=pixi_path,
        partition=partition,
        qos=qos,
        account=account,
        ntasks=ntasks,
        cpus_per_task=resolved_cpus,
        mem=resolved_mem,
        time=resolved_time,
        max_retries=max_retries,
        mail_user=mail_user,
    )


def _echo_submit_all_summary(rows: list[dict[str, str]]) -> None:
    """Print a compact submit-all summary table."""
    click.echo()
    click.echo("Submitted analyses summary:")
    click.echo("analysis\tmode\tjobs\tfinalize_job_id")
    for row in rows:
        click.echo(f"{row['analysis']}\t{row['mode']}\t{row['jobs']}\t{row['finalize_job_id']}")


def _echo_qos_tip_if_needed(partition: str | None, qos: str | None) -> None:
    """Print a gentle QoS reminder for common HPC queues."""
    if partition and not qos:
        click.echo(
            "Tip: Many HPC clusters require --qos to be set. If jobs fail to start,\n"
            "try: polyzymd compare submit <analysis> --qos <your-qos>"
        )


def _resolve_manifest_task_condition(
    manifest,
    condition_index: int,
):
    """Resolve and validate manifest condition index.

    Parameters
    ----------
    manifest : AnalysisJobManifest
        Loaded job manifest.
    condition_index : int
        Requested condition index.

    Returns
    -------
    ConditionTaskSpec
        Condition spec from the frozen manifest.

    Raises
    ------
    click.ClickException
        If condition_index is not in manifest bounds.
    """
    if condition_index < 0 or condition_index >= len(manifest.condition_specs):
        raise click.ClickException(
            f"Condition index {condition_index} is out of range for "
            f"manifest with {len(manifest.condition_specs)} conditions"
        )
    return manifest.condition_specs[condition_index]


def _settings_from_manifest(plugin, manifest):
    """Build plugin settings from the frozen manifest snapshot.

    Parameters
    ----------
    plugin : Analysis
        Analysis plugin instance.
    manifest : AnalysisJobManifest
        Loaded job manifest.

    Returns
    -------
    BaseModel
        Parsed plugin settings from manifest snapshot.

    Raises
    ------
    click.ClickException
        If snapshot settings do not validate for the plugin.
    """
    try:
        return plugin.Settings.model_validate(manifest.settings_snapshot)
    except (ValueError, KeyError) as exc:
        raise click.ClickException(
            f"Invalid settings_snapshot in manifest for analysis '{manifest.analysis_name}': {exc}"
        ) from exc


def _echo_branding() -> None:
    """Print the PolyzyMD ASCII logo for top-level comparison commands."""
    echo_logo()


@click.group()
def compare():
    """Compare analysis results across simulation conditions.

    The compare module allows you to statistically compare RMSF and other
    metrics across multiple simulation conditions (e.g., different polymer
    compositions, temperatures, etc.).

    \b
    Local workflow (interactive):
    1. polyzymd compare init -n <name>   # Create project with template
    2. Edit comparison.yaml              # Add your conditions
    3. polyzymd compare run <analysis>   # Run a single comparison
    4. polyzymd compare run-all          # Run all enabled comparisons

    \b
    HPC workflow (SLURM):
    1. polyzymd compare init -n <name>   # Create project with template
    2. Edit comparison.yaml              # Add your conditions
    3. polyzymd compare submit <analysis> --partition <part>  # Submit DAG
    4. polyzymd compare status <analysis>   # Monitor progress
    5. polyzymd compare finalize <analysis> # (if needed) re-run compare+plot

    \b
    Example (local):
        polyzymd compare run rmsf --eq-time 10ns

    \b
    Example (HPC):
        polyzymd compare submit sasa --partition aa100 --mem 8G --time 02:00:00
    """
    pass


@compare.command()
@click.option(
    "-n",
    "--name",
    required=True,
    help="Name for the comparison project (will create a directory with this name).",
)
@click.option(
    "--eq-time",
    default="10ns",
    help="Default equilibration time for analysis (e.g., '10ns', '5000ps').",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="Parent directory for the comparison project. Defaults to current directory.",
)
def init(name: str, eq_time: str, output_dir: Path | None):
    """Initialize a new comparison project.

    Creates a new directory NAME containing:
      - comparison.yaml: Template configuration file to edit
      - structures/: Directory for shared structure files (enzyme PDB)
      - comparison/: Directory for comparison result JSON files
      - figures/: Directory for comparison plots

    \b
    Example:
        polyzymd compare init -n polymer_stability_study
        polyzymd compare init -n my_study --eq-time 20ns
        polyzymd compare init -n my_study -o /path/to/projects
    """
    # Determine output location
    if output_dir is None:
        output_dir = Path.cwd()
    else:
        output_dir = Path(output_dir).resolve()

    project_dir = output_dir / name

    # Check if already exists
    if project_dir.exists():
        click.echo(f"Error: Directory already exists: {project_dir}", err=True)
        sys.exit(1)

    try:
        _echo_branding()

        # Create directory structure
        project_dir.mkdir(parents=True)
        (project_dir / "comparison").mkdir()
        (project_dir / "figures").mkdir()
        (project_dir / "structures").mkdir()

        # Create README in structures directory
        structures_readme = """\
# Structures Directory

Place shared structure files here for use in comparison analyses.

Keep condition-specific input structures with their simulation projects and
place only comparison-wide references in this directory.
"""
        (project_dir / "structures" / "README.md").write_text(
            prepend_file_header(structures_readme, comment_prefix="#")
        )

        # Generate and write template
        template_content = generate_comparison_template(name, eq_time)
        config_path = project_dir / "comparison.yaml"
        config_path.write_text(template_content)

        click.echo(f"Created comparison project: {project_dir}")
        click.echo()
        click.echo("Next steps:")
        click.echo(f"  1. Edit {_display_path(config_path)}")
        click.echo("     - Add your simulation conditions (paths to config.yaml files)")
        click.echo("     - Define catalytic_triad for active site analysis")
        click.echo()
        click.echo(f"  2. cd {_display_path(project_dir)}")
        click.echo("  3. Run comparisons:")
        click.echo("     polyzymd compare run rmsf      # Compare flexibility")
        click.echo("     polyzymd compare run triad     # Compare triad geometry")
        click.echo("     polyzymd compare run contacts  # Compare polymer-protein contacts")
        click.echo()
        click.echo("  On an HPC cluster, submit as SLURM jobs instead:")
        click.echo("     polyzymd compare submit sasa --partition <part> --mem 8G")
        click.echo()

    except (OSError, IOError) as e:
        click.echo(f"Error creating project: {e}", err=True)
        sys.exit(1)


@compare.command()
@click.option(
    "-f",
    "--file",
    "config_file",
    type=click.Path(path_type=Path),
    default="comparison.yaml",
    help="Path to comparison.yaml config file.",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["table", "json"]),
    default="table",
    help="Output format: table (default) or json.",
)
def validate(config_file: Path, output_format: str):
    """Validate a comparison.yaml configuration file.

    Checks the configuration for errors without running any analyses.
    Useful for catching configuration problems before running expensive
    computations.

    \b
    Validates:
      - YAML syntax and structure
      - Required fields present
      - At least 1 condition defined
      - Condition labels are unique
      - Control label matches a condition (if specified)
      - Config files exist for each condition

    \b
    Example:
        polyzymd compare validate
        polyzymd compare validate -f my_comparison.yaml
        polyzymd compare validate --format json
    """
    config_file = Path(config_file).resolve()

    # Prepare result structure
    result = {
        "file": str(config_file),
        "valid": False,
        "errors": [],
        "summary": {},
    }

    # Check file exists
    if not config_file.exists():
        result["errors"].append(f"File not found: {config_file}")
        _output_validation_result(result, output_format)
        sys.exit(1)

    # Try to load and validate
    try:
        config = ComparisonConfig.from_yaml(config_file)

        # Run validation
        errors = config.validate_config()
        result["errors"] = errors
        result["valid"] = len(errors) == 0

        # Build summary
        result["summary"] = {
            "name": config.name,
            "conditions_count": len(config.conditions),
            "condition_labels": [c.label for c in config.conditions],
            "control": config.control,
            "sections_configured": config.plugins.get_enabled_plugins(),
        }

    except yaml.YAMLError as e:
        result["errors"].append(f"YAML syntax error: {e}")
    except (ValidationError, ValueError, OSError) as e:
        result["errors"].append(f"Validation error: {e}")

    _output_validation_result(result, output_format)

    if not result["valid"]:
        sys.exit(1)


def _output_validation_result(result: dict, output_format: str) -> None:
    """Output validation result in the specified format.

    Parameters
    ----------
    result : dict
        Validation result dictionary
    output_format : str
        Output format: 'table' or 'json'
    """
    import json as json_module

    if output_format == "json":
        click.echo(json_module.dumps(result, indent=2))
        return

    # Human-readable table format
    click.echo(f"Validating: {result['file']}")
    click.echo()

    if result["valid"]:
        click.secho("✓ Configuration is valid", fg="green")
        click.echo()

        summary = result.get("summary", {})
        if summary:
            click.echo(f"  Name: {summary.get('name', 'N/A')}")
            click.echo(f"  Conditions: {summary.get('conditions_count', 0)}")
            labels = summary.get("condition_labels", [])
            if labels:
                click.echo(f"    - {', '.join(labels)}")
            control = summary.get("control")
            if control:
                click.echo(f"  Control: {control}")
            sections = summary.get("sections_configured", [])
            if sections:
                click.echo(f"  Analysis sections: {', '.join(sections)}")
    else:
        click.secho("✗ Configuration has errors", fg="red")
        click.echo()
        for error in result["errors"]:
            click.echo(f"  • {error}", err=True)


@compare.command("run")
@click.argument(
    "comparison_type",
    type=str,
    required=False,
    default=None,
)
@common_compare_options
@click.option(
    "--recompute",
    is_flag=True,
    help="Force recompute even if cached results exist.",
)
@click.option(
    "--list",
    "list_types",
    is_flag=True,
    help="List available comparison types and exit.",
)
def run_comparison(
    comparison_type: str | None,
    config_file: Path,
    eq_time: str | None,
    recompute: bool,
    output_format: str,
    output_path: Path | None,
    quiet: bool,
    debug: bool,
    list_types: bool,
):
    """Run a comparison using the analysis plugin system.

    Runs any discovered analysis plugin by name or alias. Use --list to see
    all available comparison types.

    This command runs analysis interactively in the current process. For
    large studies on HPC clusters, use 'polyzymd compare submit' to dispatch
    jobs to SLURM instead.

    \b
    Examples:
        polyzymd compare run rmsf
        polyzymd compare run triad --eq-time 10ns
        polyzymd compare run contacts --format markdown
        polyzymd compare run --list
    """
    from polyzymd.analyses.discovery import get_analysis, list_all_names, list_analyses
    from polyzymd.analyses.orchestrator import run_comparison as _run_pipeline
    from polyzymd.cli.logging_utils import setup_logging

    # Handle --list flag
    if list_types:
        analyses = list_analyses()
        click.echo("Available comparison types:")
        for name, cls in analyses.items():
            aliases = ", ".join(cls.aliases) if cls.aliases else ""
            suffix = f" (aliases: {aliases})" if aliases else ""
            click.echo(f"  - {name}{suffix}: {cls.__name__}")
        return

    # Require comparison_type if not listing
    if comparison_type is None:
        click.echo("Error: Missing argument 'COMPARISON_TYPE'.", err=True)
        click.echo("Use 'polyzymd compare run --list' to see available types.", err=True)
        sys.exit(1)

    # Set up logging
    setup_logging(quiet=quiet, debug=debug)

    # Look up the analysis plugin
    try:
        analysis_cls = get_analysis(comparison_type)
    except KeyError:
        if get_archived_analysis_plugin(comparison_type) is not None:
            raise _archived_analysis_click_exception(
                comparison_type,
                context="CLI comparison lookup",
            ) from None

        available = list_all_names()
        click.echo(f"Error: Unknown comparison type '{comparison_type}'", err=True)
        click.echo(f"Available types: {', '.join(available)}", err=True)
        click.echo("", err=True)
        click.echo("Use 'polyzymd compare run --list' to see all available types.", err=True)
        sys.exit(1)

    config = load_comparison_config(config_file)
    validate_and_report(config)

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Type: {analysis_cls.name}")
    click.echo(f"Conditions: {len(config.conditions)}")

    equilibration = eq_time or config.defaults.equilibration_time
    click.echo(f"Equilibration: {equilibration}")
    click.echo()

    # Run the full pipeline (compute replicates -> aggregate -> compare -> plot)
    analysis = analysis_cls()
    try:
        pipeline_result = _run_pipeline(
            analysis, config, recompute=recompute, equilibration=equilibration
        )
        result = pipeline_result["comparison"]
    except PluginContractError as e:
        raise _plugin_contract_click_exception(analysis_cls.name, e) from e
    except (FileNotFoundError, ValueError, OSError) as e:
        raise click.ClickException(f"Error during comparison: {e}") from e
    except Exception as e:
        # Unexpected error — show traceback hint
        click.echo(f"Unexpected error during comparison: {e}", err=True)
        if debug:
            import traceback

            traceback.print_exc()
        else:
            click.echo("Re-run with --debug for full traceback.", err=True)
        sys.exit(1)

    if result is None:
        click.echo("Warning: comparison returned no result (not enough data?).", err=True)
        sys.exit(1)

    json_path = pipeline_result.get("comparison_path")
    if json_path is not None:
        click.echo(f"Saved result: {json_path}")
    click.echo()

    # Format and display
    if output_format == "json":
        # JSON output: serialize directly, skip plugin formatter
        if hasattr(result, "model_dump_json"):
            formatted = result.model_dump_json(indent=2)
        else:
            import json as json_module

            formatted = json_module.dumps(
                result if isinstance(result, dict) else str(result), indent=2
            )
    else:
        # Table/markdown: delegate to the plugin's format() method
        fmt = "text" if output_format == "table" else output_format
        formatted = analysis.format(result, output_format=fmt)

    click.echo(formatted)

    # Save formatted output if requested
    if output_path:
        try:
            output_path = Path(output_path)
            output_path.write_text(formatted)
            click.echo(f"Saved output: {output_path}")
        except OSError as e:
            raise click.ClickException(f"Could not write output file: {e}") from e


@compare.command("plot-all")
@click.option(
    "-f",
    "--file",
    "config_file",
    type=click.Path(exists=True, path_type=Path),
    default="comparison.yaml",
    help="Path to comparison.yaml config file.",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Override output directory (default: from plot_settings in config).",
)
@click.option(
    "-a",
    "--analysis",
    "analysis_type",
    type=str,
    default=None,
    help="Generate plots for specific analysis type only (e.g., 'rmsf', 'catalytic_triad').",
)
@click.option(
    "--list-available",
    is_flag=True,
    help="List enabled analysis types that can produce plots.",
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help="Suppress INFO messages.",
)
@click.option(
    "--debug",
    is_flag=True,
    help="Enable DEBUG logging.",
)
def plot_all(
    config_file: Path,
    output_dir: Path | None,
    analysis_type: str | None,
    list_available: bool,
    quiet: bool,
    debug: bool,
):
    """Generate all plots from comparison.yaml configuration.

    This is the config-driven plotting command that reads plot_settings
    from comparison.yaml and generates all configured plots automatically.
    Each enabled analysis plugin's ``plot()`` method is called to produce
    its figures.

    \b
    Examples:
        polyzymd compare plot-all -f comparison.yaml
        polyzymd compare plot-all -f comparison.yaml -a catalytic_triad
        polyzymd compare plot-all --list-available
    """
    _echo_branding()

    # Configure logging
    log_level = logging.WARNING if quiet else (logging.DEBUG if debug else logging.INFO)
    logging.basicConfig(
        level=log_level,
        format="%(levelname)s: %(message)s",
    )

    # Load config
    try:
        config = ComparisonConfig.from_yaml(config_file)
    except (FileNotFoundError, yaml.YAMLError, ValidationError, ValueError) as e:
        raise click.ClickException(f"Error loading config: {e}") from e

    # Override output directory if specified
    if output_dir:
        config.plot_settings.output_dir = output_dir

    enabled = config.plugins.get_enabled_plugins()

    # List available analyses
    if list_available:
        click.echo("Enabled analysis types:")
        for atype in enabled:
            comparison_type = atype
            analysis_settings = config.plugins.get(atype)
            suffix = format_experimental_suffix(
                experimental_features_for_comparison_type(comparison_type, analysis_settings)
            )
            click.echo(f"  - {atype}{suffix}")
        return

    # Determine which analyses to plot
    if analysis_type:
        # Validate the requested analysis type exists
        from polyzymd.analyses.discovery import get_analysis

        try:
            get_analysis(analysis_type)
        except KeyError:
            if get_archived_analysis_plugin(analysis_type) is not None:
                raise _archived_analysis_click_exception(
                    analysis_type,
                    context="CLI plot lookup",
                ) from None

            from polyzymd.analyses.discovery import list_all_names

            available = list_all_names()
            raise click.BadParameter(
                f"Unknown analysis type '{analysis_type}'. "
                f"Available: {', '.join(sorted(available))}",
                param_hint="'-a' / '--analysis'",
            )
        target_analyses = [analysis_type]
    else:
        target_analyses = list(enabled)

    # Experimental warnings
    experimental_features: tuple[str, ...]
    if analysis_type:
        comparison_type = analysis_type
        analysis_settings = config.plugins.get(analysis_type)
        experimental_features = experimental_features_for_comparison_type(
            comparison_type, analysis_settings
        )
    else:
        feature_list: list[str] = []
        for settings_key in enabled:
            comparison_type = settings_key
            analysis_settings = config.plugins.get(settings_key)
            feature_list.extend(
                experimental_features_for_comparison_type(comparison_type, analysis_settings)
            )
        experimental_features = normalize_experimental_features(feature_list)

    echo_experimental_warning(experimental_features)

    # Resolve output dir for display
    plot_settings = config.plot_settings
    display_dir = plot_settings.output_dir
    if not display_dir.is_absolute() and config.source_path is not None:
        display_dir = (config.source_path.parent / display_dir).resolve()

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Conditions: {len(config.conditions)}")
    click.echo(f"Output directory: {display_dir}")
    click.echo()

    # Generate plots via orchestrator
    from polyzymd.analyses.orchestrator import run_all_plots

    click.echo(f"Generating plots for: {', '.join(target_analyses)}...")
    generated, failures = run_all_plots(config, target_analyses)

    click.echo()
    if generated:
        click.echo(f"Generated {len(generated)} plots:")
        for path in generated:
            click.echo(f"  - {path}")
    else:
        click.echo("No plots generated. Check that analyses are enabled in config.")

    if failures:
        click.echo()
        click.echo(f"Failed ({len(failures)}):", err=True)
        for name, msg in failures:
            click.echo(f"  - {name}: {msg}", err=True)
        sys.exit(1)


# ---------------------------------------------------------------------------
# run-all: batch runner for every comparison defined in comparison.yaml
# ---------------------------------------------------------------------------


@compare.command("run-all")
@click.option(
    "-f",
    "--file",
    "config_file",
    type=click.Path(path_type=Path),
    default="comparison.yaml",
    help="Path to comparison.yaml config file.",
)
@click.option(
    "--eq-time",
    default=None,
    help="Override equilibration time (e.g., '10ns', '5000ps').",
)
@click.option(
    "--recompute",
    is_flag=True,
    help="Force recompute even if cached results exist.",
)
@click.option(
    "--plot/--no-plot",
    "run_plots",
    default=False,
    help="Run plot-all after all comparisons complete.",
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help="Suppress INFO messages, show warnings/errors only.",
)
@click.option(
    "--debug",
    is_flag=True,
    help="Enable DEBUG logging for troubleshooting.",
)
def run_all(
    config_file: Path,
    eq_time: str | None,
    recompute: bool,
    run_plots: bool,
    quiet: bool,
    debug: bool,
):
    """Run ALL comparisons defined in comparison.yaml.

    Iterates over every enabled analysis in the config and runs the
    corresponding analysis plugin. Results are saved under
    ``comparison/<analysis>/result.json`` next to the config file.

    \b
    Workflow:
        1. Edit comparison.yaml (enable the analyses you want)
        2. polyzymd compare run-all -f comparison.yaml
        3. (optional) polyzymd compare run-all --plot  # also generate plots

    \b
    Examples:
        polyzymd compare run-all
        polyzymd compare run-all -f comparison.yaml --eq-time 10ns
        polyzymd compare run-all --recompute --plot
    """
    from polyzymd.analyses.orchestrator import run_all_comparisons
    from polyzymd.cli.logging_utils import setup_logging

    _echo_branding()

    setup_logging(quiet=quiet, debug=debug)

    # --- Load and validate config -------------------------------------------
    config = load_comparison_config(config_file)
    validate_and_report(config)

    equilibration = eq_time or config.defaults.equilibration_time
    # --- Discover enabled analyses ------------------------------------------
    enabled_analyses = config.plugins.get_enabled_plugins()
    if not enabled_analyses:
        click.echo("No analyses are enabled in comparison.yaml.", err=True)
        sys.exit(1)

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Conditions: {len(config.conditions)}")
    click.echo(f"Equilibration: {equilibration}")
    click.echo(f"Enabled analyses: {', '.join(enabled_analyses)}")
    click.echo()

    # --- Run all analyses in dependency order --------------------------------
    try:
        all_results = run_all_comparisons(
            config,
            analysis_names=None,  # run all enabled
            recompute=recompute,
            equilibration=equilibration,
        )
    except PluginContractError as e:
        analysis_name = _contract_error_analysis_name(e, "run-all")
        raise _plugin_contract_click_exception(analysis_name, e) from e

    # Categorize results for summary
    succeeded: list[str] = []
    failed: list[tuple[str, str]] = []
    skipped: list[str] = []

    for name, result in all_results.items():
        if "error" in result:
            click.echo(f"  [{name}] FAILED: {result['error']}", err=True)
            failed.append((name, result["error"]))
        elif result.get("comparison") is None:
            click.echo(f"  [{name}] skipped — no comparison result")
            skipped.append(name)
        else:
            json_path = result.get("comparison_path")
            click.echo(f"  [{name}] saved -> {json_path}")
            succeeded.append(name)

    # --- Summary ------------------------------------------------------------
    click.echo()
    click.echo("=" * 60)
    click.echo(f"  Succeeded: {len(succeeded)}")
    if succeeded:
        click.echo(f"    {', '.join(succeeded)}")
    if skipped:
        click.echo(f"  Skipped:   {len(skipped)}")
        click.echo(f"    {', '.join(skipped)}")
    if failed:
        click.echo(f"  Failed:    {len(failed)}")
        for name, msg in failed:
            click.echo(f"    {name}: {msg}")
    click.echo("=" * 60)

    # --- Optional plotting --------------------------------------------------
    if run_plots and succeeded:
        from polyzymd.analyses.orchestrator import run_all_plots

        click.echo()
        click.echo("Generating plots ...")

        generated, plot_failures = run_all_plots(
            config,
            succeeded,
            equilibration=equilibration,
        )
        click.echo(f"Generated {len(generated)} plots.")
        if plot_failures:
            click.echo(f"Plot failures ({len(plot_failures)}):", err=True)
            for name, msg in plot_failures:
                click.echo(f"  - {name}: {msg}", err=True)
            # Plot failures also trigger non-zero exit
            if not failed:
                failed.append(("plot", "some plots failed"))

    # Exit with error code if any comparison failed
    if failed:
        sys.exit(1)


@compare.command("submit")
@click.pass_context
@click.argument("analysis", type=str)
@click.option(
    "-f",
    "--file",
    "config_file",
    type=click.Path(path_type=Path),
    default="comparison.yaml",
    help="Path to comparison.yaml config file.",
)
@click.option("--partition", default=None, help="SLURM partition.")
@click.option("--qos", default=None, help="SLURM QoS.")
@click.option("--account", default=None, help="SLURM account/allocation.")
@click.option("--pixi-path", default="pixi", show_default=True, help="Path to pixi executable.")
@click.option("--ntasks", default=1, type=int, show_default=True)
@click.option("--cpus-per-task", default=1, type=int, show_default=True)
@click.option("--mem", default="4G", show_default=True, help="SLURM memory request.")
@click.option("--time", "time_limit", default="01:00:00", show_default=True, help="SLURM walltime.")
@click.option("--max-retries", default=3, type=int, show_default=True)
@click.option("--mail-user", default=None, help="Email for failure notifications.")
@click.option("--recompute", is_flag=True, help="Force recomputation in workers.")
@click.option(
    "--allow-partial",
    is_flag=True,
    help="Allow finalize to proceed when some conditions are missing aggregated results.",
)
@click.option("--equilibration", default=None, help="Override equilibration time.")
@click.option("--dry-run", is_flag=True, help="Generate scripts without submitting jobs.")
@click.option(
    "--job-arrays",
    is_flag=True,
    help="Submit one SLURM array job per condition for replicate workers.",
)
def submit_analysis_hpc(
    ctx: click.Context,
    analysis: str,
    config_file: Path,
    partition: str | None,
    qos: str | None,
    account: str | None,
    pixi_path: str,
    ntasks: int,
    cpus_per_task: int,
    mem: str,
    time_limit: str,
    max_retries: int,
    mail_user: str | None,
    recompute: bool,
    allow_partial: bool,
    equilibration: str | None,
    dry_run: bool,
    job_arrays: bool,
):
    """Submit replicate-level SLURM analysis DAG for one plugin."""
    from polyzymd.analyses.discovery import get_analysis
    from polyzymd.workflow.analysis_slurm import (
        build_manifest,
        generate_aggregate_script,
        generate_array_script,
        generate_finalize_script,
        generate_replicate_script,
        submit_analysis_graph,
        submit_analysis_graph_with_arrays,
    )

    # Suppress noisy third-party INFO logs during submission planning
    for logger_name in ("MDAnalysis", "numexpr"):
        logging.getLogger(logger_name).setLevel(logging.WARNING)

    _echo_qos_tip_if_needed(partition, qos)

    if get_archived_analysis_plugin(analysis) is not None:
        raise _archived_analysis_click_exception(
            analysis,
            context="CLI submission lookup",
        )

    if not dry_run and shutil.which("sbatch") is None:
        raise click.ClickException(
            "SLURM is not available: 'sbatch' not found on PATH. The HPC submission "
            "commands require a SLURM cluster. Run analysis locally with "
            "'polyzymd compare run' instead."
        )

    try:
        config = ComparisonConfig.from_yaml(config_file)
    except (FileNotFoundError, yaml.YAMLError, ValidationError, ValueError) as e:
        raise click.ClickException(f"Error loading config: {e}") from e

    try:
        analysis_cls = get_analysis(analysis)
    except KeyError:
        from polyzymd.analyses.discovery import list_all_names

        available = list_all_names()
        raise click.ClickException(
            f"Unknown analysis type '{analysis}'. Available: {', '.join(sorted(available))}"
        )
    plugin = analysis_cls()

    dependencies = tuple(getattr(analysis_cls, "dependencies", ()))
    if dependencies:
        source = config.source_path
        comparison_root = source.parent / "comparison" if source is not None else Path("comparison")
        for dep_name in dependencies:
            expected_path = comparison_root / dep_name / "result.json"
            if not expected_path.exists():
                raise click.UsageError(
                    (
                        f"Error: '{plugin.name}' depends on '{dep_name}', but '{dep_name}' comparison "
                        "results were not found at: "
                        f"{expected_path}\n\n"
                        f"Run '{dep_name}' first:\n"
                        f"  polyzymd compare submit {dep_name} [--partition ...] [--qos ...]\n\n"
                        "Or use 'compare submit-all' to submit all analyses with correct "
                        "dependency ordering."
                    )
                )

    resources = _resolve_submit_resources_with_hints(
        plugin=plugin,
        pixi_path=pixi_path,
        partition=partition,
        qos=qos,
        account=account,
        ntasks=ntasks,
        cpus_per_task=cpus_per_task,
        mem=mem,
        time_limit=time_limit,
        max_retries=max_retries,
        mail_user=mail_user,
        cpus_from_cli=(
            ctx.get_parameter_source("cpus_per_task") != click.core.ParameterSource.DEFAULT
        ),
        mem_from_cli=ctx.get_parameter_source("mem") != click.core.ParameterSource.DEFAULT,
        time_from_cli=ctx.get_parameter_source("time_limit") != click.core.ParameterSource.DEFAULT,
    )
    manifest = build_manifest(
        plugin,
        config,
        resources,
        recompute,
        equilibration,
        allow_partial=allow_partial,
    )
    hpc_dir = _resolve_hpc_dir(config, plugin.name)
    hpc_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = hpc_dir / "manifest.json"
    manifest.save(manifest_path)

    replicate_count = sum(len(cond.replicate_specs) for cond in manifest.condition_specs)
    array_count = len(manifest.condition_specs)
    aggregate_count = len(manifest.condition_specs)
    finalize_only = getattr(manifest, "pipeline_mode", "full") == "finalize_only"

    if dry_run:
        if not finalize_only:
            for cond in manifest.condition_specs:
                if job_arrays:
                    generate_array_script(
                        cond,
                        manifest,
                        resources,
                        [rep.replicate for rep in cond.replicate_specs],
                        hpc_dir,
                    )
                else:
                    for rep in cond.replicate_specs:
                        generate_replicate_script(manifest, rep, resources, hpc_dir)
                generate_aggregate_script(manifest, cond, resources, hpc_dir)
        generate_finalize_script(manifest, resources, hpc_dir)
        if finalize_only:
            click.echo("Would submit 1 finalize job (compare-only plugin)")
        elif job_arrays:
            total = array_count + aggregate_count + 1
            click.echo(
                "Would submit "
                f"{array_count} array jobs + {aggregate_count} aggregate + 1 finalize = {total} total"
            )
            click.echo("Mode: job arrays")
        else:
            total = replicate_count + aggregate_count + 1
            click.echo(
                "Would submit "
                f"{total} jobs ({replicate_count} replicate + {aggregate_count} aggregate + 1 finalize)"
            )
        click.echo("Dry run only: no jobs were submitted")
        return

    if finalize_only:
        graph = submit_analysis_graph(manifest, resources, hpc_dir)
    elif job_arrays:
        graph = submit_analysis_graph_with_arrays(manifest, resources, hpc_dir)
    else:
        graph = submit_analysis_graph(manifest, resources, hpc_dir)
    graph.save(hpc_dir / "job_graph.json")
    if finalize_only:
        click.echo("Submitted 1 finalize job (compare-only plugin)")
    elif job_arrays:
        total = array_count + aggregate_count + 1
        click.echo(
            "Submitted "
            f"{array_count} array jobs + {aggregate_count} aggregate + 1 finalize = {total} total"
        )
    else:
        total = replicate_count + aggregate_count + 1
        click.echo(
            "Submitted "
            f"{total} jobs ({replicate_count} replicate + {aggregate_count} aggregate + 1 finalize)"
        )


@compare.command("submit-all")
@click.pass_context
@click.option(
    "-f",
    "--file",
    "config_file",
    type=click.Path(path_type=Path),
    default="comparison.yaml",
    help="Path to comparison.yaml config file.",
)
@click.option("--partition", default=None, help="SLURM partition.")
@click.option("--qos", default=None, help="SLURM QoS.")
@click.option("--account", default=None, help="SLURM account/allocation.")
@click.option("--pixi-path", default="pixi", show_default=True, help="Path to pixi executable.")
@click.option("--ntasks", default=1, type=int, show_default=True)
@click.option("--cpus-per-task", default=1, type=int, show_default=True)
@click.option("--mem", default="4G", show_default=True, help="SLURM memory request.")
@click.option("--time", "time_limit", default="01:00:00", show_default=True, help="SLURM walltime.")
@click.option("--max-retries", default=3, type=int, show_default=True)
@click.option("--mail-user", default=None, help="Email for failure notifications.")
@click.option("--recompute", is_flag=True, help="Force recomputation in workers.")
@click.option(
    "--allow-partial",
    is_flag=True,
    help="Allow finalize to proceed when some conditions are missing aggregated results.",
)
@click.option("--equilibration", default=None, help="Override equilibration time.")
@click.option("--dry-run", is_flag=True, help="Generate scripts without submitting jobs.")
@click.option(
    "--exclude",
    "excluded_analyses",
    multiple=True,
    help="Exclude one analysis plugin name from submission.",
)
def submit_all_analyses_hpc(
    ctx: click.Context,
    config_file: Path,
    partition: str | None,
    qos: str | None,
    account: str | None,
    pixi_path: str,
    ntasks: int,
    cpus_per_task: int,
    mem: str,
    time_limit: str,
    max_retries: int,
    mail_user: str | None,
    recompute: bool,
    allow_partial: bool,
    equilibration: str | None,
    dry_run: bool,
    excluded_analyses: tuple[str, ...],
):
    """Submit all enabled analyses with dependency-ordered SLURM DAGs."""
    from polyzymd.analyses.discovery import get_analysis
    from polyzymd.analyses.orchestrator import order_analyses_for_execution
    from polyzymd.workflow.analysis_slurm import (
        build_manifest,
        generate_aggregate_script,
        generate_finalize_script,
        generate_replicate_script,
        submit_analysis_graph,
    )

    # Suppress noisy third-party INFO logs during submission planning
    for logger_name in ("MDAnalysis", "numexpr"):
        logging.getLogger(logger_name).setLevel(logging.WARNING)

    _echo_qos_tip_if_needed(partition, qos)

    try:
        config = ComparisonConfig.from_yaml(config_file)
    except (FileNotFoundError, yaml.YAMLError, ValidationError, ValueError) as e:
        raise click.ClickException(f"Error loading config: {e}") from e

    if not dry_run and shutil.which("sbatch") is None:
        raise click.ClickException(
            "SLURM is not available: 'sbatch' not found on PATH. The HPC submission "
            "commands require a SLURM cluster. Run analysis locally with "
            "'polyzymd compare run' instead."
        )

    enabled = list(config.plugins.get_enabled_plugins())
    excluded_set = set(excluded_analyses)
    filtered = [name for name in enabled if name not in excluded_set]
    if not filtered:
        raise click.ClickException("No enabled analyses remain after applying --exclude filters.")

    source = config.source_path
    comparison_root = source.parent / "comparison" if source is not None else Path("comparison")
    satisfied: set[str] = set()
    for excluded_name in excluded_set:
        result_path = comparison_root / excluded_name / "result.json"
        if result_path.exists():
            satisfied.add(excluded_name)

    ordered = order_analyses_for_execution(filtered, satisfied=satisfied)

    finalize_ids: dict[str, str] = {}
    summary_rows: list[dict[str, str]] = []
    for analysis_name in ordered:
        analysis_cls = get_analysis(analysis_name)
        plugin = analysis_cls()

        resources = _resolve_submit_resources_with_hints(
            plugin=plugin,
            pixi_path=pixi_path,
            partition=partition,
            qos=qos,
            account=account,
            ntasks=ntasks,
            cpus_per_task=cpus_per_task,
            mem=mem,
            time_limit=time_limit,
            max_retries=max_retries,
            mail_user=mail_user,
            cpus_from_cli=(
                ctx.get_parameter_source("cpus_per_task") != click.core.ParameterSource.DEFAULT
            ),
            mem_from_cli=ctx.get_parameter_source("mem") != click.core.ParameterSource.DEFAULT,
            time_from_cli=ctx.get_parameter_source("time_limit")
            != click.core.ParameterSource.DEFAULT,
        )

        manifest = build_manifest(
            plugin,
            config,
            resources,
            recompute,
            equilibration,
            allow_partial=allow_partial,
        )
        hpc_dir = _resolve_hpc_dir(config, plugin.name)
        hpc_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = hpc_dir / "manifest.json"
        manifest.save(manifest_path)

        dependencies = tuple(getattr(analysis_cls, "dependencies", ()))
        root_dependencies = [finalize_ids[dep] for dep in dependencies if dep in finalize_ids]

        replicate_count = sum(len(cond.replicate_specs) for cond in manifest.condition_specs)
        aggregate_count = len(manifest.condition_specs)
        finalize_only = getattr(manifest, "pipeline_mode", "full") == "finalize_only"

        if dry_run:
            if not finalize_only:
                for cond in manifest.condition_specs:
                    for rep in cond.replicate_specs:
                        generate_replicate_script(manifest, rep, resources, hpc_dir)
                    generate_aggregate_script(manifest, cond, resources, hpc_dir)
            generate_finalize_script(manifest, resources, hpc_dir)

            fake_finalize_id = f"dry-run:{plugin.name}:finalize"
            finalize_ids[plugin.name] = fake_finalize_id
            total_jobs = 1 if finalize_only else (replicate_count + aggregate_count + 1)
            summary_rows.append(
                {
                    "analysis": plugin.name,
                    "mode": manifest.pipeline_mode,
                    "jobs": str(total_jobs),
                    "finalize_job_id": fake_finalize_id,
                }
            )
            continue

        graph = submit_analysis_graph(
            manifest,
            resources,
            hpc_dir,
            root_dependencies=tuple(root_dependencies),
        )
        graph.save(hpc_dir / "job_graph.json")
        finalize_ids[plugin.name] = graph.finalizer_job_id
        total_jobs = 1 if finalize_only else (replicate_count + aggregate_count + 1)
        summary_rows.append(
            {
                "analysis": plugin.name,
                "mode": manifest.pipeline_mode,
                "jobs": str(total_jobs),
                "finalize_job_id": graph.finalizer_job_id,
            }
        )

    _echo_submit_all_summary(summary_rows)
    if dry_run:
        click.echo("Dry run only: no jobs were submitted")


@compare.command("status")
@click.argument("analysis", type=str)
@click.option(
    "-f",
    "--file",
    "config_file",
    type=click.Path(path_type=Path),
    default="comparison.yaml",
    help="Path to comparison.yaml config file.",
)
@click.option(
    "--reconcile",
    is_flag=True,
    help="Reconcile pending/running status files with sacct before reporting.",
)
@click.option("--json", "as_json", is_flag=True, help="Print machine-readable JSON status.")
def analysis_hpc_status(analysis: str, config_file: Path, reconcile: bool, as_json: bool):
    """Show status for submitted analysis SLURM DAG."""
    from polyzymd.analyses.discovery import get_analysis
    from polyzymd.workflow.analysis_slurm import read_analysis_status, reconcile_status_with_slurm

    try:
        config = ComparisonConfig.from_yaml(config_file)
    except (FileNotFoundError, yaml.YAMLError, ValidationError, ValueError) as e:
        raise click.ClickException(f"Error loading config: {e}") from e

    try:
        analysis_cls = get_analysis(analysis)
    except KeyError:
        if get_archived_analysis_plugin(analysis) is not None:
            raise _archived_analysis_click_exception(
                analysis,
                context="CLI status lookup",
            ) from None

        from polyzymd.analyses.discovery import list_all_names

        available = list_all_names()
        raise click.ClickException(
            f"Unknown analysis type '{analysis}'. Available: {', '.join(sorted(available))}"
        )
    hpc_dir = _resolve_hpc_dir(config, analysis_cls.name)

    reconciliation_summary: dict[str, Any] | None = None
    if reconcile:
        reconciliation_summary = reconcile_status_with_slurm(hpc_dir)

    status = read_analysis_status(hpc_dir)
    if as_json:
        if reconcile and reconciliation_summary is not None:
            status = {**status, "reconciliation": reconciliation_summary}
        click.echo(json.dumps(status, indent=2))
        return

    if reconcile and reconciliation_summary is not None:
        change_counts: dict[tuple[str, str], int] = {}
        for change in reconciliation_summary.get("changes", []):
            key = (change["new_state"], change["slurm_state"])
            change_counts[key] = change_counts.get(key, 0) + 1

        if change_counts:
            fragments = [
                f"{count} marked {new_state} (SLURM: {slurm_state})"
                for (new_state, slurm_state), count in sorted(change_counts.items())
            ]
            click.echo(
                f"Reconciled {reconciliation_summary.get('updated', 0)} jobs: "
                + ", ".join(fragments)
            )
        else:
            click.echo(
                f"Reconciled 0 jobs (checked {reconciliation_summary.get('checked', 0)} "
                "pending/running tasks)"
            )

    counts = status.get("counts", {})
    click.echo(f"Analysis: {analysis_cls.name}")
    click.echo(f"HPC dir: {hpc_dir}")
    click.echo(
        "States: "
        f"pending={counts.get('pending', 0)} "
        f"running={counts.get('running', 0)} "
        f"retrying={counts.get('retrying', 0)} "
        f"succeeded={counts.get('succeeded', 0)} "
        f"failed={counts.get('failed', 0)} "
        f"unknown={counts.get('unknown', 0)}"
    )
    warnings = status.get("warnings", [])
    if warnings:
        click.echo("⚠ Warnings:")
        for warning in warnings:
            click.echo(f"  - {warning}")


@compare.command("finalize")
@click.argument("analysis", type=str)
@click.option(
    "-f",
    "--file",
    "config_file",
    type=click.Path(path_type=Path),
    default="comparison.yaml",
    help="Path to comparison.yaml config file.",
)
@click.option("--recompute", is_flag=True, help="Regenerate comparison and plot outputs.")
@click.option(
    "--allow-partial",
    is_flag=True,
    help="Allow finalize to proceed when some conditions are missing aggregated results.",
)
def finalize_analysis_hpc(
    analysis: str,
    config_file: Path,
    recompute: bool,
    allow_partial: bool,
):
    """Run comparison + plotting from aggregated on-disk results."""
    from polyzymd.analyses.discovery import get_analysis
    from polyzymd.analyses.orchestrator import (
        finalize_comparison_from_disk,
        prepare_comparison_run,
    )
    from polyzymd.analyses.shared.paths import sanitize_label

    try:
        config = ComparisonConfig.from_yaml(config_file)
    except (FileNotFoundError, yaml.YAMLError, ValidationError, ValueError) as e:
        raise click.ClickException(f"Error loading config: {e}") from e

    try:
        analysis_cls = get_analysis(analysis)
    except KeyError:
        if get_archived_analysis_plugin(analysis) is not None:
            raise _archived_analysis_click_exception(
                analysis,
                context="CLI finalize lookup",
            ) from None

        from polyzymd.analyses.discovery import list_all_names

        available = list_all_names()
        raise click.ClickException(
            f"Unknown analysis type '{analysis}'. Available: {', '.join(sorted(available))}"
        )
    plugin = analysis_cls()
    prepared = prepare_comparison_run(
        plugin,
        config,
        None,
    )
    valid_conditions = prepared["valid_conditions"]
    settings = prepared["settings"]
    equilibration = prepared["equilibration"]
    analysis_root = prepared["analysis_root"]
    analysis_dirs: dict[str, Path] = {}
    aggregated_results: dict[str, object] = {}
    expects_aggregated_results = bool(getattr(plugin, "has_compute_stage", True))
    if not getattr(plugin, "has_aggregate_stage", True):
        expects_aggregated_results = False

    for condition in valid_conditions:
        cond_dir = analysis_root / sanitize_label(condition.label) / plugin.name
        analysis_dirs[condition.label] = cond_dir
        if expects_aggregated_results:
            aggregated = plugin._load_aggregated_result(cond_dir / "aggregated")
            if aggregated is not None:
                aggregated_results[condition.label] = aggregated

    missing_conditions: list[str] = []
    if expects_aggregated_results:
        missing_conditions = [
            condition.label
            for condition in valid_conditions
            if condition.label not in aggregated_results
        ]

    if missing_conditions:
        expected_paths = []
        for condition in valid_conditions:
            if condition.label in missing_conditions:
                cond_dir = analysis_root / sanitize_label(condition.label) / plugin.name
                aggregate_result_path = getattr(plugin, "aggregate_result_path", None)
                if callable(aggregate_result_path):
                    expected_path = aggregate_result_path(cond_dir / "aggregated")
                else:
                    expected_path = cond_dir / "aggregated" / "result.json"
                expected_paths.append(str(expected_path))
        click.echo(
            "Warning: missing aggregated results for condition(s): "
            f"{', '.join(missing_conditions)}\n"
            "Expected aggregate result path(s): "
            f"{', '.join(expected_paths)}",
            err=True,
        )
        if not allow_partial:
            raise click.ClickException(
                "Finalize aborted due to missing aggregated results. "
                "Re-run failed aggregate jobs or pass --allow-partial to continue "
                "(--allow-partial (CLI) / allow_partial=True (API))."
            )

    effective_control = (
        config.control
        if config.control and any(cond.label == config.control for cond in valid_conditions)
        else None
    )

    plot_settings = config.plot_settings
    figures_base = plot_settings.output_dir
    if not figures_base.is_absolute() and config.source_path is not None:
        figures_base = config.source_path.parent / figures_base
    figures_dir = plugin.figures_output_dir(figures_base)
    results_dir = analysis_root.parent / "comparison" / plugin.name

    config_for_finalize = config.model_copy(deep=True)
    config_for_finalize.defaults.equilibration_time = equilibration
    try:
        result = finalize_comparison_from_disk(
            analysis=plugin,
            config=config_for_finalize,
            analysis_dirs=analysis_dirs,
            aggregated_results=aggregated_results,
            results_dir=results_dir,
            figures_dir=figures_dir,
            settings=settings,
            effective_control=effective_control,
            allow_partial=allow_partial,
            recompute=recompute,
        )
    except PluginContractError as e:
        raise _plugin_contract_click_exception(plugin.name, e) from e
    click.echo(f"Saved result: {result['comparison_path']}")


@compare.command("worker-replicate", hidden=True)
@click.option("--manifest", "manifest_path", type=click.Path(path_type=Path), required=True)
@click.option("--condition-index", type=int, required=True)
@click.option("--replicate", type=int, required=True)
def worker_replicate(
    manifest_path: Path,
    condition_index: int,
    replicate: int,
):
    """Internal worker command for one replicate compute task."""
    from polyzymd.analyses.discovery import get_analysis
    from polyzymd.analyses.orchestrator import run_replicate_once
    from polyzymd.analyses.shared.paths import sanitize_label
    from polyzymd.workflow.analysis_slurm import AnalysisJobManifest, validate_manifest_snapshot

    manifest = AnalysisJobManifest.load(manifest_path)
    config = ComparisonConfig.from_yaml(Path(manifest.comparison_yaml))
    try:
        analysis_cls = get_analysis(manifest.analysis_name)
    except KeyError:
        if get_archived_analysis_plugin(manifest.analysis_name) is not None:
            raise _archived_analysis_click_exception(
                manifest.analysis_name,
                context="CLI worker lookup",
            ) from None
        raise
    plugin = analysis_cls()
    valid_conditions, equilibration, analysis_root = validate_manifest_snapshot(
        manifest,
        plugin,
        config,
    )
    settings = _settings_from_manifest(plugin, manifest)
    cond_spec = _resolve_manifest_task_condition(manifest, condition_index)
    condition = valid_conditions[condition_index]
    if condition.label != cond_spec.condition_label:
        raise click.ClickException(
            "Manifest/config drift detected: condition labels no longer align with submission"
        )
    if replicate not in [spec.replicate for spec in cond_spec.replicate_specs]:
        raise click.ClickException(
            f"Replicate {replicate} was not present in submitted manifest for "
            f"condition '{cond_spec.condition_label}'"
        )
    cond_dir = analysis_root / sanitize_label(condition.label) / plugin.name
    run_dir = cond_dir / f"run_{replicate}"
    result = run_replicate_once(
        plugin,
        condition,
        settings,
        equilibration,
        run_dir,
        replicate,
        getattr(manifest, "recompute", False),
    )
    if result is None:
        raise click.ClickException("Replicate computation produced no result")


@compare.command("worker-aggregate", hidden=True)
@click.option("--manifest", "manifest_path", type=click.Path(path_type=Path), required=True)
@click.option("--condition-index", type=int, required=True)
def worker_aggregate(manifest_path: Path, condition_index: int):
    """Internal worker command for one condition aggregation task."""
    from polyzymd.analyses.discovery import get_analysis
    from polyzymd.analyses.orchestrator import aggregate_condition_from_disk
    from polyzymd.analyses.shared.paths import sanitize_label
    from polyzymd.workflow.analysis_slurm import AnalysisJobManifest, validate_manifest_snapshot

    manifest = AnalysisJobManifest.load(manifest_path)
    config = ComparisonConfig.from_yaml(Path(manifest.comparison_yaml))
    try:
        analysis_cls = get_analysis(manifest.analysis_name)
    except KeyError:
        if get_archived_analysis_plugin(manifest.analysis_name) is not None:
            raise _archived_analysis_click_exception(
                manifest.analysis_name,
                context="CLI worker lookup",
            ) from None
        raise
    plugin = analysis_cls()
    valid_conditions, equilibration, analysis_root = validate_manifest_snapshot(
        manifest,
        plugin,
        config,
    )
    settings = _settings_from_manifest(plugin, manifest)
    cond_spec = _resolve_manifest_task_condition(manifest, condition_index)
    condition = valid_conditions[condition_index]
    if condition.label != cond_spec.condition_label:
        raise click.ClickException(
            "Manifest/config drift detected: condition labels no longer align with submission"
        )
    cond_dir = analysis_root / sanitize_label(condition.label) / plugin.name
    aggregate_condition_from_disk(
        plugin,
        condition,
        settings,
        equilibration,
        cond_dir,
        tuple(spec.replicate for spec in cond_spec.replicate_specs),
        recompute=getattr(manifest, "recompute", False),
    )


@compare.command("worker-finalize", hidden=True)
@click.option("--manifest", "manifest_path", type=click.Path(path_type=Path), required=True)
def worker_finalize(manifest_path: Path):
    """Internal worker command for comparison finalization."""
    from polyzymd.analyses.discovery import get_analysis
    from polyzymd.analyses.orchestrator import finalize_comparison_from_disk
    from polyzymd.analyses.shared.paths import sanitize_label
    from polyzymd.workflow.analysis_slurm import AnalysisJobManifest, validate_manifest_snapshot

    manifest = AnalysisJobManifest.load(manifest_path)
    config = ComparisonConfig.from_yaml(Path(manifest.comparison_yaml))
    try:
        analysis_cls = get_analysis(manifest.analysis_name)
    except KeyError:
        if get_archived_analysis_plugin(manifest.analysis_name) is not None:
            raise _archived_analysis_click_exception(
                manifest.analysis_name,
                context="CLI worker lookup",
            ) from None
        raise
    plugin = analysis_cls()
    valid_conditions, equilibration, analysis_root = validate_manifest_snapshot(
        manifest,
        plugin,
        config,
    )
    settings = _settings_from_manifest(plugin, manifest)

    analysis_dirs: dict[str, Path] = {}
    aggregated_results: dict[str, object] = {}
    missing_conditions: list[str] = []
    expects_aggregated_results = bool(getattr(plugin, "has_compute_stage", True))
    if getattr(manifest, "pipeline_mode", "full") == "finalize_only":
        expects_aggregated_results = False

    if expects_aggregated_results:
        for cond_idx, condition in enumerate(valid_conditions):
            cond_spec = _resolve_manifest_task_condition(manifest, cond_idx)
            if condition.label != cond_spec.condition_label:
                raise click.ClickException(
                    "Manifest/config drift detected: condition labels no longer align with submission"
                )
            cond_dir = analysis_root / sanitize_label(condition.label) / plugin.name
            aggregated = plugin._load_aggregated_result(cond_dir / "aggregated")
            if aggregated is not None:
                analysis_dirs[condition.label] = cond_dir
                aggregated_results[condition.label] = aggregated
            else:
                missing_conditions.append(condition.label)
    else:
        for cond_idx, condition in enumerate(valid_conditions):
            cond_spec = _resolve_manifest_task_condition(manifest, cond_idx)
            if condition.label != cond_spec.condition_label:
                raise click.ClickException(
                    "Manifest/config drift detected: condition labels no longer align with submission"
                )
            cond_dir = analysis_root / sanitize_label(condition.label) / plugin.name
            analysis_dirs[condition.label] = cond_dir

    partial_policy = getattr(manifest, "partial_policy", "strict")
    allow_partial = partial_policy == "allow_partial"
    if expects_aggregated_results and missing_conditions:
        if allow_partial:
            message = "Proceeding with partial finalize"
        else:
            message = "Strict policy enabled; finalize will fail"
        click.echo(
            "Warning: missing aggregated results for condition(s): "
            f"{', '.join(missing_conditions)}. {message}.",
            err=True,
        )
        if not allow_partial:
            raise click.ClickException(
                "Finalize aborted due to missing aggregated results under strict policy. "
                "Re-run failed jobs or submit with --allow-partial."
            )

    effective_control = (
        config.control
        if config.control and any(cond.label == config.control for cond in valid_conditions)
        else None
    )

    plot_settings = config.plot_settings
    figures_base = plot_settings.output_dir
    if not figures_base.is_absolute() and config.source_path is not None:
        figures_base = config.source_path.parent / figures_base

    config_for_finalize = config.model_copy(deep=True)
    config_for_finalize.defaults.equilibration_time = equilibration
    finalize_comparison_from_disk(
        analysis=plugin,
        config=config_for_finalize,
        analysis_dirs=analysis_dirs,
        aggregated_results=aggregated_results,
        results_dir=analysis_root.parent / "comparison" / plugin.name,
        figures_dir=plugin.figures_output_dir(figures_base),
        settings=settings,
        effective_control=effective_control,
        allow_partial=allow_partial,
        recompute=getattr(manifest, "recompute", False),
    )
