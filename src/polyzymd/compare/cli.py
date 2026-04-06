"""CLI commands for the compare module.

This module provides the `polyzymd compare` command group with subcommands
for initializing comparison projects and running comparisons.
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

from polyzymd.cli.colors import echo_logo
from polyzymd.compare.cli_utils import (
    common_compare_options,
    load_comparison_config,
    validate_and_report,
)
from polyzymd.compare.config import (
    ComparisonConfig,
    generate_comparison_template,
)
from polyzymd.core.branding import prepend_file_header
from polyzymd.core.experimental import (
    echo_experimental_warning,
    experimental_features_for_comparison_type,
    format_experimental_suffix,
    normalize_experimental_features,
)

LOGGER = logging.getLogger("polyzymd.compare.cli")


def _resolve_hpc_dir(config: ComparisonConfig, analysis_name: str) -> Path:
    """Return the analysis HPC artifact directory."""
    source = config.source_path
    if source is None:
        return Path("comparison") / analysis_name / "_hpc"
    return source.parent / "comparison" / analysis_name / "_hpc"


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
    except Exception as exc:
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
    Workflow:
    1. polyzymd compare init -n <name>   # Create project with template
    2. Edit comparison.yaml              # Add your conditions
    3. polyzymd compare run <analysis>   # Run a single comparison
    4. polyzymd compare run-all          # Run all enabled comparisons

    \b
    Example:
        polyzymd compare init -n polymer_study
        cd polymer_study
        vim comparison.yaml  # Add your conditions
        polyzymd compare run rmsf --eq-time 10ns
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

## For Binding Preference Analysis

Copy your enzyme PDB file here for SASA (solvent-accessible surface area)
calculation used in binding preference analysis:

    cp /path/to/your/enzyme.pdb structures/

Then reference it in comparison.yaml:

    plugins:
      contacts:
        compute_binding_preference: true
        enzyme_pdb_for_sasa: "structures/enzyme.pdb"

The enzyme PDB should be the reference structure (e.g., from PDB or your
prepared simulation input), NOT a trajectory frame.
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
        click.echo(f"  1. Edit {config_path.relative_to(Path.cwd())}")
        click.echo("     - Add your simulation conditions (paths to config.yaml files)")
        click.echo("     - Define catalytic_triad for active site analysis")
        click.echo()
        click.echo("  2. For binding preference analysis, copy your enzyme PDB:")
        click.echo(f"     cp /path/to/enzyme.pdb {project_dir.relative_to(Path.cwd())}/structures/")
        click.echo()
        click.echo(f"  3. cd {project_dir.relative_to(Path.cwd())}")
        click.echo("  4. Run comparisons:")
        click.echo("     polyzymd compare run rmsf      # Compare flexibility")
        click.echo("     polyzymd compare run triad     # Compare triad geometry")
        click.echo("     polyzymd compare run contacts  # Compare polymer-protein contacts")
        click.echo("     polyzymd compare run exposure  # Compare chaperone-like activity")
        click.echo()

    except Exception as e:
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
      - At least 2 conditions defined
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
    except Exception as e:
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

    This is a generic command that can run any discovered analysis plugin.
    Use --list to see available comparison types.

    \b
    Available comparison types:
      - rmsf: Compare RMSF (flexibility) across conditions
      - triad: Compare catalytic triad geometry across conditions
      - contacts: Compare polymer-protein contacts across conditions
      - exposure: Compare chaperone-like polymer activity across conditions

    \b
    Example:
        polyzymd compare run rmsf
        polyzymd compare run triad --eq-time 10ns
        polyzymd compare run contacts --format markdown
        polyzymd compare run exposure --format json
        polyzymd compare run --list
    """
    from polyzymd.analyses.discovery import get_analysis, list_all_names, list_analyses
    from polyzymd.analyses.orchestrator import run_comparison as _run_pipeline
    from polyzymd.analyses.shared.logging_utils import setup_logging

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
    except Exception as e:
        click.echo(f"Error during comparison: {e}", err=True)
        if debug:
            import traceback

            traceback.print_exc()
        sys.exit(1)

    if result is None:
        click.echo("Warning: comparison returned no result (not enough data?).", err=True)
        sys.exit(1)

    json_path = pipeline_result.get("comparison_path")
    if json_path is not None:
        click.echo(f"Saved result: {json_path}")
    click.echo()

    # Format and display — delegates to the plugin's format() method
    # Map CLI "table" format to plugin "text" format
    fmt = "text" if output_format == "table" else output_format
    try:
        formatted = analysis.format(result, output_format=fmt)
    except Exception as e:
        click.echo(f"Warning: Could not format result: {e}", err=True)
        if hasattr(result, "model_dump_json"):
            formatted = result.model_dump_json(indent=2)
        else:
            formatted = str(result)

    click.echo(formatted)

    # Save formatted output if requested
    if output_path:
        output_path = Path(output_path)
        output_path.write_text(formatted)
        click.echo(f"Saved output: {output_path}")


def _generate_analysis_plots(
    config: ComparisonConfig,
    analysis_names: list[str] | None = None,
) -> list[Path]:
    """Generate plots for analysis plugins using the plugin system.

    Builds :class:`~polyzymd.analyses.base.PlotContext` objects for each
    enabled analysis and calls the plugin's ``plot()`` method.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration loaded from comparison.yaml.
    analysis_names : list[str] | None
        Analyses to plot.  ``None`` means all enabled analyses.

    Returns
    -------
    list[Path]
        Paths to all generated figure files.
    """
    from polyzymd.analyses.base import Condition, PlotContext
    from polyzymd.analyses.discovery import get_analysis
    from polyzymd.analyses.orchestrator import _resolve_settings
    from polyzymd.compare.io.paths import sanitize_label

    if analysis_names is None:
        analysis_names = config.plugins.get_enabled_plugins()

    # Build Condition objects
    conditions = [Condition.from_condition_config(c) for c in config.conditions]

    # Resolve effective control label
    control_label = config.control
    effective_control = (
        control_label
        if control_label and any(c.label == control_label for c in conditions)
        else None
    )

    # Resolve paths relative to comparison.yaml
    source_path = config.source_path
    analysis_root = source_path.parent / "analysis" if source_path else Path("analysis")

    # Resolve output directory for figures
    plot_settings = config.plot_settings
    figures_base = plot_settings.output_dir
    if not figures_base.is_absolute():
        if source_path is not None:
            figures_base = source_path.parent / figures_base
        else:
            figures_base = Path.cwd() / figures_base
    figures_base = figures_base.resolve()

    generated: list[Path] = []

    for name in analysis_names:
        try:
            analysis_cls = get_analysis(name)
        except KeyError:
            LOGGER.warning(f"Unknown analysis type {name!r} — skipping plots.")
            continue

        analysis = analysis_cls()
        settings = _resolve_settings(analysis, config)

        # Filter conditions (pass resolved settings for plugins that need them)
        valid_conditions = analysis.filter_conditions(conditions, settings=settings)

        # Build analysis_dirs mapping (mirrors orchestrator.run_comparison)
        analysis_dirs: dict[str, Path] = {}
        for cond in valid_conditions:
            analysis_dirs[cond.label] = analysis_root / sanitize_label(cond.label) / analysis.name

        # Comparison results dir
        results_dir = analysis_root.parent / "comparison" / analysis.name
        comparison_result_path = analysis.comparison_result_path(results_dir)

        # Figures dir for this analysis
        figures_dir = analysis.figures_output_dir(figures_base)
        figures_dir.mkdir(parents=True, exist_ok=True)

        plot_ctx = PlotContext(
            conditions=valid_conditions,
            analysis_dirs=analysis_dirs,
            results_dir=results_dir,
            output_dir=figures_dir,
            settings=settings,
            plot_settings=plot_settings,
            comparison_path=comparison_result_path,
            control_label=effective_control,
        )

        try:
            paths = analysis.plot(plot_ctx)
            generated.extend(paths)
        except Exception as e:
            LOGGER.error(f"Failed to generate plots for {name}: {e}")

    return generated


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
    except Exception as e:
        click.echo(f"Error loading config: {e}", err=True)
        sys.exit(1)

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

    # Generate plots
    try:
        click.echo(f"Generating plots for: {', '.join(target_analyses)}...")
        generated = _generate_analysis_plots(config, target_analyses)
    except Exception as e:
        click.echo(f"Error generating plots: {e}", err=True)
        if debug:
            import traceback

            traceback.print_exc()
        sys.exit(1)

    click.echo()
    if generated:
        click.echo(f"Generated {len(generated)} plots:")
        for path in generated:
            click.echo(f"  - {path}")
    else:
        click.echo("No plots generated. Check that analyses are enabled in config.")


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
    from polyzymd.analyses.shared.logging_utils import setup_logging

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
    all_results = run_all_comparisons(
        config,
        analysis_names=None,  # run all enabled
        recompute=recompute,
        equilibration=equilibration,
    )

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
        click.echo()
        click.echo("Generating plots ...")

        try:
            generated = _generate_analysis_plots(config, succeeded)
            click.echo(f"Generated {len(generated)} plots.")
        except Exception as e:
            click.echo(f"Error generating plots: {e}", err=True)
            if debug:
                import traceback

                traceback.print_exc()

    # Exit with error code if any comparison failed
    if failed:
        sys.exit(1)


@compare.command("submit")
@click.argument("analysis", type=str)
@click.option(
    "--comparison-yaml",
    "comparison_yaml",
    type=click.Path(path_type=Path),
    required=True,
    help="Path to comparison.yaml config file.",
)
@click.option("--partition", default="aa100", help="SLURM partition.")
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
    analysis: str,
    comparison_yaml: Path,
    partition: str,
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
        AnalysisSlurmResources,
        build_manifest,
        generate_aggregate_script,
        generate_array_script,
        generate_finalize_script,
        generate_replicate_script,
        submit_analysis_graph,
        submit_analysis_graph_with_arrays,
    )

    if not dry_run and shutil.which("sbatch") is None:
        raise click.ClickException(
            "SLURM is not available: 'sbatch' not found on PATH. The HPC submission "
            "commands require a SLURM cluster. Run analysis locally with "
            "'polyzymd compare run' instead."
        )

    config = ComparisonConfig.from_yaml(comparison_yaml)
    analysis_cls = get_analysis(analysis)
    plugin = analysis_cls()
    resources = AnalysisSlurmResources(
        pixi_path=pixi_path,
        partition=partition,
        qos=qos,
        account=account,
        ntasks=ntasks,
        cpus_per_task=cpus_per_task,
        mem=mem,
        time=time_limit,
        max_retries=max_retries,
        mail_user=mail_user,
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

    if dry_run:
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
        if job_arrays:
            total = array_count + aggregate_count + 1
            click.echo(
                "Submitted "
                f"{array_count} array jobs + {aggregate_count} aggregate + 1 finalize = {total} total"
            )
            click.echo("Submission mode: job arrays")
        else:
            total = replicate_count + aggregate_count + 1
            click.echo(
                "Submitted "
                f"{total} jobs ({replicate_count} replicate + {aggregate_count} aggregate + 1 finalize)"
            )
        click.echo("Dry run only: no jobs were submitted")
        return

    if job_arrays:
        graph = submit_analysis_graph_with_arrays(manifest, resources, hpc_dir)
    else:
        graph = submit_analysis_graph(manifest, resources, hpc_dir)
    graph.save(hpc_dir / "job_graph.json")
    if job_arrays:
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


@compare.command("status")
@click.argument("analysis", type=str)
@click.option(
    "--comparison-yaml",
    "comparison_yaml",
    type=click.Path(path_type=Path),
    required=True,
    help="Path to comparison.yaml config file.",
)
@click.option(
    "--reconcile",
    is_flag=True,
    help="Reconcile pending/running status files with sacct before reporting.",
)
@click.option("--json", "as_json", is_flag=True, help="Print machine-readable JSON status.")
def analysis_hpc_status(analysis: str, comparison_yaml: Path, reconcile: bool, as_json: bool):
    """Show status for submitted analysis SLURM DAG."""
    from polyzymd.analyses.discovery import get_analysis
    from polyzymd.workflow.analysis_slurm import read_analysis_status, reconcile_status_with_slurm

    config = ComparisonConfig.from_yaml(comparison_yaml)
    analysis_cls = get_analysis(analysis)
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
        f"completed={counts.get('completed', 0)} "
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
    "--comparison-yaml",
    "comparison_yaml",
    type=click.Path(path_type=Path),
    required=True,
    help="Path to comparison.yaml config file.",
)
@click.option("--recompute", is_flag=True, help="Retained for CLI compatibility.")
@click.option(
    "--allow-partial",
    is_flag=True,
    help="Allow finalize to proceed when some conditions are missing aggregated results.",
)
def finalize_analysis_hpc(
    analysis: str,
    comparison_yaml: Path,
    recompute: bool,
    allow_partial: bool,
):
    """Run comparison + plotting from aggregated on-disk results."""
    del recompute
    from polyzymd.analyses.discovery import get_analysis
    from polyzymd.analyses.orchestrator import (
        finalize_comparison_from_disk,
        prepare_comparison_run,
    )
    from polyzymd.compare.io.paths import sanitize_label

    config = ComparisonConfig.from_yaml(comparison_yaml)
    analysis_cls = get_analysis(analysis)
    plugin = analysis_cls()
    valid_conditions, settings, equilibration, analysis_root = prepare_comparison_run(
        plugin,
        config,
        None,
    )
    analysis_dirs: dict[str, Path] = {}
    aggregated_results: dict[str, object] = {}
    for condition in valid_conditions:
        cond_dir = analysis_root / sanitize_label(condition.label) / plugin.name
        aggregated = plugin._load_aggregated_result(cond_dir / "aggregated")
        if aggregated is not None:
            analysis_dirs[condition.label] = cond_dir
            aggregated_results[condition.label] = aggregated

    missing_conditions = [
        condition.label
        for condition in valid_conditions
        if condition.label not in aggregated_results
    ]
    if missing_conditions:
        click.echo(
            "Warning: missing aggregated results for condition(s): "
            f"{', '.join(missing_conditions)}",
            err=True,
        )
        if not allow_partial:
            raise click.ClickException(
                "Finalize aborted due to missing aggregated results. "
                "Re-run failed jobs or pass --allow-partial to continue."
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
    )
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
    from polyzymd.compare.io.paths import sanitize_label
    from polyzymd.workflow.analysis_slurm import AnalysisJobManifest, validate_manifest_snapshot

    manifest = AnalysisJobManifest.load(manifest_path)
    config = ComparisonConfig.from_yaml(Path(manifest.comparison_yaml))
    analysis_cls = get_analysis(manifest.analysis_name)
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
        manifest.recompute,
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
    from polyzymd.compare.io.paths import sanitize_label
    from polyzymd.workflow.analysis_slurm import AnalysisJobManifest, validate_manifest_snapshot

    manifest = AnalysisJobManifest.load(manifest_path)
    config = ComparisonConfig.from_yaml(Path(manifest.comparison_yaml))
    analysis_cls = get_analysis(manifest.analysis_name)
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
    )


@compare.command("worker-finalize", hidden=True)
@click.option("--manifest", "manifest_path", type=click.Path(path_type=Path), required=True)
def worker_finalize(manifest_path: Path):
    """Internal worker command for comparison finalization."""
    from polyzymd.analyses.discovery import get_analysis
    from polyzymd.analyses.orchestrator import finalize_comparison_from_disk
    from polyzymd.compare.io.paths import sanitize_label
    from polyzymd.workflow.analysis_slurm import AnalysisJobManifest, validate_manifest_snapshot

    manifest = AnalysisJobManifest.load(manifest_path)
    config = ComparisonConfig.from_yaml(Path(manifest.comparison_yaml))
    analysis_cls = get_analysis(manifest.analysis_name)
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

    partial_policy = getattr(manifest, "partial_policy", "strict")
    allow_partial = partial_policy == "allow_partial"
    if missing_conditions:
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
    )
