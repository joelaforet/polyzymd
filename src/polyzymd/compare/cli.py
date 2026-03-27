"""CLI commands for the compare module.

This module provides the `polyzymd compare` command group with subcommands
for initializing comparison projects and running comparisons.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

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
    experimental_features_for_plot_type,
    format_experimental_suffix,
    normalize_experimental_features,
)

LOGGER = logging.getLogger("polyzymd.compare.cli")


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
def init(name: str, eq_time: str, output_dir: Optional[Path]):
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
        click.echo("     polyzymd compare exposure      # Compare chaperone-like activity")
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
    import json as json_module

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
    comparison_type: Optional[str],
    config_file: Path,
    eq_time: Optional[str],
    recompute: bool,
    output_format: str,
    output_path: Optional[Path],
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
    "-p",
    "--plot-type",
    type=str,
    default=None,
    help="Generate specific plot type only (e.g., 'triad_kde_panel').",
)
@click.option(
    "--list-available",
    is_flag=True,
    help="List available plot types for enabled analyses.",
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
    output_dir: Optional[Path],
    analysis_type: Optional[str],
    plot_type: Optional[str],
    list_available: bool,
    quiet: bool,
    debug: bool,
):
    """Generate all plots from comparison.yaml configuration.

    This is the config-driven plotting command that reads plot_settings
    from comparison.yaml and generates all configured plots automatically.

    \b
    Available plot types (registered via PlotterRegistry):
      - triad_kde_panel: Multi-row KDE panel for catalytic triad distances
      - triad_threshold_bars: Grouped bar chart of triad contact fractions
      - rmsf_comparison: Bar chart comparing whole-protein RMSF
      - rmsf_profile: Per-residue RMSF line plot
      - distance_kde: KDE distribution plots for distance pairs
      - distance_threshold_bars: Contact fraction bar chart
      - bfe_heatmap: ΔG_sel heatmap (AA groups × conditions, per polymer type)
      - bfe_bars: ΔG_sel grouped bar chart with ±k_BT reference lines

    \b
    Examples:
        polyzymd compare plot-all -f comparison.yaml
        polyzymd compare plot-all -f comparison.yaml -a catalytic_triad
        polyzymd compare plot-all -f comparison.yaml -p triad_kde_panel
        polyzymd compare plot-all --list-available
    """
    from polyzymd.compare.plotter import ComparisonPlotter, PlotterRegistry

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

    # Create plotter
    plotter = ComparisonPlotter(config)

    # List available plots
    if list_available:
        click.echo("Registered plot types:")
        for ptype in PlotterRegistry.list_available():
            suffix = format_experimental_suffix(experimental_features_for_plot_type(ptype))
            click.echo(f"  - {ptype}{suffix}")
        click.echo()
        click.echo("Available plots for enabled analyses:")
        available = plotter.list_available_plots()
        for atype, ptypes in available.items():
            click.echo(f"  {atype}:")
            for pt in ptypes:
                suffix = format_experimental_suffix(experimental_features_for_plot_type(pt))
                click.echo(f"    - {pt}{suffix}")
        return

    experimental_features: tuple[str, ...]
    if plot_type:
        experimental_features = experimental_features_for_plot_type(plot_type)
    elif analysis_type:
        comparison_type = analysis_type
        analysis_settings = config.plugins.get(analysis_type)
        experimental_features = experimental_features_for_comparison_type(
            comparison_type, analysis_settings
        )
    else:
        feature_list: list[str] = []
        for settings_key in config.plugins.get_enabled_plugins():
            comparison_type = settings_key
            analysis_settings = config.plugins.get(settings_key)
            feature_list.extend(
                experimental_features_for_comparison_type(comparison_type, analysis_settings)
            )
        experimental_features = normalize_experimental_features(feature_list)

    echo_experimental_warning(experimental_features)

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Conditions: {len(config.conditions)}")
    click.echo(f"Output directory: {plotter.output_dir}")
    click.echo()

    # Generate plots
    try:
        if plot_type and analysis_type:
            # Specific plot type for specific analysis
            click.echo(f"Generating {plot_type} for {analysis_type}...")
            generated = plotter.plot_single(plot_type, analysis_type)
        elif analysis_type:
            # All plots for specific analysis
            click.echo(f"Generating plots for {analysis_type}...")
            generated = plotter.plot_analysis(analysis_type)
        else:
            # All plots for all analyses
            click.echo("Generating all plots...")
            generated = plotter.plot_all()
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
    from polyzymd.analyses.discovery import get_analysis
    from polyzymd.analyses.orchestrator import run_comparison as _run_pipeline
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

    # --- Run each analysis plugin -------------------------------------------
    succeeded: list[str] = []
    failed: list[tuple[str, str]] = []
    skipped: list[str] = []

    for settings_key in enabled_analyses:
        # Resolve plugin — try the settings_key directly as an analysis name
        try:
            analysis_cls = get_analysis(settings_key)
        except KeyError:
            click.echo(f"  [{settings_key}] skipped — no registered analysis plugin")
            skipped.append(settings_key)
            continue

        click.echo(f"  [{settings_key}] running {analysis_cls.name} comparison ...")

        try:
            analysis = analysis_cls()
            pipeline_result = _run_pipeline(
                analysis, config, recompute=recompute, equilibration=equilibration
            )
            result = pipeline_result["comparison"]
        except Exception as e:
            msg = str(e)
            click.echo(f"  [{settings_key}] FAILED: {msg}", err=True)
            if debug:
                import traceback

                traceback.print_exc()
            failed.append((settings_key, msg))
            continue

        if result is None:
            click.echo(f"  [{settings_key}] skipped — no comparison result")
            skipped.append(settings_key)
            continue

        json_path = pipeline_result.get("comparison_path")
        click.echo(f"  [{settings_key}] saved -> {json_path}")
        succeeded.append(settings_key)

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
        from polyzymd.compare.plotter import ComparisonPlotter

        try:
            plotter = ComparisonPlotter(config)
            generated = plotter.plot_all()
            click.echo(f"Generated {len(generated)} plots.")
        except Exception as e:
            click.echo(f"Error generating plots: {e}", err=True)
            if debug:
                import traceback

                traceback.print_exc()

    # Exit with error code if any comparison failed
    if failed:
        sys.exit(1)
