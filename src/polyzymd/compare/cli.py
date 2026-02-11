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

from polyzymd.compare.config import ComparisonConfig, generate_comparison_template
from polyzymd.compare.comparator import RMSFComparator
from polyzymd.compare.formatters import format_result

LOGGER = logging.getLogger("polyzymd.compare.cli")


@click.group()
def compare():
    """Compare analysis results across simulation conditions.

    The compare module allows you to statistically compare RMSF and other
    metrics across multiple simulation conditions (e.g., different polymer
    compositions, temperatures, etc.).

    \b
    Workflow:
    1. polyzymd compare init <name>    # Create project with template
    2. Edit comparison.yaml            # Add your conditions
    3. polyzymd compare rmsf           # Run comparison

    \b
    Example:
        polyzymd compare init polymer_study
        cd polymer_study
        vim comparison.yaml  # Add your conditions
        polyzymd compare rmsf --eq-time 10ns
    """
    pass


@compare.command()
@click.argument("name")
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
      - results/: Directory for comparison result JSON files
      - figures/: Directory for comparison plots

    \b
    Example:
        polyzymd compare init polymer_stability_study
        polyzymd compare init my_study --eq-time 20ns
        polyzymd compare init my_study -o /path/to/projects
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
        # Create directory structure
        project_dir.mkdir(parents=True)
        (project_dir / "results").mkdir()
        (project_dir / "figures").mkdir()

        # Generate and write template
        template_content = generate_comparison_template(name, eq_time)
        config_path = project_dir / "comparison.yaml"
        config_path.write_text(template_content)

        click.echo(f"Created comparison project: {project_dir}")
        click.echo()
        click.echo("Next steps:")
        click.echo(f"  1. Edit {config_path.relative_to(Path.cwd())}")
        click.echo("     - Add your simulation conditions")
        click.echo("     - Set paths to config.yaml files")
        click.echo("     - Specify replicate numbers")
        click.echo()
        click.echo(f"  2. cd {project_dir.relative_to(Path.cwd())}")
        click.echo("  3. polyzymd compare rmsf")
        click.echo()

    except Exception as e:
        click.echo(f"Error creating project: {e}", err=True)
        sys.exit(1)


@compare.command()
@click.option(
    "-f",
    "--file",
    "config_file",
    type=click.Path(exists=True, path_type=Path),
    default="comparison.yaml",
    help="Path to comparison.yaml config file.",
)
@click.option(
    "--eq-time",
    default=None,
    help="Override equilibration time (e.g., '10ns', '5000ps').",
)
@click.option(
    "--selection",
    default=None,
    help="Override atom selection (e.g., 'protein and name CA').",
)
@click.option(
    "--recompute",
    is_flag=True,
    help="Force recompute RMSF even if cached results exist.",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["table", "markdown", "json"]),
    default="table",
    help="Output format: table (default), markdown, or json.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Save output to file. Also saves JSON result to results/ directory.",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Show verbose logging output.",
)
def rmsf(
    config_file: Path,
    eq_time: Optional[str],
    selection: Optional[str],
    recompute: bool,
    output_format: str,
    output_path: Optional[Path],
    verbose: bool,
):
    """Compare RMSF across conditions defined in comparison.yaml.

    Loads RMSF results for each condition (computing them if necessary),
    then performs statistical comparisons including t-tests, effect sizes,
    and ANOVA.

    \b
    Example:
        polyzymd compare rmsf
        polyzymd compare rmsf --eq-time 10ns --format markdown
        polyzymd compare rmsf -f my_comparison.yaml -o report.md
    """
    # Set up logging
    if verbose:
        logging.basicConfig(level=logging.INFO, format="%(message)s")
    else:
        logging.basicConfig(level=logging.WARNING, format="%(message)s")

    config_file = Path(config_file).resolve()

    # Load and validate config
    click.echo(f"Loading config: {config_file}")
    try:
        config = ComparisonConfig.from_yaml(config_file)
    except FileNotFoundError as e:
        click.echo(f"Error: {e}", err=True)
        click.echo("Run 'polyzymd compare init <name>' to create a comparison project.", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Error loading config: {e}", err=True)
        sys.exit(1)

    # Validate config
    errors = config.validate_config()
    if errors:
        click.echo("Configuration errors:", err=True)
        for error in errors:
            click.echo(f"  - {error}", err=True)
        sys.exit(1)

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Conditions: {len(config.conditions)}")

    # Apply overrides
    equilibration = eq_time or config.defaults.equilibration_time
    atom_selection = selection or config.defaults.selection

    click.echo(f"Equilibration: {equilibration}")
    click.echo(f"Selection: {atom_selection}")
    click.echo()

    # Run comparison
    try:
        comparator = RMSFComparator(
            config=config,
            equilibration=equilibration,
            selection=atom_selection,
        )
        result = comparator.compare(recompute=recompute)
    except Exception as e:
        click.echo(f"Error during comparison: {e}", err=True)
        if verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)

    # Format and display output
    formatted = format_result(result, format=output_format)
    click.echo(formatted)

    # Save results
    results_dir = config_file.parent / "results"
    if results_dir.exists():
        # Always save JSON result
        json_path = results_dir / f"rmsf_comparison_{result.name}.json"
        result.save(json_path)
        click.echo(f"Results saved: {json_path}")

    # Save formatted output if requested
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(formatted)
        click.echo(f"Output saved: {output_path}")


@compare.command()
@click.option(
    "-f",
    "--file",
    "config_file",
    type=click.Path(exists=True, path_type=Path),
    default="comparison.yaml",
    help="Path to comparison.yaml config file.",
)
def validate(config_file: Path):
    """Validate a comparison.yaml configuration file.

    Checks that:
      - At least 2 conditions are defined
      - No duplicate labels
      - Control condition exists (if specified)
      - Config paths exist for each condition

    \b
    Example:
        polyzymd compare validate
        polyzymd compare validate -f my_comparison.yaml
    """
    config_file = Path(config_file).resolve()

    try:
        config = ComparisonConfig.from_yaml(config_file)
    except FileNotFoundError:
        click.echo(f"Error: Config file not found: {config_file}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Error parsing config: {e}", err=True)
        sys.exit(1)

    click.echo(f"Validating: {config_file}")
    click.echo()

    errors = config.validate_config()

    if errors:
        click.echo("Validation FAILED:", err=True)
        for error in errors:
            click.echo(f"  - {error}", err=True)
        sys.exit(1)
    else:
        click.echo("Validation PASSED")
        click.echo()
        click.echo(f"  Name: {config.name}")
        click.echo(f"  Conditions: {len(config.conditions)}")
        for cond in config.conditions:
            marker = " (control)" if cond.label == config.control else ""
            click.echo(f"    - {cond.label}{marker}: {len(cond.replicates)} replicates")
        click.echo(f"  Equilibration: {config.defaults.equilibration_time}")
        click.echo(f"  Selection: {config.defaults.selection}")
        click.echo()


@compare.command()
@click.argument(
    "result_file",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["table", "markdown", "json"]),
    default="table",
    help="Output format.",
)
def show(result_file: Path, output_format: str):
    """Display a saved comparison result.

    \b
    Example:
        polyzymd compare show results/rmsf_comparison_my_study.json
        polyzymd compare show results/rmsf_comparison_my_study.json --format markdown
    """
    from polyzymd.compare.results import ComparisonResult

    try:
        result = ComparisonResult.load(result_file)
    except Exception as e:
        click.echo(f"Error loading result: {e}", err=True)
        sys.exit(1)

    formatted = format_result(result, format=output_format)
    click.echo(formatted)
