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


@compare.command()
@click.argument(
    "result_file",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(path_type=Path),
    default="figures",
    help="Output directory for plots. Default: figures/",
)
@click.option(
    "--format",
    "img_format",
    type=click.Choice(["png", "pdf", "svg"]),
    default="png",
    help="Image format for saved plots.",
)
@click.option(
    "--dpi",
    default=150,
    help="Resolution for PNG output (default: 150).",
)
@click.option(
    "--summary/--no-summary",
    default=True,
    help="Generate combined summary panel (default: yes).",
)
@click.option(
    "--show/--no-show",
    default=False,
    help="Display plots interactively after saving.",
)
def plot(
    result_file: Path,
    output_dir: Path,
    img_format: str,
    dpi: int,
    summary: bool,
    show: bool,
):
    """Generate comparison plots from saved results.

    Creates publication-ready figures from a comparison result JSON file.

    \b
    Generated plots:
      - rmsf_comparison: Bar chart of RMSF by condition
      - percent_change: Bar chart of % change vs control
      - effect_sizes: Forest plot of Cohen's d values
      - summary_panel: Combined multi-panel figure (if --summary)

    \b
    Color coding:
      - Green: Significant improvement (p<0.05)
      - Blue: Large effect (d>0.8) but not significant
      - Gray: Control or no effect
      - Red: Worse than control

    \b
    Example:
        polyzymd compare plot results/rmsf_comparison_my_study.json
        polyzymd compare plot results/rmsf_comparison_my_study.json -o figures/ --dpi 300
        polyzymd compare plot results/rmsf_comparison_my_study.json --format pdf --show
    """
    from polyzymd.compare.results import ComparisonResult
    from polyzymd.compare.plotting import (
        plot_rmsf_comparison,
        plot_percent_change,
        plot_effect_sizes,
        plot_summary_panel,
    )

    # Load result
    try:
        result = ComparisonResult.load(result_file)
    except Exception as e:
        click.echo(f"Error loading result: {e}", err=True)
        sys.exit(1)

    click.echo(f"Loaded comparison: {result.name}")
    click.echo(f"Conditions: {len(result.conditions)}")
    click.echo()

    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    generated = []

    # Plot 1: RMSF comparison
    rmsf_path = output_dir / f"rmsf_comparison.{img_format}"
    click.echo(f"Generating RMSF comparison plot...")
    plot_rmsf_comparison(result, save_path=rmsf_path, dpi=dpi)
    generated.append(rmsf_path)

    # Plot 2: Percent change (requires control)
    if result.control_label and result.pairwise_comparisons:
        pct_path = output_dir / f"percent_change.{img_format}"
        click.echo(f"Generating percent change plot...")
        plot_percent_change(result, save_path=pct_path, dpi=dpi)
        generated.append(pct_path)

        # Plot 3: Effect sizes
        effect_path = output_dir / f"effect_sizes.{img_format}"
        click.echo(f"Generating effect sizes plot...")
        plot_effect_sizes(result, save_path=effect_path, dpi=dpi)
        generated.append(effect_path)
    else:
        click.echo("Skipping percent change and effect size plots (no control condition)")

    # Plot 4: Summary panel
    if summary and result.control_label and result.pairwise_comparisons:
        summary_path = output_dir / f"summary_panel.{img_format}"
        click.echo(f"Generating summary panel...")
        plot_summary_panel(result, save_path=summary_path, dpi=dpi)
        generated.append(summary_path)

    click.echo()
    click.echo("Generated plots:")
    for path in generated:
        click.echo(f"  - {path}")

    if show:
        import matplotlib.pyplot as plt

        plt.show()
