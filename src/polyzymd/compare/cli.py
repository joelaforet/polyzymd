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

from polyzymd.compare.comparator import RMSFComparator
from polyzymd.compare.config import (
    ComparisonConfig,
    RMSFComparisonConfig,
    generate_comparison_template,
)
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
        click.echo("     - Add your simulation conditions (paths to config.yaml files)")
        click.echo("     - Define catalytic_triad for active site analysis")
        click.echo()
        click.echo(f"  2. cd {project_dir.relative_to(Path.cwd())}")
        click.echo("  3. Run comparisons:")
        click.echo("     polyzymd compare rmsf      # Compare flexibility")
        click.echo("     polyzymd compare triad     # Compare triad geometry")
        click.echo("     polyzymd compare contacts  # Compare polymer-protein contacts")
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
    "--override",
    is_flag=True,
    help="Enable CLI overrides for RMSF settings. Without this flag, YAML is the single source of truth.",
)
@click.option(
    "--selection",
    default=None,
    help="Override atom selection (requires --override flag).",
)
@click.option(
    "--reference-mode",
    default=None,
    type=click.Choice(["centroid", "average", "frame"]),
    help="Override reference mode (requires --override flag).",
)
@click.option(
    "--reference-frame",
    default=None,
    type=int,
    help="Override reference frame (requires --override flag, used with --reference-mode frame).",
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
def rmsf(
    config_file: Path,
    eq_time: Optional[str],
    override: bool,
    selection: Optional[str],
    reference_mode: Optional[str],
    reference_frame: Optional[int],
    recompute: bool,
    output_format: str,
    output_path: Optional[Path],
    quiet: bool,
    debug: bool,
):
    """Compare RMSF across conditions defined in comparison.yaml.

    Requires an `rmsf:` section in comparison.yaml. This is the single source
    of truth for RMSF settings. Use --override to allow CLI options to take
    precedence.

    \b
    Example:
        polyzymd compare rmsf
        polyzymd compare rmsf --eq-time 10ns --format markdown
        polyzymd compare rmsf --override --selection "protein and name CA"
        polyzymd compare rmsf -f my_comparison.yaml -o report.md
    """
    # Set up logging with colored output
    from polyzymd.analysis.core.logging_utils import setup_logging

    setup_logging(quiet=quiet, debug=debug)

    # Check if override options used without --override flag
    if not override and any([selection, reference_mode, reference_frame]):
        click.echo("Error: CLI overrides require the --override flag.", err=True)
        click.echo("", err=True)
        click.echo("Either:", err=True)
        click.echo("  1. Add --override to use CLI options", err=True)
        click.echo("  2. Set values in comparison.yaml (recommended)", err=True)
        sys.exit(1)

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

    # Check for rmsf: section (required)
    if config.rmsf is None:
        click.echo("Error: No 'rmsf:' section in comparison.yaml", err=True)
        click.echo("", err=True)
        click.echo(
            "RMSF comparison requires an 'rmsf:' section. Add to your comparison.yaml:", err=True
        )
        click.echo("", err=True)
        click.echo("  rmsf:", err=True)
        click.echo('    selection: "protein and name CA"', err=True)
        click.echo('    reference_mode: "centroid"', err=True)
        click.echo("", err=True)
        click.echo("For frame-based alignment:", err=True)
        click.echo("", err=True)
        click.echo("  rmsf:", err=True)
        click.echo('    selection: "protein and name CA"', err=True)
        click.echo('    reference_mode: "frame"', err=True)
        click.echo("    reference_frame: 500", err=True)
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

    # Build effective settings (YAML + optional overrides)
    equilibration = eq_time or config.defaults.equilibration_time
    effective_selection = selection if override and selection else config.rmsf.selection
    effective_ref_mode = (
        reference_mode if override and reference_mode else config.rmsf.reference_mode
    )
    effective_ref_frame = (
        reference_frame if override and reference_frame else config.rmsf.reference_frame
    )

    click.echo(f"Equilibration: {equilibration}")
    click.echo(f"Selection: {effective_selection}")
    click.echo(f"Reference mode: {effective_ref_mode}")
    if effective_ref_frame:
        click.echo(f"Reference frame: {effective_ref_frame}")
    click.echo()

    # Run comparison
    try:
        comparator = RMSFComparator(
            config=config,
            rmsf_config=config.rmsf,
            equilibration=equilibration,
            selection_override=selection if override else None,
            reference_mode_override=reference_mode if override else None,
            reference_frame_override=reference_frame if override else None,
        )
        result = comparator.compare(recompute=recompute)
    except Exception as e:
        click.echo(f"Error during comparison: {e}", err=True)
        if debug:
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
    from polyzymd.compare.plotting import (
        plot_effect_sizes,
        plot_percent_change,
        plot_rmsf_comparison,
        plot_summary_panel,
    )
    from polyzymd.compare.results import ComparisonResult

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
    click.echo("Generating RMSF comparison plot...")
    plot_rmsf_comparison(result, save_path=rmsf_path, dpi=dpi)
    generated.append(rmsf_path)

    # Plot 2: Percent change (requires control)
    if result.control_label and result.pairwise_comparisons:
        pct_path = output_dir / f"percent_change.{img_format}"
        click.echo("Generating percent change plot...")
        plot_percent_change(result, save_path=pct_path, dpi=dpi)
        generated.append(pct_path)

        # Plot 3: Effect sizes
        effect_path = output_dir / f"effect_sizes.{img_format}"
        click.echo("Generating effect sizes plot...")
        plot_effect_sizes(result, save_path=effect_path, dpi=dpi)
        generated.append(effect_path)
    else:
        click.echo("Skipping percent change and effect size plots (no control condition)")

    # Plot 4: Summary panel
    if summary and result.control_label and result.pairwise_comparisons:
        summary_path = output_dir / f"summary_panel.{img_format}"
        click.echo("Generating summary panel...")
        plot_summary_panel(result, save_path=summary_path, dpi=dpi)
        generated.append(summary_path)

    click.echo()
    click.echo("Generated plots:")
    for path in generated:
        click.echo(f"  - {path}")

    if show:
        import matplotlib.pyplot as plt

        plt.show()


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
    "--recompute",
    is_flag=True,
    help="Force recompute triad analysis even if cached results exist.",
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
def triad(
    config_file: Path,
    eq_time: Optional[str],
    recompute: bool,
    output_format: str,
    output_path: Optional[Path],
    quiet: bool,
    debug: bool,
):
    """Compare catalytic triad geometry across conditions.

    Analyzes the catalytic triad/active site defined in comparison.yaml,
    comparing the simultaneous contact fraction across all conditions.
    Higher contact fraction indicates better triad integrity.

    Requires a catalytic_triad section in comparison.yaml defining:
    - name: Name of the triad
    - threshold: Contact distance threshold (Angstroms)
    - pairs: List of distance pairs with selections

    \b
    Example:
        polyzymd compare triad
        polyzymd compare triad --eq-time 10ns --format markdown
        polyzymd compare triad -f my_comparison.yaml -o report.md
    """
    from polyzymd.analysis.core.logging_utils import setup_logging
    from polyzymd.compare.triad_comparator import TriadComparator
    from polyzymd.compare.triad_formatters import format_triad_result

    # Set up logging with colored output
    setup_logging(quiet=quiet, debug=debug)

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

    # Check for catalytic_triad section
    if config.catalytic_triad is None:
        click.echo("Error: No catalytic_triad section in comparison.yaml", err=True)
        click.echo("", err=True)
        click.echo("Add a catalytic_triad section to your comparison.yaml:", err=True)
        click.echo("", err=True)
        click.echo("  catalytic_triad:", err=True)
        click.echo('    name: "enzyme_triad"', err=True)
        click.echo("    threshold: 3.5", err=True)
        click.echo("    pairs:", err=True)
        click.echo('      - label: "Asp-His"', err=True)
        click.echo('        selection_a: "resid 133 and name OD1"', err=True)
        click.echo('        selection_b: "resid 156 and name ND1"', err=True)
        sys.exit(1)

    # Validate config
    errors = config.validate_config()
    if errors:
        click.echo("Configuration errors:", err=True)
        for error in errors:
            click.echo(f"  - {error}", err=True)
        sys.exit(1)

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Triad: {config.catalytic_triad.name}")
    click.echo(f"Pairs: {', '.join(config.catalytic_triad.get_pair_labels())}")
    click.echo(f"Conditions: {len(config.conditions)}")

    # Apply overrides
    equilibration = eq_time or config.defaults.equilibration_time

    click.echo(f"Equilibration: {equilibration}")
    click.echo(f"Threshold: {config.catalytic_triad.threshold} A")
    click.echo()

    # Run comparison
    try:
        comparator = TriadComparator(
            config=config,
            equilibration=equilibration,
        )
        result = comparator.compare(recompute=recompute)
    except Exception as e:
        click.echo(f"Error during comparison: {e}", err=True)
        if debug:
            import traceback

            traceback.print_exc()
        sys.exit(1)

    # Save JSON result
    if output_path or True:  # Always save JSON to results/
        results_dir = config_file.parent / "results"
        results_dir.mkdir(exist_ok=True)
        json_path = results_dir / f"triad_comparison_{config.name.replace(' ', '_')}.json"
        result.save(json_path)
        click.echo(f"Saved result: {json_path}")
        click.echo()

    # Format and display output
    formatted = format_triad_result(result, format=output_format)
    click.echo(formatted)

    # Save formatted output if requested
    if output_path:
        output_path = Path(output_path)
        output_path.write_text(formatted)
        click.echo(f"Saved output: {output_path}")


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
    "--polymer-selection",
    default=None,
    help="Override polymer selection (e.g., 'resname SBM EGM').",
)
@click.option(
    "--cutoff",
    default=None,
    type=float,
    help="Override contact cutoff distance (Angstroms).",
)
@click.option(
    "--fdr-alpha",
    default=None,
    type=float,
    help="Override FDR alpha for Benjamini-Hochberg correction (default: 0.05).",
)
@click.option(
    "--recompute",
    is_flag=True,
    help="Force recompute contacts analysis even if cached results exist.",
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
def contacts(
    config_file: Path,
    eq_time: Optional[str],
    polymer_selection: Optional[str],
    cutoff: Optional[float],
    fdr_alpha: Optional[float],
    recompute: bool,
    output_format: str,
    output_path: Optional[Path],
    quiet: bool,
    debug: bool,
):
    """Compare polymer-protein contacts across conditions.

    Analyzes polymer-protein contact statistics across all conditions,
    comparing aggregate metrics (coverage, mean contact fraction) between
    conditions with statistical tests and effect sizes (Cohen's d).

    Optionally define a contacts section in comparison.yaml for custom settings:
    - polymer_selection: MDAnalysis selection for polymer (default: "resname SBM EGM")
    - cutoff: Contact distance threshold (default: 4.5 A)

    Note: For per-residue contact-RMSF correlations (mechanistic insights),
    use 'polyzymd compare report' after running both RMSF and contacts analyses.

    \b
    Example:
        polyzymd compare contacts
        polyzymd compare contacts --eq-time 10ns --format markdown
        polyzymd compare contacts --polymer-selection "resname SBM" --cutoff 5.0
        polyzymd compare contacts -f my_comparison.yaml -o report.md
    """
    from polyzymd.analysis.core.logging_utils import setup_logging
    from polyzymd.compare.config import ContactsComparisonConfig
    from polyzymd.compare.contacts_comparator import ContactsComparator
    from polyzymd.compare.contacts_formatters import format_contacts_result

    # Set up logging with colored output
    setup_logging(quiet=quiet, debug=debug)

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

    # Build contacts config with overrides
    if config.contacts is not None:
        contacts_config = config.contacts
    else:
        contacts_config = ContactsComparisonConfig()

    # Apply CLI overrides
    if polymer_selection:
        contacts_config = ContactsComparisonConfig(
            **{**contacts_config.model_dump(), "polymer_selection": polymer_selection}
        )
    if cutoff is not None:
        contacts_config = ContactsComparisonConfig(
            **{**contacts_config.model_dump(), "cutoff": cutoff}
        )
    if fdr_alpha is not None:
        contacts_config = ContactsComparisonConfig(
            **{**contacts_config.model_dump(), "fdr_alpha": fdr_alpha}
        )

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Conditions: {len(config.conditions)}")

    # Apply overrides
    equilibration = eq_time or config.defaults.equilibration_time

    click.echo(f"Equilibration: {equilibration}")
    click.echo(f"Polymer selection: {contacts_config.polymer_selection}")
    click.echo(f"Contact cutoff: {contacts_config.cutoff} A")
    click.echo(f"FDR alpha: {contacts_config.fdr_alpha}")
    if config.control:
        click.echo(f"Control: {config.control}")
    click.echo()

    # Run comparison
    try:
        comparator = ContactsComparator(
            config=config,
            contacts_config=contacts_config,
            equilibration=equilibration,
        )
        result = comparator.compare(recompute=recompute)
    except Exception as e:
        click.echo(f"Error during comparison: {e}", err=True)
        if debug:
            import traceback

            traceback.print_exc()
        sys.exit(1)

    # Show auto-excluded conditions warning
    if result.excluded_conditions:
        click.secho(
            f"Note: {len(result.excluded_conditions)} condition(s) auto-excluded (no polymer atoms): "
            f"{', '.join(result.excluded_conditions)}",
            fg="yellow",
        )
        click.echo()

    # Save JSON result
    results_dir = config_file.parent / "results"
    results_dir.mkdir(exist_ok=True)
    json_path = results_dir / f"contacts_comparison_{config.name.replace(' ', '_')}.json"
    result.save(json_path)
    click.echo(f"Saved result: {json_path}")
    click.echo()

    # Format and display output
    formatted = format_contacts_result(
        result,
        format=output_format,
    )
    click.echo(formatted)

    # Save formatted output if requested
    if output_path:
        output_path = Path(output_path)
        output_path.write_text(formatted)
        click.echo(f"Saved output: {output_path}")
