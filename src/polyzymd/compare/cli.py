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

from polyzymd.compare.cli_utils import (
    common_compare_options,
    load_comparison_config,
    validate_and_report,
)
from polyzymd.compare.comparators.rmsf import RMSFComparator
from polyzymd.compare.config import (
    ComparisonConfig,
    generate_comparison_template,
)
from polyzymd.compare.formatters import format_result
from polyzymd.compare.settings import (
    BindingFreeEnergyAnalysisSettings,
    BindingFreeEnergyComparisonSettings,
    CatalyticTriadAnalysisSettings,
    ContactsAnalysisSettings,
    ContactsComparisonSettings,
    RMSFAnalysisSettings,
)

LOGGER = logging.getLogger("polyzymd.compare.cli")


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
    3. polyzymd compare rmsf             # Run comparison

    \b
    Example:
        polyzymd compare init -n polymer_study
        cd polymer_study
        vim comparison.yaml  # Add your conditions
        polyzymd compare rmsf --eq-time 10ns
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
      - results/: Directory for comparison result JSON files
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
        # Create directory structure
        project_dir.mkdir(parents=True)
        (project_dir / "results").mkdir()
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

    analysis_settings:
      contacts:
        compute_binding_preference: true
        enzyme_pdb_for_sasa: "structures/enzyme.pdb"

The enzyme PDB should be the reference structure (e.g., from PDB or your
prepared simulation input), NOT a trajectory frame.
"""
        (project_dir / "structures" / "README.md").write_text(structures_readme)

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
            "sections_configured": config.analysis_settings.get_enabled_analyses(),
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


@compare.command()
@common_compare_options
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

    config = load_comparison_config(config_file)

    # Get RMSF settings from analysis_settings
    rmsf_settings = config.analysis_settings.get("rmsf")
    if rmsf_settings is None:
        click.echo("Error: No 'rmsf' in analysis_settings section", err=True)
        click.echo("", err=True)
        click.echo("RMSF comparison requires an 'analysis_settings.rmsf' section:", err=True)
        click.echo("", err=True)
        click.echo("  analysis_settings:", err=True)
        click.echo("    rmsf:", err=True)
        click.echo('      selection: "protein and name CA"', err=True)
        click.echo('      reference_mode: "centroid"', err=True)
        click.echo("", err=True)
        click.echo("And a corresponding comparison_settings entry:", err=True)
        click.echo("", err=True)
        click.echo("  comparison_settings:", err=True)
        click.echo("    rmsf: {}", err=True)
        sys.exit(1)

    # Cast to concrete type for attribute access
    rmsf_settings = RMSFAnalysisSettings.model_validate(rmsf_settings.model_dump())

    validate_and_report(config)

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Conditions: {len(config.conditions)}")

    # Build effective settings (YAML + optional overrides)
    equilibration = eq_time or config.defaults.equilibration_time
    effective_selection = selection if override and selection else rmsf_settings.selection
    effective_ref_mode = (
        reference_mode if override and reference_mode else rmsf_settings.reference_mode
    )
    effective_ref_frame = (
        reference_frame if override and reference_frame else rmsf_settings.reference_frame
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
            analysis_settings=rmsf_settings,
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
@common_compare_options
@click.option(
    "--recompute",
    is_flag=True,
    help="Force recompute triad analysis even if cached results exist.",
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
    from polyzymd.compare.comparators.triad import TriadComparator
    from polyzymd.compare.triad_formatters import format_triad_result

    # Set up logging with colored output
    setup_logging(quiet=quiet, debug=debug)

    config = load_comparison_config(config_file)

    # Check for catalytic_triad section in analysis_settings
    triad_settings = config.analysis_settings.get("catalytic_triad")
    if triad_settings is None:
        click.echo("Error: No 'catalytic_triad' in analysis_settings section", err=True)
        click.echo("", err=True)
        click.echo("Add a catalytic_triad section to your comparison.yaml:", err=True)
        click.echo("", err=True)
        click.echo("  analysis_settings:", err=True)
        click.echo("    catalytic_triad:", err=True)
        click.echo('      name: "enzyme_triad"', err=True)
        click.echo("      threshold: 3.5", err=True)
        click.echo("      pairs:", err=True)
        click.echo('        - label: "Asp-His"', err=True)
        click.echo('          selection_a: "resid 133 and name OD1"', err=True)
        click.echo('          selection_b: "resid 156 and name ND1"', err=True)
        click.echo("", err=True)
        click.echo("And a corresponding comparison_settings entry:", err=True)
        click.echo("", err=True)
        click.echo("  comparison_settings:", err=True)
        click.echo("    catalytic_triad: {}", err=True)
        sys.exit(1)

    # Cast to concrete type for attribute access
    triad_settings = CatalyticTriadAnalysisSettings.model_validate(triad_settings.model_dump())

    validate_and_report(config)

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Triad: {triad_settings.name}")
    click.echo(f"Pairs: {', '.join(triad_settings.get_pair_labels())}")
    click.echo(f"Conditions: {len(config.conditions)}")

    # Apply overrides
    equilibration = eq_time or config.defaults.equilibration_time

    click.echo(f"Equilibration: {equilibration}")
    click.echo(f"Threshold: {triad_settings.threshold} A")
    click.echo()

    # Run comparison
    try:
        comparator = TriadComparator(
            config=config,
            analysis_settings=triad_settings,
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
@common_compare_options
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
    from polyzymd.compare.comparators.contacts import ContactsComparator
    from polyzymd.compare.contacts_formatters import format_contacts_result

    # Set up logging with colored output
    setup_logging(quiet=quiet, debug=debug)

    config = load_comparison_config(config_file)

    validate_and_report(config)

    # Get contacts settings from analysis_settings
    contacts_analysis = config.analysis_settings.get("contacts")
    if contacts_analysis is None:
        click.echo("Error: No 'contacts' in analysis_settings section", err=True)
        click.echo("", err=True)
        click.echo(
            "Contacts comparison requires an 'analysis_settings.contacts' section:", err=True
        )
        click.echo("", err=True)
        click.echo("  analysis_settings:", err=True)
        click.echo("    contacts:", err=True)
        click.echo('      polymer_selection: "chainID C"', err=True)
        click.echo("      cutoff: 4.5", err=True)
        click.echo("", err=True)
        click.echo("And a corresponding comparison_settings entry:", err=True)
        click.echo("", err=True)
        click.echo("  comparison_settings:", err=True)
        click.echo("    contacts:", err=True)
        click.echo("      fdr_alpha: 0.05", err=True)
        sys.exit(1)

    # Cast to concrete types
    contacts_analysis = ContactsAnalysisSettings.model_validate(contacts_analysis.model_dump())

    # Get comparison settings (optional, defaults will be used if not present)
    contacts_comparison_raw = config.comparison_settings.get("contacts")
    if contacts_comparison_raw is not None:
        contacts_comparison = ContactsComparisonSettings.model_validate(
            contacts_comparison_raw.model_dump()
        )
    else:
        contacts_comparison = ContactsComparisonSettings()

    # Apply CLI overrides to analysis settings
    if polymer_selection:
        contacts_analysis = ContactsAnalysisSettings(
            **{**contacts_analysis.model_dump(), "polymer_selection": polymer_selection}
        )
    if cutoff is not None:
        contacts_analysis = ContactsAnalysisSettings(
            **{**contacts_analysis.model_dump(), "cutoff": cutoff}
        )

    # Apply CLI overrides to comparison settings
    if fdr_alpha is not None:
        contacts_comparison = ContactsComparisonSettings(
            **{**contacts_comparison.model_dump(), "fdr_alpha": fdr_alpha}
        )

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Conditions: {len(config.conditions)}")

    # Apply overrides
    equilibration = eq_time or config.defaults.equilibration_time

    click.echo(f"Equilibration: {equilibration}")
    click.echo(f"Polymer selection: {contacts_analysis.polymer_selection}")
    click.echo(f"Contact cutoff: {contacts_analysis.cutoff} A")
    click.echo(f"FDR alpha: {contacts_comparison.fdr_alpha}")
    click.echo(f"Compute binding preference: {contacts_analysis.compute_binding_preference}")
    if contacts_analysis.compute_binding_preference:
        click.echo(f"  Surface threshold: {contacts_analysis.surface_exposure_threshold}")
        click.echo(f"  Enzyme PDB: {contacts_analysis.enzyme_pdb_for_sasa}")
    if config.control:
        click.echo(f"Control: {config.control}")
    click.echo()

    # Run comparison
    try:
        comparator = ContactsComparator(
            config=config,
            analysis_settings=contacts_analysis,
            comparison_settings=contacts_comparison,
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


@compare.command()
@common_compare_options
@click.option(
    "--exposure-threshold",
    default=None,
    type=float,
    help="Override exposure threshold (fraction SASA; default: 0.20).",
)
@click.option(
    "--polymer-resnames",
    default=None,
    help="Override polymer residue names as comma-separated list (e.g., 'SBM,EGM').",
)
@click.option(
    "--recompute-sasa",
    is_flag=True,
    help="Force recompute SASA even if cached results exist.",
)
@click.option(
    "--recompute-exposure",
    is_flag=True,
    help="Force recompute exposure dynamics even if cached results exist.",
)
def exposure(
    config_file: Path,
    eq_time: Optional[str],
    exposure_threshold: Optional[float],
    polymer_resnames: Optional[str],
    recompute_sasa: bool,
    recompute_exposure: bool,
    output_format: str,
    output_path: Optional[Path],
    quiet: bool,
    debug: bool,
):
    """Compare chaperone-like polymer activity across conditions.

    Analyzes exposure dynamics: classifies protein residues as stably
    exposed, stably buried, or transiently exposed, then detects chaperone
    events (buried -> exposed -> polymer contact -> re-buried) across all
    conditions.

    Requires contacts analysis to have been run first for each condition
    (contacts_rep{N}.json must exist in analysis/contacts/).

    Optionally define an exposure section in comparison.yaml for custom settings:
    - exposure_threshold: fraction SASA defining 'exposed' (default: 0.20)
    - polymer_resnames: list of polymer residue names for enrichment
    - transient_lower / transient_upper: thresholds for classifying residues

    \b
    Example:
        polyzymd compare exposure
        polyzymd compare exposure --eq-time 10ns --format markdown
        polyzymd compare exposure --recompute-sasa --recompute-exposure
        polyzymd compare exposure --exposure-threshold 0.25 -o report.md
    """
    from polyzymd.analysis.core.logging_utils import setup_logging
    from polyzymd.compare.comparators.exposure import ExposureDynamicsComparator
    from polyzymd.compare.exposure_formatters import format_exposure_result
    from polyzymd.compare.settings import ExposureAnalysisSettings, ExposureComparisonSettings

    # Set up logging with colored output
    setup_logging(quiet=quiet, debug=debug)

    config = load_comparison_config(config_file)

    validate_and_report(config)

    # Get exposure settings from analysis_settings (optional — use defaults if absent)
    exposure_analysis_raw = config.analysis_settings.get("exposure")
    if exposure_analysis_raw is not None:
        exposure_analysis = ExposureAnalysisSettings.model_validate(
            exposure_analysis_raw.model_dump()
        )
    else:
        exposure_analysis = ExposureAnalysisSettings()

    # Get comparison settings (optional)
    exposure_comparison_raw = config.comparison_settings.get("exposure")
    if exposure_comparison_raw is not None:
        exposure_comparison = ExposureComparisonSettings.model_validate(
            exposure_comparison_raw.model_dump()
        )
    else:
        exposure_comparison = ExposureComparisonSettings()

    # Apply CLI overrides to analysis settings
    overrides: dict = {}
    if exposure_threshold is not None:
        overrides["exposure_threshold"] = exposure_threshold
    if polymer_resnames is not None:
        overrides["polymer_resnames"] = [r.strip() for r in polymer_resnames.split(",")]
    if overrides:
        exposure_analysis = ExposureAnalysisSettings(
            **{**exposure_analysis.model_dump(), **overrides}
        )

    # Recompute flags: --recompute-sasa or --recompute-exposure both force a full recompute
    recompute = recompute_sasa or recompute_exposure

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Conditions: {len(config.conditions)}")

    equilibration = eq_time or config.defaults.equilibration_time
    click.echo(f"Equilibration: {equilibration}")
    click.echo(f"Exposure threshold: {exposure_analysis.exposure_threshold}")
    click.echo(f"Probe radius: {exposure_analysis.probe_radius_nm} nm")
    if exposure_analysis.polymer_resnames:
        click.echo(f"Polymer residue names: {', '.join(exposure_analysis.polymer_resnames)}")
    if config.control:
        click.echo(f"Control: {config.control}")
    if recompute_sasa:
        click.echo("Recompute SASA: yes")
    if recompute_exposure:
        click.echo("Recompute exposure: yes")
    click.echo()

    # Run comparison
    try:
        comparator = ExposureDynamicsComparator(
            config=config,
            analysis_settings=exposure_analysis,
            comparison_settings=exposure_comparison,
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
            f"Note: {len(result.excluded_conditions)} condition(s) auto-excluded (no polymer): "
            f"{', '.join(result.excluded_conditions)}",
            fg="yellow",
        )
        click.echo()

    # Save JSON result
    results_dir = config_file.parent / "results"
    results_dir.mkdir(exist_ok=True)
    json_path = results_dir / f"exposure_comparison_{config.name.replace(' ', '_')}.json"
    result.save(json_path)
    click.echo(f"Saved result: {json_path}")
    click.echo()

    # Format and display output
    formatted = format_exposure_result(
        result,
        format=output_format,
    )
    click.echo(formatted)

    # Save formatted output if requested
    if output_path:
        output_path = Path(output_path)
        output_path.write_text(formatted)
        click.echo(f"Saved output: {output_path}")


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
    """Run a comparison using the registry pattern.

    This is a generic command that can run any registered comparison type.
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
    from polyzymd.analysis.core.logging_utils import setup_logging
    from polyzymd.compare.core.registry import ComparatorRegistry

    # Handle --list flag
    if list_types:
        available = ComparatorRegistry.list_available()
        click.echo("Available comparison types:")
        for comp_type in available:
            comparator_cls = ComparatorRegistry.get(comp_type)
            click.echo(f"  - {comp_type}: {comparator_cls.__name__}")
        return

    # Require comparison_type if not listing
    if comparison_type is None:
        click.echo("Error: Missing argument 'COMPARISON_TYPE'.", err=True)
        click.echo("Use 'polyzymd compare run --list' to see available types.", err=True)
        sys.exit(1)

    # Set up logging
    setup_logging(quiet=quiet, debug=debug)

    # Check if comparison type is registered
    available = ComparatorRegistry.list_available()
    if comparison_type not in available:
        click.echo(f"Error: Unknown comparison type '{comparison_type}'", err=True)
        click.echo(f"Available types: {', '.join(available)}", err=True)
        click.echo("", err=True)
        click.echo("Use 'polyzymd compare run --list' to see all available types.", err=True)
        sys.exit(1)

    config = load_comparison_config(config_file)

    validate_and_report(config)

    # Get the comparator class from registry
    comparator_cls = ComparatorRegistry.get(comparison_type)

    # Get analysis settings for this comparison type
    # Map comparison_type to analysis_settings key
    settings_key_map = {
        "rmsf": "rmsf",
        "triad": "catalytic_triad",
        "contacts": "contacts",
        "distances": "distances",
        "exposure": "exposure",
        "binding_free_energy": "binding_free_energy",
    }
    settings_key = settings_key_map.get(comparison_type, comparison_type)

    analysis_settings = config.analysis_settings.get(settings_key)
    if analysis_settings is None:
        click.echo(f"Error: No '{settings_key}' in analysis_settings section", err=True)
        click.echo("", err=True)
        click.echo(
            f"Add an analysis_settings.{settings_key} section to your comparison.yaml", err=True
        )
        sys.exit(1)

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Type: {comparison_type}")
    click.echo(f"Conditions: {len(config.conditions)}")

    # Build equilibration
    equilibration = eq_time or config.defaults.equilibration_time
    click.echo(f"Equilibration: {equilibration}")
    click.echo()

    # Create comparator and run
    try:
        comparator = comparator_cls(
            config=config,
            analysis_settings=analysis_settings,
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
    results_dir = config_file.parent / "results"
    results_dir.mkdir(exist_ok=True)
    json_path = results_dir / f"{comparison_type}_comparison_{config.name.replace(' ', '_')}.json"
    result.save(json_path)
    click.echo(f"Saved result: {json_path}")
    click.echo()

    # Format output based on comparison type
    try:
        if comparison_type == "rmsf":
            formatted = format_result(result, format=output_format)
        elif comparison_type == "triad":
            from polyzymd.compare.triad_formatters import format_triad_result

            formatted = format_triad_result(result, format=output_format)
        elif comparison_type == "contacts":
            from polyzymd.compare.contacts_formatters import format_contacts_result

            formatted = format_contacts_result(result, format=output_format)
        elif comparison_type == "distances":
            from polyzymd.compare.distances_formatters import format_distances_result

            formatted = format_distances_result(result, format=output_format)
        elif comparison_type == "exposure":
            from polyzymd.compare.exposure_formatters import format_exposure_result

            formatted = format_exposure_result(result, format=output_format)
        elif comparison_type == "binding_free_energy":
            from polyzymd.compare.binding_free_energy_formatters import format_bfe_result

            formatted = format_bfe_result(result, format=output_format)
        else:
            # Generic JSON output for unknown types
            formatted = result.model_dump_json(indent=2)
    except Exception as e:
        click.echo(f"Warning: Could not format result: {e}", err=True)
        formatted = result.model_dump_json(indent=2)

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
      - bfe_heatmap: ΔΔG heatmap (AA groups × conditions, per polymer type)
      - bfe_bars: ΔΔG grouped bar chart with ±k_BT reference lines

    \b
    Examples:
        polyzymd compare plot-all -f comparison.yaml
        polyzymd compare plot-all -f comparison.yaml -a catalytic_triad
        polyzymd compare plot-all -f comparison.yaml -p triad_kde_panel
        polyzymd compare plot-all --list-available
    """
    from polyzymd.compare.plotter import ComparisonPlotter, PlotterRegistry

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
            click.echo(f"  - {ptype}")
        click.echo()
        click.echo("Available plots for enabled analyses:")
        available = plotter.list_available_plots()
        for atype, ptypes in available.items():
            click.echo(f"  {atype}:")
            for pt in ptypes:
                click.echo(f"    - {pt}")
        return

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


@compare.command(name="binding-free-energy")
@common_compare_options
@click.option(
    "--units",
    type=click.Choice(["kcal/mol", "kJ/mol"]),
    default=None,
    help="Override energy units (default: kcal/mol).",
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
    help="Force recompute binding preference analysis even if cached results exist.",
)
def binding_free_energy(
    config_file: Path,
    eq_time: Optional[str],
    units: Optional[str],
    fdr_alpha: Optional[float],
    recompute: bool,
    output_format: str,
    output_path: Optional[Path],
    quiet: bool,
    debug: bool,
):
    """Compute ΔΔG (Gibbs free energy differences) from binding preference data.

    Converts existing binding preference probability data into physically grounded
    Gibbs free energy differences via Boltzmann inversion:

    \b
        ΔΔG = -k_B·T · ln(contact_share / expected_share)

    This answers: which residue groups does each polymer condition preferentially
    contact, and by how much in real energy units?

    Requires a 'binding_free_energy' section in analysis_settings (comparison.yaml),
    with contacts analysis already cached for each condition.

    \b
    Example:
        polyzymd compare binding-free-energy
        polyzymd compare binding-free-energy --units kJ/mol --format markdown
        polyzymd compare binding-free-energy --fdr-alpha 0.01 -o bfe_report.md
    """
    from polyzymd.analysis.core.logging_utils import setup_logging
    from polyzymd.compare.binding_free_energy_formatters import format_bfe_result
    from polyzymd.compare.comparators.binding_free_energy import BindingFreeEnergyComparator

    setup_logging(quiet=quiet, debug=debug)

    config = load_comparison_config(config_file)

    validate_and_report(config)

    # Get binding_free_energy settings from analysis_settings
    bfe_analysis_raw = config.analysis_settings.get("binding_free_energy")
    if bfe_analysis_raw is None:
        click.echo("Error: No 'binding_free_energy' in analysis_settings section", err=True)
        click.echo("", err=True)
        click.echo(
            "Binding free energy comparison requires an "
            "'analysis_settings.binding_free_energy' section:",
            err=True,
        )
        click.echo("", err=True)
        click.echo("  analysis_settings:", err=True)
        click.echo("    binding_free_energy:", err=True)
        click.echo("      units: kcal/mol", err=True)
        click.echo("      surface_exposure_threshold: 0.2", err=True)
        click.echo("", err=True)
        click.echo("And a corresponding comparison_settings entry:", err=True)
        click.echo("", err=True)
        click.echo("  comparison_settings:", err=True)
        click.echo("    binding_free_energy:", err=True)
        click.echo("      fdr_alpha: 0.05", err=True)
        sys.exit(1)

    bfe_analysis = BindingFreeEnergyAnalysisSettings.model_validate(bfe_analysis_raw.model_dump())

    bfe_comparison_raw = config.comparison_settings.get("binding_free_energy")
    if bfe_comparison_raw is not None:
        bfe_comparison = BindingFreeEnergyComparisonSettings.model_validate(
            bfe_comparison_raw.model_dump()
        )
    else:
        bfe_comparison = BindingFreeEnergyComparisonSettings()

    # Apply CLI overrides
    if units is not None:
        bfe_analysis = BindingFreeEnergyAnalysisSettings(
            **{**bfe_analysis.model_dump(), "units": units}
        )
    if fdr_alpha is not None:
        bfe_comparison = BindingFreeEnergyComparisonSettings(
            **{**bfe_comparison.model_dump(), "fdr_alpha": fdr_alpha}
        )

    equilibration = eq_time or config.defaults.equilibration_time

    click.echo(f"Comparison: {config.name}")
    click.echo(f"Conditions: {len(config.conditions)}")
    click.echo(f"Equilibration: {equilibration}")
    click.echo(f"Units: {bfe_analysis.units}")
    click.echo(f"Surface exposure threshold: {bfe_analysis.surface_exposure_threshold}")
    click.echo(f"FDR alpha: {bfe_comparison.fdr_alpha}")
    if config.control:
        click.echo(f"Control: {config.control}")
    click.echo()

    # Run comparison
    try:
        comparator = BindingFreeEnergyComparator(
            config=config,
            analysis_settings=bfe_analysis,
            comparison_settings=bfe_comparison,
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
    results_dir = config_file.parent / "results"
    results_dir.mkdir(exist_ok=True)
    json_path = results_dir / f"binding_free_energy_comparison_{config.name.replace(' ', '_')}.json"
    result.save(json_path)
    click.echo(f"Saved result: {json_path}")
    click.echo()

    # Format and display output
    formatted = format_bfe_result(result, format=output_format)
    click.echo(formatted)

    # Save formatted output if requested
    if output_path:
        output_path = Path(output_path)
        output_path.write_text(formatted)
        click.echo(f"Saved output: {output_path}")
