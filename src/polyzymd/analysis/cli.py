"""Command-line interface for PolyzyMD analysis module.

This module provides CLI commands for trajectory analysis:
- `polyzymd analyze rmsf` - RMSF analysis
- `polyzymd analyze distances` - Distance analysis
- `polyzymd analyze triad` - Catalytic triad/active site analysis
- `polyzymd plot rmsf` - Plot RMSF results
- `polyzymd plot distances` - Plot distance results
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

import click

LOGGER = logging.getLogger("polyzymd.analysis")


# =============================================================================
# Helper Functions
# =============================================================================


def parse_replicates(rep_str: str) -> list[int]:
    """Parse replicate specification string.

    Examples
    --------
    >>> parse_replicates("1-5")
    [1, 2, 3, 4, 5]
    >>> parse_replicates("1,3,5")
    [1, 3, 5]
    >>> parse_replicates("1")
    [1]
    >>> parse_replicates("1-3,7,9-10")
    [1, 2, 3, 7, 9, 10]
    """
    result = set()
    for part in rep_str.split(","):
        part = part.strip()
        if "-" in part:
            start, end = part.split("-", 1)
            result.update(range(int(start.strip()), int(end.strip()) + 1))
        else:
            result.add(int(part))
    return sorted(result)


def parse_distance_pairs(pair_strs: tuple[str, ...]) -> list[tuple[str, str]]:
    """Parse distance pair specifications.

    Each pair should be formatted as "selection1 : selection2"

    Examples
    --------
    >>> parse_distance_pairs(("resid 77 and name OG : resid 133 and name NE2",))
    [("resid 77 and name OG", "resid 133 and name NE2")]
    """
    pairs = []
    for pair_str in pair_strs:
        if ":" not in pair_str:
            raise click.BadParameter(
                f"Invalid pair format: '{pair_str}'. Use 'selection1 : selection2' format."
            )
        sel1, sel2 = pair_str.split(":", 1)
        pairs.append((sel1.strip(), sel2.strip()))
    return pairs


def require_analysis_deps() -> None:
    """Check that analysis dependencies are installed."""
    try:
        import MDAnalysis  # noqa: F401
    except ImportError:
        click.echo(
            click.style(
                "Error: MDAnalysis is required for analysis.\n"
                "Install with: pip install polyzymd[analysis]",
                fg="red",
            ),
            err=True,
        )
        sys.exit(1)


def require_matplotlib() -> None:
    """Check that matplotlib is installed."""
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        click.echo(
            click.style(
                "Error: matplotlib is required for plotting.\n"
                "Install with: pip install polyzymd[analysis]",
                fg="red",
            ),
            err=True,
        )
        sys.exit(1)


# =============================================================================
# Analyze Command Group
# =============================================================================


@click.group()
def analyze() -> None:
    """Analyze MD trajectories.

    Run various analysis types on simulation trajectories, including RMSF
    (Root Mean Square Fluctuation) and inter-atomic distance calculations.

    Analysis results are saved as JSON files for later plotting and comparison.

    \b
    Workflow:
    1. cd your_project/        # Directory with config.yaml
    2. polyzymd analyze init   # Create analysis.yaml template
    3. Edit analysis.yaml      # Configure analyses
    4. polyzymd analyze run    # Run all enabled analyses

    \b
    Or run individual analyses:
        polyzymd analyze rmsf -c config.yaml -r 1-5 --eq-time 100ns
        polyzymd analyze distances -c config.yaml -r 1 --pair "sel1 : sel2"
    """
    pass


# -----------------------------------------------------------------------------
# Init Command - Create analysis.yaml scaffold
# -----------------------------------------------------------------------------


@analyze.command()
@click.option(
    "--eq-time",
    default="10ns",
    help="Default equilibration time (e.g., '10ns', '5000ps').",
)
def init(eq_time: str) -> None:
    """Initialize analysis configuration for this simulation.

    Must be run from a directory containing config.yaml.

    \b
    Creates:
      - analysis.yaml: Configuration for which analyses to run
      - analysis/: Directory structure for results
        ├── rmsf/
        ├── distances/
        └── triad/

    \b
    Example:
        cd my_project
        polyzymd analyze init
        polyzymd analyze init --eq-time 20ns
    """
    from polyzymd.analysis.config import generate_analysis_template

    cwd = Path.cwd()
    config_yaml = cwd / "config.yaml"

    # Validate we're in the right place
    if not config_yaml.exists():
        click.echo(
            click.style(
                "Error: config.yaml not found in current directory.",
                fg="red",
            ),
            err=True,
        )
        click.echo("Run this command from your simulation project directory.", err=True)
        sys.exit(1)

    analysis_yaml = cwd / "analysis.yaml"
    if analysis_yaml.exists():
        click.echo(
            click.style(
                f"Error: analysis.yaml already exists at {analysis_yaml}",
                fg="red",
            ),
            err=True,
        )
        click.echo("Delete it first if you want to regenerate.", err=True)
        sys.exit(1)

    # Create directory structure
    analysis_dir = cwd / "analysis"
    analysis_dir.mkdir(exist_ok=True)
    (analysis_dir / "rmsf").mkdir(exist_ok=True)
    (analysis_dir / "distances").mkdir(exist_ok=True)
    (analysis_dir / "catalytic_triad").mkdir(exist_ok=True)

    # Generate template
    template = generate_analysis_template(eq_time)
    analysis_yaml.write_text(template)

    click.echo(f"Created: {analysis_yaml.name}")
    click.echo(f"Created: {analysis_dir.name}/")
    click.echo()
    click.echo("Next steps:")
    click.echo("  1. Edit analysis.yaml to configure your analyses")
    click.echo("  2. Run: polyzymd analyze run")


# -----------------------------------------------------------------------------
# Validate Command - Check analysis.yaml for errors
# -----------------------------------------------------------------------------


@analyze.command()
@click.option(
    "-f",
    "--file",
    "config_file",
    type=click.Path(path_type=Path),
    default="analysis.yaml",
    help="Path to analysis.yaml config file.",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["table", "json"]),
    default="table",
    help="Output format: table (default) or json.",
)
def validate(config_file: Path, output_format: str):
    """Validate an analysis.yaml configuration file.

    Checks the configuration for errors without running any analyses.
    Useful for catching configuration problems before running expensive
    computations.

    \b
    Validates:
      - YAML syntax and structure
      - Required fields present
      - Replicates list is non-empty
      - Distance pairs defined if distances enabled
      - Triad pairs defined if catalytic_triad enabled
      - Contact selections defined if contacts enabled

    \b
    Example:
        polyzymd analyze validate
        polyzymd analyze validate -f my_analysis.yaml
        polyzymd analyze validate --format json
    """
    import json as json_module

    import yaml

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
        from polyzymd.analysis.config import AnalysisConfig

        config = AnalysisConfig.from_yaml(config_file)

        # Run validation
        errors = config.validate_config()
        result["errors"] = errors
        result["valid"] = len(errors) == 0

        # Build summary
        result["summary"] = {
            "replicates": config.replicates,
            "equilibration_time": config.defaults.equilibration_time,
            "enabled_analyses": config.get_enabled_analyses(),
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
            replicates = summary.get("replicates", [])
            if replicates:
                click.echo(f"  Replicates: {replicates}")
            eq_time = summary.get("equilibration_time")
            if eq_time:
                click.echo(f"  Equilibration time: {eq_time}")
            enabled = summary.get("enabled_analyses", [])
            if enabled:
                click.echo(f"  Enabled analyses: {', '.join(enabled)}")
            else:
                click.echo("  Enabled analyses: (none)")
    else:
        click.secho("✗ Configuration has errors", fg="red")
        click.echo()
        for error in result["errors"]:
            click.echo(f"  • {error}", err=True)


# -----------------------------------------------------------------------------
# Run Command - Execute all enabled analyses
# -----------------------------------------------------------------------------


@analyze.command("run")
@click.option(
    "-c",
    "--config",
    "analysis_config",
    type=click.Path(exists=True, path_type=Path),
    default="analysis.yaml",
    help="Path to analysis.yaml file.",
)
@click.option(
    "--recompute",
    is_flag=True,
    help="Force recompute even if cached results exist.",
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
def run_analyses(analysis_config: Path, recompute: bool, quiet: bool, debug: bool) -> None:
    """Run all enabled analyses defined in analysis.yaml.

    Reads analysis.yaml and runs each enabled analysis type
    (RMSF, distances, catalytic_triad) for all specified replicates.

    \b
    Example:
        polyzymd analyze run
        polyzymd analyze run -c analysis.yaml --recompute
        polyzymd analyze run -q
    """
    require_analysis_deps()

    from polyzymd.analysis.config import AnalysisConfig
    from polyzymd.analysis import RMSFCalculator, DistanceCalculator, CatalyticTriadAnalyzer
    from polyzymd.analysis.core.logging_utils import setup_logging
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.compare.config import CatalyticTriadConfig as CompareTriadConfig
    from polyzymd.compare.config import TriadPairConfig as CompareTriadPairConfig

    # Set up logging with colored output
    setup_logging(quiet=quiet, debug=debug)

    analysis_config = Path(analysis_config).resolve()

    # Load analysis config
    click.echo(f"Loading: {analysis_config.name}")
    try:
        config = AnalysisConfig.from_yaml(analysis_config)
    except FileNotFoundError:
        click.echo(
            click.style(
                f"Error: {analysis_config} not found.",
                fg="red",
            ),
            err=True,
        )
        click.echo("Run 'polyzymd analyze init' first.", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(
            click.style(f"Error loading config: {e}", fg="red"),
            err=True,
        )
        sys.exit(1)

    # Find simulation config.yaml (same directory)
    sim_config_path = analysis_config.parent / "config.yaml"
    if not sim_config_path.exists():
        click.echo(
            click.style(
                "Error: config.yaml not found alongside analysis.yaml",
                fg="red",
            ),
            err=True,
        )
        click.echo("analysis.yaml must be in the same directory as config.yaml.", err=True)
        sys.exit(1)

    # Load simulation config
    sim_config = SimulationConfig.from_yaml(sim_config_path)

    # Validate analysis config
    issues = config.validate_config()
    if issues:
        click.echo(click.style("Configuration issues:", fg="yellow"))
        for issue in issues:
            click.echo(f"  - {issue}")
        click.echo()

    # Get enabled analyses
    enabled = config.get_enabled_analyses()
    if not enabled:
        click.echo(click.style("No analyses enabled in analysis.yaml", fg="yellow"))
        sys.exit(0)

    click.echo(f"Replicates: {config.replicates}")
    click.echo(f"Equilibration: {config.defaults.equilibration_time}")
    click.echo(f"Enabled analyses: {', '.join(enabled)}")
    click.echo()

    # Track results for summary
    completed = []
    failed = []

    # Run RMSF analysis
    if config.rmsf.enabled:
        click.echo(click.style("=" * 60, fg="blue"))
        click.echo(click.style("Running RMSF analysis...", fg="blue"))
        click.echo(click.style("=" * 60, fg="blue"))
        try:
            calculator = RMSFCalculator(
                config=sim_config,
                selection=config.rmsf.selection,
                equilibration=config.defaults.equilibration_time,
                reference_mode=config.rmsf.reference_mode,
                reference_frame=config.rmsf.reference_frame,
            )
            for rep in config.replicates:
                click.echo(f"  Replicate {rep}...", nl=False)
                try:
                    result = calculator.compute(replicate=rep, recompute=recompute)
                    click.echo(
                        click.style(f" done (mean RMSF: {result.mean_rmsf:.2f} Å)", fg="green")
                    )
                except Exception as e:
                    click.echo(click.style(f" FAILED: {e}", fg="red"))
                    failed.append(f"RMSF rep {rep}")
            completed.append(f"RMSF ({len(config.replicates)} replicates)")
        except Exception as e:
            click.echo(click.style(f"RMSF analysis failed: {e}", fg="red"))
            failed.append("RMSF")

    # Run distance analysis
    if config.distances.enabled:
        click.echo()
        click.echo(click.style("=" * 60, fg="blue"))
        click.echo(click.style("Running distance analysis...", fg="blue"))
        click.echo(click.style("=" * 60, fg="blue"))
        try:
            pairs = [(p.selection_a, p.selection_b) for p in config.distances.pairs]
            calculator = DistanceCalculator(
                config=sim_config,
                pairs=pairs,
                equilibration=config.defaults.equilibration_time,
            )
            for rep in config.replicates:
                click.echo(f"  Replicate {rep}...", nl=False)
                try:
                    result = calculator.compute(replicate=rep, recompute=recompute)
                    click.echo(click.style(f" done ({len(pairs)} pairs)", fg="green"))
                except Exception as e:
                    click.echo(click.style(f" FAILED: {e}", fg="red"))
                    failed.append(f"Distances rep {rep}")
            completed.append(f"Distances ({len(config.replicates)} replicates)")
        except Exception as e:
            click.echo(click.style(f"Distance analysis failed: {e}", fg="red"))
            failed.append("Distances")

    # Run catalytic triad analysis
    if config.catalytic_triad.enabled:
        click.echo()
        click.echo(click.style("=" * 60, fg="blue"))
        click.echo(click.style("Running catalytic triad analysis...", fg="blue"))
        click.echo(click.style("=" * 60, fg="blue"))
        try:
            # Convert to compare module's TriadConfig format
            triad_pairs = [
                CompareTriadPairConfig(
                    label=p.label,
                    selection_a=p.selection_a,
                    selection_b=p.selection_b,
                )
                for p in config.catalytic_triad.pairs
            ]
            triad_config = CompareTriadConfig(
                name=config.catalytic_triad.name,
                threshold=config.catalytic_triad.threshold,
                pairs=triad_pairs,
            )
            analyzer = CatalyticTriadAnalyzer(
                config=sim_config,
                triad_config=triad_config,
                equilibration=config.defaults.equilibration_time,
            )
            for rep in config.replicates:
                click.echo(f"  Replicate {rep}...", nl=False)
                try:
                    result = analyzer.compute(replicate=rep, recompute=recompute)
                    contact_pct = result.simultaneous_contact_fraction * 100
                    click.echo(
                        click.style(
                            f" done (simultaneous contact: {contact_pct:.1f}%)",
                            fg="green",
                        )
                    )
                except Exception as e:
                    click.echo(click.style(f" FAILED: {e}", fg="red"))
                    failed.append(f"Triad rep {rep}")
            completed.append(f"Catalytic triad ({len(config.replicates)} replicates)")
        except Exception as e:
            click.echo(click.style(f"Catalytic triad analysis failed: {e}", fg="red"))
            failed.append("Catalytic triad")

    # Run contacts analysis
    if config.contacts.enabled:
        click.echo()
        click.echo(click.style("=" * 60, fg="blue"))
        click.echo(click.style("Running contact analysis...", fg="blue"))
        click.echo(click.style("=" * 60, fg="blue"))
        try:
            from polyzymd.analysis.contacts import ParallelContactAnalyzer
            from polyzymd.analysis.contacts.aggregator import aggregate_contact_results
            from polyzymd.analysis.common.selectors import MDAnalysisSelector
            from polyzymd.analysis.core.loader import (
                TrajectoryLoader,
                parse_time_string,
                time_to_frame,
            )

            # Parse equilibration time
            eq_value, eq_unit = parse_time_string(config.defaults.equilibration_time)

            # Create selectors
            target_selector = MDAnalysisSelector(config.contacts.protein_selection)
            query_selector = MDAnalysisSelector(config.contacts.polymer_selection)

            # Use TrajectoryLoader for consistent path resolution
            loader = TrajectoryLoader(sim_config)

            contact_results = []
            binding_pref_results = []

            for rep in config.replicates:
                click.echo(f"  Replicate {rep}...", nl=False)
                try:
                    # Load universe
                    universe = loader.load_universe(rep)

                    # Convert equilibration time to start frame
                    timestep = loader.get_timestep(rep)
                    if eq_unit == "ns":
                        eq_time_ps = eq_value * 1000
                    elif eq_unit == "us":
                        eq_time_ps = eq_value * 1e6
                    else:
                        eq_time_ps = eq_value
                    start_frame = time_to_frame(eq_time_ps, "ps", timestep, "ps")

                    # Create analyzer
                    analyzer = ParallelContactAnalyzer(
                        target_selector=target_selector,
                        query_selector=query_selector,
                        cutoff=config.contacts.cutoff,
                    )

                    # Run analysis
                    result = analyzer.run(universe, start=start_frame)
                    contact_results.append(result)

                    # Report basic stats
                    n_contacted = result.n_contacted_residues
                    coverage = result.coverage_fraction()
                    click.echo(
                        click.style(
                            f" done ({n_contacted}/{result.n_protein_residues} residues, "
                            f"{coverage:.1%} coverage)",
                            fg="green",
                        )
                    )

                    # Compute binding preference if enabled
                    if config.contacts.compute_binding_preference:
                        from polyzymd.analysis.contacts import (
                            compute_binding_preference_from_config,
                        )

                        # Get enzyme PDB path
                        enzyme_pdb = config.contacts.enzyme_pdb_for_sasa
                        if enzyme_pdb is None:
                            enzyme_pdb = sim_config.system.enzyme_pdb
                        enzyme_pdb = Path(enzyme_pdb)
                        if not enzyme_pdb.is_absolute():
                            enzyme_pdb = sim_config_path.parent / enzyme_pdb

                        click.echo(f"    Computing binding preference...", nl=False)
                        bp_result = compute_binding_preference_from_config(
                            contact_result=result,
                            universe=universe,
                            enzyme_pdb_path=enzyme_pdb,
                            config=config.contacts,
                        )
                        binding_pref_results.append(bp_result)
                        n_entries = len(bp_result.entries)
                        click.echo(click.style(f" done ({n_entries} entries)", fg="green"))

                except Exception as e:
                    click.echo(click.style(f" FAILED: {e}", fg="red"))
                    failed.append(f"Contacts rep {rep}")

            # Aggregate and save results
            if contact_results:
                # Save individual results
                analysis_dir = sim_config.output.projects_directory / "analysis" / "contacts"
                analysis_dir.mkdir(parents=True, exist_ok=True)

                for i, (rep, result) in enumerate(zip(config.replicates, contact_results)):
                    output_file = analysis_dir / f"contacts_rep{rep}.json"
                    result.save(output_file)

                    if binding_pref_results and i < len(binding_pref_results):
                        bp_file = analysis_dir / f"binding_preference_rep{rep}.json"
                        binding_pref_results[i].save(bp_file)

                # Aggregate if multiple replicates
                if len(contact_results) > 1:
                    click.echo("  Aggregating results...", nl=False)
                    agg_result = aggregate_contact_results(contact_results)
                    agg_file = (
                        analysis_dir
                        / f"contacts_aggregated_reps{config.replicates[0]}-{config.replicates[-1]}.json"
                    )
                    import json

                    with open(agg_file, "w") as f:
                        json.dump(agg_result.to_dict(), f, indent=2)

                    # Aggregate binding preference
                    if binding_pref_results:
                        from polyzymd.analysis.contacts import aggregate_binding_preference

                        agg_bp = aggregate_binding_preference(binding_pref_results)
                        agg_bp_file = (
                            analysis_dir
                            / f"binding_preference_aggregated_reps{config.replicates[0]}-{config.replicates[-1]}.json"
                        )
                        agg_bp.save(agg_bp_file)
                        click.echo(click.style(" done (contacts + binding preference)", fg="green"))
                    else:
                        click.echo(click.style(" done", fg="green"))

                click.echo(f"  Results saved: {analysis_dir}")

            completed.append(f"Contacts ({len(config.replicates)} replicates)")
        except Exception as e:
            click.echo(click.style(f"Contact analysis failed: {e}", fg="red"))
            if LOGGER.level == logging.DEBUG:
                import traceback

                traceback.print_exc()
            failed.append("Contacts")

    # Print summary
    click.echo()
    click.echo(click.style("=" * 60, fg="green" if not failed else "yellow"))
    click.echo(click.style("Summary", fg="green" if not failed else "yellow"))
    click.echo(click.style("=" * 60, fg="green" if not failed else "yellow"))
    if completed:
        click.echo("Completed:")
        for item in completed:
            click.echo(f"  - {item}")
    if failed:
        click.echo(click.style("Failed:", fg="red"))
        for item in failed:
            click.echo(click.style(f"  - {item}", fg="red"))

    if not failed:
        click.echo()
        click.echo(click.style("All analyses completed successfully!", fg="green"))
    else:
        sys.exit(1)


@analyze.command()
@click.option(
    "-c",
    "--config",
    required=True,
    type=click.Path(exists=True),
    help="Path to YAML configuration file",
)
@click.option(
    "-r",
    "--replicates",
    default="1",
    help="Replicate specification: '1-5', '1,3,5', or '1'",
)
@click.option(
    "--eq-time",
    default="0ns",
    help="Equilibration time to skip: '100ns', '5000ps' [default: 0ns]",
)
@click.option(
    "--selection",
    default="protein and name CA",
    help="MDAnalysis selection string [default: 'protein and name CA']",
)
@click.option(
    "--reference-mode",
    type=click.Choice(["centroid", "average", "frame"]),
    default="centroid",
    help=(
        "Reference structure for alignment: "
        "'centroid' (most populated state), "
        "'average' (mean structure), "
        "'frame' (user-specified) [default: centroid]"
    ),
)
@click.option(
    "--reference-frame",
    type=int,
    default=None,
    help="Reference frame when --reference-mode=frame (1-indexed)",
)
@click.option(
    "--alignment-selection",
    default="protein and name CA",
    help="Selection for trajectory alignment [default: 'protein and name CA']",
)
@click.option(
    "--centroid-selection",
    default="protein",
    help="Selection for centroid finding [default: 'protein']",
)
@click.option(
    "--plot",
    is_flag=True,
    help="Generate plot after analysis",
)
@click.option(
    "--recompute",
    is_flag=True,
    help="Force recompute even if cached results exist",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    default=None,
    help="Custom output directory for results",
)
def rmsf(
    config: str,
    replicates: str,
    eq_time: str,
    selection: str,
    reference_mode: str,
    reference_frame: Optional[int],
    alignment_selection: str,
    centroid_selection: str,
    plot: bool,
    recompute: bool,
    output_dir: Optional[str],
) -> None:
    """Compute RMSF (Root Mean Square Fluctuation) analysis.

    Calculates per-residue flexibility from MD trajectories with proper
    trajectory alignment and statistical handling.

    The trajectory is aligned to a reference structure before RMSF computation.
    Choose the reference mode based on your scientific question:

    \b
    - centroid: Most populated conformational state (K-Means clustering).
                Best for equilibrium flexibility analysis.
    - average:  Mathematical mean structure. Pure thermal fluctuations.
    - frame:    User-specified frame. Best for functional state analysis
                (e.g., catalytically competent conformation).

    \b
    Examples:
        # Default: align to most populated state
        polyzymd analyze rmsf -c config.yaml -r 1 --eq-time 100ns

        # Use average structure as reference
        polyzymd analyze rmsf -c config.yaml -r 1 --reference-mode average

        # Align to a specific frame (e.g., catalytically competent)
        polyzymd analyze rmsf -c config.yaml -r 1 --reference-mode frame --reference-frame 500
    """
    require_analysis_deps()

    from polyzymd.analysis import RMSFCalculator
    from polyzymd.config.schema import SimulationConfig

    # Validate reference-mode and reference-frame combination
    if reference_mode == "frame" and reference_frame is None:
        click.echo(
            click.style(
                "Error: --reference-frame is required when --reference-mode=frame",
                fg="red",
            ),
            err=True,
        )
        sys.exit(1)

    # Parse inputs
    rep_list = parse_replicates(replicates)
    output_path = Path(output_dir) if output_dir else None

    click.echo(f"Loading configuration from: {config}")
    try:
        sim_config = SimulationConfig.from_yaml(config)
    except Exception as e:
        click.echo(click.style(f"Error loading config: {e}", fg="red"), err=True)
        sys.exit(1)

    click.echo(f"RMSF Analysis: {sim_config.name}")
    click.echo(f"  Replicates: {replicates}")
    click.echo(f"  Equilibration: {eq_time}")
    if eq_time in ("0ns", "0ps", "0us"):
        click.echo(
            click.style(
                "  Warning: No equilibration time specified. RMSF will include "
                "potentially non-equilibrated frames. Consider using --eq-time.",
                fg="yellow",
            ),
            err=True,
        )
    click.echo(f"  Selection: {selection}")
    click.echo(f"  Alignment: {reference_mode}", nl=False)
    if reference_mode == "frame":
        click.echo(f" (frame {reference_frame})")
    else:
        click.echo()

    try:
        # Create calculator
        calc = RMSFCalculator(
            config=sim_config,
            selection=selection,
            equilibration=eq_time,
            reference_mode=reference_mode,
            reference_frame=reference_frame,
            alignment_selection=alignment_selection,
            centroid_selection=centroid_selection,
        )

        if len(rep_list) == 1:
            # Single replicate
            result = calc.compute(
                replicate=rep_list[0],
                save=True,
                output_dir=output_path,
                recompute=recompute,
            )
            click.echo()
            click.echo(click.style("RMSF Analysis Complete", fg="green"))
            click.echo(f"  Mean RMSF: {result.mean_rmsf:.3f} Å")
            click.echo(f"  Min RMSF:  {result.min_rmsf:.3f} Å")
            click.echo(f"  Max RMSF:  {result.max_rmsf:.3f} Å")
            click.echo(f"  Frames used: {result.n_frames_used}")

            if result.correlation_time:
                click.echo(
                    f"  Correlation time: {result.correlation_time:.2f} {result.correlation_time_unit}"
                )

            # Generate plot if requested
            if plot:
                require_matplotlib()
                from polyzymd.analysis.rmsf import plot_rmsf

                fig, ax = plot_rmsf(result)
                plot_dir = sim_config.output.projects_directory / "plots" / "rmsf"
                plot_dir.mkdir(parents=True, exist_ok=True)
                plot_path = plot_dir / f"rmsf_run{rep_list[0]}.png"
                fig.savefig(plot_path, dpi=150, bbox_inches="tight")
                click.echo(f"  Plot saved: {plot_path}")

        else:
            # Multiple replicates - compute aggregated
            result = calc.compute_aggregated(
                replicates=rep_list,
                save=True,
                output_dir=output_path,
                recompute=recompute,
            )
            click.echo()
            click.echo(click.style("RMSF Analysis Complete (Aggregated)", fg="green"))
            click.echo(f"  Replicates: {result.replicate_range}")
            click.echo(
                f"  Mean RMSF: {result.overall_mean_rmsf:.3f} ± {result.overall_sem_rmsf:.3f} Å"
            )
            click.echo(f"  Min RMSF:  {result.overall_min_rmsf:.3f} Å")
            click.echo(f"  Max RMSF:  {result.overall_max_rmsf:.3f} Å")

            # Generate plot if requested
            if plot:
                require_matplotlib()
                from polyzymd.analysis.rmsf import plot_rmsf

                fig, ax = plot_rmsf(result, show_error=True)
                plot_dir = sim_config.output.projects_directory / "plots" / "rmsf"
                plot_dir.mkdir(parents=True, exist_ok=True)
                plot_path = plot_dir / f"rmsf_aggregated_reps{rep_list[0]}-{rep_list[-1]}.png"
                fig.savefig(plot_path, dpi=150, bbox_inches="tight")
                click.echo(f"  Plot saved: {plot_path}")

    except Exception as e:
        click.echo(click.style(f"Analysis failed: {e}", fg="red"), err=True)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


# -----------------------------------------------------------------------------
# Distance Analysis
# -----------------------------------------------------------------------------


@analyze.command()
@click.option(
    "-c",
    "--config",
    required=True,
    type=click.Path(exists=True),
    help="Path to YAML configuration file",
)
@click.option(
    "-r",
    "--replicates",
    default="1",
    help="Replicate specification: '1-5', '1,3,5', or '1'",
)
@click.option(
    "--eq-time",
    default="0ns",
    help="Equilibration time to skip: '100ns', '5000ps' [default: 0ns]",
)
@click.option(
    "--pair",
    "pairs",
    multiple=True,
    required=True,
    help="Distance pair as 'selection1 : selection2'. Can specify multiple times.",
)
@click.option(
    "--threshold",
    type=float,
    default=None,
    help="Distance threshold for contact analysis (Angstroms)",
)
@click.option(
    "--plot",
    is_flag=True,
    help="Generate plots after analysis",
)
@click.option(
    "--recompute",
    is_flag=True,
    help="Force recompute even if cached results exist",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    default=None,
    help="Custom output directory for results",
)
def distances(
    config: str,
    replicates: str,
    eq_time: str,
    pairs: tuple[str, ...],
    threshold: Optional[float],
    plot: bool,
    recompute: bool,
    output_dir: Optional[str],
) -> None:
    """Compute inter-atomic distance analysis.

    Calculates distances between specified atom pairs across trajectories,
    with optional contact analysis (fraction of frames below threshold).

    \b
    Examples:
        # Single distance pair
        polyzymd analyze distances -c config.yaml -r 1 \\
            --pair "resid 77 and name OG : resid 133 and name NE2"

        # Multiple pairs with contact threshold
        polyzymd analyze distances -c config.yaml -r 1-5 \\
            --pair "resid 77 and name OG : resid 133 and name NE2" \\
            --pair "resid 133 and name NE2 : resid 156 and name OD1" \\
            --threshold 3.5 --plot
    """
    require_analysis_deps()

    from polyzymd.analysis import DistanceCalculator
    from polyzymd.config.schema import SimulationConfig

    # Parse inputs
    rep_list = parse_replicates(replicates)
    output_path = Path(output_dir) if output_dir else None

    try:
        distance_pairs = parse_distance_pairs(pairs)
    except click.BadParameter as e:
        click.echo(click.style(str(e), fg="red"), err=True)
        sys.exit(1)

    click.echo(f"Loading configuration from: {config}")
    try:
        sim_config = SimulationConfig.from_yaml(config)
    except Exception as e:
        click.echo(click.style(f"Error loading config: {e}", fg="red"), err=True)
        sys.exit(1)

    click.echo(f"Distance Analysis: {sim_config.name}")
    click.echo(f"  Replicates: {replicates}")
    click.echo(f"  Equilibration: {eq_time}")
    if eq_time in ("0ns", "0ps", "0us"):
        click.echo(
            click.style(
                "  Warning: No equilibration time specified. Analysis will include "
                "potentially non-equilibrated frames. Consider using --eq-time.",
                fg="yellow",
            ),
            err=True,
        )
    click.echo(f"  Distance pairs: {len(distance_pairs)}")
    for i, (sel1, sel2) in enumerate(distance_pairs, 1):
        click.echo(f"    {i}. {sel1} <-> {sel2}")
    if threshold:
        click.echo(f"  Contact threshold: {threshold} Å")

    try:
        # Create calculator
        calc = DistanceCalculator(
            config=sim_config,
            pairs=distance_pairs,
            equilibration=eq_time,
            thresholds=threshold,
        )

        if len(rep_list) == 1:
            # Single replicate
            result = calc.compute(
                replicate=rep_list[0],
                save=True,
                output_dir=output_path,
                recompute=recompute,
            )
            click.echo()
            click.echo(click.style("Distance Analysis Complete", fg="green"))
            for pr in result.pair_results:
                click.echo(f"  {pr.pair_label}:")
                click.echo(f"    Mean: {pr.mean_distance:.2f} ± {pr.std_distance:.2f} Å")
                click.echo(f"    Min:  {pr.min_distance:.2f} Å")
                click.echo(f"    Max:  {pr.max_distance:.2f} Å")
                if pr.fraction_below_threshold is not None:
                    click.echo(
                        f"    Contact fraction (<{threshold}Å): {pr.fraction_below_threshold:.1%}"
                    )

            # Generate plot if requested
            if plot:
                require_matplotlib()
                from polyzymd.analysis.distances import plot_distance_histogram

                plot_dir = sim_config.output.projects_directory / "plots" / "distances"
                plot_dir.mkdir(parents=True, exist_ok=True)
                for pr in result.pair_results:
                    fig, ax = plot_distance_histogram(pr)
                    plot_path = plot_dir / f"dist_{pr.pair_label}_run{rep_list[0]}.png"
                    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
                    click.echo(f"  Plot saved: {plot_path}")

        else:
            # Multiple replicates - compute aggregated
            result = calc.compute_aggregated(
                replicates=rep_list,
                save=True,
                output_dir=output_path,
                recompute=recompute,
            )
            click.echo()
            click.echo(click.style("Distance Analysis Complete (Aggregated)", fg="green"))
            click.echo(f"  Replicates: {result.replicate_range}")
            for pr in result.pair_results:
                click.echo(f"  {pr.pair_label}:")
                click.echo(f"    Mean: {pr.overall_mean:.2f} ± {pr.overall_sem:.2f} Å")
                if pr.overall_fraction_below is not None:
                    click.echo(
                        f"    Contact fraction: {pr.overall_fraction_below:.1%} ± {pr.sem_fraction_below:.1%}"
                    )

            # Generate plot if requested
            if plot:
                require_matplotlib()
                from polyzymd.analysis.distances import plot_contact_fraction_bar

                plot_dir = sim_config.output.projects_directory / "plots" / "distances"
                plot_dir.mkdir(parents=True, exist_ok=True)
                if threshold:
                    # Plot contact fraction bar chart for all pairs
                    fig, ax = plot_contact_fraction_bar(result.pair_results)
                    plot_path = (
                        plot_dir
                        / f"contact_fraction_aggregated_reps{rep_list[0]}-{rep_list[-1]}.png"
                    )
                    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
                    click.echo(f"  Plot saved: {plot_path}")

    except Exception as e:
        click.echo(click.style(f"Analysis failed: {e}", fg="red"), err=True)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


# -----------------------------------------------------------------------------
# Contact Analysis
# -----------------------------------------------------------------------------


@analyze.command()
@click.option(
    "-c",
    "--config",
    required=True,
    type=click.Path(exists=True),
    help="Path to YAML configuration file (config.yaml or analysis.yaml)",
)
@click.option(
    "-r",
    "--replicates",
    default="1",
    help="Replicate specification: '1-5', '1,3,5', or '1'",
)
@click.option(
    "--eq-time",
    default="0ns",
    help="Equilibration time to skip: '100ns', '5000ps' [default: 0ns]",
)
@click.option(
    "--cutoff",
    type=float,
    default=4.0,
    help="Contact distance cutoff in Angstroms [default: 4.0]",
)
@click.option(
    "--polymer-selection",
    default="segid C",
    help="MDAnalysis selection for polymer atoms [default: 'segid C']",
)
@click.option(
    "--protein-selection",
    default="protein",
    help="MDAnalysis selection for protein atoms [default: 'protein']",
)
@click.option(
    "--grouping",
    type=click.Choice(["residue", "segment", "atom"]),
    default="residue",
    help="How to group contacts [default: residue]",
)
@click.option(
    "--residence-times",
    is_flag=True,
    help="Compute residence time statistics",
)
@click.option(
    "--plot",
    is_flag=True,
    help="Generate contact map plot after analysis",
)
@click.option(
    "--recompute",
    is_flag=True,
    help="Force recompute even if cached results exist",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    default=None,
    help="Custom output directory for results",
)
@click.option(
    "--binding-preference",
    is_flag=True,
    help="Compute binding preference analysis with enrichment ratios (requires rust-sasa-python)",
)
@click.option(
    "--enzyme-pdb",
    type=click.Path(exists=True),
    default=None,
    help="Path to enzyme PDB for SASA calculation (default: uses config enzyme_pdb)",
)
@click.option(
    "--surface-threshold",
    type=float,
    default=0.2,
    help="Relative SASA threshold for surface exposure (default: 0.2 = 20%)",
)
def contacts(
    config: str,
    replicates: str,
    eq_time: str,
    cutoff: float,
    polymer_selection: str,
    protein_selection: str,
    grouping: str,
    residence_times: bool,
    plot: bool,
    recompute: bool,
    output_dir: Optional[str],
    binding_preference: bool,
    enzyme_pdb: Optional[str],
    surface_threshold: float,
) -> None:
    """Analyze polymer-protein contacts.

    Computes contact frequencies between polymer and protein residues,
    with statistical analysis following LiveCoMS best practices.

    \b
    PolyzyMD Chain Convention:
      Chain A: Protein/Enzyme
      Chain B: Substrate/Ligand
      Chain C: Polymers
      Chain D+: Solvent (water, ions, co-solvents)

    IMPORTANT: This command uses solvated_system.pdb as the topology
    to ensure correct chain assignments (NOT production_N_topology.pdb).

    \b
    Examples:
        # Basic contact analysis
        polyzymd analyze contacts -c config.yaml -r 1 --eq-time 100ns

        # Multi-replicate analysis
        polyzymd analyze contacts -c config.yaml -r 1-5 --eq-time 100ns

        # Custom selections with residence time analysis
        polyzymd analyze contacts -c config.yaml -r 1-3 \\
            --polymer-selection "segid C and resname SBM" \\
            --protein-selection "protein and (resname TRP PHE TYR)" \\
            --residence-times --plot

        # With binding preference analysis (enrichment ratios)
        polyzymd analyze contacts -c config.yaml -r 1-3 \\
            --binding-preference --surface-threshold 0.2
    """
    require_analysis_deps()

    from polyzymd.analysis.contacts import ParallelContactAnalyzer
    from polyzymd.analysis.contacts.aggregator import aggregate_contact_results
    from polyzymd.analysis.common.selectors import MDAnalysisSelector
    from polyzymd.analysis.core.loader import TrajectoryLoader, parse_time_string, time_to_frame
    from polyzymd.config.schema import SimulationConfig

    # Parse inputs
    rep_list = parse_replicates(replicates)
    output_path = Path(output_dir) if output_dir else None

    click.echo(f"Loading configuration from: {config}")
    try:
        sim_config = SimulationConfig.from_yaml(config)
    except Exception as e:
        click.echo(click.style(f"Error loading config: {e}", fg="red"), err=True)
        sys.exit(1)

    click.echo(f"Contact Analysis: {sim_config.name}")
    click.echo(f"  Replicates: {replicates}")
    click.echo(f"  Equilibration: {eq_time}")
    if eq_time in ("0ns", "0ps", "0us"):
        click.echo(
            click.style(
                "  Warning: No equilibration time specified. Analysis will include "
                "potentially non-equilibrated frames. Consider using --eq-time.",
                fg="yellow",
            ),
            err=True,
        )
    click.echo(f"  Cutoff: {cutoff} Å")
    click.echo(f"  Polymer selection: {polymer_selection}")
    click.echo(f"  Protein selection: {protein_selection}")
    click.echo(f"  Grouping: {grouping}")
    if binding_preference:
        click.echo(f"  Binding preference: enabled (threshold={surface_threshold})")

    # Parse equilibration time
    eq_value, eq_unit = parse_time_string(eq_time)

    try:
        # Create selectors
        target_selector = MDAnalysisSelector(protein_selection)
        query_selector = MDAnalysisSelector(polymer_selection)

        # Use TrajectoryLoader for consistent path resolution
        loader = TrajectoryLoader(sim_config)

        # Process each replicate
        results = []
        binding_pref_results = []
        universes = []  # Keep universe references for binding preference
        for rep in rep_list:
            click.echo(f"  Processing replicate {rep}...", nl=False)

            try:
                # Use TrajectoryLoader for consistent path resolution
                # This correctly uses scratch_directory for trajectories
                universe = loader.load_universe(rep)
            except FileNotFoundError as e:
                click.echo(click.style(f" FAILED: {e}", fg="red"))
                continue

            # Convert equilibration time to start frame
            timestep = loader.get_timestep(rep)
            # Convert eq_time to ps for consistent units
            if eq_unit == "ns":
                eq_time_ps = eq_value * 1000
            elif eq_unit == "us":
                eq_time_ps = eq_value * 1e6
            else:
                eq_time_ps = eq_value  # already ps
            start_frame = time_to_frame(eq_time_ps, "ps", timestep, "ps")

            # Create analyzer
            analyzer = ParallelContactAnalyzer(
                target_selector=target_selector,
                query_selector=query_selector,
                cutoff=cutoff,
            )

            # Run analysis
            result = analyzer.run(
                universe,
                start=start_frame,
            )

            results.append(result)
            universes.append(universe)

            # Report basic stats
            n_contacted = result.n_contacted_residues
            coverage = result.coverage_fraction()
            mean_frac = result.mean_contact_fraction()
            click.echo(
                click.style(
                    f" done ({n_contacted}/{result.n_protein_residues} residues contacted, "
                    f"{coverage:.1%} coverage, {mean_frac:.1%} mean contact)",
                    fg="green",
                )
            )

        if not results:
            click.echo(click.style("No results generated!", fg="red"), err=True)
            sys.exit(1)

        # Compute binding preference if enabled
        if binding_preference and results:
            click.echo()
            click.echo("Computing binding preference...")

            from polyzymd.analysis.contacts import (
                SurfaceExposureFilter,
                compute_binding_preference,
                resolve_protein_group_selections,
                aggregate_binding_preference,
            )

            # Determine enzyme PDB path
            if enzyme_pdb:
                enzyme_pdb_path = Path(enzyme_pdb)
            else:
                enzyme_pdb_path = sim_config.system.enzyme_pdb
                if not Path(enzyme_pdb_path).is_absolute():
                    enzyme_pdb_path = Path(config).parent / enzyme_pdb_path

            # Calculate surface exposure (once - same for all replicates)
            click.echo(f"  Calculating surface exposure from: {enzyme_pdb_path}")
            exposure_filter = SurfaceExposureFilter(threshold=surface_threshold)
            try:
                surface_exposure = exposure_filter.calculate(enzyme_pdb_path)
                click.echo(
                    f"  Found {surface_exposure.exposed_count}/{surface_exposure.total_count} "
                    f"surface-exposed residues"
                )
            except Exception as e:
                click.echo(
                    click.style(f"  SASA calculation failed: {e}", fg="red"),
                    err=True,
                )
                click.echo(
                    "  Install rust-sasa-python: pip install rust-sasa-python",
                    err=True,
                )
                binding_preference = False  # Disable for rest of processing

            if binding_preference:
                # Resolve protein groups (use default AA class selections)
                universe = universes[0]  # Use first universe for selection resolution
                protein_groups = resolve_protein_group_selections(universe, None)
                click.echo(f"  Protein groups: {', '.join(protein_groups.keys())}")

                # Compute binding preference for each replicate
                for i, (result, universe) in enumerate(zip(results, universes)):
                    try:
                        bp_result = compute_binding_preference(
                            contact_result=result,
                            surface_exposure=surface_exposure,
                            protein_groups=protein_groups,
                        )
                        binding_pref_results.append(bp_result)

                        # Save per-replicate binding preference file
                        if output_path:
                            bp_rep_file = output_path / f"binding_preference_rep{rep_list[i]}.json"
                        else:
                            bp_rep_file = (
                                sim_config.output.projects_directory
                                / "analysis"
                                / "contacts"
                                / f"binding_preference_rep{rep_list[i]}.json"
                            )
                        bp_rep_file.parent.mkdir(parents=True, exist_ok=True)
                        bp_result.save(bp_rep_file)

                        click.echo(
                            f"  Replicate {rep_list[i]}: "
                            f"{len(bp_result.polymer_types())} polymer types × "
                            f"{len(bp_result.protein_groups())} groups"
                        )
                    except Exception as e:
                        click.echo(
                            click.style(
                                f"  Replicate {rep_list[i]} binding preference FAILED: {e}",
                                fg="red",
                            )
                        )

        # Aggregate if multiple replicates
        if len(results) > 1:
            click.echo()
            click.echo("Aggregating results across replicates...")
            agg_result = aggregate_contact_results(results)

            click.echo(click.style("Aggregated Contact Analysis Complete", fg="green"))
            click.echo(
                f"  Contact fraction: {agg_result.mean_contact_fraction:.1%} "
                f"± {agg_result.mean_contact_fraction_sem:.1%}"
            )

            if residence_times and agg_result.residence_time_by_polymer_type:
                click.echo("  Residence time by polymer type:")
                for ptype, (rt_mean, rt_sem) in sorted(
                    agg_result.residence_time_by_polymer_type.items()
                ):
                    click.echo(f"    {ptype}: {rt_mean:.2f} ± {rt_sem:.2f} frames")

            # Save aggregated result
            if output_path:
                output_file = output_path / "contacts_aggregated.json"
            else:
                analysis_dir = sim_config.output.projects_directory / "analysis" / "contacts"
                analysis_dir.mkdir(parents=True, exist_ok=True)
                output_file = (
                    analysis_dir / f"contacts_aggregated_reps{rep_list[0]}-{rep_list[-1]}.json"
                )

            output_file.parent.mkdir(parents=True, exist_ok=True)
            import json

            with open(output_file, "w") as f:
                json.dump(agg_result.to_dict(), f, indent=2)
            click.echo(f"  Results saved: {output_file}")

            # Aggregate and save binding preference if computed
            if binding_pref_results:
                from polyzymd.analysis.contacts import aggregate_binding_preference

                agg_bp = aggregate_binding_preference(binding_pref_results)
                bp_file = (
                    output_file.parent
                    / f"binding_preference_aggregated_reps{rep_list[0]}-{rep_list[-1]}.json"
                )
                agg_bp.save(bp_file)
                click.echo(f"  Binding preference saved: {bp_file}")

                # Print enrichment summary grouped by polymer type
                click.echo("  Enrichment summary (mean ± SEM):")
                for ptype in agg_bp.polymer_types():
                    click.echo(f"    {ptype}:")
                    for entry in agg_bp.entries:
                        if entry.polymer_type == ptype and entry.mean_enrichment is not None:
                            click.echo(
                                f"      {entry.protein_group}: "
                                f"{entry.mean_enrichment:.2f} ± {entry.sem_enrichment or 0:.2f}"
                            )

                # Print surface-exposed coverage
                exposed_resids = surface_exposure.exposed_resids
                contacted_exposed = sum(
                    1
                    for rs in agg_result.residue_stats
                    if rs.protein_resid in exposed_resids and rs.contact_fraction_mean > 0
                )
                click.echo(
                    f"  Surface-exposed coverage: {contacted_exposed}/{len(exposed_resids)} "
                    f"({contacted_exposed / len(exposed_resids):.1%})"
                )

        else:
            # Single replicate
            result = results[0]
            click.echo()
            click.echo(click.style("Contact Analysis Complete", fg="green"))
            click.echo(
                f"  Contacted residues: {result.n_contacted_residues}/{result.n_protein_residues}"
            )
            click.echo(f"  Coverage: {result.coverage_fraction():.1%}")
            click.echo(f"  Mean contact fraction: {result.mean_contact_fraction():.1%}")

            if residence_times:
                residence_summary = result.residence_time_summary()
                if residence_summary:
                    click.echo("  Residence time by polymer type:")
                    for ptype, stats in sorted(residence_summary.items()):
                        if stats["total_events"] > 0:
                            click.echo(
                                f"    {ptype}: mean={stats['mean_frames']:.2f} frames, "
                                f"max={stats['max_frames']:.0f} frames "
                                f"({stats['total_events']} events)"
                            )

            # Save result
            if output_path:
                output_file = output_path / f"contacts_rep{rep_list[0]}.json"
            else:
                analysis_dir = sim_config.output.projects_directory / "analysis" / "contacts"
                analysis_dir.mkdir(parents=True, exist_ok=True)
                output_file = analysis_dir / f"contacts_rep{rep_list[0]}.json"

            output_file.parent.mkdir(parents=True, exist_ok=True)
            result.save(output_file)
            click.echo(f"  Results saved: {output_file}")

            # Save binding preference if computed
            if binding_pref_results:
                bp_file = output_file.parent / f"binding_preference_rep{rep_list[0]}.json"
                binding_pref_results[0].save(bp_file)
                click.echo(f"  Binding preference saved: {bp_file}")

                # Print enrichment summary grouped by polymer type
                bp_result = binding_pref_results[0]
                click.echo("  Enrichment summary:")
                for ptype in bp_result.polymer_types():
                    click.echo(f"    {ptype}:")
                    for entry in bp_result.entries:
                        if entry.polymer_type == ptype and entry.enrichment_ratio is not None:
                            click.echo(f"      {entry.protein_group}: {entry.enrichment_ratio:.2f}")

                # Print surface-exposed coverage
                exposed_resids = surface_exposure.exposed_resids
                contacted_exposed = sum(
                    1
                    for rc in result.residue_contacts
                    if rc.protein_resid in exposed_resids and rc.total_contact_events > 0
                )
                click.echo(
                    f"  Surface-exposed coverage: {contacted_exposed}/{len(exposed_resids)} "
                    f"({contacted_exposed / len(exposed_resids):.1%})"
                )

        # Generate plot if requested
        if plot:
            require_matplotlib()
            click.echo()
            click.echo("Generating contact map plot...")
            # TODO: Implement contact map plotting
            click.echo(
                click.style(
                    "  Note: Contact map plotting not yet implemented",
                    fg="yellow",
                )
            )

    except Exception as e:
        click.echo(click.style(f"Analysis failed: {e}", fg="red"), err=True)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


# -----------------------------------------------------------------------------
# Catalytic Triad Analysis
# -----------------------------------------------------------------------------


@analyze.command()
@click.option(
    "-c",
    "--comparison",
    required=True,
    type=click.Path(exists=True),
    help="Path to comparison.yaml file with catalytic_triad configuration",
)
@click.option(
    "--condition",
    default=None,
    help="Condition label to analyze (default: analyze all conditions)",
)
@click.option(
    "-r",
    "--replicates",
    default=None,
    help="Replicate specification: '1-5', '1,3,5', or '1'. Overrides condition replicates.",
)
@click.option(
    "--eq-time",
    default=None,
    help="Equilibration time to skip (overrides comparison.yaml defaults)",
)
@click.option(
    "--recompute",
    is_flag=True,
    help="Force recompute even if cached results exist",
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(),
    default=None,
    help="Custom output directory for results",
)
def triad(
    comparison: str,
    condition: Optional[str],
    replicates: Optional[str],
    eq_time: Optional[str],
    recompute: bool,
    output_dir: Optional[str],
) -> None:
    """Analyze catalytic triad/active site geometry.

    Computes per-pair distances and simultaneous contact fraction for the
    catalytic triad defined in the comparison.yaml file.

    The key metric is "simultaneous contact fraction" - the percentage of
    frames where ALL pairs are below the contact threshold at the same time.

    \b
    Examples:
        # Analyze all conditions in comparison.yaml
        polyzymd analyze triad -c comparison.yaml

        # Analyze specific condition
        polyzymd analyze triad -c comparison.yaml --condition "No Polymer"

        # Override replicates and equilibration time
        polyzymd analyze triad -c comparison.yaml --condition "No Polymer" \\
            -r 1-3 --eq-time 100ns
    """
    require_analysis_deps()

    from polyzymd.analysis.triad import CatalyticTriadAnalyzer
    from polyzymd.compare.config import ComparisonConfig
    from polyzymd.config.loader import load_config

    # Load comparison config
    click.echo(f"Loading comparison config from: {comparison}")
    try:
        comp_config = ComparisonConfig.from_yaml(comparison)
    except Exception as e:
        click.echo(click.style(f"Error loading comparison config: {e}", fg="red"), err=True)
        sys.exit(1)

    # Check catalytic triad is defined
    if comp_config.catalytic_triad is None:
        click.echo(
            click.style(
                "Error: No catalytic_triad section found in comparison.yaml.\n"
                "Add a catalytic_triad section with pairs to analyze.\n"
                "See: polyzymd compare init --help",
                fg="red",
            ),
            err=True,
        )
        sys.exit(1)

    triad_config = comp_config.catalytic_triad
    click.echo(f"Triad: {triad_config.name}")
    if triad_config.description:
        click.echo(f"  Description: {triad_config.description}")
    click.echo(f"  Pairs: {triad_config.n_pairs}")
    for pair in triad_config.pairs:
        click.echo(f"    - {pair.label}")
    click.echo(f"  Threshold: {triad_config.threshold} Å")

    # Determine equilibration time
    if eq_time is None:
        eq_time = comp_config.defaults.equilibration_time
    click.echo(f"  Equilibration: {eq_time}")
    if eq_time in ("0ns", "0ps", "0us"):
        click.echo(
            click.style(
                "  Warning: No equilibration time specified. Analysis will include "
                "potentially non-equilibrated frames. Consider using --eq-time.",
                fg="yellow",
            ),
            err=True,
        )

    # Select conditions to analyze
    if condition:
        try:
            conditions = [comp_config.get_condition(condition)]
        except KeyError as e:
            click.echo(click.style(str(e), fg="red"), err=True)
            sys.exit(1)
    else:
        conditions = comp_config.conditions

    output_path = Path(output_dir) if output_dir else None

    # Analyze each condition
    for cond in conditions:
        click.echo()
        click.echo(click.style(f"=== {cond.label} ===", fg="cyan", bold=True))

        # Load simulation config
        try:
            sim_config = load_config(cond.config)
        except Exception as e:
            click.echo(click.style(f"Error loading {cond.config}: {e}", fg="red"), err=True)
            continue

        # Determine replicates
        if replicates:
            rep_list = parse_replicates(replicates)
        else:
            rep_list = cond.replicates

        click.echo(f"  Replicates: {rep_list}")

        try:
            # Create analyzer
            analyzer = CatalyticTriadAnalyzer(
                config=sim_config,
                triad_config=triad_config,
                equilibration=eq_time,
            )

            if len(rep_list) == 1:
                # Single replicate
                result = analyzer.compute(
                    replicate=rep_list[0],
                    save=True,
                    output_dir=output_path,
                    recompute=recompute,
                )
                click.echo()
                click.echo(click.style("Triad Analysis Complete", fg="green"))
                for pr in result.pair_results:
                    frac_str = ""
                    if pr.fraction_below_threshold is not None:
                        frac_str = f" ({pr.fraction_below_threshold * 100:.1f}% below threshold)"
                    click.echo(f"  {pr.pair_label}: {pr.mean_distance:.2f} Å{frac_str}")
                click.echo()
                click.echo(
                    f"  Simultaneous contact: {result.simultaneous_contact_fraction * 100:.1f}%"
                )
                if result.sim_contact_sem is not None:
                    click.echo(f"    (SEM: ±{result.sim_contact_sem * 100:.1f}%)")

            else:
                # Multiple replicates - compute aggregated
                result = analyzer.compute_aggregated(
                    replicates=rep_list,
                    save=True,
                    output_dir=output_path,
                    recompute=recompute,
                )
                click.echo()
                click.echo(click.style("Triad Analysis Complete (Aggregated)", fg="green"))
                for pr in result.pair_results:
                    frac_str = ""
                    if pr.overall_fraction_below is not None:
                        sem_f = pr.sem_fraction_below or 0
                        frac_str = (
                            f" ({pr.overall_fraction_below * 100:.1f} ± {sem_f * 100:.1f}% below)"
                        )
                    click.echo(
                        f"  {pr.pair_label}: {pr.overall_mean:.2f} ± {pr.overall_sem:.2f} Å{frac_str}"
                    )
                click.echo()
                click.echo(
                    f"  Simultaneous contact: {result.overall_simultaneous_contact * 100:.1f} "
                    f"± {result.sem_simultaneous_contact * 100:.1f}%"
                )
                click.echo(f"  Per-replicate:")
                for rep, frac in zip(result.replicates, result.per_replicate_simultaneous):
                    click.echo(f"    Rep {rep}: {frac * 100:.1f}%")

        except Exception as e:
            click.echo(click.style(f"Analysis failed: {e}", fg="red"), err=True)
            if LOGGER.level == logging.DEBUG:
                import traceback

                traceback.print_exc()
            continue


# =============================================================================
# Plot Command Group
# =============================================================================


@click.group()
def plot() -> None:
    """Plot and compare analysis results.

    Load saved analysis results (JSON files) and create publication-quality
    figures. Useful for comparing results across different conditions.

    \b
    Examples:
        polyzymd plot rmsf --inputs result1.json result2.json --labels "Apo" "Holo"
        polyzymd plot distances --inputs dist1.json dist2.json --compare histogram
    """
    pass


# -----------------------------------------------------------------------------
# RMSF Plotting
# -----------------------------------------------------------------------------


@plot.command("rmsf")
@click.option(
    "--inputs",
    "-i",
    multiple=True,
    required=True,
    type=click.Path(exists=True),
    help="Result JSON files to plot (can specify multiple)",
)
@click.option(
    "--labels",
    "-l",
    multiple=True,
    help="Labels for each input (one per input file)",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    default=None,
    help="Output image path (default: rmsf_comparison.png)",
)
@click.option(
    "--title",
    default=None,
    help="Plot title",
)
@click.option(
    "--dpi",
    type=int,
    default=150,
    help="Resolution in DPI [default: 150]",
)
@click.option(
    "--show-error/--no-show-error",
    default=True,
    help="Show error bars/bands [default: show]",
)
def plot_rmsf_cmd(
    inputs: tuple[str, ...],
    labels: tuple[str, ...],
    output: Optional[str],
    title: Optional[str],
    dpi: int,
    show_error: bool,
) -> None:
    """Plot and compare RMSF results.

    Load one or more RMSF result JSON files and create comparison plots.

    \b
    Examples:
        # Single result
        polyzymd plot rmsf -i rmsf_aggregated.json

        # Compare multiple conditions
        polyzymd plot rmsf -i apo.json holo.json -l "Apo" "Holo" -o comparison.png
    """
    require_matplotlib()

    from polyzymd.analysis.results.rmsf import RMSFAggregatedResult, RMSFResult

    # Load results
    results = []
    for input_path in inputs:
        path = Path(input_path)
        click.echo(f"Loading: {path}")
        try:
            # Try loading as aggregated first
            result = RMSFAggregatedResult.load(path)
        except Exception:
            # Fall back to single result
            result = RMSFResult.load(path)
        results.append(result)

    # Validate labels
    if labels:
        if len(labels) != len(results):
            click.echo(
                click.style(
                    f"Error: Number of labels ({len(labels)}) must match "
                    f"number of inputs ({len(results)})",
                    fg="red",
                ),
                err=True,
            )
            sys.exit(1)
        label_list = list(labels)
    else:
        label_list = [Path(p).stem for p in inputs]

    # Create plot
    if len(results) == 1:
        from polyzymd.analysis.rmsf import plot_rmsf

        fig, ax = plot_rmsf(results[0], show_error=show_error, label=label_list[0])
    else:
        from polyzymd.analysis.rmsf import plot_rmsf_comparison

        fig, ax = plot_rmsf_comparison(results, labels=label_list, show_error=show_error)

    if title:
        ax.set_title(title)

    # Save or show
    output_path = Path(output) if output else Path("rmsf_comparison.png")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    click.echo(click.style(f"Plot saved: {output_path}", fg="green"))


# -----------------------------------------------------------------------------
# Distance Plotting
# -----------------------------------------------------------------------------


@plot.command("distances")
@click.option(
    "--inputs",
    "-i",
    multiple=True,
    required=True,
    type=click.Path(exists=True),
    help="Result JSON files to plot (can specify multiple)",
)
@click.option(
    "--labels",
    "-l",
    multiple=True,
    help="Labels for each input (one per input file)",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    default=None,
    help="Output image path (default: distance_comparison.png)",
)
@click.option(
    "--title",
    default=None,
    help="Plot title",
)
@click.option(
    "--dpi",
    type=int,
    default=150,
    help="Resolution in DPI [default: 150]",
)
@click.option(
    "--threshold",
    type=float,
    default=None,
    help="Distance threshold to mark on plot (Angstroms)",
)
@click.option(
    "--compare",
    type=click.Choice(["histogram", "timeseries", "bar"]),
    default="histogram",
    help="Comparison plot type [default: histogram]",
)
def plot_distances_cmd(
    inputs: tuple[str, ...],
    labels: tuple[str, ...],
    output: Optional[str],
    title: Optional[str],
    dpi: int,
    threshold: Optional[float],
    compare: str,
) -> None:
    """Plot and compare distance analysis results.

    Load one or more distance result JSON files and create comparison plots.

    \b
    Examples:
        # Single result histogram
        polyzymd plot distances -i distances.json --threshold 3.5

        # Compare conditions
        polyzymd plot distances -i apo.json holo.json -l "Apo" "Holo"
    """
    require_matplotlib()

    from polyzymd.analysis.results.distances import DistanceAggregatedResult, DistanceResult

    # Load results
    results = []
    for input_path in inputs:
        path = Path(input_path)
        click.echo(f"Loading: {path}")
        try:
            result = DistanceAggregatedResult.load(path)
        except Exception:
            result = DistanceResult.load(path)
        results.append(result)

    # Validate labels
    if labels:
        if len(labels) != len(results):
            click.echo(
                click.style(
                    f"Error: Number of labels ({len(labels)}) must match "
                    f"number of inputs ({len(results)})",
                    fg="red",
                ),
                err=True,
            )
            sys.exit(1)
        label_list = list(labels)
    else:
        label_list = [Path(p).stem for p in inputs]

    # Create plot based on comparison type
    output_path = Path(output) if output else Path("distance_comparison.png")

    if len(results) == 1:
        # Single result - plot histogram for each pair
        from polyzymd.analysis.distances import plot_distance_histogram

        result = results[0]
        if hasattr(result, "pair_results"):
            # Has multiple pairs
            for pr in result.pair_results:
                fig, ax = plot_distance_histogram(pr, threshold=threshold)
                if title:
                    ax.set_title(title)
                pair_output = (
                    output_path.parent / f"{output_path.stem}_{pr.pair_label}{output_path.suffix}"
                )
                pair_output.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(pair_output, dpi=dpi, bbox_inches="tight")
                click.echo(click.style(f"Plot saved: {pair_output}", fg="green"))
        else:
            fig, ax = plot_distance_histogram(result, threshold=threshold)
            if title:
                ax.set_title(title)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
            click.echo(click.style(f"Plot saved: {output_path}", fg="green"))

    else:
        # Multiple results - comparison plot
        from polyzymd.analysis.distances import plot_distance_comparison

        fig, ax = plot_distance_comparison(results, labels=label_list, threshold=threshold)
        if title:
            ax.set_title(title)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        click.echo(click.style(f"Plot saved: {output_path}", fg="green"))
