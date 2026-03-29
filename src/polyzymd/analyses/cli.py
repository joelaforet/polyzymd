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

import click

from polyzymd.cli.colors import echo_logo
from polyzymd.core.experimental import echo_experimental_warning

LOGGER = logging.getLogger("polyzymd.analyses")


def _echo_branding() -> None:
    """Print the PolyzyMD ASCII logo for top-level analysis commands."""
    echo_logo()


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

    \b
    Note: Analysis currently supports OpenMM-produced trajectories only
    (DCD format in PolyzyMD's standard directory layout). GROMACS XTC
    trajectory support is planned for v1.2.1 (see GitHub issue #47).
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
    from polyzymd.analyses._analysis_config import generate_analysis_template

    _echo_branding()

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
        from polyzymd.analyses._analysis_config import AnalysisConfig

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
    output_dir: str | None,
    binding_preference: bool,
    enzyme_pdb: str | None,
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

    from polyzymd.analyses.contacts import ParallelContactAnalyzer
    from polyzymd.analyses.contacts._aggregator import aggregate_contact_results
    from polyzymd.analyses.shared.loader import TrajectoryLoader, parse_time_string, time_to_frame
    from polyzymd.analyses.shared.selectors import MDAnalysisSelector
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
        echo_experimental_warning(("contacts_binding_preference",))

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

            from polyzymd.analyses.shared.binding_preference import (
                aggregate_binding_preference,
                compute_binding_preference,
                extract_polymer_composition,
                resolve_protein_group_selections,
            )
            from polyzymd.analyses.shared.surface_exposure import SurfaceExposureFilter

            # Determine enzyme PDB path
            if enzyme_pdb:
                enzyme_pdb_path = Path(enzyme_pdb)
            else:
                enzyme_pdb_path = sim_config.enzyme.pdb_path
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

                # Extract polymer composition (once - same for all replicates)
                polymer_composition = extract_polymer_composition(universe, None)
                click.echo(
                    f"  Polymer composition: {polymer_composition.total_residues} residues, "
                    f"{polymer_composition.total_heavy_atoms} heavy atoms"
                )

                # Compute binding preference for each replicate
                for i, (result, universe) in enumerate(zip(results, universes)):
                    try:
                        bp_result = compute_binding_preference(
                            contact_result=result,
                            surface_exposure=surface_exposure,
                            protein_groups=protein_groups,
                            polymer_composition=polymer_composition,
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
            agg_result.save(output_file)
            click.echo(f"  Results saved: {output_file}")

            # Aggregate and save binding preference if computed
            if binding_pref_results:
                from polyzymd.analyses.shared.binding_preference import (
                    aggregate_binding_preference,
                )

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
