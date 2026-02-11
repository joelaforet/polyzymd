"""Command-line interface for PolyzyMD analysis module.

This module provides CLI commands for trajectory analysis:
- `polyzymd analyze rmsf` - RMSF analysis
- `polyzymd analyze distances` - Distance analysis
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
    Examples:
        polyzymd analyze rmsf -c config.yaml -r 1-5 --eq-time 100ns
        polyzymd analyze distances -c config.yaml -r 1 --pair "resid 77 and name OG : resid 133 and name NE2"
    """
    pass


# -----------------------------------------------------------------------------
# RMSF Analysis
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
            threshold=threshold,
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
