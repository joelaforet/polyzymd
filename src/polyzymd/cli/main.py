"""
PolyzyMD Command Line Interface.

This module provides the main CLI entry point for PolyzyMD, using Click
for argument parsing and command organization.

Usage:
    polyzymd --help
    polyzymd build --config simulation.yaml
    polyzymd submit --config simulation.yaml --replicates 1-5
    polyzymd run-segment --config simulation.yaml --replicate 1
    polyzymd run-gromacs --config simulation.yaml --replicate 1
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

import click

from polyzymd.cli.colors import colored_echo, setup_colored_logging

# Bootstrap a minimal root handler so suppress_openff_logs() works at import
# time.  setup_colored_logging() replaces this handler when the CLI runs.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
LOGGER = logging.getLogger("polyzymd")


def suppress_openff_logs() -> None:
    """Suppress verbose OpenFF Interchange and Toolkit log messages.

    OpenFF libraries generate extensive INFO-level logs during system building,
    including per-atom messages for charge assignment and parameter collisions.
    For large systems (>10,000 atoms), this can produce millions of log lines.

    This function sets OpenFF-related loggers to WARNING level to suppress
    these messages while still showing important warnings and errors.
    """
    # Suppress OpenFF Interchange logs (e.g., "Preset charges applied to atom index X")
    logging.getLogger("openff.interchange").setLevel(logging.WARNING)

    # Suppress OpenFF Toolkit logs
    logging.getLogger("openff.toolkit").setLevel(logging.WARNING)

    # Suppress root logger INFO messages from OpenFF (e.g., "Key collision" messages)
    # These come from OpenFF but are logged to root logger
    # We set a filter to suppress INFO messages that look like OpenFF messages
    root_logger = logging.getLogger()

    class OpenFFFilter(logging.Filter):
        """Filter out verbose OpenFF messages from the root logger."""

        def filter(self, record: logging.LogRecord) -> bool:
            # Allow all non-INFO messages
            if record.levelno != logging.INFO:
                return True
            # Suppress known verbose OpenFF message patterns
            msg = record.getMessage()
            if "Key collision" in msg:
                return False
            if "associated_handler" in msg:
                return False
            return True

    root_logger.addFilter(OpenFFFilter())


def enable_openff_logs() -> None:
    """Re-enable OpenFF log messages for debugging.

    Sets OpenFF loggers back to INFO level and removes the root logger filter.
    """
    logging.getLogger("openff.interchange").setLevel(logging.INFO)
    logging.getLogger("openff.toolkit").setLevel(logging.INFO)

    # Remove OpenFF filter from root logger
    root_logger = logging.getLogger()
    for f in root_logger.filters[:]:
        if f.__class__.__name__ == "OpenFFFilter":
            root_logger.removeFilter(f)


# Suppress OpenFF logs by default
suppress_openff_logs()


@click.group()
@click.version_option(prog_name="polyzymd")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output")
@click.option(
    "--openff-logs",
    is_flag=True,
    help="Enable verbose OpenFF Interchange/Toolkit logs (suppressed by default)",
)
@click.option(
    "--no-color",
    is_flag=True,
    help="Disable colored output (also respects NO_COLOR env var)",
)
def cli(verbose: bool, openff_logs: bool, no_color: bool) -> None:
    """PolyzyMD: MD simulations for enzyme-polymer systems.

    A toolkit for building, running, and analyzing molecular dynamics
    simulations of enzymes with co-polymers.
    """
    # Replace the bootstrap handler with the colored formatter
    setup_colored_logging(verbose=verbose, no_color=no_color)

    if openff_logs:
        enable_openff_logs()
        LOGGER.info("OpenFF verbose logging enabled")


# =============================================================================
# Build Command
# =============================================================================


@cli.command()
@click.option(
    "-c",
    "--config",
    required=True,
    type=click.Path(exists=True),
    help="Path to YAML configuration file",
)
@click.option(
    "-r",
    "--replicate",
    default=1,
    type=int,
    help="Replicate number (default: 1)",
)
@click.option(
    "-o",
    "--output-dir",
    default=None,
    type=click.Path(),
    help="Output directory (default: from config). Alias for --scratch-dir.",
)
@click.option(
    "--scratch-dir",
    default=None,
    type=click.Path(),
    help="Scratch directory for simulation output (trajectories, checkpoints)",
)
@click.option(
    "--projects-dir",
    default=None,
    type=click.Path(),
    help="Projects directory for scripts and logs (long-term storage)",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Validate config without building",
)
@click.option(
    "--gromacs",
    is_flag=True,
    help="Export to GROMACS format (.gro, .top, .mdp) instead of preparing for OpenMM",
)
def build(
    config: str,
    replicate: int,
    output_dir: Optional[str],
    scratch_dir: Optional[str],
    projects_dir: Optional[str],
    dry_run: bool,
    gromacs: bool,
) -> None:
    """Build a simulation system from configuration.

    Loads the YAML configuration, constructs the molecular system
    (enzyme, substrate, polymers, solvent), and prepares it for simulation.

    By default, prepares the system for OpenMM simulation. Use --gromacs to
    export GROMACS-compatible files instead (.gro, .top, .mdp).

    GROMACS Export Notes:
        - Output files are placed in {projects_dir}/{replicate}/gromacs/
        - Filenames are derived from config: {enzyme_name}_{polymer_prefix}.*
        - The .mdp file is a stub for single-point energy; modify for production
        - Topology is split into .itp files for cleaner multi-component systems
    """
    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.config.schema import SimulationConfig

    colored_echo(f"Loading configuration from: {config}", phase="build")

    try:
        sim_config = SimulationConfig.from_yaml(config)
        colored_echo(f"Configuration validated: {sim_config.name}", phase="build")

        # Override directories if provided via CLI
        effective_scratch = scratch_dir or output_dir  # output_dir is alias for scratch_dir
        if effective_scratch:
            sim_config.output.scratch_directory = Path(effective_scratch)
        if projects_dir:
            sim_config.output.projects_directory = Path(projects_dir)

        if dry_run:
            colored_echo("Dry run - configuration is valid", phase="build")
            colored_echo(f"  Enzyme: {sim_config.enzyme.name}", phase="build")
            if sim_config.substrate:
                colored_echo(f"  Substrate: {sim_config.substrate.name}", phase="build")
            if sim_config.polymers and sim_config.polymers.enabled:
                colored_echo(f"  Polymers: {sim_config.polymers.type_prefix}", phase="build")
                colored_echo(f"  Polymer count: {sim_config.polymers.count}", phase="build")
            colored_echo(f"  Temperature: {sim_config.thermodynamics.temperature} K", phase="build")
            colored_echo(
                f"  Production time: {sim_config.simulation_phases.production.duration} ns",
                phase="build",
            )
            colored_echo(phase="build")
            colored_echo("Directories:", phase="build")
            colored_echo(f"  Projects: {sim_config.output.projects_directory}", phase="build")
            colored_echo(
                f"  Scratch: {sim_config.output.effective_scratch_directory}", phase="build"
            )
            if gromacs:
                colored_echo(phase="build")
                colored_echo("GROMACS export enabled:", phase="build")
                colored_echo("  Output: {projects_dir}/{replicate}/gromacs/", phase="build")
            return

        colored_echo(f"Building system for replicate {replicate}...", phase="build")
        working_dir = sim_config.get_working_directory(replicate)
        builder = SystemBuilder.from_config(sim_config)
        interchange = builder.build_from_config(
            config=sim_config,
            working_dir=working_dir,
            polymer_seed=replicate,
        )

        # Branch based on export format
        if gromacs:
            # Export to GROMACS format
            colored_echo("Exporting to GROMACS format...", phase="export")
            gromacs_dir = (
                sim_config.output.projects_directory / f"replicate_{replicate}" / "gromacs"
            )
            export_result = builder.export_to_gromacs(gromacs_dir)

            colored_echo("GROMACS export successful!", phase="export")
            colored_echo(f"Output directory: {gromacs_dir}", phase="export")
            colored_echo("Files generated:", phase="export")
            colored_echo(f"  - {export_result['gro'].name} (coordinates)", phase="export")
            colored_echo(f"  - {export_result['top'].name} (topology)", phase="export")
            colored_echo(
                f"  - {export_result['em_mdp'].name} (energy minimization)", phase="export"
            )
            for eq_mdp in export_result.get("eq_mdps", []):
                colored_echo(f"  - {eq_mdp.name} (equilibration)", phase="export")
            colored_echo(f"  - {export_result['prod_mdp'].name} (production)", phase="export")
            if export_result.get("posres_defines"):
                colored_echo("Position restraints added to molecule ITP files:", phase="export")
                for component, define in export_result["posres_defines"].items():
                    colored_echo(f"  - {component}: #ifdef {define}", phase="export")
            colored_echo(f"  - {export_result['run_script'].name} (run script)", phase="export")
            colored_echo(phase="export")
            colored_echo(
                f"To run: cd {gromacs_dir} && ./{export_result['run_script'].name}", phase="export"
            )

        else:
            # Default: prepare for OpenMM simulation
            colored_echo("Extracting OpenMM components...", phase="build")
            omm_topology, omm_system, omm_positions = builder.get_openmm_components()

            # Apply restraints if configured
            if sim_config.restraints:
                from polyzymd.core.restraints import RestraintFactory, apply_restraints

                colored_echo(
                    f"Applying {len(sim_config.restraints)} restraint(s)...", phase="build"
                )
                restraint_defs = []
                for r in sim_config.restraints:
                    if not r.enabled:
                        colored_echo(f"  - {r.name}: DISABLED (skipping)", phase="build")
                        continue

                    # Create restraint definition from config
                    restraint_def = RestraintFactory.from_config(r.model_dump())

                    # Validate the selection resolves to exactly one atom each
                    try:
                        indices1 = restraint_def.atom1.resolve(omm_topology)
                        indices2 = restraint_def.atom2.resolve(omm_topology)

                        if len(indices1) != 1:
                            colored_echo(
                                f"Error: Restraint '{r.name}' atom1 selection matched "
                                f"{len(indices1)} atoms (need exactly 1)",
                                err=True,
                                level=logging.ERROR,
                            )
                            sys.exit(1)
                        if len(indices2) != 1:
                            colored_echo(
                                f"Error: Restraint '{r.name}' atom2 selection matched "
                                f"{len(indices2)} atoms (need exactly 1)",
                                err=True,
                                level=logging.ERROR,
                            )
                            sys.exit(1)

                        colored_echo(
                            f"  - {r.name}: atom {indices1[0]} <-> atom {indices2[0]} "
                            f"(type={r.type.value}, d={r.distance} A, "
                            f"k={r.force_constant} kJ/mol/nm^2)",
                            phase="build",
                        )
                        restraint_defs.append(restraint_def)

                    except ValueError as e:
                        colored_echo(
                            f"Error: Restraint '{r.name}' invalid: {e}",
                            err=True,
                            level=logging.ERROR,
                        )
                        sys.exit(1)

                # Apply all validated restraints to the system
                if restraint_defs:
                    apply_restraints(restraint_defs, omm_topology, omm_system)
                    colored_echo(
                        f"Successfully applied {len(restraint_defs)} restraint(s)", phase="build"
                    )

            # Save OpenMM system to XML for --skip-build support
            from openmm import XmlSerializer

            system_xml_path = working_dir / "system.xml"
            colored_echo(f"Saving OpenMM system to {system_xml_path}...", phase="build")
            with open(system_xml_path, "w") as f:
                f.write(XmlSerializer.serialize(omm_system))

            colored_echo("System built successfully!", phase="build")
            colored_echo(f"Output directory: {working_dir}", phase="build")
            colored_echo("Files saved:", phase="build")
            colored_echo("  - solvated_system.pdb (topology + positions)", phase="build")
            colored_echo("  - system.xml (OpenMM system with restraints)", phase="build")
            colored_echo(
                "Use 'polyzymd run --skip-build' to run without rebuilding.", phase="build"
            )

    except FileNotFoundError as e:
        colored_echo(f"Error: {e}", err=True, level=logging.ERROR)
        sys.exit(1)
    except Exception as e:
        colored_echo(f"Build failed: {e}", err=True, level=logging.ERROR)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


# =============================================================================
# Run-GROMACS Command
# =============================================================================


@cli.command("run-gromacs")
@click.option(
    "-c",
    "--config",
    required=True,
    type=click.Path(exists=True),
    help="Path to YAML configuration file",
)
@click.option(
    "-r",
    "--replicate",
    default=1,
    type=int,
    help="Replicate number (default: 1)",
)
@click.option(
    "--scratch-dir",
    default=None,
    type=click.Path(),
    help="Scratch directory for simulation output (trajectories, checkpoints)",
)
@click.option(
    "--projects-dir",
    default=None,
    type=click.Path(),
    help="Projects directory for scripts and logs (long-term storage)",
)
@click.option(
    "--gmx-path",
    default="gmx",
    help="Path to GROMACS executable (default: gmx)",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Export files but don't run simulation",
)
def run_gromacs(
    config: str,
    replicate: int,
    scratch_dir: Optional[str],
    projects_dir: Optional[str],
    gmx_path: str,
    dry_run: bool,
) -> None:
    """Run a simulation using GROMACS.

    Builds the system, exports to GROMACS format (.gro, .top, .mdp),
    and executes the full GROMACS workflow locally (EM, equilibration,
    production, and trajectory post-processing).

    \b
    Notes:
        - Requires GROMACS to be installed and accessible
        - Use --gmx-path to specify a custom GROMACS executable
        - Use --dry-run to export files without running the simulation
    """
    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.config.schema import SimulationConfig

    colored_echo(f"Loading configuration from: {config}", phase="export")

    try:
        sim_config = SimulationConfig.from_yaml(config)
        colored_echo(f"Running GROMACS simulation: {sim_config.name}", phase="export")

        # Override directories if provided via CLI
        if scratch_dir:
            sim_config.output.scratch_directory = Path(scratch_dir)
        if projects_dir:
            sim_config.output.projects_directory = Path(projects_dir)

        _run_gromacs_impl(
            sim_config=sim_config,
            replicate=replicate,
            gmx_path=gmx_path,
            dry_run=dry_run,
        )

    except Exception as e:
        colored_echo(f"Simulation failed: {e}", err=True, level=logging.ERROR)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


def _run_gromacs_impl(
    sim_config: "SimulationConfig",
    replicate: int,
    gmx_path: str,
    dry_run: bool,
) -> None:
    """Run simulation using GROMACS.

    Parameters
    ----------
    sim_config : SimulationConfig
        Validated simulation configuration.
    replicate : int
        Replicate number.
    gmx_path : str
        Path to GROMACS executable.
    dry_run : bool
        If True, export files but don't run simulation.
    """
    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.exporters.gromacs import GromacsError, GromacsExporter, GromacsRunner

    # Determine output directory for GROMACS files
    gromacs_dir = sim_config.output.projects_directory / f"replicate_{replicate}" / "gromacs"
    working_dir = sim_config.get_working_directory(replicate)

    click.echo(f"Building system for replicate {replicate}...")
    builder = SystemBuilder.from_config(sim_config)
    interchange = builder.build_from_config(
        config=sim_config,
        working_dir=working_dir,
        polymer_seed=replicate,
    )

    # Get component info for position restraints
    component_info = builder.get_component_info()

    # Export to GROMACS format
    colored_echo("Exporting to GROMACS format...", phase="export")
    exporter = GromacsExporter(
        interchange=interchange,
        config=sim_config,
        component_info=component_info,
    )
    export_result = exporter.export(
        output_dir=gromacs_dir,
        gmx_command=gmx_path,
    )

    colored_echo(f"\nGROMACS files exported to: {gromacs_dir}", phase="export")
    colored_echo("Files generated:", phase="export")
    colored_echo(f"  - {export_result['gro'].name} (coordinates)", phase="export")
    colored_echo(f"  - {export_result['top'].name} (topology)", phase="export")
    colored_echo(f"  - {export_result['em_mdp'].name} (energy minimization)", phase="export")
    for eq_mdp in export_result["eq_mdps"]:
        colored_echo(f"  - {eq_mdp.name} (equilibration)", phase="export")
    colored_echo(f"  - {export_result['prod_mdp'].name} (production)", phase="export")
    if export_result.get("posres_defines"):
        colored_echo("Position restraints added to molecule ITP files:", phase="export")
        for component, define in export_result["posres_defines"].items():
            colored_echo(f"  - {component}: #ifdef {define}", phase="export")
    colored_echo(f"  - {export_result['run_script'].name} (run script)", phase="export")

    if dry_run:
        colored_echo(
            "\n--dry-run specified: Files exported but simulation not started.", phase="export"
        )
        colored_echo(
            f"To run manually: cd {gromacs_dir} && ./{export_result['run_script'].name}",
            phase="export",
        )
        return

    # Run GROMACS workflow
    colored_echo("\nStarting GROMACS simulation...", phase="export")
    colored_echo(f"Using GROMACS executable: {gmx_path}", phase="export")

    # Get equilibration MDP filenames
    eq_mdp_names = [p.name for p in export_result["eq_mdps"]]

    # Generate prefix from export result
    prefix = export_result["gro"].stem  # e.g., "lysozyme_PEG" from "lysozyme_PEG.gro"

    try:
        runner = GromacsRunner(
            working_dir=gromacs_dir,
            prefix=prefix,
            equilibration_mdps=eq_mdp_names,
            gmx_command=gmx_path,
        )
        runner.run_full_workflow()

        colored_echo("\nGROMACS simulation completed successfully!", phase="export")
        colored_echo(f"Output directory: {gromacs_dir}", phase="export")

    except GromacsError as e:
        colored_echo(f"\nGROMACS simulation failed: {e}", err=True, level=logging.ERROR)
        colored_echo(f"Check log files in: {gromacs_dir}", err=True, level=logging.ERROR)
        sys.exit(1)

    except FileNotFoundError as e:
        colored_echo(f"\nError: {e}", err=True, level=logging.ERROR)
        colored_echo(
            "Ensure GROMACS is installed and in your PATH, or use --gmx-path.",
            err=True,
            level=logging.ERROR,
        )
        sys.exit(1)


# =============================================================================
# Submit Command (SLURM)
# =============================================================================


@cli.command()
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
    help="Replicate range (e.g., '1-5', '1,3,5')",
)
@click.option(
    "--scratch-dir",
    default=None,
    type=click.Path(),
    help="Scratch directory for simulation output (high-performance storage)",
)
@click.option(
    "--projects-dir",
    default=None,
    type=click.Path(),
    help="Projects directory for scripts/logs (long-term storage)",
)
@click.option(
    "--preset",
    type=click.Choice(["aa100", "al40", "blanca-shirts", "bridges2", "testing"]),
    default="aa100",
    help="SLURM partition preset (default: aa100)",
)
@click.option(
    "--email",
    default="",
    help="Email for job notifications",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Generate scripts without submitting",
)
@click.option(
    "--output-dir",
    default=None,
    help="Directory for job scripts (default: {projects_dir}/job_scripts)",
)
@click.option(
    "--time-limit",
    default=None,
    help="Override SLURM time limit (format: HH:MM:SS or M:SS)",
)
@click.option(
    "--memory",
    default=None,
    help="Override SLURM memory allocation (e.g. '4G', '8G'). Not needed for bridges2 (allocated per GPU).",
)
@click.option(
    "--account",
    default=None,
    help=(
        "Override SLURM account / allocation ID. Required for bridges2 "
        "(find yours at https://www.psc.edu/resources/bridges-2/user-guide)."
    ),
)
@click.option(
    "--gpu-type",
    default=None,
    type=click.Choice(["v100-16", "v100-32", "l40s-48", "h100-80"]),
    help=(
        "Override GPU type for presets that use the --gpus directive (e.g. bridges2). "
        "Valid types: v100-16, v100-32, l40s-48, h100-80. Default for bridges2: v100-32."
    ),
)
@click.option(
    "--openff-logs",
    "submit_openff_logs",
    is_flag=True,
    help="Enable verbose OpenFF logs in generated job scripts (suppressed by default)",
)
@click.option(
    "--skip-build",
    is_flag=True,
    help="Skip system building in generated jobs (use pre-built system from 'polyzymd build')",
)
@click.option(
    "--pixi-env",
    default=None,
    type=click.Choice(["cuda-12-4", "cuda-12-6"]),
    help=(
        "Pixi environment for SLURM jobs. If omitted, inferred from --preset "
        "(blanca/alpine presets → cuda-12-4, bridges2 → cuda-12-6)."
    ),
)
def submit(
    config: str,
    replicates: str,
    scratch_dir: Optional[str],
    projects_dir: Optional[str],
    preset: str,
    email: str,
    dry_run: bool,
    output_dir: Optional[str],
    time_limit: Optional[str],
    memory: Optional[str],
    account: Optional[str],
    gpu_type: Optional[str],
    submit_openff_logs: bool,
    skip_build: bool,
    pixi_env: Optional[str],
) -> None:
    """Submit simulation jobs to SLURM.

    Creates and submits self-resubmitting jobs (one per replicate)
    that automatically checkpoint and resume until the full
    production duration is complete.

    Directory structure:
    - projects_dir: Where job scripts and SLURM logs are stored (long-term storage)
    - scratch_dir: Where simulation data is written (high-performance storage)
    """
    from polyzymd.workflow.daisy_chain import submit_daisy_chain
    from polyzymd.workflow.slurm import PRESET_DEFAULT_PIXI_ENV

    # Resolve pixi environment: explicit flag > preset default
    resolved_pixi_env = pixi_env or PRESET_DEFAULT_PIXI_ENV.get(preset, "cuda-12-4")

    colored_echo(f"Loading configuration from: {config}", phase="workflow")
    colored_echo(f"Submitting jobs with preset: {preset}", phase="workflow")
    colored_echo(f"Pixi environment: {resolved_pixi_env}", phase="workflow")
    colored_echo(f"Replicates: {replicates}", phase="workflow")
    if scratch_dir:
        colored_echo(f"Scratch directory: {scratch_dir}", phase="workflow")
    if projects_dir:
        colored_echo(f"Projects directory: {projects_dir}", phase="workflow")
    if account:
        colored_echo(f"Account: {account}", phase="workflow")
    if memory:
        colored_echo(f"Memory allocation: {memory}", phase="workflow")
    if gpu_type:
        colored_echo(f"GPU type override: {gpu_type}", phase="workflow")
    if skip_build:
        colored_echo("Skip-build mode: using pre-built systems", phase="workflow")

    if dry_run:
        colored_echo("DRY RUN MODE - scripts will be created but not submitted", phase="workflow")

    try:
        results = submit_daisy_chain(
            config_path=config,
            slurm_preset=preset,
            replicates=replicates,
            email=email,
            dry_run=dry_run,
            pixi_env=resolved_pixi_env,
            output_dir=output_dir,
            scratch_dir=scratch_dir,
            projects_dir=projects_dir,
            time_limit=time_limit,
            memory=memory,
            account=account,
            gpu_type=gpu_type,
            openff_logs=submit_openff_logs,
            skip_build=skip_build,
        )

        if not dry_run:
            colored_echo("\nJob submission complete!", phase="workflow")
            colored_echo("Monitor with: squeue -u $USER", phase="workflow")

    except Exception as e:
        colored_echo(f"Submission failed: {e}", err=True, level=logging.ERROR)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


# =============================================================================
# Run-segment Command — unified entry point for self-resubmitting jobs
# =============================================================================


@cli.command("run-segment")
@click.option(
    "-c",
    "--config",
    required=True,
    type=click.Path(exists=True),
    help="Path to YAML configuration file",
)
@click.option(
    "-r",
    "--replicate",
    default=1,
    type=int,
    help="Replicate number (default: 1)",
)
@click.option(
    "--scratch-dir",
    default=None,
    type=click.Path(),
    help="Override scratch directory for simulation output",
)
@click.option(
    "--skip-build",
    is_flag=True,
    help="Skip system building (use existing) for initial segment",
)
def run_segment(
    config: str,
    replicate: int,
    scratch_dir: Optional[str],
    skip_build: bool,
) -> None:
    """Run the next simulation segment (self-resubmitting job entry point).

    This is the unified command called by every SLURM job in the
    self-resubmitting architecture. It determines what work remains
    by loading progress state, then either:

    \b
    - Builds, equilibrates, and runs the first production segment (segment 0)
    - Continues from the last completed segment
    - Exits cleanly if the simulation is already complete

    Exit codes:
        0  - Segment completed (check progress to decide resubmission)
        1  - Error
        99 - Graceful interruption (wall-time signal, should resubmit)
    """
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.simulation.progress import (
        SegmentStatus,
        SimulationProgress,
        SimulationStatus,
        get_next_segment_info,
        load_or_scan_progress,
        save_progress,
    )

    colored_echo(f"Loading configuration from: {config}", phase="simulation")

    try:
        sim_config = SimulationConfig.from_yaml(config)
    except Exception as e:
        colored_echo(f"Failed to load config: {e}", err=True, level=logging.ERROR)
        sys.exit(1)

    # Determine working directory
    if scratch_dir:
        working_dir = Path(scratch_dir)
    else:
        working_dir = sim_config.get_working_directory(replicate)

    working_dir.mkdir(parents=True, exist_ok=True)

    # Calculate total steps and samples from config
    prod = sim_config.simulation_phases.production
    timestep_fs = prod.time_step
    total_steps = int(prod.duration * 1e6 / timestep_fs)
    total_samples = prod.samples

    colored_echo(f"Working directory: {working_dir}", phase="simulation")
    colored_echo(f"Total production: {prod.duration} ns = {total_steps} steps", phase="simulation")

    # Load or create progress
    progress = load_or_scan_progress(
        working_dir=working_dir,
        config_path=str(Path(config).resolve()),
        total_steps=total_steps,
        total_samples=total_samples,
        timestep_fs=timestep_fs,
        replicate=replicate,
    )
    save_progress(working_dir, progress)

    # Check if simulation is already complete
    if progress.is_complete:
        colored_echo("Simulation already complete — nothing to do.", phase="simulation")
        sys.exit(0)

    # ---- Handle FAILED segments before determining next work ----
    # A FAILED segment (no state.xml, no INTERRUPTED, no checkpoint) means
    # no recoverable state exists.  Clean up the directory and remove the
    # record so the segment can be retried from the appropriate starting
    # point (initial build if segment 0, or continuation from the last good
    # segment otherwise).
    import shutil

    failed_segments = [s for s in progress.segments if s.status == SegmentStatus.FAILED]
    for failed in failed_segments:
        failed_dir = working_dir / f"production_{failed.index}"
        if failed_dir.exists():
            colored_echo(
                f"Cleaning up failed segment {failed.index} (no recoverable state) — will retry",
                phase="simulation",
                level=logging.WARNING,
            )
            shutil.rmtree(failed_dir)
        progress.segments = [s for s in progress.segments if s.index != failed.index]
    if failed_segments:
        save_progress(working_dir, progress)

    # Determine what to run next
    seg_info = get_next_segment_info(progress, total_steps, total_samples)
    if seg_info is None:
        colored_echo("No remaining work — simulation complete.", phase="simulation")
        sys.exit(0)

    seg_idx = seg_info["segment_index"]
    steps_to_run = seg_info["steps_to_run"]
    samples_to_write = seg_info["samples_to_write"]
    report_interval = seg_info["report_interval"]
    duration_ns = (steps_to_run * timestep_fs) / 1e6

    colored_echo(
        f"Next segment: {seg_idx} "
        f"({duration_ns:.3f} ns, {steps_to_run} steps, {samples_to_write} frames)",
        phase="simulation",
    )

    try:
        if seg_idx == 0 and not progress.segments:
            # ---- FIRST RUN: build, equilibrate, run segment 0 ----
            _run_initial_segment(
                sim_config=sim_config,
                working_dir=working_dir,
                replicate=replicate,
                skip_build=skip_build,
                duration_ns=duration_ns,
                num_samples=samples_to_write,
                timestep_fs=timestep_fs,
                report_interval=report_interval,
                checkpoint_interval_s=prod.checkpoint_interval,
            )
        else:
            # ---- CONTINUATION: load previous state, run next segment ----
            _run_continuation_segment(
                working_dir=working_dir,
                segment_index=seg_idx,
                duration_ns=duration_ns,
                num_samples=samples_to_write,
                timestep_fs=timestep_fs,
                report_interval=report_interval,
                checkpoint_interval_s=prod.checkpoint_interval,
            )

        colored_echo(f"Segment {seg_idx} completed successfully.", phase="simulation")

    except Exception as e:
        # Check if this is a GracefulExit (signal-based interruption)
        from polyzymd.simulation.signals import EXIT_CODE_INTERRUPTED, GracefulExit

        if isinstance(e, GracefulExit):
            colored_echo(
                f"Segment {seg_idx} interrupted (graceful shutdown): {e}",
                phase="simulation",
                level=logging.WARNING,
            )
            sys.exit(EXIT_CODE_INTERRUPTED)

        colored_echo(f"Segment {seg_idx} failed: {e}", err=True, level=logging.ERROR)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


def _run_initial_segment(
    sim_config: "SimulationConfig",
    working_dir: Path,
    replicate: int,
    skip_build: bool,
    duration_ns: float,
    num_samples: int,
    timestep_fs: float,
    report_interval: int | None = None,
    checkpoint_interval_s: float = 60.0,
) -> None:
    """Build system, equilibrate, and run the first production segment.

    Parameters
    ----------
    sim_config : SimulationConfig
        Validated simulation configuration.
    working_dir : Path
        Working directory for output.
    replicate : int
        Replicate number.
    skip_build : bool
        If True, load pre-built system from disk.
    duration_ns : float
        Production segment duration in nanoseconds.
    num_samples : int
        Number of trajectory frames to save.
    timestep_fs : float
        Integration timestep in femtoseconds.
    report_interval : int or None
        Fixed reporter interval in steps. Overrides per-segment
        interval calculation when provided.
    checkpoint_interval_s : float
        Wall-time interval in seconds between restart checkpoints.
    """
    from polyzymd.simulation.runner import SimulationRunner

    temperature = sim_config.thermodynamics.temperature
    pressure = sim_config.thermodynamics.pressure

    if not skip_build:
        from polyzymd.builders.system_builder import SystemBuilder

        colored_echo(f"Building system for replicate {replicate}...", phase="build")
        builder = SystemBuilder.from_config(sim_config)
        builder.build_from_config(
            config=sim_config,
            working_dir=working_dir,
            polymer_seed=replicate,
        )

        colored_echo("Extracting OpenMM components...", phase="build")
        omm_topology, omm_system, omm_positions = builder.get_openmm_components()

        # Apply restraints if configured
        if sim_config.restraints:
            from polyzymd.core.restraints import RestraintFactory, apply_restraints

            colored_echo(f"Applying {len(sim_config.restraints)} restraint(s)...", phase="build")
            restraint_defs = []
            for r in sim_config.restraints:
                if not r.enabled:
                    continue
                restraint_def = RestraintFactory.from_config(r.model_dump())
                restraint_defs.append(restraint_def)
            if restraint_defs:
                apply_restraints(restraint_defs, omm_topology, omm_system)
                colored_echo(f"Applied {len(restraint_defs)} restraint(s)", phase="build")
    else:
        from openmm import XmlSerializer
        from openmm.app import PDBFile

        colored_echo("Loading pre-built system...", phase="simulation")
        pdb_path = working_dir / "solvated_system.pdb"
        system_path = working_dir / "system.xml"
        if not pdb_path.exists() or not system_path.exists():
            raise FileNotFoundError(
                f"Pre-built system not found in {working_dir}. "
                "Run 'polyzymd build' first or remove --skip-build."
            )
        pdb = PDBFile(str(pdb_path))
        omm_topology = pdb.topology
        omm_positions = pdb.positions
        with open(system_path, "r") as f:
            omm_system = XmlSerializer.deserialize(f.read())

    # Create runner
    runner = SimulationRunner(
        topology=omm_topology,
        system=omm_system,
        positions=omm_positions,
        working_dir=working_dir,
    )

    # Minimize
    colored_echo("Running energy minimization...", phase="simulation")
    runner.minimize()

    # Equilibrate
    phases = sim_config.simulation_phases
    eq_duration = phases.total_equilibration_duration
    eq_mode = "multi-stage" if phases.uses_staged_equilibration else "simple"
    colored_echo(f"Running equilibration: {eq_duration:.3f} ns ({eq_mode})...", phase="simulation")
    eq_result = runner.run_equilibration(temperature=temperature, config=phases)

    # Save equilibration progress so a resubmitted job knows eq is done
    from polyzymd.simulation.progress import (
        EquilibrationStageRecord,
        SegmentStatus,
        load_progress,
        save_progress,
    )

    progress = load_progress(working_dir)
    if progress is not None:
        eq_stages = []
        if eq_result.get("type") == "staged_equilibration":
            for stage_info in eq_result.get("stages", []):
                eq_stages.append(
                    EquilibrationStageRecord(
                        index=stage_info["stage_index"],
                        name=stage_info["stage_name"],
                        status=SegmentStatus.COMPLETED,
                        duration_ns=stage_info["duration_ns"],
                        ensemble=stage_info.get("ensemble", "NVT"),
                    )
                )
        progress.equilibration_stages = eq_stages
        save_progress(working_dir, progress)
        colored_echo(f"Saved equilibration progress ({len(eq_stages)} stages)", phase="simulation")

    # Run first production segment
    colored_echo(
        f"Running production segment 0: {duration_ns:.3f} ns, {num_samples} frames...",
        phase="simulation",
    )
    runner.run_production(
        temperature=temperature,
        duration_ns=duration_ns,
        num_samples=num_samples,
        timestep_fs=timestep_fs,
        pressure=pressure,
        segment_index=0,
        report_interval=report_interval,
        checkpoint_interval_s=checkpoint_interval_s,
    )


def _run_continuation_segment(
    working_dir: Path,
    segment_index: int,
    duration_ns: float,
    num_samples: int,
    timestep_fs: float,
    report_interval: int | None = None,
    checkpoint_interval_s: float = 60.0,
) -> None:
    """Continue from the last completed segment.

    Parameters
    ----------
    working_dir : Path
        Working directory containing simulation outputs.
    segment_index : int
        Segment index to run.
    duration_ns : float
        Duration of this segment in nanoseconds.
    num_samples : int
        Number of trajectory frames to save.
    timestep_fs : float
        Integration timestep in femtoseconds.
    report_interval : int or None
        Fixed reporter interval in steps. Overrides per-segment
        interval calculation when provided.
    checkpoint_interval_s : float
        Wall-time interval in seconds between restart checkpoints.
    """
    from polyzymd.simulation.continuation import ContinuationManager

    colored_echo(f"Loading state from segment {segment_index - 1}...", phase="simulation")
    manager = ContinuationManager(
        working_dir=working_dir,
        segment_index=segment_index,
    )
    manager.load_previous_state()

    colored_echo(
        f"Running segment {segment_index}: {duration_ns:.3f} ns, {num_samples} frames...",
        phase="simulation",
    )
    manager.run_segment(
        duration_ns=duration_ns,
        num_samples=num_samples,
        timestep_fs=timestep_fs,
        report_interval=report_interval,
        checkpoint_interval_s=checkpoint_interval_s,
    )


# =============================================================================
# Check-progress Command
# =============================================================================


@cli.command("check-progress")
@click.option(
    "-c",
    "--config",
    required=True,
    type=click.Path(exists=True),
    help="Path to YAML configuration file",
)
@click.option(
    "-r",
    "--replicate",
    default=1,
    type=int,
    help="Replicate number (default: 1)",
)
@click.option(
    "--scratch-dir",
    default=None,
    type=click.Path(),
    help="Override scratch directory",
)
def check_progress(
    config: str,
    replicate: int,
    scratch_dir: Optional[str],
) -> None:
    """Check whether a simulation is complete.

    Loads progress state and exits with code 0 if the simulation is
    complete, or code 1 if work remains. Used by SLURM resubmission
    logic to decide whether to resubmit.

    \b
    Exit codes:
        0 - Simulation complete (do NOT resubmit)
        1 - Work remains (resubmit)
    """
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.simulation.progress import load_or_scan_progress

    try:
        sim_config = SimulationConfig.from_yaml(config)
    except Exception as e:
        colored_echo(f"Failed to load config: {e}", err=True, level=logging.ERROR)
        sys.exit(1)

    if scratch_dir:
        working_dir = Path(scratch_dir)
    else:
        working_dir = sim_config.get_working_directory(replicate)

    prod = sim_config.simulation_phases.production
    timestep_fs = prod.time_step
    total_steps = int(prod.duration * 1e6 / timestep_fs)
    total_samples = prod.samples

    progress = load_or_scan_progress(
        working_dir=working_dir,
        config_path=str(Path(config).resolve()),
        total_steps=total_steps,
        total_samples=total_samples,
        timestep_fs=timestep_fs,
        replicate=replicate,
    )

    pct = progress.fraction_complete() * 100
    colored_echo(
        f"Progress: {progress.total_steps_completed}/{progress.total_steps_requested} steps "
        f"({pct:.1f}%), {len(progress.segments)} segment(s)",
        phase="progress",
    )

    if progress.is_complete:
        colored_echo("Status: COMPLETE", phase="progress")
        sys.exit(0)
    else:
        remaining_ns = (progress.steps_remaining * timestep_fs) / 1e6
        colored_echo(
            f"Status: {progress.status.value} — {remaining_ns:.3f} ns remaining",
            phase="progress",
        )
        sys.exit(1)


# =============================================================================
# Status Command
# =============================================================================


@cli.command()
@click.option(
    "-c",
    "--config",
    required=True,
    type=click.Path(exists=True),
    help="Path to YAML configuration file",
)
def status(config: str) -> None:
    """Show progress overview for all replicates.

    Auto-detects replicate directories and displays a compact progress
    summary with colored bars, completion percentage, ns progress, and
    simulation status for each replicate.

    \b
    Example:
        polyzymd status -c config.yaml
    """
    from polyzymd.cli.colors import render_progress_bar
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.simulation.progress import SimulationStatus, load_progress

    try:
        sim_config = SimulationConfig.from_yaml(config)
    except Exception as e:
        click.echo(click.style(f"Error: Failed to load config: {e}", fg="red"), err=True)
        sys.exit(1)

    # Total production time from config
    total_ns = sim_config.simulation_phases.production.duration

    # Build a human-readable system name from the directory template
    # (format with replicate=1, then strip the trailing "_run1")
    dir_name = sim_config._format_run_directory_name(1)
    system_name = dir_name.rsplit("_run", 1)[0] if "_run" in dir_name else dir_name

    # Discover replicate directories
    replicates = sim_config.discover_replicate_dirs()

    # Also include configured replicates that don't exist on disk yet
    # (they show as "not found")
    found_nums = {num for num, _ in replicates}
    rep_map: dict[int, Path | None] = dict(replicates)

    # If no replicates found on disk, show a message
    if not replicates:
        click.echo(click.style(f"  polyzymd status — {system_name}", bold=True))
        click.echo(f"  {'─' * 50}")
        click.echo()
        click.echo("  No replicate directories found.")
        click.echo(
            f"  Expected pattern: {sim_config._format_run_directory_name(1).rsplit('1', 1)[0]}*"
        )
        click.echo(f"  In: {sim_config.output.effective_scratch_directory}")
        sys.exit(0)

    # Header
    click.echo()
    click.echo(click.style(f"  polyzymd status — {system_name}", bold=True))
    click.echo(f"  {'─' * 50}")
    click.echo()

    # Determine the widest replicate label for alignment
    max_rep = max(found_nums)
    label_width = len(f"run{max_rep}")

    need_attention = 0

    for rep_num, rep_path in sorted(rep_map.items()):
        label = f"run{rep_num}"

        if rep_path is None:
            # Directory not found on disk
            frac = 0.0
            completed_ns = 0.0
            status_str = "not_found"
            status_display = "not found"
        else:
            progress = load_progress(rep_path)
            if progress is None:
                frac = 0.0
                completed_ns = 0.0
                status_str = "not_started"
                status_display = SimulationStatus.NOT_STARTED.value
            else:
                frac = progress.fraction_complete()
                # Compute ns from total steps (not time_completed_ns which
                # only counts COMPLETED segments, ignoring INTERRUPTED/RUNNING).
                completed_ns = (progress.total_steps_completed * progress.timestep_fs) / 1e6
                status_val = progress.status
                status_str = status_val.value
                status_display = status_val.value

        bar = render_progress_bar(frac, status_str)
        pct = frac * 100

        # Count replicates needing attention
        if status_str in ("interrupted", "failed", "not_started", "not_found"):
            need_attention += 1

        # Format: "  run1  ████░░░░  100.0%  100.0/100.0 ns  completed"
        click.echo(
            f"  {label:<{label_width}}  {bar}  {pct:5.1f}%  "
            f"{completed_ns:6.1f}/{total_ns:.1f} ns  {status_display}"
        )

    click.echo()
    total_reps = len(rep_map)
    if need_attention > 0:
        click.echo(
            f"  {need_attention}/{total_reps} need attention "
            f"(recover with: polyzymd recover -c {config} -r <N> --submit)"
        )
    else:
        click.echo(click.style(f"  All {total_reps} replicates completed!", fg="green"))
    click.echo()


# =============================================================================
# Validate Command
# =============================================================================


@cli.command()
@click.option(
    "-c",
    "--config",
    required=True,
    type=click.Path(exists=True),
    help="Path to YAML configuration file",
)
def validate(config: str) -> None:
    """Validate a configuration file.

    Checks that the configuration is valid and all referenced
    files exist.
    """
    from polyzymd.config.schema import SimulationConfig

    colored_echo(f"Validating configuration: {config}")

    try:
        sim_config = SimulationConfig.from_yaml(config)

        click.echo(click.style("Configuration is valid!", fg="green"))
        colored_echo()
        colored_echo("Summary:")
        colored_echo(f"  Name: {sim_config.name}")
        colored_echo(f"  Enzyme: {sim_config.enzyme.name}")

        if sim_config.substrate:
            colored_echo(f"  Substrate: {sim_config.substrate.name}")
        else:
            colored_echo("  Substrate: None (apo simulation)")

        if sim_config.polymers and sim_config.polymers.enabled:
            colored_echo(f"  Polymers: {sim_config.polymers.type_prefix}")
            colored_echo(f"    Count: {sim_config.polymers.count}")
            colored_echo(f"    Length: {sim_config.polymers.length}")
            for m in sim_config.polymers.monomers:
                colored_echo(f"    Monomer {m.label}: {m.probability * 100:.1f}%")
        else:
            colored_echo("  Polymers: Disabled")

        colored_echo(f"  Temperature: {sim_config.thermodynamics.temperature} K")
        colored_echo(f"  Pressure: {sim_config.thermodynamics.pressure} atm")
        colored_echo()
        colored_echo("Simulation phases:")
        eq = sim_config.simulation_phases.equilibration
        colored_echo(f"  Equilibration: {eq.duration} ns ({eq.ensemble.value})")
        prod = sim_config.simulation_phases.production
        colored_echo(f"  Production: {prod.duration} ns ({prod.ensemble.value})")

        if sim_config.restraints:
            colored_echo()
            colored_echo(f"Restraints: {len(sim_config.restraints)}")
            for r in sim_config.restraints:
                status = "enabled" if r.enabled else "disabled"
                colored_echo(f"  - {r.name} ({r.type.value}): {status}")

    except FileNotFoundError as e:
        click.echo(click.style(f"Error: {e}", fg="red"), err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Validation failed: {e}", fg="red"), err=True)
        sys.exit(1)


# =============================================================================
# Init Command
# =============================================================================


@cli.command()
@click.option(
    "-n",
    "--name",
    required=True,
    help="Name of the project directory to create",
)
def init(name: str) -> None:
    """Initialize a new PolyzyMD project directory.

    Creates a scaffold with a template configuration file and placeholder
    structure files to help new users get started.

    \b
    Example:
        polyzymd init --name my_simulation

    This creates:
        my_simulation/
        ├── config.yaml              <- Edit this file
        ├── structures/              <- Add your PDB/SDF files
        ├── job_scripts/
        └── slurm_logs/
    """
    import shutil
    from importlib import resources

    project_dir = Path(name)

    # Check if directory already exists
    if project_dir.exists():
        click.echo(
            click.style(f"Error: Directory '{name}' already exists.", fg="red"),
            err=True,
        )
        colored_echo("Choose a different name or remove the existing directory.")
        sys.exit(1)

    try:
        # Create directory structure
        colored_echo(f"Creating project directory: {name}/")
        project_dir.mkdir(parents=True)
        (project_dir / "structures").mkdir()
        (project_dir / "job_scripts").mkdir()
        (project_dir / "slurm_logs").mkdir()

        # Copy template configuration
        template_path = resources.files("polyzymd.configs.templates").joinpath(
            "config_template.yaml"
        )
        config_dest = project_dir / "config.yaml"
        with resources.as_file(template_path) as template_file:
            shutil.copy(template_file, config_dest)

        # Create placeholder files
        protein_placeholder = project_dir / "structures" / "place_protein_here.placeholder.txt"
        protein_placeholder.write_text("""\
# ============================================================================
# PLACEHOLDER: Place your protein PDB file here
# ============================================================================
#
# This directory should contain your input structure files.
#
# PROTEIN PDB FILE (required):
#   - Standard amino acid residue names (no nonstandard residues)
#   - Properly protonated at your simulation pH
#   - No missing heavy atoms in regions of interest
#   - Rename to match your config.yaml enzyme.pdb_path setting
#
# PREPARATION:
#   Use `polyzymd clean-pdb -i your_protein.pdb` or another PDB
#   preparation tool to replace nonstandard residues and add
#   missing hydrogens before building your system.
#
# Delete this placeholder file after adding your protein structure.
# ============================================================================
""")

        ligand_placeholder = project_dir / "structures" / "place_ligand_here.placeholder.txt"
        ligand_placeholder.write_text("""\
# ============================================================================
# PLACEHOLDER: Place your ligand SDF file here (if using substrate)
# ============================================================================
#
# If your simulation includes a docked substrate or small molecule,
# place the SDF file in this directory.
#
# LIGAND SDF FILE (optional):
#   - 3D coordinates (from docking, crystal structure, or conformer generation)
#   - Correct protonation state for simulation pH
#   - Single conformer recommended (or specify conformer_index in config)
#   - Rename to match your config.yaml substrate.sdf_path setting
#
# SUPPORTED FORMATS:
#   - SDF (recommended)
#   - MOL2
#
# If you're not using a substrate, you can delete this placeholder
# and comment out the 'substrate' section in config.yaml.
# ============================================================================
""")

        # Success message
        colored_echo()
        click.echo(click.style("Project created successfully!", fg="green"))
        colored_echo()
        colored_echo("Directory structure:")
        colored_echo(f"  {name}/")
        colored_echo("  ├── config.yaml              <- Edit this file")
        colored_echo("  ├── structures/              <- Add your PDB/SDF files")
        colored_echo("  ├── job_scripts/")
        colored_echo("  └── slurm_logs/")
        colored_echo()
        colored_echo("Next steps:")
        colored_echo(f"  1. Add structure files to {name}/structures/")
        colored_echo(f"  2. Edit {name}/config.yaml (uncomment and customize sections)")
        colored_echo(f"  3. Validate: polyzymd validate -c {name}/config.yaml")
        colored_echo(f"  4. Build:    polyzymd build -c {name}/config.yaml -r 1")
        colored_echo()
        colored_echo(
            "Documentation: https://polyzymd.readthedocs.io/en/latest/tutorials/quickstart.html"
        )

    except Exception as e:
        # Clean up on failure
        if project_dir.exists():
            shutil.rmtree(project_dir)
        click.echo(click.style(f"Error creating project: {e}", fg="red"), err=True)
        sys.exit(1)


# =============================================================================
# Clean-PDB Command
# =============================================================================


@cli.command("clean-pdb")
@click.option(
    "-i",
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the input PDB file.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    default=None,
    type=click.Path(dir_okay=False),
    help="Path for the cleaned PDB file. Defaults to <name>_clean.pdb.",
)
@click.option(
    "--ph",
    default=7.4,
    type=float,
    show_default=True,
    help="pH for hydrogen addition.",
)
def clean_pdb(input_path: str, output_path: str | None, ph: float) -> None:
    """Clean a PDB file for use with PolyzyMD.

    Replaces nonstandard residues with their standard equivalents and adds
    missing hydrogens at the specified pH.  Chain IDs and residue numbers
    are preserved in the output.

    Requires the ``pdbfixer`` package (available via conda-forge).

    \b
    Examples:
        polyzymd clean-pdb -i structures/my_protein.pdb
        polyzymd clean-pdb -i raw.pdb -o cleaned.pdb --ph 7.0
    """
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    input_file = Path(input_path)
    if output_path is None:
        output_file = input_file.with_name(f"{input_file.stem}_clean.pdb")
    else:
        output_file = Path(output_path)

    colored_echo(f"Cleaning PDB: {input_file}")
    colored_echo(f"  pH: {ph}")

    fixer = PDBFixer(filename=str(input_file))

    fixer.findNonstandardResidues()
    n_nonstandard = len(fixer.nonstandardResidues)
    if n_nonstandard > 0:
        colored_echo(f"  Replacing {n_nonstandard} nonstandard residue(s)...")
    fixer.replaceNonstandardResidues()

    colored_echo("  Adding missing hydrogens...")
    fixer.addMissingHydrogens(ph)

    with open(output_file, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    colored_echo()
    click.echo(click.style(f"Cleaned PDB written to: {output_file}", fg="green"))


# =============================================================================
# Recover Command — resume a stalled simulation
# =============================================================================


@cli.command()
@click.option(
    "-c",
    "--config",
    required=True,
    type=click.Path(exists=True),
    help="Path to YAML configuration file",
)
@click.option(
    "-r",
    "--replicate",
    default=1,
    type=int,
    help="Replicate number (default: 1)",
)
@click.option(
    "--scratch-dir",
    default=None,
    type=click.Path(),
    help="Override scratch directory",
)
@click.option(
    "--preset",
    type=click.Choice(["aa100", "al40", "blanca-shirts", "bridges2", "testing"]),
    default="aa100",
    help="SLURM partition preset (default: aa100)",
)
@click.option(
    "--submit/--no-submit",
    default=False,
    help="Submit a self-resubmitting job to resume (default: status only)",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show status and what would be submitted without actually submitting",
)
@click.option(
    "--memory",
    default=None,
    help="Override SLURM memory allocation (e.g. '4G', '8G'). Not needed for bridges2 (allocated per GPU).",
)
@click.option(
    "--pixi-env",
    default=None,
    type=click.Choice(["cuda-12-4", "cuda-12-6"]),
    help=("Pixi environment for the recovery SLURM job. If omitted, inferred from --preset."),
)
def recover(
    config: str,
    replicate: int,
    scratch_dir: Optional[str],
    preset: str,
    submit: bool,
    dry_run: bool,
    memory: Optional[str],
    pixi_env: Optional[str],
) -> None:
    """Resume a stalled or interrupted simulation.

    Scans the working directory, loads progress state, and reports how
    much work remains. With ``--submit``, generates and submits a
    self-resubmitting SLURM job that will automatically continue from
    the last completed segment.

    \b
    Examples:
        # Check status only
        polyzymd recover -c config.yaml -r 1

        # Submit a recovery job
        polyzymd recover -c config.yaml -r 1 --submit --preset blanca-shirts

        # Dry-run (show what would be submitted)
        polyzymd recover -c config.yaml -r 1 --submit --dry-run
    """
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.simulation.progress import load_or_scan_progress, save_progress

    try:
        sim_config = SimulationConfig.from_yaml(config)
    except Exception as e:
        colored_echo(f"Failed to load config: {e}", err=True, phase="workflow", level=logging.ERROR)
        sys.exit(1)

    if scratch_dir:
        working_dir = Path(scratch_dir)
    else:
        working_dir = sim_config.get_working_directory(replicate)

    if not working_dir.exists():
        colored_echo(
            f"Working directory not found: {working_dir}",
            err=True,
            phase="workflow",
            level=logging.ERROR,
        )
        sys.exit(1)

    # Calculate total steps from config
    prod = sim_config.simulation_phases.production
    timestep_fs = prod.time_step
    total_steps = int(prod.duration * 1e6 / timestep_fs)
    total_samples = prod.samples

    # Load progress
    progress = load_or_scan_progress(
        working_dir=working_dir,
        config_path=str(Path(config).resolve()),
        total_steps=total_steps,
        total_samples=total_samples,
        timestep_fs=timestep_fs,
        replicate=replicate,
    )
    save_progress(working_dir, progress)

    # Report status
    pct = progress.fraction_complete() * 100
    remaining_ns = (progress.steps_remaining * timestep_fs) / 1e6

    colored_echo(f"Working directory: {working_dir}", phase="workflow")
    colored_echo(
        f"Progress: {progress.total_steps_completed}/{progress.total_steps_requested} steps "
        f"({pct:.1f}%)",
        phase="workflow",
    )
    colored_echo(f"Status: {progress.status.value}", phase="workflow")
    colored_echo(f"Segments: {len(progress.segments)}", phase="workflow")

    for seg in progress.segments:
        colored_echo(
            f"  segment {seg.index}: {seg.status.value} "
            f"({seg.duration_ns:.3f} ns, {seg.steps_completed} steps)",
            phase="workflow",
        )

    if progress.is_complete:
        click.echo(click.style("\nSimulation is complete — nothing to recover.", fg="green"))
        return

    colored_echo(
        f"\nRemaining: {remaining_ns:.3f} ns ({progress.steps_remaining} steps)",
        phase="workflow",
    )

    if not submit:
        colored_echo(
            "\nTo resume, run:\n"
            f"  polyzymd recover -c {config} -r {replicate} --submit --preset <preset>",
            phase="workflow",
        )
        return

    # Generate and submit a self-resubmitting SLURM job
    from polyzymd.workflow.daisy_chain import create_job_name
    from polyzymd.workflow.slurm import PRESET_DEFAULT_PIXI_ENV, SlurmConfig, SlurmScriptGenerator

    # Resolve pixi environment: explicit flag > preset default
    resolved_pixi_env = pixi_env or PRESET_DEFAULT_PIXI_ENV.get(preset, "cuda-12-4")

    colored_echo(
        f"\nGenerating recovery job (preset: {preset}, pixi env: {resolved_pixi_env})...",
        phase="workflow",
    )

    slurm_config = SlurmConfig.from_preset(preset)
    if memory:
        slurm_config.memory = memory
    generator = SlurmScriptGenerator(slurm_config, pixi_env=resolved_pixi_env)

    # Use the same descriptive job naming as `polyzymd submit`
    job_name = create_job_name(sim_config, replicate)
    logs_subdir = sim_config.output.slurm_logs_subdir
    output_file = f"{logs_subdir}/{job_name}.%j.out"

    config_path_abs = str(Path(config).resolve())
    script_content = generator.generate_job_script(
        config_path=config_path_abs,
        replicate=replicate,
        working_dir=str(working_dir),
        job_name=job_name,
        output_file=output_file,
    )

    # Write script
    script_dir = working_dir / "recovery_scripts"
    script_dir.mkdir(exist_ok=True)
    script_path = script_dir / f"recover_rep{replicate}.sh"
    script_path.write_text(script_content)
    script_path.chmod(0o755)

    colored_echo(f"Script: {script_path}", phase="workflow")

    if dry_run:
        colored_echo("\n[DRY RUN] Would submit:", phase="workflow")
        colored_echo(f"  sbatch {script_path}", phase="workflow")
        return

    # Submit
    import subprocess

    result = subprocess.run(
        ["sbatch", str(script_path)],
        capture_output=True,
        text=True,
    )

    if result.returncode == 0:
        colored_echo(f"Submitted: {result.stdout.strip()}", phase="workflow")
        colored_echo("Monitor with: squeue -u $USER", phase="workflow")
    else:
        colored_echo(
            f"Submission failed: {result.stderr.strip()}",
            err=True,
            phase="workflow",
            level=logging.ERROR,
        )
        sys.exit(1)


def _find_topology_pdb(working_dir: Path) -> Path:
    """Find a suitable topology PDB in the working directory.

    Parameters
    ----------
    working_dir : Path
        Working directory to search.

    Returns
    -------
    Path
        Path to the PDB file.

    Raises
    ------
    FileNotFoundError
        If no suitable PDB is found.
    """
    patterns = [
        "*solvated*.pdb",
        "*_solvated.pdb",
        "solvated_*.pdb",
        "equilibration/*_topology.pdb",
        "production_0/*_topology.pdb",
    ]
    for pattern in patterns:
        pdb_files = list(working_dir.glob(pattern))
        if pdb_files:
            return pdb_files[0]

    pdb_files = list(working_dir.glob("**/*.pdb"))
    if pdb_files:
        return pdb_files[0]

    raise FileNotFoundError(f"Could not find topology PDB in {working_dir}")


# =============================================================================
# Info Command
# =============================================================================


@cli.command()
def info() -> None:
    """Show PolyzyMD installation information."""
    from polyzymd import __version__

    colored_echo("PolyzyMD - Molecular Dynamics for Enzyme-Polymer Systems", phase="cli")
    colored_echo(f"Version: {__version__}", phase="cli")
    colored_echo("", phase="cli")

    # Check dependencies
    colored_echo("Dependencies:", phase="cli")

    try:
        import openmm

        colored_echo(f"  OpenMM: {openmm.__version__}", phase="cli")
    except ImportError:
        colored_echo("  OpenMM: NOT INSTALLED", phase="cli")

    try:
        from openff.toolkit import __version__ as off_version

        colored_echo(f"  OpenFF Toolkit: {off_version}", phase="cli")
    except ImportError:
        colored_echo("  OpenFF Toolkit: NOT INSTALLED", phase="cli")

    try:
        from openff.interchange import __version__ as int_version

        colored_echo(f"  OpenFF Interchange: {int_version}", phase="cli")
    except ImportError:
        colored_echo("  OpenFF Interchange: NOT INSTALLED", phase="cli")

    try:
        import pydantic

        colored_echo(f"  Pydantic: {pydantic.__version__}", phase="cli")
    except ImportError:
        colored_echo("  Pydantic: NOT INSTALLED", phase="cli")

    colored_echo("", phase="cli")
    colored_echo("Example configs: polyzymd/configs/examples/", phase="cli")


def main() -> int:
    """Main entry point."""
    try:
        cli()
        return 0
    except SystemExit as e:
        return e.code if isinstance(e.code, int) else 1
    except Exception as e:
        colored_echo(f"Error: {e}", err=True, phase="cli", level=logging.ERROR)
        return 1


if __name__ == "__main__":
    sys.exit(main())
