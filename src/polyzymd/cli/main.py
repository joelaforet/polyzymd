"""
PolyzyMD Command Line Interface.

This module provides the main CLI entry point for PolyzyMD, using Click
for argument parsing and command organization.

Usage:
    polyzymd --help
    polyzymd build --config simulation.yaml
    polyzymd run --config simulation.yaml --replicate 1
    polyzymd submit --config simulation.yaml --replicates 1-5
    polyzymd continue --working-dir path/to/sim --segment 2
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

import click

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


def _apply_restraints(sim_config, omm_topology, omm_system) -> None:
    """Validate and apply configured restraints to an OpenMM system.

    Iterates over restraints defined in sim_config, validates that each
    atom selection resolves to exactly one atom, and applies all enabled
    restraints to the system.

    Parameters
    ----------
    sim_config : SimulationConfig
        Simulation configuration containing restraint definitions.
    omm_topology : openmm.app.Topology
        OpenMM topology for resolving atom selections.
    omm_system : openmm.System
        OpenMM system to apply restraints to.
    """
    if not sim_config.restraints:
        return

    from polyzymd.core.restraints import RestraintFactory, apply_restraints

    click.echo(f"Applying {len(sim_config.restraints)} restraint(s)...")
    restraint_defs = []
    for r in sim_config.restraints:
        if not r.enabled:
            click.echo(f"  - {r.name}: DISABLED (skipping)")
            continue

        # Create restraint definition from config
        restraint_def = RestraintFactory.from_config(r.model_dump())

        # Validate the selection resolves to exactly one atom each
        try:
            indices1 = restraint_def.atom1.resolve(omm_topology)
            indices2 = restraint_def.atom2.resolve(omm_topology)

            if len(indices1) != 1:
                click.echo(
                    f"Error: Restraint '{r.name}' atom1 selection matched "
                    f"{len(indices1)} atoms (need exactly 1)",
                    err=True,
                )
                sys.exit(1)
            if len(indices2) != 1:
                click.echo(
                    f"Error: Restraint '{r.name}' atom2 selection matched "
                    f"{len(indices2)} atoms (need exactly 1)",
                    err=True,
                )
                sys.exit(1)

            click.echo(
                f"  - {r.name}: atom {indices1[0]} <-> atom {indices2[0]} "
                f"(type={r.type.value}, d={r.distance} A, "
                f"k={r.force_constant} kJ/mol/nm^2)"
            )
            restraint_defs.append(restraint_def)

        except ValueError as e:
            click.echo(f"Error: Restraint '{r.name}' invalid: {e}", err=True)
            sys.exit(1)

    # Apply all validated restraints to the system
    if restraint_defs:
        apply_restraints(restraint_defs, omm_topology, omm_system)
        click.echo(f"Successfully applied {len(restraint_defs)} restraint(s)")


@click.group()
@click.version_option(prog_name="polyzymd")
@click.option(
    "-q", "--quiet", is_flag=True, help="Suppress INFO messages, show warnings/errors only"
)
@click.option("--debug", is_flag=True, help="Enable DEBUG logging for troubleshooting")
@click.option(
    "--openff-logs",
    is_flag=True,
    help="Enable verbose OpenFF Interchange/Toolkit logs (suppressed by default)",
)
def cli(quiet: bool, debug: bool, openff_logs: bool) -> None:
    """PolyzyMD: MD simulations for enzyme-polymer systems.

    A toolkit for building, running, and analyzing molecular dynamics
    simulations of enzymes with co-polymers.
    """
    from polyzymd.analysis.core.logging_utils import setup_logging

    setup_logging(quiet=quiet, debug=debug)

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

    click.echo(f"Loading configuration from: {config}")

    try:
        sim_config = SimulationConfig.from_yaml(config)
        click.echo(f"Configuration validated: {sim_config.name}")

        # Override directories if provided via CLI
        effective_scratch = scratch_dir or output_dir  # output_dir is alias for scratch_dir
        if effective_scratch:
            sim_config.output.scratch_directory = Path(effective_scratch)
        if projects_dir:
            sim_config.output.projects_directory = Path(projects_dir)

        if dry_run:
            click.echo("Dry run - configuration is valid")
            click.echo(f"  Enzyme: {sim_config.enzyme.name}")
            if sim_config.substrate:
                click.echo(f"  Substrate: {sim_config.substrate.name}")
            if sim_config.polymers and sim_config.polymers.enabled:
                click.echo(f"  Polymers: {sim_config.polymers.type_prefix}")
                click.echo(f"  Polymer count: {sim_config.polymers.count}")
            click.echo(f"  Temperature: {sim_config.thermodynamics.temperature} K")
            click.echo(f"  Production time: {sim_config.simulation_phases.production.duration} ns")
            click.echo()
            click.echo("Directories:")
            click.echo(f"  Projects: {sim_config.output.projects_directory}")
            click.echo(f"  Scratch: {sim_config.output.effective_scratch_directory}")
            if gromacs:
                click.echo()
                click.echo("GROMACS export enabled:")
                click.echo("  Output: {projects_dir}/{replicate}/gromacs/")
            return

        click.echo(f"Building system for replicate {replicate}...")
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
            click.echo("Exporting to GROMACS format...")
            gromacs_dir = (
                sim_config.output.projects_directory / f"replicate_{replicate}" / "gromacs"
            )
            export_result = builder.export_to_gromacs(gromacs_dir)

            click.echo("GROMACS export successful!")
            click.echo(f"Output directory: {gromacs_dir}")
            click.echo("Files generated:")
            click.echo(f"  - {export_result['gro'].name} (coordinates)")
            click.echo(f"  - {export_result['top'].name} (topology)")
            click.echo(f"  - {export_result['em_mdp'].name} (energy minimization)")
            for eq_mdp in export_result.get("eq_mdps", []):
                click.echo(f"  - {eq_mdp.name} (equilibration)")
            click.echo(f"  - {export_result['prod_mdp'].name} (production)")
            if export_result.get("posres_defines"):
                click.echo("Position restraints added to molecule ITP files:")
                for component, define in export_result["posres_defines"].items():
                    click.echo(f"  - {component}: #ifdef {define}")
            click.echo(f"  - {export_result['run_script'].name} (run script)")
            click.echo()
            click.echo(f"To run: cd {gromacs_dir} && ./{export_result['run_script'].name}")

        else:
            # Default: prepare for OpenMM simulation
            click.echo("Extracting OpenMM components...")
            omm_topology, omm_system, omm_positions = builder.get_openmm_components()

            # Apply restraints if configured
            _apply_restraints(sim_config, omm_topology, omm_system)

            # Save OpenMM system to XML for --skip-build support
            from openmm import XmlSerializer

            system_xml_path = working_dir / "system.xml"
            click.echo(f"Saving OpenMM system to {system_xml_path}...")
            with open(system_xml_path, "w") as f:
                f.write(XmlSerializer.serialize(omm_system))

            click.echo("System built successfully!")
            click.echo(f"Output directory: {working_dir}")
            click.echo("Files saved:")
            click.echo("  - solvated_system.pdb (topology + positions)")
            click.echo("  - system.xml (OpenMM system with restraints)")
            click.echo("Use 'polyzymd run --skip-build' to run without rebuilding.")

    except FileNotFoundError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Build failed: {e}", err=True)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


# =============================================================================
# Run Command
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
    help="Scratch directory for simulation output (trajectories, checkpoints)",
)
@click.option(
    "--projects-dir",
    default=None,
    type=click.Path(),
    help="Projects directory for scripts and logs (long-term storage)",
)
@click.option(
    "--segment-time",
    default=None,
    type=float,
    help="Override production time per segment (ns) [OpenMM only]",
)
@click.option(
    "--segment-frames",
    default=None,
    type=int,
    help="Override frames per segment [OpenMM only]",
)
@click.option(
    "--skip-build",
    is_flag=True,
    help="Skip system building (use existing) [OpenMM only]",
)
@click.option(
    "--gromacs",
    is_flag=True,
    help="Run simulation using GROMACS instead of OpenMM",
)
@click.option(
    "--gmx-path",
    default="gmx",
    help="Path to GROMACS executable (default: gmx) [GROMACS only]",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Export files but don't run simulation [GROMACS only]",
)
def run(
    config: str,
    replicate: int,
    scratch_dir: Optional[str],
    projects_dir: Optional[str],
    segment_time: Optional[float],
    segment_frames: Optional[int],
    skip_build: bool,
    gromacs: bool,
    gmx_path: str,
    dry_run: bool,
) -> None:
    """Run a simulation from configuration.

    By default, runs using OpenMM. Use --gromacs to run using GROMACS instead.

    OpenMM Mode (default):
        Builds the system (unless --skip-build), runs equilibration,
        and then runs production simulation using OpenMM.

    GROMACS Mode (--gromacs):
        Builds the system, exports to GROMACS format (.gro, .top, .mdp),
        and executes the full GROMACS workflow locally (EM, equilibration,
        production, and trajectory post-processing).

    GROMACS Notes:
        - Requires GROMACS to be installed and accessible
        - Use --gmx-path to specify a custom GROMACS executable
        - Use --dry-run to export files without running the simulation
        - --skip-build is not supported for GROMACS (always rebuilds)
    """
    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.config.schema import SimulationConfig

    click.echo(f"Loading configuration from: {config}")

    try:
        sim_config = SimulationConfig.from_yaml(config)
        click.echo(f"Running simulation: {sim_config.name}")

        # Override directories if provided via CLI
        if scratch_dir:
            sim_config.output.scratch_directory = Path(scratch_dir)
        if projects_dir:
            sim_config.output.projects_directory = Path(projects_dir)

        # Branch based on simulation engine
        if gromacs:
            _run_gromacs(
                sim_config=sim_config,
                replicate=replicate,
                gmx_path=gmx_path,
                dry_run=dry_run,
                skip_build=skip_build,
            )
        else:
            _run_openmm(
                sim_config=sim_config,
                replicate=replicate,
                scratch_dir=scratch_dir,
                segment_time=segment_time,
                segment_frames=segment_frames,
                skip_build=skip_build,
            )

    except Exception as e:
        click.echo(f"Simulation failed: {e}", err=True)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


def _run_gromacs(
    sim_config: "SimulationConfig",
    replicate: int,
    gmx_path: str,
    dry_run: bool,
    skip_build: bool,
) -> None:
    """Run simulation using GROMACS.

    Args:
        sim_config: Validated simulation configuration.
        replicate: Replicate number.
        gmx_path: Path to GROMACS executable.
        dry_run: If True, export files but don't run simulation.
        skip_build: If True, skip system building (not supported for GROMACS MVP).
    """
    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.exporters.gromacs import GromacsError, GromacsExporter, GromacsRunner

    # Warn about unsupported options
    if skip_build:
        click.echo(
            "Warning: --skip-build is not supported for GROMACS mode. System will be rebuilt.",
            err=True,
        )

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
    click.echo("Exporting to GROMACS format...")
    exporter = GromacsExporter(
        interchange=interchange,
        config=sim_config,
        component_info=component_info,
    )
    export_result = exporter.export(
        output_dir=gromacs_dir,
        gmx_command=gmx_path,
    )

    click.echo(f"\nGROMACS files exported to: {gromacs_dir}")
    click.echo("Files generated:")
    click.echo(f"  - {export_result['gro'].name} (coordinates)")
    click.echo(f"  - {export_result['top'].name} (topology)")
    click.echo(f"  - {export_result['em_mdp'].name} (energy minimization)")
    for eq_mdp in export_result["eq_mdps"]:
        click.echo(f"  - {eq_mdp.name} (equilibration)")
    click.echo(f"  - {export_result['prod_mdp'].name} (production)")
    if export_result.get("posres_defines"):
        click.echo("Position restraints added to molecule ITP files:")
        for component, define in export_result["posres_defines"].items():
            click.echo(f"  - {component}: #ifdef {define}")
    click.echo(f"  - {export_result['run_script'].name} (run script)")

    if dry_run:
        click.echo("\n--dry-run specified: Files exported but simulation not started.")
        click.echo(f"To run manually: cd {gromacs_dir} && ./{export_result['run_script'].name}")
        return

    # Run GROMACS workflow
    click.echo("\nStarting GROMACS simulation...")
    click.echo(f"Using GROMACS executable: {gmx_path}")

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

        click.echo("\nGROMACS simulation completed successfully!")
        click.echo(f"Output directory: {gromacs_dir}")

    except GromacsError as e:
        click.echo(f"\nGROMACS simulation failed: {e}", err=True)
        click.echo(f"Check log files in: {gromacs_dir}", err=True)
        sys.exit(1)

    except FileNotFoundError as e:
        click.echo(f"\nError: {e}", err=True)
        click.echo("Ensure GROMACS is installed and in your PATH, or use --gmx-path.", err=True)
        sys.exit(1)


def _run_openmm(
    sim_config: "SimulationConfig",
    replicate: int,
    scratch_dir: Optional[str],
    segment_time: Optional[float],
    segment_frames: Optional[int],
    skip_build: bool,
) -> None:
    """Run simulation using OpenMM.

    Args:
        sim_config: Validated simulation configuration.
        replicate: Replicate number.
        scratch_dir: Override for scratch directory.
        segment_time: Override for production time per segment.
        segment_frames: Override for frames per segment.
        skip_build: If True, load pre-built system from disk.
    """
    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.simulation.runner import SimulationRunner

    # Determine working directory
    if scratch_dir:
        working_dir = Path(scratch_dir)
    else:
        working_dir = sim_config.get_working_directory(replicate)

    if not skip_build:
        click.echo(f"Building system for replicate {replicate}...")
        builder = SystemBuilder.from_config(sim_config)
        interchange = builder.build_from_config(
            config=sim_config,
            working_dir=working_dir,
            polymer_seed=replicate,
        )

        # Extract OpenMM components from Interchange
        click.echo("Extracting OpenMM components...")
        omm_topology, omm_system, omm_positions = builder.get_openmm_components()

        # Apply restraints if configured
        _apply_restraints(sim_config, omm_topology, omm_system)

    else:
        # --skip-build: Load pre-built system from disk
        click.echo("Skipping build, loading pre-built system...")
        from openmm import XmlSerializer
        from openmm.app import PDBFile

        # Check that required files exist
        pdb_path = working_dir / "solvated_system.pdb"
        system_path = working_dir / "system.xml"

        if not pdb_path.exists():
            click.echo(f"Error: {pdb_path} not found. Run 'polyzymd build' first.", err=True)
            sys.exit(1)
        if not system_path.exists():
            click.echo(f"Error: {system_path} not found. Run 'polyzymd build' first.", err=True)
            sys.exit(1)

        # Load topology and positions from PDB
        click.echo(f"Loading topology and positions from {pdb_path}...")
        pdb = PDBFile(str(pdb_path))
        omm_topology = pdb.topology
        omm_positions = pdb.positions

        # Load system from XML (already includes restraints from build)
        click.echo(f"Loading OpenMM system from {system_path}...")
        with open(system_path, "r") as f:
            omm_system = XmlSerializer.deserialize(f.read())

        click.echo("Pre-built system loaded successfully")

    # Create runner
    runner = SimulationRunner(
        topology=omm_topology,
        system=omm_system,
        positions=omm_positions,
        working_dir=working_dir,
    )

    # Run energy minimization first
    click.echo("Running energy minimization...")
    runner.minimize()

    # Get thermodynamic parameters
    temperature = sim_config.thermodynamics.temperature
    pressure = sim_config.thermodynamics.pressure

    # Run equilibration
    phases = sim_config.simulation_phases
    eq_duration = phases.total_equilibration_duration
    eq_mode = "multi-stage" if phases.uses_staged_equilibration else "simple"
    click.echo(f"Running equilibration: {eq_duration:.3f} ns at {temperature} K ({eq_mode})...")
    runner.run_equilibration(
        temperature=temperature,
        config=phases,
    )

    # Calculate segment parameters
    total_time = sim_config.simulation_phases.production.duration
    num_segments = sim_config.simulation_phases.segments
    seg_time = segment_time or (total_time / num_segments)
    seg_frames = segment_frames or (sim_config.simulation_phases.production.samples // num_segments)

    # Run first production segment
    click.echo(f"Running production segment 0: {seg_time} ns, {seg_frames} frames (NPT)...")
    runner.run_production(
        temperature=temperature,
        duration_ns=seg_time,
        num_samples=seg_frames,
        pressure=pressure,
        segment_index=0,
    )

    click.echo("Simulation completed successfully!")
    click.echo(f"Output: {working_dir}")


# =============================================================================
# Submit Command (Daisy-chain)
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
    type=click.Choice(["aa100", "al40", "blanca-shirts", "testing"]),
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
    help="Override SLURM memory allocation (default: 3G). Increase to 4-8G for larger systems or if you encounter OOM errors.",
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
    submit_openff_logs: bool,
    skip_build: bool,
) -> None:
    """Submit daisy-chain simulation jobs to SLURM.

    Creates and submits a chain of dependent jobs that will run
    sequentially, allowing long simulations to be broken into
    smaller segments that fit within HPC time limits.

    Directory structure:
    - projects_dir: Where job scripts and SLURM logs are stored (long-term storage)
    - scratch_dir: Where simulation data is written (high-performance storage)
    """
    from polyzymd.workflow.daisy_chain import submit_daisy_chain

    click.echo(f"Loading configuration from: {config}")
    click.echo(f"Submitting jobs with preset: {preset}")
    click.echo(f"Replicates: {replicates}")
    if scratch_dir:
        click.echo(f"Scratch directory: {scratch_dir}")
    if projects_dir:
        click.echo(f"Projects directory: {projects_dir}")
    if memory:
        click.echo(f"Memory allocation: {memory}")
    if skip_build:
        click.echo("Skip-build mode: using pre-built systems")

    if dry_run:
        click.echo("DRY RUN MODE - scripts will be created but not submitted")

    try:
        results = submit_daisy_chain(
            config_path=config,
            slurm_preset=preset,
            replicates=replicates,
            email=email,
            dry_run=dry_run,
            output_dir=output_dir,
            scratch_dir=scratch_dir,
            projects_dir=projects_dir,
            time_limit=time_limit,
            memory=memory,
            openff_logs=submit_openff_logs,
            skip_build=skip_build,
        )

        if not dry_run:
            click.echo("\nJob submission complete!")
            click.echo("Monitor with: squeue -u $USER")

    except Exception as e:
        click.echo(f"Submission failed: {e}", err=True)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


# =============================================================================
# Continue Command
# =============================================================================


@cli.command("continue")
@click.option(
    "-w",
    "--working-dir",
    required=True,
    type=click.Path(exists=True),
    help="Working directory with previous segment",
)
@click.option(
    "-s",
    "--segment",
    required=True,
    type=int,
    help="Segment index to run (1-based)",
)
@click.option(
    "-t",
    "--segment-time",
    required=True,
    type=float,
    help="Duration of this segment (ns)",
)
@click.option(
    "-n",
    "--num-samples",
    default=250,
    type=int,
    help="Number of frames to save (default: 250)",
)
def continue_sim(
    working_dir: str,
    segment: int,
    segment_time: float,
    num_samples: int,
) -> None:
    """Continue a simulation from a previous segment.

    Loads state from the previous production segment and continues
    the simulation. Used by daisy-chain continuation jobs.
    """
    from polyzymd.simulation.continuation import ContinuationManager

    click.echo(f"Continuing simulation from segment {segment - 1}")
    click.echo(f"Working directory: {working_dir}")
    click.echo(f"Duration: {segment_time} ns, Frames: {num_samples}")

    try:
        manager = ContinuationManager(
            working_dir=working_dir,
            segment_index=segment,
        )

        click.echo("Loading previous state...")
        manager.load_previous_state()

        click.echo(f"Running segment {segment}...")
        results = manager.run_segment(
            duration_ns=segment_time,
            num_samples=num_samples,
        )

        click.echo(f"Segment {segment} completed!")
        click.echo(f"Output: {results['output_dir']}")

    except FileNotFoundError as e:
        click.echo(f"Error: Previous segment not found: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Continuation failed: {e}", err=True)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


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

    click.echo(f"Validating configuration: {config}")

    try:
        sim_config = SimulationConfig.from_yaml(config)

        click.echo(click.style("Configuration is valid!", fg="green"))
        click.echo()
        click.echo("Summary:")
        click.echo(f"  Name: {sim_config.name}")
        click.echo(f"  Enzyme: {sim_config.enzyme.name}")

        if sim_config.substrate:
            click.echo(f"  Substrate: {sim_config.substrate.name}")
        else:
            click.echo("  Substrate: None (apo simulation)")

        if sim_config.polymers and sim_config.polymers.enabled:
            click.echo(f"  Polymers: {sim_config.polymers.type_prefix}")
            click.echo(f"    Count: {sim_config.polymers.count}")
            click.echo(f"    Length: {sim_config.polymers.length}")
            for m in sim_config.polymers.monomers:
                click.echo(f"    Monomer {m.label}: {m.probability * 100:.1f}%")
        else:
            click.echo("  Polymers: Disabled")

        click.echo(f"  Temperature: {sim_config.thermodynamics.temperature} K")
        click.echo(f"  Pressure: {sim_config.thermodynamics.pressure} atm")
        click.echo()
        click.echo("Simulation phases:")
        eq = sim_config.simulation_phases.equilibration
        click.echo(f"  Equilibration: {eq.duration} ns ({eq.ensemble.value})")
        prod = sim_config.simulation_phases.production
        click.echo(f"  Production: {prod.duration} ns ({prod.ensemble.value})")
        click.echo(f"  Segments: {sim_config.simulation_phases.segments}")

        if sim_config.restraints:
            click.echo()
            click.echo(f"Restraints: {len(sim_config.restraints)}")
            for r in sim_config.restraints:
                status = "enabled" if r.enabled else "disabled"
                click.echo(f"  - {r.name} ({r.type.value}): {status}")

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
        click.echo("Choose a different name or remove the existing directory.")
        sys.exit(1)

    try:
        # Create directory structure
        click.echo(f"Creating project directory: {name}/")
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
#   - Properly protonated (use PDB2PQR, Reduce, or similar)
#   - No missing residues in regions of interest
#   - Standard amino acid residue names
#   - Rename to match your config.yaml enzyme.pdb_path setting
#
# PREPARATION TIPS:
#   1. Download structure from PDB or use your own model
#   2. Remove waters, ligands, and alternate conformations
#   3. Add hydrogens at your simulation pH
#   4. Check for missing loops/residues
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
        click.echo()
        click.echo(click.style("Project created successfully!", fg="green"))
        click.echo()
        click.echo("Directory structure:")
        click.echo(f"  {name}/")
        click.echo("  ├── config.yaml              <- Edit this file")
        click.echo("  ├── structures/              <- Add your PDB/SDF files")
        click.echo("  ├── job_scripts/")
        click.echo("  └── slurm_logs/")
        click.echo()
        click.echo("Next steps:")
        click.echo(f"  1. Add structure files to {name}/structures/")
        click.echo(f"  2. Edit {name}/config.yaml (uncomment and customize sections)")
        click.echo(f"  3. Validate: polyzymd validate -c {name}/config.yaml")
        click.echo(f"  4. Build:    polyzymd build -c {name}/config.yaml -r 1")
        click.echo()
        click.echo(
            "Documentation: https://polyzymd.readthedocs.io/en/latest/tutorials/quickstart.html"
        )

    except Exception as e:
        # Clean up on failure
        if project_dir.exists():
            shutil.rmtree(project_dir)
        click.echo(click.style(f"Error creating project: {e}", fg="red"), err=True)
        sys.exit(1)


# =============================================================================
# Info Command
# =============================================================================


@cli.command()
def info() -> None:
    """Show PolyzyMD installation information."""
    from polyzymd import __version__

    click.echo("PolyzyMD - Molecular Dynamics for Enzyme-Polymer Systems")
    click.echo(f"Version: {__version__}")
    click.echo()

    # Check dependencies
    click.echo("Dependencies:")

    try:
        import openmm

        click.echo(f"  OpenMM: {openmm.__version__}")
    except ImportError:
        click.echo("  OpenMM: NOT INSTALLED")

    try:
        from openff.toolkit import __version__ as off_version

        click.echo(f"  OpenFF Toolkit: {off_version}")
    except ImportError:
        click.echo("  OpenFF Toolkit: NOT INSTALLED")

    try:
        from openff.interchange import __version__ as int_version

        click.echo(f"  OpenFF Interchange: {int_version}")
    except ImportError:
        click.echo("  OpenFF Interchange: NOT INSTALLED")

    try:
        import pydantic

        click.echo(f"  Pydantic: {pydantic.__version__}")
    except ImportError:
        click.echo("  Pydantic: NOT INSTALLED")

    click.echo()
    click.echo("Example configs: polyzymd/configs/examples/")


# =============================================================================
# Analysis Commands (from analysis module)
# =============================================================================

# Register analysis command groups
from polyzymd.analysis.cli import analyze, plot

cli.add_command(analyze)
cli.add_command(plot)

# Register compare command group
from polyzymd.compare.cli import compare

cli.add_command(compare)


def main() -> int:
    """Main entry point."""
    try:
        cli()
        return 0
    except SystemExit as e:
        return e.code if isinstance(e.code, int) else 1
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
