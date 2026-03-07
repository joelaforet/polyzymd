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

# Setup logging
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
def cli(verbose: bool, openff_logs: bool) -> None:
    """PolyzyMD: MD simulations for enzyme-polymer systems.

    A toolkit for building, running, and analyzing molecular dynamics
    simulations of enzymes with co-polymers.
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

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
            if sim_config.restraints:
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
                            f"(type={r.type.value}, d={r.distance} A, k={r.force_constant} kJ/mol/nm^2)"
                        )
                        restraint_defs.append(restraint_def)

                    except ValueError as e:
                        click.echo(f"Error: Restraint '{r.name}' invalid: {e}", err=True)
                        sys.exit(1)

                # Apply all validated restraints to the system
                if restraint_defs:
                    apply_restraints(restraint_defs, omm_topology, omm_system)
                    click.echo(f"Successfully applied {len(restraint_defs)} restraint(s)")

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
        if sim_config.restraints:
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
                        f"(type={r.type.value}, d={r.distance} A, k={r.force_constant} kJ/mol/nm^2)"
                    )
                    restraint_defs.append(restraint_def)

                except ValueError as e:
                    click.echo(f"Error: Restraint '{r.name}' invalid: {e}", err=True)
                    sys.exit(1)

            # Apply all validated restraints to the system
            if restraint_defs:
                apply_restraints(restraint_defs, omm_topology, omm_system)
                click.echo(f"Successfully applied {len(restraint_defs)} restraint(s)")

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
    if account:
        click.echo(f"Account: {account}")
    if memory:
        click.echo(f"Memory allocation: {memory}")
    if gpu_type:
        click.echo(f"GPU type override: {gpu_type}")
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
            account=account,
            gpu_type=gpu_type,
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
        # Check if this is a GracefulExit (signal-based interruption)
        from polyzymd.simulation.signals import EXIT_CODE_INTERRUPTED, GracefulExit

        if isinstance(e, GracefulExit):
            click.echo(f"Segment {segment} interrupted (graceful shutdown): {e}", err=True)
            sys.exit(EXIT_CODE_INTERRUPTED)
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

    click.echo(f"Cleaning PDB: {input_file}")
    click.echo(f"  pH: {ph}")

    fixer = PDBFixer(filename=str(input_file))

    fixer.findNonstandardResidues()
    n_nonstandard = len(fixer.nonstandardResidues)
    if n_nonstandard > 0:
        click.echo(f"  Replacing {n_nonstandard} nonstandard residue(s)...")
    fixer.replaceNonstandardResidues()

    click.echo("  Adding missing hydrogens...")
    fixer.addMissingHydrogens(ph)

    with open(output_file, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    click.echo()
    click.echo(click.style(f"Cleaned PDB written to: {output_file}", fg="green"))


# =============================================================================
# Recover Command — finish an interrupted segment
# =============================================================================


@cli.command()
@click.option(
    "-w",
    "--working-dir",
    required=True,
    type=click.Path(exists=True),
    help="Working directory containing simulation output",
)
@click.option(
    "-s",
    "--segment",
    required=True,
    type=int,
    help="Segment index that was interrupted",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show what would be recovered without running",
)
def recover(working_dir: str, segment: int, dry_run: bool) -> None:
    """Recover an interrupted simulation segment.

    Reads the INTERRUPTED marker and emergency state files from a
    previously interrupted segment, then runs the remaining steps
    in a ``production_N_resume/`` subdirectory.

    After recovery, the normal end-of-segment files (state.xml,
    system.xml, parameters.json) are written to the original
    ``production_N/`` directory so that the next ContinuationManager
    can pick up seamlessly.
    """
    import json

    from openmm import XmlSerializer
    from openmm.app import PDBFile, Simulation

    from polyzymd.simulation.continuation import ContinuationManager, quantity_from_dict

    work = Path(working_dir)
    seg_dir = work / f"production_{segment}"
    marker = seg_dir / "INTERRUPTED"

    if not marker.exists():
        click.echo(f"No INTERRUPTED marker found in {seg_dir}", err=True)
        sys.exit(1)

    # Parse INTERRUPTED marker
    meta = {}
    for line in marker.read_text().strip().splitlines():
        key, val = line.split("=", 1)
        meta[key.strip()] = int(val.strip())

    steps_completed = meta["steps_completed"]
    total_steps = meta["total_steps"]
    remaining_steps = meta["remaining_steps"]

    click.echo(f"Segment {segment}: {steps_completed}/{total_steps} steps completed")
    click.echo(f"  Remaining: {remaining_steps} steps")

    # Load parameters to get reporter settings
    # Try segment's own parameters first, then look in the segment directory
    param_path = seg_dir / f"production_{segment}_parameters.json"
    if not param_path.exists():
        click.echo(f"Parameters file not found: {param_path}", err=True)
        sys.exit(1)

    with open(param_path) as f:
        param_dict = json.load(f)

    integ_raw = param_dict["__values__"]["integ_params"]["__values__"]
    num_samples = integ_raw.get("num_samples", 250)
    report_interval = max(1, total_steps // num_samples)

    # Calculate remaining frames
    remaining_frames = remaining_steps // report_interval
    click.echo(f"  Remaining frames: ~{remaining_frames}")

    if dry_run:
        click.echo()
        click.echo("[DRY RUN] Would recover segment, no simulation run.")
        return

    # Verify emergency state files exist
    emergency_state = seg_dir / "emergency_state.xml"
    emergency_system = seg_dir / "emergency_system.xml"
    if not emergency_state.exists():
        click.echo(f"Emergency state not found: {emergency_state}", err=True)
        sys.exit(1)
    if not emergency_system.exists():
        click.echo(f"Emergency system not found: {emergency_system}", err=True)
        sys.exit(1)

    # Create resume subdirectory
    resume_dir = work / f"production_{segment}_resume"
    resume_dir.mkdir(exist_ok=True)
    click.echo(f"  Resume directory: {resume_dir}")

    # Load system and state
    click.echo("Loading emergency state...")
    import openmm

    with open(emergency_system) as f:
        system = XmlSerializer.deserialize(f.read())

    # Load topology (same as ContinuationManager)
    pdb_path = _find_topology_pdb(work)
    topology = PDBFile(str(pdb_path)).topology

    # Create integrator from parameters
    thermo_raw = param_dict["__values__"]["thermo_params"]["__values__"]
    thermostat_raw = thermo_raw["thermostat_params"]["__values__"]
    time_step = quantity_from_dict(integ_raw["time_step"])
    temperature = quantity_from_dict(thermostat_raw["temperature"])
    friction_coeff = quantity_from_dict(thermostat_raw["timescale"])

    integrator = openmm.LangevinMiddleIntegrator(temperature, friction_coeff, time_step)

    # Create simulation and load state
    simulation = Simulation(topology, system, integrator)
    simulation.loadState(str(emergency_state))

    # Setup reporters for resume directory
    from openmm.app import DCDReporter, StateDataReporter

    resume_report_interval = max(1, remaining_steps // max(1, remaining_frames))

    traj_path = resume_dir / f"production_{segment}_resume_trajectory.dcd"
    state_data_path = resume_dir / f"production_{segment}_resume_state_data.csv"

    simulation.reporters.append(DCDReporter(str(traj_path), resume_report_interval))
    simulation.reporters.append(
        StateDataReporter(
            str(state_data_path),
            resume_report_interval,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
        )
    )

    # Install signal handlers (recovery can itself be interrupted)
    from polyzymd.simulation.signals import (
        GracefulExit,
        install_handlers,
        is_interrupted,
        save_emergency_state,
    )

    install_handlers()

    # Run remaining steps (chunked for signal checking)
    click.echo(f"Running {remaining_steps} remaining steps...")
    chunk_size = min(resume_report_interval, remaining_steps)
    steps_done = 0
    try:
        while steps_done < remaining_steps:
            this_chunk = min(chunk_size, remaining_steps - steps_done)
            simulation.step(this_chunk)
            steps_done += this_chunk
            if is_interrupted():
                click.echo(
                    f"Recovery interrupted at step {steps_done}/{remaining_steps}",
                    err=True,
                )
                save_emergency_state(
                    simulation=simulation,
                    output_dir=seg_dir,
                    segment_index=segment,
                    steps_completed=steps_completed + steps_done,
                    total_steps=total_steps,
                )
                sys.exit(99)
    except GracefulExit:
        sys.exit(99)

    click.echo("Recovery simulation complete — writing end-of-segment files")

    # Write normal end-of-segment files to the ORIGINAL segment directory
    # so that ContinuationManager for segment+1 can find them
    state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getForces=True,
        getEnergy=True,
        getParameters=True,
    )

    state_xml_path = seg_dir / f"production_{segment}_state.xml"
    with open(state_xml_path, "w") as f:
        f.write(XmlSerializer.serialize(state))

    system_xml_path = seg_dir / f"production_{segment}_system.xml"
    with open(system_xml_path, "w") as f:
        f.write(XmlSerializer.serialize(system))

    # Copy parameters file (already exists, but confirm it's there)
    # param_path already verified above

    # Save checkpoint
    chk_path = seg_dir / f"production_{segment}_checkpoint.chk"
    simulation.saveCheckpoint(str(chk_path))

    # Remove INTERRUPTED marker — segment is now complete
    marker.unlink()
    click.echo(f"Removed INTERRUPTED marker from {seg_dir}")

    click.echo(click.style(f"Segment {segment} recovery complete!", fg="green"))


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
# Recover-chain Command — scan and resubmit a broken chain
# =============================================================================


@cli.command("recover-chain")
@click.option(
    "-w",
    "--working-dir",
    required=True,
    type=click.Path(exists=True),
    help="Working directory containing simulation output",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show chain status without taking action",
)
def recover_chain(working_dir: str, dry_run: bool) -> None:
    """Scan a daisy-chain working directory and report segment status.

    For each ``production_N/`` directory, reports whether the segment
    completed, was interrupted (recoverable), or is missing (crashed).

    With ``--dry-run`` (default behavior), only prints the status report.
    Without ``--dry-run``, recovers any interrupted segments in-place.
    """
    work = Path(working_dir)

    # Find all production_N directories
    seg_dirs = sorted(work.glob("production_*"))
    seg_dirs = [d for d in seg_dirs if d.is_dir() and "_resume" not in d.name]

    if not seg_dirs:
        click.echo("No production segments found.")
        return

    click.echo(f"Scanning {work}")
    click.echo(f"Found {len(seg_dirs)} segment(s):\n")

    statuses = []
    for seg_dir in seg_dirs:
        seg_name = seg_dir.name
        # Extract segment index from directory name
        try:
            seg_idx = int(seg_name.split("_")[1])
        except (IndexError, ValueError):
            continue

        marker = seg_dir / "INTERRUPTED"
        state_xml = seg_dir / f"production_{seg_idx}_state.xml"

        if marker.exists():
            # Parse marker for details
            meta = {}
            for line in marker.read_text().strip().splitlines():
                k, v = line.split("=", 1)
                meta[k.strip()] = int(v.strip())
            pct = 100 * meta.get("steps_completed", 0) / max(meta.get("total_steps", 1), 1)
            status = f"INTERRUPTED ({pct:.0f}% done, {meta.get('remaining_steps', '?')} remaining)"
            statuses.append(("interrupted", seg_idx))
        elif state_xml.exists():
            status = "COMPLETED"
            statuses.append(("completed", seg_idx))
        else:
            status = "MISSING (no state.xml, no INTERRUPTED marker — likely crashed)"
            statuses.append(("missing", seg_idx))

        # Check for resume directory
        resume_dir = work / f"production_{seg_idx}_resume"
        resume_note = " [has resume/]" if resume_dir.exists() else ""

        click.echo(f"  segment {seg_idx:>3d}: {status}{resume_note}")

    click.echo()

    interrupted = [idx for st, idx in statuses if st == "interrupted"]
    missing = [idx for st, idx in statuses if st == "missing"]
    completed = [idx for st, idx in statuses if st == "completed"]

    click.echo(
        f"Summary: {len(completed)} completed, "
        f"{len(interrupted)} interrupted, {len(missing)} missing"
    )

    if not interrupted and not missing:
        click.echo(click.style("All segments healthy!", fg="green"))
        return

    if missing:
        click.echo(
            click.style(
                f"\nSegments {missing} have no recoverable state — "
                "these must be re-run from the last completed segment.",
                fg="red",
            )
        )

    if interrupted and not dry_run:
        click.echo(f"\nRecovering {len(interrupted)} interrupted segment(s)...")
        from click.testing import CliRunner

        runner = CliRunner()
        for seg_idx in interrupted:
            click.echo(f"\n--- Recovering segment {seg_idx} ---")
            result = runner.invoke(
                recover,
                ["-w", working_dir, "-s", str(seg_idx)],
            )
            click.echo(result.output)
            if result.exit_code != 0:
                click.echo(
                    click.style(f"Recovery of segment {seg_idx} failed!", fg="red"),
                    err=True,
                )
                sys.exit(1)

        click.echo(click.style("\nAll interrupted segments recovered!", fg="green"))
    elif interrupted and dry_run:
        click.echo(f"\n[DRY RUN] Would recover segments: {interrupted}")


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
