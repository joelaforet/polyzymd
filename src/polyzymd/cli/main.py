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


@click.group()
@click.version_option(prog_name="polyzymd")
@click.option("-v", "--verbose", is_flag=True, help="Enable verbose output")
def cli(verbose: bool) -> None:
    """PolyzyMD: MD simulations for enzyme-polymer systems.

    A toolkit for building, running, and analyzing molecular dynamics
    simulations of enzymes with co-polymers.
    """
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)


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
def build(
    config: str,
    replicate: int,
    output_dir: Optional[str],
    scratch_dir: Optional[str],
    projects_dir: Optional[str],
    dry_run: bool,
) -> None:
    """Build a simulation system from configuration.

    Loads the YAML configuration, constructs the molecular system
    (enzyme, substrate, polymers, solvent), and prepares it for simulation.
    """
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.builders.system_builder import SystemBuilder

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
            return

        click.echo(f"Building system for replicate {replicate}...")
        working_dir = sim_config.get_working_directory(replicate)
        builder = SystemBuilder.from_config(sim_config)
        interchange = builder.build_from_config(
            config=sim_config,
            working_dir=working_dir,
            polymer_seed=replicate,
        )

        click.echo(f"System built successfully!")
        click.echo(f"Output directory: {working_dir}")

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
    help="Override production time per segment (ns)",
)
@click.option(
    "--segment-frames",
    default=None,
    type=int,
    help="Override frames per segment",
)
@click.option(
    "--skip-build",
    is_flag=True,
    help="Skip system building (use existing)",
)
def run(
    config: str,
    replicate: int,
    scratch_dir: Optional[str],
    projects_dir: Optional[str],
    segment_time: Optional[float],
    segment_frames: Optional[int],
    skip_build: bool,
) -> None:
    """Run a simulation from configuration.

    Builds the system (unless --skip-build), runs equilibration,
    and then runs production simulation.
    """
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.simulation.runner import SimulationRunner

    click.echo(f"Loading configuration from: {config}")

    try:
        sim_config = SimulationConfig.from_yaml(config)
        click.echo(f"Running simulation: {sim_config.name}")

        # Override directories if provided via CLI
        if scratch_dir:
            sim_config.output.scratch_directory = Path(scratch_dir)
        if projects_dir:
            sim_config.output.projects_directory = Path(projects_dir)

        working_dir = sim_config.get_working_directory(replicate)

        if not skip_build:
            click.echo(f"Building system for replicate {replicate}...")
            builder = SystemBuilder.from_config(sim_config)
            interchange = builder.build_from_config(
                config=sim_config,
                working_dir=working_dir,
                polymer_seed=replicate,
            )
        else:
            click.echo("Skipping build, loading existing system...")
            # Load existing interchange would go here
            raise NotImplementedError("--skip-build requires pre-built system")

        # Extract OpenMM components from Interchange
        click.echo("Extracting OpenMM components...")
        omm_topology, omm_system, omm_positions = builder.get_openmm_components()

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
        eq_config = sim_config.simulation_phases.equilibration
        click.echo(f"Running equilibration: {eq_config.duration} ns at {temperature} K (NVT)...")
        runner.run_equilibration(
            temperature=temperature,
            duration_ns=eq_config.duration,
            num_samples=eq_config.samples,
            timestep_fs=eq_config.time_step,
        )

        # Calculate segment parameters
        total_time = sim_config.simulation_phases.production.duration
        num_segments = sim_config.simulation_phases.segments
        seg_time = segment_time or (total_time / num_segments)
        seg_frames = segment_frames or (
            sim_config.simulation_phases.production.samples // num_segments
        )

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

    except Exception as e:
        click.echo(f"Simulation failed: {e}", err=True)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


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
