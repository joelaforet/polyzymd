"""
PolyzyMD Command Line Interface.

This module provides the main CLI entry point for PolyzyMD, using Click
for argument parsing and command organization.

Usage:
    polyzymd --help
    polyzymd build --config simulation.yaml --replicates 1-3
    polyzymd run --config simulation.yaml --replicates 1-3 --engine openmm
    polyzymd submit --config simulation.yaml --replicates 1-5
    polyzymd run-segment --config simulation.yaml --replicate 1
"""

from __future__ import annotations

import json
import logging
import sys
import tempfile
from pathlib import Path
from typing import Any

import click
import yaml
from pydantic import ValidationError

from polyzymd.cli.colors import colored_echo, echo_logo, setup_colored_logging
from polyzymd.core.branding import prepend_file_header

# Bootstrap a minimal root handler so suppress_openff_logs() works at import
# time.  setup_colored_logging() replaces this handler when the CLI runs.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
LOGGER = logging.getLogger("polyzymd")


def _echo_branding() -> None:
    """Print the PolyzyMD ASCII logo for top-level user-facing commands."""
    echo_logo()


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


def _resolve_replicates_option(
    replicates: str | None,
    replicate: int | None,
    command_name: str,
) -> list[int]:
    """Resolve --replicates / --replicate into a list of replicate numbers.

    Parameters
    ----------
    replicates : str or None
        Value from ``--replicates`` (range syntax, e.g. ``"1-3"``)
    replicate : int or None
        Value from deprecated ``--replicate`` (single integer)
    command_name : str
        Name of the CLI command (for error messages)

    Returns
    -------
    list[int]
        Sorted, deduplicated list of replicate numbers

    Raises
    ------
    click.UsageError
        If both ``--replicates`` and ``--replicate`` are given
    """
    from polyzymd.utils.replicates import parse_replicate_range

    if replicates is not None and replicate is not None:
        raise click.UsageError(
            f"Cannot use both --replicates and --replicate in '{command_name}'. "
            "Use --replicates (e.g., --replicates 1-3)."
        )

    if replicate is not None:
        if replicate <= 0:
            raise click.BadParameter(
                f"Replicate must be a positive integer, got {replicate}.",
                param_hint="'--replicate'",
            )
        click.echo(
            f"Warning: --replicate is deprecated in '{command_name}', use --replicates instead.",
            err=True,
        )
        return [replicate]

    if replicates is not None:
        try:
            return parse_replicate_range(replicates)
        except ValueError as e:
            raise click.BadParameter(str(e), param_hint="'--replicates'") from e

    # Default to a single replicate
    return [1]


def _resolve_engine_name(sim_config: object, override: str | None = None) -> str:
    """Resolve the simulation engine from CLI override or config.

    Parameters
    ----------
    sim_config : object
        Simulation configuration with optional ``engine`` attribute.
    override : str or None
        CLI-provided engine override (takes priority).

    Returns
    -------
    str
        Resolved engine name (lowercase).
    """
    if override:
        return override.lower()
    return getattr(sim_config, "engine", "openmm") or "openmm"


def _generate_system_prefix(sim_config: object) -> str:
    """Generate a system filename prefix from simulation config.

    Replicates ``GromacsExporter._generate_prefix`` so CLI checks use the
    same naming convention as build and submit workflows.

    Parameters
    ----------
    sim_config : object
        Simulation configuration object.

    Returns
    -------
    str
        System prefix (e.g. ``"CALB_SBMA-EGMA"``).
    """
    parts: list[str] = []

    enzyme = getattr(sim_config, "enzyme", None)
    enzyme_name = getattr(enzyme, "name", None)
    if isinstance(enzyme_name, str) and enzyme_name:
        parts.append(enzyme_name)

    polymers = getattr(sim_config, "polymers", None)
    polymers_enabled = getattr(polymers, "enabled", False)
    polymer_prefix = getattr(polymers, "type_prefix", None)
    if isinstance(polymers_enabled, bool) and polymers_enabled and isinstance(polymer_prefix, str):
        parts.append(polymer_prefix)

    return "_".join(parts) if parts else "system"


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


@cli.command(
    help=(
        "Build simulation input files (OpenMM, GROMACS, LAMMPS, AMBER) without "
        "running. Use --format to select the output engine. "
        "Use run --engine <gromacs|openmm> to execute locally, or submit for "
        "OpenMM SLURM jobs."
    )
)
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
    default=None,
    type=str,
    help="Replicate range (e.g., '1', '1-3', '1,3,5'). Default: 1",
)
@click.option(
    "--replicate",
    default=None,
    type=int,
    hidden=True,
    help="[Deprecated] Use --replicates instead.",
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
    "--format",
    "export_format",
    default=None,
    type=click.Choice(["gromacs", "lammps", "amber"], case_sensitive=False),
    help=(
        "Export format: gromacs, lammps (planned), amber (planned). Default: OpenMM (no export)."
    ),
)
@click.option(
    "--gromacs",
    is_flag=True,
    hidden=True,
    help="[Deprecated] Use --format gromacs instead.",
)
def build(
    config: str,
    replicates: str | None,
    replicate: int | None,
    output_dir: str | None,
    scratch_dir: str | None,
    projects_dir: str | None,
    dry_run: bool,
    export_format: str | None,
    gromacs: bool,
) -> None:
    """Build simulation input files from configuration.

    Loads the YAML configuration, constructs the molecular system
    (enzyme, substrate, polymers, solvent), and writes engine-ready input
    artifacts for one or more replicates. No simulation is executed.

    By default, this prepares OpenMM inputs in the working directory. Use
    ``--format gromacs`` to export GROMACS files (``.gro``, ``.top``,
    ``.itp``, ``.mdp``). AMBER and LAMMPS export are not yet supported.

    Use ``run --engine gromacs`` if you want PolyzyMD to build and then
    execute the full local GROMACS workflow. Use ``run --engine openmm`` for
    local OpenMM execution, or ``submit`` for OpenMM self-resubmitting SLURM
    workflows.

    The ``--replicates`` option accepts range syntax (for example ``1-3`` or
    ``1,3,5``). Each replicate is built independently with a different polymer
    random seed.

    \b
    Export Notes:
        - Output files are placed in replicate_<n>/<format>/
        - Filenames are derived from config: {enzyme_name}_{polymer_prefix}.*
        - MDP files include energy minimization, equilibration, and production stages
        - Topology is split into .itp files for cleaner multi-component systems
    """
    from pydantic import ValidationError as PydanticValidationError

    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.config.schema import SimulationConfig

    # Resolve export format: --format takes priority, --gromacs is deprecated alias
    if gromacs and export_format is not None:
        raise click.UsageError(
            "Cannot use both --gromacs and --format. Use --format gromacs instead."
        )
    if gromacs:
        click.echo(
            "Warning: --gromacs is deprecated, use --format gromacs instead.",
            err=True,
        )
        export_format = "gromacs"

    replicate_list = _resolve_replicates_option(replicates, replicate, "build")

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
            colored_echo("=" * 60, phase="build")
            colored_echo("DRY RUN — Validation Report", phase="build")
            colored_echo("=" * 60, phase="build")
            colored_echo(phase="build")

            colored_echo("Configuration Summary:", phase="build")
            colored_echo(f"  Name: {sim_config.name}", phase="build")
            if sim_config.description:
                colored_echo(f"  Description: {sim_config.description}", phase="build")
            colored_echo(f"  Config file: {config}", phase="build")
            colored_echo(phase="build")

            colored_echo("Replicates:", phase="build")
            colored_echo(f"  Count: {len(replicate_list)}", phase="build")
            colored_echo(f"  IDs: {replicate_list}", phase="build")
            colored_echo(f"  Polymer seeds: {replicate_list} (one per replicate)", phase="build")
            colored_echo(phase="build")

            colored_echo("System Components:", phase="build")
            colored_echo(f"  Chain A (Protein): {sim_config.enzyme.name}", phase="build")
            colored_echo(f"    PDB: {sim_config.enzyme.pdb_path}", phase="build")
            if sim_config.substrate:
                colored_echo(f"  Chain B (Substrate): {sim_config.substrate.name}", phase="build")
                colored_echo(f"    SDF: {sim_config.substrate.sdf_path}", phase="build")
            else:
                colored_echo("  Chain B (Substrate): none (apo system)", phase="build")
            if sim_config.polymers and sim_config.polymers.enabled:
                colored_echo(
                    f"  Chain C (Polymer): {sim_config.polymers.type_prefix}", phase="build"
                )
                colored_echo(f"    Count: {sim_config.polymers.count}", phase="build")
                colored_echo(f"    Length: {sim_config.polymers.length} monomers", phase="build")
                for monomer in sim_config.polymers.monomers:
                    colored_echo(
                        f"    Monomer: {monomer.label} ({monomer.probability * 100:.0f}%)",
                        phase="build",
                    )
            else:
                colored_echo("  Chain C (Polymer): none (no polymer)", phase="build")
            colored_echo(f"  Chain D+ (Solvent): {sim_config.solvent.primary.model}", phase="build")
            colored_echo(f"    Box padding: {sim_config.solvent.box.padding} nm", phase="build")
            colored_echo(
                f"    NaCl concentration: {sim_config.solvent.ions.nacl_concentration} M",
                phase="build",
            )
            colored_echo(phase="build")

            colored_echo("Parameterization Plan:", phase="build")
            colored_echo(f"  Protein FF: {sim_config.force_field.protein}", phase="build")
            colored_echo(
                f"  Small molecule FF: {sim_config.force_field.small_molecule}", phase="build"
            )
            colored_echo(f"  Water model: {sim_config.solvent.primary.model}", phase="build")
            colored_echo(phase="build")

            colored_echo("Thermodynamics:", phase="build")
            colored_echo(f"  Temperature: {sim_config.thermodynamics.temperature} K", phase="build")
            pressure = sim_config.thermodynamics.pressure
            if pressure is not None:
                colored_echo(f"  Pressure: {pressure} atm", phase="build")
            colored_echo(phase="build")

            colored_echo("Simulation Phases:", phase="build")
            eq_stages = sim_config.simulation_phases.equilibration_stages
            if eq_stages:
                colored_echo(f"  Equilibration: {len(eq_stages)} stage(s)", phase="build")
                for i, stage in enumerate(eq_stages, 1):
                    colored_echo(
                        f"    Stage {i}: {stage.duration} ns, {stage.ensemble}",
                        phase="build",
                    )
            colored_echo(
                f"  Production: {sim_config.simulation_phases.production.duration} ns",
                phase="build",
            )
            colored_echo(
                f"  Samples: {sim_config.simulation_phases.production.samples}",
                phase="build",
            )
            colored_echo(phase="build")

            colored_echo("Directories:", phase="build")
            colored_echo(f"  Projects: {sim_config.output.projects_directory}", phase="build")
            colored_echo(
                f"  Scratch: {sim_config.output.effective_scratch_directory}", phase="build"
            )
            colored_echo(phase="build")

            colored_echo("Per-Replicate Output:", phase="build")
            for rep in replicate_list:
                working_dir = sim_config.get_working_directory(rep)
                colored_echo(f"  Replicate {rep}:", phase="build")
                colored_echo(f"    Working dir: {working_dir}", phase="build")
                if export_format:
                    export_dir = sim_config.get_working_directory(rep) / export_format
                    colored_echo(f"    Export dir:  {export_dir}", phase="build")
            colored_echo(phase="build")

            if export_format:
                colored_echo(f"Files to Generate ({export_format.upper()}):", phase="build")
                colored_echo("  Per replicate:", phase="build")
                if export_format == "gromacs":
                    colored_echo("    - *.gro (coordinates)", phase="build")
                    colored_echo("    - *.top (topology)", phase="build")
                    colored_echo("    - *.itp (molecule parameters)", phase="build")
                    colored_echo("    - em.mdp (energy minimization)", phase="build")
                    if eq_stages:
                        for i, stage in enumerate(eq_stages, 1):
                            colored_echo(
                                f"    - eq_{i:02d}_{stage.name}.mdp (equilibration)",
                                phase="build",
                            )
                    colored_echo("    - prod.mdp (production)", phase="build")
                    colored_echo(
                        "    - Position restraints appended to molecule *.itp files",
                        phase="build",
                    )
                    colored_echo("    - run_*_gromacs.sh (run script)", phase="build")
                elif export_format in ("lammps", "amber"):
                    colored_echo(
                        f"    ({export_format.upper()} export is not yet supported)",
                        phase="build",
                    )
            else:
                colored_echo("Files to Generate (OpenMM):", phase="build")
                colored_echo("  Per replicate:", phase="build")
                colored_echo("    - solvated_system.pdb (topology + positions)", phase="build")
                colored_echo("    - system.xml (OpenMM system with restraints)", phase="build")
            colored_echo(phase="build")

            if sim_config.restraints:
                colored_echo("Restraints:", phase="build")
                for r in sim_config.restraints:
                    status = "ENABLED" if r.enabled else "DISABLED"
                    colored_echo(
                        f"  - {r.name}: {r.type.value}, d={r.distance} Å, "
                        f"k={r.force_constant} kJ/mol/nm², {status}",
                        phase="build",
                    )
                colored_echo(phase="build")

            colored_echo("=" * 60, phase="build")
            if export_format in ("lammps", "amber"):
                colored_echo(
                    f"Validation passed. {export_format.upper()} export is not yet implemented.",
                    phase="build",
                )
            else:
                colored_echo("Validation passed. Ready to build.", phase="build")
            colored_echo("=" * 60, phase="build")
            return

        for rep in replicate_list:
            colored_echo(f"Building system for replicate {rep}...", phase="build")
            working_dir = sim_config.get_working_directory(rep)
            builder = SystemBuilder.from_config(sim_config)
            interchange = builder.build_from_config(
                config=sim_config,
                working_dir=working_dir,
                polymer_seed=rep,
            )

            # Branch based on export format
            if export_format:
                # Export to requested engine format
                from polyzymd.exporters.interchange import export_system

                colored_echo(f"Exporting to {export_format.upper()} format...", phase="export")
                export_dir = sim_config.get_working_directory(rep) / export_format
                export_result = export_system(
                    interchange=interchange,
                    config=sim_config,
                    output_dir=export_dir,
                    fmt=export_format,
                    component_info=builder.get_component_info(),
                )

                colored_echo(f"{export_format.upper()} export successful!", phase="export")
                colored_echo(f"Output directory: {export_dir}", phase="export")
                colored_echo("Files generated:", phase="export")
                colored_echo(f"  - {export_result['gro'].name} (coordinates)", phase="export")
                colored_echo(f"  - {export_result['top'].name} (topology)", phase="export")
                colored_echo("  - *.itp (molecule parameters)", phase="export")
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
                    f"To run: cd {export_dir} && ./{export_result['run_script'].name}",
                    phase="export",
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
                            f"Successfully applied {len(restraint_defs)} restraint(s)",
                            phase="build",
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
                    "Use 'polyzymd submit' to submit for HPC execution,",
                    phase="build",
                )
                colored_echo(
                    "or 'polyzymd run-segment' to run a single segment locally.",
                    phase="build",
                )

    except PydanticValidationError as e:
        colored_echo("Configuration error:", err=True, level=logging.ERROR)
        for error in e.errors():
            loc = " → ".join(str(x) for x in error["loc"])
            colored_echo(f"  {loc}: {error['msg']}", err=True, level=logging.ERROR)
        sys.exit(1)

    except FileNotFoundError as e:
        colored_echo(f"File not found: {e}", err=True, level=logging.ERROR)
        sys.exit(1)

    except ValueError as e:
        colored_echo(f"Validation error: {e}", err=True, level=logging.ERROR)
        sys.exit(1)

    except NotImplementedError as e:
        colored_echo(f"Not yet supported: {e}", err=True, level=logging.ERROR)
        sys.exit(1)

    except Exception as e:
        colored_echo(f"Unexpected error: {e}", err=True, level=logging.ERROR)
        colored_echo(
            "This may be a bug. Re-run with --verbose for details.",
            err=True,
            level=logging.ERROR,
        )
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


# =============================================================================
# Run Command
# =============================================================================


def _print_run_dry_run_report(
    sim_config: "SimulationConfig",
    config_path: str,
    replicate_list: list[int],
    engine: str,
    gmx_path: str | None,
) -> None:
    """Print a preview report for ``run --dry-run``.

    Parameters
    ----------
    sim_config : SimulationConfig
        Validated simulation configuration.
    config_path : str
        Path to the YAML configuration file.
    replicate_list : list[int]
        Replicates that would be run.
    engine : str
        Execution engine (``"gromacs"`` or ``"openmm"``).
    gmx_path : str or None
        Optional GROMACS executable path.
    """
    phase = "simulation"
    colored_echo("=" * 60, phase=phase)
    colored_echo("DRY RUN — Run Command Preview", phase=phase)
    colored_echo("=" * 60, phase=phase)
    colored_echo(phase=phase)

    colored_echo("Configuration Summary:", phase=phase)
    colored_echo(f"  Name: {sim_config.name}", phase=phase)
    if sim_config.description:
        colored_echo(f"  Description: {sim_config.description}", phase=phase)
    colored_echo(f"  Config file: {config_path}", phase=phase)
    colored_echo(f"  Engine: {engine}", phase=phase)
    if engine == "gromacs":
        colored_echo(f"  GROMACS executable: {gmx_path or 'gmx'}", phase=phase)
    colored_echo(phase=phase)

    colored_echo("Replicates:", phase=phase)
    colored_echo(f"  Count: {len(replicate_list)}", phase=phase)
    colored_echo(f"  IDs: {replicate_list}", phase=phase)
    colored_echo(phase=phase)

    colored_echo("Planned output:", phase=phase)
    for rep in replicate_list:
        working_dir = sim_config.get_working_directory(rep)
        colored_echo(f"  Replicate {rep}:", phase=phase)
        colored_echo(f"    Working dir: {working_dir}", phase=phase)
        if engine == "gromacs":
            gromacs_dir = sim_config.get_working_directory(rep) / "gromacs"
            colored_echo(f"    GROMACS dir: {gromacs_dir}", phase=phase)
        else:
            colored_echo(
                "    Workflow: build -> minimize -> equilibrate -> production", phase=phase
            )
            colored_echo(
                "    Outputs: solvated_system.pdb, system.xml, production_0/*", phase=phase
            )

    colored_echo(phase=phase)
    colored_echo("Dry run complete. No files were written.", phase=phase)
    colored_echo("=" * 60, phase=phase)


@cli.command(
    help=(
        "Build and run simulations locally with a selected engine. "
        "Use --engine gromacs for full local GROMACS workflow, or "
        "--engine openmm for local OpenMM execution. Use --dry-run for "
        "preview-only validation without writing files."
    ),
)
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
    default=None,
    type=str,
    help="Replicate range (e.g., '1', '1-3', '1,3,5'). Default: 1",
)
@click.option(
    "--replicate",
    default=None,
    type=int,
    hidden=True,
    help="[Deprecated] Use --replicates instead.",
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
    default=None,
    help="Path to GROMACS executable (only valid with --engine gromacs)",
)
@click.option(
    "--engine",
    required=True,
    type=click.Choice(["gromacs", "openmm"], case_sensitive=False),
    help="Simulation engine to run locally",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Validate config and preview planned execution without writing files",
)
def run(
    config: str,
    replicates: str | None,
    replicate: int | None,
    scratch_dir: str | None,
    projects_dir: str | None,
    gmx_path: str | None,
    engine: str,
    dry_run: bool,
) -> None:
    """Build and run a simulation locally.

    For each replicate, this command builds the system and executes the full
    local simulation workflow using the selected engine.

    Use ``--dry-run`` to validate the configuration and preview planned output
    without building or running.

    \b
    Notes:
        - ``--engine gromacs`` requires GROMACS to be installed and accessible
        - ``--gmx-path`` is valid only with ``--engine gromacs``
        - ``--dry-run`` validates and previews without writing files
    """
    from pydantic import ValidationError as PydanticValidationError

    from polyzymd.config.schema import SimulationConfig

    engine = engine.lower()
    if engine == "openmm" and gmx_path is not None:
        raise click.UsageError("--gmx-path can only be used with --engine gromacs")

    replicate_list = _resolve_replicates_option(replicates, replicate, "run")

    colored_echo(f"Loading configuration from: {config}", phase="simulation")

    try:
        sim_config = SimulationConfig.from_yaml(config)
        colored_echo(f"Running local simulation: {sim_config.name}", phase="simulation")
        colored_echo(f"Engine: {engine}", phase="simulation")

        # Override directories if provided via CLI
        if scratch_dir:
            sim_config.output.scratch_directory = Path(scratch_dir)
        if projects_dir:
            sim_config.output.projects_directory = Path(projects_dir)

        if dry_run:
            _print_run_dry_run_report(
                sim_config=sim_config,
                config_path=config,
                replicate_list=replicate_list,
                engine=engine,
                gmx_path=gmx_path,
            )
            return

        resolved_gmx_path = gmx_path or "gmx"

        succeeded = 0
        for rep in replicate_list:
            colored_echo(
                f"\n--- Replicate {rep} ({succeeded + 1}/{len(replicate_list)}) ---",
                phase="simulation",
            )

            if engine == "gromacs":
                _run_gromacs_impl(
                    sim_config=sim_config,
                    replicate=rep,
                    gmx_path=resolved_gmx_path,
                )
            else:
                _run_openmm_impl(
                    sim_config=sim_config,
                    replicate=rep,
                )

            succeeded += 1

        if len(replicate_list) > 1:
            colored_echo(
                f"\nAll {succeeded} replicate(s) completed successfully.", phase="simulation"
            )

    except PydanticValidationError as e:
        colored_echo("Configuration error:", err=True, level=logging.ERROR)
        for error in e.errors():
            loc = " → ".join(str(x) for x in error["loc"])
            colored_echo(f"  {loc}: {error['msg']}", err=True, level=logging.ERROR)
        sys.exit(1)

    except FileNotFoundError as e:
        colored_echo(f"File not found: {e}", err=True, level=logging.ERROR)
        sys.exit(1)

    except ValueError as e:
        colored_echo(f"Validation error: {e}", err=True, level=logging.ERROR)
        sys.exit(1)

    except (RuntimeError, OSError) as e:
        colored_echo(f"Unexpected error: {e}", err=True, level=logging.ERROR)
        colored_echo(
            "This may be a bug. Re-run with --verbose for details.",
            err=True,
            level=logging.ERROR,
        )
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


def _run_gromacs_impl(
    sim_config: "SimulationConfig",
    replicate: int,
    gmx_path: str,
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
    """
    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.exporters.gromacs import GromacsError, GromacsExporter, GromacsRunner

    working_dir = sim_config.get_working_directory(replicate)
    # Determine output directory for GROMACS files
    gromacs_dir = working_dir / "gromacs"

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
    colored_echo("  - *.itp (molecule parameters)", phase="export")
    colored_echo(f"  - {export_result['em_mdp'].name} (energy minimization)", phase="export")
    for eq_mdp in export_result["eq_mdps"]:
        colored_echo(f"  - {eq_mdp.name} (equilibration)", phase="export")
    colored_echo(f"  - {export_result['prod_mdp'].name} (production)", phase="export")
    if export_result.get("posres_defines"):
        colored_echo("Position restraints added to molecule ITP files:", phase="export")
        for component, define in export_result["posres_defines"].items():
            colored_echo(f"  - {component}: #ifdef {define}", phase="export")
    colored_echo(f"  - {export_result['run_script'].name} (run script)", phase="export")

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


def _run_openmm_impl(
    sim_config: "SimulationConfig",
    replicate: int,
) -> None:
    """Build and run a full local OpenMM simulation.

    Parameters
    ----------
    sim_config : SimulationConfig
        Validated simulation configuration.
    replicate : int
        Replicate number.
    """
    production = sim_config.simulation_phases.production
    working_dir = sim_config.get_working_directory(replicate)

    colored_echo(f"Building and running OpenMM in {working_dir}", phase="simulation")
    _run_initial_segment(
        sim_config=sim_config,
        working_dir=working_dir,
        replicate=replicate,
        skip_build=False,
        duration_ns=production.duration,
        num_samples=production.samples,
        timestep_fs=production.time_step,
    )
    colored_echo("OpenMM simulation completed successfully.", phase="simulation")
    colored_echo(f"Output directory: {working_dir}", phase="simulation")


# =============================================================================
# Submit Command (SLURM)
# =============================================================================


def _print_gromacs_dry_run_details(
    sim_config: "SimulationConfig",
    preset: str,
    replicate_list: list[int],
    time_limit: str | None,
    memory: str | None,
    account: str | None,
    gpu_type: str | None,
    constraint: str | None,
    nodelist: str | None,
    partition: str | None = None,
    qos: str | None = None,
    email: str = "",
) -> None:
    """Print GROMACS-specific dry-run details with effective SLURM configuration.

    Parameters
    ----------
    sim_config : SimulationConfig
        Validated simulation configuration.
    preset : str
        SLURM partition preset name.
    replicate_list : list[int]
        Replicates that would be submitted.
    time_limit : str or None
        CLI time-limit override.
    memory : str or None
        CLI memory override.
    account : str or None
        CLI account override.
    gpu_type : str or None
        CLI GPU type override.
    constraint : str or None
        CLI constraint override.
    nodelist : str or None
        CLI nodelist override.
    partition : str or None, optional
        CLI partition override.
    qos : str or None, optional
        CLI QoS override.
    email : str, optional
        CLI email override.
    """
    from polyzymd.engines import create_engine
    from polyzymd.workflow.slurm import SlurmConfig

    phase = "workflow"
    engine_impl = create_engine(sim_config, override="gromacs", defer_binary=True)
    gromacs_cfg = sim_config.gromacs

    base_slurm = SlurmConfig.from_preset(preset)
    if time_limit:
        base_slurm.time_limit = time_limit
    if memory:
        base_slurm.memory = memory
    if account:
        base_slurm.account = account
    if partition:
        base_slurm.partition = partition
    if qos:
        base_slurm.qos = qos
    if email:
        base_slurm.email = email
    if gpu_type:
        base_slurm.gpu_type = gpu_type
    if constraint:
        base_slurm.constraint = constraint
    if nodelist:
        base_slurm.nodelist = nodelist

    effective = engine_impl._resolve_slurm_config(base_slurm)
    effective_flags = engine_impl._resolve_mdrun_flags(effective)

    colored_echo("  SLURM configuration (effective):", phase=phase)
    colored_echo(f"    Partition:      {effective.partition}", phase=phase)
    colored_echo(f"    Time limit:     {effective.time_limit}", phase=phase)
    colored_echo(f"    Memory:         {effective.memory}", phase=phase)
    colored_echo(f"    Account:        {effective.account or '(none)'}", phase=phase)
    colored_echo(f"    Email:          {effective.email or '(none)'}", phase=phase)
    colored_echo(f"    Nodes:          {effective.nodes}", phase=phase)
    colored_echo(f"    Tasks:          {effective.ntasks}", phase=phase)
    colored_echo(f"    CPUs/task:      {effective.cpus_per_task}", phase=phase)
    colored_echo(f"    GPUs:           {effective.gpus}", phase=phase)
    if effective.constraint:
        colored_echo(f"    Constraint:     {effective.constraint}", phase=phase)
    if effective.nodelist:
        colored_echo(f"    Nodelist:       {effective.nodelist}", phase=phase)
    if effective.qos:
        colored_echo(f"    QoS:            {effective.qos}", phase=phase)
    colored_echo(phase=phase)

    colored_echo("  GROMACS configuration:", phase=phase)
    colored_echo(f"    Binary:         {gromacs_cfg.gmx_binary or '(auto-detect)'}", phase=phase)
    colored_echo(f"    GPU mode:       {'yes' if gromacs_cfg.gpu else 'no'}", phase=phase)
    colored_echo(f"    ntmpi:          {gromacs_cfg.ntmpi}", phase=phase)
    colored_echo(f"    ntomp:          {gromacs_cfg.ntomp}", phase=phase)
    colored_echo(f"    mdrun flags:    {effective_flags or '(none)'}", phase=phase)
    if gromacs_cfg.module_load:
        colored_echo(f"    Module load:    {gromacs_cfg.module_load}", phase=phase)
    colored_echo(phase=phase)

    colored_echo("  Per-replicate output:", phase=phase)
    for rep in replicate_list:
        working_dir = sim_config.get_working_directory(rep) / "gromacs"
        script_dir = working_dir / "daisy_chain_scripts"
        script_name = f"run_rep{rep}.sh"
        colored_echo(f"    Rep {rep}:", phase=phase)
        colored_echo(f"      Working dir: {working_dir}", phase=phase)
        colored_echo(f"      Script:      {script_dir / script_name}", phase=phase)


@cli.command(
    help=(
        "Submit self-resubmitting SLURM jobs for OpenMM or GROMACS. "
        "Use --engine to override the engine configured in YAML."
    )
)
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
    help="Preview submission plan only (no files written, no submission)",
)
@click.option(
    "--generate-only",
    is_flag=True,
    help="Generate job scripts without submitting to SLURM",
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
    "--partition",
    default=None,
    help="Override SLURM partition",
)
@click.option(
    "--qos",
    default=None,
    help="Override SLURM QoS",
)
@click.option(
    "--gpu-type",
    default=None,
    help="GPU type for GRES (e.g., 'a100', 'a40', 'v100-32', 'mi100').",
)
@click.option(
    "--constraint",
    default=None,
    help=(
        "SLURM --constraint expression for node feature selection. "
        "Supports boolean expressions: 'A40' (single), 'A40|A100' (OR), 'avx2&rh8' (AND)."
    ),
)
@click.option(
    "--nodelist",
    default=None,
    help="SLURM --nodelist override (e.g., 'node01' or 'node[01-04]').",
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
@click.option(
    "--engine",
    default=None,
    type=click.Choice(["gromacs", "openmm"], case_sensitive=False),
    help="Engine override for submission backend (default: from config or openmm)",
)
@click.option(
    "--force",
    is_flag=True,
    help="Skip duplicate-job check and submit even if a SLURM job is already running for the replicate",
)
def submit(
    config: str,
    replicates: str,
    scratch_dir: str | None,
    projects_dir: str | None,
    preset: str,
    email: str,
    dry_run: bool,
    generate_only: bool,
    output_dir: str | None,
    time_limit: str | None,
    memory: str | None,
    account: str | None,
    partition: str | None,
    qos: str | None,
    gpu_type: str | None,
    constraint: str | None,
    nodelist: str | None,
    submit_openff_logs: bool,
    skip_build: bool,
    pixi_env: str | None,
    engine: str | None,
    force: bool,
) -> None:
    """Submit OpenMM or GROMACS simulation jobs to SLURM.

    Creates and optionally submits one self-resubmitting job per replicate.
    OpenMM submission uses the existing daisy-chain flow, while GROMACS
    submission uses the engine submission interface.

    \b
    Directory structure:
        - projects_dir: Where job scripts and SLURM logs are stored (long-term storage)
        - scratch_dir: Where simulation data is written (high-performance storage)
    """
    from polyzymd.workflow.daisy_chain import submit_daisy_chain
    from polyzymd.workflow.slurm import PRESET_DEFAULT_PIXI_ENV

    if dry_run and generate_only:
        raise click.UsageError("Cannot use both --dry-run and --generate-only")

    # Resolve pixi environment: explicit flag > preset default
    resolved_pixi_env = pixi_env or PRESET_DEFAULT_PIXI_ENV.get(preset, "cuda-12-4")

    _echo_branding()
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
    if partition:
        colored_echo(f"Partition override: {partition}", phase="workflow")
    if qos:
        colored_echo(f"QoS override: {qos}", phase="workflow")
    if memory:
        colored_echo(f"Memory allocation: {memory}", phase="workflow")
    if gpu_type:
        colored_echo(f"GPU type override: {gpu_type}", phase="workflow")
    if constraint:
        colored_echo(f"Constraint: {constraint}", phase="workflow")
    if skip_build:
        colored_echo("Skip-build mode: using pre-built systems", phase="workflow")

    if dry_run:
        from polyzymd.config.schema import SimulationConfig

        replicate_list = _resolve_replicates_option(replicates, None, "submit")
        sim_config = SimulationConfig.from_yaml(config)
        engine_name = _resolve_engine_name(sim_config, override=engine)
        if scratch_dir:
            sim_config.output.scratch_directory = Path(scratch_dir)
        if projects_dir:
            sim_config.output.projects_directory = Path(projects_dir)

        if engine_name == "gromacs" and submit_openff_logs:
            colored_echo(
                "Warning: --openff-logs has no effect with --engine gromacs",
                phase="workflow",
                level=logging.WARNING,
            )

        colored_echo("=" * 60, phase="workflow")
        colored_echo("DRY RUN — Submission Preview", phase="workflow")
        colored_echo("=" * 60, phase="workflow")
        colored_echo(phase="workflow")
        colored_echo(f"  Simulation:  {sim_config.name}", phase="workflow")
        colored_echo(f"  Engine:      {engine_name}", phase="workflow")
        colored_echo(f"  Preset:      {preset}", phase="workflow")
        colored_echo(f"  Pixi env:    {resolved_pixi_env}", phase="workflow")
        colored_echo(f"  Replicates:  {replicate_list}", phase="workflow")
        colored_echo(phase="workflow")

        if engine_name == "gromacs":
            _print_gromacs_dry_run_details(
                sim_config=sim_config,
                preset=preset,
                replicate_list=replicate_list,
                time_limit=time_limit,
                memory=memory,
                account=account,
                gpu_type=gpu_type,
                constraint=constraint,
                nodelist=nodelist,
                partition=partition,
                qos=qos,
                email=email,
            )
        else:
            script_dir = (
                Path(output_dir) if output_dir else sim_config.output.get_job_scripts_directory()
            )
            colored_echo(f"  Script dir:  {script_dir}", phase="workflow")

        colored_echo(phase="workflow")
        colored_echo("Dry run complete. No files were written.", phase="workflow")
        colored_echo("=" * 60, phase="workflow")
        return

    from polyzymd.config.schema import SimulationConfig

    sim_config = SimulationConfig.from_yaml(config)
    engine_name = _resolve_engine_name(sim_config, override=engine)

    if scratch_dir:
        sim_config.output.scratch_directory = Path(scratch_dir)
    if projects_dir:
        sim_config.output.projects_directory = Path(projects_dir)

    if engine_name == "gromacs":
        from polyzymd.engines import create_engine
        from polyzymd.engines.base import EngineSubmitRequest
        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs, create_job_name
        from polyzymd.workflow.slurm import SlurmConfig

        if submit_openff_logs:
            colored_echo(
                "Warning: --openff-logs has no effect with --engine gromacs",
                phase="workflow",
                level=logging.WARNING,
            )

        engine_impl = create_engine(sim_config, override="gromacs", defer_binary=True)
        replicate_list = _resolve_replicates_option(replicates, None, "submit")
        config_path_abs = Path(config).resolve()

        colored_echo("Using GROMACS submission backend", phase="workflow")

        for rep in replicate_list:
            slurm_config = SlurmConfig.from_preset(preset)
            if email:
                slurm_config.email = email
            if time_limit:
                slurm_config.time_limit = time_limit
            if memory:
                slurm_config.memory = memory
            if account:
                slurm_config.account = account
            if partition:
                slurm_config.partition = partition
            if qos:
                slurm_config.qos = qos
            if gpu_type:
                slurm_config.gpu_type = gpu_type
            if constraint:
                slurm_config.constraint = constraint
            if nodelist:
                slurm_config.nodelist = nodelist

            working_dir = sim_config.get_working_directory(rep) / "gromacs"
            job_name = create_job_name(sim_config, rep)

            if not force:
                existing = check_existing_slurm_jobs(job_name)
                if existing:
                    ids = ", ".join(existing)
                    colored_echo(
                        f"Replicate {rep} already has RUNNING/PENDING SLURM "
                        f"job(s): {ids}. Use --force to submit anyway.",
                        err=True,
                        phase="workflow",
                        level=logging.ERROR,
                    )
                    continue

            request = EngineSubmitRequest(
                replicate=rep,
                config_path=config_path_abs,
                working_dir=working_dir,
                slurm_config=slurm_config,
                job_name=job_name,
                extra={"pixi_env": resolved_pixi_env, "skip_build": skip_build},
            )

            if generate_only:
                script_path = engine_impl.prepare_submission(request)
                colored_echo(f"  Rep {rep}: script at {script_path}", phase="workflow")
            else:
                result = engine_impl.submit(request)
                if result.get("submitted"):
                    colored_echo(f"  Rep {rep}: {result['stdout']}", phase="workflow")
                else:
                    colored_echo(
                        f"  Rep {rep}: script at {result['script_path']} (sbatch not available)",
                        phase="workflow",
                    )

        if not generate_only:
            colored_echo("\nGROMACS job submission complete!", phase="workflow")
            colored_echo("Monitor with: squeue -u $USER", phase="workflow")
        return

    if generate_only:
        colored_echo(
            "GENERATE-ONLY MODE - scripts will be created but not submitted", phase="workflow"
        )

    try:
        submit_daisy_chain(
            config_path=config,
            slurm_preset=preset,
            replicates=replicates,
            email=email,
            dry_run=dry_run,
            generate_only=generate_only,
            force=force,
            pixi_env=resolved_pixi_env,
            output_dir=output_dir,
            scratch_dir=scratch_dir,
            projects_dir=projects_dir,
            time_limit=time_limit,
            memory=memory,
            account=account,
            partition=partition,
            qos=qos,
            gpu_type=gpu_type,
            constraint=constraint,
            nodelist=nodelist,
            openff_logs=submit_openff_logs,
            skip_build=skip_build,
        )

        if not generate_only:
            colored_echo("\nJob submission complete!", phase="workflow")
            colored_echo("Monitor with: squeue -u $USER", phase="workflow")

    except (FileNotFoundError, ValueError, RuntimeError, OSError) as e:
        colored_echo(f"Submission failed: {e}", err=True, level=logging.ERROR)
        if LOGGER.level == logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)


# =============================================================================
# Run-segment Command — unified entry point for self-resubmitting jobs
# =============================================================================


@cli.command("run-segment", hidden=True)
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
    scratch_dir: str | None,
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
        _derive_overall_status,
        get_next_segment_info,
        load_or_scan_progress,
        save_progress,
    )

    colored_echo(f"Loading configuration from: {config}", phase="simulation")

    try:
        sim_config = SimulationConfig.from_yaml(config)
    except (FileNotFoundError, yaml.YAMLError, ValidationError, ValueError) as e:
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
        progress.status = _derive_overall_status(
            progress.segments, is_complete=progress.is_complete
        )
        save_progress(working_dir, progress)

    # ---- Handle RUNNING segments (concurrency guard) ----
    # If any segment appears to still be executing (recent checkpoint),
    # refuse to start a new segment to prevent concurrent execution and
    # the associated data-loss / overlap bugs.  Exit with a dedicated
    # code (EXIT_CODE_CONCURRENT = 2) so the SLURM bash wrapper knows
    # this is a duplicate chain and terminates WITHOUT resubmitting.
    #
    # Previously this exited with code 0, which caused the resubmission
    # logic to call check-progress (work remains → exit 1) and resubmit,
    # creating an infinite submit-cancel-resubmit loop.
    from polyzymd.simulation.signals import EXIT_CODE_CONCURRENT

    running_segments = [s for s in progress.segments if s.status == SegmentStatus.RUNNING]
    if running_segments:
        indices = ", ".join(str(s.index) for s in running_segments)
        colored_echo(
            f"Segment(s) {indices} appear(s) to still be running "
            f"(checkpoint written recently). Refusing to start a new "
            f"segment to avoid concurrent execution — this duplicate "
            f"chain will terminate without resubmitting.",
            phase="simulation",
            level=logging.WARNING,
        )
        sys.exit(EXIT_CODE_CONCURRENT)

    # ---- Handle hard-killed segments (no INTERRUPTED marker) ----
    # When SLURM preempts a job with SIGKILL (no grace period) the
    # graceful shutdown handler never runs: no ``INTERRUPTED`` marker,
    # no ``interrupted_state.xml``.  The filesystem scanner classifies
    # these segments as INTERRUPTED via the stale-checkpoint heuristic,
    # but the only recovery file is ``restart_state.xml`` which may be
    # at an arbitrarily early position within the segment's run.
    #
    # Advancing to a new segment index in this situation causes massive
    # data loss: the new segment loads from ``restart_state.xml`` (an
    # early wall-time checkpoint) and re-simulates work that the killed
    # segment had already completed beyond that point.  The step
    # accounting also becomes corrupted because the filesystem scanner
    # estimates steps from the CSV file rather than an authoritative
    # INTERRUPTED marker.
    #
    # Fix: clean up the hard-killed segment's directory and remove it
    # from progress so that ``get_next_segment_info()`` assigns the
    # *same* segment index.  The segment will be retried from the
    # previous good state (the last segment with a proper completion
    # or graceful interruption).
    if progress.segments:
        last_seg = max(progress.segments, key=lambda s: s.index)
        if last_seg.status == SegmentStatus.INTERRUPTED:
            last_seg_dir = working_dir / f"production_{last_seg.index}"
            interrupted_marker = last_seg_dir / "INTERRUPTED"
            if last_seg_dir.exists() and not interrupted_marker.exists():
                colored_echo(
                    f"Segment {last_seg.index} was hard-killed (no INTERRUPTED "
                    f"marker — only stale checkpoint found). Cleaning up "
                    f"directory to retry from previous good state.",
                    phase="simulation",
                    level=logging.WARNING,
                )
                shutil.rmtree(last_seg_dir)
                progress.segments = [s for s in progress.segments if s.index != last_seg.index]
                progress.status = _derive_overall_status(
                    progress.segments, is_complete=progress.is_complete
                )
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

    # ------------------------------------------------------------------
    # Fast-path: skip minimize + equilibration when recovering a job that
    # already completed equilibration but was interrupted before starting
    # production.  Without this, the runner would re-run all equilibration
    # stages (or worse, rebuild the system with a different atom count).
    # ------------------------------------------------------------------
    if skip_build:
        from polyzymd.simulation.progress import load_progress

        progress = load_progress(working_dir)
        if progress is not None and progress.equilibration_complete:
            phases = sim_config.simulation_phases
            eq_stages_cfg = phases.equilibration_stages
            last_idx = len(eq_stages_cfg) - 1
            last_name = eq_stages_cfg[last_idx].name

            colored_echo(
                f"Equilibration already complete "
                f"({progress.num_eq_stages_completed} stage(s)) "
                "— skipping minimize + equilibrate, jumping to production",
                phase="simulation",
            )
            runner._load_eq_stage_state(last_idx, last_name)

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
            return

    # Minimize
    colored_echo("Running energy minimization...", phase="simulation")
    runner.minimize()

    # Equilibrate
    phases = sim_config.simulation_phases
    eq_duration = phases.total_equilibration_duration
    colored_echo(f"Running staged equilibration: {eq_duration:.3f} ns...", phase="simulation")
    eq_result = runner.run_equilibration(temperature=temperature, config=phases)

    # Save equilibration progress so a resubmitted job knows eq is done
    from datetime import datetime, timezone

    from polyzymd.simulation.progress import (
        EquilibrationStageRecord,
        SegmentStatus,
        save_progress,
    )
    from polyzymd.simulation.progress import (
        load_progress as _load_progress,
    )

    progress = _load_progress(working_dir)
    if progress is not None:
        eq_stages = []
        if eq_result.get("type") == "staged_equilibration":
            now_iso = datetime.now(timezone.utc).isoformat()
            for stage_info in eq_result.get("stages", []):
                eq_stages.append(
                    EquilibrationStageRecord(
                        index=stage_info["stage_index"],
                        name=stage_info["stage_name"],
                        status=SegmentStatus.COMPLETED,
                        duration_ns=stage_info["duration_ns"],
                        ensemble=stage_info.get("ensemble", "NVT"),
                        finished_at=now_iso,
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


@cli.command("check-progress", hidden=True)
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
    "--engine",
    default=None,
    type=click.Choice(["gromacs", "openmm"], case_sensitive=False),
    hidden=True,
)
def check_progress(
    config: str,
    replicate: int,
    scratch_dir: str | None,
    engine: str | None,
) -> None:
    """Check whether a simulation is complete.

    Loads progress state and exits with code 0 if the simulation is
    complete, or code 1 if work remains. Used by SLURM resubmission
    logic to decide whether to resubmit.

    \b
    Exit codes:
        0 - Simulation complete (do NOT resubmit)
        1 - Work remains (resubmit)
        3 - Error (do NOT resubmit)
    """
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.engines import create_engine
    from polyzymd.simulation.progress import save_progress
    from polyzymd.simulation.signals import EXIT_CODE_CHECK_ERROR

    try:
        sim_config = SimulationConfig.from_yaml(config)
    except (FileNotFoundError, yaml.YAMLError, ValidationError, ValueError) as e:
        colored_echo(f"Failed to load config: {e}", err=True, level=logging.ERROR)
        sys.exit(EXIT_CODE_CHECK_ERROR)

    engine_name = _resolve_engine_name(sim_config, override=engine)
    engine_inst = create_engine(sim_config, override=engine_name, defer_binary=True)

    if scratch_dir:
        working_dir = Path(scratch_dir)
    else:
        working_dir = engine_inst.get_engine_working_directory(sim_config, replicate)

    prod = sim_config.simulation_phases.production
    timestep_fs = prod.time_step
    try:
        progress = engine_inst.load_or_scan_progress(working_dir, replicate)
        save_progress(working_dir, progress)
    except (FileNotFoundError, ValueError, OSError) as e:
        colored_echo(f"Failed to load progress: {e}", err=True, level=logging.ERROR)
        sys.exit(EXIT_CODE_CHECK_ERROR)

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
    from polyzymd.engines import create_engine
    from polyzymd.simulation.progress import SimulationStatus, save_progress

    try:
        sim_config = SimulationConfig.from_yaml(config)
    except (FileNotFoundError, yaml.YAMLError, ValidationError, ValueError) as e:
        click.echo(click.style(f"Error: Failed to load config: {e}", fg="red"), err=True)
        sys.exit(1)

    engine_name = _resolve_engine_name(sim_config, override=None)
    engine_inst = create_engine(sim_config, override=engine_name, defer_binary=True)

    # Total production metadata from config
    prod = sim_config.simulation_phases.production
    total_ns = prod.duration

    # Build a human-readable system name from the directory template
    # (format with replicate=1, then strip the trailing "_run1")
    dir_name = sim_config._format_run_directory_name(1)
    system_name = dir_name.rsplit("_run", 1)[0] if "_run" in dir_name else dir_name

    # Discover replicate directories (scratch-based, works for both engines)
    replicates = sim_config.discover_replicate_dirs()

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
    max_rep = max(rep_map)
    label_width = len(f"run{max_rep}")

    need_attention = 0
    completed_count = 0
    running_count = 0

    for rep_num, rep_path in sorted(rep_map.items()):
        label = f"run{rep_num}"

        if rep_path is None:
            # Directory not found on disk
            frac = 0.0
            completed_ns = 0.0
            status_str = "not_found"
            status_display = "not found"
        else:
            engine_dir = engine_inst.resolve_engine_working_directory(rep_path)
            progress = engine_inst.load_or_scan_progress(engine_dir, rep_num)
            save_progress(engine_dir, progress)

            frac = progress.fraction_complete()
            # Compute ns from total steps (not time_completed_ns which
            # only counts COMPLETED segments, ignoring INTERRUPTED/RUNNING).
            completed_ns = (progress.total_steps_completed * progress.timestep_fs) / 1e6
            status_val = progress.status
            status_str = status_val.value
            status_display = status_val.value

        bar = render_progress_bar(frac, status_str)
        pct = frac * 100

        # Count replicates by category
        if status_str == "completed":
            completed_count += 1
        elif status_str == "running":
            running_count += 1
        else:
            need_attention += 1

        # Format: "  run1  ████░░░░  100.0%  100.0/100.0 ns  completed"
        click.echo(
            f"  {label:<{label_width}}  {bar}  {pct:5.1f}%  "
            f"{completed_ns:6.1f}/{total_ns:.1f} ns  {status_display}"
        )

    click.echo()
    total_reps = len(rep_map)
    if completed_count == total_reps:
        click.echo(click.style(f"  All {total_reps} replicates completed!", fg="green"))
    else:
        if need_attention > 0:
            click.echo(
                f"  {need_attention}/{total_reps} need attention "
                f"(recover with: polyzymd recover -c {config} -r <N> --submit)"
            )
        if running_count > 0:
            click.echo(f"  {completed_count}/{total_reps} completed, {running_count} still running")
        if need_attention == 0 and running_count == 0:
            click.echo(f"  {completed_count}/{total_reps} completed")
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
        eq_stages = sim_config.simulation_phases.equilibration_stages
        colored_echo(
            f"  Equilibration: {sim_config.simulation_phases.total_equilibration_duration} ns "
            f"across {len(eq_stages)} stage(s)"
        )
        for stage in eq_stages:
            colored_echo(f"    - {stage.name}: {stage.duration} ns ({stage.ensemble.value})")
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
    except (yaml.YAMLError, ValidationError, ValueError) as e:
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
        _echo_branding()
        colored_echo(f"Creating project directory: {name}/")
        project_dir.mkdir(parents=True)
        (project_dir / "structures").mkdir()
        (project_dir / "job_scripts").mkdir()
        (project_dir / "slurm_logs").mkdir()

        # Copy template configuration
        template_path = resources.files("polyzymd.templates.templates").joinpath(
            "config_template.yaml"
        )
        config_dest = project_dir / "config.yaml"
        with resources.as_file(template_path) as template_file:
            shutil.copy(template_file, config_dest)

        config_dest.write_text(prepend_file_header(config_dest.read_text(), comment_prefix="#"))

        # Create placeholder files
        protein_placeholder = project_dir / "structures" / "place_protein_here.placeholder.txt"
        protein_placeholder.write_text(
            prepend_file_header(
                """\
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
""",
                comment_prefix="#",
            )
        )

        ligand_placeholder = project_dir / "structures" / "place_ligand_here.placeholder.txt"
        ligand_placeholder.write_text(
            prepend_file_header(
                """\
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
""",
                comment_prefix="#",
            )
        )

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
        # Broad catch is intentional — must clean up partially-created
        # directory regardless of what went wrong
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
    help="Path for the cleaned PDB file. Defaults to <stem>_cleaned.pdb.",
)
@click.option(
    "--ph",
    default=7.4,
    type=float,
    show_default=True,
    help="Accepted for compatibility; pure-Python normalization does not add hydrogens.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    help="Plan normalization and print diagnostics without writing a PDB file.",
)
@click.option(
    "--report-json",
    "report_json",
    default=None,
    type=click.Path(dir_okay=False),
    help="Write a JSON-safe normalization plan report to this path.",
)
@click.option(
    "--check-pablo",
    is_flag=True,
    help="Validate the normalized PDB with OpenFF Pablo after clean-PDB planning succeeds.",
)
def clean_pdb(
    input_path: str,
    output_path: str | None,
    ph: float,
    dry_run: bool,
    report_json: str | None,
    check_pablo: bool,
) -> None:
    """Normalize a PDB file for PolyzyMD conjugation ingestion.

    The command performs pure-Python chain and residue normalization only. It
    accepts canonical protein residues plus explicitly linked PTM/glycan/polymer
    moieties, rejects waters, ions, solvents, and free ligands, and never
    mutates the input file.

    \b
    Examples:
        polyzymd clean-pdb -i structures/my_protein.pdb
        polyzymd clean-pdb -i raw.pdb -o cleaned.pdb --dry-run
    """
    from polyzymd.builders.conjugation.structure_normalization import (
        default_cleaned_pdb_path,
        plan_pdb_chain_normalization,
        write_normalized_pdb,
    )

    input_file = Path(input_path)
    output_file = (
        Path(output_path) if output_path is not None else default_cleaned_pdb_path(input_file)
    )
    report_path = Path(report_json) if report_json is not None else None

    colored_echo(f"Planning clean-PDB normalization: {input_file}")
    colored_echo("  Mode: pure-Python chain and residue normalization")
    colored_echo(f"  Compatibility pH option ignored: {ph:g}")

    plan = plan_pdb_chain_normalization(input_file)

    _echo_clean_pdb_plan_summary(plan, output_file, dry_run)
    if not plan.valid:
        _echo_clean_pdb_issues(plan.issues)
        _write_clean_pdb_report(
            report_path,
            plan,
            _skipped_pablo_validation("Clean-PDB structural validation failed")
            if check_pablo
            else None,
        )
        raise click.ClickException(
            "Clean-PDB validation failed; strip waters, ions, solvents, and unsupported free "
            "components before normalization. Intentional bound cofactors/metals are not yet "
            "supported by this command."
        )

    if dry_run:
        pablo_validation = None
        if check_pablo:
            temporary_path = _temporary_clean_pdb_path(output_file)
            try:
                write_normalized_pdb(input_file, temporary_path, plan)
                pablo_validation = _validate_clean_pdb_with_pablo(
                    temporary_path,
                    temporary_check_file=True,
                )
            except ValueError as exc:
                pablo_validation = _failed_pablo_validation(
                    temporary_path,
                    temporary_check_file=True,
                    error=str(exc),
                )
            finally:
                temporary_path.unlink(missing_ok=True)
        _write_clean_pdb_report(report_path, plan, pablo_validation)
        if pablo_validation is not None and pablo_validation["status"] != "success":
            raise click.ClickException(_pablo_failure_message(pablo_validation))
        colored_echo()
        click.echo(click.style("Dry run complete; no PDB file written.", fg="yellow"))
        return

    try:
        written = write_normalized_pdb(input_file, output_file, plan)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    pablo_validation = None
    if check_pablo:
        pablo_validation = _validate_clean_pdb_with_pablo(
            written,
            temporary_check_file=False,
        )
    _write_clean_pdb_report(report_path, plan, pablo_validation)
    if pablo_validation is not None and pablo_validation["status"] != "success":
        raise click.ClickException(_pablo_failure_message(pablo_validation))

    colored_echo()
    click.echo(click.style(f"Cleaned PDB written to: {written}", fg="green"))


def _write_clean_pdb_report(
    report_path: Path | None,
    plan: Any,
    pablo_validation: dict[str, Any] | None = None,
) -> None:
    """Write a clean-PDB JSON report when requested.

    Parameters
    ----------
    report_path : Path or None
        Destination report path, or ``None`` when reporting is disabled.
    plan : Any
        PDB normalization plan with Pydantic ``model_dump`` support.
    pablo_validation : dict[str, Any] or None, optional
        Optional Pablo validation object to append to the report, by default
        ``None``.
    """
    if report_path is None:
        return
    payload = plan.model_dump(mode="json")
    if pablo_validation is not None:
        payload["pablo_validation"] = pablo_validation
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(payload, indent=2) + "\n")
    colored_echo(f"  Report JSON: {report_path}")


def _temporary_clean_pdb_path(output_file: Path) -> Path:
    """Create a temporary path for dry-run Pablo validation.

    Parameters
    ----------
    output_file : Path
        User-facing output path whose parent is preferred for the temporary
        cleaned PDB.

    Returns
    -------
    Path
        Temporary PDB path that the caller must remove.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{output_file.stem}_pablo_",
        suffix=".pdb",
        dir=output_file.parent,
        delete=False,
    ) as handle:
        return Path(handle.name)


def _validate_clean_pdb_with_pablo(
    checked_path: Path,
    *,
    temporary_check_file: bool,
) -> dict[str, Any]:
    """Validate a normalized PDB with Pablo and return a JSON-safe report.

    Parameters
    ----------
    checked_path : Path
        Normalized PDB path to validate.
    temporary_check_file : bool
        Whether ``checked_path`` is a temporary dry-run file.

    Returns
    -------
    dict[str, Any]
        JSON-safe validation status and result diagnostics.
    """
    if not _pdb_has_explicit_hydrogens(checked_path):
        return _failed_pablo_validation(
            checked_path,
            temporary_check_file=temporary_check_file,
            error=(
                "OpenFF Pablo requires explicit hydrogens, but the normalized PDB contains no "
                "hydrogen atoms. Add hydrogens with a preparation backend before using "
                "clean-pdb --check-pablo."
            ),
        )

    from polyzymd.builders.conjugation.pablo_adapter import PabloIngestor
    from polyzymd.config.schema import (
        ConjugationCcdPabloPolicyConfig,
        ConjugationChainPolicyConfig,
    )

    try:
        result = PabloIngestor(ConjugationCcdPabloPolicyConfig()).ingest_structure(
            checked_path,
            chain_policy=ConjugationChainPolicyConfig(),
        )
    except Exception as exc:  # noqa: BLE001 - CLI must normalize third-party parser failures
        return _failed_pablo_validation(
            checked_path,
            temporary_check_file=temporary_check_file,
            error=str(exc),
        )

    payload: dict[str, Any] = {
        "attempted": True,
        "status": "success" if getattr(result, "success", False) else "failed",
        "checked_path": str(checked_path),
        "temporary_check_file": temporary_check_file,
        "result": _json_safe_pablo_payload(result),
    }
    if payload["status"] != "success":
        payload["diagnostics"] = _json_safe_pablo_payload(getattr(result, "diagnostics", []))
    return payload


def _pdb_has_explicit_hydrogens(path: Path) -> bool:
    """Return whether a PDB file contains explicit hydrogen atom records.

    Parameters
    ----------
    path : Path
        PDB file to inspect.

    Returns
    -------
    bool
        ``True`` when at least one ATOM/HETATM line has element H or a hydrogen
        atom name.
    """
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = line[76:78].strip().upper()
        atom_name = line[12:16].strip().upper()
        if element == "H" or atom_name.startswith("H"):
            return True
    return False


def _json_safe_pablo_payload(value: Any) -> Any:
    """Convert Pablo adapter objects to JSON-safe data.

    Parameters
    ----------
    value : Any
        Pydantic model, sequence, mapping, or scalar value to serialize.

    Returns
    -------
    Any
        JSON-compatible representation.
    """
    if hasattr(value, "model_dump"):
        return value.model_dump(mode="json")
    if isinstance(value, dict):
        return {str(key): _json_safe_pablo_payload(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe_pablo_payload(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    try:
        json.dumps(value)
    except TypeError:
        return str(value)
    return value


def _failed_pablo_validation(
    checked_path: Path,
    *,
    temporary_check_file: bool,
    error: str,
) -> dict[str, Any]:
    """Build a failed Pablo validation report.

    Parameters
    ----------
    checked_path : Path
        Path that was checked or intended for checking.
    temporary_check_file : bool
        Whether the path is a temporary dry-run artifact.
    error : str
        Normalized parser or import error message.

    Returns
    -------
    dict[str, Any]
        JSON-safe failed validation object.
    """
    return {
        "attempted": True,
        "status": "failed",
        "checked_path": str(checked_path),
        "temporary_check_file": temporary_check_file,
        "diagnostics": [{"severity": "error", "message": error}],
    }


def _skipped_pablo_validation(reason: str) -> dict[str, Any]:
    """Build a skipped Pablo validation report.

    Parameters
    ----------
    reason : str
        Reason Pablo validation did not run.

    Returns
    -------
    dict[str, Any]
        JSON-safe skipped validation object.
    """
    return {
        "attempted": False,
        "status": "skipped",
        "checked_path": None,
        "temporary_check_file": False,
        "diagnostics": [{"severity": "info", "message": reason}],
    }


def _pablo_failure_message(pablo_validation: dict[str, Any]) -> str:
    """Return a concise Click error message for Pablo validation failures.

    Parameters
    ----------
    pablo_validation : dict[str, Any]
        Pablo validation report produced by clean-PDB.

    Returns
    -------
    str
        User-facing failure message.
    """
    diagnostics = pablo_validation.get("diagnostics") or []
    if diagnostics:
        first = diagnostics[0]
        if isinstance(first, dict) and first.get("message"):
            return f"Pablo validation failed: {first['message']}"
    return "Pablo validation failed; review the report JSON or parser diagnostics."


def _echo_clean_pdb_plan_summary(plan: Any, output_file: Path, dry_run: bool) -> None:
    """Print a concise clean-PDB normalization plan summary.

    Parameters
    ----------
    plan
        Normalization plan returned by the pure-Python PDB planner.
    output_file : Path
        Destination path that would be used for a valid non-dry run.
    dry_run : bool
        Whether the CLI is planning only.
    """
    status = "valid" if plan.valid else "invalid"
    colored_echo(f"  Validation: {status}")
    colored_echo(f"  Protein residues: {plan.protein_residue_count}")
    colored_echo(f"  Attached moiety residues: {plan.moiety_residue_count}")
    colored_echo(f"  Chain ID changes: {plan.chain_id_change_count}")
    colored_echo(f"  Residue renumbering changes: {plan.residue_number_change_count}")
    if plan.source_chains_collapsed:
        collapsed = "; ".join(
            f"{target} <= {', '.join(sources)}"
            for target, sources in plan.source_chains_collapsed.items()
        )
        colored_echo(f"  Chain collapse: {collapsed}")
    if plan.warnings:
        colored_echo("  Warnings:")
        for warning in plan.warnings[:5]:
            colored_echo(f"    - {warning}")
        if len(plan.warnings) > 5:
            colored_echo(f"    - ... {len(plan.warnings) - 5} more warning(s)")

    output_label = "Would write" if dry_run else "Output path"
    colored_echo(f"  {output_label}: {output_file}")
    _echo_clean_pdb_action_examples(plan.actions)


def _echo_clean_pdb_action_examples(actions: list[Any]) -> None:
    """Print the first few residue normalization action mappings.

    Parameters
    ----------
    actions : list
        Residue-level action mappings from the normalization plan.
    """
    if not actions:
        return
    colored_echo("  First action mappings:")
    for action in actions[:5]:
        source_chain = action.source_chain or "blank"
        source_number = action.res_seq or "?"
        source_icode = action.insertion_code or ""
        target_icode = action.target_insertion_code or ""
        colored_echo(
            "    - "
            f"{source_chain}:{action.residue_name}{source_number}{source_icode} -> "
            f"{action.target_chain}:{action.residue_name}{action.target_res_seq}{target_icode}"
        )
    if len(actions) > 5:
        colored_echo(f"    - ... {len(actions) - 5} more mapping(s)")


def _echo_clean_pdb_issues(issues: list[Any]) -> None:
    """Print clean-PDB validation issues.

    Parameters
    ----------
    issues : list
        Validation issues from the normalization plan.
    """
    if not issues:
        return
    colored_echo("  Issues:")
    for issue in issues[:8]:
        residue = f" ({issue.residue_id})" if issue.residue_id else ""
        colored_echo(f"    - [{issue.code}]{residue} {issue.message}")
    if len(issues) > 8:
        colored_echo(f"    - ... {len(issues) - 8} more issue(s)")


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
    "--email",
    default="",
    help="Email for job notifications",
)
@click.option(
    "--memory",
    default=None,
    help="Override SLURM memory allocation (e.g. '4G', '8G'). Not needed for bridges2 (allocated per GPU).",
)
@click.option(
    "--partition",
    default=None,
    help="Override SLURM partition",
)
@click.option(
    "--qos",
    default=None,
    help="Override SLURM QoS",
)
@click.option(
    "--constraint",
    default=None,
    help=(
        "SLURM --constraint expression for node feature selection. "
        "Supports boolean expressions: 'A40' (single), 'A40|A100' (OR)."
    ),
)
@click.option(
    "--nodelist",
    default=None,
    help="SLURM --nodelist override (e.g., 'node01' or 'node[01-04]').",
)
@click.option(
    "--pixi-env",
    default=None,
    type=click.Choice(["cuda-12-4", "cuda-12-6"]),
    help=("Pixi environment for the recovery SLURM job. If omitted, inferred from --preset."),
)
@click.option(
    "--force",
    is_flag=True,
    help="Skip duplicate-job check and submit even if a SLURM job is already running for the replicate",
)
@click.option(
    "--engine",
    default=None,
    help="Override simulation engine",
)
def recover(
    config: str,
    replicate: int,
    scratch_dir: str | None,
    preset: str,
    submit: bool,
    dry_run: bool,
    email: str,
    memory: str | None,
    partition: str | None,
    qos: str | None,
    constraint: str | None,
    nodelist: str | None,
    pixi_env: str | None,
    force: bool,
    engine: str | None,
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
    from polyzymd.engines import create_engine
    from polyzymd.simulation.progress import save_progress

    _echo_branding()

    try:
        sim_config = SimulationConfig.from_yaml(config)
    except (FileNotFoundError, yaml.YAMLError, ValidationError, ValueError) as e:
        colored_echo(f"Failed to load config: {e}", err=True, phase="workflow", level=logging.ERROR)
        sys.exit(1)

    engine_name = _resolve_engine_name(sim_config, override=engine)
    engine_impl = create_engine(sim_config, override=engine_name, defer_binary=True)

    if scratch_dir:
        working_dir = Path(scratch_dir)
    else:
        working_dir = engine_impl.get_engine_working_directory(sim_config, replicate)

    if not working_dir.exists():
        colored_echo(
            f"Working directory not found: {working_dir}",
            err=True,
            phase="workflow",
            level=logging.ERROR,
        )
        sys.exit(1)

    # Calculate progress metadata from config
    prod = sim_config.simulation_phases.production
    timestep_fs = prod.time_step

    # Load progress
    progress = engine_impl.load_or_scan_progress(working_dir, replicate)
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
            f"  polyzymd recover -c {config} -r {replicate} --submit --preset <preset>"
            f" [engine={engine_name}]",
            phase="workflow",
        )
        return

    # Generate and submit a self-resubmitting SLURM job
    from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs, create_job_name
    from polyzymd.workflow.slurm import PRESET_DEFAULT_PIXI_ENV, SlurmConfig, SlurmScriptGenerator

    # Resolve pixi environment: explicit flag > preset default
    resolved_pixi_env = pixi_env or PRESET_DEFAULT_PIXI_ENV.get(preset, "cuda-12-4")

    colored_echo(
        f"\nGenerating recovery job (preset: {preset}, engine: {engine_name}, "
        f"pixi env: {resolved_pixi_env})...",
        phase="workflow",
    )

    # Use the same descriptive job naming as `polyzymd submit`
    job_name = create_job_name(sim_config, replicate)

    # Best-effort duplicate guard
    if not force:
        existing = check_existing_slurm_jobs(job_name)
        if existing:
            ids = ", ".join(existing)
            colored_echo(
                f"Replicate {replicate} already has RUNNING/PENDING SLURM "
                f"job(s): {ids} (job name '{job_name}'). "
                "Use --force to submit anyway.",
                err=True,
                phase="workflow",
                level=logging.ERROR,
            )
            sys.exit(1)

    if engine_name == "gromacs":
        import shutil

        from polyzymd.engines import create_engine
        from polyzymd.engines.base import EngineSubmitRequest

        slurm_config = SlurmConfig.from_preset(preset)
        if email:
            slurm_config.email = email
        if memory:
            slurm_config.memory = memory
        if partition:
            slurm_config.partition = partition
        if qos:
            slurm_config.qos = qos
        if constraint:
            slurm_config.constraint = constraint
        if nodelist:
            slurm_config.nodelist = nodelist

        prefix = _generate_system_prefix(sim_config)
        gromacs_inputs_exist = all(
            (working_dir / f).exists()
            for f in [f"{prefix}.top", f"{prefix}.gro", "em.mdp", "prod.mdp"]
        )
        if gromacs_inputs_exist:
            colored_echo(
                "Detected existing GROMACS inputs — recovery will reuse them",
                phase="workflow",
            )
            if not list(working_dir.glob("eq_*.mdp")):
                colored_echo(
                    "Warning: core GROMACS inputs found but no equilibration MDPs "
                    "(eq_*.mdp). The recovery job will skip equilibration and go "
                    "straight to production.",
                    phase="workflow",
                    level=logging.WARNING,
                )

        config_path_abs = str(Path(config).resolve())
        request = EngineSubmitRequest(
            replicate=replicate,
            config_path=Path(config_path_abs),
            working_dir=working_dir,
            slurm_config=slurm_config,
            job_name=job_name,
            extra={"pixi_env": resolved_pixi_env, "skip_build": gromacs_inputs_exist},
        )

        engine_impl = create_engine(sim_config, override="gromacs", defer_binary=True)
        script_path = engine_impl.prepare_submission(request)

        recovery_dir = working_dir / "recovery_scripts"
        recovery_dir.mkdir(exist_ok=True)
        recovery_path = recovery_dir / f"recover_rep{replicate}.sh"
        shutil.copy2(script_path, recovery_path)
        recovery_path.chmod(0o755)

        colored_echo(f"Script: {recovery_path}", phase="workflow")

        if dry_run:
            colored_echo("\n[DRY RUN] Would submit:", phase="workflow")
            colored_echo(f"  sbatch {recovery_path}", phase="workflow")
            return

        from polyzymd.workflow.slurm_submit import run_sbatch

        module_load = (
            getattr(sim_config.gromacs, "module_load", None) if sim_config.gromacs else None
        )
        result = run_sbatch(recovery_path, module_load=module_load)
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
    else:
        # OpenMM recovery path
        slurm_config = SlurmConfig.from_preset(preset)
        if email:
            slurm_config.email = email
        if memory:
            slurm_config.memory = memory
        if partition:
            slurm_config.partition = partition
        if qos:
            slurm_config.qos = qos
        if constraint:
            slurm_config.constraint = constraint
        if nodelist:
            slurm_config.nodelist = nodelist

        # Detect pre-built system files so the recovery job loads the existing
        # topology instead of rebuilding (non-deterministic packing would produce
        # a different atom count and crash on checkpoint reload).
        system_already_built = (working_dir / "solvated_system.pdb").exists() and (
            working_dir / "system.xml"
        ).exists()
        if system_already_built:
            colored_echo(
                "Detected pre-built system — recovery job will use --skip-build",
                phase="workflow",
            )

        generator = SlurmScriptGenerator(
            slurm_config, pixi_env=resolved_pixi_env, skip_build=system_already_built
        )

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

    _echo_branding()
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
    colored_echo("Example configs: polyzymd/templates/examples/", phase="cli")


# =============================================================================
# Scaffold Command
# =============================================================================


@cli.command("new-analysis")
@click.argument("name")
@click.option(
    "--class-name",
    default=None,
    help="PascalCase class prefix (default: auto-derived from NAME).",
)
@click.option(
    "--style",
    type=click.Choice(["dict", "pydantic"], case_sensitive=False),
    default="dict",
    help="Template style: 'dict' (default) for plain dicts, 'pydantic' for typed result models.",
)
@click.option(
    "--project-root",
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
    default=None,
    help="Repository root. Default: auto-detected from this file's location.",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Overwrite existing files.",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help="Print what would be created without writing files.",
)
def new_analysis(
    name: str,
    class_name: str | None,
    style: str,
    project_root: str | None,
    force: bool,
    dry_run: bool,
) -> None:
    """Scaffold a new analysis plugin.

    NAME is the snake_case plugin name (e.g. 'solvent_shell').

    Creates:

    \b
      src/polyzymd/analyses/<NAME>/__init__.py    — plugin class
      tests/analyses/plugins/test_<NAME>.py       — smoke tests

    Styles:

    \b
      dict      Plain dicts for results (default, simplest)
      pydantic  Typed Pydantic result models (better for complex analyses)

    Run the generated tests with:

    \b
      pixi run -e build pytest tests/analyses/plugins/test_<NAME>.py -v
    """
    from polyzymd.cli.scaffold import generate_scaffold, validate_class_name, validate_name

    # Validate name
    error = validate_name(name)
    if error:
        raise click.BadParameter(error, param_hint="'NAME'")

    # Validate class name if provided
    if class_name is not None:
        cls_error = validate_class_name(class_name)
        if cls_error:
            raise click.BadParameter(cls_error, param_hint="'--class-name'")

    # Resolve project root
    if project_root is None:
        # Walk up from this file to find pyproject.toml
        root = Path(__file__).resolve().parent
        for _ in range(10):
            if (root / "pyproject.toml").exists():
                break
            root = root.parent
        else:
            raise click.UsageError(
                "Could not auto-detect project root. Pass --project-root explicitly."
            )
    else:
        root = Path(project_root)

    try:
        created = generate_scaffold(
            name,
            root,
            class_name=class_name,
            style=style,
            force=force,
            dry_run=dry_run,
        )
    except FileExistsError as exc:
        raise click.ClickException(str(exc)) from exc

    verb = "Would create" if dry_run else "Created"
    for p in created:
        try:
            display = p.relative_to(root)
        except ValueError:
            display = p
        colored_echo(f"  {verb}: {display}", phase="cli")

    if not dry_run:
        colored_echo(f"\nPlugin '{name}' scaffolded successfully!", phase="cli")
        colored_echo(
            f"Run tests: pixi run -e build pytest tests/analyses/plugins/test_{name}.py -v",
            phase="cli",
        )


# =============================================================================
# Analysis Commands (from analysis module)
# =============================================================================


def _register_optional_command_groups() -> None:
    """Register optional command groups when deps are importable."""
    from polyzymd.cli.compare import compare

    cli.add_command(compare)


# =============================================================================
# Hidden: GROMACS progress update (called by GROMACS SLURM scripts)
# =============================================================================


@cli.command("_update-gromacs-progress", hidden=True)
@click.option(
    "--working-dir",
    required=True,
    type=click.Path(exists=True),
    help="GROMACS working directory",
)
@click.option(
    "--config-path",
    default="",
    help="Simulation config path for metadata",
)
@click.option(
    "--replicate",
    default=1,
    type=int,
    help="Replicate index",
)
@click.option(
    "--mark-complete",
    is_flag=True,
    help="Force completion status after post-processing",
)
def update_gromacs_progress_cmd(
    working_dir: str,
    config_path: str,
    replicate: int,
    mark_complete: bool,
) -> None:
    """Update GROMACS progress.json from prod.log scan.

    Called by the GROMACS self-resubmitting SLURM wrapper after each
    mdrun invocation. Merges the latest log scan into progress.json.
    """
    from polyzymd.engines.gromacs.progress import update_gromacs_progress

    progress = update_gromacs_progress(
        working_dir=Path(working_dir),
        config_path=config_path,
        replicate=replicate,
        mark_complete=mark_complete,
    )
    pct = progress.fraction_complete() * 100
    click.echo(
        f"Progress updated: {progress.total_steps_completed}/{progress.total_steps_requested} "
        f"steps ({pct:.1f}%) — status: {progress.status.value}"
    )


_register_optional_command_groups()


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
