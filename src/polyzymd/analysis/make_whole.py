"""
Make trajectories visually whole for visualization.

This module provides tools to post-process MD trajectories by:
1. Unwrapping molecules split across periodic boundaries
2. Centering protein and substrate in the simulation box

The output trajectories are suitable for visualization in PyMOL, VMD, etc.
without visual artifacts from periodic boundary conditions.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional, Tuple, Union

import MDAnalysis as mda
from MDAnalysis import transformations

if TYPE_CHECKING:
    from polyzymd.config.schema import SimulationConfig

LOGGER = logging.getLogger(__name__)


def find_production_trajectories(
    working_dir: Union[str, Path],
    segment_range: Optional[Tuple[int, int]] = None,
) -> Tuple[Path, List[Path]]:
    """
    Find topology and trajectory files in a simulation directory.

    Automatically discovers all production_N directories and their DCD files,
    sorted in correct numerical order.

    Args:
        working_dir: Path to simulation working directory (in scratch)
        segment_range: Optional (start, end) to load only segments start through end
                       (inclusive). E.g., (3, 7) loads production_3 through production_7.

    Returns:
        Tuple of (topology_path, [trajectory_paths...])

    Raises:
        FileNotFoundError: If no topology or trajectories found
    """
    working_dir = Path(working_dir)

    if not working_dir.exists():
        raise FileNotFoundError(f"Working directory not found: {working_dir}")

    # Find all production_N directories
    production_dirs = list(working_dir.glob("production_*"))

    if not production_dirs:
        raise FileNotFoundError(f"No production_N directories found in {working_dir}")

    # Sort directories numerically
    def extract_dir_number(dirpath: Path) -> int:
        match = re.search(r"production_(\d+)$", dirpath.name)
        return int(match.group(1)) if match else 0

    production_dirs.sort(key=extract_dir_number)

    # Filter by segment range if specified
    if segment_range is not None:
        start, end = segment_range
        production_dirs = [d for d in production_dirs if start <= extract_dir_number(d) <= end]
        if not production_dirs:
            raise FileNotFoundError(f"No production directories found in range {start}-{end}")

    LOGGER.info(f"Found {len(production_dirs)} production directories")

    # Find topology file (from first production directory)
    topology_file: Optional[Path] = None
    for prod_dir in production_dirs:
        # Try production_N_topology.pdb pattern
        candidates = list(prod_dir.glob("*_topology.pdb"))
        if candidates:
            topology_file = candidates[0]
            break

    if topology_file is None:
        # Try solvated_system.pdb in working directory
        solvated = working_dir / "solvated_system.pdb"
        if solvated.exists():
            topology_file = solvated
        else:
            raise FileNotFoundError(
                f"No topology file found in {working_dir}. "
                "Expected *_topology.pdb in production_N/ or solvated_system.pdb"
            )

    LOGGER.info(f"Using topology: {topology_file}")

    # Collect all DCD files from production directories
    all_trajectories: List[Path] = []

    for prod_dir in production_dirs:
        # Find trajectory files (production_N_trajectory.dcd pattern)
        dcd_files = list(prod_dir.glob("*_trajectory.dcd"))

        # Sort by number in filename
        def extract_traj_number(filepath: Path) -> int:
            match = re.search(r"production_(\d+)_trajectory\.dcd$", filepath.name)
            return int(match.group(1)) if match else 0

        dcd_files.sort(key=extract_traj_number)
        all_trajectories.extend(dcd_files)

    if not all_trajectories:
        raise FileNotFoundError(
            f"No trajectory files (*_trajectory.dcd) found in production directories"
        )

    LOGGER.info(f"Found {len(all_trajectories)} trajectory files")

    return topology_file, all_trajectories


def make_whole_trajectory(
    config: Union[str, Path, "SimulationConfig"],
    replicate: int = 1,
    output: Optional[Union[str, Path]] = None,
    segment_range: Optional[Tuple[int, int]] = None,
    strip_solvent: bool = False,
    strip_cosolvent: bool = False,
    water_resnames: str = "HOH SOL WAT TIP3",
) -> Path:
    """
    Make trajectory visually whole by unwrapping PBC and centering protein+ligand.

    Applies MDAnalysis transformations to:
    1. Unwrap molecules split across periodic boundaries
    2. Center protein and substrate in the simulation box

    The ligand residue name is extracted from config.substrate.residue_name
    (single source of truth from config.yaml).

    Output is written to the scratch directory derived from config.output.

    Args:
        config: Path to config.yaml or SimulationConfig object
        replicate: Replicate number (used to determine working directory)
        output: Output trajectory path. If None, uses default naming:
                - Full trajectory: merged_madeWhole_production.dcd
                - Partial: merged_madeWhole_production(start-end).dcd
        segment_range: Optional (start, end) to process only certain segments (inclusive).
                       E.g., (3, 7) loads production_3 through production_7.
        strip_solvent: If True, exclude water from output (smaller file, no water analysis)
        strip_cosolvent: If True, exclude co-solvents from output
        water_resnames: Space-separated water residue names for stripping
                        (default: "HOH SOL WAT TIP3")

    Returns:
        Path to output trajectory file

    Raises:
        FileNotFoundError: If config, working directory, or trajectories not found

    Example:
        >>> # Process all segments for replicate 1
        >>> make_whole_trajectory("config.yaml", replicate=1)

        >>> # Process only segments 3-7
        >>> make_whole_trajectory("config.yaml", replicate=1, segment_range=(3, 7))

        >>> # Strip solvent for smaller visualization file
        >>> make_whole_trajectory("config.yaml", replicate=1, strip_solvent=True)
    """
    from polyzymd.config.schema import SimulationConfig

    # Load config if path provided
    if isinstance(config, (str, Path)):
        config_path = Path(config)
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        config = SimulationConfig.from_yaml(config_path)

    # Get working directory from config (single source of truth)
    working_dir = config.get_working_directory(replicate)
    LOGGER.info(f"Working directory: {working_dir}")

    # Determine output path
    if output is None:
        if segment_range:
            filename = f"merged_madeWhole_production({segment_range[0]}-{segment_range[1]}).dcd"
        else:
            filename = "merged_madeWhole_production.dcd"
        output = working_dir / filename
    else:
        output = Path(output)

    LOGGER.info(f"Output will be written to: {output}")

    # Find topology and trajectory files
    topology_file, trajectory_files = find_production_trajectories(working_dir, segment_range)

    # Build center selection from config
    # Protein is always included; add ligand if present
    ligand_resname = config.substrate.residue_name if config.substrate else None
    if ligand_resname:
        center_selection = f"protein or resname {ligand_resname}"
        LOGGER.info(f"Centering on: protein and {ligand_resname}")
    else:
        center_selection = "protein"
        LOGGER.info("Centering on: protein only (no substrate)")

    # Load universe with all trajectories
    LOGGER.info(f"Loading {len(trajectory_files)} trajectory files...")
    trajectory_strs = [str(t) for t in trajectory_files]
    u = mda.Universe(str(topology_file), *trajectory_strs)

    LOGGER.info(f"Loaded universe: {u.atoms.n_atoms} atoms, {len(u.trajectory)} frames")

    # Build transformation pipeline
    # 1. Unwrap - make molecules whole across PBC
    unwrap_transform = transformations.unwrap(u.atoms)

    # 2. Center protein+ligand in box
    center_ag = u.select_atoms(center_selection)
    if len(center_ag) == 0:
        raise ValueError(
            f"Center selection '{center_selection}' matched no atoms. "
            "Check that the topology contains protein atoms."
        )
    center_transform = transformations.center_in_box(center_ag, center="mass")

    # Apply transformations in order
    u.trajectory.add_transformations(unwrap_transform, center_transform)

    LOGGER.info("Applied transformations: unwrap, center_in_box")

    # Determine which atoms to write
    if strip_solvent or strip_cosolvent:
        exclude_parts = []

        if strip_solvent:
            # Build exclusion for water
            water_sel = " ".join(f"resname {r}" for r in water_resnames.split())
            exclude_parts.append(f"({water_sel})")
            LOGGER.info(f"Stripping solvent: {water_resnames}")

        if strip_cosolvent:
            # Get co-solvent resnames from config
            for cs in config.solvent.co_solvents:
                exclude_parts.append(f"resname {cs.residue_name}")
                LOGGER.info(f"Stripping co-solvent: {cs.residue_name}")

        exclusion = " or ".join(exclude_parts)
        selection = f"not ({exclusion})"
        output_atoms = u.select_atoms(selection)
        LOGGER.info(f"Writing {output_atoms.n_atoms} atoms (after stripping)")
    else:
        output_atoms = u.atoms
        LOGGER.info(f"Writing all {output_atoms.n_atoms} atoms")

    # Ensure output directory exists
    output.parent.mkdir(parents=True, exist_ok=True)

    # Write transformed trajectory
    LOGGER.info(f"Writing trajectory to {output}...")
    with mda.Writer(str(output), output_atoms.n_atoms) as writer:
        for ts in u.trajectory:
            writer.write(output_atoms)

    LOGGER.info(f"Successfully wrote {len(u.trajectory)} frames to {output}")

    return output
