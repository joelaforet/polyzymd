#!/usr/bin/env python
"""Convert pre-PolyzyMD legacy simulation directories for PolyzyMD analysis.

This script takes legacy simulation directories (pre-PolyzyMD format) and generates
the files needed for PolyzyMD's current ``compare`` workflow:

1. ``solvated_system.pdb`` — topology with corrected chain assignments
   (A=protein, B=substrate, C=polymer, D=solvent/ions)
2. ``config.yaml`` — valid PolyzyMD SimulationConfig for TrajectoryLoader

The script creates an output directory that mirrors the legacy structure but uses
symlinks for the (multi-GB) trajectory data, so no data duplication occurs.

Dependencies
------------
- Python 3.10+
- MDAnalysis
- PyYAML
- biopandas (for protein atom name fixing)

All available in the PolyzyMD ``build`` pixi environment.

Usage
-----
::

    # Single simulation
    pixi run -e build python convert_legacy.py \\
        /path/to/10A_RESTRAINT_LipA_..._run1 \\
        --output-dir /path/to/polyzymd_converted \\
        --reference-pdb /path/to/NH3_terminal_His_proton_updated.pdb

    # Batch: all simulations in a parent directory
    pixi run -e build python convert_legacy.py \\
        /path/to/PREPOLYZYMD_SIMULATIONS \\
        --output-dir /path/to/polyzymd_converted \\
        --reference-pdb /path/to/NH3_terminal_His_proton_updated.pdb \\
        --batch

    # Dry-run (preview only)
    pixi run -e build python convert_legacy.py \\
        /path/to/PREPOLYZYMD_SIMULATIONS \\
        --output-dir /path/to/polyzymd_converted \\
        --reference-pdb /path/to/NH3_terminal_His_proton_updated.pdb \\
        --batch --dry-run

Legacy Naming Convention
------------------------
``<restraint>_<Enzyme>_<Substrate>_[conf<N>_][<Polymer>_<count>x_]<Temp>_<Equil>-<Ens>_<Prod>-<Ens>_<Run>``

Examples::

    10A_RESTRAINT_LipA_Resorufin-Butyrate_363.0K_0.5ns-NVT_1000.0ns-NPT_run1
    10A_RESTRAINT_LipA_Resorufin-Butyrate_SBMA-EGMA-50%_38x_363.0K_0.5ns-NVT_1000.0ns-NPT_run1
    10A_RESTRAINT_CALB_Resorufin-Butyrate_conf2_300.0K_0.5ns-NVT_1000.0ns-NPT_run1
    10A_RESTRAINT_CALB_Resorufin-Butyrate_conf2_SBMA-EGMA-50%_77x_363.0K_0.5ns-NVT_1000.0ns-NPT_run1
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)
logger = logging.getLogger(__name__)

# Known polymer residue names for chain assignment
KNOWN_POLYMER_RESNAMES = {"SBM", "EGM", "EGP"}

# Known solvent / ion residue names
KNOWN_SOLVENT_RESNAMES = {"HOH", "WAT", "TIP3", "Na+", "Cl-", "NA", "CL", "K+", "SOL"}

# Polymer name mapping: folder substring -> (type_prefix, monomer_A_name, monomer_B_name)
POLYMER_NAME_MAP = {
    "SBMA-EGMA": ("SBMA-EGMA", "SBMA", "EGMA"),
    "SBMA-EGPMA": ("SBMA-EGPMA", "SBMA", "EGPMA"),
}


# =============================================================================
# Data classes for parsed simulation metadata
# =============================================================================


@dataclass
class PolymerInfo:
    """Parsed polymer information from folder name."""

    type_prefix: str  # e.g., "SBMA-EGMA"
    monomer_a_name: str  # e.g., "SBMA"
    monomer_b_name: str  # e.g., "EGPMA"
    percent_b: float  # e.g., 50.0 (percentage of second component)
    chain_count: int  # e.g., 38


@dataclass
class PhaseParams:
    """Simulation phase parameters extracted from *_parameters.json."""

    ensemble: str  # "NVT" or "NPT"
    duration_ns: float  # total duration in nanoseconds
    samples: int  # number of trajectory frames
    time_step_fs: float  # time step in femtoseconds
    temperature_K: float
    pressure_atm: float
    thermostat: str = "LangevinMiddle"
    barostat: str = "MC"
    barostat_frequency: int = 25


@dataclass
class SimMetadata:
    """Complete parsed simulation metadata."""

    folder_name: str
    restraint_distance: str  # e.g., "10A"
    enzyme_name: str  # e.g., "LipA"
    substrate_name: str  # e.g., "Resorufin-Butyrate"
    temperature_K: float  # e.g., 363.0
    replicate: int  # e.g., 1
    conformation: int | None = None  # e.g., 2 (from "conf2" in folder name)
    polymer: PolymerInfo | None = None
    equilibration: PhaseParams | None = None
    production: PhaseParams | None = None
    n_segments: int = 1
    topology_path: Path | None = None
    trajectory_paths: list[Path] = field(default_factory=list)
    production_dirs: list[Path] = field(default_factory=list)


# =============================================================================
# Folder name parsing
# =============================================================================


def parse_folder_name(folder_name: str) -> SimMetadata:
    """Parse a legacy simulation folder name into structured metadata.

    Parameters
    ----------
    folder_name : str
        The directory name (not full path), e.g.,
        ``10A_RESTRAINT_LipA_Resorufin-Butyrate_SBMA-EGMA-50%_38x_363.0K_0.5ns-NVT_1000.0ns-NPT_run1``
        ``10A_RESTRAINT_CALB_Resorufin-Butyrate_conf2_SBMA-EGMA-50%_77x_363.0K_0.5ns-NVT_1000.0ns-NPT_run1``

    Returns
    -------
    SimMetadata
        Parsed simulation metadata.

    Raises
    ------
    ValueError
        If the folder name does not match the expected pattern.
    """
    # Pattern for the full folder name with optional conformation and polymer fields
    # Group structure:
    #   1: restraint distance (e.g., "10A")
    #   2: enzyme name (e.g., "LipA", "CALB")
    #   3: substrate name (e.g., "Resorufin-Butyrate")
    #   4: optional conformation number (e.g., "2" from "conf2")
    #   5: optional polymer composition (e.g., "SBMA-EGMA-50%")
    #   6: optional polymer chain count (e.g., "38", "77")
    #   7: temperature (e.g., "363.0")
    #   8: equil duration (e.g., "0.5")
    #   9: equil ensemble (e.g., "NVT")
    #  10: prod duration (e.g., "1000.0")
    #  11: prod ensemble (e.g., "NPT")
    #  12: replicate (e.g., "1")
    pattern = re.compile(
        r"^(\w+)_RESTRAINT_"  # restraint distance
        r"(\w+)_"  # enzyme name
        r"([\w-]+?)"  # substrate name (non-greedy)
        r"(?:_conf(\d+))?"  # optional: conformation number
        r"(?:_([\w-]+?-\d+%)_(\d+)x)?"  # optional: polymer_comp_count
        r"_([\d.]+)K_"  # temperature
        r"([\d.]+)ns-(\w+)_"  # equil duration-ensemble
        r"([\d.]+)ns-(\w+)_"  # prod duration-ensemble
        r"run(\d+)$"  # replicate
    )

    match = pattern.match(folder_name)
    if not match:
        raise ValueError(
            f"Folder name does not match expected legacy naming convention:\n"
            f"  {folder_name}\n"
            f"Expected pattern: <restraint>_RESTRAINT_<Enzyme>_<Substrate>"
            f"[_conf<N>][_<Polymer>_<count>x]_<Temp>K_<Equil>ns-<Ens>_<Prod>ns-<Ens>_run<N>"
        )

    restraint = match.group(1)
    enzyme = match.group(2)
    substrate = match.group(3)
    conf_str = match.group(4)  # None if no conformation tag
    polymer_str = match.group(5)  # None if no polymer
    polymer_count_str = match.group(6)  # None if no polymer
    temp = float(match.group(7))
    equil_dur = float(match.group(8))
    equil_ens = match.group(9)
    prod_dur = float(match.group(10))
    prod_ens = match.group(11)
    replicate = int(match.group(12))

    # Parse polymer info if present
    polymer = None
    if polymer_str is not None:
        polymer = _parse_polymer_string(polymer_str, int(polymer_count_str))

    metadata = SimMetadata(
        folder_name=folder_name,
        restraint_distance=restraint,
        enzyme_name=enzyme,
        substrate_name=substrate,
        temperature_K=temp,
        replicate=replicate,
        conformation=int(conf_str) if conf_str else None,
        polymer=polymer,
    )

    # Store folder-name-derived phase info (will be overwritten by JSON if available)
    metadata.equilibration = PhaseParams(
        ensemble=equil_ens,
        duration_ns=equil_dur,
        samples=10,  # placeholder, updated from JSON
        time_step_fs=2.0,
        temperature_K=temp,
        pressure_atm=1.0,
    )
    metadata.production = PhaseParams(
        ensemble=prod_ens,
        duration_ns=prod_dur,
        samples=2500,  # placeholder, updated from JSON
        time_step_fs=2.0,
        temperature_K=temp,
        pressure_atm=1.0,
    )

    return metadata


def _parse_polymer_string(polymer_str: str, chain_count: int) -> PolymerInfo:
    """Parse polymer composition string like 'SBMA-EGMA-50%' or 'SBMA-EGPMA-5%'.

    Parameters
    ----------
    polymer_str : str
        Polymer composition string from folder name.
    chain_count : int
        Number of polymer chains (e.g., 38).

    Returns
    -------
    PolymerInfo
        Parsed polymer information.
    """
    # Pattern: <typeA>-<typeB>-<percent>%
    match = re.match(r"^([\w]+)-([\w]+)-(\d+)%$", polymer_str)
    if not match:
        raise ValueError(f"Cannot parse polymer string: {polymer_str}")

    monomer_a = match.group(1)  # e.g., "SBMA"
    monomer_b = match.group(2)  # e.g., "EGMA" or "EGPMA"
    percent_b = float(match.group(3))

    type_prefix = f"{monomer_a}-{monomer_b}"

    return PolymerInfo(
        type_prefix=type_prefix,
        monomer_a_name=monomer_a,
        monomer_b_name=monomer_b,
        percent_b=percent_b,
        chain_count=chain_count,
    )


# =============================================================================
# File discovery
# =============================================================================


def discover_files(sim_dir: Path, metadata: SimMetadata) -> None:
    """Discover topology, trajectory, and production directories in a legacy sim folder.

    Updates ``metadata`` in-place with discovered file paths.

    Parameters
    ----------
    sim_dir : Path
        Path to the legacy simulation directory.
    metadata : SimMetadata
        Metadata object to update with discovered paths.
    """
    # --- Discover production directories ---
    # Check for daisy-chain segments: production_0/, production_1/, ...
    segment_dirs = sorted(
        [d for d in sim_dir.iterdir() if d.is_dir() and re.match(r"production_\d+$", d.name)],
        key=lambda d: int(d.name.split("_")[1]),
    )

    if segment_dirs:
        metadata.production_dirs = segment_dirs
        metadata.n_segments = len(segment_dirs)
    else:
        # Single production directory
        prod_dir = sim_dir / "production"
        if prod_dir.is_dir():
            metadata.production_dirs = [prod_dir]
            metadata.n_segments = 1
        else:
            raise FileNotFoundError(f"No production directory found in {sim_dir}")

    # --- Discover topology ---
    # Priority order matching TrajectoryLoader._find_topology()
    candidates = []

    # Check for updated_topology.pdb (cleanest chain assignments)
    for prod_dir in metadata.production_dirs:
        updated = prod_dir / "updated_topology.pdb"
        if updated.exists():
            candidates.append(("updated_topology", updated))

    # Check daisy-chain first segment topology
    if segment_dirs:
        seg0_topo = segment_dirs[0] / f"{segment_dirs[0].name}_topology.pdb"
        if seg0_topo.exists():
            candidates.append(("segment_0_topology", seg0_topo))

    # Check single production topology
    single_topo = sim_dir / "production" / "production_topology.pdb"
    if single_topo.exists():
        candidates.append(("production_topology", single_topo))

    if not candidates:
        raise FileNotFoundError(f"No topology PDB found in {sim_dir}")

    chosen_label, chosen_path = candidates[0]
    metadata.topology_path = chosen_path
    logger.info(f"  Topology: {chosen_label} -> {chosen_path.name}")

    # --- Discover trajectories ---
    trajectories = []

    if segment_dirs:
        # Daisy-chain: production_N/production_N_trajectory.dcd
        for seg_dir in segment_dirs:
            traj = seg_dir / f"{seg_dir.name}_trajectory.dcd"
            if traj.exists():
                trajectories.append(traj)
            else:
                logger.warning(f"  Missing trajectory: {traj}")
    else:
        # Single production: production/production_trajectory.dcd
        traj = sim_dir / "production" / "production_trajectory.dcd"
        if traj.exists():
            trajectories.append(traj)
        else:
            raise FileNotFoundError(f"No production trajectory found in {sim_dir}")

    metadata.trajectory_paths = trajectories
    logger.info(f"  Trajectories: {len(trajectories)} segment(s) found")


# =============================================================================
# Parameters JSON parsing
# =============================================================================


def read_parameters_json(sim_dir: Path, metadata: SimMetadata) -> None:
    """Read simulation parameters from legacy JSON files.

    Handles both older (flat thermo_params) and newer (nested thermostat_params)
    schema formats.

    Parameters
    ----------
    sim_dir : Path
        Path to the legacy simulation directory.
    metadata : SimMetadata
        Metadata object to update with parsed parameters.
    """
    # --- Equilibration parameters ---
    equil_json = sim_dir / "equilibration" / "equilibration_parameters.json"
    if equil_json.exists():
        metadata.equilibration = _parse_phase_json(equil_json)
        logger.info(
            f"  Equilibration: {metadata.equilibration.duration_ns} ns "
            f"{metadata.equilibration.ensemble}, "
            f"{metadata.equilibration.samples} samples"
        )

    # --- Production parameters ---
    # For segmented sims, read from the first segment
    prod_json = None
    if metadata.n_segments > 1 and metadata.production_dirs:
        seg0 = metadata.production_dirs[0]
        prod_json = seg0 / f"{seg0.name}_parameters.json"

    if prod_json is None or not prod_json.exists():
        prod_json = sim_dir / "production" / "production_parameters.json"

    if prod_json.exists():
        segment_params = _parse_phase_json(prod_json)

        # For segmented sims, the JSON shows per-segment duration.
        # We want the total across all segments.
        if metadata.n_segments > 1:
            total_duration = segment_params.duration_ns * metadata.n_segments
            total_samples = segment_params.samples * metadata.n_segments
            logger.info(
                f"  Production: {metadata.n_segments} segments x "
                f"{segment_params.duration_ns:.1f} ns = {total_duration:.1f} ns total, "
                f"{total_samples} total samples"
            )
            segment_params.duration_ns = total_duration
            segment_params.samples = total_samples
        else:
            logger.info(
                f"  Production: {segment_params.duration_ns} ns "
                f"{segment_params.ensemble}, "
                f"{segment_params.samples} samples"
            )

        metadata.production = segment_params


def _parse_phase_json(json_path: Path) -> PhaseParams:
    """Parse a single *_parameters.json file.

    Parameters
    ----------
    json_path : Path
        Path to the JSON file.

    Returns
    -------
    PhaseParams
        Parsed phase parameters.
    """
    with open(json_path) as f:
        data = json.load(f)

    vals = data["__values__"]
    integ = vals["integ_params"]["__values__"]
    thermo = vals["thermo_params"]["__values__"]

    time_step = integ["time_step"]["__values__"]["value"]
    duration = integ["total_time"]["__values__"]["value"]
    duration_unit = integ["total_time"]["__values__"]["unit"]
    samples = integ["num_samples"]

    # Convert duration to nanoseconds if needed
    if duration_unit == "nanosecond":
        duration_ns = duration
    elif duration_unit == "picosecond":
        duration_ns = duration / 1000.0
    else:
        duration_ns = duration  # assume ns

    # Handle both old and new thermo_params schema
    if "thermostat_params" in thermo:
        # New schema: nested thermostat_params / barostat_params
        thermo_p = thermo["thermostat_params"]["__values__"]
        temperature = thermo_p["temperature"]["__values__"]["value"]
        thermostat_type = thermo_p.get("thermostat", {}).get("__values__", "LangevinMiddle")
        thermostat = _normalize_thermostat(thermostat_type)

        pressure = 1.0
        barostat = "MC"
        barostat_freq = 25

        if "barostat_params" in thermo and thermo["barostat_params"] is not None:
            baro_p = thermo["barostat_params"]["__values__"]
            pressure = baro_p["pressure"]["__values__"]["value"]
            barostat_freq = baro_p.get("update_frequency", 25)
            baro_type = baro_p.get("barostat", {}).get("__values__", "MC")
            barostat = _normalize_barostat(baro_type)
            ensemble = "NPT"
        else:
            ensemble = "NVT"  # No barostat means NVT
    else:
        # Old schema: flat thermo_params
        ensemble = thermo.get("ensemble", "NPT")
        temperature = thermo["temperature"]["__values__"]["value"]
        pressure = thermo.get("pressure", {}).get("__values__", {}).get("value", 1.0)
        thermostat = "LangevinMiddle"
        barostat = "MC"
        barostat_freq = thermo.get("barostat_freq", 25)

    return PhaseParams(
        ensemble=ensemble,
        duration_ns=duration_ns,
        samples=samples,
        time_step_fs=time_step,
        temperature_K=temperature,
        pressure_atm=pressure,
        thermostat=thermostat,
        barostat=barostat,
        barostat_frequency=barostat_freq,
    )


def _normalize_thermostat(raw: str) -> str:
    """Normalize thermostat name to PolyzyMD convention."""
    mapping = {
        "LANGEVIN_MIDDLE": "LangevinMiddle",
        "LANGEVIN": "Langevin",
        "LangevinMiddle": "LangevinMiddle",
        "Langevin": "Langevin",
    }
    return mapping.get(raw, "LangevinMiddle")


def _normalize_barostat(raw: str) -> str:
    """Normalize barostat name to PolyzyMD convention."""
    mapping = {
        "MONTE_CARLO": "MC",
        "MC": "MC",
        "MonteCarlo": "MC",
    }
    return mapping.get(raw, "MC")


# =============================================================================
# Topology rewriting (chain ID + atom name fixing)
# =============================================================================


def rewrite_topology(
    metadata: SimMetadata,
    reference_pdb: Path,
    output_path: Path,
) -> dict[str, int]:
    """Rewrite topology PDB with corrected chain IDs and protein atom names.

    Chain assignment convention:
    - A: protein (standard amino acid residues)
    - B: substrate (resname RBY)
    - C: polymer (resnames SBM, EGM, EGP, or other known polymer residues)
    - D: solvent, ions, everything else (HOH, Na+, Cl-, etc.)

    Protein atom names are fixed using a reference PDB where atom names are
    correct. The mapping is 1:1 by atom index within the protein selection,
    validated by matching resid and resname.

    Parameters
    ----------
    metadata : SimMetadata
        Parsed simulation metadata with topology_path set.
    reference_pdb : Path
        Path to the reference PDB with correct protein atom names.
    output_path : Path
        Where to write the corrected solvated_system.pdb.

    Returns
    -------
    dict[str, int]
        Chain assignment summary: {"A": n_protein, "B": n_substrate, ...}
    """
    import MDAnalysis as mda
    from biopandas.pdb import PandasPdb

    logger.info(f"  Loading topology: {metadata.topology_path}")
    u = mda.Universe(str(metadata.topology_path))

    # --- Step 1: Fix protein atom names using reference PDB ---
    logger.info(f"  Fixing protein atom names from reference: {reference_pdb.name}")
    ppdb = PandasPdb().read_pdb(str(reference_pdb))
    ref_df = ppdb.df["ATOM"]

    prot_atoms = u.select_atoms("protein")

    if len(ref_df) != len(prot_atoms):
        raise ValueError(
            f"Reference PDB has {len(ref_df)} protein atoms, but topology has "
            f"{len(prot_atoms)} protein atoms. These must match exactly."
        )

    mismatches = []
    for k in range(len(prot_atoms)):
        atom = prot_atoms[k]
        ref_row = ref_df.iloc[k]
        if ref_row["residue_number"] != atom.resid or ref_row["residue_name"] != atom.resname:
            mismatches.append(
                f"  Atom {k}: ref=({ref_row['residue_name']} {ref_row['residue_number']}) "
                f"vs topo=({atom.resname} {atom.resid})"
            )

    if mismatches:
        raise ValueError(
            f"Reference PDB / topology protein mismatch at {len(mismatches)} atoms:\n"
            + "\n".join(mismatches[:10])
            + ("\n  ..." if len(mismatches) > 10 else "")
        )

    for k in range(len(prot_atoms)):
        atom = prot_atoms[k]
        ref_name = ref_df.iloc[k]["atom_name"]
        atom.name = ref_name
        atom.type = ref_name

    logger.info(f"  Fixed {len(prot_atoms)} protein atom names")

    # --- Step 2: Detect polymer residue names present in topology ---
    all_resnames = set(u.atoms.resnames)
    polymer_resnames = all_resnames & KNOWN_POLYMER_RESNAMES
    if metadata.polymer and not polymer_resnames:
        logger.warning(
            "  Folder name indicates polymers but no known polymer residues "
            f"({KNOWN_POLYMER_RESNAMES}) found in topology. "
            "Chain C will be empty."
        )
    elif polymer_resnames:
        logger.info(f"  Polymer residue names detected: {sorted(polymer_resnames)}")

    # --- Step 3: Assign chain IDs ---
    # Build selection strings for each chain
    chain_summary = {}

    # Chain A: protein
    protein_sel = u.select_atoms("protein")
    protein_sel.atoms.chainIDs = ["A"] * len(protein_sel)
    chain_summary["A (protein)"] = len(protein_sel)

    # Chain B: substrate (RBY)
    substrate_sel = u.select_atoms("resname RBY")
    if len(substrate_sel) > 0:
        substrate_sel.atoms.chainIDs = ["B"] * len(substrate_sel)
        chain_summary["B (substrate)"] = len(substrate_sel)
    else:
        logger.warning("  No substrate (resname RBY) found in topology")

    # Chain C: polymer residues
    if polymer_resnames:
        polymer_sel_str = " or ".join(f"resname {rn}" for rn in sorted(polymer_resnames))
        # Exclude atoms already assigned to protein or substrate
        polymer_sel = u.select_atoms(f"({polymer_sel_str}) and not protein and not resname RBY")
        if len(polymer_sel) > 0:
            polymer_sel.atoms.chainIDs = ["C"] * len(polymer_sel)
            chain_summary["C (polymer)"] = len(polymer_sel)

    # Chain D: everything else (solvent, ions)
    remaining = u.select_atoms("not protein and not resname RBY")
    if polymer_resnames:
        polymer_sel_str = " or ".join(f"resname {rn}" for rn in sorted(polymer_resnames))
        remaining = remaining.select_atoms(f"not ({polymer_sel_str})")

    if len(remaining) > 0:
        remaining.atoms.chainIDs = ["D"] * len(remaining)
        chain_summary["D (solvent/ions)"] = len(remaining)

        # Renumber solvent/ion residues sequentially.
        # Legacy PDBs have all water/ions as resid 1 (one molecule per chain).
        # After collapsing to chain D, we need unique resids for each residue.
        for i, res in enumerate(remaining.residues, start=1):
            res.atoms.residues.resids = i
        logger.info(
            f"  Renumbered {len(remaining.residues)} solvent/ion residues (1-{len(remaining.residues)})"
        )

    # --- Step 4: Write the corrected PDB ---
    logger.info(f"  Writing: {output_path}")
    u.atoms.write(str(output_path))

    for chain_label, count in chain_summary.items():
        logger.info(f"    {chain_label}: {count} atoms")

    return chain_summary


# =============================================================================
# config.yaml generation
# =============================================================================


def generate_config_yaml(
    metadata: SimMetadata,
    output_dir: Path,
    config_path: Path,
) -> None:
    """Generate a PolyzyMD-compatible config.yaml for a legacy simulation.

    The key trick is the naming_template: it is set to the literal legacy folder
    name with only {replicate} as a variable. This makes
    ``config.get_working_directory(replicate=N)`` resolve to the correct directory
    in the output folder.

    Parameters
    ----------
    metadata : SimMetadata
        Parsed simulation metadata.
    output_dir : Path
        Base output directory (parent of the simulation folder in the converted
        structure). This becomes ``scratch_directory``.
    config_path : Path
        Path where config.yaml will be written.
    """
    import yaml

    # Build the naming template: the literal folder name with {replicate} substituted
    # Replace the actual replicate number with {replicate} placeholder
    naming_template = re.sub(
        r"_run\d+$",
        "_run{replicate}",
        metadata.folder_name,
    )

    # Build config dict
    config = {
        "name": _build_sim_name(metadata),
        "description": (
            f"Legacy simulation converted for PolyzyMD analysis. Original: {metadata.folder_name}"
        ),
        # Enzyme
        "enzyme": {
            "name": metadata.enzyme_name,
            "pdb_path": "legacy_placeholder.pdb",  # Not used by analysis pipeline
        },
        # Substrate
        "substrate": {
            "name": metadata.substrate_name,
            "sdf_path": "legacy_placeholder.sdf",  # Not used by analysis pipeline
            "residue_name": "RBY",
        },
        # Thermodynamics
        "thermodynamics": {
            "temperature": metadata.temperature_K,
            "pressure": metadata.production.pressure_atm if metadata.production else 1.0,
        },
        # Simulation phases
        "simulation_phases": _build_phases_config(metadata),
        # Output — this is the critical section
        "output": {
            "projects_directory": str(output_dir),
            "scratch_directory": str(output_dir),
            "naming_template": naming_template,
        },
    }

    # Polymer config (if applicable)
    if metadata.polymer:
        prob_b = metadata.polymer.percent_b / 100.0
        prob_a = 1.0 - prob_b
        config["polymers"] = {
            "enabled": True,
            "generation_mode": "cached",
            "type_prefix": metadata.polymer.type_prefix,
            "sdf_directory": ".",  # Placeholder — not used by analysis pipeline
            "monomers": [
                {
                    "label": "A",
                    "name": metadata.polymer.monomer_a_name,
                    "probability": round(prob_a, 4),
                },
                {
                    "label": "B",
                    "name": metadata.polymer.monomer_b_name,
                    "probability": round(prob_b, 4),
                },
            ],
            "length": 5,  # 5-mers (from naming convention context)
            "count": metadata.polymer.chain_count,
        }
    else:
        config["polymers"] = None

    # Write YAML
    logger.info(f"  Writing: {config_path}")
    with open(config_path, "w") as f:
        f.write("# ============================================================================\n")
        f.write("# PolyzyMD Configuration — Auto-generated from legacy simulation\n")
        f.write("# ============================================================================\n")
        f.write(f"# Source: {metadata.folder_name}\n")
        f.write("# Generated by: scripts/convert_legacy.py\n")
        f.write("#\n")
        f.write("# NOTE: pdb_path, sdf_path, and sdf_directory are placeholders because\n")
        f.write("# this config is only used with current compare commands, not builds.\n")
        f.write(
            "# ============================================================================\n\n"
        )
        yaml.dump(config, f, default_flow_style=False, sort_keys=False, allow_unicode=True)


def _build_sim_name(metadata: SimMetadata) -> str:
    """Build a clean simulation name from metadata."""
    substrate_clean = metadata.substrate_name.replace("-", "")
    if metadata.polymer:
        return f"{metadata.enzyme_name}_{substrate_clean}_{metadata.polymer.type_prefix}"
    return f"{metadata.enzyme_name}_{substrate_clean}_control"


def _build_phases_config(metadata: SimMetadata) -> dict:
    """Build simulation_phases config section from metadata."""
    phases = {}

    if metadata.equilibration:
        eq = metadata.equilibration
        phases["equilibration_stages"] = [
            {
                "name": "legacy_equilibration",
                "ensemble": eq.ensemble,
                "duration": eq.duration_ns,
                "samples": eq.samples,
                "time_step": eq.time_step_fs,
                "temperature": eq.temperature_K,
                "thermostat": eq.thermostat,
                "thermostat_timescale": 1.0,
            }
        ]

    if metadata.production:
        prod = metadata.production
        phases["production"] = {
            "ensemble": prod.ensemble,
            "duration": prod.duration_ns,
            "samples": prod.samples,
            "report_interval": max(
                1,
                int(prod.duration_ns * 1e6 / prod.time_step_fs) // prod.samples,
            ),
            "checkpoint_interval": 60.0,
            "time_step": prod.time_step_fs,
            "thermostat": prod.thermostat,
            "thermostat_timescale": 1.0,
            "barostat": prod.barostat,
            "barostat_frequency": prod.barostat_frequency,
        }

    return phases


# =============================================================================
# Symlink creation
# =============================================================================


def create_symlinks(
    sim_dir: Path,
    output_sim_dir: Path,
    metadata: SimMetadata,
) -> None:
    """Create symlinks from output directory to original production directories.

    Parameters
    ----------
    sim_dir : Path
        Path to the original legacy simulation directory.
    output_sim_dir : Path
        Path to the output simulation directory (where solvated_system.pdb lives).
    metadata : SimMetadata
        Metadata with production directory info.
    """
    for prod_dir in metadata.production_dirs:
        link_path = output_sim_dir / prod_dir.name
        target = prod_dir.resolve()

        if link_path.exists() or link_path.is_symlink():
            if link_path.is_symlink() and link_path.resolve() == target:
                logger.info(f"  Symlink already exists: {prod_dir.name} -> {target}")
                continue
            else:
                logger.warning(
                    f"  {link_path} already exists (not a symlink to expected target). Skipping."
                )
                continue

        link_path.symlink_to(target)
        logger.info(f"  Symlink: {prod_dir.name} -> {target}")


# =============================================================================
# Validation
# =============================================================================


def validate_output(output_sim_dir: Path, metadata: SimMetadata) -> bool:
    """Validate the converted output directory.

    Checks:
    1. solvated_system.pdb exists and is loadable
    2. config.yaml exists and is parseable
    3. Trajectory files are discoverable via symlinks
    4. Chain assignments are correct

    Parameters
    ----------
    output_sim_dir : Path
        Path to the converted simulation directory.
    metadata : SimMetadata
        Original simulation metadata.

    Returns
    -------
    bool
        True if validation passes.
    """
    import MDAnalysis as mda

    success = True

    # Check solvated_system.pdb
    topo_path = output_sim_dir / "solvated_system.pdb"
    if not topo_path.exists():
        logger.error(f"  FAIL: solvated_system.pdb not found at {topo_path}")
        return False

    # Load and verify chain assignments
    try:
        u = mda.Universe(str(topo_path))
        chains = set(u.atoms.chainIDs)

        # Must have chain A (protein)
        prot = u.select_atoms("protein and chainID A")
        if len(prot) == 0:
            logger.error("  FAIL: No protein atoms on chain A")
            success = False
        else:
            logger.info(f"  OK: Chain A (protein): {len(prot)} atoms")

        # Check substrate
        sub = u.select_atoms("resname RBY and chainID B")
        if len(sub) == 0:
            logger.warning("  WARN: No substrate (RBY) atoms on chain B")
        else:
            logger.info(f"  OK: Chain B (substrate): {len(sub)} atoms")

        # Check polymer (if expected)
        if metadata.polymer:
            poly = u.select_atoms("chainID C")
            if len(poly) == 0:
                logger.error("  FAIL: No polymer atoms on chain C (expected from folder name)")
                success = False
            else:
                logger.info(f"  OK: Chain C (polymer): {len(poly)} atoms")

        # Check atom names on protein
        ca_atoms = u.select_atoms("protein and name CA")
        if len(ca_atoms) == 0:
            logger.warning(
                "  WARN: No 'CA' atoms found in protein selection. "
                "Atom name fixing may have failed. This could affect RMSF backbone selections."
            )
        else:
            logger.info(f"  OK: Protein atom names fixed ({len(ca_atoms)} CA atoms found)")

    except Exception as e:
        logger.error(f"  FAIL: Could not load solvated_system.pdb: {e}")
        success = False

    # Check config.yaml
    config_path = output_sim_dir / "config.yaml"
    if not config_path.exists():
        logger.error(f"  FAIL: config.yaml not found at {config_path}")
        success = False
    else:
        logger.info("  OK: config.yaml exists")

    # Check trajectory discovery via symlinks
    for prod_dir in metadata.production_dirs:
        link = output_sim_dir / prod_dir.name
        if not link.exists():
            logger.error(f"  FAIL: Production symlink missing: {link}")
            success = False
        else:
            # Check that trajectory files are accessible
            if metadata.n_segments > 1:
                traj_name = f"{prod_dir.name}_trajectory.dcd"
            else:
                traj_name = "production_trajectory.dcd"
            traj_path = link / traj_name
            if traj_path.exists():
                logger.info(f"  OK: Trajectory accessible: {prod_dir.name}/{traj_name}")
            else:
                logger.warning(f"  WARN: Trajectory not found at expected path: {traj_path}")

    # Try loading config with PolyzyMD (optional — may not be installed)
    try:
        from polyzymd.config.schema import SimulationConfig

        cfg = SimulationConfig.from_yaml(config_path)
        wd = cfg.get_working_directory(metadata.replicate)
        if wd.exists():
            logger.info(f"  OK: PolyzyMD config resolves working directory: {wd}")
        else:
            logger.warning(f"  WARN: Working directory does not exist: {wd}")
    except ImportError:
        logger.info("  SKIP: PolyzyMD not importable, skipping config validation")
    except Exception as e:
        logger.warning(f"  WARN: PolyzyMD config validation failed: {e}")

    return success


# =============================================================================
# Comparison YAML generation
# =============================================================================


def _condition_key(metadata: SimMetadata) -> tuple[str, str, float, str | None, float | None]:
    """Return a hashable key that groups replicates into a single condition.

    The key captures everything about the simulation *except* replicate number:
    (enzyme, substrate, temperature, polymer_type_prefix, polymer_percent_b).

    Parameters
    ----------
    metadata : SimMetadata
        Parsed simulation metadata.

    Returns
    -------
    tuple
        Hashable condition key.
    """
    poly_prefix = metadata.polymer.type_prefix if metadata.polymer else None
    poly_pct = metadata.polymer.percent_b if metadata.polymer else None
    return (
        metadata.enzyme_name,
        metadata.substrate_name,
        metadata.temperature_K,
        poly_prefix,
        poly_pct,
    )


def _condition_label(metadata: SimMetadata) -> str:
    """Generate a human-readable label for a condition.

    Parameters
    ----------
    metadata : SimMetadata
        Parsed simulation metadata (any replicate from the condition).

    Returns
    -------
    str
        Display label, e.g., "No Polymer (Control)", "SBMA-EGMA 50%",
        "pSBMA Homopolymer".
    """
    if metadata.polymer is None:
        return "No Polymer (Control)"

    poly = metadata.polymer
    # Special case: 0% second component = pSBMA homopolymer (pure SBMA)
    if poly.percent_b == 0.0 and poly.type_prefix in ("SBMA-EGPMA", "SBMA-EGMA"):
        return "pSBMA Homopolymer"

    pct = int(poly.percent_b) if poly.percent_b == int(poly.percent_b) else poly.percent_b
    return f"{poly.type_prefix} {pct}%"


def group_conditions(
    output_dir: Path,
    temperature_filter: float | None = None,
) -> dict[str, dict]:
    """Scan converted output directory and group simulations by condition.

    Parameters
    ----------
    output_dir : Path
        Base output directory containing converted simulation subdirectories.
    temperature_filter : float, optional
        If provided, only include simulations at this temperature (K).

    Returns
    -------
    dict[str, dict]
        Mapping of condition label to:
        - ``"config"`` : Path to run1's config.yaml (canonical config)
        - ``"replicates"`` : sorted list of replicate ints
        - ``"metadata"`` : SimMetadata of the first replicate (for reference)
    """
    # Find all subdirectories that have config.yaml (i.e., successfully converted)
    sim_dirs = sorted(
        [
            d
            for d in output_dir.iterdir()
            if d.is_dir() and (d / "config.yaml").exists() and "_RESTRAINT_" in d.name
        ]
    )

    if not sim_dirs:
        raise FileNotFoundError(
            f"No converted simulation directories with config.yaml found in {output_dir}"
        )

    # Group by condition
    conditions: dict[tuple, dict] = {}
    skipped = 0

    for sim_dir in sim_dirs:
        try:
            metadata = parse_folder_name(sim_dir.name)
        except ValueError:
            skipped += 1
            continue

        # Temperature filter
        if temperature_filter is not None and metadata.temperature_K != temperature_filter:
            skipped += 1
            continue

        key = _condition_key(metadata)
        if key not in conditions:
            label = _condition_label(metadata)
            conditions[key] = {
                "label": label,
                "replicates": [],
                "config_paths": {},  # replicate -> config path
                "metadata": metadata,
            }

        conditions[key]["replicates"].append(metadata.replicate)
        conditions[key]["config_paths"][metadata.replicate] = sim_dir / "config.yaml"

    if not conditions:
        temp_msg = f" at {temperature_filter} K" if temperature_filter else ""
        raise FileNotFoundError(
            f"No conditions found{temp_msg} in {output_dir}. "
            f"Scanned {len(sim_dirs)} directories, skipped {skipped}."
        )

    # Finalize: sort replicates, pick run1 (or lowest) as canonical config
    result = {}
    for key, cond in conditions.items():
        cond["replicates"].sort()
        # Use the lowest replicate's config as the canonical one
        lowest_rep = cond["replicates"][0]
        cond["config"] = cond["config_paths"][lowest_rep]
        del cond["config_paths"]
        result[cond["label"]] = cond

    if skipped > 0:
        logger.info(f"  Skipped {skipped} directories (parse error or temperature filter)")

    return result


def generate_comparison_yaml(
    output_dir: Path,
    comparison_dir: Path,
    project_name: str = "legacy_polymer_comparison",
    control_label: str = "No Polymer (Control)",
    temperature_filter: float | None = None,
    enzyme_pdb_path: Path | None = None,
) -> Path:
    """Generate a complete comparison.yaml and project directory.

    Creates the comparison project directory structure::

        <comparison_dir>/
        ├── comparison.yaml
        ├── structures/
        │   └── enzyme.pdb  (copied if enzyme_pdb_path provided)
        ├── results/
        └── figures/

    Parameters
    ----------
    output_dir : Path
        Base output directory containing converted simulation subdirectories.
    comparison_dir : Path
        Path where the comparison project directory will be created.
    project_name : str
        Name for the comparison project.
    control_label : str
        Label of the control condition (must match a generated condition label).
    temperature_filter : float, optional
        Only include simulations at this temperature (K).
    enzyme_pdb_path : Path, optional
        Path to enzyme PDB for SASA calculations. Will be copied to
        ``structures/enzyme.pdb`` in the comparison directory.

    Returns
    -------
    Path
        Path to the generated comparison.yaml file.
    """
    import shutil

    import yaml

    logger.info(f"\n{'=' * 80}")
    logger.info("GENERATING COMPARISON PROJECT")
    logger.info(f"{'=' * 80}")

    # --- Group conditions ---
    conditions = group_conditions(output_dir, temperature_filter)
    logger.info(f"  Found {len(conditions)} conditions:")
    for label, cond in sorted(conditions.items()):
        reps = cond["replicates"]
        logger.info(f"    {label}: replicates {reps}")

    # --- Validate control ---
    if control_label not in conditions:
        available = ", ".join(sorted(conditions.keys()))
        raise ValueError(
            f"Control label '{control_label}' not found in conditions.\nAvailable: {available}"
        )

    # --- Create project directory structure ---
    comparison_dir.mkdir(parents=True, exist_ok=True)
    (comparison_dir / "structures").mkdir(exist_ok=True)
    (comparison_dir / "results").mkdir(exist_ok=True)
    (comparison_dir / "figures").mkdir(exist_ok=True)

    # --- Copy enzyme PDB if provided ---
    enzyme_pdb_rel = None
    if enzyme_pdb_path is not None:
        dest = comparison_dir / "structures" / "enzyme.pdb"
        if not dest.exists():
            shutil.copy2(enzyme_pdb_path, dest)
            logger.info(f"  Copied enzyme PDB to {dest}")
        else:
            logger.info(f"  Enzyme PDB already exists at {dest}")
        enzyme_pdb_rel = "structures/enzyme.pdb"

    # --- Build conditions list for YAML ---
    # Config paths are relative to comparison_dir
    conditions_yaml = []

    # Control first, then alphabetical
    sorted_labels = sorted(conditions.keys())
    ordered_labels = [control_label] + [lb for lb in sorted_labels if lb != control_label]

    for label in ordered_labels:
        cond = conditions[label]
        config_abs = cond["config"].resolve()
        # Make path relative to comparison_dir
        try:
            config_rel = config_abs.relative_to(comparison_dir.resolve())
            config_str = str(config_rel)
        except ValueError:
            # Not a subpath — use relative path with ../
            try:
                config_str = str(
                    Path("..") / config_abs.relative_to(comparison_dir.resolve().parent)
                )
            except ValueError:
                # Fallback to absolute path
                config_str = str(config_abs)

        conditions_yaml.append(
            {
                "label": label,
                "config": config_str,
                "replicates": cond["replicates"],
            }
        )

    # --- Build the full comparison config ---
    comparison_config = {
        "name": project_name,
        "description": (
            f"Comparison of enzyme-polymer conjugates"
            f"{f' at {temperature_filter:.0f} K' if temperature_filter else ''}"
            f" — {len(conditions)} conditions, generated by convert_legacy.py"
        ),
        "control": control_label,
        "conditions": conditions_yaml,
        "defaults": {
            "equilibration_time": "10ns",
        },
        "plugins": _build_plugin_settings(enzyme_pdb_rel),
        "plot_settings": _build_plot_settings(),
    }

    # --- Write comparison.yaml ---
    comparison_yaml_path = comparison_dir / "comparison.yaml"
    logger.info(f"\n  Writing: {comparison_yaml_path}")

    with open(comparison_yaml_path, "w") as f:
        f.write("# ============================================================================\n")
        f.write(f"# PolyzyMD Comparison — {project_name}\n")
        f.write("# ============================================================================\n")
        f.write("# Generated by: scripts/convert_legacy.py --generate-comparison\n")
        if temperature_filter:
            f.write(f"# Temperature filter: {temperature_filter:.0f} K\n")
        f.write(f"# Control: {control_label}\n")
        f.write(f"# Conditions: {len(conditions)}\n")
        f.write(
            "# ============================================================================\n\n"
        )
        yaml.dump(
            comparison_config,
            f,
            default_flow_style=False,
            sort_keys=False,
            allow_unicode=True,
        )

    # --- Summary ---
    logger.info(f"\n  Comparison project created at: {comparison_dir}")
    logger.info(f"  comparison.yaml: {comparison_yaml_path}")
    logger.info(f"\n  Conditions ({len(conditions)}):")
    for entry in conditions_yaml:
        logger.info(f"    - {entry['label']}: replicates {entry['replicates']}")
        logger.info(f"      config: {entry['config']}")

    logger.info(
        f"\n  Next steps:"
        f"\n    1. Review and edit {comparison_yaml_path} if needed"
        f"\n    2. Ensure enzyme PDB is at {comparison_dir}/structures/enzyme.pdb"
        f"\n    3. Validate the comparison configuration:"
        f"\n       polyzymd compare validate -f {comparison_yaml_path}"
        f"\n    4. Run analyses:"
        f"\n       polyzymd compare run rmsf -f {comparison_yaml_path}"
        f"\n       polyzymd compare run catalytic_triad -f {comparison_yaml_path}"
        f"\n       polyzymd compare run contacts -f {comparison_yaml_path}"
        f"\n    5. Generate all plots:"
        f"\n       polyzymd compare plot-all -f {comparison_yaml_path}"
    )

    return comparison_yaml_path


def _build_plugin_settings(enzyme_pdb_rel: str | None) -> dict:
    """Build the plugins section of comparison.yaml.

    Parameters
    ----------
    enzyme_pdb_rel : str or None
        Relative path to enzyme PDB retained for backward-compatible calls.

    Returns
    -------
    dict
        Plugin settings configuration.
    """
    del enzyme_pdb_rel

    settings: dict = {
        # RMSF: per-residue fluctuations (backbone CA atoms)
        "rmsf": {
            "selection": "protein and name CA",
            "reference_mode": "average",
        },
        # Catalytic triad: Ser77-His156-Asp133 distances
        "catalytic_triad": {
            "name": "CalB Catalytic Triad",
            "description": "Ser77(OG)-His156(NE2) and His156(ND1)-Asp133(OD2) distances",
            "threshold": 3.5,
            "pairs": [
                {
                    "label": "Ser77-His156",
                    "selection_a": "resid 77 and name OG",
                    "selection_b": "resid 156 and name NE2",
                },
                {
                    "label": "His156-Asp133",
                    "selection_a": "resid 156 and name ND1",
                    "selection_b": "resid 133 and name OD2",
                },
            ],
        },
        # Contacts: polymer-protein contacts
        "contacts": {
            "polymer_selection": "chainid C",
            "protein_selection": "chainid A",
            "cutoff": 4.5,
            "grouping": "aa_class",
            "compute_residence_times": True,
            "protein_groups": {
                "catalytic_triad": [77, 133, 156],
            },
        },
    }

    return settings


def _build_plot_settings() -> dict:
    """Build the plot_settings section of comparison.yaml.

    Returns
    -------
    dict
        Plot settings configuration.
    """
    return {
        "output_dir": "figures/",
        "format": "png",
        "dpi": 300,
        "style": "compact",
        "color_palette": "tab10",
        "rmsf": {
            "show_error": True,
            "highlight_residues": [77, 133, 156],
            "figsize_profile": [14, 4],
            "figsize_comparison": [8, 6],
        },
        "triad": {
            "generate_kde_panel": True,
            "generate_bars": True,
            "figsize_bars": [10, 6],
        },
        "contacts": {
            "figsize": [10, 8],
        },
    }


# =============================================================================
# Main conversion logic
# =============================================================================


def convert_simulation(
    sim_dir: Path,
    output_dir: Path,
    reference_pdb: Path,
    dry_run: bool = False,
) -> bool:
    """Convert a single legacy simulation directory.

    Parameters
    ----------
    sim_dir : Path
        Path to the legacy simulation directory.
    output_dir : Path
        Base output directory for converted files.
    reference_pdb : Path
        Path to the reference PDB with correct protein atom names.
    dry_run : bool
        If True, only print what would be done without writing files.

    Returns
    -------
    bool
        True if conversion (or dry-run) succeeded.
    """
    folder_name = sim_dir.name
    logger.info(f"\n{'=' * 80}")
    logger.info(f"Converting: {folder_name}")
    logger.info(f"{'=' * 80}")

    # Step 1: Parse folder name
    try:
        metadata = parse_folder_name(folder_name)
    except ValueError as e:
        logger.error(f"  SKIP: {e}")
        return False

    logger.info(f"  Enzyme: {metadata.enzyme_name}")
    logger.info(f"  Substrate: {metadata.substrate_name}")
    logger.info(f"  Temperature: {metadata.temperature_K} K")
    logger.info(f"  Replicate: {metadata.replicate}")
    if metadata.polymer:
        logger.info(
            f"  Polymer: {metadata.polymer.type_prefix} "
            f"({metadata.polymer.percent_b}% {metadata.polymer.monomer_b_name}), "
            f"{metadata.polymer.chain_count} chains"
        )
    else:
        logger.info("  Polymer: none (control)")

    # Step 2: Discover files
    try:
        discover_files(sim_dir, metadata)
    except FileNotFoundError as e:
        logger.error(f"  SKIP: {e}")
        return False

    # Step 3: Read parameters JSON
    read_parameters_json(sim_dir, metadata)

    # Output paths
    output_sim_dir = output_dir / folder_name
    solvated_pdb = output_sim_dir / "solvated_system.pdb"
    config_yaml = output_sim_dir / "config.yaml"

    if dry_run:
        logger.info("\n  [DRY RUN] Would create:")
        logger.info(f"    {output_sim_dir}/")
        logger.info(f"    {solvated_pdb}")
        logger.info(f"    {config_yaml}")
        for prod_dir in metadata.production_dirs:
            logger.info(f"    {output_sim_dir / prod_dir.name} -> {prod_dir.resolve()}")
        return True

    # Step 4: Create output directory
    output_sim_dir.mkdir(parents=True, exist_ok=True)

    # Step 5: Rewrite topology
    try:
        rewrite_topology(metadata, reference_pdb, solvated_pdb)
    except Exception as e:
        logger.error(f"  FAIL: Topology rewrite failed: {e}")
        return False

    # Step 6: Generate config.yaml
    try:
        generate_config_yaml(metadata, output_dir, config_yaml)
    except Exception as e:
        logger.error(f"  FAIL: Config generation failed: {e}")
        return False

    # Step 7: Create symlinks
    create_symlinks(sim_dir, output_sim_dir, metadata)

    # Step 8: Validate
    logger.info("\n  Validation:")
    valid = validate_output(output_sim_dir, metadata)

    if valid:
        logger.info(f"\n  SUCCESS: {folder_name}")
    else:
        logger.warning(f"\n  WARNINGS during conversion of {folder_name}")

    return True


def find_legacy_sim_dirs(parent_dir: Path) -> list[Path]:
    """Find all legacy simulation directories in a parent directory.

    Identifies directories matching the legacy naming convention by
    checking for the ``_RESTRAINT_`` pattern and presence of production
    subdirectories.

    Parameters
    ----------
    parent_dir : Path
        Parent directory to search.

    Returns
    -------
    list[Path]
        Sorted list of legacy simulation directory paths.
    """
    sim_dirs = []
    for d in sorted(parent_dir.iterdir()):
        if not d.is_dir():
            continue
        # Must contain _RESTRAINT_ and _run pattern
        if "_RESTRAINT_" in d.name and re.search(r"_run\d+$", d.name):
            # Must have a production directory
            has_production = (d / "production").is_dir() or any(d.glob("production_[0-9]*"))
            if has_production:
                sim_dirs.append(d)
    return sim_dirs


# =============================================================================
# CLI
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Convert pre-PolyzyMD legacy simulation directories for PolyzyMD analysis.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "input_path",
        type=Path,
        nargs="?",
        default=None,
        help=(
            "Path to a single legacy simulation directory, or (with --batch) "
            "the parent directory containing multiple legacy simulations. "
            "Not required when using --generate-comparison alone."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Base output directory for converted files (will be created if needed).",
    )
    parser.add_argument(
        "--reference-pdb",
        type=Path,
        default=None,
        help="Path to reference PDB with correct protein atom names (required for conversion).",
    )
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Treat input_path as a parent directory and convert all legacy sims within.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview what would be done without writing any files.",
    )

    # Comparison generation options
    parser.add_argument(
        "--generate-comparison",
        type=Path,
        metavar="DIR",
        help=(
            "After conversion, generate a comparison project at DIR. "
            "Scans --output-dir for converted simulations and groups them "
            "into conditions. Creates comparison.yaml with all analysis types."
        ),
    )
    parser.add_argument(
        "--comparison-name",
        type=str,
        default="legacy_polymer_comparison",
        help="Name for the comparison project (default: legacy_polymer_comparison).",
    )
    parser.add_argument(
        "--control",
        type=str,
        default="No Polymer (Control)",
        help='Label of the control condition (default: "No Polymer (Control)").',
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=None,
        help="Only include simulations at this temperature (K) in the comparison.",
    )
    parser.add_argument(
        "--enzyme-pdb",
        type=Path,
        default=None,
        help="Path to enzyme PDB for SASA calculations. Copied to comparison structures/.",
    )

    args = parser.parse_args()

    # Determine mode: conversion, comparison-only, or both
    do_conversion = args.input_path is not None
    do_comparison = args.generate_comparison is not None

    if not do_conversion and not do_comparison:
        parser.error(
            "Either provide input_path for conversion, or use "
            "--generate-comparison for comparison-only mode."
        )

    # Validate conversion inputs
    if do_conversion:
        if not args.input_path.exists():
            logger.error(f"Input path does not exist: {args.input_path}")
            sys.exit(1)

        if args.reference_pdb is None:
            parser.error("--reference-pdb is required for conversion.")

        if not args.reference_pdb.exists():
            logger.error(f"Reference PDB does not exist: {args.reference_pdb}")
            sys.exit(1)

    # Resolve paths
    output_dir = args.output_dir.resolve()
    input_path = args.input_path.resolve() if do_conversion else None
    reference_pdb = args.reference_pdb.resolve() if do_conversion else None

    if args.dry_run:
        logger.info("[DRY RUN MODE] No files will be written.\n")

    # --- Conversion phase ---
    n_success = 0
    n_total = 0

    if do_conversion:
        # Find simulation directories
        if args.batch:
            sim_dirs = find_legacy_sim_dirs(input_path)
            if not sim_dirs:
                logger.error(f"No legacy simulation directories found in {input_path}")
                sys.exit(1)
            logger.info(f"Found {len(sim_dirs)} legacy simulation(s):")
            for d in sim_dirs:
                logger.info(f"  - {d.name}")
        else:
            if not input_path.is_dir():
                logger.error(f"Input path is not a directory: {input_path}")
                sys.exit(1)
            sim_dirs = [input_path]

        # Create output directory
        if not args.dry_run:
            output_dir.mkdir(parents=True, exist_ok=True)

        # Convert each simulation
        results = {}
        for sim_dir in sim_dirs:
            success = convert_simulation(sim_dir, output_dir, reference_pdb, args.dry_run)
            results[sim_dir.name] = success

        # Summary
        logger.info(f"\n{'=' * 80}")
        logger.info("CONVERSION SUMMARY")
        logger.info(f"{'=' * 80}")
        n_success = sum(1 for v in results.values() if v)
        n_total = len(results)
        for name, success in results.items():
            status = "OK" if success else "FAILED"
            logger.info(f"  [{status}] {name}")
        logger.info(f"\n  {n_success}/{n_total} simulations converted successfully")

    if not args.dry_run:
        if do_conversion:
            logger.info(f"\n  Output directory: {output_dir}")

        # Generate comparison project if requested
        if do_comparison:
            comparison_dir = args.generate_comparison.resolve()
            enzyme_pdb = args.enzyme_pdb.resolve() if args.enzyme_pdb else None

            if enzyme_pdb and not enzyme_pdb.exists():
                logger.error(f"Enzyme PDB does not exist: {enzyme_pdb}")
                sys.exit(1)

            try:
                generate_comparison_yaml(
                    output_dir=output_dir,
                    comparison_dir=comparison_dir,
                    project_name=args.comparison_name,
                    control_label=args.control,
                    temperature_filter=args.temperature,
                    enzyme_pdb_path=enzyme_pdb,
                )
            except (FileNotFoundError, ValueError) as e:
                logger.error(f"Comparison generation failed: {e}")
                sys.exit(1)
        elif do_conversion:
            logger.info(
                "\n  Next steps:"
                "\n    1. Re-run with --generate-comparison <dir> to create comparison.yaml"
                "\n    2. Review comparison.yaml and point each condition to its config.yaml"
                "\n    3. Validate the comparison configuration:"
                "\n       polyzymd compare validate -f comparison.yaml"
                "\n    4. Run analyses by canonical plugin name:"
                "\n       polyzymd compare run <analysis_type> -f comparison.yaml"
                "\n    5. Generate figures:"
                "\n       polyzymd compare plot-all -f comparison.yaml"
            )

    # Exit code: success if no conversion was done, or all conversions passed
    if do_conversion:
        sys.exit(0 if n_success == n_total else 1)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()
