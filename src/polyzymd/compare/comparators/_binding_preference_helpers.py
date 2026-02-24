"""Shared helpers for computing binding preference from contacts data.

Extracted from ContactsComparator so that both ContactsComparator and
BindingFreeEnergyComparator can compute binding preference without
duplicating logic.

Both comparators follow the same contract:
1. Resolve enzyme PDB → compute SASA-based surface exposure
2. Resolve protein groups from surface exposure
3. Extract polymer composition from topology
4. For each replicate: load contacts_rep{N}.json → compute_binding_preference
5. Aggregate across replicates → save aggregated result

This module provides three public functions:
- compute_condition_binding_preference()
- find_enzyme_pdb()
- try_load_cached_binding_preference()
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from polyzymd.analysis.contacts.binding_preference import (
        AggregatedBindingPreferenceResult,
        BindingPreferenceResult,
    )
    from polyzymd.compare.config import ConditionConfig

logger = logging.getLogger("polyzymd.compare")


def find_enzyme_pdb(sim_config: Any) -> Path | None:
    """Find enzyme PDB file from simulation config.

    Searches common locations relative to the project directory.

    Parameters
    ----------
    sim_config : SimulationConfig
        Simulation configuration (must have ``output.projects_directory``).

    Returns
    -------
    Path or None
        Path to enzyme PDB, or None if not found.
    """
    import glob as glob_module

    project_dir = sim_config.output.projects_directory
    possible_paths = [
        project_dir / "structures" / "enzyme.pdb",
        project_dir / "input" / "enzyme.pdb",
        project_dir.parent / "structures" / "enzyme.pdb",
        project_dir.parent / "enzyme.pdb",
    ]

    for path in possible_paths:
        if path.exists():
            return path

    # Try glob for any PDB with "enzyme" in name
    patterns = [
        str(project_dir / "**" / "*enzyme*.pdb"),
        str(project_dir.parent / "*enzyme*.pdb"),
    ]
    for pattern in patterns:
        matches = glob_module.glob(pattern, recursive=True)
        if matches:
            return Path(matches[0])

    return None


def resolve_enzyme_pdb(
    enzyme_pdb_setting: str | None,
    source_path: Path | None,
    sim_config: Any,
) -> Path | None:
    """Resolve the enzyme PDB path from settings or auto-discovery.

    Parameters
    ----------
    enzyme_pdb_setting : str or None
        Explicit enzyme PDB path from analysis settings (e.g.,
        ``enzyme_pdb_for_sasa``). If relative, resolved against
        *source_path*'s parent directory.
    source_path : Path or None
        Path to the comparison.yaml file (used to resolve relative paths).
    sim_config : Any
        Simulation configuration for auto-discovery fallback.

    Returns
    -------
    Path or None
        Resolved enzyme PDB path, or None if not found.
    """
    if enzyme_pdb_setting:
        if source_path:
            enzyme_pdb = source_path.parent / enzyme_pdb_setting
        else:
            enzyme_pdb = Path(enzyme_pdb_setting)
        if enzyme_pdb.exists():
            return enzyme_pdb
        logger.warning(f"Explicit enzyme_pdb_for_sasa not found at {enzyme_pdb}")
        return None

    return find_enzyme_pdb(sim_config)


def compute_condition_binding_preference(
    cond: "ConditionConfig",
    sim_config: Any,
    analysis_dir: Path,
    *,
    enzyme_pdb: Path,
    threshold: float = 0.2,
    include_default_aa_groups: bool = True,
    custom_protein_groups: dict[str, list[int]] | None = None,
    protein_partitions: dict[str, list[str]] | None = None,
    polymer_type_selections: dict[str, str] | None = None,
) -> "AggregatedBindingPreferenceResult | None":
    """Compute binding preference for a condition from contacts data.

    Computes surface exposure from enzyme PDB, resolves protein groups,
    extracts polymer composition, and calculates binding preference
    enrichment for each replicate.

    This function requires that per-replicate contacts files
    (``contacts_rep{N}.json``) already exist in *analysis_dir*. It does
    **not** trigger contacts computation.

    Parameters
    ----------
    cond : ConditionConfig
        Condition to compute (provides label and replicate list).
    sim_config : Any
        Simulation configuration (provides working directory for topology).
    analysis_dir : Path
        Contacts analysis directory containing ``contacts_rep{N}.json``
        files and where output will be saved.
    enzyme_pdb : Path
        Path to enzyme PDB for SASA calculation.
    threshold : float, optional
        Relative SASA threshold for surface-exposed residues. Default 0.2.
    include_default_aa_groups : bool, optional
        Include default AA class groupings. Default True.
    custom_protein_groups : dict or None, optional
        Custom protein groups as ``{name: [resid, ...]}``.
    protein_partitions : dict or None, optional
        Custom partitions as ``{partition_name: [group1, ...]}``.
    polymer_type_selections : dict or None, optional
        Custom polymer type selections as ``{name: "MDAnalysis selection"}``.

    Returns
    -------
    AggregatedBindingPreferenceResult or None
        Computed and aggregated result, or None on failure.
    """
    import MDAnalysis as mda

    from polyzymd.analysis.contacts.binding_preference import (
        PolymerComposition,
        aggregate_binding_preference,
        compute_binding_preference,
        extract_polymer_composition,
        resolve_protein_groups_from_surface_exposure,
    )
    from polyzymd.analysis.contacts.results import ContactResult
    from polyzymd.analysis.contacts.surface_exposure import SurfaceExposureFilter

    # --- Step 1: Compute surface exposure ---
    try:
        exposure_filter = SurfaceExposureFilter(threshold=threshold)
        surface_exposure = exposure_filter.calculate(str(enzyme_pdb))
        logger.debug(
            f"Computed surface exposure for {cond.label}: "
            f"{surface_exposure.exposed_count}/{surface_exposure.total_count} residues exposed"
        )
    except Exception as e:
        logger.warning(f"Failed to compute surface exposure for {cond.label}: {e}")
        return None

    # --- Step 2: Resolve protein groups ---
    protein_groups = resolve_protein_groups_from_surface_exposure(
        surface_exposure,
        include_default_aa_groups=include_default_aa_groups,
        custom_protein_groups=custom_protein_groups,
    )

    if not protein_groups:
        logger.warning(f"No protein groups resolved for {cond.label}")
        return None

    # --- Step 3: Extract polymer composition from first replicate topology ---
    polymer_composition = None
    first_rep = cond.replicates[0] if cond.replicates else 1
    run_dir = sim_config.get_working_directory(first_rep)
    topology_path = run_dir / "solvated_system.pdb"

    if topology_path.exists():
        try:
            universe = mda.Universe(str(topology_path))
            polymer_composition = extract_polymer_composition(universe, polymer_type_selections)
            logger.debug(
                f"Extracted polymer composition for {cond.label}: "
                f"{polymer_composition.total_residues} residues, "
                f"{polymer_composition.total_heavy_atoms} heavy atoms"
            )
        except Exception as e:
            logger.warning(f"Failed to extract polymer composition for {cond.label}: {e}")
    else:
        logger.warning(
            f"Cannot extract polymer composition for {cond.label}: "
            f"topology not found at {topology_path}"
        )

    if polymer_composition is None:
        polymer_composition = PolymerComposition()
        logger.warning(
            f"Using empty polymer composition for {cond.label} — enrichment ratios will be NaN"
        )

    # --- Step 4: Compute binding preference per replicate ---
    rep_results = []
    for rep in cond.replicates:
        contact_path = analysis_dir / f"contacts_rep{rep}.json"
        if not contact_path.exists():
            logger.warning(f"Contacts file not found: {contact_path}")
            continue

        try:
            contact_result = ContactResult.load(contact_path)
            bp_result = compute_binding_preference(
                contact_result=contact_result,
                surface_exposure=surface_exposure,
                protein_groups=protein_groups,
                polymer_composition=polymer_composition,
                protein_partitions=protein_partitions,
            )
            rep_results.append(bp_result)

            # Save per-replicate result
            rep_bp_path = analysis_dir / f"binding_preference_rep{rep}.json"
            bp_result.save(rep_bp_path)
            logger.debug(f"Computed and saved binding preference for {cond.label} rep{rep}")

        except Exception as e:
            logger.warning(f"Failed to compute binding preference for {cond.label} rep{rep}: {e}")
            continue

    if not rep_results:
        logger.warning(f"No binding preference results computed for {cond.label}")
        return None

    # --- Step 5: Aggregate and save ---
    agg_result = aggregate_binding_preference(rep_results)
    rep_range = f"{min(cond.replicates)}-{max(cond.replicates)}"
    agg_path = analysis_dir / f"binding_preference_aggregated_reps{rep_range}.json"
    agg_result.save(agg_path)
    logger.info(
        f"Computed binding preference for {cond.label}: "
        f"{len(rep_results)} replicates, {len(protein_groups)} protein groups"
    )

    return agg_result


def try_load_cached_binding_preference(
    cond: "ConditionConfig",
    analysis_dir: Path,
) -> "AggregatedBindingPreferenceResult | BindingPreferenceResult | None":
    """Try to load cached binding preference results for a condition.

    Searches for binding preference files in order of preference:
    1. binding_preference_aggregated.json
    2. binding_preference_aggregated_reps*.json (glob pattern)
    3. binding_preference.json (single replicate)
    4. Per-replicate files (binding_preference_rep{N}.json)

    Parameters
    ----------
    cond : ConditionConfig
        Condition to load.
    analysis_dir : Path
        Analysis directory for this condition.

    Returns
    -------
    AggregatedBindingPreferenceResult | BindingPreferenceResult | None
        Loaded result, or None if not found.
    """
    import glob as glob_module

    from polyzymd.analysis.contacts.binding_preference import (
        AggregatedBindingPreferenceResult,
        BindingPreferenceResult,
        aggregate_binding_preference,
    )

    # Try aggregated result first (multi-replicate)
    agg_path = analysis_dir / "binding_preference_aggregated.json"
    if agg_path.exists():
        result = AggregatedBindingPreferenceResult.load(agg_path)
        logger.debug(f"Loaded aggregated binding preference for {cond.label}")
        return result

    # Try aggregated result with rep range in name (e.g., _reps1-3.json)
    agg_pattern = str(analysis_dir / "binding_preference_aggregated_reps*.json")
    agg_matches = sorted(glob_module.glob(agg_pattern))
    if agg_matches:
        result = AggregatedBindingPreferenceResult.load(agg_matches[-1])
        logger.debug(f"Loaded aggregated binding preference for {cond.label}")
        return result

    # Try single replicate result
    single_path = analysis_dir / "binding_preference.json"
    if single_path.exists():
        result = BindingPreferenceResult.load(single_path)
        logger.debug(f"Loaded single binding preference for {cond.label}")
        return result

    # Try per-replicate results and aggregate them
    rep_results = []
    for rep in cond.replicates:
        rep_path = analysis_dir / f"binding_preference_rep{rep}.json"
        if rep_path.exists():
            rep_results.append(BindingPreferenceResult.load(rep_path))

    if rep_results:
        agg_result = aggregate_binding_preference(rep_results)
        logger.debug(
            f"Aggregated {len(rep_results)} replicate binding preference results for {cond.label}"
        )
        return agg_result

    return None
