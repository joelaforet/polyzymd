"""Config-driven orchestration for binding preference workflows."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Mapping, Protocol, Sequence

from ._compute import compute_binding_preference
from ._models import PolymerComposition
from ._polymer import extract_polymer_composition
from ._resolution import (
    resolve_polymer_type_selections,
    resolve_protein_group_selections,
    resolve_protein_groups_from_surface_exposure,
)

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

    from ._models import AggregatedBindingPreferenceResult, BindingPreferenceResult

logger = logging.getLogger(__name__)


class ConditionLike(Protocol):
    """Minimal condition protocol required by orchestration helpers.

    Attributes
    ----------
    label : str
        Condition label.
    replicates : Sequence[int]
        Replicate IDs associated with this condition.
    """

    label: str
    replicates: Sequence[int]


def compute_condition_binding_preference(
    cond: ConditionLike,
    sim_config: Any,
    analysis_dir: Path,
    *,
    enzyme_pdb: Path,
    contact_results_by_replicate: Mapping[int, Path],
    load_contact_result: Callable[[Path], Any],
    threshold: float,
    include_default_aa_groups: bool,
    custom_protein_groups: dict[str, list[int]] | None,
    protein_partitions: dict[str, list[str]] | None,
    polymer_type_selections: dict[str, str] | None,
    polymer_chain: str,
    settings_fp: str | None,
) -> "AggregatedBindingPreferenceResult | None":
    """Compute condition-level binding preference from contacts data.

    Parameters
    ----------
    cond : ConditionLike
        Condition to compute.
    sim_config : Any
        Simulation config used for topology lookup.
    analysis_dir : Path
        Contacts analysis directory where outputs are saved.
    enzyme_pdb : Path
        Enzyme PDB path used for SASA filtering.
    contact_results_by_replicate : Mapping[int, Path]
        Mapping from replicate index to contacts JSON path.
    load_contact_result : callable
        Callable used to deserialize contacts results.
    threshold : float
        Relative SASA threshold for exposed residues.
    include_default_aa_groups : bool
        Whether default AA class groups are included.
    custom_protein_groups : dict[str, list[int]] | None
        Optional custom protein groups.
    protein_partitions : dict[str, list[str]] | None
        Optional user-defined protein partitions.
    polymer_type_selections : dict[str, str] | None
        Optional explicit polymer selection map.
    polymer_chain : str
        Polymer chain ID for auto-detection.
    settings_fp : str | None
        Optional settings fingerprint used for cache filenames.

    Returns
    -------
    AggregatedBindingPreferenceResult | None
        Aggregated result, or None if computation failed.
    """
    import MDAnalysis as mda

    from polyzymd.analyses.shared.loader import TrajectoryLoader
    from polyzymd.analyses.shared.surface_exposure import SurfaceExposureFilter

    from ._aggregate import aggregate_binding_preference

    # Compute surface exposure from enzyme structure
    try:
        exposure_filter = SurfaceExposureFilter(threshold=threshold)
        surface_exposure = exposure_filter.calculate(str(enzyme_pdb))
        logger.debug(
            f"Computed surface exposure for {cond.label}: "
            f"{surface_exposure.exposed_count}/{surface_exposure.total_count} residues exposed"
        )
    except (FileNotFoundError, ValueError, OSError, ImportError) as exc:
        logger.warning(f"Failed to compute surface exposure for {cond.label}: {exc}")
        return None

    protein_groups = resolve_protein_groups_from_surface_exposure(
        surface_exposure,
        include_default_aa_groups=include_default_aa_groups,
        custom_protein_groups=custom_protein_groups,
    )
    if not protein_groups:
        logger.warning(f"No protein groups resolved for {cond.label}")
        return None

    polymer_composition: PolymerComposition | None = None
    first_rep = cond.replicates[0] if cond.replicates else 1
    run_dir = sim_config.get_working_directory(first_rep)

    topology_path = None
    try:
        loader = TrajectoryLoader(sim_config)
        topology_path = loader.find_topology(run_dir)
    except (FileNotFoundError, ImportError):
        logger.warning(
            f"Cannot extract polymer composition for {cond.label}: "
            f"no topology file found in {run_dir}"
        )

    if topology_path is not None:
        try:
            universe = mda.Universe(str(topology_path))
            polymer_composition = extract_polymer_composition(
                universe, polymer_type_selections, polymer_chain=polymer_chain
            )
            logger.debug(
                f"Extracted polymer composition for {cond.label}: "
                f"{polymer_composition.total_residues} residues, "
                f"{polymer_composition.total_heavy_atoms} heavy atoms"
            )
        except (ValueError, OSError, AttributeError) as exc:
            logger.warning(f"Failed to extract polymer composition for {cond.label}: {exc}")

    if polymer_composition is None:
        polymer_composition = PolymerComposition()
        logger.warning(
            f"Using empty polymer composition for {cond.label} — enrichment ratios will be NaN"
        )

    rep_results = []
    for rep in cond.replicates:
        contact_path = contact_results_by_replicate.get(rep)
        if contact_path is None:
            logger.warning(f"Contacts file not found for {cond.label} rep{rep} in {analysis_dir}")
            continue

        try:
            contact_result = load_contact_result(contact_path)
            bp_result = compute_binding_preference(
                contact_result=contact_result,
                surface_exposure=surface_exposure,
                protein_groups=protein_groups,
                polymer_composition=polymer_composition,
                protein_partitions=protein_partitions,
            )
            rep_results.append(bp_result)

            if settings_fp is not None:
                rep_bp_path = analysis_dir / f"binding_preference_s{settings_fp}_rep{rep}.json"
            else:
                rep_bp_path = analysis_dir / f"binding_preference_rep{rep}.json"
            bp_result.save(rep_bp_path)
            logger.debug(f"Computed and saved binding preference for {cond.label} rep{rep}")

        except (FileNotFoundError, json.JSONDecodeError, ValueError, KeyError, OSError) as exc:
            logger.warning(f"Failed to compute binding preference for {cond.label} rep{rep}: {exc}")

    if not rep_results:
        logger.warning(f"No binding preference results computed for {cond.label}")
        return None

    agg_result = aggregate_binding_preference(rep_results)
    rep_range = f"{min(cond.replicates)}-{max(cond.replicates)}"
    if settings_fp is not None:
        agg_path = (
            analysis_dir / f"binding_preference_aggregated_s{settings_fp}_reps{rep_range}.json"
        )
    else:
        agg_path = analysis_dir / f"binding_preference_aggregated_reps{rep_range}.json"
    agg_result.save(agg_path)
    logger.info(
        f"Computed binding preference for {cond.label}: "
        f"{len(rep_results)} replicates, {len(protein_groups)} protein groups"
    )

    return agg_result


def compute_binding_preference_from_config(
    contact_result: Any,
    universe: "Universe",
    enzyme_pdb_path: Path | str,
    config: Any,
) -> BindingPreferenceResult:
    """Compute binding preference using contacts plugin settings.

    This is a convenience function that orchestrates the full binding
    preference computation using configuration from analysis.yaml.
    It:
    1. Calculates surface exposure from enzyme PDB
    2. Resolves protein group selections to residue IDs
    3. Computes binding preference with enrichment ratios

    Parameters
    ----------
    contact_result : Any
        Contact analysis results from trajectory
    universe : Universe
        MDAnalysis Universe (used for resolving selections)
    enzyme_pdb_path : Path or str
        Path to enzyme PDB file for SASA calculation
    config : Any
        Settings object with binding preference fields

    Returns
    -------
    BindingPreferenceResult
        Binding preference with enrichment ratios

    Notes
    -----
    This function requires rust_sasa_python to be installed for
    surface exposure calculation. Install with:
        pip install rust-sasa-python

    Examples
    --------
    >>> from polyzymd.analyses.contacts import ContactsSettings
    >>> config = ContactsSettings(
    ...     compute_binding_preference=True,
    ...     surface_exposure_threshold=0.2,
    ... )
    >>> result = compute_binding_preference_from_config(
    ...     contact_result, universe, "enzyme.pdb", config
    ... )
    >>> print(result.enrichment_matrix())
    {'SBM': {'aromatic': 1.45, 'polar': 0.82, ...}, ...}
    """
    from polyzymd.analyses.shared.surface_exposure import SurfaceExposureFilter

    logger.info("Computing binding preference from config...")

    # Step 1: Calculate surface exposure
    exposure_filter = SurfaceExposureFilter(threshold=config.surface_exposure_threshold)
    surface_exposure = exposure_filter.calculate(enzyme_pdb_path)

    logger.info(
        f"Surface exposure: {surface_exposure.exposed_count}/{surface_exposure.total_count} "
        f"residues exposed (threshold={config.surface_exposure_threshold})"
    )

    # Step 2: Resolve protein group selections to residue IDs
    protein_groups = resolve_protein_group_selections(universe, config.protein_group_selections)

    for group_name, resids in protein_groups.items():
        exposed_in_group = surface_exposure.exposed_in_selection(resids)
        logger.debug(f"  {group_name}: {len(resids)} residues, {len(exposed_in_group)} exposed")

    # Step 3: Get polymer types
    # Extract polymer chain ID from config.polymer_selection (e.g. "chainID C" → "C").
    # Falls back to default "C" if the selection doesn't match the "chainID X" pattern.
    _polymer_chain = "C"
    _ps = getattr(config, "polymer_selection", "")
    if _ps.startswith("chainID ") and len(_ps.split()) == 2:
        _polymer_chain = _ps.split()[1]

    polymer_types = resolve_polymer_type_selections(
        universe, config.polymer_type_selections, polymer_chain=_polymer_chain
    )
    logger.info(f"Polymer types for binding preference: {polymer_types}")

    # Step 4: Extract polymer composition for dual normalization
    polymer_composition = extract_polymer_composition(
        universe, config.polymer_type_selections, polymer_chain=_polymer_chain
    )

    # Step 5: Compute binding preference
    result = compute_binding_preference(
        contact_result=contact_result,
        surface_exposure=surface_exposure,
        protein_groups=protein_groups,
        polymer_composition=polymer_composition,
        polymer_types=polymer_types if polymer_types else None,
        protein_group_selections=config.protein_group_selections,
        polymer_type_selections=config.polymer_type_selections,
    )

    return result
