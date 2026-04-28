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


class BindingPreferenceContactIdentityError(RuntimeError):
    """Raised when a contacts artifact does not match the requested BP input."""


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
    contact_settings_fp: str | None = None,
    equilibration: str | None = None,
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
        Optional binding-preference settings fingerprint used for cache filenames.
    contact_settings_fp : str or None, optional
        Upstream contacts settings fingerprint expected for all contacts inputs.
    equilibration : str or None, optional
        Requested contacts equilibration window. Stored in cache metadata so
        downstream consumers can reject cross-window artifacts.

    Returns
    -------
    AggregatedBindingPreferenceResult | None
        Aggregated result, or None if computation failed.
    """
    import MDAnalysis as mda

    from polyzymd.analyses.shared.loader import TrajectoryLoader
    from polyzymd.analyses.shared.paths import format_replicate_cache_token
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
    successful_replicates: list[int] = []
    contact_window_label: str | None = None
    for rep in cond.replicates:
        contact_path = contact_results_by_replicate.get(rep)
        if contact_path is None:
            logger.warning(f"Contacts file not found for {cond.label} rep{rep} in {analysis_dir}")
            continue

        try:
            contact_result = load_contact_result(contact_path)
            _validate_contact_result_contract(
                contact_result,
                contact_path,
                expected_replicate=rep,
                contact_settings_fp=contact_settings_fp,
            )
            result_window = _contact_result_window_label(contact_result)
            if equilibration is not None and not _window_labels_match(result_window, equilibration):
                logger.warning(
                    "Contacts file for %s rep%d has window %s, expected %s: %s",
                    cond.label,
                    rep,
                    result_window or "<missing>",
                    equilibration,
                    contact_path,
                )
                continue
            if contact_window_label is None:
                contact_window_label = result_window or equilibration
            bp_result = compute_binding_preference(
                contact_result=contact_result,
                surface_exposure=surface_exposure,
                protein_groups=protein_groups,
                polymer_composition=polymer_composition,
                protein_partitions=protein_partitions,
            )
            bp_metadata = dict(getattr(bp_result, "metadata", {}) or {})
            if settings_fp is not None:
                bp_metadata["binding_preference_settings_fingerprint"] = settings_fp
                bp_metadata["settings_fingerprint"] = settings_fp
            if contact_settings_fp is not None:
                bp_metadata["contacts_settings_fingerprint"] = contact_settings_fp
            if contact_window_label is not None:
                bp_metadata["equilibration"] = contact_window_label
            bp_metadata["replicate"] = rep
            bp_result.metadata = bp_metadata
            rep_results.append(bp_result)
            successful_replicates.append(rep)

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
    agg_metadata = dict(getattr(agg_result, "metadata", {}) or {})
    if settings_fp is not None:
        agg_metadata["binding_preference_settings_fingerprint"] = settings_fp
        agg_metadata["settings_fingerprint"] = settings_fp
    if contact_settings_fp is not None:
        agg_metadata["contacts_settings_fingerprint"] = contact_settings_fp
    if contact_window_label is not None:
        agg_metadata["equilibration"] = contact_window_label
    agg_metadata["replicates"] = successful_replicates
    agg_result.metadata = agg_metadata
    rep_token = format_replicate_cache_token(successful_replicates)
    if settings_fp is not None:
        agg_path = analysis_dir / f"binding_preference_aggregated_s{settings_fp}_{rep_token}.json"
    else:
        agg_path = analysis_dir / f"binding_preference_aggregated_{rep_token}.json"
    agg_result.save(agg_path)
    logger.info(
        f"Computed binding preference for {cond.label}: "
        f"{len(rep_results)} replicates, {len(protein_groups)} protein groups"
    )

    return agg_result


def _validate_contact_result_contract(
    contact_result: Any,
    contact_path: Path,
    *,
    expected_replicate: int,
    contact_settings_fp: str | None,
) -> None:
    """Validate that a contacts result matches the requested BP input contract.

    Parameters
    ----------
    contact_result : Any
        Loaded contacts result.
    contact_path : Path
        Contacts artifact path used for metadata fallback.
    expected_replicate : int
        Replicate ID requested by the condition workflow.
    contact_settings_fp : str or None
        Expected contacts settings fingerprint, when known.

    Raises
    ------
    BindingPreferenceContactIdentityError
        If replicate or settings metadata conflicts with the requested input.
    """

    stored_replicate = _contact_result_replicate(contact_result)
    if stored_replicate is not None and stored_replicate != int(expected_replicate):
        raise BindingPreferenceContactIdentityError(
            "Contacts artifact replicate mismatch for binding preference: "
            f"stored={stored_replicate}, expected={expected_replicate}, path={contact_path}"
        )

    if contact_settings_fp is None:
        return

    stored_settings_fp = _contact_result_settings_fingerprint(contact_result, contact_path)
    if stored_settings_fp != contact_settings_fp:
        raise BindingPreferenceContactIdentityError(
            "Contacts artifact settings fingerprint mismatch for binding preference: "
            f"stored={stored_settings_fp or '<missing>'}, expected={contact_settings_fp}, "
            f"path={contact_path}"
        )


def _contact_result_replicate(contact_result: Any) -> int | None:
    """Return the replicate declared by a contacts result, if present."""

    candidates = [getattr(contact_result, "replicate", None)]
    metadata = getattr(contact_result, "metadata", None)
    if isinstance(metadata, dict):
        candidates.append(metadata.get("replicate"))
    for candidate in candidates:
        if candidate is None:
            continue
        try:
            return int(candidate)
        except (TypeError, ValueError) as exc:
            raise BindingPreferenceContactIdentityError(
                f"Contacts artifact has invalid replicate metadata: {candidate!r}"
            ) from exc
    return None


def _contact_result_settings_fingerprint(contact_result: Any, contact_path: Path) -> str | None:
    """Return the settings fingerprint declared by a contacts result."""

    metadata = getattr(contact_result, "metadata", None)
    if isinstance(metadata, dict):
        stored = (
            metadata.get("contacts_settings_fingerprint")
            or metadata.get("contact_settings_fingerprint")
            or metadata.get("settings_fingerprint")
            or metadata.get("settings_fp")
        )
        if stored is not None:
            return str(stored)

    stored = getattr(contact_result, "contacts_settings_fingerprint", None) or getattr(
        contact_result,
        "contact_settings_fingerprint",
        None,
    )
    stored = (
        stored
        or getattr(contact_result, "settings_fingerprint", None)
        or getattr(
            contact_result,
            "settings_fp",
            None,
        )
    )
    if stored is not None:
        return str(stored)

    from polyzymd.analyses.contacts._paths import extract_contact_settings_fingerprint

    return extract_contact_settings_fingerprint(contact_path)


def _contact_result_window_label(contact_result: Any) -> str | None:
    """Return a compact equilibration label recorded on a contacts result.

    Parameters
    ----------
    contact_result : Any
        Loaded contacts result.

    Returns
    -------
    str or None
        Window label such as ``"10ns"`` when metadata is available.
    """

    stored_time = getattr(contact_result, "equilibration_time", None)
    stored_unit = getattr(contact_result, "equilibration_unit", None)
    if stored_time is not None and stored_unit is not None:
        try:
            return f"{float(stored_time):g}{stored_unit}"
        except (TypeError, ValueError):
            return None
    metadata = getattr(contact_result, "metadata", {}) or {}
    if isinstance(metadata, dict):
        stored = metadata.get("equilibration") or metadata.get("equilibration_time_label")
        if isinstance(stored, str):
            return stored
    return None


def _window_labels_match(stored: str | None, requested: str) -> bool:
    """Return whether two equilibration labels refer to the same time."""

    if stored is None:
        return False

    import math

    from polyzymd.analyses.shared.loader import convert_time, parse_time_string

    try:
        stored_time, stored_unit = parse_time_string(stored)
        requested_time, requested_unit = parse_time_string(requested)
        stored_ps = convert_time(float(stored_time), stored_unit, "ps")
        requested_ps = convert_time(float(requested_time), requested_unit, "ps")
    except (TypeError, ValueError):
        return False
    return math.isclose(stored_ps, requested_ps, rel_tol=0.0, abs_tol=1.0e-9)


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
