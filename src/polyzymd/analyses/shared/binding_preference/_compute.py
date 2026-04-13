"""Binding preference compute kernels and partition logic."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, Literal

from polyzymd.analyses.shared.aa_classification import DEFAULT_AA_CLASS_SELECTIONS

from ._models import (
    BindingPreferenceEntry,
    BindingPreferenceResult,
    PartitionBindingEntry,
    PartitionBindingResult,
    PartitionCoverageEntry,
    PartitionCoverageResult,
    PolymerBindingPreferenceResult,
    PolymerComposition,
    SystemCoverageResult,
)

if TYPE_CHECKING:
    from polyzymd.analyses.shared.surface_exposure import SurfaceExposureResult

logger = logging.getLogger(__name__)


def _detect_overlapping_groups(groups: dict[str, set[int]]) -> list[tuple[str, str]]:
    """Detect pairs of groups that share residue IDs.

    Parameters
    ----------
    groups : dict[str, set[int]]
        Mapping of group name to residue IDs

    Returns
    -------
    list[tuple[str, str]]
        List of (group_a, group_b) pairs that have overlapping residue IDs
    """
    group_names = list(groups.keys())
    overlaps = []

    for i, g1 in enumerate(group_names):
        for g2 in group_names[i + 1 :]:
            if groups[g1] & groups[g2]:  # Non-empty intersection
                overlaps.append((g1, g2))

    return overlaps


def _compute_enrichment(coverage_share: float, expected_share: float) -> float | None:
    """Calculate zero-centered enrichment: (observed / expected) - 1.

    Parameters
    ----------
    coverage_share : float
        Fraction of contacts to this element
    expected_share : float
        Expected fraction based on surface availability

    Returns
    -------
    float | None
        Enrichment value, or None if expected_share is 0
    """
    if expected_share > 0 and coverage_share > 0:
        return (coverage_share / expected_share) - 1
    elif expected_share > 0 and coverage_share == 0:
        return -1.0  # Complete avoidance
    else:
        return None  # Cannot compute (no expected share)


def _compute_partition_binding(
    partition_name: str,
    partition_type: Literal["aa_class", "user_defined"],
    partition_groups: dict[str, set[int]],
    exposed_partition: dict[str, set[int]],
    polymer_type: str,
    contact_data_for_polymer: dict[str, dict[str, Any]],
    total_exposed: int,
    total_contact_frames_for_polymer: int,
) -> "PartitionBindingResult":
    """Compute binding preference for a single partition and single polymer type.

    This is the per-polymer version of _compute_partition_coverage.
    It computes contact_share and enrichment for one polymer type across
    all elements in a partition.

    Parameters
    ----------
    partition_name : str
        Name of the partition (e.g., "aa_class", "lid_helices")
    partition_type : str
        One of: "aa_class", "user_defined"
    partition_groups : dict[str, set[int]]
        Mapping of element name to ALL residue IDs in that element
    exposed_partition : dict[str, set[int]]
        Mapping of element name to EXPOSED residue IDs only
    polymer_type : str
        Polymer type this binding is for (e.g., "SBM")
    contact_data_for_polymer : dict[str, dict[str, Any]]
        Contact data for this polymer: {element_name: {"total_frames": int, "residues_contacted": set}}
    total_exposed : int
        Total number of exposed residues across all elements in the partition
    total_contact_frames_for_polymer : int
        Total contact frames for this polymer type (for contact_share denominator)

    Returns
    -------
    PartitionBindingResult
        Binding result for this polymer type on this partition
    """
    # The classes below are defined earlier in this same module -- no import needed.
    # (Previously used a self-import to avoid forward references.)

    # Build partition entries
    binding_entries = []
    total_contact_share = 0.0
    total_expected_share = 0.0

    for element_name in sorted(partition_groups.keys()):
        n_total = len(partition_groups.get(element_name, set()))
        n_exposed = len(exposed_partition.get(element_name, set()))

        # Calculate expected share based on surface availability
        expected_share = n_exposed / total_exposed if total_exposed > 0 else 0.0
        total_expected_share += expected_share

        # Get contact data for this element
        edata = contact_data_for_polymer.get(
            element_name, {"total_frames": 0, "residues_contacted": set()}
        )
        total_frames = edata.get("total_frames", 0)
        residues_contacted = edata.get("residues_contacted", set())
        n_residues_contacted = len(residues_contacted)

        # Calculate contact share (fraction of this polymer's contacts to this element)
        if total_contact_frames_for_polymer > 0:
            contact_share = total_frames / total_contact_frames_for_polymer
        else:
            contact_share = 0.0
        total_contact_share += contact_share

        # Calculate enrichment
        enrichment = _compute_enrichment(contact_share, expected_share)

        binding_entries.append(
            PartitionBindingEntry(
                partition_element=element_name,
                polymer_type=polymer_type,
                total_contact_frames=total_frames,
                contact_share=contact_share,
                expected_share=expected_share,
                enrichment=enrichment,
                n_exposed_in_element=n_exposed,
                n_residues_in_element=n_total,
                n_residues_contacted=n_residues_contacted,
            )
        )

    return PartitionBindingResult(
        partition_name=partition_name,
        partition_type=partition_type,
        polymer_type=polymer_type,
        entries=binding_entries,
        total_contact_share=total_contact_share,
        total_expected_share=total_expected_share,
        total_contact_frames=total_contact_frames_for_polymer,
    )


def _compute_polymer_binding_preference(
    contact_data: dict[str, dict[str, dict[str, Any]]],
    total_contacts_by_polymer: dict[str, int],
    protein_groups: dict[str, set[int]],
    exposed_groups: dict[str, set[int]],
    protein_partitions: dict[str, list[str]] | None,
    total_exposed: int,
    n_frames: int,
    surface_exposure_threshold: float | None,
    polymer_composition: "PolymerComposition | None",
    protein_groups_used: dict[str, str] | None,
) -> "PolymerBindingPreferenceResult":
    """Compute per-polymer binding preference using partition-based analysis.

    Parameters
    ----------
    contact_data : dict[str, dict[str, dict[str, Any]]]
        Contact data: {polymer_type: {group_name: {"total_frames": int, "residues_contacted": set}}}
    total_contacts_by_polymer : dict[str, int]
        Total contact frames per polymer type
    protein_groups : dict[str, set[int]]
        All protein groups (AA classes + custom)
    exposed_groups : dict[str, set[int]]
        Exposed residues per group
    protein_partitions : dict[str, list[str]] | None
        User-defined partitions
    total_exposed : int
        Total exposed residues
    n_frames : int
        Number of frames analyzed
    surface_exposure_threshold : float | None
        SASA threshold used
    polymer_composition : PolymerComposition | None
        Polymer composition metadata
    protein_groups_used : dict[str, str] | None
        Selection strings used (for metadata)

    Returns
    -------
    PolymerBindingPreferenceResult
        Per-polymer partition-based binding preference
    """
    # PolymerBindingPreferenceResult is defined later in this same module.
    # (Previously used a self-import to avoid forward references.)

    all_polymer_types = sorted(contact_data.keys())

    # Separate AA class groups from custom groups
    aa_class_names = set(DEFAULT_AA_CLASS_SELECTIONS.keys())
    aa_class_groups: dict[str, set[int]] = {}
    aa_class_exposed: dict[str, set[int]] = {}

    for group_name, resids in protein_groups.items():
        if group_name in aa_class_names:
            aa_class_groups[group_name] = resids
            aa_class_exposed[group_name] = exposed_groups.get(group_name, set())

    # Get all exposed residue IDs (for computing "rest_of_protein")
    all_exposed_resids: set[int] = set()
    for resids in exposed_groups.values():
        all_exposed_resids.update(resids)

    # Compute total exposed for AA class partition
    aa_class_total_exposed = sum(len(resids) for resids in aa_class_exposed.values())

    # ---------------------------------------------------------------------
    # 1. Compute AA Class Partition for each polymer type
    # ---------------------------------------------------------------------
    aa_class_binding: dict[str, "PartitionBindingResult"] = {}

    for poly_type in all_polymer_types:
        poly_contact_data = contact_data.get(poly_type, {})
        total_frames = total_contacts_by_polymer.get(poly_type, 0)

        # Filter contact data to AA class groups only
        aa_contact_data: dict[str, dict[str, Any]] = {}
        for group_name in aa_class_groups:
            if group_name in poly_contact_data:
                aa_contact_data[group_name] = poly_contact_data[group_name]
            else:
                aa_contact_data[group_name] = {"total_frames": 0, "residues_contacted": set()}

        # Compute total contact frames for AA class groups only
        aa_total_frames = sum(d.get("total_frames", 0) for d in aa_contact_data.values())

        aa_class_binding[poly_type] = _compute_partition_binding(
            partition_name="aa_class",
            partition_type="aa_class",
            partition_groups=aa_class_groups,
            exposed_partition=aa_class_exposed,
            polymer_type=poly_type,
            contact_data_for_polymer=aa_contact_data,
            total_exposed=aa_class_total_exposed,
            total_contact_frames_for_polymer=aa_total_frames,
        )

    # ---------------------------------------------------------------------
    # 2. Compute User-Defined Partitions for each polymer type
    # ---------------------------------------------------------------------
    user_defined_partitions: dict[str, dict[str, "PartitionBindingResult"]] = {}

    if protein_partitions:
        for partition_name, group_names in protein_partitions.items():
            # Build partition groups from referenced group names
            partition_groups_map: dict[str, set[int]] = {}
            partition_exposed_map: dict[str, set[int]] = {}
            all_partition_exposed: set[int] = set()

            for group_name in group_names:
                if group_name not in protein_groups:
                    logger.warning(
                        f"Partition '{partition_name}' references undefined group "
                        f"'{group_name}' - skipping this group"
                    )
                    continue

                partition_groups_map[group_name] = protein_groups[group_name]
                partition_exposed_map[group_name] = exposed_groups.get(group_name, set())
                all_partition_exposed.update(partition_exposed_map[group_name])

            if not partition_groups_map:
                logger.warning(f"Partition '{partition_name}' has no valid groups - skipping")
                continue

            # Check if we need to add rest_of_protein
            rest_exposed = all_exposed_resids - all_partition_exposed
            if rest_exposed:
                # Partition doesn't cover all residues - add rest_of_protein
                rest_all = set().union(*exposed_groups.values()) - all_partition_exposed
                partition_groups_map["rest_of_protein"] = rest_all
                partition_exposed_map["rest_of_protein"] = rest_exposed
                logger.debug(
                    f"Partition '{partition_name}': auto-added 'rest_of_protein' "
                    f"with {len(rest_exposed)} exposed residues"
                )

            # Compute total exposed for this partition
            user_partition_total_exposed = sum(len(r) for r in partition_exposed_map.values())

            # Compute for each polymer type
            user_defined_partitions[partition_name] = {}
            for poly_type in all_polymer_types:
                poly_contact_data = contact_data.get(poly_type, {})

                # Build contact data for this partition
                partition_contact_data: dict[str, dict[str, Any]] = {}
                for group_name in partition_groups_map:
                    if group_name == "rest_of_protein":
                        # Aggregate contacts from groups NOT in this partition
                        rest_frames = 0
                        rest_residues: set[int] = set()
                        for gname, gdata in poly_contact_data.items():
                            if gname not in group_names:
                                rest_frames += gdata.get("total_frames", 0)
                                rest_residues.update(gdata.get("residues_contacted", set()))
                        partition_contact_data["rest_of_protein"] = {
                            "total_frames": rest_frames,
                            "residues_contacted": rest_residues,
                        }
                    elif group_name in poly_contact_data:
                        partition_contact_data[group_name] = poly_contact_data[group_name]
                    else:
                        partition_contact_data[group_name] = {
                            "total_frames": 0,
                            "residues_contacted": set(),
                        }

                # Compute total frames for this partition
                partition_total_frames = sum(
                    d.get("total_frames", 0) for d in partition_contact_data.values()
                )

                user_defined_partitions[partition_name][poly_type] = _compute_partition_binding(
                    partition_name=partition_name,
                    partition_type="user_defined",
                    partition_groups=partition_groups_map,
                    exposed_partition=partition_exposed_map,
                    polymer_type=poly_type,
                    contact_data_for_polymer=partition_contact_data,
                    total_exposed=user_partition_total_exposed,
                    total_contact_frames_for_polymer=partition_total_frames,
                )

            logger.info(
                f"Computed user-defined partition '{partition_name}' binding for "
                f"{len(all_polymer_types)} polymer types"
            )

    return PolymerBindingPreferenceResult(
        aa_class_binding=aa_class_binding,
        user_defined_partitions=user_defined_partitions,
        n_frames=n_frames,
        total_exposed_residues=total_exposed,
        surface_exposure_threshold=surface_exposure_threshold,
        polymer_types=all_polymer_types,
        polymer_composition=polymer_composition,
        protein_groups_used=protein_groups_used or {},
    )


def _compute_partition_coverage(
    partition_name: str,
    partition_type: Literal["aa_class", "binary_custom", "combined_custom", "user_defined"],
    partition_groups: dict[str, set[int]],
    exposed_partition: dict[str, set[int]],
    entries: list["BindingPreferenceEntry"],
    total_exposed: int,
    all_polymer_types: list[str],
) -> PartitionCoverageResult:
    """Compute coverage for a single partition.

    A partition divides the protein surface into mutually exclusive elements.
    This function validates the partition and computes coverage metrics.

    Parameters
    ----------
    partition_name : str
        Name of the partition (e.g., "aa_class", "lid_helix_5_vs_rest")
    partition_type : str
        One of: "aa_class", "binary_custom", "combined_custom", "user_defined"
    partition_groups : dict[str, set[int]]
        Mapping of element name to ALL residue IDs in that element
    exposed_partition : dict[str, set[int]]
        Mapping of element name to EXPOSED residue IDs only
    entries : list[BindingPreferenceEntry]
        Binding preference entries to aggregate
    total_exposed : int
        Total number of exposed residues across all elements
    all_polymer_types : list[str]
        List of all polymer types

    Returns
    -------
    PartitionCoverageResult
        Coverage result for this partition
    """
    # Build a mapping: resid -> partition_element
    resid_to_element: dict[int, str] = {}
    for element_name, resids in partition_groups.items():
        for resid in resids:
            resid_to_element[resid] = element_name

    # Aggregate contact frames by partition element
    # Structure: {element: {"total_frames": int, "by_polymer": {poly: frames}}}
    element_totals: dict[str, dict[str, Any]] = {
        element_name: {"total_frames": 0, "by_polymer": dict.fromkeys(all_polymer_types, 0)}
        for element_name in partition_groups.keys()
    }

    for entry in entries:
        # Determine which element this entry's protein_group maps to
        element_name = entry.protein_group
        if element_name in element_totals:
            element_totals[element_name]["total_frames"] += entry.total_contact_frames
            element_totals[element_name]["by_polymer"][entry.polymer_type] = (
                element_totals[element_name]["by_polymer"].get(entry.polymer_type, 0)
                + entry.total_contact_frames
            )

    # Calculate grand total of contacts across all elements
    grand_total = sum(et["total_frames"] for et in element_totals.values())

    # Build partition entries
    coverage_entries = []
    total_coverage_share = 0.0
    total_expected_share = 0.0

    for element_name in sorted(partition_groups.keys()):
        n_total = len(partition_groups.get(element_name, set()))
        n_exposed = len(exposed_partition.get(element_name, set()))

        # Calculate expected share
        expected_share = n_exposed / total_exposed if total_exposed > 0 else 0.0
        total_expected_share += expected_share

        # Get contact data
        edata = element_totals.get(element_name, {"total_frames": 0, "by_polymer": {}})
        total_frames = edata["total_frames"]
        by_polymer = edata.get("by_polymer", {})

        # Calculate coverage share
        coverage_share = total_frames / grand_total if grand_total > 0 else 0.0
        total_coverage_share += coverage_share

        # Calculate enrichment
        enrichment = _compute_enrichment(coverage_share, expected_share)

        # Calculate polymer contributions
        polymer_contributions: dict[str, float] = {}
        if total_frames > 0:
            for pt, pf in by_polymer.items():
                polymer_contributions[pt] = pf / total_frames
        else:
            for pt in all_polymer_types:
                polymer_contributions[pt] = 0.0

        coverage_entries.append(
            PartitionCoverageEntry(
                partition_element=element_name,
                total_contact_frames=total_frames,
                coverage_share=coverage_share,
                expected_share=expected_share,
                coverage_enrichment=enrichment,
                n_exposed_in_element=n_exposed,
                n_residues_in_element=n_total,
                polymer_contributions=polymer_contributions,
            )
        )

    return PartitionCoverageResult(
        partition_name=partition_name,
        partition_type=partition_type,
        entries=coverage_entries,
        total_coverage_share=total_coverage_share,
        total_expected_share=total_expected_share,
    )


def _compute_system_coverage(
    entries: list[BindingPreferenceEntry],
    protein_groups: dict[str, set[int]],
    exposed_groups: dict[str, set[int]],
    total_exposed: int,
    n_frames: int,
    surface_exposure_threshold: float | None,
    protein_group_selections: dict[str, str] | None,
    protein_partitions: dict[str, list[str]] | None = None,
) -> SystemCoverageResult:
    """Compute system-level coverage using partition-based analysis.

    This function computes coverage enrichments using proper partitions to
    avoid the overlap bug where custom groups and AA classes can inflate
    the expected_share denominator.

    Partition Strategy
    ------------------
    1. **AA Class Partition**: 5-way partition by amino acid class.
       Always computed, uses only the 5 default AA classes.

    2. **Binary Custom Partitions**: For each custom group, compute a
       binary partition (group vs rest_of_protein).

    3. **Combined Custom Partition**: If custom groups don't overlap,
       combine them all + rest_of_protein into a single partition.

    4. **User-Defined Partitions**: Custom partitions from config.
       Each references groups from protein_groups and must be mutually
       exclusive. 'rest_of_protein' is auto-added if needed.

    Parameters
    ----------
    entries : list[BindingPreferenceEntry]
        Binding preference entries (per polymer type × protein group)
    protein_groups : dict[str, set[int]]
        Mapping of group name to ALL residue IDs in that group
    exposed_groups : dict[str, set[int]]
        Mapping of group name to EXPOSED residue IDs only
    total_exposed : int
        Total number of exposed residues
    n_frames : int
        Number of trajectory frames analyzed
    surface_exposure_threshold : float | None
        SASA threshold used for surface filtering
    protein_group_selections : dict[str, str] | None
        Original MDAnalysis selections (for metadata)
    protein_partitions : dict[str, list[str]] | None
        User-defined partitions: {partition_name: [group1, group2, ...]}
        Groups must exist in protein_groups.

    Returns
    -------
    SystemCoverageResult
        Partition-based coverage metrics (schema v2)
    """
    from polyzymd.analyses.shared.aa_classification import DEFAULT_AA_CLASS_SELECTIONS

    # Collect all polymer types
    all_polymer_types = sorted({e.polymer_type for e in entries})

    # Separate AA class groups from custom groups
    aa_class_names = set(DEFAULT_AA_CLASS_SELECTIONS.keys())
    aa_class_groups: dict[str, set[int]] = {}
    aa_class_exposed: dict[str, set[int]] = {}
    custom_groups: dict[str, set[int]] = {}
    custom_exposed: dict[str, set[int]] = {}
    custom_selections: dict[str, str] = {}

    for group_name, resids in protein_groups.items():
        if group_name in aa_class_names:
            aa_class_groups[group_name] = resids
            aa_class_exposed[group_name] = exposed_groups.get(group_name, set())
        else:
            custom_groups[group_name] = resids
            custom_exposed[group_name] = exposed_groups.get(group_name, set())
            if protein_group_selections and group_name in protein_group_selections:
                custom_selections[group_name] = protein_group_selections[group_name]

    # Get all exposed residue IDs (for computing "rest_of_protein")
    all_exposed_resids: set[int] = set()
    for resids in exposed_groups.values():
        all_exposed_resids.update(resids)

    # Filter entries by group type for proper partition computation
    aa_class_entries = [e for e in entries if e.protein_group in aa_class_names]
    custom_entries = [e for e in entries if e.protein_group not in aa_class_names]

    # ---------------------------------------------------------------------
    # 1. Compute AA Class Partition (always)
    # ---------------------------------------------------------------------
    # Compute total exposed for AA class partition
    aa_class_total_exposed = sum(len(resids) for resids in aa_class_exposed.values())

    aa_class_coverage = _compute_partition_coverage(
        partition_name="aa_class",
        partition_type="aa_class",
        partition_groups=aa_class_groups,
        exposed_partition=aa_class_exposed,
        entries=aa_class_entries,
        total_exposed=aa_class_total_exposed,
        all_polymer_types=all_polymer_types,
    )

    # ---------------------------------------------------------------------
    # 2. Compute Binary Custom Partitions (one per custom group)
    # ---------------------------------------------------------------------
    custom_group_coverages: dict[str, PartitionCoverageResult] = {}

    for group_name, group_resids in custom_groups.items():
        group_exposed = custom_exposed.get(group_name, set())

        # Compute "rest_of_protein" as all exposed residues NOT in this group
        rest_exposed = all_exposed_resids - group_exposed
        rest_all = set()
        for gname, gresids in protein_groups.items():
            if gname != group_name:
                rest_all.update(gresids)
        # Actually, rest_all should be all residues NOT in this custom group
        # We need the full protein residue set - but we only have groups
        # For now, use the exposed residues as the proxy

        # Build binary partition
        binary_partition_groups = {
            group_name: group_resids,
            "rest_of_protein": rest_all - group_resids,
        }
        binary_partition_exposed = {
            group_name: group_exposed,
            "rest_of_protein": rest_exposed,
        }

        # Create synthetic entries for "rest_of_protein" by aggregating all other entries
        # We need to compute the contact frames to "rest_of_protein"
        rest_contact_frames: dict[str, int] = dict.fromkeys(all_polymer_types, 0)
        group_contact_frames: dict[str, int] = dict.fromkeys(all_polymer_types, 0)

        for entry in entries:
            if entry.protein_group == group_name:
                group_contact_frames[entry.polymer_type] = entry.total_contact_frames
            elif entry.protein_group not in custom_groups:
                # It's an AA class group - add to rest
                rest_contact_frames[entry.polymer_type] = (
                    rest_contact_frames.get(entry.polymer_type, 0) + entry.total_contact_frames
                )
            # Note: Other custom groups are NOT added to rest - they'll have their own partition

        # Build synthetic entries for the binary partition
        binary_entries: list[BindingPreferenceEntry] = []

        # Add entries for the custom group
        for entry in entries:
            if entry.protein_group == group_name:
                binary_entries.append(entry)

        # Create synthetic entries for rest_of_protein
        for pt, frames in rest_contact_frames.items():
            binary_entries.append(
                BindingPreferenceEntry(
                    polymer_type=pt,
                    protein_group="rest_of_protein",
                    total_contact_frames=frames,
                    mean_contact_fraction=0.0,  # Not used for partition coverage
                    n_residues_contacted=0,  # Not used for partition coverage
                    contact_share=0.0,
                    expected_share=0.0,
                    enrichment=None,
                    n_exposed_in_group=len(rest_exposed),
                    n_residues_in_group=len(binary_partition_groups["rest_of_protein"]),
                )
            )

        binary_total_exposed = len(group_exposed) + len(rest_exposed)

        binary_coverage = _compute_partition_coverage(
            partition_name=f"{group_name}_vs_rest",
            partition_type="binary_custom",
            partition_groups=binary_partition_groups,
            exposed_partition=binary_partition_exposed,
            entries=binary_entries,
            total_exposed=binary_total_exposed,
            all_polymer_types=all_polymer_types,
        )

        custom_group_coverages[group_name] = binary_coverage

    # ---------------------------------------------------------------------
    # 3. Check for overlaps among custom groups
    # ---------------------------------------------------------------------
    overlapping_pairs = _detect_overlapping_groups(custom_exposed)
    has_overlaps = len(overlapping_pairs) > 0

    if has_overlaps:
        logger.warning(
            f"Custom protein groups have overlapping residues: {overlapping_pairs}. "
            f"Combined custom partition will not be computed."
        )

    # ---------------------------------------------------------------------
    # 4. Compute Combined Custom Partition (if no overlaps)
    # ---------------------------------------------------------------------
    combined_custom_coverage: PartitionCoverageResult | None = None

    if custom_groups and not has_overlaps:
        # Build combined partition: all custom groups + rest_of_protein
        combined_partition_groups: dict[str, set[int]] = dict(custom_groups)
        combined_partition_exposed: dict[str, set[int]] = dict(custom_exposed)

        # Compute rest_of_protein for combined partition
        all_custom_exposed: set[int] = set()
        all_custom_resids: set[int] = set()
        for resids in custom_exposed.values():
            all_custom_exposed.update(resids)
        for resids in custom_groups.values():
            all_custom_resids.update(resids)

        rest_exposed_combined = all_exposed_resids - all_custom_exposed
        rest_all_combined: set[int] = set()
        for gname, gresids in protein_groups.items():
            if gname not in custom_groups:
                rest_all_combined.update(gresids)
        rest_all_combined = rest_all_combined - all_custom_resids

        combined_partition_groups["rest_of_protein"] = rest_all_combined
        combined_partition_exposed["rest_of_protein"] = rest_exposed_combined

        # Build entries for combined partition
        combined_entries: list[BindingPreferenceEntry] = list(custom_entries)

        # Add synthetic entries for rest_of_protein
        rest_contact_frames_combined: dict[str, int] = dict.fromkeys(all_polymer_types, 0)
        for entry in aa_class_entries:
            rest_contact_frames_combined[entry.polymer_type] = (
                rest_contact_frames_combined.get(entry.polymer_type, 0) + entry.total_contact_frames
            )

        for pt, frames in rest_contact_frames_combined.items():
            combined_entries.append(
                BindingPreferenceEntry(
                    polymer_type=pt,
                    protein_group="rest_of_protein",
                    total_contact_frames=frames,
                    mean_contact_fraction=0.0,  # Not used for partition coverage
                    n_residues_contacted=0,  # Not used for partition coverage
                    contact_share=0.0,
                    expected_share=0.0,
                    enrichment=None,
                    n_exposed_in_group=len(rest_exposed_combined),
                    n_residues_in_group=len(rest_all_combined),
                )
            )

        combined_total_exposed = sum(len(r) for r in combined_partition_exposed.values())

        combined_custom_coverage = _compute_partition_coverage(
            partition_name="combined_custom",
            partition_type="combined_custom",
            partition_groups=combined_partition_groups,
            exposed_partition=combined_partition_exposed,
            entries=combined_entries,
            total_exposed=combined_total_exposed,
            all_polymer_types=all_polymer_types,
        )

    # ---------------------------------------------------------------------
    # 5. Compute User-Defined Partitions (from protein_partitions config)
    # ---------------------------------------------------------------------
    user_defined_partitions: dict[str, PartitionCoverageResult] = {}

    if protein_partitions:
        for partition_name, group_names in protein_partitions.items():
            # Build partition groups from referenced group names
            partition_groups_map: dict[str, set[int]] = {}
            partition_exposed_map: dict[str, set[int]] = {}
            partition_entries_list: list[BindingPreferenceEntry] = []

            # Collect residues from all specified groups
            all_partition_exposed: set[int] = set()

            for group_name in group_names:
                if group_name not in protein_groups:
                    logger.warning(
                        f"Partition '{partition_name}' references undefined group "
                        f"'{group_name}' - skipping this group"
                    )
                    continue

                partition_groups_map[group_name] = protein_groups[group_name]
                partition_exposed_map[group_name] = exposed_groups.get(group_name, set())
                all_partition_exposed.update(partition_exposed_map[group_name])

                # Add entries for this group
                for entry in entries:
                    if entry.protein_group == group_name:
                        partition_entries_list.append(entry)

            if not partition_groups_map:
                logger.warning(f"Partition '{partition_name}' has no valid groups - skipping")
                continue

            # Check if we need to add rest_of_protein
            # (if partition doesn't cover all exposed residues)
            rest_exposed = all_exposed_resids - all_partition_exposed
            if rest_exposed:
                # Partition doesn't cover all residues - add rest_of_protein
                partition_groups_map["rest_of_protein"] = (
                    set().union(*exposed_groups.values()) - all_partition_exposed
                )
                partition_exposed_map["rest_of_protein"] = rest_exposed

                # Create synthetic entries for rest_of_protein
                rest_contact_frames_user: dict[str, int] = dict.fromkeys(all_polymer_types, 0)
                for entry in entries:
                    if entry.protein_group not in group_names:
                        rest_contact_frames_user[entry.polymer_type] = (
                            rest_contact_frames_user.get(entry.polymer_type, 0)
                            + entry.total_contact_frames
                        )

                for pt, frames in rest_contact_frames_user.items():
                    partition_entries_list.append(
                        BindingPreferenceEntry(
                            polymer_type=pt,
                            protein_group="rest_of_protein",
                            total_contact_frames=frames,
                            mean_contact_fraction=0.0,
                            n_residues_contacted=0,
                            contact_share=0.0,
                            expected_share=0.0,
                            enrichment=None,
                            n_exposed_in_group=len(rest_exposed),
                            n_residues_in_group=len(partition_groups_map["rest_of_protein"]),
                        )
                    )

                logger.debug(
                    f"Partition '{partition_name}': auto-added 'rest_of_protein' "
                    f"with {len(rest_exposed)} exposed residues"
                )
            else:
                logger.debug(
                    f"Partition '{partition_name}': covers all exposed residues, "
                    f"no 'rest_of_protein' needed"
                )

            # Compute total exposed for this partition
            user_partition_total_exposed = sum(len(r) for r in partition_exposed_map.values())

            # Compute the partition coverage
            user_partition_coverage = _compute_partition_coverage(
                partition_name=partition_name,
                partition_type="user_defined",
                partition_groups=partition_groups_map,
                exposed_partition=partition_exposed_map,
                entries=partition_entries_list,
                total_exposed=user_partition_total_exposed,
                all_polymer_types=all_polymer_types,
            )

            user_defined_partitions[partition_name] = user_partition_coverage
            logger.info(
                f"Computed user-defined partition '{partition_name}' with "
                f"{len(partition_groups_map)} groups"
            )

    # ---------------------------------------------------------------------
    # 6. Compute total contact frames
    # ---------------------------------------------------------------------
    total_contact_frames = sum(e.total_contact_frames for e in entries)

    return SystemCoverageResult(
        aa_class_coverage=aa_class_coverage,
        custom_group_coverages=custom_group_coverages,
        combined_custom_coverage=combined_custom_coverage,
        user_defined_partitions=user_defined_partitions,
        n_frames=n_frames,
        total_contact_frames=total_contact_frames,
        total_exposed_residues=total_exposed,
        surface_exposure_threshold=surface_exposure_threshold,
        custom_group_selections=custom_selections,
        polymer_types_included=all_polymer_types,
        has_overlapping_custom_groups=has_overlaps,
        overlapping_group_pairs=overlapping_pairs,
    )


def compute_binding_preference(
    contact_result: Any,
    surface_exposure: "SurfaceExposureResult",
    protein_groups: dict[str, set[int]],
    polymer_composition: PolymerComposition,
    polymer_types: list[str] | None = None,
    protein_group_selections: dict[str, str] | None = None,
    polymer_type_selections: dict[str, str] | None = None,
    protein_partitions: dict[str, list[str]] | None = None,
) -> BindingPreferenceResult:
    """Compute binding preference from contact results.

    This function computes enrichment ratios for each (polymer_type, protein_group)
    combination, answering: "Does this polymer type preferentially bind this
    protein group compared to random chance?"

    The enrichment calculation accounts for:
    1. Surface exposure (only exposed residues are considered)
    2. Contact duration (contact frames are summed, not binary counts)
    3. Polymer composition (normalization by residue count or heavy atom count)

    Parameters
    ----------
    contact_result : Any
        Raw contact analysis results from trajectory
    surface_exposure : SurfaceExposureResult
        Surface exposure data for filtering buried residues
    protein_groups : dict[str, set[int]]
        Mapping of group name to residue IDs in that group.
        Only surface-exposed residues in each group are counted.
        Example: {"aromatic": {12, 45, 67}, "charged": {23, 34}}
    polymer_composition : PolymerComposition
        Polymer composition data (residue and heavy atom counts per type).
        Used for dual normalization of enrichment ratios.
    polymer_types : list[str], optional
        If provided, only compute for these polymer residue types.
        If None, all polymer types found in contacts are used.
    protein_group_selections : dict[str, str], optional
        Original MDAnalysis selections (for metadata/reproducibility)
    polymer_type_selections : dict[str, str], optional
        Original MDAnalysis selections (for metadata/reproducibility)
    protein_partitions : dict[str, list[str]], optional
        User-defined partitions for system coverage plots.
        Each partition maps a name to a list of group names from protein_groups.
        Groups within a partition must be mutually exclusive.
        Example: {"lid_helices": ["lid_helix_5", "lid_helix_10"]}

    Returns
    -------
    BindingPreferenceResult
        Binding preference metrics with dual enrichment normalization

    Notes
    -----
    Enrichment calculation is centered at zero

    For each (polymer_type, protein_group) pair:

        contact_share = polymer_contacts_to_group / polymer_total_contacts

    Two normalization methods are computed:

    1. **Residue-based** (matches experimental concentration ratios):

        expected_by_residue = polymer_residue_count / total_polymer_residues
        enrichment_by_residue = (contact_share / expected_by_residue) - 1

    2. **Atom-based** (accounts for monomer size via heavy atoms):

        expected_by_atoms = polymer_heavy_atoms / total_polymer_heavy_atoms
        enrichment_by_atoms = (contact_share / expected_by_atoms) - 1

    Interpretation (both methods):
    - enrichment > 0: Preferential binding (more contacts than expected)
    - enrichment = 0: Neutral (matches random chance)
    - enrichment < 0: Avoidance (fewer contacts than expected)
    - enrichment = -1: Complete avoidance (no contacts at all)

    When to use which metric:
    - enrichment_by_residue: Direct comparison to experimental concentration ratios
    - enrichment_by_atoms: Reveals true chemical affinity vs. geometric/steric effects
    """
    exposed_resids = surface_exposure.exposed_resids
    n_frames = contact_result.n_frames
    total_exposed = len(exposed_resids)

    logger.info(
        f"Computing binding preference: {total_exposed} exposed residues, "
        f"{len(protein_groups)} protein groups"
    )
    logger.info(
        f"Polymer composition: {polymer_composition.residue_counts} residues, "
        f"{polymer_composition.heavy_atom_counts} heavy atoms"
    )

    # Filter protein_groups to only exposed residues
    exposed_groups: dict[str, set[int]] = {}
    for group_name, resids in protein_groups.items():
        exposed_groups[group_name] = resids & exposed_resids
        n_orig = len(resids)
        n_exposed = len(exposed_groups[group_name])
        logger.debug(f"  {group_name}: {n_exposed}/{n_orig} residues are exposed")

    # Collect contact frame counts per (polymer_type, protein_group)
    # Structure: {polymer_type: {protein_group: {"frames": int, "residues": set}}}
    contact_data: dict[str, dict[str, dict[str, Any]]] = {}

    # Also track total contacts per polymer type (for contact_share denominator)
    total_contacts_by_polymer: dict[str, int] = {}

    for rc in contact_result.residue_contacts:
        resid = rc.protein_resid
        if resid not in exposed_resids:
            continue  # Skip buried residues

        # Determine which protein groups this residue belongs to
        residue_groups = [gname for gname, gresids in exposed_groups.items() if resid in gresids]

        if not residue_groups:
            # Residue is exposed but not in any defined group
            continue

        # Get contacts by polymer type for this residue
        contacts_by_type = rc.contacts_by_polymer_type(n_frames)

        for poly_type, frac in contacts_by_type.items():
            if polymer_types and poly_type not in polymer_types:
                continue

            # Convert fraction to frame count
            contact_frames = int(round(frac * n_frames))

            if contact_frames == 0:
                continue

            # Initialize polymer type if needed
            if poly_type not in contact_data:
                contact_data[poly_type] = {}
                total_contacts_by_polymer[poly_type] = 0

            # Add to total contacts for this polymer
            total_contacts_by_polymer[poly_type] += contact_frames

            # Add to each group this residue belongs to
            for gname in residue_groups:
                if gname not in contact_data[poly_type]:
                    contact_data[poly_type][gname] = {
                        "total_frames": 0,
                        "residues_contacted": set(),
                    }
                contact_data[poly_type][gname]["total_frames"] += contact_frames
                contact_data[poly_type][gname]["residues_contacted"].add(resid)

    # Get polymer composition totals (for metadata/secondary analysis)
    total_poly_residues = polymer_composition.total_residues
    total_poly_atoms = polymer_composition.total_heavy_atoms

    # Build result entries with enrichment calculations
    entries = []

    for poly_type in sorted(contact_data.keys()):
        total_poly_contacts = total_contacts_by_polymer.get(poly_type, 0)

        # Get polymer composition for this type (stored as metadata)
        poly_res_count = polymer_composition.residue_counts.get(poly_type, 0)
        poly_atom_count = polymer_composition.heavy_atom_counts.get(poly_type, 0)

        for group_name in sorted(protein_groups.keys()):
            n_total_in_group = len(protein_groups.get(group_name, set()))
            n_exposed_in_group = len(exposed_groups.get(group_name, set()))

            # Calculate expected share based on PROTEIN SURFACE AVAILABILITY
            # This is the correct normalization: how much of the exposed surface
            # is this protein group?
            if total_exposed > 0:
                expected_share = n_exposed_in_group / total_exposed
            else:
                expected_share = 0.0

            # Get contact data for this (polymer, group) pair
            gdata = contact_data[poly_type].get(group_name, {})
            contact_frames = gdata.get("total_frames", 0)
            residues_contacted = gdata.get("residues_contacted", set())
            n_residues_contacted = len(residues_contacted)

            # Calculate mean contact fraction
            if n_exposed_in_group > 0 and n_frames > 0:
                # Per-residue average: total_frames / (n_frames * n_exposed)
                mean_frac = contact_frames / (n_frames * n_exposed_in_group)
            else:
                mean_frac = 0.0

            # Calculate contact share (what fraction of this polymer's contacts
            # went to this protein group?)
            if total_poly_contacts > 0:
                contact_share = contact_frames / total_poly_contacts
            else:
                contact_share = 0.0

            # Calculate enrichment (centered at zero)
            # enrichment = (contact_share / expected_share) - 1
            # Positive = prefers this group, Negative = avoids this group
            enrichment = _compute_enrichment(contact_share, expected_share)

            entries.append(
                BindingPreferenceEntry(
                    polymer_type=poly_type,
                    protein_group=group_name,
                    total_contact_frames=contact_frames,
                    mean_contact_fraction=mean_frac,
                    n_residues_in_group=n_total_in_group,
                    n_exposed_in_group=n_exposed_in_group,
                    n_residues_contacted=n_residues_contacted,
                    contact_share=contact_share,
                    expected_share=expected_share,
                    enrichment=enrichment,
                    # Polymer composition metadata (for secondary analysis)
                    polymer_residue_count=poly_res_count,
                    total_polymer_residues=total_poly_residues,
                    polymer_heavy_atom_count=poly_atom_count,
                    total_polymer_heavy_atoms=total_poly_atoms,
                )
            )

    # Compute system-level coverage (collapsed across polymer types)
    system_coverage = _compute_system_coverage(
        entries=entries,
        protein_groups=protein_groups,
        exposed_groups=exposed_groups,
        total_exposed=total_exposed,
        n_frames=n_frames,
        surface_exposure_threshold=surface_exposure.threshold,
        protein_group_selections=protein_group_selections,
        protein_partitions=protein_partitions,
    )

    # Compute per-polymer partition-based binding preference (NEW in v5)
    # This is the primary output - contact_share sums to 1.0 within each partition
    binding_preference = _compute_polymer_binding_preference(
        contact_data=contact_data,
        total_contacts_by_polymer=total_contacts_by_polymer,
        protein_groups=protein_groups,
        exposed_groups=exposed_groups,
        protein_partitions=protein_partitions,
        total_exposed=total_exposed,
        n_frames=n_frames,
        surface_exposure_threshold=surface_exposure.threshold,
        polymer_composition=polymer_composition,
        protein_groups_used=protein_group_selections,
    )

    result = BindingPreferenceResult(
        entries=entries,  # DEPRECATED: kept for backward compat
        n_frames=n_frames,
        total_exposed_residues=total_exposed,
        surface_exposure_threshold=surface_exposure.threshold,
        protein_groups_used=protein_group_selections or {},
        polymer_types_used=polymer_type_selections or {},
        polymer_composition=polymer_composition,
        system_coverage=system_coverage,
        binding_preference=binding_preference,  # NEW: partition-based per-polymer
    )

    # Log summary
    polymer_types_found = result.polymer_types()
    logger.info(
        f"Binding preference computed: {len(entries)} entries for "
        f"{len(polymer_types_found)} polymer types × {len(protein_groups)} groups"
    )
    n_aa_classes = len(system_coverage.aa_class_coverage.entries)
    n_custom_groups = len(system_coverage.custom_group_coverages)
    logger.info(
        f"System coverage computed: {n_aa_classes} AA classes, "
        f"{n_custom_groups} custom groups, "
        f"{system_coverage.total_contact_frames} total contact frames"
    )

    # Log partition-based binding preference summary
    if binding_preference:
        n_user_partitions = len(binding_preference.user_defined_partitions)
        logger.info(
            f"Partition-based binding preference computed: "
            f"{len(binding_preference.polymer_types)} polymer types, "
            f"AA class partition + {n_user_partitions} user partitions"
        )
        # Validate contact_share sums to ~1.0
        for poly_type, aa_result in binding_preference.aa_class_binding.items():
            total_share = aa_result.total_contact_share
            if abs(total_share - 1.0) > 0.01:
                logger.warning(
                    f"AA class partition for {poly_type}: contact_share sums to "
                    f"{total_share:.4f} (expected ~1.0)"
                )

    return result
