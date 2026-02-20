"""Binding preference metrics for polymer-protein contacts.

This module computes enrichment ratios that answer the scientific question:
"Does polymer type X preferentially bind amino acid class Y?"

Enrichment Calculation (Zero-Centered, Surface-Normalized)
-----------------------------------------------------------
For each (polymer_type, protein_group) pair:

    contact_share = Σ(contact_frames for exposed residues in group) /
                    Σ(contact_frames for all exposed residues)

    expected_share = n_exposed_residues_in_group / n_total_exposed_residues

    enrichment = (contact_share / expected_share) - 1

The expected share is based on **protein surface availability**, answering:
"Given that X% of the exposed surface is aromatic, does this polymer type
contact aromatic residues X% of the time, or more/less?"

Interpretation (zero-centered scale):
- enrichment > 0: Preferential binding (more contacts than expected)
    - +0.5 means "50% more contacts than expected"
    - +1.0 means "2× as many contacts as expected"
- enrichment = 0: Neutral (contact frequency matches surface availability)
- enrichment < 0: Avoidance (fewer contacts than expected)
    - -0.3 means "30% fewer contacts than expected"
- enrichment = -1: Complete avoidance (no contacts at all)

The contact-frame weighting ensures that residues contacted for longer durations
contribute proportionally more to the enrichment calculation. A residue contacted
for 60% of the simulation contributes 60× more than one contacted for 1 frame.

Polymer Composition Metadata
----------------------------
The polymer composition (residue counts and heavy atom counts per polymer type)
is stored as metadata for secondary analysis but is NOT used in the enrichment
calculation. The enrichment formula is purely based on protein surface availability.

Examples
--------
>>> from polyzymd.analysis.contacts.binding_preference import (
...     compute_binding_preference,
...     extract_polymer_composition,
... )
>>> from polyzymd.analysis.contacts.surface_exposure import SurfaceExposureFilter
>>>
>>> # Load contact results and compute surface exposure
>>> contact_result = ContactResult.load("contacts.json")
>>> exposure_filter = SurfaceExposureFilter(threshold=0.2)
>>> surface_exposure = exposure_filter.calculate("enzyme.pdb")
>>>
>>> # Define protein groups (resolved to residue IDs)
>>> protein_groups = {
...     "aromatic": {12, 45, 67, 89},  # resids of aromatic residues
...     "charged": {23, 34, 56, 78},
... }
>>>
>>> # Extract polymer composition from trajectory (stored as metadata)
>>> polymer_composition = extract_polymer_composition(universe)
>>>
>>> # Compute binding preference
>>> result = compute_binding_preference(
...     contact_result,
...     surface_exposure,
...     protein_groups,
...     polymer_composition,
... )
>>>
>>> # Check enrichment (surface-normalized)
>>> print(result.get_enrichment("SBM", "aromatic"))
0.45  # SBMA has 45% more contacts with aromatics than surface availability predicts
>>>
>>> print(result.enrichment_matrix())
{'SBM': {'aromatic': 0.45, 'charged': -0.18}, 'EGM': {...}}
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
from pydantic import BaseModel, Field

from polyzymd.analysis.common.aa_classification import DEFAULT_AA_CLASS_SELECTIONS

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

    from polyzymd.analysis.config import ContactsConfig
    from polyzymd.analysis.contacts.results import ContactResult
    from polyzymd.analysis.contacts.surface_exposure import SurfaceExposureResult

logger = logging.getLogger(__name__)


class PolymerComposition(BaseModel):
    """Polymer composition data extracted from trajectory.

    Contains counts of residues and heavy atoms for each polymer type,
    enabling normalization of enrichment ratios by polymer availability.

    This data is used for dual normalization of binding preference:

    - **Residue-based**: Normalizes by number of polymer residues per type.
      Matches the experimental viewpoint where concentrations are specified
      in terms of monomer units.

    - **Atom-based**: Normalizes by number of heavy atoms (non-hydrogen) per type.
      Accounts for differences in monomer size, since larger monomers have
      more surface area and thus more opportunity for contacts.

    Attributes
    ----------
    residue_counts : dict[str, int]
        Number of residues per polymer type (e.g., {"SBM": 50, "EGM": 50})
    heavy_atom_counts : dict[str, int]
        Number of heavy atoms (non-hydrogen) per polymer type.
        Heavy atoms are defined as all atoms with element != 'H'.

    Examples
    --------
    >>> composition = PolymerComposition(
    ...     residue_counts={"SBM": 50, "EGM": 50},
    ...     heavy_atom_counts={"SBM": 750, "EGM": 400},
    ... )
    >>> composition.total_residues
    100
    >>> composition.residue_fraction("SBM")
    0.5
    >>> composition.heavy_atom_fraction("SBM")
    0.652  # SBM has larger monomers
    """

    residue_counts: dict[str, int] = Field(
        default_factory=dict,
        description="Number of residues per polymer type",
    )
    heavy_atom_counts: dict[str, int] = Field(
        default_factory=dict,
        description="Number of heavy atoms (non-H) per polymer type",
    )

    @property
    def total_residues(self) -> int:
        """Total polymer residues across all types."""
        return sum(self.residue_counts.values())

    @property
    def total_heavy_atoms(self) -> int:
        """Total heavy atoms across all polymer types."""
        return sum(self.heavy_atom_counts.values())

    def residue_fraction(self, polymer_type: str) -> float:
        """Fraction of residues that are this polymer type.

        Parameters
        ----------
        polymer_type : str
            Polymer type name (e.g., "SBM")

        Returns
        -------
        float
            Fraction in range [0, 1], or 0.0 if type not found
        """
        total = self.total_residues
        if total == 0:
            return 0.0
        return self.residue_counts.get(polymer_type, 0) / total

    def heavy_atom_fraction(self, polymer_type: str) -> float:
        """Fraction of heavy atoms that are this polymer type.

        Parameters
        ----------
        polymer_type : str
            Polymer type name (e.g., "SBM")

        Returns
        -------
        float
            Fraction in range [0, 1], or 0.0 if type not found
        """
        total = self.total_heavy_atoms
        if total == 0:
            return 0.0
        return self.heavy_atom_counts.get(polymer_type, 0) / total

    def polymer_types(self) -> list[str]:
        """Get sorted list of polymer types in this composition."""
        return sorted(set(self.residue_counts.keys()) | set(self.heavy_atom_counts.keys()))


class BindingPreferenceEntry(BaseModel):
    """Single entry in the binding preference matrix.

    Represents the binding preference metrics for one
    (polymer_type, protein_group) combination.

    Enrichment Interpretation (centered at zero)
    --------------------------------------------
    The enrichment ratio measures whether a polymer type contacts a protein
    group more or less than expected based on **surface availability**.

    - enrichment > 0: Preferential binding (more contacts than expected)
        - +0.5 means "50% more contacts than expected"
        - +1.0 means "2× as many contacts as expected"
    - enrichment = 0: Neutral (contact frequency matches surface availability)
    - enrichment < 0: Avoidance (fewer contacts than expected)
        - -0.3 means "30% fewer contacts than expected"
    - enrichment = -1: Complete avoidance (no contacts at all)

    The expected share is based on protein surface availability:
        expected_share = n_exposed_in_group / total_exposed_residues

    This normalization answers: "Given how much of the protein surface is
    aromatic/charged/etc., does this polymer type contact that surface
    proportionally, more than proportionally, or less?"

    Attributes
    ----------
    polymer_type : str
        Polymer residue type (e.g., "SBM", "EGM")
    protein_group : str
        Protein group label (e.g., "aromatic", "charged_positive")
    total_contact_frames : int
        Sum of contact frames across all exposed residues in this group.
    mean_contact_fraction : float
        Average per-residue contact fraction within this group.
    n_residues_in_group : int
        Total residues in this protein group (exposed + buried)
    n_exposed_in_group : int
        Surface-exposed residues in this group (used for enrichment)
    n_residues_contacted : int
        Number of exposed residues that had at least one contact
    contact_share : float
        Fraction of this polymer's total contacts that went to this group.
    expected_share : float
        Expected contact share based on surface availability
        (n_exposed_in_group / total_exposed_residues)
    enrichment : float | None
        Zero-centered enrichment: (contact_share / expected_share) - 1
    polymer_residue_count : int
        Number of residues of this polymer type (metadata)
    total_polymer_residues : int
        Total polymer residues across all types (metadata)
    polymer_heavy_atom_count : int
        Number of heavy atoms for this polymer type (metadata)
    total_polymer_heavy_atoms : int
        Total polymer heavy atoms across all types (metadata)
    """

    polymer_type: str
    protein_group: str
    total_contact_frames: int = Field(
        description="Sum of contact frames for all exposed residues in group"
    )
    mean_contact_fraction: float = Field(
        description="Average per-residue contact fraction in group"
    )
    n_residues_in_group: int = Field(description="Total residues in this protein group")
    n_exposed_in_group: int = Field(description="Surface-exposed residues in group")
    n_residues_contacted: int = Field(
        default=0, description="Exposed residues with at least one contact"
    )
    contact_share: float = Field(
        default=0.0, description="Fraction of polymer's contacts to this group"
    )
    expected_share: float = Field(
        default=0.0,
        description="Expected share based on protein surface availability",
    )
    enrichment: float | None = Field(
        default=None,
        description="Zero-centered enrichment: (contact_share / expected_share) - 1",
    )

    # Polymer composition metadata (for secondary analysis)
    polymer_residue_count: int = Field(
        default=0,
        description="Number of residues of this polymer type (metadata)",
    )
    total_polymer_residues: int = Field(
        default=0,
        description="Total polymer residues across all types (metadata)",
    )
    polymer_heavy_atom_count: int = Field(
        default=0,
        description="Number of heavy atoms (non-H) for this polymer type (metadata)",
    )
    total_polymer_heavy_atoms: int = Field(
        default=0,
        description="Total heavy atoms across all polymer types (metadata)",
    )


class BindingPreferenceResult(BaseModel):
    """Complete binding preference analysis result.

    Provides enrichment-normalized metrics for polymer-protein
    binding preferences, answering questions like:

    - "Does SBMA preferentially bind aromatic residues?"
    - "How does EGMA's preference for charged residues compare to SBMA?"
    - "Which amino acid class does this polymer type prefer?"

    Enrichment values are centered at zero:
    - enrichment > 0: Preferential binding
    - enrichment = 0: Neutral (random chance)
    - enrichment < 0: Avoidance

    Attributes
    ----------
    entries : list[BindingPreferenceEntry]
        All (polymer_type × protein_group) combinations
    n_frames : int
        Total frames analyzed
    total_exposed_residues : int
        Number of surface-exposed residues considered
    surface_exposure_threshold : float
        SASA threshold used for surface filtering
    protein_groups_used : dict[str, str]
        Mapping of group name to MDAnalysis selection string
    polymer_types_used : dict[str, str]
        Mapping of polymer type name to MDAnalysis selection string
    polymer_composition : PolymerComposition
        Polymer composition data (residue/atom counts per type, metadata only)
    system_coverage : SystemCoverageResult, optional
        System-level coverage metrics collapsed across polymer types.
        Answers: "What does this polymer mixture collectively cover?"
    schema_version : int
        Version for forward compatibility. Version 4 adds system_coverage.
    """

    entries: list[BindingPreferenceEntry] = Field(default_factory=list)
    n_frames: int = 0
    total_exposed_residues: int = 0
    surface_exposure_threshold: float | None = None
    protein_groups_used: dict[str, str] = Field(default_factory=dict)
    polymer_types_used: dict[str, str] = Field(default_factory=dict)
    polymer_composition: PolymerComposition | None = Field(
        default=None,
        description="Polymer composition data (residue/atom counts per type, metadata only)",
    )
    system_coverage: "SystemCoverageResult | None" = Field(
        default=None,
        description="System-level coverage metrics collapsed across polymer types",
    )
    metadata: dict[str, Any] = Field(default_factory=dict)
    schema_version: int = 4  # Version 4: adds system_coverage field

    def to_dataframe(self) -> "pd.DataFrame":
        """Convert to pandas DataFrame for analysis/plotting.

        Returns
        -------
        pd.DataFrame
            Columns: polymer_type, protein_group, total_contact_frames,
            mean_contact_fraction, n_residues_in_group, n_exposed_in_group,
            n_residues_contacted, contact_share, expected_share, enrichment
        """
        import pandas as pd

        return pd.DataFrame([e.model_dump() for e in self.entries])

    def enrichment_matrix(self) -> dict[str, dict[str, float]]:
        """Get enrichment as nested dict: {polymer_type: {protein_group: value}}.

        Enrichment values are centered at zero and normalized by protein
        surface availability:
        - > 0: Preferential binding (more contacts than expected)
        - = 0: Neutral (matches surface availability)
        - < 0: Avoidance (fewer contacts than expected)

        Returns
        -------
        dict[str, dict[str, float]]
            Nested mapping of enrichment values.
            Missing/invalid values are returned as 0.0.

        Examples
        --------
        >>> matrix = result.enrichment_matrix()
        >>> print(matrix["SBM"]["aromatic"])
        0.45  # 45% more contacts than expected based on surface availability
        """
        result: dict[str, dict[str, float]] = {}
        for entry in self.entries:
            if entry.polymer_type not in result:
                result[entry.polymer_type] = {}
            value = entry.enrichment
            result[entry.polymer_type][entry.protein_group] = value if value is not None else 0.0
        return result

    def contact_fraction_matrix(self) -> dict[str, dict[str, float]]:
        """Get mean contact fractions as nested dict.

        Returns
        -------
        dict[str, dict[str, float]]
            Nested mapping: {polymer_type: {protein_group: mean_frac}}
        """
        result: dict[str, dict[str, float]] = {}
        for entry in self.entries:
            if entry.polymer_type not in result:
                result[entry.polymer_type] = {}
            result[entry.polymer_type][entry.protein_group] = entry.mean_contact_fraction
        return result

    def contact_share_matrix(self) -> dict[str, dict[str, float]]:
        """Get contact shares as nested dict.

        Returns
        -------
        dict[str, dict[str, float]]
            Nested mapping: {polymer_type: {protein_group: contact_share}}
        """
        result: dict[str, dict[str, float]] = {}
        for entry in self.entries:
            if entry.polymer_type not in result:
                result[entry.polymer_type] = {}
            result[entry.polymer_type][entry.protein_group] = entry.contact_share
        return result

    def get_enrichment(self, polymer_type: str, protein_group: str) -> float | None:
        """Get enrichment for a specific (polymer_type, protein_group) pair.

        Parameters
        ----------
        polymer_type : str
            Polymer type name
        protein_group : str
            Protein group name

        Returns
        -------
        float or None
            Enrichment value (centered at zero), or None if pair not found.
            Enrichment is based on protein surface availability:
            (contact_share / expected_share) - 1
        """
        for entry in self.entries:
            if entry.polymer_type == polymer_type and entry.protein_group == protein_group:
                return entry.enrichment
        return None

    def get_entry(self, polymer_type: str, protein_group: str) -> BindingPreferenceEntry | None:
        """Get the full entry for a (polymer_type, protein_group) pair.

        Parameters
        ----------
        polymer_type : str
            Polymer type name
        protein_group : str
            Protein group name

        Returns
        -------
        BindingPreferenceEntry or None
            Full entry, or None if not found
        """
        for entry in self.entries:
            if entry.polymer_type == polymer_type and entry.protein_group == protein_group:
                return entry
        return None

    def polymer_types(self) -> list[str]:
        """Get list of polymer types in this result."""
        return sorted(set(e.polymer_type for e in self.entries))

    def protein_groups(self) -> list[str]:
        """Get list of protein groups in this result."""
        return sorted(set(e.protein_group for e in self.entries))

    def save(self, path: str | Path) -> None:
        """Save to JSON file.

        Parameters
        ----------
        path : str or Path
            Output path for JSON file
        """
        Path(path).write_text(json.dumps(self.model_dump(), indent=2))
        logger.info(f"Saved binding preference result to {path}")

    @classmethod
    def load(cls, path: str | Path) -> "BindingPreferenceResult":
        """Load from JSON file.

        Parameters
        ----------
        path : str or Path
            Path to JSON file

        Returns
        -------
        BindingPreferenceResult
            Loaded result
        """
        data = json.loads(Path(path).read_text())
        return cls.model_validate(data)


# =============================================================================
# System Coverage Helper Functions
# =============================================================================


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


def _compute_partition_coverage(
    partition_name: str,
    partition_type: Literal["aa_class", "binary_custom", "combined_custom"],
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
        One of: "aa_class", "binary_custom", "combined_custom"
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

    Returns
    -------
    SystemCoverageResult
        Partition-based coverage metrics (schema v2)
    """
    from polyzymd.analysis.common.aa_classification import DEFAULT_AA_CLASS_SELECTIONS

    # Collect all polymer types
    all_polymer_types = sorted(set(e.polymer_type for e in entries))

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
    # 5. Compute total contact frames
    # ---------------------------------------------------------------------
    total_contact_frames = sum(e.total_contact_frames for e in entries)

    return SystemCoverageResult(
        aa_class_coverage=aa_class_coverage,
        custom_group_coverages=custom_group_coverages,
        combined_custom_coverage=combined_custom_coverage,
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
    contact_result: "ContactResult",
    surface_exposure: "SurfaceExposureResult",
    protein_groups: dict[str, set[int]],
    polymer_composition: PolymerComposition,
    polymer_types: list[str] | None = None,
    protein_group_selections: dict[str, str] | None = None,
    polymer_type_selections: dict[str, str] | None = None,
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
    contact_result : ContactResult
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

    Returns
    -------
    BindingPreferenceResult
        Binding preference metrics with dual enrichment normalization

    Notes
    -----
    Enrichment Calculation (centered at zero)
    -----------------------------------------
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

    # Helper function for enrichment calculation (centered at zero)
    def calc_enrichment(contact_share: float, expected: float) -> float | None:
        """Calculate enrichment: (observed / expected) - 1, centered at zero."""
        if expected > 0 and contact_share > 0:
            return (contact_share / expected) - 1
        elif expected > 0 and contact_share == 0:
            return -1.0  # Complete avoidance
        else:
            return None  # Cannot compute (no expected share)

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
            enrichment = calc_enrichment(contact_share, expected_share)

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
    )

    result = BindingPreferenceResult(
        entries=entries,
        n_frames=n_frames,
        total_exposed_residues=total_exposed,
        surface_exposure_threshold=surface_exposure.threshold,
        protein_groups_used=protein_group_selections or {},
        polymer_types_used=polymer_type_selections or {},
        polymer_composition=polymer_composition,
        system_coverage=system_coverage,
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

    return result


def aggregate_binding_preference(
    results: list[BindingPreferenceResult],
) -> "AggregatedBindingPreferenceResult":
    """Aggregate binding preference across replicates.

    Computes mean ± SEM for both residue-based and atom-based enrichment
    ratios across multiple replicates.

    Parameters
    ----------
    results : list[BindingPreferenceResult]
        Binding preference results from multiple replicates

    Returns
    -------
    AggregatedBindingPreferenceResult
        Aggregated results with mean and SEM for both normalization methods
    """
    if not results:
        raise ValueError("No results to aggregate")

    # Helper function for computing mean and SEM
    def _compute_stats(values: list[float]) -> tuple[float | None, float | None]:
        """Compute mean and SEM from a list of values."""
        n = len(values)
        if n == 0:
            return None, None
        mean_val = float(np.mean(values))
        sem_val = float(np.std(values, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        return mean_val, sem_val

    # Collect all (polymer_type, protein_group) pairs
    all_pairs: set[tuple[str, str]] = set()
    for r in results:
        for e in r.entries:
            all_pairs.add((e.polymer_type, e.protein_group))

    entries = []
    for poly_type, prot_group in sorted(all_pairs):
        # Collect values from each replicate
        enrichments = []
        contact_fractions = []
        contact_shares = []

        for r in results:
            entry = r.get_entry(poly_type, prot_group)
            if entry is not None:
                if entry.enrichment is not None:
                    enrichments.append(entry.enrichment)
                contact_fractions.append(entry.mean_contact_fraction)
                contact_shares.append(entry.contact_share)

        # Compute statistics
        mean_enrichment, sem_enrichment = _compute_stats(enrichments)
        mean_contact_fraction, sem_contact_fraction = _compute_stats(contact_fractions)
        mean_contact_share = float(np.mean(contact_shares)) if contact_shares else 0.0

        # Get group metadata from first result
        first_entry = results[0].get_entry(poly_type, prot_group)
        n_exposed = first_entry.n_exposed_in_group if first_entry else 0
        n_total = first_entry.n_residues_in_group if first_entry else 0
        expected_share = first_entry.expected_share if first_entry else 0.0

        entries.append(
            AggregatedBindingPreferenceEntry(
                polymer_type=poly_type,
                protein_group=prot_group,
                # Enrichment (surface-normalized)
                mean_enrichment=mean_enrichment,
                sem_enrichment=sem_enrichment,
                per_replicate_enrichments=enrichments,
                # Contact metrics
                mean_contact_fraction=mean_contact_fraction if mean_contact_fraction else 0.0,
                sem_contact_fraction=sem_contact_fraction if sem_contact_fraction else 0.0,
                mean_contact_share=mean_contact_share,
                # Expected share (from protein surface)
                expected_share=expected_share,
                # Group metadata
                n_exposed_in_group=n_exposed,
                n_residues_in_group=n_total,
                n_replicates=len(enrichments),
            )
        )

    # Aggregate system coverage if present in all results
    aggregated_system_coverage = None
    system_coverages = [r.system_coverage for r in results if r.system_coverage is not None]
    if len(system_coverages) == len(results) and len(system_coverages) > 0:
        aggregated_system_coverage = aggregate_system_coverage(system_coverages)
        logger.debug(
            f"Aggregated system coverage: {len(aggregated_system_coverage.entries)} groups "
            f"from {len(system_coverages)} replicates"
        )

    return AggregatedBindingPreferenceResult(
        entries=entries,
        n_replicates=len(results),
        total_exposed_residues=results[0].total_exposed_residues if results else 0,
        surface_exposure_threshold=results[0].surface_exposure_threshold if results else None,
        protein_groups_used=results[0].protein_groups_used if results else {},
        polymer_types_used=results[0].polymer_types_used if results else {},
        polymer_composition=results[0].polymer_composition if results else None,
        system_coverage=aggregated_system_coverage,
    )


class AggregatedBindingPreferenceEntry(BaseModel):
    """Aggregated binding preference for one (polymer_type, protein_group) pair.

    Contains mean ± SEM across replicates for enrichment based on
    protein surface availability.

    Enrichment values are centered at zero:
    - > 0: Preferential binding (more contacts than expected by surface area)
    - = 0: Neutral (contact frequency matches surface availability)
    - < 0: Avoidance (fewer contacts than expected by surface area)

    The expected share is based on protein surface availability:
        expected_share = n_exposed_in_group / total_exposed_residues

    This normalization answers: "Given how much of the protein surface is
    aromatic/charged/etc., does this polymer type contact that surface
    proportionally, more than proportionally, or less?"
    """

    polymer_type: str
    protein_group: str

    # Enrichment (surface-availability normalized)
    mean_enrichment: float | None = Field(
        default=None,
        description="Mean enrichment across replicates (surface-normalized)",
    )
    sem_enrichment: float | None = Field(
        default=None,
        description="Standard error of enrichment",
    )
    per_replicate_enrichments: list[float] = Field(
        default_factory=list,
        description="Enrichment values from each replicate",
    )

    # Contact metrics
    mean_contact_fraction: float = Field(
        default=0.0,
        description="Mean per-residue contact fraction",
    )
    sem_contact_fraction: float = Field(
        default=0.0,
        description="Standard error of contact fraction",
    )
    mean_contact_share: float = Field(
        default=0.0,
        description="Mean contact share",
    )

    # Expected share (from protein surface availability)
    expected_share: float = Field(
        default=0.0,
        description="Expected contact share based on protein surface availability",
    )

    # Group metadata
    n_exposed_in_group: int = Field(
        default=0,
        description="Surface-exposed residues in group",
    )
    n_residues_in_group: int = Field(
        default=0,
        description="Total residues in group",
    )
    n_replicates: int = Field(
        default=0,
        description="Number of replicates with valid data",
    )


class AggregatedBindingPreferenceResult(BaseModel):
    """Binding preference aggregated across replicates.

    Contains mean ± SEM for all metrics across multiple replicates.
    Enrichment is normalized by protein surface availability.

    Enrichment values are centered at zero:
    - > 0: Preferential binding (more contacts than expected by surface area)
    - = 0: Neutral (contact frequency matches surface availability)
    - < 0: Avoidance (fewer contacts than expected by surface area)
    """

    entries: list[AggregatedBindingPreferenceEntry] = Field(default_factory=list)
    n_replicates: int = 0
    total_exposed_residues: int = 0
    surface_exposure_threshold: float | None = None
    protein_groups_used: dict[str, str] = Field(default_factory=dict)
    polymer_types_used: dict[str, str] = Field(default_factory=dict)
    polymer_composition: PolymerComposition | None = Field(
        default=None,
        description="Polymer composition data (residue/atom counts per type, metadata only)",
    )
    system_coverage: "AggregatedSystemCoverageResult | None" = Field(
        default=None,
        description="Aggregated system-level coverage metrics",
    )
    schema_version: int = 4  # Version 4: adds system_coverage field

    def to_dataframe(self) -> "pd.DataFrame":
        """Convert to pandas DataFrame."""
        import pandas as pd

        return pd.DataFrame([e.model_dump() for e in self.entries])

    def enrichment_matrix(self) -> dict[str, dict[str, float]]:
        """Get mean enrichment as nested dict.

        Enrichment values are centered at zero and normalized by protein
        surface availability:
        - > 0: Preferential binding (more contacts than expected)
        - = 0: Neutral (matches surface availability)
        - < 0: Avoidance (fewer contacts than expected)

        Returns
        -------
        dict[str, dict[str, float]]
            Nested mapping: {polymer_type: {protein_group: mean_enrichment}}.
            Missing/invalid values are returned as 0.0.
        """
        result: dict[str, dict[str, float]] = {}
        for entry in self.entries:
            if entry.polymer_type not in result:
                result[entry.polymer_type] = {}

            value = entry.mean_enrichment
            result[entry.polymer_type][entry.protein_group] = value if value is not None else 0.0
        return result

    def get_entry(
        self, polymer_type: str, protein_group: str
    ) -> AggregatedBindingPreferenceEntry | None:
        """Get entry for a (polymer_type, protein_group) pair."""
        for entry in self.entries:
            if entry.polymer_type == polymer_type and entry.protein_group == protein_group:
                return entry
        return None

    def polymer_types(self) -> list[str]:
        """Get list of polymer types."""
        return sorted(set(e.polymer_type for e in self.entries))

    def protein_groups(self) -> list[str]:
        """Get list of protein groups."""
        return sorted(set(e.protein_group for e in self.entries))

    def save(self, path: str | Path) -> None:
        """Save to JSON file."""
        Path(path).write_text(json.dumps(self.model_dump(), indent=2))
        logger.info(f"Saved aggregated binding preference to {path}")

    @classmethod
    def load(cls, path: str | Path) -> "AggregatedBindingPreferenceResult":
        """Load from JSON file."""
        data = json.loads(Path(path).read_text())
        return cls.model_validate(data)


# =============================================================================
# System Coverage Models
# =============================================================================


class SystemCoverageEntry(BaseModel):
    """System-level coverage entry for one protein group.

    While BindingPreferenceEntry answers "What does SBMA prefer?", this entry
    answers "What fraction of ALL polymer contacts in this system go to this
    protein group?" — collapsing across polymer types for condition-level analysis.

    Use Case
    --------
    Compare copolymer compositions (conditions) with each other. For example:
    "Does a 70:30 SBMA:EGMA mixture cover aromatic residues differently than
    a 30:70 mixture?"

    Coverage Enrichment Calculation
    -------------------------------
    For each protein group:

        coverage_share = Σ(all polymer contacts to group) / Σ(all polymer contacts)
        expected_share = n_exposed_in_group / total_exposed_residues
        coverage_enrichment = (coverage_share / expected_share) - 1

    Interpretation (centered at zero):
    - coverage_enrichment > 0: Preferential coverage (more than surface predicts)
        - +0.5 means "50% more coverage than expected"
    - coverage_enrichment = 0: Neutral (coverage matches surface availability)
    - coverage_enrichment < 0: Under-coverage (less than surface predicts)
        - -0.3 means "30% less coverage than expected"
    - coverage_enrichment = -1: Complete avoidance (no coverage at all)

    Attributes
    ----------
    protein_group : str
        Protein group label (e.g., "aromatic", "charged_positive")
    total_contact_frames : int
        Sum of contact frames from ALL polymer types to this group.
    coverage_share : float
        Fraction of all polymer contacts that went to this group.
    expected_share : float
        Expected coverage based on protein surface availability
        (n_exposed_in_group / total_exposed_residues)
    coverage_enrichment : float | None
        Zero-centered enrichment: (coverage_share / expected_share) - 1
    n_exposed_in_group : int
        Surface-exposed residues in this group
    n_residues_in_group : int
        Total residues in this group (exposed + buried)
    polymer_contributions : dict[str, float]
        Breakdown of coverage by polymer type: {"SBMA": 0.35, "EGMA": 0.65}
        Values sum to 1.0 (fraction of contacts to this group from each polymer)
    """

    protein_group: str
    total_contact_frames: int = Field(
        default=0,
        description="Sum of contact frames from ALL polymer types to this group",
    )
    coverage_share: float = Field(
        default=0.0,
        description="Fraction of all polymer contacts that went to this group",
    )
    expected_share: float = Field(
        default=0.0,
        description="Expected coverage based on protein surface availability",
    )
    coverage_enrichment: float | None = Field(
        default=None,
        description="Zero-centered enrichment: (coverage_share / expected_share) - 1",
    )
    n_exposed_in_group: int = Field(
        default=0,
        description="Surface-exposed residues in this group",
    )
    n_residues_in_group: int = Field(
        default=0,
        description="Total residues in this group",
    )
    polymer_contributions: dict[str, float] = Field(
        default_factory=dict,
        description="Fraction of contacts to this group from each polymer type (sums to 1.0)",
    )


# =============================================================================
# Partition-Based Coverage Models (Schema v2)
# =============================================================================
# These models fix the overlap bug in the original SystemCoverageResult by
# enforcing proper partition semantics: elements are mutually exclusive and
# collectively exhaustive, so coverage_share and expected_share both sum to 1.0.


class PartitionCoverageEntry(BaseModel):
    """Coverage metrics for one element in a partition.

    A partition element is a mutually exclusive subset of protein residues.
    Within a partition, all elements together cover the entire protein surface
    exactly once (no residue is counted in multiple elements).

    This ensures that:
    - coverage_share sums to 1.0 across all elements in the partition
    - expected_share sums to 1.0 across all elements in the partition
    - enrichment is mathematically valid (no inflated denominators)

    Attributes
    ----------
    partition_element : str
        Name of this partition element (e.g., "aromatic", "lid_helix_5", "rest_of_protein")
    total_contact_frames : int
        Sum of contact frames from ALL polymer types to residues in this element
    coverage_share : float
        Fraction of all polymer contacts that went to this element.
        Sums to 1.0 across all elements in the partition.
    expected_share : float
        Expected coverage based on surface availability (n_exposed / total_exposed).
        Sums to 1.0 across all elements in the partition.
    coverage_enrichment : float | None
        Zero-centered enrichment: (coverage_share / expected_share) - 1
    n_exposed_in_element : int
        Number of surface-exposed residues in this element
    n_residues_in_element : int
        Total residues in this element (exposed + buried)
    polymer_contributions : dict[str, float]
        Breakdown of coverage by polymer type (sums to 1.0 for this element)
    """

    partition_element: str
    total_contact_frames: int = Field(
        default=0,
        description="Sum of contact frames from ALL polymer types to this element",
    )
    coverage_share: float = Field(
        default=0.0,
        description="Fraction of all polymer contacts that went to this element",
    )
    expected_share: float = Field(
        default=0.0,
        description="Expected coverage based on protein surface availability",
    )
    coverage_enrichment: float | None = Field(
        default=None,
        description="Zero-centered enrichment: (coverage_share / expected_share) - 1",
    )
    n_exposed_in_element: int = Field(
        default=0,
        description="Surface-exposed residues in this element",
    )
    n_residues_in_element: int = Field(
        default=0,
        description="Total residues in this element",
    )
    polymer_contributions: dict[str, float] = Field(
        default_factory=dict,
        description="Fraction of contacts to this element from each polymer type",
    )


class PartitionCoverageResult(BaseModel):
    """Coverage analysis for a complete partition of the protein surface.

    A partition divides the protein surface into mutually exclusive regions
    that together cover the entire surface. This ensures valid enrichment
    calculations where both coverage_share and expected_share sum to 1.0.

    Partition Types
    ---------------
    - **aa_class**: 5-way partition by amino acid class
      (aromatic, polar, nonpolar, charged_positive, charged_negative)
    - **binary_custom**: 2-way partition for a custom group vs rest_of_protein
      (e.g., "lid_helix_5" vs "rest_of_protein")
    - **combined_custom**: N+1 way partition with all non-overlapping custom groups
      plus "rest_of_protein"

    Attributes
    ----------
    partition_name : str
        Descriptive name (e.g., "aa_class", "lid_helix_5_vs_rest")
    partition_type : str
        One of: "aa_class", "binary_custom", "combined_custom"
    entries : list[PartitionCoverageEntry]
        Coverage metrics for each element in the partition
    total_coverage_share : float
        Validation check: should be ~1.0
    total_expected_share : float
        Validation check: should be ~1.0
    """

    partition_name: str
    partition_type: Literal["aa_class", "binary_custom", "combined_custom"]
    entries: list[PartitionCoverageEntry] = Field(default_factory=list)
    total_coverage_share: float = Field(
        default=1.0,
        description="Sum of coverage_share across elements (validation: should be ~1.0)",
    )
    total_expected_share: float = Field(
        default=1.0,
        description="Sum of expected_share across elements (validation: should be ~1.0)",
    )

    def to_dataframe(self) -> "pd.DataFrame":
        """Convert to pandas DataFrame for analysis/plotting."""
        import pandas as pd

        return pd.DataFrame([e.model_dump() for e in self.entries])

    def coverage_enrichment_dict(self) -> dict[str, float]:
        """Get coverage enrichment as dict: {element: enrichment}."""
        return {
            e.partition_element: (
                e.coverage_enrichment if e.coverage_enrichment is not None else 0.0
            )
            for e in self.entries
        }

    def coverage_share_dict(self) -> dict[str, float]:
        """Get coverage shares as dict: {element: share}."""
        return {e.partition_element: e.coverage_share for e in self.entries}

    def expected_share_dict(self) -> dict[str, float]:
        """Get expected shares as dict: {element: share}."""
        return {e.partition_element: e.expected_share for e in self.entries}

    def get_entry(self, element_name: str) -> PartitionCoverageEntry | None:
        """Get the entry for a specific partition element."""
        for entry in self.entries:
            if entry.partition_element == element_name:
                return entry
        return None

    def element_names(self) -> list[str]:
        """Get list of partition element names."""
        return [e.partition_element for e in self.entries]


class SystemCoverageResult(BaseModel):
    """System-level coverage analysis with proper partition structure.

    This result uses partitions to ensure mathematically valid enrichment
    calculations. A partition divides the protein surface into mutually
    exclusive regions, avoiding the overlap bug where custom groups and
    AA class groups can inflate the expected_share denominator.

    Partition Strategy
    ------------------
    1. **AA Class Partition** (always computed):
       5-way partition by amino acid class. Every surface residue belongs
       to exactly one class.

    2. **Binary Custom Partitions** (per custom group):
       Each custom group is compared to "rest_of_protein". This answers:
       "Does my lid_helix_5 have enriched polymer contacts vs non-lid regions?"

    3. **Combined Custom Partition** (optional):
       If custom groups don't overlap, all custom groups + rest_of_protein
       form a single partition. If groups overlap, this is not computed
       and an error is raised if explicitly requested.

    Attributes
    ----------
    aa_class_coverage : PartitionCoverageResult
        5-way partition by amino acid class. Always computed.
    custom_group_coverages : dict[str, PartitionCoverageResult]
        Binary partitions for each custom group vs rest_of_protein.
        Keys are custom group names.
    combined_custom_coverage : PartitionCoverageResult | None
        All custom groups + rest_of_protein as a single partition.
        Only computed if custom groups don't overlap.
    n_frames : int
        Total frames analyzed
    total_contact_frames : int
        Sum of all polymer contacts across all groups
    total_exposed_residues : int
        Number of surface-exposed protein residues
    surface_exposure_threshold : float | None
        SASA threshold used for surface filtering
    custom_group_selections : dict[str, str]
        Custom group name to MDAnalysis selection (for metadata)
    polymer_types_included : list[str]
        Polymer types that contributed to coverage
    has_overlapping_custom_groups : bool
        True if custom groups share residues (combined partition not computed)
    overlapping_group_pairs : list[tuple[str, str]]
        Pairs of custom groups that overlap (for diagnostics)
    schema_version : int
        Schema version (2 = partition-based)
    """

    aa_class_coverage: PartitionCoverageResult
    custom_group_coverages: dict[str, PartitionCoverageResult] = Field(default_factory=dict)
    combined_custom_coverage: PartitionCoverageResult | None = None

    # Metadata
    n_frames: int = 0
    total_contact_frames: int = Field(
        default=0,
        description="Sum of all polymer contacts across all groups",
    )
    total_exposed_residues: int = 0
    surface_exposure_threshold: float | None = None
    custom_group_selections: dict[str, str] = Field(default_factory=dict)
    polymer_types_included: list[str] = Field(default_factory=list)
    has_overlapping_custom_groups: bool = False
    overlapping_group_pairs: list[tuple[str, str]] = Field(default_factory=list)
    schema_version: int = 2

    def get_aa_class_enrichment(self, aa_class: str) -> float | None:
        """Get coverage enrichment for an AA class.

        Parameters
        ----------
        aa_class : str
            One of: aromatic, polar, nonpolar, charged_positive, charged_negative

        Returns
        -------
        float | None
            Coverage enrichment, or None if not found
        """
        entry = self.aa_class_coverage.get_entry(aa_class)
        return entry.coverage_enrichment if entry else None

    def get_custom_group_enrichment(self, group_name: str) -> float | None:
        """Get coverage enrichment for a custom group (vs rest_of_protein).

        Parameters
        ----------
        group_name : str
            Custom group name (e.g., "lid_helix_5")

        Returns
        -------
        float | None
            Coverage enrichment for the custom group, or None if not found
        """
        if group_name not in self.custom_group_coverages:
            return None
        partition = self.custom_group_coverages[group_name]
        entry = partition.get_entry(group_name)
        return entry.coverage_enrichment if entry else None

    def aa_class_enrichment_dict(self) -> dict[str, float]:
        """Get AA class enrichments as dict: {aa_class: enrichment}."""
        return self.aa_class_coverage.coverage_enrichment_dict()

    def custom_group_enrichment_dict(self) -> dict[str, float]:
        """Get custom group enrichments as dict: {group_name: enrichment}.

        Each custom group's enrichment is relative to rest_of_protein.
        """
        result = {}
        for group_name, partition in self.custom_group_coverages.items():
            entry = partition.get_entry(group_name)
            if entry and entry.coverage_enrichment is not None:
                result[group_name] = entry.coverage_enrichment
            else:
                result[group_name] = 0.0
        return result

    def aa_class_names(self) -> list[str]:
        """Get list of AA class names in canonical order."""
        canonical_order = ["aromatic", "polar", "nonpolar", "charged_positive", "charged_negative"]
        names = self.aa_class_coverage.element_names()
        return [n for n in canonical_order if n in names]

    def custom_group_names(self) -> list[str]:
        """Get list of custom group names."""
        return sorted(self.custom_group_coverages.keys())

    def save(self, path: str | Path) -> None:
        """Save to JSON file.

        Parameters
        ----------
        path : str or Path
            Output path for JSON file
        """
        Path(path).write_text(json.dumps(self.model_dump(), indent=2))
        logger.info(f"Saved system coverage result to {path}")

    @classmethod
    def load(cls, path: str | Path) -> "SystemCoverageResult":
        """Load from JSON file.

        Parameters
        ----------
        path : str or Path
            Path to JSON file

        Returns
        -------
        SystemCoverageResult
            Loaded result
        """
        data = json.loads(Path(path).read_text())
        return cls.model_validate(data)


# =============================================================================
# Aggregated Partition Coverage Models (Schema v2)
# =============================================================================


class AggregatedPartitionCoverageEntry(BaseModel):
    """Aggregated coverage for one partition element across replicates.

    Contains mean ± SEM for coverage metrics, enabling statistical comparison
    of coverage enrichment across conditions.

    Attributes
    ----------
    partition_element : str
        Element name (e.g., "aromatic", "lid_helix_5", "rest_of_protein")
    mean_coverage_share : float
        Mean coverage share across replicates
    sem_coverage_share : float
        Standard error of coverage share
    mean_coverage_enrichment : float | None
        Mean coverage enrichment across replicates
    sem_coverage_enrichment : float | None
        Standard error of coverage enrichment
    per_replicate_enrichments : list[float]
        Coverage enrichment values from each replicate
    expected_share : float
        Expected coverage based on surface availability
    n_exposed_in_element : int
        Surface-exposed residues in this element
    n_residues_in_element : int
        Total residues in this element
    n_replicates : int
        Number of replicates with valid data
    mean_polymer_contributions : dict[str, float]
        Mean polymer contributions across replicates
    """

    partition_element: str
    mean_coverage_share: float = Field(
        default=0.0,
        description="Mean coverage share across replicates",
    )
    sem_coverage_share: float = Field(
        default=0.0,
        description="Standard error of coverage share",
    )
    mean_coverage_enrichment: float | None = Field(
        default=None,
        description="Mean coverage enrichment across replicates",
    )
    sem_coverage_enrichment: float | None = Field(
        default=None,
        description="Standard error of coverage enrichment",
    )
    per_replicate_enrichments: list[float] = Field(
        default_factory=list,
        description="Coverage enrichment values from each replicate",
    )
    expected_share: float = Field(
        default=0.0,
        description="Expected coverage based on surface availability",
    )
    n_exposed_in_element: int = Field(
        default=0,
        description="Surface-exposed residues in this element",
    )
    n_residues_in_element: int = Field(
        default=0,
        description="Total residues in this element",
    )
    n_replicates: int = Field(
        default=0,
        description="Number of replicates with valid data",
    )
    mean_polymer_contributions: dict[str, float] = Field(
        default_factory=dict,
        description="Mean polymer contributions across replicates",
    )


class AggregatedPartitionCoverageResult(BaseModel):
    """Aggregated coverage for a partition across replicates.

    Contains mean ± SEM for all elements in the partition.

    Attributes
    ----------
    partition_name : str
        Name of the partition
    partition_type : str
        One of: "aa_class", "binary_custom", "combined_custom"
    entries : list[AggregatedPartitionCoverageEntry]
        Aggregated coverage for each element
    n_replicates : int
        Number of replicates aggregated
    """

    partition_name: str
    partition_type: Literal["aa_class", "binary_custom", "combined_custom"]
    entries: list[AggregatedPartitionCoverageEntry] = Field(default_factory=list)
    n_replicates: int = 0

    def to_dataframe(self) -> "pd.DataFrame":
        """Convert to pandas DataFrame."""
        import pandas as pd

        return pd.DataFrame([e.model_dump() for e in self.entries])

    def coverage_enrichment_dict(self) -> dict[str, float]:
        """Get mean coverage enrichment as dict: {element: enrichment}."""
        return {
            e.partition_element: (
                e.mean_coverage_enrichment if e.mean_coverage_enrichment is not None else 0.0
            )
            for e in self.entries
        }

    def get_entry(self, element_name: str) -> AggregatedPartitionCoverageEntry | None:
        """Get entry for a specific partition element."""
        for entry in self.entries:
            if entry.partition_element == element_name:
                return entry
        return None

    def element_names(self) -> list[str]:
        """Get list of partition element names."""
        return [e.partition_element for e in self.entries]


class AggregatedSystemCoverageResult(BaseModel):
    """System coverage aggregated across replicates (schema v2).

    Contains aggregated partition coverages with mean ± SEM statistics
    for statistical comparison between conditions.

    Attributes
    ----------
    aa_class_coverage : AggregatedPartitionCoverageResult
        Aggregated 5-way AA class partition
    custom_group_coverages : dict[str, AggregatedPartitionCoverageResult]
        Aggregated binary partitions for each custom group
    combined_custom_coverage : AggregatedPartitionCoverageResult | None
        Aggregated combined custom partition (if applicable)
    n_replicates : int
        Number of replicates aggregated
    total_exposed_residues : int
        Number of surface-exposed protein residues
    surface_exposure_threshold : float | None
        SASA threshold used for surface filtering
    custom_group_selections : dict[str, str]
        Custom group name to MDAnalysis selection
    polymer_types_included : list[str]
        Polymer types that contributed to coverage
    has_overlapping_custom_groups : bool
        True if custom groups share residues
    schema_version : int
        Schema version (2 = partition-based)
    """

    aa_class_coverage: AggregatedPartitionCoverageResult
    custom_group_coverages: dict[str, AggregatedPartitionCoverageResult] = Field(
        default_factory=dict
    )
    combined_custom_coverage: AggregatedPartitionCoverageResult | None = None

    # Metadata
    n_replicates: int = 0
    total_exposed_residues: int = 0
    surface_exposure_threshold: float | None = None
    custom_group_selections: dict[str, str] = Field(default_factory=dict)
    polymer_types_included: list[str] = Field(default_factory=list)
    has_overlapping_custom_groups: bool = False
    schema_version: int = 2

    def get_aa_class_enrichment(self, aa_class: str) -> float | None:
        """Get mean coverage enrichment for an AA class."""
        entry = self.aa_class_coverage.get_entry(aa_class)
        return entry.mean_coverage_enrichment if entry else None

    def get_custom_group_enrichment(self, group_name: str) -> float | None:
        """Get mean coverage enrichment for a custom group (vs rest_of_protein)."""
        if group_name not in self.custom_group_coverages:
            return None
        partition = self.custom_group_coverages[group_name]
        entry = partition.get_entry(group_name)
        return entry.mean_coverage_enrichment if entry else None

    def aa_class_enrichment_dict(self) -> dict[str, float]:
        """Get AA class mean enrichments as dict: {aa_class: enrichment}."""
        return self.aa_class_coverage.coverage_enrichment_dict()

    def custom_group_enrichment_dict(self) -> dict[str, float]:
        """Get custom group mean enrichments as dict: {group_name: enrichment}."""
        result = {}
        for group_name, partition in self.custom_group_coverages.items():
            entry = partition.get_entry(group_name)
            if entry and entry.mean_coverage_enrichment is not None:
                result[group_name] = entry.mean_coverage_enrichment
            else:
                result[group_name] = 0.0
        return result

    def aa_class_names(self) -> list[str]:
        """Get list of AA class names in canonical order."""
        canonical_order = ["aromatic", "polar", "nonpolar", "charged_positive", "charged_negative"]
        names = self.aa_class_coverage.element_names()
        return [n for n in canonical_order if n in names]

    def custom_group_names(self) -> list[str]:
        """Get list of custom group names."""
        return sorted(self.custom_group_coverages.keys())

    def save(self, path: str | Path) -> None:
        """Save to JSON file."""
        Path(path).write_text(json.dumps(self.model_dump(), indent=2))
        logger.info(f"Saved aggregated system coverage to {path}")

    @classmethod
    def load(cls, path: str | Path) -> "AggregatedSystemCoverageResult":
        """Load from JSON file."""
        data = json.loads(Path(path).read_text())
        return cls.model_validate(data)


def _aggregate_partition_coverage(
    partitions: list[PartitionCoverageResult],
) -> AggregatedPartitionCoverageResult:
    """Aggregate a partition's coverage across replicates.

    Parameters
    ----------
    partitions : list[PartitionCoverageResult]
        Same partition from multiple replicates

    Returns
    -------
    AggregatedPartitionCoverageResult
        Aggregated partition coverage
    """
    if not partitions:
        raise ValueError("No partitions to aggregate")

    # Helper function for computing mean and SEM
    def _compute_stats(values: list[float]) -> tuple[float | None, float | None]:
        """Compute mean and SEM from a list of values."""
        n = len(values)
        if n == 0:
            return None, None
        mean_val = float(np.mean(values))
        sem_val = float(np.std(values, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        return mean_val, sem_val

    # Get partition metadata from first result
    first_partition = partitions[0]

    # Collect all element names
    all_elements: set[str] = set()
    for p in partitions:
        for e in p.entries:
            all_elements.add(e.partition_element)

    # Collect all polymer types
    all_polymer_types: set[str] = set()
    for p in partitions:
        for e in p.entries:
            all_polymer_types.update(e.polymer_contributions.keys())

    entries = []
    for element_name in sorted(all_elements):
        # Collect values from each replicate
        coverage_shares = []
        enrichments = []
        polymer_contributions_all: dict[str, list[float]] = {pt: [] for pt in all_polymer_types}

        for p in partitions:
            entry = p.get_entry(element_name)
            if entry is not None:
                coverage_shares.append(entry.coverage_share)
                if entry.coverage_enrichment is not None:
                    enrichments.append(entry.coverage_enrichment)
                for poly_type in all_polymer_types:
                    contrib = entry.polymer_contributions.get(poly_type, 0.0)
                    polymer_contributions_all[poly_type].append(contrib)

        # Compute statistics
        mean_coverage_share, sem_coverage_share = _compute_stats(coverage_shares)
        mean_enrichment, sem_enrichment = _compute_stats(enrichments)

        # Compute mean polymer contributions
        mean_polymer_contributions = {}
        for poly_type, contribs in polymer_contributions_all.items():
            if contribs:
                mean_polymer_contributions[poly_type] = float(np.mean(contribs))
            else:
                mean_polymer_contributions[poly_type] = 0.0

        # Get metadata from first partition
        first_entry = first_partition.get_entry(element_name)
        n_exposed = first_entry.n_exposed_in_element if first_entry else 0
        n_total = first_entry.n_residues_in_element if first_entry else 0
        expected_share = first_entry.expected_share if first_entry else 0.0

        entries.append(
            AggregatedPartitionCoverageEntry(
                partition_element=element_name,
                mean_coverage_share=mean_coverage_share if mean_coverage_share else 0.0,
                sem_coverage_share=sem_coverage_share if sem_coverage_share else 0.0,
                mean_coverage_enrichment=mean_enrichment,
                sem_coverage_enrichment=sem_enrichment,
                per_replicate_enrichments=enrichments,
                expected_share=expected_share,
                n_exposed_in_element=n_exposed,
                n_residues_in_element=n_total,
                n_replicates=len(enrichments),
                mean_polymer_contributions=mean_polymer_contributions,
            )
        )

    return AggregatedPartitionCoverageResult(
        partition_name=first_partition.partition_name,
        partition_type=first_partition.partition_type,
        entries=entries,
        n_replicates=len(partitions),
    )


def aggregate_system_coverage(
    results: list[SystemCoverageResult],
) -> AggregatedSystemCoverageResult:
    """Aggregate system coverage across replicates.

    Computes mean ± SEM for coverage metrics across multiple replicates
    for all partitions (AA class, custom groups, combined).

    Parameters
    ----------
    results : list[SystemCoverageResult]
        System coverage results from multiple replicates

    Returns
    -------
    AggregatedSystemCoverageResult
        Aggregated results with mean and SEM for all partitions
    """
    if not results:
        raise ValueError("No results to aggregate")

    # Aggregate AA class partition
    aa_class_partitions = [r.aa_class_coverage for r in results]
    aggregated_aa_class = _aggregate_partition_coverage(aa_class_partitions)

    # Aggregate custom group partitions
    # First, find all custom group names across all results
    all_custom_groups: set[str] = set()
    for r in results:
        all_custom_groups.update(r.custom_group_coverages.keys())

    aggregated_custom_groups: dict[str, AggregatedPartitionCoverageResult] = {}
    for group_name in sorted(all_custom_groups):
        # Collect this partition from all results that have it
        group_partitions = [
            r.custom_group_coverages[group_name]
            for r in results
            if group_name in r.custom_group_coverages
        ]
        if group_partitions:
            aggregated_custom_groups[group_name] = _aggregate_partition_coverage(group_partitions)

    # Aggregate combined custom partition (if all results have it)
    aggregated_combined: AggregatedPartitionCoverageResult | None = None
    combined_partitions = [
        r.combined_custom_coverage for r in results if r.combined_custom_coverage is not None
    ]
    if combined_partitions and len(combined_partitions) == len(results):
        aggregated_combined = _aggregate_partition_coverage(combined_partitions)

    # Collect polymer types
    all_polymer_types: set[str] = set()
    for r in results:
        all_polymer_types.update(r.polymer_types_included)

    # Check for overlaps
    has_overlaps = any(r.has_overlapping_custom_groups for r in results)

    return AggregatedSystemCoverageResult(
        aa_class_coverage=aggregated_aa_class,
        custom_group_coverages=aggregated_custom_groups,
        combined_custom_coverage=aggregated_combined,
        n_replicates=len(results),
        total_exposed_residues=results[0].total_exposed_residues if results else 0,
        surface_exposure_threshold=results[0].surface_exposure_threshold if results else None,
        custom_group_selections=results[0].custom_group_selections if results else {},
        polymer_types_included=sorted(all_polymer_types),
        has_overlapping_custom_groups=has_overlaps,
    )


# =============================================================================
# Selection Resolution Helpers
# =============================================================================


def resolve_protein_group_selections(
    universe: "Universe",
    protein_group_selections: dict[str, str] | None = None,
) -> dict[str, set[int]]:
    """Resolve protein group MDAnalysis selections to residue IDs.

    Converts MDAnalysis selection strings into sets of residue IDs at
    analysis time. This allows the user to define groups with flexible
    selections while we work with concrete residue IDs internally.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe with loaded topology
    protein_group_selections : dict[str, str], optional
        Mapping of group name to MDAnalysis selection string.
        If None, uses DEFAULT_AA_CLASS_SELECTIONS.

    Returns
    -------
    dict[str, set[int]]
        Mapping of group name to set of residue IDs

    Examples
    --------
    >>> selections = {"aromatic": "protein and resname PHE TRP TYR HIS"}
    >>> groups = resolve_protein_group_selections(universe, selections)
    >>> print(groups["aromatic"])
    {12, 45, 67, 89, ...}  # Set of aromatic residue IDs
    """
    if protein_group_selections is None:
        protein_group_selections = DEFAULT_AA_CLASS_SELECTIONS.copy()

    resolved: dict[str, set[int]] = {}

    for group_name, selection in protein_group_selections.items():
        try:
            atoms = universe.select_atoms(selection)
            resids = set(atoms.residues.resids)
            resolved[group_name] = resids
            logger.debug(
                f"Resolved protein group '{group_name}': {len(resids)} residues ({selection})"
            )
        except Exception as e:
            logger.warning(
                f"Failed to resolve protein group '{group_name}' with selection '{selection}': {e}"
            )
            resolved[group_name] = set()

    return resolved


def resolve_protein_groups_from_surface_exposure(
    surface_exposure: "SurfaceExposureResult",
    include_default_aa_groups: bool = True,
    custom_protein_groups: dict[str, list[int]] | None = None,
) -> dict[str, set[int]]:
    """Resolve protein groups from surface exposure data without Universe.

    This function derives protein group → residue ID mappings using only the
    surface exposure result (which already contains resid, resname, and aa_class
    for each residue). This allows binding preference computation without
    requiring an MDAnalysis Universe at comparison time.

    The function supports:
    - Default AA class groups (aromatic, polar, nonpolar, charged_+, charged_-)
      derived from the aa_class field in surface exposure data
    - Custom user-defined groups specified as resid lists
    - Override behavior: if a custom group has the same name as a default,
      the custom definition takes precedence

    Parameters
    ----------
    surface_exposure : SurfaceExposureResult
        Surface exposure analysis result containing residue data
    include_default_aa_groups : bool, default True
        If True, include default AA class groupings (aromatic, polar, etc.)
        derived from surface exposure data
    custom_protein_groups : dict[str, list[int]], optional
        User-defined protein groups as {group_name: [resid1, resid2, ...]}.
        If a group name matches a default AA class, it overrides that default.

    Returns
    -------
    dict[str, set[int]]
        Mapping of group name to set of residue IDs

    Examples
    --------
    >>> from polyzymd.analysis.contacts.surface_exposure import SurfaceExposureFilter
    >>> filter = SurfaceExposureFilter(threshold=0.2)
    >>> surface_result = filter.calculate("enzyme.pdb")
    >>> # Get default AA groups + custom active_site group
    >>> groups = resolve_protein_groups_from_surface_exposure(
    ...     surface_result,
    ...     include_default_aa_groups=True,
    ...     custom_protein_groups={"active_site": [77, 133, 156]}
    ... )
    >>> print(groups.keys())
    dict_keys(['aromatic', 'polar', 'nonpolar', 'charged_positive', 'charged_negative', 'active_site'])
    """
    resolved: dict[str, set[int]] = {}
    all_valid_resids = surface_exposure.all_resids

    # Step 1: Build default AA class groups from surface exposure data
    if include_default_aa_groups:
        # Group residues by their aa_class from surface exposure
        aa_class_groups: dict[str, set[int]] = {}
        for res_exp in surface_exposure.residue_exposures:
            aa_class = res_exp.aa_class
            if aa_class and aa_class != "unknown":
                if aa_class not in aa_class_groups:
                    aa_class_groups[aa_class] = set()
                aa_class_groups[aa_class].add(res_exp.resid)

        for group_name, resids in aa_class_groups.items():
            resolved[group_name] = resids
            logger.debug(
                f"Default AA group '{group_name}': {len(resids)} residues from surface exposure"
            )

    # Step 2: Add custom protein groups (with override behavior)
    if custom_protein_groups:
        for group_name, resid_list in custom_protein_groups.items():
            # Convert to set and validate
            requested_resids = set(resid_list)
            valid_resids = requested_resids & all_valid_resids
            invalid_resids = requested_resids - all_valid_resids

            # Warn about invalid resids
            if invalid_resids:
                logger.warning(
                    f"Custom protein group '{group_name}': {len(invalid_resids)} resids not found "
                    f"in enzyme structure and will be ignored: {sorted(invalid_resids)}"
                )

            # Override if same name as default (user intent takes precedence)
            if group_name in resolved:
                logger.info(f"Custom protein group '{group_name}' overrides default AA class group")

            resolved[group_name] = valid_resids
            logger.debug(f"Custom protein group '{group_name}': {len(valid_resids)} valid residues")

    return resolved


def resolve_polymer_type_selections(
    universe: "Universe",
    polymer_type_selections: dict[str, str] | None = None,
) -> list[str]:
    """Resolve polymer type selections and return list of polymer types.

    If explicit selections are provided, validates them and returns the keys.
    If None, auto-detects polymer types from the standard polymer chain (C).

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe with loaded topology
    polymer_type_selections : dict[str, str], optional
        Mapping of type name to MDAnalysis selection string.
        If None, auto-detects from "chainID C".

    Returns
    -------
    list[str]
        List of polymer type names (resnames or selection keys)

    Examples
    --------
    >>> # Auto-detect from chain C
    >>> polymer_types = resolve_polymer_type_selections(universe, None)
    >>> print(polymer_types)
    ['SBM', 'EGM']

    >>> # Explicit selections
    >>> selections = {"SBMA": "chainID C and resname SBM"}
    >>> polymer_types = resolve_polymer_type_selections(universe, selections)
    >>> print(polymer_types)
    ['SBMA']
    """
    if polymer_type_selections is not None:
        # Validate selections and return keys
        valid_types = []
        for type_name, selection in polymer_type_selections.items():
            try:
                atoms = universe.select_atoms(selection)
                if len(atoms) > 0:
                    valid_types.append(type_name)
                    logger.debug(f"Polymer type '{type_name}': {len(atoms)} atoms ({selection})")
                else:
                    logger.warning(
                        f"Polymer type '{type_name}' selection matched no atoms: {selection}"
                    )
            except Exception as e:
                logger.warning(f"Invalid selection for polymer type '{type_name}': {e}")
        return valid_types

    # Auto-detect from chain C
    try:
        polymer_atoms = universe.select_atoms("chainID C")
        if len(polymer_atoms) == 0:
            logger.warning("No atoms found in chainID C for polymer auto-detection")
            return []

        # Get unique resnames
        resnames = set(polymer_atoms.residues.resnames)
        logger.debug(f"Auto-detected polymer types from chain C: {resnames}")
        return sorted(resnames)

    except Exception as e:
        logger.warning(f"Failed to auto-detect polymer types: {e}")
        return []


def extract_polymer_composition(
    universe: "Universe",
    polymer_type_selections: dict[str, str] | None = None,
) -> PolymerComposition:
    """Extract polymer composition (residue and heavy atom counts) from trajectory.

    This function counts the number of residues and heavy atoms for each
    polymer type in the system. The composition data is used for dual
    normalization of binding preference enrichment ratios.

    Parameters
    ----------
    universe : Universe
        MDAnalysis Universe with loaded topology
    polymer_type_selections : dict[str, str], optional
        Mapping of type name to MDAnalysis selection string.
        If None, auto-detects from "chainID C" using unique resnames.

    Returns
    -------
    PolymerComposition
        Composition data with residue_counts and heavy_atom_counts per type.
        Heavy atoms are defined as all atoms with element != 'H'.

    Notes
    -----
    For auto-detection (polymer_type_selections=None), the function:
    1. Selects all atoms in chainID C
    2. Groups residues by resname
    3. Counts residues and heavy atoms per resname

    For explicit selections, each selection string defines which atoms
    belong to that polymer type. The function counts residues (unique
    resids) and heavy atoms within each selection.

    Examples
    --------
    >>> # Auto-detect from chain C
    >>> composition = extract_polymer_composition(universe)
    >>> print(composition.residue_counts)
    {'SBM': 50, 'EGM': 50}
    >>> print(composition.heavy_atom_counts)
    {'SBM': 750, 'EGM': 400}

    >>> # Explicit selections
    >>> selections = {"SBMA": "chainID C and resname SBM"}
    >>> composition = extract_polymer_composition(universe, selections)
    """
    residue_counts: dict[str, int] = {}
    heavy_atom_counts: dict[str, int] = {}

    if polymer_type_selections is not None:
        # Use explicit selections
        for type_name, selection in polymer_type_selections.items():
            try:
                atoms = universe.select_atoms(selection)
                if len(atoms) == 0:
                    logger.warning(
                        f"Polymer type '{type_name}' selection matched no atoms: {selection}"
                    )
                    continue

                # Count unique residues
                n_residues = len(atoms.residues)
                residue_counts[type_name] = n_residues

                # Count heavy atoms (non-hydrogen)
                # MDAnalysis stores element info - filter out hydrogens
                try:
                    heavy_atoms = atoms.select_atoms("not element H")
                    n_heavy = len(heavy_atoms)
                except Exception:
                    # Fallback: count atoms with mass > 1.1 (heavier than H)
                    n_heavy = sum(1 for a in atoms if a.mass > 1.1)

                heavy_atom_counts[type_name] = n_heavy

                logger.debug(
                    f"Polymer type '{type_name}': {n_residues} residues, {n_heavy} heavy atoms"
                )

            except Exception as e:
                logger.warning(f"Failed to extract composition for '{type_name}': {e}")

    else:
        # Auto-detect from chain C
        try:
            polymer_atoms = universe.select_atoms("chainID C")
            if len(polymer_atoms) == 0:
                logger.warning("No atoms found in chainID C for polymer composition")
                return PolymerComposition()

            # Group by resname
            resnames = set(polymer_atoms.residues.resnames)

            for resname in resnames:
                # Select atoms of this resname within chain C
                type_atoms = universe.select_atoms(f"chainID C and resname {resname}")

                # Count residues
                n_residues = len(type_atoms.residues)
                residue_counts[resname] = n_residues

                # Count heavy atoms
                try:
                    heavy_atoms = type_atoms.select_atoms("not element H")
                    n_heavy = len(heavy_atoms)
                except Exception:
                    # Fallback: count atoms with mass > 1.1
                    n_heavy = sum(1 for a in type_atoms if a.mass > 1.1)

                heavy_atom_counts[resname] = n_heavy

                logger.debug(
                    f"Auto-detected polymer '{resname}': {n_residues} residues, "
                    f"{n_heavy} heavy atoms"
                )

        except Exception as e:
            logger.warning(f"Failed to extract polymer composition: {e}")
            return PolymerComposition()

    composition = PolymerComposition(
        residue_counts=residue_counts,
        heavy_atom_counts=heavy_atom_counts,
    )

    logger.info(
        f"Polymer composition extracted: {composition.total_residues} total residues, "
        f"{composition.total_heavy_atoms} total heavy atoms across {len(residue_counts)} types"
    )

    return composition


def compute_binding_preference_from_config(
    contact_result: "ContactResult",
    universe: "Universe",
    enzyme_pdb_path: Path | str,
    config: "ContactsConfig",
) -> BindingPreferenceResult:
    """Compute binding preference using ContactsConfig settings.

    This is a convenience function that orchestrates the full binding
    preference computation using configuration from analysis.yaml.
    It:
    1. Calculates surface exposure from enzyme PDB
    2. Resolves protein group selections to residue IDs
    3. Computes binding preference with enrichment ratios

    Parameters
    ----------
    contact_result : ContactResult
        Contact analysis results from trajectory
    universe : Universe
        MDAnalysis Universe (used for resolving selections)
    enzyme_pdb_path : Path or str
        Path to enzyme PDB file for SASA calculation
    config : ContactsConfig
        Configuration from analysis.yaml with binding preference settings

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
    >>> from polyzymd.analysis.config import ContactsConfig
    >>> config = ContactsConfig(
    ...     compute_binding_preference=True,
    ...     surface_exposure_threshold=0.2,
    ... )
    >>> result = compute_binding_preference_from_config(
    ...     contact_result, universe, "enzyme.pdb", config
    ... )
    >>> print(result.enrichment_matrix())
    {'SBM': {'aromatic': 1.45, 'polar': 0.82, ...}, ...}
    """
    from polyzymd.analysis.contacts.surface_exposure import SurfaceExposureFilter

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
    polymer_types = resolve_polymer_type_selections(universe, config.polymer_type_selections)
    logger.info(f"Polymer types for binding preference: {polymer_types}")

    # Step 4: Extract polymer composition for dual normalization
    polymer_composition = extract_polymer_composition(universe, config.polymer_type_selections)

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
