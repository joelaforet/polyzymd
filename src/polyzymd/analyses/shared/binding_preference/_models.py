"""Binding preference result models and data containers."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field

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

    entries: list[BindingPreferenceEntry] = Field(
        default_factory=list,
        description="DEPRECATED: Overlapping-groups entries. Use binding_preference instead.",
    )
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
    binding_preference: "PolymerBindingPreferenceResult | None" = Field(
        default=None,
        description=(
            "Partition-based per-polymer binding preference. "
            "contact_share sums to 1.0 within each partition for each polymer. "
            "This is the primary binding preference output (v5+)."
        ),
    )
    metadata: dict[str, Any] = Field(default_factory=dict)
    schema_version: int = 5  # Version 5: adds partition-based binding_preference

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
        return sorted({e.polymer_type for e in self.entries})

    def protein_groups(self) -> list[str]:
        """Get list of protein groups in this result."""
        return sorted({e.protein_group for e in self.entries})

    def save(self, path: str | Path) -> None:
        """Save to JSON file."""
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
    binding_preference: "AggregatedPolymerBindingPreferenceResult | None" = Field(
        default=None,
        description=(
            "Aggregated partition-based per-polymer binding preference. "
            "contact_share sums to 1.0 within each partition for each polymer. "
            "This is the primary binding preference output (v5+)."
        ),
    )
    schema_version: int = 5  # Version 5: adds partition-based binding_preference

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
        return sorted({e.polymer_type for e in self.entries})

    def protein_groups(self) -> list[str]:
        """Get list of protein groups."""
        return sorted({e.protein_group for e in self.entries})

    def save(self, path: str | Path) -> None:
        """Save to JSON file."""
        Path(path).write_text(json.dumps(self.model_dump(), indent=2))
        logger.info(f"Saved aggregated binding preference to {path}")

    @classmethod
    def load(cls, path: str | Path) -> "AggregatedBindingPreferenceResult":
        """Load from JSON file."""
        data = json.loads(Path(path).read_text())
        return cls.model_validate(data)


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
    partition_type: Literal["aa_class", "binary_custom", "combined_custom", "user_defined"]
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

    4. **User-Defined Partitions** (from protein_partitions config):
       Custom partitions specified by the user in the YAML config. Each
       partition references groups from protein_groups and must be mutually
       exclusive. 'rest_of_protein' is auto-added if the groups don't cover
       all exposed protein residues. One plot per partition is generated.

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
    user_defined_partitions : dict[str, PartitionCoverageResult]
        User-defined partitions from protein_partitions config.
        Keys are partition names, values are the computed coverage partitions.
        'rest_of_protein' is auto-added if groups don't fully cover the protein.
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
    user_defined_partitions: dict[str, PartitionCoverageResult] = Field(
        default_factory=dict,
        description=(
            "User-defined partitions from protein_partitions config. "
            "Each partition contains mutually exclusive groups defined by the user, "
            "with 'rest_of_protein' auto-added if groups don't cover all protein residues."
        ),
    )

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

    def user_partition_names(self) -> list[str]:
        """Get list of user-defined partition names."""
        return sorted(self.user_defined_partitions.keys())

    def get_user_partition(self, partition_name: str) -> PartitionCoverageResult | None:
        """Get a user-defined partition by name.

        Parameters
        ----------
        partition_name : str
            Name of the partition (e.g., "lid_helices")

        Returns
        -------
        PartitionCoverageResult | None
            The partition result, or None if not found
        """
        return self.user_defined_partitions.get(partition_name)

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


class PartitionBindingEntry(BaseModel):
    """Binding metrics for one partition element for a specific polymer type.

    A partition element is a mutually exclusive subset of protein residues.
    Within a partition, all elements together cover the entire protein surface
    exactly once (no residue is counted in multiple elements).

    This entry is for a SINGLE polymer type, answering:
    "What fraction of SBMA's contacts go to aromatic residues?"

    This ensures that:
    - contact_share sums to 1.0 across all elements in the partition (for this polymer)
    - expected_share sums to 1.0 across all elements in the partition
    - enrichment is mathematically valid (no inflated denominators)

    Attributes
    ----------
    partition_element : str
        Name of this partition element (e.g., "aromatic", "lid_helix_5", "rest_of_protein")
    polymer_type : str
        Polymer type this entry is for (e.g., "SBM", "EGM")
    total_contact_frames : int
        Sum of contact frames from THIS polymer type to residues in this element
    contact_share : float
        Fraction of this polymer's contacts that went to this element.
        Sums to 1.0 across all elements in the partition (for this polymer).
    expected_share : float
        Expected share based on surface availability (n_exposed / total_exposed).
        Sums to 1.0 across all elements in the partition.
    enrichment : float | None
        Zero-centered enrichment: (contact_share / expected_share) - 1
    n_exposed_in_element : int
        Number of surface-exposed residues in this element
    n_residues_in_element : int
        Total residues in this element (exposed + buried)
    n_residues_contacted : int
        Number of exposed residues that had at least one contact from this polymer
    """

    partition_element: str
    polymer_type: str
    total_contact_frames: int = Field(
        default=0,
        description="Sum of contact frames from THIS polymer type to this element",
    )
    contact_share: float = Field(
        default=0.0,
        description="Fraction of this polymer's contacts that went to this element",
    )
    expected_share: float = Field(
        default=0.0,
        description="Expected share based on protein surface availability",
    )
    enrichment: float | None = Field(
        default=None,
        description="Zero-centered enrichment: (contact_share / expected_share) - 1",
    )
    n_exposed_in_element: int = Field(
        default=0,
        description="Surface-exposed residues in this element",
    )
    n_residues_in_element: int = Field(
        default=0,
        description="Total residues in this element",
    )
    n_residues_contacted: int = Field(
        default=0,
        description="Exposed residues contacted by this polymer type",
    )


class PartitionBindingResult(BaseModel):
    """Binding preference for a complete partition for ONE polymer type.

    A partition divides the protein surface into mutually exclusive regions
    that together cover the entire surface. This class stores the binding
    preference of a single polymer type across all partition elements.

    This ensures valid enrichment calculations where both contact_share
    and expected_share sum to 1.0.

    Partition Types
    ---------------
    - **aa_class**: 5-way partition by amino acid class
      (aromatic, polar, nonpolar, charged_positive, charged_negative)
    - **user_defined**: N+1 way partition with user-specified groups
      plus "rest_of_protein" (auto-added if groups don't cover all residues)

    Attributes
    ----------
    partition_name : str
        Descriptive name (e.g., "aa_class", "lid_helices")
    partition_type : str
        One of: "aa_class", "user_defined"
    polymer_type : str
        Polymer type this result is for (e.g., "SBM", "EGM")
    entries : list[PartitionBindingEntry]
        Binding metrics for each element in the partition
    total_contact_share : float
        Validation check: should be ~1.0
    total_expected_share : float
        Validation check: should be ~1.0
    total_contact_frames : int
        Total contact frames from this polymer type (across all elements)
    """

    partition_name: str
    partition_type: Literal["aa_class", "user_defined"]
    polymer_type: str
    entries: list[PartitionBindingEntry] = Field(default_factory=list)
    total_contact_share: float = Field(
        default=1.0,
        description="Sum of contact_share across elements (validation: should be ~1.0)",
    )
    total_expected_share: float = Field(
        default=1.0,
        description="Sum of expected_share across elements (validation: should be ~1.0)",
    )
    total_contact_frames: int = Field(
        default=0,
        description="Total contact frames from this polymer type",
    )

    def to_dataframe(self) -> "pd.DataFrame":
        """Convert to pandas DataFrame for analysis/plotting."""
        import pandas as pd

        return pd.DataFrame([e.model_dump() for e in self.entries])

    def enrichment_dict(self) -> dict[str, float]:
        """Get enrichment as dict: {element: enrichment}."""
        return {
            e.partition_element: (e.enrichment if e.enrichment is not None else 0.0)
            for e in self.entries
        }

    def contact_share_dict(self) -> dict[str, float]:
        """Get contact shares as dict: {element: share}."""
        return {e.partition_element: e.contact_share for e in self.entries}

    def expected_share_dict(self) -> dict[str, float]:
        """Get expected shares as dict: {element: share}."""
        return {e.partition_element: e.expected_share for e in self.entries}

    def get_entry(self, element_name: str) -> PartitionBindingEntry | None:
        """Get the entry for a specific partition element."""
        for entry in self.entries:
            if entry.partition_element == element_name:
                return entry
        return None

    def element_names(self) -> list[str]:
        """Get list of partition element names."""
        return [e.partition_element for e in self.entries]


class PolymerBindingPreferenceResult(BaseModel):
    """Per-polymer binding preference using proper partition structure.

    This result stores binding preference for ALL polymer types, with each
    polymer having its own partition-based enrichment calculations.

    Unlike SystemCoverageResult (which collapses all polymer contacts), this
    maintains per-polymer data to answer: "Does SBMA prefer aromatic residues
    more than EGMA does?"

    Partition Strategy (per polymer type)
    -------------------------------------
    1. **AA Class Partition** (always computed):
       5-way partition by amino acid class. Every surface residue belongs
       to exactly one class. Each polymer type gets its own enrichment values.

    2. **User-Defined Partitions** (from protein_partitions config):
       Custom partitions specified by the user. Each partition references groups
       from protein_groups. 'rest_of_protein' is auto-added if groups don't
       cover all exposed protein residues. Each polymer type gets its own
       enrichment values per partition.

    Attributes
    ----------
    aa_class_binding : dict[str, PartitionBindingResult]
        AA class partition binding for each polymer type.
        Keys are polymer type names (e.g., "SBM", "EGM").
    user_defined_partitions : dict[str, dict[str, PartitionBindingResult]]
        User-defined partitions for each polymer type.
        Outer keys are partition names, inner keys are polymer types.
        Example: {"lid_helices": {"SBM": ..., "EGM": ...}}
    n_frames : int
        Total frames analyzed
    total_exposed_residues : int
        Number of surface-exposed protein residues
    surface_exposure_threshold : float | None
        SASA threshold used for surface filtering
    polymer_types : list[str]
        Polymer types included in this result
    polymer_composition : PolymerComposition | None
        Polymer composition metadata
    schema_version : int
        Schema version (5 = partition-based binding preference)
    """

    aa_class_binding: dict[str, PartitionBindingResult] = Field(
        default_factory=dict,
        description="AA class partition binding for each polymer type",
    )
    user_defined_partitions: dict[str, dict[str, PartitionBindingResult]] = Field(
        default_factory=dict,
        description=(
            "User-defined partitions for each polymer type. "
            "Outer key: partition name, inner key: polymer type."
        ),
    )

    # Metadata
    n_frames: int = 0
    total_exposed_residues: int = 0
    surface_exposure_threshold: float | None = None
    polymer_types: list[str] = Field(default_factory=list)
    polymer_composition: PolymerComposition | None = None
    protein_groups_used: dict[str, str] = Field(default_factory=dict)
    schema_version: int = 5  # Version 5: partition-based per-polymer binding

    def get_aa_class_enrichment(self, polymer_type: str, aa_class: str) -> float | None:
        """Get binding enrichment for an AA class for a specific polymer.

        Parameters
        ----------
        polymer_type : str
            Polymer type (e.g., "SBM")
        aa_class : str
            One of: aromatic, polar, nonpolar, charged_positive, charged_negative

        Returns
        -------
        float | None
            Binding enrichment, or None if not found
        """
        if polymer_type not in self.aa_class_binding:
            return None
        entry = self.aa_class_binding[polymer_type].get_entry(aa_class)
        return entry.enrichment if entry else None

    def get_user_partition_enrichment(
        self, partition_name: str, polymer_type: str, element_name: str
    ) -> float | None:
        """Get binding enrichment for a user partition element for a specific polymer.

        Parameters
        ----------
        partition_name : str
            Name of the user-defined partition (e.g., "lid_helices")
        polymer_type : str
            Polymer type (e.g., "SBM")
        element_name : str
            Element within the partition (e.g., "lid_helix_5")

        Returns
        -------
        float | None
            Binding enrichment, or None if not found
        """
        if partition_name not in self.user_defined_partitions:
            return None
        if polymer_type not in self.user_defined_partitions[partition_name]:
            return None
        entry = self.user_defined_partitions[partition_name][polymer_type].get_entry(element_name)
        return entry.enrichment if entry else None

    def aa_class_enrichment_matrix(self) -> dict[str, dict[str, float]]:
        """Get AA class enrichments as nested dict: {polymer_type: {aa_class: enrichment}}.

        Returns
        -------
        dict[str, dict[str, float]]
            Nested mapping of enrichment values.
        """
        result: dict[str, dict[str, float]] = {}
        for poly_type, partition_result in self.aa_class_binding.items():
            result[poly_type] = partition_result.enrichment_dict()
        return result

    def user_partition_enrichment_matrix(self, partition_name: str) -> dict[str, dict[str, float]]:
        """Get user partition enrichments as nested dict.

        Parameters
        ----------
        partition_name : str
            Name of the user-defined partition

        Returns
        -------
        dict[str, dict[str, float]]
            {polymer_type: {element_name: enrichment}}
        """
        if partition_name not in self.user_defined_partitions:
            return {}
        result: dict[str, dict[str, float]] = {}
        for poly_type, partition_result in self.user_defined_partitions[partition_name].items():
            result[poly_type] = partition_result.enrichment_dict()
        return result

    def aa_class_names(self) -> list[str]:
        """Get list of AA class names in canonical order."""
        canonical_order = ["aromatic", "polar", "nonpolar", "charged_positive", "charged_negative"]
        if not self.aa_class_binding:
            return []
        # Get from first polymer type
        first_poly = next(iter(self.aa_class_binding.values()))
        names = first_poly.element_names()
        return [n for n in canonical_order if n in names]

    def user_partition_names(self) -> list[str]:
        """Get list of user-defined partition names."""
        return sorted(self.user_defined_partitions.keys())

    def save(self, path: str | Path) -> None:
        """Save to JSON file."""
        Path(path).write_text(json.dumps(self.model_dump(), indent=2))
        logger.info(f"Saved polymer binding preference result to {path}")

    @classmethod
    def load(cls, path: str | Path) -> "PolymerBindingPreferenceResult":
        """Load from JSON file."""
        data = json.loads(Path(path).read_text())
        return cls.model_validate(data)


class AggregatedPartitionBindingEntry(BaseModel):
    """Aggregated binding metrics for one partition element for a specific polymer type.

    Contains mean ± SEM across replicates for binding preference, enabling
    statistical comparison of binding enrichment across conditions.

    Attributes
    ----------
    partition_element : str
        Element name (e.g., "aromatic", "lid_helix_5", "rest_of_protein")
    polymer_type : str
        Polymer type this entry is for (e.g., "SBM", "EGM")
    mean_contact_share : float
        Mean contact share across replicates
    sem_contact_share : float
        Standard error of contact share
    mean_enrichment : float | None
        Mean enrichment across replicates
    sem_enrichment : float | None
        Standard error of enrichment
    per_replicate_enrichments : list[float]
        Enrichment values from each replicate
    expected_share : float
        Expected share based on surface availability
    n_exposed_in_element : int
        Surface-exposed residues in this element
    n_residues_in_element : int
        Total residues in this element
    n_replicates : int
        Number of replicates with valid data
    """

    partition_element: str
    polymer_type: str
    mean_contact_share: float = Field(
        default=0.0,
        description="Mean contact share across replicates",
    )
    sem_contact_share: float = Field(
        default=0.0,
        description="Standard error of contact share",
    )
    mean_enrichment: float | None = Field(
        default=None,
        description="Mean enrichment across replicates",
    )
    sem_enrichment: float | None = Field(
        default=None,
        description="Standard error of enrichment",
    )
    per_replicate_enrichments: list[float] = Field(
        default_factory=list,
        description="Enrichment values from each replicate",
    )
    expected_share: float = Field(
        default=0.0,
        description="Expected share based on surface availability",
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


class AggregatedPartitionBindingResult(BaseModel):
    """Aggregated binding preference for a partition for ONE polymer type.

    Contains aggregated statistics across replicates for all partition elements.

    Attributes
    ----------
    partition_name : str
        Descriptive name (e.g., "aa_class", "lid_helices")
    partition_type : str
        One of: "aa_class", "user_defined"
    polymer_type : str
        Polymer type this result is for
    entries : list[AggregatedPartitionBindingEntry]
        Aggregated binding metrics for each element
    mean_total_contact_share : float
        Mean of total_contact_share across replicates (validation: should be ~1.0)
    n_replicates : int
        Number of replicates
    """

    partition_name: str
    partition_type: Literal["aa_class", "user_defined"]
    polymer_type: str
    entries: list[AggregatedPartitionBindingEntry] = Field(default_factory=list)
    mean_total_contact_share: float = Field(
        default=1.0,
        description="Mean sum of contact_share across elements (should be ~1.0)",
    )
    n_replicates: int = Field(default=0)

    def enrichment_dict(self) -> dict[str, float]:
        """Get mean enrichment as dict: {element: mean_enrichment}."""
        return {
            e.partition_element: (e.mean_enrichment if e.mean_enrichment is not None else 0.0)
            for e in self.entries
        }

    def element_names(self) -> list[str]:
        """Get list of partition element names."""
        return [e.partition_element for e in self.entries]


class AggregatedPolymerBindingPreferenceResult(BaseModel):
    """Aggregated per-polymer binding preference across replicates.

    Contains mean ± SEM for all partition-based binding metrics.

    Attributes
    ----------
    aa_class_binding : dict[str, AggregatedPartitionBindingResult]
        Aggregated AA class partition binding for each polymer type.
    user_defined_partitions : dict[str, dict[str, AggregatedPartitionBindingResult]]
        Aggregated user-defined partitions for each polymer type.
    n_replicates : int
        Number of replicates
    total_exposed_residues : int
        Number of surface-exposed protein residues
    surface_exposure_threshold : float | None
        SASA threshold used
    polymer_types : list[str]
        Polymer types included
    schema_version : int
        Schema version
    """

    aa_class_binding: dict[str, AggregatedPartitionBindingResult] = Field(
        default_factory=dict,
        description="Aggregated AA class partition binding for each polymer type",
    )
    user_defined_partitions: dict[str, dict[str, AggregatedPartitionBindingResult]] = Field(
        default_factory=dict,
        description="Aggregated user-defined partitions for each polymer type",
    )

    n_replicates: int = 0
    total_exposed_residues: int = 0
    surface_exposure_threshold: float | None = None
    polymer_types: list[str] = Field(default_factory=list)
    schema_version: int = 5

    def aa_class_enrichment_matrix(self) -> dict[str, dict[str, float]]:
        """Get AA class enrichments as nested dict: {polymer_type: {aa_class: mean_enrichment}}."""
        result: dict[str, dict[str, float]] = {}
        for poly_type, partition_result in self.aa_class_binding.items():
            result[poly_type] = partition_result.enrichment_dict()
        return result

    def aa_class_names(self) -> list[str]:
        """Get list of AA class names in canonical order."""
        canonical_order = ["aromatic", "polar", "nonpolar", "charged_positive", "charged_negative"]
        if not self.aa_class_binding:
            return []
        first_poly = next(iter(self.aa_class_binding.values()))
        names = first_poly.element_names()
        return [n for n in canonical_order if n in names]

    def user_partition_names(self) -> list[str]:
        """Get list of user-defined partition names."""
        return sorted(self.user_defined_partitions.keys())


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
    partition_type: Literal["aa_class", "binary_custom", "combined_custom", "user_defined"]
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
    user_defined_partitions: dict[str, AggregatedPartitionCoverageResult] = Field(
        default_factory=dict,
        description=(
            "Aggregated user-defined partitions from protein_partitions config. "
            "Keys are partition names, values are aggregated coverage results."
        ),
    )

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


