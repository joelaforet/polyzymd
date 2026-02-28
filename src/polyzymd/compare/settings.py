"""Analysis and comparison settings for the comparison workflow.

This module defines the concrete settings classes for each analysis type,
registered via the AnalysisSettingsRegistry and ComparisonSettingsRegistry.

Analysis Settings (WHAT to analyze):
- RMSFAnalysisSettings: RMSF calculation parameters
- DistancesAnalysisSettings: Distance pair monitoring parameters
- CatalyticTriadAnalysisSettings: Active site distance analysis
- ContactsAnalysisSettings: Polymer-protein contact parameters

Comparison Settings (HOW to compare):
- RMSFComparisonSettings: (no comparison-specific params)
- DistancesComparisonSettings: (no comparison-specific params)
- CatalyticTriadComparisonSettings: (no comparison-specific params)
- ContactsComparisonSettings: FDR, effect size, top residues

All settings classes are auto-registered on module import.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional

from pydantic import Field, field_validator, model_validator

from polyzymd.analysis.core.constants import (
    DEFAULT_CONTACT_CUTOFF,
    DEFAULT_DISTANCE_THRESHOLD,
    DEFAULT_SURFACE_EXPOSURE_THRESHOLD,
)
from polyzymd.analysis.core.registry import (
    AnalysisSettingsRegistry,
    BaseAnalysisSettings,
    BaseComparisonSettings,
    ComparisonSettingsRegistry,
)

if TYPE_CHECKING:
    from polyzymd.analysis.core.alignment import AlignmentConfig

# ============================================================================
# RMSF Settings
# ============================================================================


@AnalysisSettingsRegistry.register("rmsf")
class RMSFAnalysisSettings(BaseAnalysisSettings):
    """RMSF analysis settings.

    Attributes
    ----------
    selection : str
        MDAnalysis selection string for RMSF calculation.
    reference_mode : str
        Reference structure mode: centroid, average, frame, or external.
    reference_frame : int, optional
        Frame number if reference_mode is 'frame' (1-indexed).
    reference_file : str, optional
        Path to external PDB file if reference_mode is 'external'.
    """

    selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection string for RMSF calculation",
    )
    reference_mode: str = Field(
        default="centroid",
        description="Reference structure mode: centroid, average, frame, or external",
    )
    reference_frame: Optional[int] = Field(
        default=None,
        description="Frame number if reference_mode is 'frame' (1-indexed)",
    )
    reference_file: Optional[str] = Field(
        default=None,
        description=(
            "Path to external PDB file if reference_mode is 'external'. "
            "The PDB must contain protein atoms matching the simulation topology."
        ),
    )

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "rmsf"

    @field_validator("reference_mode", mode="after")
    @classmethod
    def validate_reference_mode(cls, v: str) -> str:
        """Validate reference mode is one of the allowed values."""
        valid = {"centroid", "average", "frame", "external"}
        if v not in valid:
            raise ValueError(f"reference_mode must be one of {valid}, got '{v}'")
        return v

    @model_validator(mode="after")
    def validate_reference_params(self) -> "RMSFAnalysisSettings":
        """Validate reference_frame and reference_file for their modes."""
        if self.reference_mode == "frame" and self.reference_frame is None:
            raise ValueError("reference_frame is required when reference_mode is 'frame'")
        if self.reference_mode == "external" and self.reference_file is None:
            raise ValueError(
                "reference_file is required when reference_mode is 'external'. "
                "Provide a path to the external PDB reference structure."
            )
        return self

    def to_analysis_yaml_dict(self) -> dict[str, Any]:
        """Convert to analysis.yaml-compatible dictionary."""
        result = {
            "enabled": True,
            "selection": self.selection,
            "reference_mode": self.reference_mode,
        }
        if self.reference_frame is not None:
            result["reference_frame"] = self.reference_frame
        if self.reference_file is not None:
            result["reference_file"] = self.reference_file
        return result


@ComparisonSettingsRegistry.register("rmsf")
class RMSFComparisonSettings(BaseComparisonSettings):
    """Comparison settings for RMSF analysis.

    Currently empty — all RMSF comparison behavior uses defaults from
    ``BaseComparisonSettings``.  This class exists as an extension point:
    add fields here when RMSF-specific comparison parameters are needed
    (e.g., a per-residue significance threshold) without modifying the
    orchestrator or other comparison types.
    """

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "rmsf"


# ============================================================================
# Distance Analysis Settings
# ============================================================================


class DistancePairSettings(BaseAnalysisSettings):
    """Configuration for a single distance pair.

    Attributes
    ----------
    label : str
        Human-readable label for this pair.
    selection_a : str
        First atom/point selection.
    selection_b : str
        Second atom/point selection.
    threshold : float, optional
        Per-pair distance threshold (Angstroms). If None, uses the global
        threshold from DistancesAnalysisSettings.
    """

    label: str = Field(..., description="Human-readable label for this pair")
    selection_a: str = Field(..., description="First atom/point selection")
    selection_b: str = Field(..., description="Second atom/point selection")
    threshold: Optional[float] = Field(
        default=None,
        description="Per-pair distance threshold (Angstroms). If None, uses global threshold.",
    )

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "distance_pair"

    def to_analysis_yaml_dict(self) -> dict[str, Any]:
        """Convert to analysis.yaml-compatible dictionary."""
        result: dict[str, Any] = {
            "label": self.label,
            "selection_a": self.selection_a,
            "selection_b": self.selection_b,
        }
        if self.threshold is not None:
            result["threshold"] = self.threshold
        return result


@AnalysisSettingsRegistry.register("distances")
class DistancesAnalysisSettings(BaseAnalysisSettings):
    """Distance analysis settings.

    Attributes
    ----------
    threshold : float, optional
        Distance threshold for contact analysis (Angstroms).
    pairs : list[DistancePairSettings]
        List of atom pairs to measure distances between.
    use_pbc : bool
        Use PBC-aware minimum image distances. Default True.
    align_trajectory : bool
        Align trajectory before distance calculation. Default True.
        When enabled, removes rotational drift and COM motion that
        can add noise to inter-domain distance measurements.
    alignment_selection : str
        MDAnalysis selection for trajectory alignment.
        Default: "protein and name CA".
    alignment_mode : str
        Reference mode for alignment: "centroid", "average", or "frame".
        Default: "centroid".
    alignment_frame : int, optional
        Reference frame (1-indexed) when alignment_mode="frame".
    """

    threshold: Optional[float] = Field(
        default=DEFAULT_DISTANCE_THRESHOLD,
        description="Distance threshold for contact analysis (Angstroms)",
    )
    pairs: list[DistancePairSettings] = Field(
        default_factory=list, description="Distance pairs to monitor"
    )

    # PBC and alignment settings (new)
    use_pbc: bool = Field(
        default=True,
        description="Use PBC-aware minimum image distances",
    )
    align_trajectory: bool = Field(
        default=True,
        description="Align trajectory before distance calculation (removes drift)",
    )
    alignment_selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection for trajectory alignment",
    )
    alignment_mode: str = Field(
        default="centroid",
        description="Reference mode: centroid, average, or frame",
    )
    alignment_frame: Optional[int] = Field(
        default=None,
        description="Reference frame (1-indexed) when alignment_mode='frame'",
    )

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "distances"

    @field_validator("pairs", mode="after")
    @classmethod
    def validate_pairs(cls, v: list[DistancePairSettings]) -> list[DistancePairSettings]:
        """Ensure at least one pair is defined."""
        if len(v) == 0:
            raise ValueError("At least one distance pair must be defined")
        return v

    @field_validator("alignment_mode", mode="after")
    @classmethod
    def validate_alignment_mode(cls, v: str) -> str:
        """Validate alignment mode is one of the allowed values."""
        valid = {"centroid", "average", "frame"}
        if v not in valid:
            raise ValueError(f"alignment_mode must be one of {valid}, got '{v}'")
        return v

    @model_validator(mode="after")
    def validate_alignment_frame_required(self) -> "DistancesAnalysisSettings":
        """Ensure alignment_frame is provided when alignment_mode is 'frame'."""
        if (
            self.align_trajectory
            and self.alignment_mode == "frame"
            and self.alignment_frame is None
        ):
            raise ValueError("alignment_frame is required when alignment_mode is 'frame'")
        return self

    def to_analysis_yaml_dict(self) -> dict[str, Any]:
        """Convert to analysis.yaml-compatible dictionary."""
        result: dict[str, Any] = {
            "enabled": True,
            "pairs": [p.to_analysis_yaml_dict() for p in self.pairs],
            "use_pbc": self.use_pbc,
            "align_trajectory": self.align_trajectory,
        }
        if self.align_trajectory:
            result["alignment_selection"] = self.alignment_selection
            result["alignment_mode"] = self.alignment_mode
            if self.alignment_frame is not None:
                result["alignment_frame"] = self.alignment_frame
        return result

    def get_pair_selections(self) -> list[tuple[str, str]]:
        """Get list of (selection_a, selection_b) tuples."""
        return [(p.selection_a, p.selection_b) for p in self.pairs]

    def get_pair_labels(self) -> list[str]:
        """Get list of pair labels."""
        return [p.label for p in self.pairs]

    def get_pair_thresholds(self) -> list[float | None]:
        """Get list of thresholds per pair, using global threshold as fallback.

        Returns
        -------
        list[float | None]
            List of thresholds, one per pair. If a pair has no explicit threshold,
            the global threshold is used. If neither is set, None is returned.
        """
        return [p.threshold if p.threshold is not None else self.threshold for p in self.pairs]

    def get_alignment_config(self) -> "AlignmentConfig":
        """Build an AlignmentConfig from these settings.

        Returns
        -------
        AlignmentConfig
            Configuration for trajectory alignment, ready to pass to
            align_trajectory() or DistanceCalculator.

        Notes
        -----
        Import is done inside the method to avoid circular imports.
        """
        from polyzymd.analysis.core.alignment import AlignmentConfig

        return AlignmentConfig(
            enabled=self.align_trajectory,
            reference_mode=self.alignment_mode,  # type: ignore[arg-type]
            reference_frame=self.alignment_frame,
            selection=self.alignment_selection,
        )


@ComparisonSettingsRegistry.register("distances")
class DistancesComparisonSettings(BaseComparisonSettings):
    """Comparison settings for distance analysis.

    Currently empty — all distance comparison behavior uses defaults from
    ``BaseComparisonSettings``.  This class exists as an extension point:
    add fields here when distance-specific comparison parameters are needed
    (e.g., per-pair significance thresholds) without modifying the
    orchestrator or other comparison types.
    """

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "distances"


# ============================================================================
# Catalytic Triad Settings
# ============================================================================


class TriadPairSettings(BaseAnalysisSettings):
    """Configuration for one distance pair in a catalytic triad/active site.

    Attributes
    ----------
    label : str
        Human-readable label for this pair (e.g., "Asp133-His156").
    selection_a : str
        First atom/point selection.
    selection_b : str
        Second atom/point selection.
    """

    label: str = Field(..., description="Human-readable label for this pair")
    selection_a: str = Field(..., description="First atom/point selection")
    selection_b: str = Field(..., description="Second atom/point selection")

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "triad_pair"

    def to_analysis_yaml_dict(self) -> dict[str, Any]:
        """Convert to analysis.yaml-compatible dictionary."""
        return {
            "label": self.label,
            "selection_a": self.selection_a,
            "selection_b": self.selection_b,
        }


@AnalysisSettingsRegistry.register("catalytic_triad")
class CatalyticTriadAnalysisSettings(BaseAnalysisSettings):
    """Catalytic triad/active site analysis settings.

    Attributes
    ----------
    name : str
        Name of the triad/active site (e.g., "LipA_catalytic_triad").
    pairs : list[TriadPairSettings]
        Distance pairs to monitor.
    threshold : float
        Distance threshold for contact/H-bond analysis (Angstroms).
    description : str, optional
        Description of the active site.
    """

    name: str = Field(..., description="Name of the catalytic triad/active site")
    pairs: list[TriadPairSettings] = Field(..., description="Distance pairs to monitor")
    threshold: float = Field(
        default=DEFAULT_DISTANCE_THRESHOLD,
        description="Distance threshold for contact analysis (Angstroms)",
    )
    description: Optional[str] = Field(default=None, description="Description of the active site")

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "catalytic_triad"

    @field_validator("pairs", mode="after")
    @classmethod
    def validate_pairs(cls, v: list[TriadPairSettings]) -> list[TriadPairSettings]:
        """Ensure at least one pair is defined."""
        if len(v) == 0:
            raise ValueError("At least one distance pair must be defined")
        return v

    @property
    def n_pairs(self) -> int:
        """Number of distance pairs."""
        return len(self.pairs)

    def get_pair_selections(self) -> list[tuple[str, str]]:
        """Get list of (selection_a, selection_b) tuples."""
        return [(p.selection_a, p.selection_b) for p in self.pairs]

    def get_pair_labels(self) -> list[str]:
        """Get list of pair labels."""
        return [p.label for p in self.pairs]

    def to_analysis_yaml_dict(self) -> dict[str, Any]:
        """Convert to analysis.yaml-compatible dictionary."""
        result: dict[str, Any] = {
            "enabled": True,
            "name": self.name,
            "threshold": self.threshold,
            "pairs": [p.to_analysis_yaml_dict() for p in self.pairs],
        }
        if self.description:
            result["description"] = self.description
        return result


@ComparisonSettingsRegistry.register("catalytic_triad")
class CatalyticTriadComparisonSettings(BaseComparisonSettings):
    """Comparison settings for catalytic triad analysis.

    Currently empty — all triad comparison behavior uses defaults from
    ``BaseComparisonSettings``.  This class exists as an extension point:
    add fields here when triad-specific comparison parameters are needed
    (e.g., functional distance thresholds) without modifying the
    orchestrator or other comparison types.
    """

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "catalytic_triad"


# ============================================================================
# Polymer-Protein Contacts Settings
# ============================================================================


class BindingPreferenceFieldsMixin(BaseAnalysisSettings):
    """Shared fields for binding preference analysis.

    Both ``ContactsAnalysisSettings`` and ``BindingFreeEnergyAnalysisSettings``
    need identical fields for surface exposure, protein grouping, and polymer
    type selection. This mixin provides them once, keeping defaults in sync.

    Attributes
    ----------
    surface_exposure_threshold : float
        Relative SASA threshold for surface exposure (0.0-1.0).
    enzyme_pdb_for_sasa : str, optional
        Path to enzyme PDB for SASA calculation.
    include_default_aa_groups : bool
        Include default AA class groupings (aromatic, polar, etc.).
    protein_groups : dict[str, list[int]], optional
        Custom protein groups as {name: [resid1, resid2, ...]}.
    protein_partitions : dict[str, list[str]], optional
        Custom partitions for system coverage comparison.
    polymer_type_selections : dict[str, str], optional
        Custom polymer type selections as {name: "MDAnalysis selection"}.
    """

    surface_exposure_threshold: float = Field(
        default=DEFAULT_SURFACE_EXPOSURE_THRESHOLD,
        ge=0.0,
        le=1.0,
        description="Relative SASA threshold for surface exposure (0.2 = 20%)",
    )
    enzyme_pdb_for_sasa: Optional[str] = Field(
        default=None,
        description="Path to enzyme PDB for SASA calculation (relative to comparison.yaml)",
    )
    include_default_aa_groups: bool = Field(
        default=True,
        description="Include default AA class groupings (aromatic, polar, nonpolar, charged)",
    )
    protein_groups: Optional[dict[str, list[int]]] = Field(
        default=None,
        description="Custom protein groups as {name: [resid1, resid2, ...]}",
    )
    protein_partitions: Optional[dict[str, list[str]]] = Field(
        default=None,
        description=(
            "Custom partitions for system coverage comparison. "
            "Each partition defines a mutually exclusive set of protein groups "
            "that will generate one comparison plot. Format: {partition_name: [group1, group2, ...]}. "
            "Groups must be defined in protein_groups. If groups don't cover all protein residues, "
            "'rest_of_protein' is auto-added. Overlapping groups within a partition cause validation error."
        ),
    )
    polymer_type_selections: Optional[dict[str, str]] = Field(
        default=None,
        description="Custom polymer type selections as {name: 'MDAnalysis selection'}",
    )

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier (override in subclass)."""
        raise NotImplementedError


@AnalysisSettingsRegistry.register("contacts")
class ContactsAnalysisSettings(BindingPreferenceFieldsMixin):
    """Polymer-protein contact analysis settings.

    Inherits binding preference fields (surface_exposure_threshold,
    enzyme_pdb_for_sasa, include_default_aa_groups, protein_groups,
    protein_partitions, polymer_type_selections) from
    ``BindingPreferenceFieldsMixin``.

    Attributes
    ----------
    polymer_selection : str
        MDAnalysis selection for polymer atoms.
    protein_selection : str
        MDAnalysis selection for protein atoms.
    cutoff : float
        Distance cutoff for contacts in Angstroms.
    polymer_types : list[str], optional
        Filter contacts by polymer residue names.
    grouping : str
        How to group protein residues: aa_class, secondary_structure, or none.
    compute_residence_times : bool
        If True, compute residence time statistics.
    compute_binding_preference : bool
        If True, compute binding preference enrichment analysis.
    enrichment_normalization : str
        **DEPRECATED** (kept for backward compatibility).
        Enrichment is now always normalized by protein surface availability.
        This field is ignored.
    """

    polymer_selection: str = Field(
        default="chainID C", description="MDAnalysis selection for polymer atoms"
    )
    protein_selection: str = Field(
        default="protein", description="MDAnalysis selection for protein atoms"
    )
    cutoff: float = Field(
        default=DEFAULT_CONTACT_CUTOFF, description="Contact distance cutoff in Angstroms"
    )
    polymer_types: Optional[list[str]] = Field(
        default=None, description="Filter by polymer residue names"
    )
    grouping: str = Field(
        default="aa_class", description="Group by: aa_class, secondary_structure, or none"
    )
    compute_residence_times: bool = Field(
        default=True, description="Compute residence time statistics"
    )

    # Binding preference settings
    compute_binding_preference: bool = Field(
        default=False, description="Compute binding preference enrichment analysis"
    )
    enrichment_normalization: str = Field(
        default="residue",
        description="DEPRECATED: Enrichment is now always normalized by protein surface availability. This field is ignored.",
    )

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "contacts"

    @field_validator("grouping", mode="after")
    @classmethod
    def validate_grouping(cls, v: str) -> str:
        """Validate grouping mode."""
        valid = {"aa_class", "secondary_structure", "none"}
        if v not in valid:
            raise ValueError(f"grouping must be one of {valid}, got '{v}'")
        return v

    @field_validator("enrichment_normalization", mode="after")
    @classmethod
    def validate_enrichment_normalization(cls, v: str) -> str:
        """Validate enrichment normalization method."""
        valid = {"residue", "atoms"}
        if v not in valid:
            raise ValueError(f"enrichment_normalization must be one of {valid}, got '{v}'")
        return v

    @model_validator(mode="after")
    def validate_protein_partitions(self) -> "ContactsAnalysisSettings":
        """Validate protein_partitions references and mutual exclusivity.

        Validates:
        1. All groups referenced in partitions exist in protein_groups
        2. Groups within each partition don't overlap (mutually exclusive)
        """
        if not self.protein_partitions:
            return self

        # protein_groups must exist if partitions are defined
        if not self.protein_groups:
            raise ValueError(
                "protein_partitions requires protein_groups to be defined. "
                "Define the groups first, then reference them in partitions."
            )

        protein_groups = self.protein_groups

        for partition_name, group_names in self.protein_partitions.items():
            if not group_names:
                raise ValueError(
                    f"Partition '{partition_name}' is empty. "
                    "Each partition must contain at least one group."
                )

            # Check all referenced groups exist
            for group_name in group_names:
                if group_name not in protein_groups:
                    available = ", ".join(sorted(protein_groups.keys()))
                    raise ValueError(
                        f"Partition '{partition_name}' references undefined group '{group_name}'. "
                        f"Available groups: {available}"
                    )

            # Check for overlapping groups within this partition
            seen_resids: dict[int, str] = {}  # resid -> first group that contains it
            for group_name in group_names:
                group_resids = protein_groups[group_name]
                for resid in group_resids:
                    if resid in seen_resids:
                        raise ValueError(
                            f"Partition '{partition_name}' has overlapping groups: "
                            f"residue {resid} is in both '{seen_resids[resid]}' and '{group_name}'. "
                            "Groups within a partition must be mutually exclusive."
                        )
                    seen_resids[resid] = group_name

        return self

    def to_analysis_yaml_dict(self) -> dict[str, Any]:
        """Convert to analysis.yaml-compatible dictionary."""
        result: dict[str, Any] = {
            "enabled": True,
            "polymer_selection": self.polymer_selection,
            "protein_selection": self.protein_selection,
            "cutoff": self.cutoff,
            "grouping": self.grouping,
            "compute_residence_times": self.compute_residence_times,
        }
        if self.polymer_types:
            result["polymer_types"] = self.polymer_types

        # Binding preference settings (only include if enabled)
        if self.compute_binding_preference:
            result["compute_binding_preference"] = True
            result["surface_exposure_threshold"] = self.surface_exposure_threshold
            result["include_default_aa_groups"] = self.include_default_aa_groups
            # Note: enrichment_normalization is deprecated and no longer included
            if self.enzyme_pdb_for_sasa:
                result["enzyme_pdb_for_sasa"] = self.enzyme_pdb_for_sasa
            if self.protein_groups:
                result["protein_groups"] = self.protein_groups
            if self.protein_partitions:
                result["protein_partitions"] = self.protein_partitions
            if self.polymer_type_selections:
                result["polymer_type_selections"] = self.polymer_type_selections

        return result


@ComparisonSettingsRegistry.register("contacts")
class ContactsComparisonSettings(BaseComparisonSettings):
    """Comparison settings for polymer-protein contacts analysis.

    Attributes
    ----------
    fdr_alpha : float
        False discovery rate alpha for Benjamini-Hochberg correction.
    min_effect_size : float
        Minimum Cohen's d effect size to highlight in reports.
    top_residues : int
        Number of top residues (by effect size) to display in console.
    """

    fdr_alpha: float = Field(
        default=0.05, description="FDR alpha for Benjamini-Hochberg correction"
    )
    min_effect_size: float = Field(
        default=0.5,
        description="Minimum Cohen's d to highlight (0.2=small, 0.5=medium, 0.8=large)",
    )
    top_residues: int = Field(
        default=10, description="Number of top residues to display in console"
    )

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "contacts"

    @field_validator("fdr_alpha", mode="after")
    @classmethod
    def validate_fdr_alpha(cls, v: float) -> float:
        """Validate FDR alpha is in valid range."""
        if not 0 < v < 1:
            raise ValueError(f"fdr_alpha must be between 0 and 1, got {v}")
        return v


# ============================================================================
# Utility Functions
# ============================================================================


# ============================================================================
# Exposure Dynamics Settings
# ============================================================================


@AnalysisSettingsRegistry.register("exposure")
class ExposureAnalysisSettings(BaseAnalysisSettings):
    """Exposure dynamics analysis settings (dynamic SASA-based chaperone analysis).

    Attributes
    ----------
    protein_selection : str
        MDAnalysis selection for protein atoms (chain A by default).
    polymer_selection : str
        MDAnalysis selection for polymer atoms (chain C by default).
    exposure_threshold : float
        Relative SASA threshold for classifying a residue as exposed.
    transient_lower : float
        Lower bound of exposure fraction for "transient" classification.
    transient_upper : float
        Upper bound of exposure fraction for "transient" classification.
    min_event_length : int
        Minimum exposed-window length (frames) to count as an event.
    probe_radius_nm : float
        Probe radius for MDTraj shrake_rupley, in nm.
    n_sphere_points : int
        Number of sphere points for shrake_rupley.
    protein_chain : str
        Chain letter for protein (default "A").
    polymer_resnames : list[str], optional
        Subset of polymer monomer resnames to include. If None, all detected.
    """

    protein_selection: str = Field(
        default="protein", description="MDAnalysis selection for protein"
    )
    polymer_selection: str = Field(
        default="chainID C", description="MDAnalysis selection for polymer"
    )
    exposure_threshold: float = Field(
        default=DEFAULT_SURFACE_EXPOSURE_THRESHOLD,
        ge=0.0,
        le=1.0,
        description="Relative SASA threshold for exposed classification",
    )
    transient_lower: float = Field(
        default=0.2,
        ge=0.0,
        le=1.0,
        description="Lower exposure fraction bound for 'transient' residues",
    )
    transient_upper: float = Field(
        default=0.8,
        ge=0.0,
        le=1.0,
        description="Upper exposure fraction bound for 'transient' residues",
    )
    min_event_length: int = Field(
        default=1,
        ge=1,
        description="Minimum exposed-window length (frames) to count as event",
    )
    probe_radius_nm: float = Field(default=0.14, description="Probe radius for SASA in nm")
    n_sphere_points: int = Field(
        default=960, description="Number of sphere points for shrake_rupley"
    )
    protein_chain: str = Field(default="A", description="Chain letter for protein")
    polymer_resnames: Optional[list[str]] = Field(
        default=None,
        description="Subset of polymer resnames to analyze. If None, all detected.",
    )

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "exposure"

    def to_analysis_yaml_dict(self) -> dict[str, Any]:
        """Convert to analysis.yaml-compatible dictionary."""
        result: dict[str, Any] = {
            "enabled": True,
            "exposure_threshold": self.exposure_threshold,
            "transient_lower": self.transient_lower,
            "transient_upper": self.transient_upper,
            "min_event_length": self.min_event_length,
            "protein_chain": self.protein_chain,
        }
        if self.polymer_resnames:
            result["polymer_resnames"] = self.polymer_resnames
        return result


@ComparisonSettingsRegistry.register("exposure")
class ExposureComparisonSettings(BaseComparisonSettings):
    """Comparison settings for exposure dynamics analysis.

    Currently empty — all exposure comparison behavior uses defaults from
    ``BaseComparisonSettings``.  This class exists as an extension point:
    add fields here when exposure-specific comparison parameters are needed
    (e.g., transient classification thresholds) without modifying the
    orchestrator or other comparison types.
    """

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "exposure"


# ============================================================================
# Binding Free Energy Settings
# ============================================================================


@AnalysisSettingsRegistry.register("binding_free_energy")
class BindingFreeEnergyAnalysisSettings(BindingPreferenceFieldsMixin):
    """Settings for binding free energy analysis via Boltzmann inversion.

    Computes the selectivity free energy difference:

        ΔΔG = -k_B·T · ln(contact_share / expected_share)

    where:
    - contact_share  = fraction of polymer contacts directed at an AA group
    - expected_share = fraction of exposed surface belonging to that AA group
    - T              = simulation temperature (from SimulationConfig)

    This is a post-processing analysis that consumes binding preference results
    from the contacts analysis layer (no new per-frame computation is needed).

    Inherits binding preference fields (surface_exposure_threshold,
    enzyme_pdb_for_sasa, include_default_aa_groups, protein_groups,
    protein_partitions, polymer_type_selections) from
    ``BindingPreferenceFieldsMixin``.

    Attributes
    ----------
    units : str
        Energy units for output. One of "kT" (dimensionless, in units of
        k_bT — the thermal energy), "kcal/mol", or "kJ/mol".
    compute_binding_preference : bool
        Compute binding preference from contacts data when cached results
        are not found.
    """

    units: str = Field(
        default="kT",
        description="Energy units: 'kT' (default, dimensionless), 'kcal/mol', or 'kJ/mol'",
    )
    compute_binding_preference: bool = Field(
        default=True,
        description=(
            "Compute binding preference from contacts data when cached results "
            "are not found. Set to False to only load pre-existing results."
        ),
    )

    @field_validator("units")
    @classmethod
    def validate_units(cls, v: str) -> str:
        """Validate energy units."""
        allowed = {"kT", "kcal/mol", "kJ/mol"}
        if v not in allowed:
            raise ValueError(f"units must be one of {sorted(allowed)}, got '{v}'")
        return v

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "binding_free_energy"

    def k_b(self) -> float:
        """Return k_B in the selected energy units.

        Returns
        -------
        float
            Boltzmann constant in kcal/(mol·K) or kJ/(mol·K).
            When units='kT', returns 0.0 — callers should use kT=1.0 directly
            instead of k_b() * T.
        """
        if self.units == "kT":
            return 0.0  # Not used; comparator sets kT=1.0 directly
        if self.units == "kJ/mol":
            return 0.0083144626  # kJ/(mol·K)
        return 0.0019872041  # kcal/(mol·K)  [default]

    def to_analysis_yaml_dict(self) -> dict:
        """Convert to analysis.yaml-compatible dictionary.

        Returns
        -------
        dict
            Dictionary suitable for writing to analysis.yaml.
        """
        result: dict = {
            "enabled": True,
            "units": self.units,
            "compute_binding_preference": self.compute_binding_preference,
            "surface_exposure_threshold": self.surface_exposure_threshold,
        }
        if self.enzyme_pdb_for_sasa is not None:
            result["enzyme_pdb_for_sasa"] = self.enzyme_pdb_for_sasa
        if self.protein_groups is not None:
            result["protein_groups"] = self.protein_groups
        if self.protein_partitions is not None:
            result["protein_partitions"] = self.protein_partitions
        if self.polymer_type_selections is not None:
            result["polymer_type_selections"] = self.polymer_type_selections
        return result


@ComparisonSettingsRegistry.register("binding_free_energy")
class BindingFreeEnergyComparisonSettings(BaseComparisonSettings):
    """Comparison settings for binding free energy analysis.

    Attributes
    ----------
    fdr_alpha : float
        False discovery rate alpha for Benjamini-Hochberg correction
        of p-values across (polymer_type, AA_group) pairs.
    """

    fdr_alpha: float = Field(
        default=0.05,
        description="FDR alpha for Benjamini-Hochberg correction",
    )

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "binding_free_energy"


def get_all_analysis_types() -> list[str]:
    """Get all registered analysis types.

    Returns
    -------
    list[str]
        Sorted list of registered analysis type names.
    """
    return AnalysisSettingsRegistry.list_available()


def get_all_comparison_types() -> list[str]:
    """Get all registered comparison settings types.

    Returns
    -------
    list[str]
        Sorted list of registered comparison type names.
    """
    return ComparisonSettingsRegistry.list_available()
