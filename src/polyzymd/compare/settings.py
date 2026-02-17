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

from typing import Any, Optional

from pydantic import Field, field_validator, model_validator

from polyzymd.analysis.core.registry import (
    AnalysisSettingsRegistry,
    BaseAnalysisSettings,
    BaseComparisonSettings,
    ComparisonSettingsRegistry,
)

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
        Reference structure mode: centroid, average, or frame.
    reference_frame : int, optional
        Frame number if reference_mode is 'frame' (1-indexed).
    """

    selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection string for RMSF calculation",
    )
    reference_mode: str = Field(
        default="centroid",
        description="Reference structure mode: centroid, average, or frame",
    )
    reference_frame: Optional[int] = Field(
        default=None,
        description="Frame number if reference_mode is 'frame' (1-indexed)",
    )

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "rmsf"

    @field_validator("reference_mode", mode="after")
    @classmethod
    def validate_reference_mode(cls, v: str) -> str:
        """Validate reference mode is one of the allowed values."""
        valid = {"centroid", "average", "frame"}
        if v not in valid:
            raise ValueError(f"reference_mode must be one of {valid}, got '{v}'")
        return v

    @model_validator(mode="after")
    def validate_reference_frame_required(self) -> "RMSFAnalysisSettings":
        """Ensure reference_frame is provided when reference_mode is 'frame'."""
        if self.reference_mode == "frame" and self.reference_frame is None:
            raise ValueError("reference_frame is required when reference_mode is 'frame'")
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
        return result


@ComparisonSettingsRegistry.register("rmsf")
class RMSFComparisonSettings(BaseComparisonSettings):
    """Comparison settings for RMSF analysis.

    Currently no comparison-specific parameters. This class exists to
    enforce the pattern: analysis_settings.rmsf + comparison_settings.rmsf.
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
    """

    threshold: Optional[float] = Field(
        default=3.5, description="Distance threshold for contact analysis (Angstroms)"
    )
    pairs: list[DistancePairSettings] = Field(
        default_factory=list, description="Distance pairs to monitor"
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

    def to_analysis_yaml_dict(self) -> dict[str, Any]:
        """Convert to analysis.yaml-compatible dictionary."""
        return {
            "enabled": True,
            "pairs": [p.to_analysis_yaml_dict() for p in self.pairs],
        }

    def get_pair_selections(self) -> list[tuple[str, str]]:
        """Get list of (selection_a, selection_b) tuples."""
        return [(p.selection_a, p.selection_b) for p in self.pairs]

    def get_pair_labels(self) -> list[str]:
        """Get list of pair labels."""
        return [p.label for p in self.pairs]


@ComparisonSettingsRegistry.register("distances")
class DistancesComparisonSettings(BaseComparisonSettings):
    """Comparison settings for distance analysis.

    Currently no comparison-specific parameters beyond analysis settings.
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
        default=3.5, description="Distance threshold for contact analysis (Angstroms)"
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

    Currently no comparison-specific parameters beyond analysis settings.
    """

    @classmethod
    def analysis_type(cls) -> str:
        """Return the analysis type identifier."""
        return "catalytic_triad"


# ============================================================================
# Polymer-Protein Contacts Settings
# ============================================================================


@AnalysisSettingsRegistry.register("contacts")
class ContactsAnalysisSettings(BaseAnalysisSettings):
    """Polymer-protein contact analysis settings.

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
    """

    polymer_selection: str = Field(
        default="chainID C", description="MDAnalysis selection for polymer atoms"
    )
    protein_selection: str = Field(
        default="protein", description="MDAnalysis selection for protein atoms"
    )
    cutoff: float = Field(default=4.5, description="Contact distance cutoff in Angstroms")
    polymer_types: Optional[list[str]] = Field(
        default=None, description="Filter by polymer residue names"
    )
    grouping: str = Field(
        default="aa_class", description="Group by: aa_class, secondary_structure, or none"
    )
    compute_residence_times: bool = Field(
        default=True, description="Compute residence time statistics"
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
