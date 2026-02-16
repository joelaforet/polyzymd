"""Configuration schema for comparison projects.

This module defines the YAML schema for comparison.yaml files that
specify which simulation conditions to compare.

Also includes catalytic triad/active site configuration for distance
and H-bond analysis.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator

# ============================================================================
# RMSF Comparison Configuration
# ============================================================================


class RMSFComparisonConfig(BaseModel):
    """Configuration for RMSF comparison analysis.

    Required to run `polyzymd compare rmsf`. If this section is not
    present in comparison.yaml, RMSF comparison will fail with an
    informative error message.

    Attributes
    ----------
    selection : str
        MDAnalysis selection string for RMSF calculation.
        Default: "protein and name CA"
    reference_mode : str
        Reference structure mode: centroid, average, or frame.
        Default: "centroid"
    reference_frame : int, optional
        Frame number if reference_mode is 'frame' (1-indexed, PyMOL convention)

    Examples
    --------
    In comparison.yaml:

    ```yaml
    rmsf:
      selection: "protein and name CA"
      reference_mode: "centroid"
    ```

    For aligning to a specific frame:

    ```yaml
    rmsf:
      selection: "protein and name CA"
      reference_mode: "frame"
      reference_frame: 500
    ```
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

    @field_validator("reference_mode", mode="after")
    @classmethod
    def validate_reference_mode(cls, v: str) -> str:
        """Validate reference mode is one of the allowed values."""
        valid = {"centroid", "average", "frame"}
        if v not in valid:
            raise ValueError(f"reference_mode must be one of {valid}, got '{v}'")
        return v

    @model_validator(mode="after")
    def validate_reference_frame_required(self) -> "RMSFComparisonConfig":
        """Ensure reference_frame is provided when reference_mode is 'frame'."""
        if self.reference_mode == "frame" and self.reference_frame is None:
            raise ValueError("reference_frame is required when reference_mode is 'frame'")
        return self


# ============================================================================
# Catalytic Triad / Active Site Configuration
# ============================================================================


class TriadPairConfig(BaseModel):
    """Configuration for one distance pair in a catalytic triad/active site.

    Each pair defines two atom selections between which distances will be
    measured. Selections support special syntax:
    - Standard MDAnalysis: "resid 77 and name OG"
    - Midpoint: "midpoint(resid 133 and name OD1 OD2)"
    - Center of mass: "com(resid 50-75)"

    Attributes
    ----------
    label : str
        Human-readable label for this pair (e.g., "Asp133-His156")
    selection_a : str
        First atom/point selection
    selection_b : str
        Second atom/point selection

    Examples
    --------
    >>> pair = TriadPairConfig(
    ...     label="Asp133-His156",
    ...     selection_a="midpoint(resid 133 and name OD1 OD2)",
    ...     selection_b="resid 156 and name ND1",
    ... )
    """

    label: str = Field(..., description="Human-readable label for this pair")
    selection_a: str = Field(..., description="First atom/point selection")
    selection_b: str = Field(..., description="Second atom/point selection")


class CatalyticTriadConfig(BaseModel):
    """Configuration for catalytic triad/active site distance analysis.

    Defines multiple distance pairs that together comprise an active site
    or catalytic machinery. Computes:
    - Per-pair distance distributions and statistics
    - Simultaneous contact fraction (all pairs below threshold at same time)

    Attributes
    ----------
    name : str
        Name of the triad/active site (e.g., "LipA_catalytic_triad")
    pairs : list[TriadPairConfig]
        Distance pairs to monitor
    threshold : float
        Distance threshold for contact/H-bond analysis (Angstroms)
    description : str, optional
        Description of the active site

    Examples
    --------
    In comparison.yaml:

    ```yaml
    catalytic_triad:
      name: "LipA_catalytic_triad"
      description: "Ser-His-Asp catalytic triad of Lipase A"
      threshold: 3.5
      pairs:
        - label: "Asp133-His156"
          selection_a: "midpoint(resid 133 and name OD1 OD2)"
          selection_b: "resid 156 and name ND1"
        - label: "His156-Ser77"
          selection_a: "resid 156 and name NE2"
          selection_b: "resid 77 and name OG"
    ```
    """

    name: str = Field(..., description="Name of the catalytic triad/active site")
    pairs: list[TriadPairConfig] = Field(..., description="Distance pairs to monitor")
    threshold: float = Field(
        default=3.5, description="Distance threshold for contact analysis (Angstroms)"
    )
    description: Optional[str] = Field(default=None, description="Description of the active site")

    @field_validator("pairs", mode="after")
    @classmethod
    def validate_pairs(cls, v: list[TriadPairConfig]) -> list[TriadPairConfig]:
        """Ensure at least one pair is defined."""
        if len(v) == 0:
            raise ValueError("At least one distance pair must be defined")
        return v

    @property
    def n_pairs(self) -> int:
        """Number of distance pairs."""
        return len(self.pairs)

    def get_pair_selections(self) -> list[tuple[str, str]]:
        """Get list of (selection_a, selection_b) tuples for DistanceCalculator."""
        return [(p.selection_a, p.selection_b) for p in self.pairs]

    def get_pair_labels(self) -> list[str]:
        """Get list of pair labels."""
        return [p.label for p in self.pairs]


# ============================================================================
# Polymer-Protein Contacts Configuration
# ============================================================================


class ContactsComparisonConfig(BaseModel):
    """Configuration for polymer-protein contacts comparison analysis.

    Defines parameters for comparing contact statistics across conditions,
    including polymer selection, contact criteria, and statistical thresholds.

    Attributes
    ----------
    name : str
        Name of the contacts analysis (e.g., "polymer_protein_contacts")
    polymer_selection : str
        MDAnalysis selection string for polymer atoms.
        Default: "resname SBM EGM" (3-character PDB standard for SBMA/EGMA)
    protein_selection : str
        MDAnalysis selection string for protein atoms.
        Default: "protein"
    cutoff : float
        Distance cutoff for contacts in Angstroms.
        Default: 4.5
    contact_criteria : str
        Contact criteria: "distance", "heavy_atom", or "any_atom".
        Default: "heavy_atom"
    fdr_alpha : float
        False discovery rate alpha for Benjamini-Hochberg correction.
        Default: 0.05
    min_effect_size : float
        Minimum Cohen's d effect size to highlight in reports.
        Default: 0.5 (medium effect)
    top_residues : int
        Number of top residues (by effect size) to display in console.
        Default: 10
    description : str, optional
        Description of the contacts analysis

    Examples
    --------
    In comparison.yaml:

    ```yaml
    contacts:
      name: "polymer_contacts"
      description: "Polymer-protein contact analysis"
      polymer_selection: "resname SBM EGM"
      protein_selection: "protein"
      cutoff: 4.5
      contact_criteria: "heavy_atom"
      fdr_alpha: 0.05
      min_effect_size: 0.5
      top_residues: 10
    ```
    """

    name: str = Field(
        default="polymer_protein_contacts", description="Name of the contacts analysis"
    )
    polymer_selection: str = Field(
        default="resname SBM EGM", description="MDAnalysis selection for polymer atoms"
    )
    protein_selection: str = Field(
        default="protein", description="MDAnalysis selection for protein atoms"
    )
    cutoff: float = Field(default=4.5, description="Contact distance cutoff in Angstroms")
    contact_criteria: str = Field(
        default="heavy_atom", description="Contact criteria: distance, heavy_atom, or any_atom"
    )
    fdr_alpha: float = Field(
        default=0.05, description="FDR alpha for Benjamini-Hochberg correction"
    )
    min_effect_size: float = Field(
        default=0.5, description="Minimum Cohen's d to highlight (0.2=small, 0.5=medium, 0.8=large)"
    )
    top_residues: int = Field(
        default=10, description="Number of top residues to display in console"
    )
    description: Optional[str] = Field(
        default=None, description="Description of the contacts analysis"
    )

    @field_validator("contact_criteria", mode="after")
    @classmethod
    def validate_criteria(cls, v: str) -> str:
        """Validate contact criteria."""
        valid = {"distance", "heavy_atom", "any_atom"}
        if v not in valid:
            raise ValueError(f"contact_criteria must be one of {valid}, got '{v}'")
        return v

    @field_validator("fdr_alpha", mode="after")
    @classmethod
    def validate_fdr_alpha(cls, v: float) -> float:
        """Validate FDR alpha is in valid range."""
        if not 0 < v < 1:
            raise ValueError(f"fdr_alpha must be between 0 and 1, got {v}")
        return v


# ============================================================================
# Comparison Configuration
# ============================================================================


class ConditionConfig(BaseModel):
    """Configuration for one condition in a comparison.

    Attributes
    ----------
    label : str
        Display name for this condition (e.g., "No Polymer", "100% SBMA")
    config : Path
        Path to the simulation's config.yaml file
    replicates : list[int]
        List of replicate numbers to include in the analysis
    """

    label: str
    config: Path
    replicates: list[int]

    @field_validator("config", mode="before")
    @classmethod
    def resolve_path(cls, v: str | Path) -> Path:
        """Convert string paths to Path objects."""
        return Path(v)

    @field_validator("replicates", mode="before")
    @classmethod
    def ensure_list(cls, v: list[int] | int) -> list[int]:
        """Ensure replicates is a list."""
        if isinstance(v, int):
            return [v]
        return list(v)


class AnalysisDefaults(BaseModel):
    """Default parameters for comparison analyses.

    These can be overridden on the command line.

    Attributes
    ----------
    equilibration_time : str
        Time to skip for equilibration (e.g., "10ns", "5000ps").
        Shared across all analyses (RMSF, contacts, triad).
    """

    equilibration_time: str = "10ns"


class ComparisonConfig(BaseModel):
    """Schema for comparison.yaml configuration files.

    A comparison config defines multiple simulation conditions to compare,
    along with default analysis parameters and optional analysis-specific
    configuration sections.

    Attributes
    ----------
    name : str
        Name of the comparison project
    description : str, optional
        Description of what is being compared
    control : str, optional
        Label of the control condition for relative comparisons
    conditions : list[ConditionConfig]
        List of conditions to compare
    defaults : AnalysisDefaults
        Default analysis parameters (equilibration_time)
    rmsf : RMSFComparisonConfig, optional
        Configuration for RMSF comparison analysis.
        Required to run `polyzymd compare rmsf`.
    catalytic_triad : CatalyticTriadConfig, optional
        Configuration for catalytic triad/active site analysis.
        Required to run `polyzymd compare triad`.
    contacts : ContactsComparisonConfig, optional
        Configuration for polymer-protein contacts analysis.
        Required to run `polyzymd compare contacts`.

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> print(config.name)
    "Polymer Stabilization Study"
    >>> for cond in config.conditions:
    ...     print(f"{cond.label}: {cond.config}")
    >>> if config.rmsf:
    ...     print(f"RMSF selection: {config.rmsf.selection}")
    >>> if config.catalytic_triad:
    ...     print(f"Triad: {config.catalytic_triad.name}")
    >>> if config.contacts:
    ...     print(f"Contacts: {config.contacts.name}")
    """

    name: str
    description: Optional[str] = None
    control: Optional[str] = None
    conditions: list[ConditionConfig]
    defaults: AnalysisDefaults = AnalysisDefaults()
    rmsf: Optional[RMSFComparisonConfig] = None
    catalytic_triad: Optional[CatalyticTriadConfig] = None
    contacts: Optional[ContactsComparisonConfig] = None

    @classmethod
    def from_yaml(cls, path: Path | str) -> "ComparisonConfig":
        """Load comparison config from YAML file.

        Parameters
        ----------
        path : Path or str
            Path to comparison.yaml file

        Returns
        -------
        ComparisonConfig
            Loaded and validated configuration

        Raises
        ------
        FileNotFoundError
            If the config file doesn't exist
        ValidationError
            If the config is invalid
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Comparison config not found: {path}")

        with open(path) as f:
            data = yaml.safe_load(f)

        # Resolve relative paths relative to the config file location
        config_dir = path.parent.resolve()
        if "conditions" in data:
            for cond in data["conditions"]:
                if "config" in cond:
                    cond_path = Path(cond["config"])
                    if not cond_path.is_absolute():
                        cond["config"] = str(config_dir / cond_path)

        return cls(**data)

    def to_yaml(self, path: Path | str) -> None:
        """Save comparison config to YAML file.

        Parameters
        ----------
        path : Path or str
            Output path for comparison.yaml
        """
        path = Path(path)

        # Convert to dict, handling Path objects
        data = self.model_dump()
        for cond in data["conditions"]:
            cond["config"] = str(cond["config"])

        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    def get_condition(self, label: str) -> ConditionConfig:
        """Get a condition by its label.

        Parameters
        ----------
        label : str
            The condition label to find

        Returns
        -------
        ConditionConfig
            The matching condition

        Raises
        ------
        KeyError
            If no condition with that label exists
        """
        for cond in self.conditions:
            if cond.label == label:
                return cond
        raise KeyError(f"Condition '{label}' not found in: {[c.label for c in self.conditions]}")

    def validate_config(self) -> list[str]:
        """Validate the comparison configuration.

        Returns
        -------
        list[str]
            List of error messages (empty if valid)
        """
        errors = []

        # Check minimum conditions
        if len(self.conditions) < 2:
            errors.append("Need at least 2 conditions to compare")

        # Check for duplicate labels
        labels = [c.label for c in self.conditions]
        if len(labels) != len(set(labels)):
            errors.append("Duplicate condition labels found")

        # Check control label exists
        if self.control and self.control not in labels:
            errors.append(f"Control '{self.control}' not in conditions: {labels}")

        # Check config files exist
        for cond in self.conditions:
            if not cond.config.exists():
                errors.append(f"Config not found for '{cond.label}': {cond.config}")

        return errors


def generate_comparison_template(name: str, eq_time: str = "10ns") -> str:
    """Generate a template comparison.yaml file.

    Parameters
    ----------
    name : str
        Project name
    eq_time : str
        Default equilibration time

    Returns
    -------
    str
        YAML template content
    """
    return f"""\
# ============================================================================
# PolyzyMD Comparison Configuration
# ============================================================================
# Compare analyses across simulation conditions (e.g., polymer, temperature).
# Docs: https://polyzymd.readthedocs.io/en/latest/
# ============================================================================

name: "{name}"
description: "Comparison of simulation conditions"

# Control condition for relative comparisons.
# Must match one of the 'label' values in 'conditions' below, or null if none.
# Example: control: "noPoly" (if you have a condition labeled "noPoly")
control: null

# ============================================================================
# Conditions
# ============================================================================
# Each condition points to a simulation's config.yaml file.
conditions:
  - label: "Condition A"
    config: "../path/to/condition_a/config.yaml"
    replicates: [1, 2, 3]

  - label: "Condition B"
    config: "../path/to/condition_b/config.yaml"
    replicates: [1, 2, 3]

# ============================================================================
# Defaults
# ============================================================================
# Shared parameters across all analyses.
# Note: If equilibration_time is 0ns, you'll get a warning - set appropriately!
defaults:
  equilibration_time: "{eq_time}"

# ============================================================================
# RMSF (for polyzymd compare rmsf)
# ============================================================================
# Compare protein flexibility across conditions. Uncomment to enable.
#
# Run: polyzymd compare rmsf -c comparison.yaml
#
# rmsf:
#   selection: "protein and name CA"    # MDAnalysis selection string
#   reference_mode: "centroid"          # centroid, average, or frame
#   # reference_frame: 500              # Required if reference_mode is "frame"

# ============================================================================
# Catalytic Triad (for polyzymd compare triad)
# ============================================================================
# Define your enzyme's catalytic machinery ONCE here - it's shared across
# all conditions since the enzyme is the same.
#
# Run: polyzymd compare triad -c comparison.yaml
#
# catalytic_triad:
#   name: "enzyme_catalytic_triad"
#   threshold: 3.5  # Angstroms (H-bond cutoff)
#   pairs:
#     - label: "Asp-His"
#       selection_a: "midpoint(resid 133 and name OD1 OD2)"
#       selection_b: "resid 156 and name ND1"
#     - label: "His-Ser"
#       selection_a: "resid 156 and name NE2"
#       selection_b: "resid 77 and name OG"

# ============================================================================
# Distance Analysis (for polyzymd compare distances)
# ============================================================================
# Compare specific inter-atomic distances across conditions.
# Useful for monitoring substrate positioning, active site geometry, etc.
#
# Run: polyzymd compare distances -c comparison.yaml
#
# distances:
#   threshold: 3.5                          # Optional: contact threshold (Angstroms)
#   pairs:
#     - label: "Ser77-Substrate"
#       selection_a: "resid 77 and name OG"
#       selection_b: "resname RBY and name C1"
#     - label: "His156-Substrate"
#       selection_a: "resid 156 and name NE2"
#       selection_b: "resname RBY and name C1"

# ============================================================================
# Polymer-Protein Contacts (for polyzymd compare contacts)
# ============================================================================
# Configure polymer-protein contact analysis with per-residue statistics.
#
# Run: polyzymd compare contacts -c comparison.yaml
#
# contacts:
#   name: "polymer_protein_contacts"
#   description: "Polymer-protein contact analysis"
#   polymer_selection: "resname SBM EGM"  # 3-char PDB names for SBMA/EGMA
#   protein_selection: "protein"
#   cutoff: 4.5                            # Contact distance (Angstroms)
#   contact_criteria: "heavy_atom"         # distance, heavy_atom, or any_atom
#   fdr_alpha: 0.05                         # FDR for per-residue tests
#   min_effect_size: 0.5                    # Cohen's d threshold (0.5 = medium)
#   top_residues: 10                        # Top residues to show in console
"""
