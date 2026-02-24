"""Configuration schema for analysis.yaml files.

This module defines the YAML schema for analysis.yaml files that
configure which analyses to run for a simulation.

The analysis.yaml must live alongside config.yaml (same directory)
to maintain the config.yaml as the single source of truth.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import yaml
from pydantic import BaseModel, Field, field_validator

from polyzymd.analysis.core.constants import (
    DEFAULT_CONTACT_CUTOFF,
    DEFAULT_DISTANCE_THRESHOLD,
    DEFAULT_SURFACE_EXPOSURE_THRESHOLD,
)


# =============================================================================
# Individual Analysis Configurations
# =============================================================================


class RMSFConfig(BaseModel):
    """Configuration for RMSF analysis.

    Attributes
    ----------
    enabled : bool
        Whether to run RMSF analysis
    selection : str
        MDAnalysis selection string for atoms to analyze
    reference_mode : str
        Reference structure mode: 'centroid', 'average', or 'frame'
    reference_frame : int, optional
        Frame number if reference_mode is 'frame'
    """

    enabled: bool = False
    selection: str = "protein and name CA"
    reference_mode: str = "centroid"
    reference_frame: Optional[int] = None


class DistancePairConfig(BaseModel):
    """Configuration for a single distance pair.

    Attributes
    ----------
    label : str
        Human-readable label for this pair
    selection_a : str
        First atom/point selection
    selection_b : str
        Second atom/point selection
    """

    label: str
    selection_a: str
    selection_b: str


class DistancesConfig(BaseModel):
    """Configuration for distance analysis.

    Attributes
    ----------
    enabled : bool
        Whether to run distance analysis
    pairs : list[DistancePairConfig]
        List of atom pairs to measure distances between
    """

    enabled: bool = False
    pairs: list[DistancePairConfig] = Field(default_factory=list)


class TriadPairConfig(BaseModel):
    """Configuration for a catalytic triad distance pair.

    Attributes
    ----------
    label : str
        Human-readable label for this pair (e.g., "Asp133-His156")
    selection_a : str
        First atom/point selection
    selection_b : str
        Second atom/point selection
    """

    label: str
    selection_a: str
    selection_b: str


class CatalyticTriadConfig(BaseModel):
    """Configuration for catalytic triad analysis.

    Attributes
    ----------
    enabled : bool
        Whether to run triad analysis
    name : str
        Name of the catalytic triad/active site
    threshold : float
        Distance threshold for contact/H-bond analysis (Angstroms)
    pairs : list[TriadPairConfig]
        List of atom pairs that define the triad
    """

    enabled: bool = False
    name: str = "catalytic_triad"
    threshold: float = DEFAULT_DISTANCE_THRESHOLD
    pairs: list[TriadPairConfig] = Field(default_factory=list)

    @field_validator("pairs", mode="after")
    @classmethod
    def validate_pairs_if_enabled(cls, v: list[TriadPairConfig], info) -> list[TriadPairConfig]:
        """Warn if enabled but no pairs defined."""
        # Note: We can't access 'enabled' here easily, so validation
        # happens at runtime in the CLI
        return v


class ContactsConfig(BaseModel):
    """Configuration for polymer-protein contact analysis.

    IMPORTANT: PolyzyMD Chain Convention
    ------------------------------------
    - Chain A: Protein/Enzyme
    - Chain B: Substrate/Ligand
    - Chain C: Polymers
    - Chain D+: Solvent (water, ions, co-solvents)

    You MUST use `solvated_system.pdb` as the topology file (NOT
    `production_N_topology.pdb`) to get correct chain assignments.

    Attributes
    ----------
    enabled : bool
        Whether to run contact analysis
    polymer_selection : str
        MDAnalysis selection for polymer atoms. Default uses chain C.
    protein_selection : str
        MDAnalysis selection for protein atoms.
    cutoff : float
        Distance cutoff for contacts in Angstroms
    polymer_types : list[str], optional
        Filter contacts by polymer residue names (e.g., ["SBM", "EGP"]).
        If None, all polymer types are included.
    grouping : str
        How to group protein residues: "aa_class" (aromatic, charged, etc.),
        "secondary_structure", or "none".
    compute_residence_times : bool
        If True, compute residence time statistics for contacts.
    compute_binding_preference : bool
        If True, compute binding preference analysis with enrichment ratios.
        Requires surface exposure calculation using rust_sasa_python.
    polymer_type_selections : dict[str, str], optional
        Define polymer types with MDAnalysis selections.
        Keys are type labels, values are MDAnalysis selection strings.
        Example: {"SBMA": "chainID C and resname SBM", "EGMA": "chainID C and resname EGM"}
    protein_group_selections : dict[str, str], optional
        Define protein groups with MDAnalysis selections.
        Keys are group labels, values are MDAnalysis selection strings.
        Example: {"aromatic": "protein and resname PHE TRP TYR HIS"}
        If not specified, defaults to standard AA class groupings.
    surface_exposure_threshold : float
        Relative SASA threshold (0-1) for surface exposure filtering.
        Residues with SASA/maxSASA > threshold are considered exposed.
        Default 0.2 (20% of max theoretical SASA).
    enzyme_pdb_for_sasa : str, optional
        Path to enzyme PDB for SASA calculation. If not specified,
        uses the enzyme_pdb from the main simulation config.
    """

    enabled: bool = False
    polymer_selection: str = "chainID C"
    protein_selection: str = "protein"
    cutoff: float = DEFAULT_CONTACT_CUTOFF
    polymer_types: Optional[list[str]] = None
    grouping: str = "aa_class"
    compute_residence_times: bool = True

    # Binding preference analysis
    compute_binding_preference: bool = False
    polymer_type_selections: Optional[dict[str, str]] = None
    protein_group_selections: Optional[dict[str, str]] = None
    surface_exposure_threshold: float = DEFAULT_SURFACE_EXPOSURE_THRESHOLD
    enzyme_pdb_for_sasa: Optional[str] = None


# =============================================================================
# Main Analysis Configuration
# =============================================================================


class AnalysisDefaults(BaseModel):
    """Default parameters applied to all analyses.

    Attributes
    ----------
    equilibration_time : str
        Time to skip for equilibration (e.g., "10ns", "5000ps")
    """

    equilibration_time: str = "10ns"


class AnalysisConfig(BaseModel):
    """Schema for analysis.yaml configuration files.

    This configuration defines which analyses to run for a simulation.
    It must be located in the same directory as the simulation's config.yaml.

    Attributes
    ----------
    replicates : list[int]
        List of replicate numbers to analyze
    defaults : AnalysisDefaults
        Default parameters for all analyses
    rmsf : RMSFConfig
        RMSF analysis configuration
    distances : DistancesConfig
        Distance analysis configuration
    catalytic_triad : CatalyticTriadConfig
        Catalytic triad analysis configuration
    contacts : ContactsConfig
        Polymer-protein contact analysis configuration

    Examples
    --------
    >>> config = AnalysisConfig.from_yaml("analysis.yaml")
    >>> if config.rmsf.enabled:
    ...     print("RMSF analysis enabled")
    """

    replicates: list[int] = Field(default_factory=lambda: [1, 2, 3])
    defaults: AnalysisDefaults = Field(default_factory=AnalysisDefaults)
    rmsf: RMSFConfig = Field(default_factory=RMSFConfig)
    distances: DistancesConfig = Field(default_factory=DistancesConfig)
    catalytic_triad: CatalyticTriadConfig = Field(default_factory=CatalyticTriadConfig)
    contacts: ContactsConfig = Field(default_factory=ContactsConfig)

    @field_validator("replicates", mode="before")
    @classmethod
    def ensure_list(cls, v: list[int] | int) -> list[int]:
        """Ensure replicates is a list."""
        if isinstance(v, int):
            return [v]
        return list(v)

    @classmethod
    def from_yaml(cls, path: Path | str) -> "AnalysisConfig":
        """Load analysis config from YAML file.

        Parameters
        ----------
        path : Path or str
            Path to analysis.yaml file

        Returns
        -------
        AnalysisConfig
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
            raise FileNotFoundError(f"Analysis config not found: {path}")

        with open(path) as f:
            data = yaml.safe_load(f)

        # Handle empty file
        if data is None:
            data = {}

        return cls(**data)

    def to_yaml(self, path: Path | str) -> None:
        """Save analysis config to YAML file.

        Parameters
        ----------
        path : Path or str
            Output path for analysis.yaml
        """
        path = Path(path)
        data = self.model_dump(exclude_none=True)

        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    def get_enabled_analyses(self) -> list[str]:
        """Get list of enabled analysis types.

        Returns
        -------
        list[str]
            Names of enabled analyses (e.g., ['rmsf', 'catalytic_triad', 'contacts'])
        """
        enabled = []
        if self.rmsf.enabled:
            enabled.append("rmsf")
        if self.distances.enabled:
            enabled.append("distances")
        if self.catalytic_triad.enabled:
            enabled.append("catalytic_triad")
        if self.contacts.enabled:
            enabled.append("contacts")
        return enabled

    def validate_config(self) -> list[str]:
        """Validate the analysis configuration.

        Returns
        -------
        list[str]
            List of error/warning messages (empty if valid)
        """
        issues = []

        # Check replicates
        if not self.replicates:
            issues.append("No replicates specified")

        # Check distances pairs if enabled
        if self.distances.enabled and not self.distances.pairs:
            issues.append("Distance analysis enabled but no pairs defined")

        # Check triad pairs if enabled
        if self.catalytic_triad.enabled and not self.catalytic_triad.pairs:
            issues.append("Catalytic triad analysis enabled but no pairs defined")

        # Check contacts config if enabled
        if self.contacts.enabled:
            if not self.contacts.polymer_selection:
                issues.append("Contact analysis enabled but no polymer_selection defined")
            if not self.contacts.protein_selection:
                issues.append("Contact analysis enabled but no protein_selection defined")

        return issues


def generate_analysis_template(eq_time: str = "10ns") -> str:
    """Generate a template analysis.yaml file.

    Parameters
    ----------
    eq_time : str
        Default equilibration time

    Returns
    -------
    str
        YAML template content
    """
    return f'''\
# ============================================================================
# PolyzyMD Analysis Configuration
# ============================================================================
# Configures which analyses to run for this simulation.
# Must be in same directory as config.yaml.
#
# Run all enabled analyses: polyzymd analyze run
# Docs: https://polyzymd.readthedocs.io/en/latest/
# ============================================================================

# Which replicates to analyze
replicates: [1, 2, 3]

# Default parameters (override per-analysis if needed)
# Note: Set equilibration_time appropriately for your simulation!
defaults:
  equilibration_time: "{eq_time}"

# ============================================================================
# RMSF Analysis
# ============================================================================
rmsf:
  enabled: true
  selection: "protein and name CA"
  reference_mode: "centroid"  # centroid, average, or frame

# ============================================================================
# Distance Analysis
# ============================================================================
# IMPORTANT: Always use "protein and resid X" for protein residues!
# Residue numbers restart per chain. Without "protein and", your selection
# may match atoms from polymer or water chains, causing incorrect distances.
#
# distances:
#   enabled: true
#   pairs:
#     - label: "Ser77-Substrate"
#       selection_a: "protein and resid 77 and name OG"
#       selection_b: "resname RBY and name C1"

# ============================================================================
# Catalytic Triad Analysis
# ============================================================================
# For cross-condition comparison, define triad in comparison.yaml instead.
#
# IMPORTANT: Always use "protein and resid X" for protein residues!
# Residue numbers restart per chain. Without "protein and", your selection
# may match atoms from polymer or water chains, causing incorrect distances.
#
# catalytic_triad:
#   enabled: true
#   name: "LipA_catalytic_triad"
#   threshold: 3.5
#   pairs:
#     - label: "Asp133-His156"
#       selection_a: "midpoint(protein and resid 133 and name OD1 OD2)"
#       selection_b: "protein and resid 156 and name ND1"
#     - label: "His156-Ser77"
#       selection_a: "protein and resid 156 and name NE2"
#       selection_b: "protein and resid 77 and name OG"

# ============================================================================
# Polymer-Protein Contact Analysis
# ============================================================================
# Analyzes contacts between polymer chains and protein residues.
#
# IMPORTANT: PolyzyMD Chain Convention
#   - Chain A: Protein/Enzyme
#   - Chain B: Substrate/Ligand
#   - Chain C: Polymers
#   - Chain D+: Solvent
# Use solvated_system.pdb as topology (NOT production_N_topology.pdb).
#
# contacts:
#   enabled: true
#   polymer_selection: "chainID C"  # Default: all polymers
#   protein_selection: "protein"
#   cutoff: 4.5                     # Distance cutoff in Angstroms
#   polymer_types: ["SBM", "EGP"]   # Optional: filter by polymer type
#   grouping: "aa_class"            # Group by: aa_class, secondary_structure, none
#   compute_residence_times: true
#
#   # -------------------------------------------------------------------------
#   # Binding Preference Analysis (Enrichment Ratios)
#   # -------------------------------------------------------------------------
#   # Computes enrichment ratios to answer: "Does polymer type X preferentially
#   # bind amino acid class Y?" Requires rust_sasa_python for surface exposure.
#   #
#   # Enrichment > 1.0 = preferential binding
#   # Enrichment = 1.0 = random/expected
#   # Enrichment < 1.0 = avoidance
#   #
#   compute_binding_preference: true
#
#   # Define polymer types with MDAnalysis selections
#   # Each key becomes a bar group in plots
#   polymer_type_selections:
#     SBMA: "chainID C and resname SBM"
#     EGMA: "chainID C and resname EGM"
#
#   # Define protein groups with MDAnalysis selections
#   # Each key becomes a category in enrichment analysis
#   # These defaults classify residues by amino acid physicochemical properties
#   protein_group_selections:
#     aromatic: "protein and resname PHE TRP TYR HIS"
#     polar: "protein and resname SER THR ASN GLN CYS"
#     nonpolar: "protein and resname ALA VAL ILE LEU MET GLY PRO"
#     charged_positive: "protein and resname ARG LYS"
#     charged_negative: "protein and resname ASP GLU"
#     # Add custom groups as needed:
#     # active_site: "protein and resid 77 133 156"
#     # hydrophobic_patch: "protein and resid 12 15 45 67"
#
#   # Surface exposure filtering
#   # Only residues with relative SASA > threshold are considered
#   surface_exposure_threshold: 0.2  # 20% of max theoretical SASA
'''
