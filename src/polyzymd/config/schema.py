"""
Configuration schema for PolyzyMD simulations.

This module defines Pydantic models for all configuration sections,
providing validation, type safety, and YAML/JSON serialization support.
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

from pydantic import BaseModel, Field, field_validator, model_validator


class ChargeMethod(str, Enum):
    """Supported charge assignment methods for small molecules."""

    NAGL = "nagl"
    ESPALOMA = "espaloma"
    AM1BCC = "am1bcc"


class WaterModel(str, Enum):
    """Supported water models."""

    TIP3P = "tip3p"
    SPCE = "spce"
    TIP4P = "tip4p"
    TIP4PEW = "tip4pew"
    OPC = "opc"


class BoxShape(str, Enum):
    """Supported simulation box shapes."""

    CUBE = "cube"
    RHOMBIC_DODECAHEDRON = "rhombic_dodecahedron"
    TRUNCATED_OCTAHEDRON = "truncated_octahedron"


class Ensemble(str, Enum):
    """Thermodynamic ensemble types."""

    NVT = "NVT"
    NPT = "NPT"
    NVE = "NVE"


class ThermostatType(str, Enum):
    """Supported thermostat types."""

    LANGEVIN_MIDDLE = "LangevinMiddle"
    LANGEVIN = "Langevin"
    ANDERSEN = "Andersen"
    NOSE_HOOVER = "NoseHoover"


class BarostatType(str, Enum):
    """Supported barostat types."""

    MONTE_CARLO = "MC"
    MONTE_CARLO_ANISOTROPIC = "MCA"


class RestraintType(str, Enum):
    """Types of restraints that can be applied."""

    FLAT_BOTTOM = "flat_bottom"
    HARMONIC = "harmonic"
    UPPER_WALL = "upper_wall"
    LOWER_WALL = "lower_wall"


# =============================================================================
# Enzyme Configuration
# =============================================================================


class EnzymeConfig(BaseModel):
    """Configuration for the enzyme/protein component.

    Attributes:
        name: Identifier for the enzyme (e.g., "LipA")
        pdb_path: Path to the PDB file containing the enzyme structure
        description: Optional description of the enzyme
    """

    name: str = Field(..., description="Enzyme identifier")
    pdb_path: Path = Field(..., description="Path to enzyme PDB file")
    description: Optional[str] = Field(None, description="Optional description")

    @field_validator("pdb_path")
    @classmethod
    def validate_pdb_path(cls, v: Path) -> Path:
        """Validate that PDB path has correct extension."""
        if v.suffix.lower() != ".pdb":
            raise ValueError(f"Expected .pdb file, got {v.suffix}")
        return v


# =============================================================================
# Substrate/Ligand Configuration
# =============================================================================


class SubstrateConfig(BaseModel):
    """Configuration for the docked substrate/ligand.

    Attributes:
        name: Identifier for the substrate (e.g., "Resorufin-Butyrate")
        sdf_path: Path to SDF file with docked conformers
        conformer_index: Which conformer to use (0-indexed)
        charge_method: Method for assigning partial charges
        residue_name: 3-letter residue name for topology (default: "LIG")
    """

    name: str = Field(..., description="Substrate identifier")
    sdf_path: Path = Field(..., description="Path to docked conformers SDF")
    conformer_index: int = Field(0, ge=0, description="Index of conformer to use")
    charge_method: ChargeMethod = Field(ChargeMethod.NAGL, description="Charge assignment method")
    residue_name: str = Field("LIG", max_length=3, description="3-letter residue name")

    @field_validator("sdf_path")
    @classmethod
    def validate_sdf_path(cls, v: Path) -> Path:
        """Validate that SDF path has correct extension."""
        if v.suffix.lower() != ".sdf":
            raise ValueError(f"Expected .sdf file, got {v.suffix}")
        return v


# =============================================================================
# Polymer Packing Configuration
# =============================================================================


class PolymerPackingConfig(BaseModel):
    """Settings for packing polymers around the solute.

    Controls the box size and PACKMOL behavior when packing polymers
    around the protein-ligand complex.

    Attributes:
        padding: Box padding around the solute in nanometers. Larger values
            give polymers more room and can speed up PACKMOL convergence.
        tolerance: Minimum molecular spacing for PACKMOL in Angstrom.
        movebadrandom: When ``True``, pass the ``movebadrandom`` keyword to
            PACKMOL.  This places badly-packed molecules at random positions
            in the box rather than near well-packed neighbours, which
            improves convergence for dense or heterogeneous systems (many
            unique chain types).  Has no effect when only a single chain
            type is present.  Default is ``False`` (PACKMOL default
            behaviour is preserved).

    Example:
        >>> PolymerPackingConfig(padding=2.5, movebadrandom=True)
    """

    padding: float = Field(2.0, gt=0.0, description="Box padding around solute (nm)")
    tolerance: float = Field(2.0, gt=0.0, description="PACKMOL tolerance (Angstrom)")
    movebadrandom: bool = Field(
        False,
        description=(
            "Pass the 'movebadrandom' keyword to PACKMOL. "
            "Improves convergence for heterogeneous polymer systems "
            "by placing badly-packed molecules at random box positions."
        ),
    )


# =============================================================================
# Polymer Configuration
# =============================================================================


class MonomerSpec(BaseModel):
    """Specification for a single monomer type in a co-polymer.

    For dynamic polymer generation, provide the raw (unactivated) monomer SMILES.
    The system will run initiation reactions to create the active fragments.

    Attributes:
        label: Single character label for this monomer (e.g., "A", "B")
        probability: Probability of selecting this monomer (0-1)
        name: Optional full name (e.g., "SBMA", "EGPMA")
        smiles: Raw monomer SMILES string (required for dynamic generation)
        residue_name: 3-character PDB residue name (auto-generated if not provided)
    """

    label: str = Field(..., min_length=1, max_length=1, description="Monomer label")
    probability: float = Field(..., ge=0.0, le=1.0, description="Selection probability")
    name: Optional[str] = Field(None, description="Full monomer name")
    smiles: Optional[str] = Field(
        None,
        description="Raw monomer SMILES (required for dynamic generation mode)",
    )
    residue_name: Optional[str] = Field(
        None,
        max_length=3,
        description="3-char PDB residue name (auto-generated from name if not provided)",
    )

    @model_validator(mode="after")
    def auto_generate_residue_name(self) -> "MonomerSpec":
        """Auto-generate residue name from monomer name if not provided."""
        if self.residue_name is None and self.name is not None:
            # Use first 3 characters of name, uppercase
            object.__setattr__(self, "residue_name", self.name[:3].upper())
        return self


class PolymerGenerationMode(str, Enum):
    """Mode for polymer generation."""

    CACHED = "cached"  # Load pre-built SDF files from disk
    DYNAMIC = "dynamic"  # Generate polymers on-the-fly using Polymerist


class ReactionConfig(BaseModel):
    """Paths to reaction templates for ATRP polymer generation.

    These .rxn files define the chemical transformations used to create
    polymer fragments from raw monomer SMILES. For ATRP, this includes:
    - Initiation: Activates the vinyl group (e.g., chlorination)
    - Polymerization: Creates chain-extending fragments
    - Termination: Restores the alkene for chain ends

    You can use "default" as a special value to use the bundled ATRP
    methacrylate reaction templates that ship with PolyzyMD.

    Example:
        reactions:
          initiation: "default"
          polymerization: "default"
          termination: "default"

    Attributes:
        initiation: Path to the initiation reaction template (.rxn) or "default"
        polymerization: Path to the polymerization reaction template (.rxn) or "default"
        termination: Path to the termination reaction template (.rxn) or "default"
    """

    initiation: Path = Field(..., description="Path to ATRP initiation .rxn file or 'default'")
    polymerization: Path = Field(
        ..., description="Path to ATRP polymerization .rxn file or 'default'"
    )
    termination: Path = Field(..., description="Path to ATRP termination .rxn file or 'default'")

    @field_validator("initiation", "polymerization", "termination", mode="before")
    @classmethod
    def resolve_default_paths(cls, v, info) -> Path:
        """Resolve 'default' to bundled ATRP reaction paths."""
        from polyzymd.data.reactions import (
            get_atrp_initiation_path,
            get_atrp_polymerization_path,
            get_atrp_termination_path,
        )

        if isinstance(v, str) and v.lower() == "default":
            # Map field name to the appropriate getter
            field_to_getter = {
                "initiation": get_atrp_initiation_path,
                "polymerization": get_atrp_polymerization_path,
                "termination": get_atrp_termination_path,
            }
            getter = field_to_getter.get(info.field_name)
            if getter:
                return getter()
            raise ValueError(f"Unknown reaction field: {info.field_name}")

        # Convert string to Path if needed
        if isinstance(v, str):
            v = Path(v)

        # Validate .rxn extension
        if isinstance(v, Path) and v.suffix.lower() != ".rxn":
            raise ValueError(f"Expected .rxn file, got {v.suffix}")

        return v


class PolymerConfig(BaseModel):
    """Configuration for polymer components.

    Supports two generation modes:
    - "cached": Load pre-built polymer SDF files from sdf_directory (legacy)
    - "dynamic": Generate polymers on-the-fly using Polymerist from SMILES

    For dynamic mode, you must provide:
    - SMILES for each monomer in monomers[].smiles
    - Reaction templates in the reactions field

    Attributes:
        enabled: Whether to include polymers in the system
        generation_mode: "cached" for pre-built SDFs, "dynamic" for on-the-fly generation
        type_prefix: Prefix for polymer type in filenames (e.g., "SBMA-EGPMA")
        monomers: List of monomer specifications with probabilities (and SMILES for dynamic)
        length: Number of monomer units per polymer chain
        count: Number of polymer chains to add
        sdf_directory: Path to pre-built polymer SDF files (for cached mode)
        reactions: Reaction templates for ATRP (required for dynamic mode)
        charger: Charge assignment method for generated polymers
        max_retries: Maximum retries for polymer generation (ring-piercing failures)
        cache_directory: Directory for caching generated polymers and fragments
        packing: Settings for packing polymers around the solute
        random_seed: Random seed for polymer sequence generation (for reproducibility)
    """

    enabled: bool = Field(True, description="Include polymers in system")
    generation_mode: PolymerGenerationMode = Field(
        PolymerGenerationMode.CACHED,
        description="Polymer generation mode: 'cached' (pre-built SDFs) or 'dynamic' (generate from SMILES)",
    )
    type_prefix: str = Field(..., description="Polymer type prefix for filenames")
    monomers: List[MonomerSpec] = Field(..., min_length=1, description="Monomer specifications")
    length: int = Field(..., ge=1, description="Monomers per chain")
    count: int = Field(..., ge=1, description="Number of polymer chains")
    sdf_directory: Optional[Path] = Field(
        None, description="Directory with pre-built polymer SDFs (for cached mode)"
    )
    reactions: Optional[ReactionConfig] = Field(
        None, description="ATRP reaction templates (required for dynamic mode)"
    )
    charger: ChargeMethod = Field(
        ChargeMethod.NAGL,
        description="Charge assignment method for generated polymers",
    )
    max_retries: int = Field(
        10,
        ge=1,
        description="Maximum retries for polymer generation (handles ring-piercing failures)",
    )
    cache_directory: Path = Field(
        Path(".polymer_cache"), description="Cache directory for generated polymers and fragments"
    )
    packing: PolymerPackingConfig = Field(
        default_factory=PolymerPackingConfig,
        description="Polymer packing settings (padding, tolerance)",
    )
    random_seed: Optional[int] = Field(
        None,
        description="Random seed for polymer sequence generation. If None, uses replicate number.",
    )

    @model_validator(mode="after")
    def validate_probabilities_sum_to_one(self) -> "PolymerConfig":
        """Ensure monomer probabilities sum to 1.0."""
        if self.enabled:
            total = sum(m.probability for m in self.monomers)
            if abs(total - 1.0) > 1e-6:
                raise ValueError(f"Monomer probabilities must sum to 1.0, got {total}")
        return self

    @model_validator(mode="after")
    def validate_generation_mode_requirements(self) -> "PolymerConfig":
        """Validate that required fields are present for the selected generation mode."""
        if not self.enabled:
            return self

        if self.generation_mode == PolymerGenerationMode.DYNAMIC:
            # Dynamic mode requires SMILES for all monomers
            missing_smiles = [m.label for m in self.monomers if m.smiles is None]
            if missing_smiles:
                raise ValueError(
                    f"Dynamic generation mode requires 'smiles' for all monomers. "
                    f"Missing SMILES for monomers: {missing_smiles}"
                )

            # Dynamic mode requires reaction templates
            if self.reactions is None:
                raise ValueError(
                    "Dynamic generation mode requires 'reactions' field with ATRP reaction templates"
                )

        elif self.generation_mode == PolymerGenerationMode.CACHED:
            # Cached mode requires sdf_directory
            if self.sdf_directory is None:
                raise ValueError(
                    "Cached generation mode requires 'sdf_directory' with pre-built polymer SDFs"
                )

        return self


# =============================================================================
# Solvent Configuration
# =============================================================================


class CoSolventSpec(BaseModel):
    """Specification for a co-solvent component.

    You must specify EITHER volume_fraction OR concentration, not both.

    For co-solvents in the built-in library (dmso, dmf, urea, ethanol, etc.),
    you can omit the smiles and density fields - they will be looked up
    automatically.

    Attributes:
        name: Identifier for the co-solvent (e.g., "dmso")
        smiles: SMILES string (optional if co-solvent is in library)
        volume_fraction: Volume fraction (0-1), e.g., 0.30 for 30% v/v
        concentration: Molar concentration (mol/L)
        density: Density in g/mL (required for volume_fraction with custom molecules)
        residue_name: 3-letter residue name (default: first 3 chars of name)

    Example (library co-solvent with volume fraction):
        >>> CoSolventSpec(name="dmso", volume_fraction=0.30)

    Example (library co-solvent with concentration):
        >>> CoSolventSpec(name="urea", concentration=2.0)

    Example (custom co-solvent):
        >>> CoSolventSpec(
        ...     name="my_solvent",
        ...     smiles="CCOC(=O)C",
        ...     density=0.902,
        ...     volume_fraction=0.15
        ... )
    """

    name: str = Field(..., description="Co-solvent identifier")
    smiles: Optional[str] = Field(None, description="SMILES string (optional if in library)")

    # Specification method 1: Volume fraction
    volume_fraction: Optional[float] = Field(
        None,
        gt=0.0,
        lt=1.0,
        description="Volume fraction (0-1), e.g., 0.30 for 30% v/v",
    )

    # Specification method 2: Molar concentration
    concentration: Optional[float] = Field(None, gt=0.0, description="Molar concentration (mol/L)")

    # Physical property for volume fraction calculation
    density: Optional[float] = Field(
        None,
        gt=0.0,
        description="Density in g/mL (required for volume_fraction with custom molecules)",
    )

    residue_name: Optional[str] = Field(None, max_length=3, description="Residue name")

    @model_validator(mode="after")
    def validate_and_populate(self) -> "CoSolventSpec":
        """Validate specification and populate missing fields from library."""
        from polyzymd.data.cosolvent_library import get_cosolvent

        # Check that exactly one of volume_fraction or concentration is specified
        has_vf = self.volume_fraction is not None
        has_conc = self.concentration is not None

        if not has_vf and not has_conc:
            raise ValueError(
                f"Co-solvent '{self.name}': Must specify either 'volume_fraction' "
                f"or 'concentration'"
            )
        if has_vf and has_conc:
            raise ValueError(
                f"Co-solvent '{self.name}': Cannot specify both 'volume_fraction' "
                f"and 'concentration' - choose one"
            )

        # Look up from library
        library_data = get_cosolvent(self.name)

        # Populate SMILES if not provided
        if self.smiles is None:
            if library_data:
                object.__setattr__(self, "smiles", library_data.smiles)
            else:
                raise ValueError(
                    f"Co-solvent '{self.name}' not in library. Please provide 'smiles' field."
                )

        # For volume_fraction, we need density
        if has_vf and self.density is None:
            if library_data:
                object.__setattr__(self, "density", library_data.density)
            else:
                raise ValueError(
                    f"Co-solvent '{self.name}' not in library. "
                    f"Please provide 'density' field (g/mL) for volume_fraction calculation."
                )

        # Set default residue name
        if self.residue_name is None:
            object.__setattr__(self, "residue_name", self.name[:3].upper())

        return self


class PrimarySolventConfig(BaseModel):
    """Configuration for the primary solvent (usually water).

    Attributes:
        type: Solvent type identifier
        model: Water model to use (if type is "water")
    """

    type: str = Field("water", description="Primary solvent type")
    model: WaterModel = Field(WaterModel.TIP3P, description="Water model")


class IonConfig(BaseModel):
    """Configuration for ions in the solvent.

    Attributes:
        neutralize: Whether to add ions for charge neutralization
        nacl_concentration: Additional NaCl concentration in mol/L
        kcl_concentration: Additional KCl concentration in mol/L
        mgcl2_concentration: Additional MgCl2 concentration in mol/L
    """

    neutralize: bool = Field(True, description="Neutralize system charge")
    nacl_concentration: float = Field(0.1, ge=0.0, description="NaCl conc. (mol/L)")
    kcl_concentration: float = Field(0.0, ge=0.0, description="KCl conc. (mol/L)")
    mgcl2_concentration: float = Field(0.0, ge=0.0, description="MgCl2 conc. (mol/L)")


class BoxConfig(BaseModel):
    """Configuration for the simulation box.

    Attributes:
        padding: Distance from solute to box edge in nm
        shape: Box geometry
        target_density: Target density in g/mL
        tolerance: Minimum molecular spacing for PACKMOL in Angstrom
    """

    padding: float = Field(1.2, gt=0.0, description="Box padding (nm)")
    shape: BoxShape = Field(BoxShape.RHOMBIC_DODECAHEDRON, description="Box shape")
    target_density: float = Field(1.0, gt=0.0, description="Target density (g/mL)")
    tolerance: float = Field(2.0, gt=0.0, description="PACKMOL tolerance (Angstrom)")


class SolventConfig(BaseModel):
    """Complete solvent configuration.

    Attributes:
        primary: Primary solvent settings
        co_solvents: List of co-solvent specifications
        ions: Ion configuration
        box: Box geometry settings
    """

    primary: PrimarySolventConfig = Field(
        default_factory=PrimarySolventConfig, description="Primary solvent"
    )
    co_solvents: List[CoSolventSpec] = Field(default_factory=list, description="Co-solvents")
    ions: IonConfig = Field(default_factory=IonConfig, description="Ion settings")
    box: BoxConfig = Field(default_factory=BoxConfig, description="Box settings")

    @model_validator(mode="after")
    def validate_volume_fractions(self) -> "SolventConfig":
        """Ensure co-solvent volume fractions don't exceed 1.0."""
        total = sum(cs.volume_fraction for cs in self.co_solvents)
        if total >= 1.0:
            raise ValueError(f"Total co-solvent volume fraction must be < 1.0, got {total}")
        return self


# =============================================================================
# Restraint Configuration
# =============================================================================


class AtomSelectionConfig(BaseModel):
    """Configuration for selecting atoms for restraints.

    Uses MDAnalysis-compatible selection syntax for flexibility.

    Attributes:
        selection: MDAnalysis selection string (e.g., "resid 77 and name OG")
        description: Optional human-readable description
    """

    selection: str = Field(..., description="MDAnalysis selection string")
    description: Optional[str] = Field(None, description="Human-readable description")


class RestraintConfig(BaseModel):
    """Configuration for a single restraint.

    Attributes:
        type: Type of restraint (flat_bottom, harmonic, etc.)
        name: Identifier for this restraint
        atom1: First atom selection
        atom2: Second atom selection
        distance: Target/threshold distance in Angstroms
        force_constant: Force constant in kJ/mol/nm^2
        enabled: Whether this restraint is active
    """

    type: RestraintType = Field(..., description="Restraint type")
    name: str = Field(..., description="Restraint identifier")
    atom1: AtomSelectionConfig = Field(..., description="First atom selection")
    atom2: AtomSelectionConfig = Field(..., description="Second atom selection")
    distance: float = Field(..., gt=0.0, description="Distance threshold (Angstrom)")
    force_constant: float = Field(10000.0, gt=0.0, description="Force constant (kJ/mol/nm^2)")
    enabled: bool = Field(True, description="Whether restraint is active")


# =============================================================================
# Thermodynamics Configuration
# =============================================================================


class ThermodynamicsConfig(BaseModel):
    """Thermodynamic conditions for the simulation.

    Attributes:
        temperature: System temperature in Kelvin
        pressure: System pressure in atmospheres (for NPT)
        salt_concentration: Ionic strength in mol/L (deprecated, use solvent.ions)
    """

    temperature: float = Field(..., gt=0.0, description="Temperature (K)")
    pressure: float = Field(1.0, gt=0.0, description="Pressure (atm)")


# =============================================================================
# Simulation Phase Configuration
# =============================================================================


class SimulationPhaseConfig(BaseModel):
    """Configuration for a single simulation phase (equilibration or production).

    Attributes:
        ensemble: Thermodynamic ensemble (NVT, NPT)
        duration: Simulation duration in nanoseconds
        samples: Number of trajectory frames to save
        time_step: Integration time step in femtoseconds
        thermostat: Thermostat type
        thermostat_timescale: Thermostat coupling timescale in ps
        barostat: Barostat type (for NPT)
        barostat_frequency: Barostat update frequency (steps)
    """

    ensemble: Ensemble = Field(..., description="Thermodynamic ensemble")
    duration: float = Field(..., gt=0.0, description="Duration (ns)")
    samples: int = Field(..., ge=1, description="Trajectory frames to save")
    time_step: float = Field(2.0, gt=0.0, description="Time step (fs)")
    thermostat: ThermostatType = Field(
        ThermostatType.LANGEVIN_MIDDLE, description="Thermostat type"
    )
    thermostat_timescale: float = Field(1.0, gt=0.0, description="Thermostat timescale (ps)")
    barostat: Optional[BarostatType] = Field(None, description="Barostat type")
    barostat_frequency: int = Field(25, ge=1, description="Barostat update frequency")

    @model_validator(mode="after")
    def validate_ensemble_barostat(self) -> "SimulationPhaseConfig":
        """Ensure NPT ensemble has a barostat configured."""
        if self.ensemble == Ensemble.NPT and self.barostat is None:
            self.barostat = BarostatType.MONTE_CARLO
        return self


# =============================================================================
# Equilibration Stage Configuration (Multi-Stage Equilibration)
# =============================================================================


class PositionRestraintConfig(BaseModel):
    """Configuration for positional restraints on an atom group.

    Position restraints apply a harmonic potential to keep atoms near their
    initial coordinates. This is commonly used during equilibration to
    prevent large structural changes while the system relaxes.

    Attributes:
        group: Predefined atom group name
        force_constant: Force constant in kJ/mol/nm^2 (4184.0 = 1.0 kcal/mol/A^2)
    """

    group: str = Field(
        ...,
        description=(
            "Atom group: protein_heavy, protein_backbone, protein_calpha, "
            "ligand_heavy, polymer_heavy, solvent, water_only, ions_only, cosolvents_only"
        ),
    )
    force_constant: float = Field(
        4184.0,  # 1.0 kcal/mol/A^2 in kJ/mol/nm^2
        gt=0.0,
        description="Force constant (kJ/mol/nm^2). Default 4184.0 = 1.0 kcal/mol/A^2",
    )

    @field_validator("group")
    @classmethod
    def validate_group_name(cls, v: str) -> str:
        """Validate that the group name is a recognized predefined group."""
        from polyzymd.core.atom_groups import PREDEFINED_GROUPS

        if v not in PREDEFINED_GROUPS:
            raise ValueError(
                f"Unknown atom group: '{v}'. Valid groups: {sorted(PREDEFINED_GROUPS)}"
            )
        return v


class EquilibrationStageConfig(BaseModel):
    """Configuration for a single equilibration stage.

    Supports two temperature modes:
    1. Constant temperature: Set 'temperature' field
    2. Temperature ramping (simulated annealing): Set 'temperature_start' and
       'temperature_end' fields

    Position restraints can be applied to hold specific atom groups in place
    during the stage.

    Attributes:
        name: Stage identifier (used in output paths)
        duration: Stage duration in nanoseconds
        samples: Number of trajectory frames to save
        ensemble: Thermodynamic ensemble (NVT or NPT)
        temperature: Constant temperature in K (mutually exclusive with ramping)
        temperature_start: Starting temperature for ramping in K
        temperature_end: Ending temperature for ramping in K
        temperature_increment: Temperature increment per update in K
        temperature_interval: Time between temperature updates in fs
        position_restraints: List of position restraints for this stage
        time_step: Optional time step override in fs
        thermostat: Optional thermostat type override
        thermostat_timescale: Optional thermostat timescale override in ps
        barostat: Optional barostat type (for NPT ensemble)
        barostat_frequency: Optional barostat update frequency
    """

    name: str = Field(..., description="Stage identifier (used in output paths)")
    duration: float = Field(..., gt=0.0, description="Duration (ns)")
    samples: int = Field(100, ge=1, description="Trajectory frames to save")
    ensemble: Ensemble = Field(Ensemble.NVT, description="Thermodynamic ensemble")

    # Constant temperature mode
    temperature: Optional[float] = Field(
        None,
        gt=0.0,
        description="Constant temperature (K). Mutually exclusive with temperature ramping.",
    )

    # Temperature ramping mode (simulated annealing)
    temperature_start: Optional[float] = Field(
        None, gt=0.0, description="Starting temperature for ramping (K)"
    )
    temperature_end: Optional[float] = Field(
        None, gt=0.0, description="Ending temperature for ramping (K)"
    )
    temperature_increment: float = Field(
        1.0, gt=0.0, description="Temperature increment per update (K)"
    )
    temperature_interval: float = Field(
        1200.0, gt=0.0, description="Time between temperature updates (fs)"
    )

    # Position restraints
    position_restraints: List[PositionRestraintConfig] = Field(
        default_factory=list, description="Position restraints for this stage"
    )

    # Optional overrides (inherit from parent/defaults if not specified)
    time_step: Optional[float] = Field(None, gt=0.0, description="Time step (fs)")
    thermostat: Optional[ThermostatType] = Field(None, description="Thermostat type")
    thermostat_timescale: Optional[float] = Field(
        None, gt=0.0, description="Thermostat timescale (ps)"
    )
    barostat: Optional[BarostatType] = Field(None, description="Barostat type (for NPT)")
    barostat_frequency: Optional[int] = Field(None, ge=1, description="Barostat update frequency")

    @model_validator(mode="after")
    def validate_temperature_mode(self) -> "EquilibrationStageConfig":
        """Ensure valid temperature specification."""
        has_constant = self.temperature is not None
        has_start = self.temperature_start is not None
        has_end = self.temperature_end is not None

        if has_constant and (has_start or has_end):
            raise ValueError(
                "Cannot specify both 'temperature' and temperature ramping "
                "('temperature_start'/'temperature_end')"
            )

        if not has_constant and not (has_start and has_end):
            raise ValueError(
                "Must specify either 'temperature' for constant temperature, "
                "or both 'temperature_start' and 'temperature_end' for ramping"
            )

        if has_start and has_end and self.temperature_start > self.temperature_end:
            raise ValueError(
                f"temperature_start ({self.temperature_start}) must be <= "
                f"temperature_end ({self.temperature_end})"
            )

        return self

    @model_validator(mode="after")
    def validate_npt_barostat(self) -> "EquilibrationStageConfig":
        """Ensure NPT ensemble has a barostat configured."""
        if self.ensemble == Ensemble.NPT and self.barostat is None:
            self.barostat = BarostatType.MONTE_CARLO
        return self

    @property
    def is_temperature_ramping(self) -> bool:
        """Check if this stage uses temperature ramping."""
        return self.temperature_start is not None

    def get_start_temperature(self) -> float:
        """Get the starting temperature of this stage."""
        if self.is_temperature_ramping:
            return self.temperature_start
        return self.temperature

    def get_final_temperature(self) -> float:
        """Get the final temperature of this stage."""
        if self.is_temperature_ramping:
            return self.temperature_end
        return self.temperature


class SimulationPhasesConfig(BaseModel):
    """Configuration for all simulation phases.

    Supports two equilibration modes (mutually exclusive):
    1. Simple mode (legacy): Single equilibration phase via 'equilibration' field
    2. Staged mode: Multi-stage protocol via 'equilibration_stages' field

    Attributes:
        equilibration_stages: Multi-stage equilibration protocol (new)
        equilibration: Simple single-stage equilibration (legacy)
        production: Production phase settings
        segments: Number of segments for daisy-chaining
    """

    # Staged equilibration (new)
    equilibration_stages: Optional[List[EquilibrationStageConfig]] = Field(
        None, description="Multi-stage equilibration protocol with position restraints"
    )

    # Simple equilibration (legacy, for backwards compatibility)
    equilibration: Optional[SimulationPhaseConfig] = Field(
        None, description="Simple single-stage equilibration (legacy)"
    )

    production: SimulationPhaseConfig = Field(..., description="Production settings")
    segments: int = Field(1, ge=1, description="Number of daisy-chain segments")

    @model_validator(mode="after")
    def validate_equilibration_mode(self) -> "SimulationPhasesConfig":
        """Ensure exactly one equilibration mode is specified."""
        has_stages = self.equilibration_stages is not None
        has_simple = self.equilibration is not None

        if has_stages and has_simple:
            raise ValueError(
                "Cannot specify both 'equilibration_stages' and 'equilibration'. "
                "Use 'equilibration_stages' for multi-stage protocols, "
                "or 'equilibration' for simple single-stage equilibration."
            )

        if not has_stages and not has_simple:
            raise ValueError("Must specify either 'equilibration_stages' or 'equilibration'")

        return self

    @property
    def uses_staged_equilibration(self) -> bool:
        """Check if using multi-stage equilibration."""
        return self.equilibration_stages is not None

    @property
    def total_equilibration_duration(self) -> float:
        """Total equilibration duration in nanoseconds.

        Works for both simple and staged equilibration modes.
        For staged mode, returns the sum of all stage durations.
        """
        if self.uses_staged_equilibration:
            return sum(stage.duration for stage in self.equilibration_stages)
        return self.equilibration.duration

    @property
    def total_equilibration_samples(self) -> int:
        """Total equilibration trajectory samples.

        Works for both simple and staged equilibration modes.
        For staged mode, returns the sum of all stage samples.
        """
        if self.uses_staged_equilibration:
            return sum(stage.samples for stage in self.equilibration_stages)
        return self.equilibration.samples


# =============================================================================
# Output Configuration
# =============================================================================


def expand_path(path: Path) -> Path:
    """Expand environment variables and user home in a path.

    Supports:
    - $VAR and ${VAR} syntax for environment variables
    - ~ for user home directory

    Args:
        path: Path that may contain environment variables

    Returns:
        Path with variables expanded
    """
    import os

    expanded = os.path.expandvars(str(path))
    expanded = os.path.expanduser(expanded)
    return Path(expanded)


class OutputConfig(BaseModel):
    """Configuration for simulation output.

    Supports separate directories for:
    - scripts/logs (projects_directory): Where job scripts and SLURM logs are written
    - simulation data (scratch_directory): Where trajectories, checkpoints go

    This separation allows running simulations on HPC systems where code lives
    in long-term storage (projects) but data is written to high-performance
    scratch storage.

    Attributes:
        projects_directory: Directory for scripts, configs, logs (long-term storage)
        scratch_directory: Directory for simulation output (high-performance storage)
        naming_template: Template for naming working directories
        job_scripts_subdir: Subdirectory name for job scripts within projects
        slurm_logs_subdir: Subdirectory name for SLURM logs within projects
        save_checkpoint: Whether to save checkpoint files
        save_state_data: Whether to save thermodynamic state data
        trajectory_format: Output trajectory format

    Example YAML:
        output:
          projects_directory: /projects/user/polyzymd
          scratch_directory: /scratch/alpine/user/simulations
          naming_template: "{enzyme}_{substrate}_{temperature}K_run{replicate}"
    """

    # Directory configuration
    projects_directory: Path = Field(
        Path("."),
        description="Base directory for scripts and logs (typically in projects/long-term storage)",
    )
    scratch_directory: Optional[Path] = Field(
        None,
        description="Base directory for simulation output (scratch storage). If None, uses projects_directory.",
    )
    naming_template: str = Field(
        "{enzyme}_{substrate}_{polymer_type}_{duration}ns_{temperature}K_run{replicate}",
        description="Directory naming template",
    )

    # Subdirectory organization
    job_scripts_subdir: str = Field(
        "job_scripts",
        description="Subdirectory for generated job scripts",
    )
    slurm_logs_subdir: str = Field(
        "slurm_logs",
        description="Subdirectory for SLURM output logs",
    )

    # Output options
    save_checkpoint: bool = Field(True, description="Save checkpoint files")
    save_state_data: bool = Field(True, description="Save state data CSV")
    trajectory_format: str = Field("dcd", description="Trajectory file format")

    # Legacy compatibility
    base_directory: Optional[Path] = Field(
        None,
        description="Deprecated: Use scratch_directory instead",
        exclude=True,
    )

    @field_validator("projects_directory", "scratch_directory", "base_directory", mode="before")
    @classmethod
    def expand_env_vars_in_paths(cls, v: Optional[Union[str, Path]]) -> Optional[Path]:
        """Expand environment variables and ~ in path fields.

        Supports $USER, ${HOME}, ~/path, etc.
        """
        if v is None:
            return None
        return expand_path(Path(v))

    @model_validator(mode="after")
    def handle_legacy_base_directory(self) -> "OutputConfig":
        """Handle legacy base_directory field for backwards compatibility."""
        if self.base_directory is not None and self.scratch_directory is None:
            self.scratch_directory = self.base_directory
        return self

    @property
    def effective_scratch_directory(self) -> Path:
        """Get the effective scratch directory (falls back to projects if not set)."""
        return self.scratch_directory if self.scratch_directory else self.projects_directory

    def format_directory_name(
        self,
        enzyme: str,
        substrate: str,
        polymer_type: str,
        temperature: float,
        replicate: int,
        duration: float = 0.0,
        **kwargs: Any,
    ) -> str:
        """Format the directory name using the template.

        Args:
            enzyme: Enzyme name
            substrate: Substrate name
            polymer_type: Polymer type (or "none")
            temperature: Temperature in K
            replicate: Replicate number
            duration: Production duration in ns
            **kwargs: Additional template variables

        Returns:
            Formatted directory name
        """
        return self.naming_template.format(
            enzyme=enzyme,
            substrate=substrate,
            polymer_type=polymer_type,
            temperature=int(temperature),
            replicate=replicate,
            duration=int(duration),
            **kwargs,
        )

    def get_job_scripts_directory(self) -> Path:
        """Get the directory for job scripts.

        Returns:
            Path to job scripts directory (within projects)
        """
        return self.projects_directory / self.job_scripts_subdir

    def get_slurm_logs_directory(self) -> Path:
        """Get the directory for SLURM log files.

        Returns:
            Path to SLURM logs directory (within projects)
        """
        return self.projects_directory / self.slurm_logs_subdir


# =============================================================================
# Force Field Configuration
# =============================================================================


class ForceFieldConfig(BaseModel):
    """Configuration for force field selection.

    Attributes:
        protein: Force field for proteins
        small_molecule: Force field for small molecules
        water: Water model force field (derived from solvent config)
    """

    protein: str = Field("ff14sb_off_impropers_0.0.4.offxml", description="Protein force field")
    small_molecule: str = Field("openff-2.0.0.offxml", description="Small molecule force field")


# =============================================================================
# Main Simulation Configuration
# =============================================================================


class SimulationConfig(BaseModel):
    """Complete simulation configuration.

    This is the top-level configuration model that contains all settings
    for a PolyzyMD simulation.

    Attributes:
        name: Simulation name/identifier
        description: Optional description
        enzyme: Enzyme configuration
        substrate: Substrate/ligand configuration
        polymers: Polymer configuration (optional)
        solvent: Solvent and box configuration
        restraints: List of restraint configurations
        thermodynamics: Temperature/pressure settings
        simulation_phases: Equilibration and production settings
        output: Output file settings
        force_field: Force field selection

    Example:
        >>> config = SimulationConfig.from_yaml("simulation.yaml")
        >>> print(config.enzyme.name)
        "LipA"
    """

    name: str = Field(..., description="Simulation identifier")
    description: Optional[str] = Field(None, description="Simulation description")

    enzyme: EnzymeConfig = Field(..., description="Enzyme configuration")
    substrate: Optional[SubstrateConfig] = Field(None, description="Substrate configuration")
    polymers: Optional[PolymerConfig] = Field(None, description="Polymer configuration")
    solvent: SolventConfig = Field(
        default_factory=SolventConfig, description="Solvent configuration"
    )
    restraints: List[RestraintConfig] = Field(
        default_factory=list, description="Restraint configurations"
    )
    thermodynamics: ThermodynamicsConfig = Field(..., description="Thermodynamic conditions")
    simulation_phases: SimulationPhasesConfig = Field(..., description="Phase settings")
    output: OutputConfig = Field(default_factory=OutputConfig, description="Output settings")
    force_field: ForceFieldConfig = Field(
        default_factory=ForceFieldConfig, description="Force field settings"
    )

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "SimulationConfig":
        """Load configuration from a YAML file.

        Args:
            path: Path to the YAML configuration file

        Returns:
            SimulationConfig instance

        Raises:
            FileNotFoundError: If file doesn't exist
            ValidationError: If configuration is invalid
        """
        from polyzymd.config.loader import load_config

        return load_config(path)

    def to_yaml(self, path: Union[str, Path]) -> None:
        """Save configuration to a YAML file.

        Args:
            path: Path to save the configuration file
        """
        from polyzymd.config.loader import save_config

        save_config(self, path)

    def get_working_directory(self, replicate: int = 1) -> Path:
        """Get the working directory path for a given replicate (in scratch).

        This returns the path where simulation output (trajectories, checkpoints)
        will be written, which is in the scratch directory.

        Args:
            replicate: Replicate number

        Returns:
            Path to the working directory (in scratch)
        """
        dir_name = self._format_run_directory_name(replicate)
        return self.output.effective_scratch_directory / dir_name

    def get_projects_directory(self) -> Path:
        """Get the projects directory path.

        This is where scripts, configs, and logs are stored.

        Returns:
            Path to the projects directory
        """
        return self.output.projects_directory

    def _format_run_directory_name(self, replicate: int) -> str:
        """Format the run directory name for a replicate.

        Args:
            replicate: Replicate number

        Returns:
            Formatted directory name
        """
        polymer_type = "none"
        if self.polymers and self.polymers.enabled:
            # Format polymer type with full composition info
            # Sort by label (A, B, C...) alphabetically - user controls which is "A"
            probs = {m.label: m.probability for m in self.polymers.monomers}
            sorted_labels = sorted(probs.keys())
            composition = "_".join(f"{lbl}{probs[lbl] * 100:.0f}" for lbl in sorted_labels)
            polymer_type = f"{self.polymers.type_prefix}_{composition}"

        substrate_name = self.substrate.name if self.substrate else "apo"

        return self.output.format_directory_name(
            enzyme=self.enzyme.name,
            substrate=substrate_name.replace("-", ""),
            polymer_type=polymer_type,
            temperature=self.thermodynamics.temperature,
            replicate=replicate,
            duration=self.simulation_phases.production.duration,
        )

    def to_signac_statepoint(self, replicate: int = 1) -> Dict[str, Any]:
        """Convert configuration to a Signac-compatible state point dictionary.

        Args:
            replicate: Replicate number

        Returns:
            Dictionary suitable for use as a Signac state point
        """
        statepoint = {
            "enzyme": self.enzyme.name,
            "temperature": self.thermodynamics.temperature,
            "replicate": replicate,
        }

        if self.substrate:
            statepoint["substrate"] = self.substrate.name
            statepoint["substrate_conformer"] = self.substrate.conformer_index

        if self.polymers and self.polymers.enabled:
            statepoint["polymer_type"] = self.polymers.type_prefix
            statepoint["polymer_length"] = self.polymers.length
            statepoint["polymer_count"] = self.polymers.count
            # Include monomer probabilities
            for monomer in self.polymers.monomers:
                statepoint[f"monomer_{monomer.label}_prob"] = monomer.probability
        else:
            statepoint["polymer_type"] = "none"

        # Include co-solvent info
        for cosolvent in self.solvent.co_solvents:
            statepoint[f"cosolvent_{cosolvent.name}_fraction"] = cosolvent.volume_fraction

        return statepoint
