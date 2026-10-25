"""
Configuration schema for PolyzyMD simulations.

This module defines Pydantic models for all configuration sections,
providing validation, type safety, and YAML/JSON serialization support.
"""

from __future__ import annotations

import functools
import logging
import math
import os
import re
import shlex
import warnings
from enum import Enum
from pathlib import Path
from string import Formatter
from typing import Annotated, Any, ClassVar, Literal, Self, Union, get_args

from pydantic import (
    BaseModel,
    ConfigDict,
    Discriminator,
    Field,
    SerializerFunctionWrapHandler,
    Tag,
    field_validator,
    model_serializer,
    model_validator,
)

LOGGER = logging.getLogger(__name__)
#: Warnings about unread keys printed in this process, so each prints once.
_WARNED_UNREAD: set[str] = set()

_TOKEN_UNSAFE_CHARS = re.compile(r"[^A-Za-z0-9_-]+")
_TOKEN_UNDERSCORES = re.compile(r"_+")


def _format_safe_token(value: object) -> str:
    """Convert a value into a safe lowercase naming token.

    Parameters
    ----------
    value : object
        Value to normalize for a directory or job-name token.

    Returns
    -------
    str
        Sanitized token using only alphanumerics, underscores, and hyphens.
    """
    token = str(value).strip().lower()
    token = _TOKEN_UNSAFE_CHARS.sub("_", token)
    token = _TOKEN_UNDERSCORES.sub("_", token)
    return token.strip("_") or "unknown"


def _format_decimal_token(value: float) -> str:
    """Format a decimal value without raw decimal points.

    Parameters
    ----------
    value : float
        Numeric value to format.

    Returns
    -------
    str
        Compact token, using ``p`` as the decimal separator when needed.
    """
    formatted = f"{value:g}"
    return formatted.replace(".", "p")


class ChargeMethod(str, Enum):
    """Supported charge assignment methods for small molecules."""

    NAGL = "nagl"
    ESPALOMA = "espaloma"
    AM1BCC = "am1bcc"


class WaterModel(str, Enum):
    """Water models a config may name.

    Only TIP3P is built: the force fields PolyzyMD loads carry the TIP3P water
    parameters. The other values stay loadable so configs of existing runs
    still load for status and analysis; new builds refuse them
    (:meth:`SimulationConfig.require_buildable`).
    """

    TIP3P = "tip3p"
    SPCE = "spce"
    TIP4P = "tip4p"
    TIP4PEW = "tip4pew"
    OPC = "opc"


class BoxShape(str, Enum):
    """Supported simulation box shapes."""

    CUBE = "cube"
    RHOMBIC_DODECAHEDRON = "rhombic_dodecahedron"


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


class _ConfigModel(BaseModel):
    """Base of every config section: a key the section does not define is an error.

    A misspelled or unsupported key (``chrage_method``) would otherwise be
    dropped without a word, and the build would use the default. Numbers must
    be finite, and ``true``/``false`` is not a number.
    """

    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)

    #: Keys a section still accepts only so the configs of existing runs load.
    #: No code reads them. Giving one prints its warning once per process, and
    #: ``save_config`` leaves out those at their default value; the config hash
    #: still includes them, so a saved config keeps its hash.
    UNREAD_FIELDS: ClassVar[dict[str, str]] = {}

    @model_validator(mode="after")
    def _warn_unread_fields(self) -> Self:
        """Warn once for each unread key the config gives."""
        for name in sorted(self.model_fields_set & self.UNREAD_FIELDS.keys()):
            message = self.UNREAD_FIELDS[name]
            if message not in _WARNED_UNREAD:
                _WARNED_UNREAD.add(message)
                LOGGER.warning(message)
        return self

    @field_validator("*", mode="before")
    @classmethod
    def _refuse_bool_as_number(cls, value: Any, info: Any) -> Any:
        """Refuse ``true``/``false`` where a number is expected (it would read as 1 or 0)."""
        if isinstance(value, bool) and info.field_name in cls.model_fields:
            annotation = cls.model_fields[info.field_name].annotation
            kinds = get_args(annotation) or (annotation,)
            if bool not in kinds and (int in kinds or float in kinds):
                raise ValueError(f"must be a number, not {str(value).lower()}")
        return value


def load_custom_substructures(path: Path) -> dict[str, dict[str, list[str]]]:
    """Read a custom substructures JSON file and check its shape.

    Parameters
    ----------
    path : Path
        JSON file that maps each residue name to SMARTS patterns, and each
        SMARTS pattern to the PDB atom names it matches, in order.

    Returns
    -------
    dict
        The file's content, ``{residue name: {SMARTS: [atom names]}}``.

    Raises
    ------
    ValueError
        When the file is missing, is not JSON or has another shape.
    """
    import json

    try:
        data = json.loads(Path(path).read_text())
    except OSError as error:
        raise ValueError(f"Cannot read custom_substructures_path {path}: {error}") from error
    except json.JSONDecodeError as error:
        raise ValueError(f"custom_substructures_path {path} is not JSON: {error}") from error
    if not isinstance(data, dict) or not all(
        isinstance(name, str)
        and isinstance(patterns, dict)
        and all(
            isinstance(smarts, str)
            and isinstance(atoms, list)
            and all(isinstance(atom, str) for atom in atoms)
            for smarts, atoms in patterns.items()
        )
        for name, patterns in data.items()
    ):
        raise ValueError(
            f"custom_substructures_path {path} must map each residue name to SMARTS "
            'patterns, and each pattern to a list of atom names: {"NCYX": {"[#7:1]...": ["N", ...]}}'
        )
    return data


#: The page that maps every key of a 1.3.0rc4 config to its place in rc5.
REWRITE_GUIDE = "https://polyzymd.readthedocs.io/en/latest/how_to/rewrite_config_rc5.html"

#: Same default as polyzymd.utils.charging.DEFAULT_NAGL_MODEL (not imported: it loads RDKit).
DEFAULT_NAGL_MODEL = "openff-gnn-am1bcc-1.0.0.pt"


class ChargeSpec(_ConfigModel):
    """Where a component's partial charges come from.

    - ``library``: the force field's library charges (ff14SB for proteins).
    - ``preset``: the charges stored in the component's SDF files.
    - ``nagl``: OpenFF NAGL; ``model`` names the network file (required).
    - ``am1bcc`` (needs AmberTools) and ``espaloma``: computed for each molecule.

    PolyzyMD assumes no charge model: the config names one for every component,
    so the choice is on record with the system it made.
    """

    method: Literal["library", "preset", "nagl", "am1bcc", "espaloma"] = Field(
        ..., description="library, preset, nagl, am1bcc or espaloma"
    )
    model: str | None = Field(
        None,
        description=(
            "NAGL model file, required for method nagl (for example "
            f"{DEFAULT_NAGL_MODEL}, the model Sage 2.3.0 names)"
        ),
    )

    @model_validator(mode="after")
    def _model_only_for_nagl(self) -> Self:
        """Require a NAGL model with nagl, and refuse one for any other method."""
        if self.method == "nagl":
            if self.model is None:
                raise ValueError(
                    "method nagl needs model: the NAGL network file, for example "
                    f"model: {DEFAULT_NAGL_MODEL} (the model Sage 2.3.0 names)"
                )
        elif self.model is not None:
            raise ValueError(
                f"model names a NAGL model, but method {self.method} does not use one; "
                "remove model or set method: nagl"
            )
        return self


#: Names of the components PolyzyMD adds itself (and their segment IDs), so no solute takes them.
_BUILT_IN_COMPONENT_NAMES = frozenset({"solvent", "ions", "SOL", "ION"})

#: Charge methods that compute charges from a molecule's graph.
COMPUTED_CHARGES = ("nagl", "am1bcc", "espaloma")


class _Solute(_ConfigModel):
    """The keys every entry of ``solutes:`` has.

    ``type`` selects the component class and is never inferred, so a typo is
    an error rather than a silently different component. ``name`` is the
    handle that restraints, ``modeling_exceptions``, provenance and analysis
    use; it defaults to the type. ``segid`` is the PSF segment ID the
    component's atoms carry; it defaults from the type.
    """

    name: str = Field(
        None,  # type: ignore[assignment]  # filled from the type below
        description="Handle for restraints, exceptions and analysis (default: the type)",
    )
    segid: str = Field(
        None,  # type: ignore[assignment]  # filled from the type below
        description="PSF segment ID, 1-4 characters (default from the type)",
    )
    force_field: list[str] = Field(
        ..., min_length=1, description="Force field files (.offxml) that parameterize it"
    )

    @model_validator(mode="after")
    def _defaults_from_type(self) -> Self:
        """Name and segment ID default from the type; check the segment ID fits a PSF."""
        from polyzymd.core.components import DEFAULT_SEGIDS, validate_segid

        type_: str = self.type  # type: ignore[attr-defined]  # each subclass defines it
        if self.name is None:
            object.__setattr__(self, "name", type_)
        if (
            not self.name.strip()
            or any(c in self.name for c in "/\\")
            or self.name.split() != [self.name]
        ):
            raise ValueError(
                f"name {self.name!r} must be non-empty, without '/', '\\' or spaces (it names "
                "run folders and GROMACS index groups)"
            )
        if self.name in _BUILT_IN_COMPONENT_NAMES:
            raise ValueError(
                f"name {self.name!r} is taken by a component PolyzyMD adds itself; choose another"
            )
        if self.segid is None:
            object.__setattr__(self, "segid", DEFAULT_SEGIDS[type_])
        validate_segid(self.segid)
        return self


class ProteinSolute(_Solute):
    """A protein (``type: protein``), read from a PDB file through residue templates.

    Its charges are the protein force field's library charges.

    Attributes:
        pdb_path: PDB file of the protein, protonated.
        description: Optional description.
        custom_substructures_path: Optional JSON file of residue templates for
            residues that OpenFF's template matcher does not know.
    """

    type: Literal["protein"]
    pdb_path: Path = Field(..., description="Path to the protein PDB file")
    description: str | None = Field(None, description="Optional description")
    custom_substructures_path: Path | None = Field(
        None,
        description=(
            "JSON file of extra residue templates, passed to OpenFF's "
            "Topology.from_pdb(_custom_substructures=...): "
            '{"RESNAME": {"SMARTS": ["ATOM", ...]}}. Relative to the config.'
        ),
    )
    charges: ChargeSpec = Field(
        ..., description="library: the protein force field's library charges"
    )

    @field_validator("custom_substructures_path")
    @classmethod
    def validate_custom_substructures(cls, v: Path | None) -> Path | None:
        """Check that the file exists and has the shape OpenFF reads."""
        if v is None:
            return v
        load_custom_substructures(v)
        return v

    @field_validator("pdb_path")
    @classmethod
    def validate_pdb_path(cls, v: Path) -> Path:
        """Validate that PDB path has correct extension."""
        if v.suffix.lower() != ".pdb":
            raise ValueError(f"Expected .pdb file, got {v.suffix}")
        return v

    @model_validator(mode="after")
    def _library_charges(self) -> Self:
        if self.charges.method != "library":
            raise ValueError(
                "a protein takes the library charges of its force field; use "
                "charges: {method: library}"
            )
        return self


class SmallMoleculeSolute(_Solute):
    """A small molecule placed as given (``type: small_molecule``), such as a docked substrate.

    Attributes:
        sdf_path: SDF file with the molecule's conformers.
        conformer_index: Which conformer to place (0-indexed).
        residue_name: 3-letter residue name.
        charges: How its partial charges are computed (nagl, am1bcc or espaloma).
    """

    type: Literal["small_molecule"]
    sdf_path: Path = Field(..., description="Path to the SDF with the docked conformers")
    conformer_index: int = Field(0, ge=0, description="Index of the conformer to place")
    residue_name: str = Field("LIG", max_length=3, description="3-letter residue name")
    charges: ChargeSpec = Field(..., description="nagl (with model), am1bcc or espaloma")

    @field_validator("sdf_path")
    @classmethod
    def validate_sdf_path(cls, v: Path) -> Path:
        """Validate that SDF path has correct extension."""
        if v.suffix.lower() != ".sdf":
            raise ValueError(f"Expected .sdf file, got {v.suffix}")
        return v

    @model_validator(mode="after")
    def _computed_charges(self) -> Self:
        if self.charges.method not in COMPUTED_CHARGES:
            raise ValueError(
                f"charges method {self.charges.method} is not built for a small molecule; "
                "use nagl, am1bcc or espaloma"
            )
        return self

    @property
    def charge_method(self) -> ChargeMethod:
        """The charge method, as the builders take it."""
        return ChargeMethod(self.charges.method)


# =============================================================================
# Polymer Packing Configuration
# =============================================================================


class PolymerPackingConfig(_ConfigModel):
    """Settings for packing polymers around the solute.

    Controls the PACKMOL behavior when packing polymers around the
    protein-ligand complex.

    Attributes:
        padding: Room reserved for the polymer chains around the solute, in
            nanometers.  It is added to ``solvent.box.padding`` when the
            periodic cell is computed (before packing, from the protein +
            substrate alone) and is the padding of the confinement sphere.
        tolerance: Minimum molecular spacing for PACKMOL in Angstrom.
        movebadrandom: When ``True``, pass the ``movebadrandom`` keyword to
            PACKMOL.  This places badly-packed molecules at random positions
            in the box rather than near well-packed neighbours, which
            improves convergence for dense or heterogeneous systems (many
            unique chain types).  Has no effect when only a single chain
            type is present.  Default is ``False`` (PACKMOL default
            behaviour is preserved).
        confine_to_sphere: When ``True`` (default), chains are packed inside a
            sphere centred on the solute with radius
            ``0.5 * |solute bbox diagonal| + padding``, in addition to the
            periodic brick.  Without it the chains spread into the corners of
            the brick and drift away from the protein.
        nloop: Maximum PACKMOL GENCAN loops per molecule type (default 200).

    Example:
        >>> PolymerPackingConfig(padding=2.5, movebadrandom=True)
    """

    padding: float = Field(
        2.0, gt=0.0, description="Room reserved for polymer chains around the solute (nm)"
    )
    tolerance: float = Field(2.0, gt=0.0, description="PACKMOL tolerance (Angstrom)")
    confine_to_sphere: bool = Field(
        True,
        description=(
            "Confine polymer chains to a sphere around the solute (radius = solute "
            "bounding-box circumradius + padding) while packing inside the final "
            "periodic brick. Set false to let chains fill the whole brick."
        ),
    )
    nloop: int = Field(200, ge=1, description="Maximum PACKMOL GENCAN loops per molecule type")
    movebadrandom: bool = Field(
        False,
        description=(
            "Pass the 'movebadrandom' keyword to PACKMOL. "
            "Improves convergence for heterogeneous polymer systems "
            "by placing badly-packed molecules at random box positions."
        ),
    )

    @model_validator(mode="before")
    @classmethod
    def reject_removed_packing_keys(cls, data: Any) -> Any:
        """Refuse the removed packing keys and name the keys that replace them."""
        if not isinstance(data, dict):
            return data
        if "exclude_solute_bbox" in data:
            raise ValueError(
                "'exclude_solute_bbox' was removed: chains are kept off the solute by the "
                "Packmol tolerance and near it by 'confine_to_sphere'. Delete the key."
            )
        if "box_vectors" in data:
            raise ValueError(
                "'box_vectors' was removed: the periodic cell is the protein + substrate "
                "bounding box plus 2 x ('polymers.packing.padding' + 'solvent.box.padding'). "
                "Delete the key and set those paddings."
            )
        return data


# =============================================================================
# Polymer Configuration
# =============================================================================


class MonomerSpec(_ConfigModel):
    """Specification for a single monomer type in a co-polymer.

    For dynamic generation, ``smiles`` is the repeat unit as it is in the chain,
    with ``*`` on the two atoms that bond to the neighbouring repeat units:
    ``[*]CC([*])(C)C(=O)OC`` for methyl methacrylate. The first ``*`` bonds
    toward the head of the chain, the second toward its tail.

    Attributes:
        label: Single character label for this monomer (e.g., "A", "B")
        probability: Probability of selecting this monomer (0-1)
        name: Optional full name (e.g., "SBMA", "EGPMA")
        smiles: Star-sited repeat-unit SMILES (required for dynamic generation)
        residue_name: 3-character PDB residue name (auto-generated if not provided)
    """

    label: str = Field(..., min_length=1, max_length=1, description="Monomer label")
    probability: float = Field(..., ge=0.0, le=1.0, description="Selection probability")
    name: str | None = Field(None, description="Full monomer name")
    smiles: str | None = Field(
        None,
        description=(
            "Repeat-unit SMILES with '*' on the two linking atoms (required for dynamic "
            "generation)"
        ),
    )
    residue_name: str | None = Field(
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
    DYNAMIC = "dynamic"  # Build chains with mBuild from star-sited SMILES


class EndGroupOption(_ConfigModel):
    """One cap a chain end can get, drawn with ``probability``.

    ``smiles`` has one ``*`` on the atom that bonds to the chain: ``[*][H]``
    caps with a hydrogen, ``[*]C`` with a methyl group, and
    ``[*]C=C(C)C(=O)OC`` with a whole methyl methacrylate unit that kept its
    double bond. ``label`` is the one character the chain's sequence shows for
    this cap. With ``residue_name`` the cap is a residue of its own; without
    it, its atoms join the residue of the repeat unit it caps.
    """

    label: str = Field(..., min_length=1, max_length=1, description="Cap label in the sequence")
    smiles: str = Field(..., min_length=1, description="Cap SMILES with one '*'")
    probability: float = Field(..., ge=0.0, le=1.0, description="Probability of this cap")
    name: str | None = Field(None, description="Cap name")
    residue_name: str | None = Field(
        None,
        max_length=3,
        description="3-char residue name; the cap is then its own residue",
    )


class PolymerEndGroups(_ConfigModel):
    """The caps of a dynamically generated chain: the options of each end.

    Each end draws one of its options per chain; the probabilities of each
    end sum to 1. A fixed cap is one option with probability 1. There is no
    default; the caps are part of the molecule the study simulates.
    """

    head: list[EndGroupOption] = Field(
        ..., min_length=1, description="Caps of the first repeat unit's head '*'"
    )
    tail: list[EndGroupOption] = Field(
        ..., min_length=1, description="Caps of the last repeat unit's tail '*'"
    )

    @field_validator("head", "tail", mode="before")
    @classmethod
    def _refuse_bare_smiles(cls, value: Any, info: Any) -> Any:
        """Refuse a cap given as a plain SMILES; each end is a list of options."""
        if isinstance(value, str):
            smiles = value if "*" in value else "[*][H]"
            raise ValueError(
                f"end_groups.{info.field_name} is a list of caps with probabilities, each "
                "SMILES with one '*', for example "
                f"{info.field_name}: [{{label: H, smiles: '{smiles}', probability: 1.0}}]"
            )
        return value


def _without_stereo(smiles: str) -> str:
    """Return ``smiles`` without its stereo marks (``@``, ``/``, ``\\``)."""
    return re.sub(r"[@/\\]", "", smiles)


def _count_stars(smiles: str) -> int:
    """Count the ``*`` atoms of a SMILES string (``*`` and ``[*]`` alike)."""
    return smiles.count("*")


class PolymerSolute(_Solute):
    """Polymer chains packed around the fixed solutes (``type: polymer``).

    Supports two generation modes:

    - "cached": Load pre-built polymer SDF files from sdf_directory
    - "dynamic": Build chains with mBuild from star-sited SMILES

    For dynamic mode, you must provide:

    - the star-sited repeat unit of each monomer in monomers[].smiles
    - the star-sited caps of each end in end_groups (head and tail), each a
      list of options with probabilities

    Attributes:
        generation_mode: "cached" for pre-built SDFs, "dynamic" for on-the-fly generation
        type_prefix: Prefix for polymer type in filenames (e.g., "SBMA-EGPMA")
        monomers: List of monomer specifications with probabilities (and SMILES for dynamic)
        length: Number of monomer units per polymer chain
        count: Number of polymer chains to add
        sdf_directory: Path to pre-built polymer SDF files (for cached mode)
        end_groups: Star-sited caps of each chain (required for dynamic mode)
        charges: preset (the charges in the cached SDFs) or, for dynamic
            generation, nagl, am1bcc or espaloma
        max_retries: Maximum retries for polymer generation (ring-piercing failures)
        cache_directory: Directory for caching generated polymers
        packing: Settings for packing polymers around the solute
        random_seed: Random seed for polymer sequence generation (for reproducibility)
    """

    type: Literal["polymer"]
    generation_mode: PolymerGenerationMode = Field(
        PolymerGenerationMode.CACHED,
        description="Polymer generation mode: 'cached' (pre-built SDFs) or 'dynamic' (generate from SMILES)",
    )
    type_prefix: str = Field(..., min_length=1, description="Polymer type prefix for filenames")
    monomers: list[MonomerSpec] = Field(..., min_length=1, description="Monomer specifications")
    length: int = Field(..., ge=1, description="Monomers per chain")
    count: int = Field(..., ge=1, description="Number of polymer chains")
    sdf_directory: Path | None = Field(
        None, description="Directory with pre-built polymer SDFs (for cached mode)"
    )
    end_groups: PolymerEndGroups | None = Field(
        None, description="Star-sited head and tail caps (required for dynamic mode)"
    )
    charges: ChargeSpec = Field(
        ..., description="preset (cached SDFs) or nagl, am1bcc, espaloma (dynamic generation)"
    )
    max_retries: int = Field(
        10,
        ge=1,
        description="Builds tried before a chain with a pierced ring or overlapping atoms is refused",
    )
    cache_directory: Path = Field(
        Path(".polymer_cache"), description="Cache directory for generated polymers"
    )
    packing: PolymerPackingConfig = Field(
        default_factory=PolymerPackingConfig,
        description="Polymer packing settings (padding, tolerance)",
    )
    random_seed: int | None = Field(
        None,
        description="Random seed for polymer sequence generation. If None, uses replicate number.",
    )

    @model_validator(mode="after")
    def validate_probabilities_sum_to_one(self) -> "PolymerSolute":
        """Ensure monomer probabilities sum to 1.0."""
        total = sum(m.probability for m in self.monomers)
        if abs(total - 1.0) > 1e-6:
            raise ValueError(f"Monomer probabilities must sum to 1.0, got {total}")
        return self

    @model_validator(mode="before")
    @classmethod
    def _refuse_reactions(cls, data: Any) -> Any:
        """Refuse ``reactions:``: chains are built from star-sited SMILES since 1.3.0rc5."""
        if isinstance(data, dict) and "reactions" in data:
            raise ValueError(
                "polymer key 'reactions' was removed in PolyzyMD 1.3.0rc5: dynamic chains are "
                "built from star-sited SMILES. Give each monomer's repeat unit with '*' on its "
                "two linking atoms ([*]CC([*])(C)C(=O)OC for methyl methacrylate) and the caps "
                "in end_groups, for example end_groups: {head: [{label: H, smiles: '[*][H]', "
                f"probability: 1.0}}], tail: [...]}}. See {REWRITE_GUIDE}"
            )
        return data

    @model_validator(mode="after")
    def validate_generation_mode_requirements(self) -> "PolymerSolute":
        """Validate that required fields are present for the selected generation mode."""
        if self.generation_mode == PolymerGenerationMode.DYNAMIC:
            # Dynamic mode requires SMILES for all monomers
            missing_smiles = [m.label for m in self.monomers if m.smiles is None]
            if missing_smiles:
                raise ValueError(
                    f"Dynamic generation mode requires 'smiles' for all monomers. "
                    f"Missing SMILES for monomers: {missing_smiles}"
                )
            for monomer in self.monomers:
                stars = _count_stars(monomer.smiles)
                if stars != 2:
                    raise ValueError(
                        f"monomer {monomer.label!r}: its smiles {monomer.smiles!r} has {stars} "
                        "'*' sites; a repeat unit needs exactly 2, on the atoms that bond to "
                        "the neighbouring units ([*]CC([*])(C)C(=O)OC for methyl "
                        "methacrylate). The notebook at "
                        "https://polyzymd.readthedocs.io/en/latest/tutorials/"
                        "star_sited_monomers.html helps you place them."
                    )
            if self.length < 3:
                raise ValueError(
                    f"length {self.length}: a dynamic chain's length counts its two ends (head "
                    "and tail caps) and its repeat units, so it is at least 3; length 5 is "
                    "2 ends + 3 repeat units."
                )
            if self.end_groups is None:
                raise ValueError(
                    "Dynamic generation mode requires 'end_groups' with the caps of the head and "
                    "tail, each a list of SMILES with one '*' and probabilities, for example "
                    "end_groups: {head: [{label: H, smiles: '[*][H]', probability: 1.0}], "
                    "tail: [{label: H, smiles: '[*][H]', probability: 1.0}]}. There is no "
                    "default cap."
                )
            for end in ("head", "tail"):
                options = getattr(self.end_groups, end)
                for option in options:
                    stars = _count_stars(option.smiles)
                    if stars != 1:
                        raise ValueError(
                            f"end_groups.{end} {option.label!r}: {option.smiles!r} has {stars} "
                            "'*' sites; a cap needs exactly 1, on the atom that bonds to the "
                            "chain ([*][H], [*]C)."
                        )
                labels = [option.label for option in options]
                repeated = sorted({label for label in labels if labels.count(label) > 1})
                if repeated:
                    raise ValueError(
                        f"end_groups.{end} options share label {repeated}; give each its own."
                    )
                total = sum(option.probability for option in options)
                if abs(total - 1.0) > 1e-6:
                    hint = (
                        "a single cap needs probability 1.0"
                        if len(options) == 1
                        else "the caps of one end are drawn from these probabilities"
                    )
                    raise ValueError(
                        f"end_groups.{end} probabilities must sum to 1.0, got {total:g}; {hint}."
                    )
            for key in ("label", "name"):
                values = [getattr(m, key) for m in self.monomers if getattr(m, key)]
                repeated = sorted({v for v in values if values.count(v) > 1})
                if repeated:
                    raise ValueError(
                        f"monomers share {key} {repeated}; give each monomer its own {key}."
                    )
            # Stereoisomers of one monomer (the letters that set tacticity) may share
            # a residue name; other monomers may not.
            constitution: dict[str, list[str]] = {}
            for m in self.monomers:
                if m.residue_name:
                    constitution.setdefault(m.residue_name, []).append(m.smiles)
            repeated = sorted(
                name
                for name, forms in constitution.items()
                if len(forms) > 1
                and (
                    len({_without_stereo(f) for f in forms}) > 1  # different monomers
                    or len(set(forms)) < len(forms)  # the same SMILES twice
                )
            )
            if repeated:
                raise ValueError(
                    f"monomers share residue_name {repeated}, so analyses could not tell them "
                    "apart. Give each monomer its own residue_name (only stereoisomers of one "
                    "monomer, which set the tacticity, may share one)."
                )

        elif self.generation_mode == PolymerGenerationMode.CACHED:
            if self.end_groups is not None and self.length < 3:
                raise ValueError(
                    f"length {self.length}: with end_groups the length counts the two ends "
                    "and the repeat units, so it is at least 3."
                )
            # Cached mode requires sdf_directory
            if self.sdf_directory is None:
                raise ValueError(
                    "Cached generation mode requires 'sdf_directory' with pre-built polymer SDFs"
                )
            if self.charges.method != "preset":
                raise ValueError(
                    "cached polymers carry the charges stored in their SDFs; use "
                    "charges: {method: preset}"
                )
        if (
            self.generation_mode == PolymerGenerationMode.DYNAMIC
            and self.charges.method not in COMPUTED_CHARGES
        ):
            raise ValueError(
                "dynamically generated polymers are charged as they are built; use charges "
                "method nagl, am1bcc or espaloma"
            )

        return self

    @property
    def enabled(self) -> bool:
        """Always True: a config without polymers has no polymer solute."""
        return True

    @property
    def charger(self) -> ChargeMethod:
        """The charger dynamic generation uses (not read for cached polymers)."""
        if self.charges.method in COMPUTED_CHARGES:
            return ChargeMethod(self.charges.method)
        return ChargeMethod.NAGL


def _solute_type(value: Any) -> str | None:
    """Return the ``type`` of a solutes entry, from YAML data or a model."""
    if isinstance(value, dict):
        return value.get("type")
    return getattr(value, "type", None)


#: One entry of ``solutes:``; ``type`` picks the class.
Solute = Annotated[
    Union[
        Annotated[ProteinSolute, Tag("protein")],
        Annotated[SmallMoleculeSolute, Tag("small_molecule")],
        Annotated[PolymerSolute, Tag("polymer")],
    ],
    Discriminator(
        _solute_type,
        custom_error_type="solute_type",
        custom_error_message="type must be protein, small_molecule or polymer",
    ),
]


# =============================================================================
# Solvent Configuration
# =============================================================================


class CoSolventSpec(_ConfigModel):
    """Specification for a co-solvent component.

    You must specify either ``mole_fraction`` or ``concentration``, not both.

    For co-solvents in the built-in library (dmso, dmf, urea, ethanol, etc.),
    you can omit the smiles and density fields - they will be looked up
    automatically.

    Example (library co-solvent with mole fraction):
        >>> CoSolventSpec(name="dmso", mole_fraction=0.10)

    Example (library co-solvent with concentration):
        >>> CoSolventSpec(name="urea", concentration=2.0)

    Example (custom co-solvent):
        >>> CoSolventSpec(
        ...     name="my_solvent",
        ...     smiles="CCOC(=O)C",
        ...     mole_fraction=0.05
        ... )
    """

    name: str = Field(..., description="Co-solvent identifier")
    smiles: str | None = Field(None, description="SMILES string (optional if in library)")

    # Specification method 1: Mole fraction
    mole_fraction: float | None = Field(
        None,
        gt=0.0,
        lt=1.0,
        description="Mole fraction (0-1), e.g., 0.10 for 10 mol%",
    )

    # Specification method 2: Molar concentration
    concentration: float | None = Field(None, gt=0.0, description="Molar concentration (mol/L)")

    # Specification method 3: number of molecules
    count: int | None = Field(None, ge=1, description="Number of molecules in the box")

    charges: ChargeSpec = Field(
        None,  # type: ignore[assignment]  # required; refused below with a fix that names preset
        description=(
            "preset (the charges PolyzyMD ships for molecules such as dmso), or nagl "
            "(with model), am1bcc or espaloma"
        ),
    )

    # Optional physical property for library metadata
    density: float | None = Field(
        None,
        gt=0.0,
        description="Optional density in g/mL for co-solvent metadata",
    )

    residue_name: str | None = Field(None, max_length=3, description="Residue name")

    @model_validator(mode="before")
    @classmethod
    def reject_volume_fraction(cls, data: Any) -> Any:
        """Reject removed volume-fraction and charge_method keys."""
        if isinstance(data, dict) and "volume_fraction" in data:
            raise ValueError(
                "Co-solvent 'volume_fraction' has been removed. Use 'mole_fraction' "
                "for solvent mole fractions or 'concentration' for molarity; values are "
                "not converted automatically."
            )
        if isinstance(data, dict) and "charge_method" in data:
            raise ValueError(
                "co-solvent 'charge_method' is now 'charges'. Releases up to 1.3.0rc4 gave a "
                "co-solvent PolyzyMD ships charges for (such as dmso) those shipped charges "
                "whatever charge_method said; to keep them, write charges: {method: preset}. "
                "Any other co-solvent: charges: {method: nagl, model: ...} (or am1bcc, "
                f"espaloma). See {REWRITE_GUIDE}"
            )
        return data

    @model_validator(mode="after")
    def validate_and_populate(self) -> "CoSolventSpec":
        """Validate specification and populate missing fields from library."""
        from polyzymd.data.cosolvent_library import get_cosolvent

        # Check that exactly one composition method is specified
        given = [
            key
            for key in ("mole_fraction", "concentration", "count")
            if getattr(self, key) is not None
        ]
        if len(given) != 1:
            raise ValueError(
                f"Co-solvent '{self.name}': give exactly one of mole_fraction, concentration "
                f"and count, not {given or 'none'}"
            )

        from polyzymd.data.solvent_molecules import is_bundled_solvent

        if self.charges is None:
            shipped = (
                "charges: {method: preset} (the charges PolyzyMD ships with it) or "
                if is_bundled_solvent(self.name)
                else ""
            )
            raise ValueError(
                f"Co-solvent '{self.name}' needs charges, for example {shipped}"
                "charges: {method: nagl, model: ...} (or am1bcc, espaloma)"
            )
        if self.charges.method == "library":
            raise ValueError(
                f"Co-solvent '{self.name}': charges method library is for proteins; use "
                "preset, nagl, am1bcc or espaloma"
            )
        if self.charges.method == "preset" and not is_bundled_solvent(self.name):
            raise ValueError(
                f"Co-solvent '{self.name}': charges preset reads the charges PolyzyMD ships, "
                "and it ships none for this molecule; use nagl, am1bcc or espaloma"
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

        # Preserve library density as metadata when available
        if self.density is None and library_data:
            object.__setattr__(self, "density", library_data.density)

        # Set default residue name: the library's, or three letters of the
        # name that no protein, nucleic, water or ion selection matches.
        if self.residue_name is None:
            if library_data:
                default = library_data.residue_name
            else:
                letters = [c for c in self.name.upper() if c.isalnum()] or ["X"]
                head = "".join(letters[:2]).ljust(2, "X")
                candidates = [head + c for c in letters[2:] + list("0123456789")]
                reserved = reserved_residue_names()
                default = next(c for c in candidates if c not in reserved)
            object.__setattr__(self, "residue_name", default)

        return self


@functools.cache
def reserved_residue_names() -> frozenset[str]:
    """Return the residue names that protein, nucleic, water or ion selections match."""
    # No MDAnalysis here: the sim-cuda environments that run simulations do not have it.
    from polyzymd.core.atom_groups import (
        ION_RESIDUE_NAMES,
        NUCLEIC_SELECTION_RESIDUE_NAMES,
        PROTEIN_SELECTION_RESIDUE_NAMES,
        WATER_RESIDUE_NAMES,
    )

    return frozenset(
        PROTEIN_SELECTION_RESIDUE_NAMES
        | NUCLEIC_SELECTION_RESIDUE_NAMES
        | {"SOL", "H2O", "OH2", "TIP", "T3P", "T4P", "T5P"}
        | WATER_RESIDUE_NAMES
        | ION_RESIDUE_NAMES
    )


class WaterConfig(_ConfigModel):
    """The water of the box.

    Attributes:
        model: Water model; only TIP3P is built.
    """

    model: WaterModel = Field(..., description="Water model (tip3p is the one built)")


class IonConfig(_ConfigModel):
    """Configuration for ions in the solvent.

    ``nacl_concentration`` is the NaCl salt concentration. When ``neutralize``
    is true, PolyzyMD adds the Na+ or Cl- ions that cancel the charge of the
    solute and co-solvents on top of the salt.

    ``ionic_strength`` instead sets the ionic strength of the whole box,
    counter-ions included: Na+ and Cl- split the solute charge between them
    (Machado and Pantano, J. Chem. Theory Comput. 2020, 16, 1367), so a buffer
    can be matched by its ionic strength.
    """

    neutralize: bool = Field(True, description="Neutralize system charge")
    nacl_concentration: float = Field(
        0.0,
        ge=0.0,
        description="NaCl salt concentration (mol/L), before neutralizing ions",
    )
    ionic_strength: float | None = Field(
        None,
        gt=0.0,
        description=(
            "Ionic strength (mol/L) of the box, counter-ions included; replaces "
            "nacl_concentration"
        ),
    )
    kcl_concentration: float = Field(0.0, ge=0.0, description="KCl conc. (mol/L); not built")
    mgcl2_concentration: float = Field(0.0, ge=0.0, description="MgCl2 conc. (mol/L); not built")

    @model_serializer(mode="wrap")
    def _serialize(self, handler: SerializerFunctionWrapHandler) -> dict[str, Any]:
        """Omit an unset ionic_strength, so configs without it keep their build hash."""
        data = handler(self)
        if self.ionic_strength is None:
            data.pop("ionic_strength", None)
        return data

    @model_validator(mode="after")
    def _one_salt_target(self) -> "IonConfig":
        """Refuse an ionic strength beside a salt concentration or without neutralizing."""
        if self.ionic_strength is None:
            return self
        if self.nacl_concentration:
            raise ValueError(
                "give either ionic_strength or nacl_concentration, not both: ionic_strength "
                "already sets the Na+ and Cl- counts"
            )
        if not self.neutralize:
            raise ValueError(
                "ionic_strength sets the ions of a neutral box; remove neutralize: false"
            )
        return self


class BoxConfig(_ConfigModel):
    """Configuration for the simulation box.

    Attributes:
        padding: Distance from solute to box edge in nm. The cell edge is the
            solute diameter plus 2 x padding.
        shape: Box geometry
        target_density: Target density in g/mL
        tolerance: Minimum molecular spacing for PACKMOL in Angstrom
    """

    padding: float = Field(1.2, gt=0.0, description="Box padding (nm)")
    shape: BoxShape = Field(BoxShape.RHOMBIC_DODECAHEDRON, description="Box shape")
    target_density: float = Field(1.0, gt=0.0, description="Target density (g/mL)")
    tolerance: float = Field(2.0, gt=0.0, description="PACKMOL tolerance (Angstrom)")


class SolventConfig(_ConfigModel):
    """Everything PolyzyMD packs around the solutes: water, co-solvents and ions.

    Water and co-solvents carry segment ID ``SOL`` (co-solvents are told apart
    by residue name) and ions ``ION``.

    Attributes:
        force_field: Force field files for water, ions and co-solvents. Sage
            carries the TIP3P water and the ion parameters.
        water: The water model.
        co_solvents: Co-solvents, each with its amount and charges.
        ions: Salt and neutralizing ions.
        box: Box shape and padding.
    """

    force_field: list[str] = Field(
        ..., min_length=1, description="Force field files for water, ions and co-solvents"
    )
    water: WaterConfig = Field(..., description="Water model")
    co_solvents: list[CoSolventSpec] = Field(default_factory=list, description="Co-solvents")
    ions: IonConfig = Field(default_factory=IonConfig, description="Ion settings")
    box: BoxConfig = Field(default_factory=BoxConfig, description="Box settings")

    @model_validator(mode="before")
    @classmethod
    def _primary_is_now_water(cls, data: Any) -> Any:
        if isinstance(data, dict) and "primary" in data:
            raise ValueError(
                "solvent.primary is now solvent.water: {model: tip3p}. " f"See {REWRITE_GUIDE}"
            )
        return data

    @model_validator(mode="after")
    def validate_mole_fractions(self) -> "SolventConfig":
        """Ensure co-solvent mole fractions leave room for water."""
        total = sum(cs.mole_fraction for cs in self.co_solvents if cs.mole_fraction is not None)
        if total >= 1.0:
            raise ValueError(f"Total co-solvent mole fraction must be < 1.0, got {total}")
        return self


# =============================================================================
# Restraint Configuration
# =============================================================================


class AtomSelectionConfig(_ConfigModel):
    """An atom of a restraint: a component, and a selection inside it.

    The selection is an MDTraj selection in which ``resid``, ``chain`` and
    ``pdbindex`` keep their MDAnalysis meaning (PDB residue number, chain
    letter, position in the built system counted from 1, ``bynum``; not the
    PDB serial, which skips a number at every TER). It is matched only among the atoms of
    ``component``, so ``resid 106`` cannot pick up a polymer residue 106.

    Attributes:
        component: Name of the solute the atom belongs to (a solute's name
            defaults to its type, so ``protein`` when the protein has no name).
        selection: Selection string (e.g., "resid 77 and name OG")
        description: Optional human-readable description
    """

    component: str = Field(..., description="Name of the solute the atom belongs to")
    selection: str = Field(..., description="Atom selection, e.g. 'resid 77 and name OG'")
    description: str | None = Field(None, description="Human-readable description")


class RestraintConfig(_ConfigModel):
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


#: The docs page listing every physics check and how to accept a failing one.
CHECKS_GUIDE = "https://polyzymd.readthedocs.io/en/latest/reference/physics_checks.html"


class ModelingException(_ConfigModel):
    """A deliberate departure from one of PolyzyMD's physics checks.

    A build runs every check; a failing check stops it unless the config
    accepts that failure here. The reason is not read by code: it is there so
    that a reader of a published config sees why a non-standard choice was
    made. The build provenance records each accepted failure with its reason
    and the components and atoms it covered.

    Attributes:
        check: The check's ID, such as ``geometry.buried_water``.
        component: The component the exception covers (a solute name, ``solvent``,
            ``ions`` or ``system``); every component the check runs on when left out.
        reason: Why the departure is deliberate (required).
    """

    check: str = Field(..., description=f"Check ID (see {CHECKS_GUIDE})")
    component: str | None = Field(None, description="Component it covers; all when left out")
    reason: str = Field(..., min_length=1, description="Why this departure is deliberate")

    @field_validator("check")
    @classmethod
    def _registered(cls, value: str) -> str:
        import polyzymd.audit.checks  # noqa: F401  # registers the checks
        from polyzymd.audit import registry

        try:
            registry.get(value)
        except KeyError as error:  # the message names the closest IDs
            raise ValueError(f"{error.args[0]}; every check: {CHECKS_GUIDE}") from None
        return value

    @field_validator("reason")
    @classmethod
    def _real_reason(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("give the reason this departure is deliberate")
        return value


# =============================================================================
# Thermodynamics Configuration
# =============================================================================


class ThermodynamicsConfig(_ConfigModel):
    """Thermodynamic conditions for the simulation.

    Attributes:
        temperature: System temperature in Kelvin
        pressure: System pressure in atmospheres (for NPT)
    """

    temperature: float = Field(..., gt=0.0, description="Temperature (K)")
    pressure: float = Field(1.0, gt=0.0, description="Pressure (atm)")


# =============================================================================
# Simulation Phase Configuration
# =============================================================================


class SimulationPhaseConfig(_ConfigModel):
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
        checkpoint_interval: Wall-time interval (seconds) between restart
            checkpoints for preemption resilience
    """

    ensemble: Ensemble = Field(..., description="Thermodynamic ensemble")
    duration: float = Field(..., gt=0.0, description="Duration (ns)")
    samples: int = Field(..., ge=1, description="Trajectory frames to save")
    time_step: float = Field(2.0, gt=0.0, description="Time step (fs)")
    thermostat: ThermostatType = Field(
        ThermostatType.LANGEVIN_MIDDLE, description="Thermostat type"
    )
    thermostat_timescale: float = Field(1.0, gt=0.0, description="Thermostat timescale (ps)")
    barostat: BarostatType | None = Field(None, description="Barostat type")
    barostat_frequency: int = Field(25, ge=1, description="Barostat update frequency")
    checkpoint_interval: float = Field(
        60.0,
        gt=0.0,
        description=(
            "Wall-time interval in seconds between restart checkpoints (default 60). "
            "Controls how frequently simulation state is saved for automatic "
            "restart on SLURM preemption or hard kill. Independent of "
            "trajectory/reporter output frequency. Must be positive so "
            "unexpected termination recovery remains enabled."
        ),
    )

    @model_validator(mode="before")
    @classmethod
    def reject_report_interval(cls, data: Any) -> Any:
        """Refuse the removed reporter interval and name the key that replaces it."""
        if isinstance(data, dict) and "report_interval" in data:
            raise ValueError(
                "'report_interval' was removed: the trajectory interval is the phase's steps "
                "divided by 'samples'. Delete 'report_interval' and set 'samples'."
            )
        return data

    @model_validator(mode="after")
    def validate_ensemble_barostat(self) -> "SimulationPhaseConfig":
        """Ensure NPT ensemble has a barostat configured."""
        if self.ensemble == Ensemble.NPT and self.barostat is None:
            self.barostat = BarostatType.MONTE_CARLO
        return self


# =============================================================================
# Equilibration Stage Configuration (Multi-Stage Equilibration)
# =============================================================================


class PositionRestraintConfig(_ConfigModel):
    """Configuration for positional restraints on an atom group.

    Position restraints apply a harmonic potential to keep atoms near their
    initial coordinates. This is commonly used during equilibration to
    prevent large structural changes while the system relaxes.

    Attributes:
        group: Predefined atom group name
        force_constant: Force constant in kJ/mol/nm^2 (4184.0 = 10 kcal/mol/A^2)
    """

    group: str = Field(
        ...,
        description=(
            "Atom group: protein_heavy, protein_backbone, protein_calpha, "
            "ligand_heavy, polymer_heavy, solvent, water_only, ions_only, cosolvents_only"
        ),
    )
    force_constant: float = Field(
        4184.0,  # 10 kcal/mol/A^2 in kJ/mol/nm^2
        gt=0.0,
        description="Force constant (kJ/mol/nm^2). Default 4184.0 = 10 kcal/mol/A^2",
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


class EquilibrationStageConfig(_ConfigModel):
    """Configuration for a single equilibration stage.

    Supports two temperature modes:

    1. Constant temperature: Set 'temperature' field.
    2. Temperature ramping (simulated annealing): Set 'temperature_start',
       'temperature_end', 'temperature_increment', and
       'temperature_interval_steps'. The duration is derived.

    Position restraints can be applied to hold specific atom groups in place
    during the stage.

    Attributes:
        name: Stage identifier (used in output paths)
        duration: Stage duration in nanoseconds (constant-temperature stages only)
        samples: Number of trajectory frames to save
        ensemble: Thermodynamic ensemble (NVT or NPT)
        temperature: Constant temperature in K (mutually exclusive with ramping)
        temperature_start: Starting temperature for ramping in K
        temperature_end: Ending temperature for ramping in K
        temperature_increment: Temperature increase per update in K
        temperature_interval_steps: MD steps between temperature updates
        position_restraints: List of position restraints for this stage
        time_step: Optional time step override in fs
        thermostat: Optional thermostat type override
        thermostat_timescale: Optional thermostat timescale override in ps
        barostat: Optional barostat type (for NPT ensemble)
        barostat_frequency: Optional barostat update frequency
    """

    name: str = Field(..., description="Stage identifier (used in output paths)")
    duration: float | None = Field(
        None,
        gt=0.0,
        description="Duration (ns). Required for constant-temperature stages; derived for ramps.",
    )
    samples: int = Field(100, ge=1, description="Trajectory frames to save")
    ensemble: Ensemble = Field(Ensemble.NVT, description="Thermodynamic ensemble")

    # Constant temperature mode
    temperature: float | None = Field(
        None,
        gt=0.0,
        description="Constant temperature (K). Mutually exclusive with temperature ramping.",
    )

    # Temperature ramping mode (simulated annealing)
    temperature_start: float | None = Field(
        None, gt=0.0, description="Starting temperature for ramping (K)"
    )
    temperature_end: float | None = Field(
        None, gt=0.0, description="Ending temperature for ramping (K)"
    )
    temperature_increment: float | None = Field(
        None,
        gt=0.0,
        description="Temperature increase per update (K)",
    )
    temperature_interval_steps: int | None = Field(
        None,
        ge=1,
        description="MD steps between temperature updates",
    )

    # Position restraints
    position_restraints: list[PositionRestraintConfig] = Field(
        default_factory=list, description="Position restraints for this stage"
    )

    # Optional overrides (inherit from parent/defaults if not specified)
    time_step: float | None = Field(None, gt=0.0, description="Time step (fs)")
    thermostat: ThermostatType | None = Field(None, description="Thermostat type")
    thermostat_timescale: float | None = Field(
        None, gt=0.0, description="Thermostat timescale (ps)"
    )
    barostat: BarostatType | None = Field(None, description="Barostat type (for NPT)")
    barostat_frequency: int | None = Field(None, ge=1, description="Barostat update frequency")

    @model_validator(mode="before")
    @classmethod
    def reject_removed_ramp_keys(cls, data: Any) -> Any:
        """Refuse the removed ramp keys and name ``temperature_interval_steps``."""
        if not isinstance(data, dict):
            return data
        if "temperature_ramp_rate" in data:
            raise ValueError(
                "'temperature_ramp_rate' is not supported; specify "
                "'temperature_increment' and 'temperature_interval_steps'"
            )
        if "temperature_interval" in data:
            raise ValueError(
                "'temperature_interval' was removed; give the MD steps between temperature "
                "updates as 'temperature_interval_steps' (1200 fs at a 2 fs time step is 600)"
            )
        ramping = any(data.get(key) is not None for key in ("temperature_start", "temperature_end"))
        if ramping and data.get("duration") is not None:
            raise ValueError(
                "Do not specify 'duration' for a temperature ramp; it is derived from "
                "'temperature_start', 'temperature_end', 'temperature_increment', "
                "and 'temperature_interval_steps'"
            )
        return data

    @model_validator(mode="after")
    def validate_temperature_mode(self) -> "EquilibrationStageConfig":
        """Ensure valid temperature specification."""
        if self.temperature is not None:
            if self.temperature_start is not None or self.temperature_end is not None:
                raise ValueError(
                    "Cannot specify both 'temperature' and temperature ramping "
                    "('temperature_start'/'temperature_end')"
                )
            if (
                self.temperature_increment is not None
                or self.temperature_interval_steps is not None
            ):
                raise ValueError(
                    "'temperature_increment' and 'temperature_interval_steps' are only "
                    "valid for temperature ramps"
                )
            if self.duration is None:
                raise ValueError("Constant-temperature stages require 'duration'")
            return self

        if self.temperature_start is None or self.temperature_end is None:
            raise ValueError(
                "Must specify either 'temperature' for constant temperature, "
                "or 'temperature_start', 'temperature_end', 'temperature_increment', "
                "and 'temperature_interval_steps' for ramping"
            )

        if self.temperature_increment is None or self.temperature_interval_steps is None:
            raise ValueError(
                "Temperature ramps require 'temperature_increment' and 'temperature_interval_steps'"
            )

        if self.temperature_start >= self.temperature_end:
            raise ValueError(
                f"temperature_start ({self.temperature_start}) must be < "
                f"temperature_end ({self.temperature_end})"
            )

        timestep_fs = self.time_step or 2.0
        total_steps = self.temperature_ramp_updates * self.temperature_interval_steps
        self.duration = total_steps * timestep_fs / 1e6

        return self

    @model_validator(mode="after")
    def validate_npt_barostat(self) -> "EquilibrationStageConfig":
        """Ensure NPT ensemble has a barostat configured."""
        if self.ensemble == Ensemble.NPT and self.barostat is None:
            self.barostat = BarostatType.MONTE_CARLO
        return self

    @model_serializer(mode="wrap")
    def serialize_stage(self, handler: SerializerFunctionWrapHandler) -> dict[str, Any]:
        """Omit internally derived ramp duration from saved configurations."""
        data = handler(self)
        if self.is_temperature_ramping:
            data.pop("duration", None)
        return data

    @property
    def is_temperature_ramping(self) -> bool:
        """Check if this stage uses temperature ramping."""
        return self.temperature_start is not None

    def get_start_temperature(self) -> float:
        """Get the starting temperature of this stage."""
        if self.is_temperature_ramping:
            assert self.temperature_start is not None
            return self.temperature_start
        assert self.temperature is not None
        return self.temperature

    def get_final_temperature(self) -> float:
        """Get the final temperature of this stage."""
        if self.is_temperature_ramping:
            assert self.temperature_end is not None
            return self.temperature_end
        assert self.temperature is not None
        return self.temperature

    @property
    def resolved_duration(self) -> float:
        """Return the validated stage duration in nanoseconds."""
        if self.duration is None:
            raise RuntimeError(f"Equilibration stage '{self.name}' has no resolved duration")
        return self.duration

    @property
    def temperature_ramp_updates(self) -> int:
        """Return the number of discrete updates needed to reach the endpoint."""
        if not self.is_temperature_ramping:
            return 0
        assert self.temperature_start is not None
        assert self.temperature_end is not None
        assert self.temperature_increment is not None
        exact_updates = (self.temperature_end - self.temperature_start) / self.temperature_increment
        nearest_updates = round(exact_updates)
        if math.isclose(exact_updates, nearest_updates, rel_tol=1e-12, abs_tol=1e-9):
            return max(1, nearest_updates)
        return max(1, math.ceil(exact_updates))

    def temperature_at_step(self, step: int, total_steps: int) -> float:
        """Return the scheduled ramp temperature at an integration step."""
        if not self.is_temperature_ramping:
            return self.get_start_temperature()
        if total_steps <= 0:
            raise ValueError("total_steps must be positive")
        assert self.temperature_start is not None
        assert self.temperature_end is not None
        assert self.temperature_increment is not None
        assert self.temperature_interval_steps is not None
        bounded_step = min(max(step, 0), total_steps)
        completed_updates = bounded_step // self.temperature_interval_steps
        return min(
            self.temperature_end,
            self.temperature_start + completed_updates * self.temperature_increment,
        )


class MinimizationConfig(_ConfigModel):
    """Energy-minimisation settings.

    Attributes:
        freeze_solute: Hold every protein and substrate atom fixed during
            minimisation so that only solvent and polymers relax. The prepared
            structure therefore enters equilibration with its coordinates
            unchanged (solute RMSD to the input structure is zero), and no
            packing artefact can be "accommodated" by protein deformation.
        max_iterations: Maximum minimiser iterations (0 = until convergence).
        tolerance: Energy tolerance in kJ/mol/nm.
    """

    model_config = ConfigDict(extra="forbid")

    freeze_solute: bool = Field(
        True, description="Freeze protein and substrate atoms during minimisation"
    )
    max_iterations: int = Field(1000, ge=0, description="Maximum minimiser iterations")
    tolerance: float = Field(10.0, gt=0, description="Energy tolerance in kJ/mol/nm")


class SimulationPhasesConfig(_ConfigModel):
    """Configuration for all simulation phases.

    Attributes:
        minimization: Energy-minimisation settings (runtime only, see below)
        equilibration_stages: Multi-stage equilibration protocol
        production: Production phase settings
    """

    model_config = ConfigDict(extra="forbid")

    # Runtime-only setting: excluded from model_dump() so that adding or changing it
    # does not alter the build-manifest config hash and invalidate pre-built bundles
    # (minimisation happens at run time, not at build time).
    minimization: MinimizationConfig = Field(
        default_factory=MinimizationConfig,
        exclude=True,
        description="Energy-minimisation settings (does not affect the built system)",
    )

    equilibration_stages: list[EquilibrationStageConfig] | None = Field(
        None, description="Multi-stage equilibration protocol with position restraints"
    )

    production: SimulationPhaseConfig = Field(..., description="Production settings")

    @model_validator(mode="after")
    def validate_equilibration_mode(self) -> "SimulationPhasesConfig":
        """Ensure staged equilibration is configured."""
        if self.equilibration_stages is None:
            raise ValueError(
                "simulation_phases.equilibration_stages is missing: list at least one stage"
            )

        if len(self.equilibration_stages) == 0:
            raise ValueError("'equilibration_stages' must contain at least one stage")

        return self

    @property
    def total_equilibration_duration(self) -> float:
        """Total equilibration duration in nanoseconds.

        Returns the sum of all stage durations.
        """
        return sum(stage.resolved_duration for stage in self.equilibration_stages or [])

    @property
    def total_equilibration_samples(self) -> int:
        """Total equilibration trajectory samples.

        Returns the sum of all stage samples.
        """
        return sum(stage.samples for stage in self.equilibration_stages or [])


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

    expanded = os.path.expandvars(str(path))
    expanded = os.path.expanduser(expanded)
    return Path(expanded)


#: Placeholders ``output.naming_template`` may use.
NAMING_PLACEHOLDERS = (
    "enzyme",
    "substrate",
    "polymer_type",
    "temperature",
    "replicate",
    "duration",
    "primary_solvent",
    "cosolvent_composition",
    "solvent_composition",
)


class OutputConfig(_ConfigModel):
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
        save_checkpoint: Not read; kept so the configs of existing runs load
        save_state_data: Not read; kept so the configs of existing runs load
        trajectory_format: Not read; the engine picks the format (OpenMM dcd, GROMACS xtc)

    Example YAML:
        output:
          projects_directory: /projects/user/polyzymd
          scratch_directory: /scratch/alpine/user/simulations
          naming_template: "{enzyme}_{substrate}_{temperature}K_run{replicate}"
    """

    model_config = ConfigDict(extra="forbid")

    UNREAD_FIELDS: ClassVar[dict[str, str]] = {
        "save_checkpoint": (
            "output.save_checkpoint is not read: both engines always write checkpoints. "
            "Leave it out of new configs."
        ),
        "save_state_data": (
            "output.save_state_data is not read: both engines always write energies. "
            "Leave it out of new configs."
        ),
        "trajectory_format": (
            "output.trajectory_format is not read: OpenMM writes dcd and GROMACS writes xtc. "
            "Leave it out of new configs."
        ),
    }

    # Directory configuration
    projects_directory: Path = Field(
        Path("."),
        description="Base directory for scripts and logs (typically in projects/long-term storage)",
    )
    scratch_directory: Path | None = Field(
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
    save_checkpoint: bool = Field(True, description="Not read: checkpoints are always written")
    save_state_data: bool = Field(True, description="Not read: energies are always written")
    trajectory_format: str = Field(
        "dcd", description="Not read: OpenMM writes dcd and GROMACS writes xtc trajectories"
    )

    @field_validator("naming_template")
    @classmethod
    def validate_naming_template(cls, v: str) -> str:
        """Refuse a placeholder the run-folder name does not fill."""
        try:
            fields = [field for _, field, _, _ in Formatter().parse(v) if field is not None]
        except ValueError as error:
            raise ValueError(f"naming_template {v!r} is not a valid template: {error}") from None
        unknown = sorted({field for field in fields if field not in NAMING_PLACEHOLDERS})
        if unknown:
            raise ValueError(
                f"naming_template {v!r} has unknown placeholders {unknown}; use only "
                + ", ".join("{" + name + "}" for name in NAMING_PLACEHOLDERS)
            )
        return v

    @field_validator("projects_directory", "scratch_directory", mode="before")
    @classmethod
    def expand_env_vars_in_paths(cls, v: str | Path | None) -> Path | None:
        """Expand environment variables and ~ in path fields.

        Supports $USER, ${HOME}, ~/path, etc.
        """
        if v is None:
            return None
        return expand_path(Path(v))

    @property
    def effective_scratch_directory(self) -> Path:
        """Get the effective scratch directory (falls back to projects if not set).

        The path is used as written: ``..`` parts are not collapsed as text,
        because behind a symlinked folder the operating system resolves them
        to another place. Printed paths are shortened where they are printed.
        """
        return self.scratch_directory if self.scratch_directory else self.projects_directory

    def format_directory_name(
        self,
        enzyme: str,
        substrate: str,
        polymer_type: str,
        temperature: float,
        replicate: int | str,
        duration: float = 0.0,
        primary_solvent: str = "water_tip3p",
        cosolvent_composition: str = "none",
        solvent_composition: str | None = None,
        **kwargs: Any,
    ) -> str:
        """Format the directory name using the template.

        Parameters
        ----------
        enzyme : str
            Enzyme name.
        substrate : str
            Substrate name.
        polymer_type : str
            Polymer type, or ``"none"``.
        temperature : float
            Temperature in K.
        replicate : int or str
            Replicate number or glob token.
        duration : float, optional
            Production duration in ns, by default 0.0. Written as whole ns
            from 1 ns up (``100``) and in full below 1 ns (``0.005``), so a
            short run does not become ``0``.
        primary_solvent : str, optional
            Primary solvent token, by default ``"water_tip3p"``.
        cosolvent_composition : str, optional
            Co-solvent composition token, by default ``"none"``.
        solvent_composition : str or None, optional
            Full solvent composition token. When not provided, the primary solvent
            and co-solvent tokens are combined.
        **kwargs : Any
            Additional template variables.

        Returns
        -------
        str
            Formatted directory name.
        """
        if solvent_composition is None:
            solvent_composition = primary_solvent
            if cosolvent_composition != "none":
                solvent_composition = f"{primary_solvent}_{cosolvent_composition}"

        return self.naming_template.format(
            enzyme=enzyme,
            substrate=substrate,
            polymer_type=polymer_type,
            temperature=int(temperature),
            replicate=replicate,
            duration=int(duration) if duration >= 1 else f"{duration:g}",
            primary_solvent=primary_solvent,
            cosolvent_composition=cosolvent_composition,
            solvent_composition=solvent_composition,
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


#: OpenMM platform names, keyed by the lowercase spelling a config may use
#: (``simulation/platform.py`` maps them; HIP must be written as is).
OPENMM_PLATFORMS = {"cuda": "CUDA", "cpu": "CPU", "opencl": "OpenCL", "reference": "Reference"}


class ForceFieldConfig(_ConfigModel):
    """The force fields and NAGL model a build applies, read from the solutes and solvent.

    Not a config section: :attr:`SimulationConfig.force_field` derives it, so the
    builders keep one place to read them from.

    Attributes:
        protein: Force field of the protein solute.
        small_molecule: Force field of every other molecule (rc5 builds them with one).
        nagl_model: NAGL model of every component charged with nagl.
    """

    protein: str = Field("ff14sb_off_impropers_0.0.4.offxml", description="Protein force field")
    small_molecule: str = Field("openff-2.3.0.offxml", description="Small molecule force field")
    nagl_model: str = Field(DEFAULT_NAGL_MODEL, description="NAGL model")


# =============================================================================
# Engine Configuration
# =============================================================================


class OpenMMEngineConfig(_ConfigModel):
    """OpenMM-specific engine settings.

    These settings control OpenMM platform selection and device configuration.
    Only relevant when ``engine`` is ``"openmm"``.

    Attributes:
        platform: OpenMM platform to use for computation.
        device_index: GPU device index (platform-specific).
        precision: Floating-point precision mode.
    """

    platform: str = Field("CUDA", description="OpenMM compute platform")
    device_index: str | None = Field(None, description="GPU device index")
    precision: str = Field("mixed", description="Floating-point precision")


class GromacsEngineConfig(_ConfigModel):
    """GROMACS-specific engine settings.

    These settings control how GROMACS binaries are located and invoked,
    and how SLURM resources are allocated for GROMACS jobs.

    GROMACS users configure parallelism via ``ntmpi`` (MPI ranks) and
    ``ntomp`` (OpenMP threads per rank). These map directly to the
    ``gmx mdrun -ntmpi`` and ``-ntomp`` flags. The corresponding SLURM
    resources (``ntasks``, ``cpus_per_task``, ``memory``) are set to
    match.

    Set ``gpu: true`` when running a CUDA-enabled GROMACS build on a
    GPU node. When ``False`` (default), the SLURM script omits GPU
    directives entirely.

    Attributes:
        gmx_binary: Path or name of the GROMACS binary.
            Resolved via config > $GMX_BIN > PATH discovery if None.
        mdrun_flags: Additional flags passed to ``gmx mdrun`` (all stages).
        grompp_flags: Additional flags passed to ``gmx grompp``.
        module_load: Module load command for HPC (e.g. ``"module load gromacs/2024"``)
        ntmpi: Number of MPI ranks for ``gmx mdrun -ntmpi``.
            Also sets SLURM ``--ntasks``.
        ntomp: Number of OpenMP threads per rank for ``gmx mdrun -ntomp``.
            Also sets SLURM ``--cpus-per-task``.
        gpu: Request a GPU via SLURM. When False, the ``--gres=gpu``
            directive is omitted entirely.
        gpus: Number of GPUs to request when ``gpu`` is True. Ignored
            when ``gpu`` is False.
        memory: SLURM ``--mem`` allocation for GROMACS jobs.
    """

    gmx_binary: str | None = Field(None, description="GROMACS binary path or name")
    analysis_topology: str | None = Field(
        None,
        description=(
            "File name of the run's GROMACS topology (.top) in the run folder. Analyses "
            "read it when MDAnalysis cannot read prod.tpr, and freeze deposits it with the "
            "files it includes. Default: <prefix>.top, the file PolyzyMD writes."
        ),
    )
    mdrun_flags: str = Field("", description="Extra flags for gmx mdrun (all stages)")
    grompp_flags: str = Field(
        "",
        description="Extra flags for gmx grompp, such as '-maxwarn 1' to accept a warning "
        "you have read; by default every grompp warning stops the run",
    )
    mdrun_flags_equilibration: str | None = Field(
        None,
        description=(
            "Override mdrun_flags for equilibration stages only. "
            "When None, falls back to mdrun_flags."
        ),
    )
    mdrun_flags_production: str | None = Field(
        None,
        description=(
            "Override mdrun_flags for production only. When None, falls back to mdrun_flags."
        ),
    )
    command_prefix: str | None = Field(
        None,
        description=(
            "Prefix prepended to all GROMACS commands. "
            "Use for container wrappers, e.g. "
            "'singularity exec --rocm --bind $PWD /path/to/gromacs.sif'"
        ),
    )
    mpi_launcher_flags: str = Field(
        "",
        description=(
            "Extra flags passed to the MPI launcher (mpirun). "
            "Only used when the binary is a real-MPI build. "
            "E.g. '-genv I_MPI_FABRICS shm:tcp' for Intel MPI."
        ),
    )
    module_load: str | None = Field(None, description="HPC module load command")
    env_exports: dict[str, str] = Field(
        default_factory=dict,
        description=("Environment variables exported before GROMACS commands (key=value pairs)"),
    )
    setup_commands: list[str] = Field(
        default_factory=list,
        description="Shell commands run after module_load and before GROMACS commands",
    )
    ntmpi: int = Field(1, ge=1, description="MPI ranks (-ntmpi); sets SLURM --ntasks")
    slurm_ntasks: int | None = Field(
        default=None,
        ge=1,
        description=(
            "Override SLURM --ntasks independently of GROMACS -ntmpi. "
            "Advanced option for multi-node MPI+GPU workflows where scheduler "
            "tasks differ from GROMACS thread-MPI ranks. When None, SLURM ntasks "
            "is set from ntmpi (default behavior)."
        ),
    )
    ntomp: int = Field(
        8,
        ge=1,
        description="OpenMP threads per rank (-ntomp); sets SLURM --cpus-per-task",
    )
    gpu: bool = Field(
        False,
        description="Request GPU via SLURM (set True for CUDA-enabled GROMACS)",
    )
    gpus: int = Field(
        1,
        ge=1,
        description="Number of GPUs to request via SLURM when gpu is True",
    )
    memory: str = Field("16G", description="SLURM --mem allocation for GROMACS jobs")

    @model_validator(mode="after")
    def _warn_gpu_ntmpi(self) -> Self:
        """Warn when GPU mode is combined with multiple thread-MPI ranks.

        Returns
        -------
        Self
            The validated config instance
        """
        if self.gpu and self.ntmpi > 1:
            LOGGER.warning(
                "GPU GROMACS uses thread-MPI (gmx mdrun -ntmpi), not real MPI ranks. "
                "With ntmpi=%d and gpus=%d, ensure enough GPUs are available or "
                "accept GPU sharing.",
                self.ntmpi,
                self.gpus,
            )
        return self

    @model_validator(mode="after")
    def _warn_mdrun_flags_conflict(self) -> Self:
        """Warn when mdrun_flags conflict with explicit ntmpi or ntomp fields.

        Returns
        -------
        Self
            The validated config instance
        """
        if not self.mdrun_flags:
            return self

        try:
            tokens = shlex.split(self.mdrun_flags)
        except ValueError:
            return self

        flag_map: dict[str, str] = {}
        for i, tok in enumerate(tokens):
            if tok in ("-ntmpi", "-ntomp") and i + 1 < len(tokens):
                flag_map[tok] = tokens[i + 1]

        if "-ntmpi" in flag_map:
            try:
                flag_val = int(flag_map["-ntmpi"])
                if flag_val != self.ntmpi:
                    LOGGER.warning(
                        "mdrun_flags contains '-ntmpi %s' but config field ntmpi=%d. "
                        "The flag string will override at runtime.",
                        flag_map["-ntmpi"],
                        self.ntmpi,
                    )
            except ValueError:
                pass

        if "-ntomp" in flag_map:
            try:
                flag_val = int(flag_map["-ntomp"])
                if flag_val != self.ntomp:
                    LOGGER.warning(
                        "mdrun_flags contains '-ntomp %s' but config field ntomp=%d. "
                        "The flag string will override at runtime.",
                        flag_map["-ntomp"],
                        self.ntomp,
                    )
            except ValueError:
                pass

        return self

    @model_validator(mode="after")
    def _refuse_noconfout(self) -> Self:
        """Refuse ``-noconfout``: without it mdrun writes the final ``.gro`` and checkpoint.

        PolyzyMD finds a finished equilibration stage by its ``.gro`` and a
        finished production by its last checkpoint, so a run with
        ``-noconfout`` would never count as finished.
        """
        for name in ("mdrun_flags", "mdrun_flags_equilibration", "mdrun_flags_production"):
            if "-noconfout" in (getattr(self, name) or "").split():
                raise ValueError(
                    f"gromacs.{name} contains -noconfout. PolyzyMD needs the final .gro "
                    "and checkpoint that mdrun writes to know a stage finished; remove it."
                )
        return self


# =============================================================================
# Main Simulation Configuration
# =============================================================================


class SimulationConfig(_ConfigModel):
    """Complete simulation configuration.

    This is the top-level configuration model that contains all settings
    for a PolyzyMD simulation.

    The system is a ``solutes:`` list (everything added to the box: a
    protein, a docked small molecule, polymer chains) and a ``solvent:`` block
    (water, co-solvents, ions, box). Each solute names its force field and
    charge source. :attr:`enzyme`, :attr:`substrate`, :attr:`polymers` and
    :attr:`force_field` are read-only views of them for the builders.

    Attributes:
        name: Simulation name/identifier
        description: Optional description
        solutes: The components added to the box, each with its ``type``
        solvent: Water, co-solvents, ions and box
        restraints: List of restraint configurations
        thermodynamics: Temperature/pressure settings
        simulation_phases: Equilibration and production settings
        output: Output file settings

    Example:
        >>> config = SimulationConfig.from_yaml("simulation.yaml")
        >>> [solute.name for solute in config.solutes]
        ['LipA', 'substrate', 'polymer']
    """

    #: Sections of 1.3.0rc4 and earlier configs, and where they went.
    REMOVED_SECTIONS: ClassVar[dict[str, str]] = {
        "enzyme": "a solutes entry with type: protein",
        "substrate": "a solutes entry with type: small_molecule",
        "polymers": "a solutes entry with type: polymer",
        "force_field": "force_field lists on each solute and on solvent",
    }

    name: str = Field(..., description="Simulation identifier")
    description: str | None = Field(None, description="Simulation description")

    solutes: list[Solute] = Field(..., min_length=1, description="Components added to the box")
    solvent: SolventConfig = Field(..., description="Water, co-solvents, ions and box")
    restraints: list[RestraintConfig] = Field(
        default_factory=list, description="Restraint configurations"
    )
    modeling_exceptions: list[ModelingException] = Field(
        default_factory=list,
        description=f"Deliberate departures from the physics checks, each with a reason ({CHECKS_GUIDE})",
    )
    thermodynamics: ThermodynamicsConfig = Field(..., description="Thermodynamic conditions")
    simulation_phases: SimulationPhasesConfig = Field(..., description="Phase settings")
    output: OutputConfig = Field(default_factory=OutputConfig, description="Output settings")
    engine: Literal["openmm", "gromacs"] = Field(..., description="Simulation engine to use")
    openmm: OpenMMEngineConfig = Field(
        default_factory=OpenMMEngineConfig, description="OpenMM engine settings"
    )
    gromacs: GromacsEngineConfig = Field(
        default_factory=GromacsEngineConfig, description="GROMACS engine settings"
    )

    @model_validator(mode="before")
    @classmethod
    def _refuse_removed_sections(cls, data: Any) -> Any:
        """Name where each section of an older config went, instead of a bare unknown-key error."""
        if not isinstance(data, dict):
            return data
        removed = [key for key in cls.REMOVED_SECTIONS if key in data]
        if removed:
            moves = "; ".join(f"{key} -> {cls.REMOVED_SECTIONS[key]}" for key in removed)
            raise ValueError(
                f"this config uses the layout of PolyzyMD 1.3.0rc4 and earlier ({moves}). "
                "Rewrite it from the template (polyzymd study add-condition LABEL --new "
                f"writes one) or see {REWRITE_GUIDE}"
            )
        return data

    @model_validator(mode="after")
    def _solutes_fit_this_release(self) -> Self:
        """Check names, segment IDs and the combinations this release builds."""
        from polyzymd.core.components import DEFAULT_SEGIDS

        types = [solute.type for solute in self.solutes]
        for type_ in sorted(set(types)):
            if types.count(type_) > 1:
                raise ValueError(
                    f"this release builds one solute of each type, and the config has "
                    f"{types.count(type_)} of type {type_}"
                )
        if "protein" not in types:
            raise ValueError(
                "this release builds systems around a protein; add a type: protein solute"
            )
        names = [solute.name for solute in self.solutes]
        repeated = sorted({name for name in names if names.count(name) > 1})
        if repeated:
            raise ValueError(
                f"solute name {', '.join(repeated)} is used twice; give each solute its own "
                "name (restraints and analysis refer to solutes by name)"
            )
        reserved = {
            DEFAULT_SEGIDS["solvent"]: "water and co-solvents",
            DEFAULT_SEGIDS["ions"]: "ions",
        }
        segids = [solute.segid for solute in self.solutes]
        for solute in self.solutes:
            if solute.segid in reserved:
                raise ValueError(
                    f"solute {solute.name!r} has segid {solute.segid}, which is the segment "
                    f"of the {reserved[solute.segid]}; choose another"
                )
            if segids.count(solute.segid) > 1:
                raise ValueError(
                    f"segid {solute.segid} is used by more than one solute; each solute "
                    "needs its own segment ID"
                )
        for solute in self.solutes:
            if len(solute.force_field) != 1:
                raise ValueError(
                    f"solute {solute.name!r}: this release applies one force field file per "
                    f"solute, and it lists {len(solute.force_field)}"
                )
        if len(self.solvent.force_field) != 1:
            raise ValueError("solvent.force_field: this release applies one force field file")
        shared = self.solvent.force_field
        for solute in self.solutes:
            if solute.type != "protein" and solute.force_field != shared:
                raise ValueError(
                    f"solute {solute.name!r} lists force_field {solute.force_field}, but "
                    f"solvent lists {shared}; this release parameterizes every molecule but the "
                    "protein with one force field, so they must be the same"
                )
        models = {charges.model for charges in self._charge_specs() if charges.method == "nagl"}
        if len(models) > 1:
            raise ValueError(
                f"components name different NAGL models ({', '.join(sorted(models))}); this "
                "release charges every nagl component with one model"
            )
        for restraint in self.restraints:
            for atom in (restraint.atom1, restraint.atom2):
                if atom.component not in names:
                    raise ValueError(
                        f"restraint {restraint.name!r} names component {atom.component!r}, "
                        f"which is not a solute; solutes are {', '.join(names)}"
                    )
        from polyzymd.audit import registry

        type_of = {solute.name: solute.type for solute in self.solutes}
        type_of.update({"solvent": "solvent", "ions": "ions"})
        for exception in self.modeling_exceptions:
            component = exception.component
            if component is None:
                continue
            spec = registry.get(exception.check)
            # A check reports on the components of the types it applies to; a
            # system check may also report on "system".
            allowed = {n for n, t in type_of.items() if spec.applies_to(t)}
            if spec.scope == "system":
                allowed.add("system")
            if component not in allowed:
                raise ValueError(
                    f"modeling_exceptions: {exception.check} cannot report on {component!r}; "
                    f"it reports on {', '.join(sorted(allowed)) or 'no component of this config'}. "
                    "Name one of those, or leave component out to cover all of them."
                )
        return self

    @property
    def nagl_model(self) -> str | None:
        """The NAGL model the components charged with nagl name, or None when none is."""
        models = sorted({c.model for c in self._charge_specs() if c.method == "nagl" and c.model})
        return models[0] if models else None

    def _charge_specs(self) -> list[ChargeSpec]:
        return [solute.charges for solute in self.solutes] + [
            cosolvent.charges for cosolvent in self.solvent.co_solvents
        ]

    def solute(self, type_: str) -> Solute | None:
        """Return the solute of ``type_``, or None (this release builds one of each)."""
        return next((solute for solute in self.solutes if solute.type == type_), None)

    @property
    def enzyme(self) -> ProteinSolute:
        """The protein solute."""
        protein = self.solute("protein")
        assert isinstance(protein, ProteinSolute)
        return protein

    @property
    def substrate(self) -> SmallMoleculeSolute | None:
        """The small-molecule solute, if any."""
        return self.solute("small_molecule")  # type: ignore[return-value]

    @property
    def polymers(self) -> PolymerSolute | None:
        """The polymer solute, if any."""
        return self.solute("polymer")  # type: ignore[return-value]

    @property
    def force_field(self) -> ForceFieldConfig:
        """The force fields and NAGL model the builders apply, read from the solutes."""
        models = [c.model for c in self._charge_specs() if c.method == "nagl" and c.model]
        return ForceFieldConfig(
            protein=self.enzyme.force_field[0],
            small_molecule=self.solvent.force_field[0],
            nagl_model=models[0] if models else DEFAULT_NAGL_MODEL,
        )

    @model_validator(mode="after")
    def validate_barostat_for_engine(self) -> "SimulationConfig":
        """Refuse barostat MCA with engine openmm, which runs only the isotropic barostat."""
        if self.engine == "openmm":
            self.require_engine_barostats("openmm")
        return self

    def require_buildable(self, engine: str | None = None, *, inputs: bool = True) -> None:
        """Raise ``ValueError`` when a new system built from this config would be wrong.

        Loading a config does not run these checks, so the configs of existing
        runs still load for status and analysis. The commands that build a new
        system (``validate``, ``build``, ``run``, ``submit``) call this. It
        also checks the input files (:func:`polyzymd.config.validation.require_inputs`).

        Parameters
        ----------
        engine : str or None
            Engine that will run the system. Default: the config's engine.
        inputs : bool
            Check the input files too. ``submit`` passes False: it runs only
            replicates that ``build`` has already built.
        """
        engine = engine or self.engine
        model = self.solvent.water.model
        if model != WaterModel.TIP3P:
            raise ValueError(
                f"PolyzyMD cannot build water model '{model.value}' correctly: the force fields "
                "it loads carry TIP3P water parameters; use model: tip3p for new systems."
            )
        ions = self.solvent.ions
        for key in ("kcl_concentration", "mgcl2_concentration"):
            if getattr(ions, key):
                raise ValueError(
                    f"{key}: the build adds only Na+ and Cl- ions. For new systems, remove "
                    f"{key} and use nacl_concentration."
                )
        for cs in self.solvent.co_solvents:
            if cs.residue_name.upper() in reserved_residue_names():
                raise ValueError(
                    f"Co-solvent '{cs.name}' has residue_name '{cs.residue_name}', an amino-acid, "
                    "nucleic-acid, water or ion name, so selections such as 'protein' would "
                    "include it. Set residue_name to another 3-letter name, such as the "
                    "molecule's PDB chemical-component code, or leave it out."
                )
        self._require_buildable_settings(engine)
        self.require_engine_restraints(engine)
        self.require_engine_barostats(engine)
        if inputs:
            from polyzymd.config.validation import require_inputs

            require_inputs(self)

    def _require_buildable_settings(self, engine: str) -> None:
        """Raise ``ValueError`` for settings a new build would not run as written."""
        phases = self.simulation_phases
        stages = phases.equilibration_stages or []
        named = [(f"equilibration stage {stage.name!r}", stage) for stage in stages]
        named.append(("production", phases.production))
        for what, phase in named:
            if phase.ensemble == Ensemble.NVE:
                raise ValueError(
                    f"{what}: ensemble NVE is not run; both engines keep the thermostat on."
                    "\nfix: use ensemble NVT or NPT."
                )
            if engine == "openmm" and phase.thermostat in (
                ThermostatType.NOSE_HOOVER,
                ThermostatType.ANDERSEN,
            ):
                raise ValueError(
                    f"{what}: thermostat {phase.thermostat.value} runs only with engine gromacs; "
                    "the OpenMM engine would run LangevinMiddle instead."
                    "\nfix: use thermostat LangevinMiddle or Langevin, or set engine: gromacs."
                )
        for what, phase in named:
            if isinstance(phase, EquilibrationStageConfig):
                time_step = phase.time_step or 2.0
                steps = int(round(phase.resolved_duration * 1e6 / time_step))
            else:
                time_step = phase.time_step
                steps = int(phase.duration * 1e6 / time_step)
            if steps < 1:
                raise ValueError(
                    f"{what}: duration {phase.duration:g} ns is shorter than one "
                    f"{time_step:g} fs time step.\nfix: give a longer duration."
                )
            if phase.samples > steps:
                raise ValueError(
                    f"{what}: samples {phase.samples} is more than its {steps} steps "
                    f"({phase.duration:g} ns / {time_step:g} fs), and at most one frame is saved "
                    f"per step.\nfix: set samples to {steps} or fewer, or give a longer duration."
                )
        has_component = {
            "ligand_heavy": self.substrate is not None,
            "polymer_heavy": bool(self.polymers and self.polymers.enabled),
        }
        for stage in stages:
            for posres in stage.position_restraints:
                if not has_component.get(posres.group, True):
                    part = "substrate" if posres.group == "ligand_heavy" else "enabled polymers"
                    raise ValueError(
                        f"equilibration stage {stage.name!r} restrains {posres.group}, but the "
                        f"config has no {part}, so nothing would be restrained."
                        f"\nfix: remove that position restraint, or add the {part}."
                    )
        fields = {f for _, f, _, _ in Formatter().parse(self.output.naming_template) if f}
        if "replicate" not in fields:
            raise ValueError(
                "output.naming_template has no {replicate}, so every replicate would write to "
                "the same folder.\nfix: add {replicate} to naming_template."
            )
        folder = self.format_run_directory_name(1)
        if not folder.strip() or folder in (".", "..") or "/" in folder or "\\" in folder:
            raise ValueError(
                f"The run folder name {folder!r} (output.naming_template filled with the enzyme, "
                "substrate and polymer names) is not a single folder name."
                "\nfix: remove '/' and '\\' from those names and from naming_template."
            )
        for stage in stages:
            if not stage.name.strip() or "/" in stage.name or "\\" in stage.name:
                raise ValueError(
                    f"Equilibration stage name {stage.name!r} is used in a folder name."
                    "\nfix: give each stage a non-empty name without '/' or '\\'."
                )
        from polyzymd.config.loader import UNSET_VARIABLE

        paths = {
            f"{self.enzyme.name}.pdb_path": self.enzyme.pdb_path,
            f"{self.enzyme.name}.custom_substructures_path": (
                self.enzyme.custom_substructures_path
            ),
            "output.projects_directory": self.output.projects_directory,
            "output.scratch_directory": self.output.scratch_directory,
        }
        if self.substrate is not None:
            paths[f"{self.substrate.name}.sdf_path"] = self.substrate.sdf_path
        if self.polymers is not None:
            paths[f"{self.polymers.name}.sdf_directory"] = self.polymers.sdf_directory
            paths[f"{self.polymers.name}.cache_directory"] = self.polymers.cache_directory
        for key, path in paths.items():
            unset = UNSET_VARIABLE.search(str(path)) if path is not None else None
            if unset:
                raise ValueError(
                    f"{key} {str(path)!r} names the environment variable ${unset.group(1)}, "
                    f"which is not set.\nfix: set {unset.group(1)} before running PolyzyMD, or "
                    "write the path out."
                )
        if engine == "openmm":
            platform = self.openmm.platform
            if platform != "HIP" and platform.lower() not in OPENMM_PLATFORMS:
                raise ValueError(
                    f"openmm.platform {platform!r} is not an OpenMM platform."
                    f"\nfix: use one of {', '.join(OPENMM_PLATFORMS.values())} or HIP."
                )
            if self.openmm.precision not in ("single", "mixed", "double"):
                raise ValueError(
                    f"openmm.precision {self.openmm.precision!r} is not an OpenMM precision."
                    "\nfix: use single, mixed or double."
                )
        if engine == "gromacs":
            if not re.fullmatch(r"\d+[KMGTkmgt]?", self.gromacs.memory):
                raise ValueError(
                    f"gromacs.memory {self.gromacs.memory!r} is not a SLURM --mem value."
                    "\nfix: give a number with an optional K, M, G or T unit, such as 16G."
                )
            from polyzymd.engines.gromacs.slurm import GromacsSlurmScriptGenerator

            prefix = self.gromacs.command_prefix
            pattern = GromacsSlurmScriptGenerator._COMMAND_PREFIX_PATTERN
            if prefix and not pattern.match(prefix):
                raise ValueError(
                    f"gromacs.command_prefix {prefix!r} has characters the job script refuses "
                    "($, quotes, ;, |, & and others).\nfix: write the command out without them."
                )

    def require_engine_restraints(self, engine: str) -> None:
        """Raise ``ValueError`` when ``engine`` would drop a restraint this config names.

        Only the OpenMM engine applies distance restraints (``restraints:``).
        The GROMACS export writes position restraints for the protein, ligand
        and polymer groups only. :meth:`require_buildable` calls this with the
        engine a command builds for, so ``build --format gromacs`` and
        ``run --engine gromacs`` on an OpenMM config are refused too.
        """
        if engine != "gromacs":
            return
        enabled = [r.name for r in self.restraints if r.enabled]
        if enabled:
            raise ValueError(
                f"Distance restraints ({', '.join(enabled)}) are applied only by the OpenMM "
                "engine; GROMACS runs would drop them. To keep them, set engine: openmm. "
                "To run on GROMACS, remove them or set enabled: false. Position restraints "
                "in equilibration stages work on both engines."
            )
        from polyzymd.exporters.gromacs import PositionRestraintGenerator

        mapping = PositionRestraintGenerator.GROUP_MAPPING
        force_constants: dict[str, set[float]] = {}
        for stage in self.simulation_phases.equilibration_stages or []:
            for posres in stage.position_restraints:
                force_constants.setdefault(posres.group, set()).add(posres.force_constant)
        unknown = sorted(group for group in force_constants if group not in mapping)
        if unknown:
            raise ValueError(
                f"GROMACS cannot apply position restraints on {', '.join(unknown)}: it restrains "
                f"only {', '.join(mapping)}. Remove those groups from position_restraints, "
                "or set engine: openmm."
            )
        # The GROMACS export writes one restraint block per group, with one
        # force constant, and the protein groups share the POSRES_PROTEIN define.
        varying = sorted(group for group, values in force_constants.items() if len(values) > 1)
        if varying:
            raise ValueError(
                f"On GROMACS a position restraint group has one force constant in every stage, "
                f"but {', '.join(varying)} has several. Use the same force_constant in "
                "every stage."
            )
        protein = sorted(group for group in force_constants if mapping[group][0] == "protein")
        if len(protein) > 1:
            raise ValueError(
                f"On GROMACS the protein groups share one restraint switch, so {', '.join(protein)} "
                "would all be on in every stage that names one of them. Use the same protein "
                "group in every stage."
            )

    def require_engine_barostats(self, engine: str) -> None:
        """Raise ``ValueError`` when ``engine`` cannot run a barostat this config names.

        The OpenMM engine runs only the isotropic Monte Carlo barostat, so
        barostat MCA (anisotropic) needs GROMACS. Commands that run OpenMM on
        a config call this, so ``run --engine openmm`` on a GROMACS config
        with MCA is refused too.
        """
        phases = self.simulation_phases
        barostats = [stage.barostat for stage in phases.equilibration_stages or []]
        barostats.append(phases.production.barostat)
        if engine == "openmm" and BarostatType.MONTE_CARLO_ANISOTROPIC in barostats:
            raise ValueError(
                "barostat MCA (anisotropic) runs only with engine gromacs. The OpenMM engine "
                "runs only the isotropic Monte Carlo barostat: set engine: gromacs, or use "
                "barostat MC."
            )

    @classmethod
    def from_yaml(cls, path: str | Path) -> "SimulationConfig":
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

    def _primary_solvent_token(self) -> str:
        """Format the primary solvent naming token.

        Returns
        -------
        str
            Primary solvent token for naming templates.
        """
        return f"water_{_format_safe_token(self.solvent.water.model.value)}"

    def _cosolvent_composition_token(self) -> str:
        """Format the co-solvent composition naming token.

        Returns
        -------
        str
            Co-solvent token or ``none`` when no co-solvents are configured.
        """
        if not self.solvent.co_solvents:
            return "none"

        tokens: list[tuple[str, str]] = []
        for cosolvent in self.solvent.co_solvents:
            name = _format_safe_token(cosolvent.name)
            if cosolvent.mole_fraction is not None:
                percent = _format_decimal_token(cosolvent.mole_fraction * 100)
                token = f"{name}_{percent}molpct"
            elif cosolvent.concentration is not None:
                concentration = _format_decimal_token(cosolvent.concentration)
                token = f"{name}_{concentration}M"
            elif cosolvent.count is not None:
                token = f"{name}_{cosolvent.count}mol"
            else:
                # Validation guarantees one composition mode, but keep errors explicit
                raise ValueError(f"Co-solvent '{cosolvent.name}' has no composition value")
            tokens.append((name, token))

        return "_".join(token for _, token in sorted(tokens, key=lambda item: item[0]))

    def _run_directory_template_values(self, replicate: int | str = 1) -> dict[str, Any]:
        """Build centralized values for run-directory naming templates.

        Parameters
        ----------
        replicate : int or str, optional
            Replicate token to insert into the template, by default 1.

        Returns
        -------
        dict of str to Any
            Template values shared by directory naming and replicate discovery.
        """
        polymer_type = "none"
        if self.polymers and self.polymers.enabled:
            probs = {m.label: m.probability for m in self.polymers.monomers}
            sorted_labels = sorted(probs.keys())
            composition = "_".join(f"{lbl}{probs[lbl] * 100:.0f}" for lbl in sorted_labels)
            polymer_type = f"{self.polymers.type_prefix}_{composition}"

        substrate_name = self.substrate.name if self.substrate else "apo"
        primary_solvent = self._primary_solvent_token()
        cosolvent_composition = self._cosolvent_composition_token()
        solvent_composition = primary_solvent
        if cosolvent_composition != "none":
            solvent_composition = f"{primary_solvent}_{cosolvent_composition}"

        return {
            "enzyme": self.enzyme.name,
            "substrate": substrate_name.replace("-", ""),
            "polymer_type": polymer_type,
            "temperature": int(self.thermodynamics.temperature),
            "replicate": replicate,
            "duration": self.simulation_phases.production.duration,
            "primary_solvent": primary_solvent,
            "cosolvent_composition": cosolvent_composition,
            "solvent_composition": solvent_composition,
        }

    def format_run_directory_name(self, replicate: int | str = 1) -> str:
        """Format the run directory name for a replicate.

        Parameters
        ----------
        replicate : int or str, optional
            Replicate number or glob token, by default 1.

        Returns
        -------
        str
            Formatted run directory name.
        """
        return self.output.format_directory_name(**self._run_directory_template_values(replicate))

    def get_working_directory(self, replicate: int = 1) -> Path:
        """Get the working directory path for a given replicate (in scratch).

        This returns the path where simulation output (trajectories, checkpoints)
        will be written, which is in the scratch directory.

        Args:
            replicate: Replicate number

        Returns:
            Path to the working directory (in scratch)
        """
        scratch = self.output.effective_scratch_directory
        path = scratch / self.format_run_directory_name(replicate)
        old_path = self._zero_duration_run_directory(replicate)
        if not path.exists() and old_path is not None and old_path.is_dir():
            return old_path
        return path

    def _zero_duration_run_directory(self, replicate: int | str) -> Path | None:
        """The run directory as named before durations below 1 ns kept their value.

        Earlier versions wrote a production shorter than 1 ns as ``0``
        (``..._0ns_...``). Returns that path, or None when the duration is
        1 ns or longer.
        """
        values = self._run_directory_template_values(replicate)
        if values["duration"] >= 1:
            return None
        values["duration"] = 0
        return self.output.effective_scratch_directory / self.output.format_directory_name(**values)

    def discover_replicate_dirs(self) -> list[tuple[int, Path]]:
        """Auto-detect all replicate directories on disk.

        Builds a glob pattern from the naming template with
        ``replicate="*"`` and scans the effective scratch directory. For a
        production shorter than 1 ns it also finds folders with the earlier
        ``0`` duration name; a folder with the current name wins.

        Returns:
            Sorted list of ``(replicate_number, directory_path)`` tuples.
        """
        fields = {
            field for _, field, _, _ in Formatter().parse(self.output.naming_template) if field
        }
        if "replicate" not in fields:
            raise ValueError(
                "Cannot discover replicate directories because output.naming_template does not "
                "include the {replicate} placeholder."
            )

        glob_patterns = [self.format_run_directory_name(replicate="*")]
        old_path = self._zero_duration_run_directory("*")
        if old_path is not None:
            glob_patterns.insert(0, old_path.name)

        scratch = self.output.effective_scratch_directory
        found: dict[int, Path] = {}
        for glob_pattern in glob_patterns:
            regex_pattern = (
                "^" + re.escape(glob_pattern).replace(r"\*", r"(?P<replicate>\d+)") + "$"
            )
            replicate_regex = re.compile(regex_pattern)
            for path in scratch.glob(glob_pattern):
                if not path.is_dir():
                    continue
                match = replicate_regex.fullmatch(path.name)
                if match:
                    found[int(match.group("replicate"))] = path

        return sorted(found.items())

    def _format_run_directory_name(self, replicate: int | str) -> str:
        """Format the run directory name for a replicate.

        Args:
            replicate: Replicate number

        Returns:
            Formatted directory name
        """
        return self.format_run_directory_name(replicate)
