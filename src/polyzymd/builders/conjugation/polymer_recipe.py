"""Polymer recipe models and Polymerist generation boundary."""

from __future__ import annotations

import random
from pathlib import Path

from pydantic import AliasChoices, BaseModel, ConfigDict, Field, field_validator, model_validator

POC_SBMA_SMILES = (
    "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])[N+]"
    "(C([H])([H])[H])(C([H])([H])[H])C([H])([H])C([H])([H])"
    "C([H])([H])S(=O)(=O)[O-])C([H])([H])[H]"
)
POC_EGPMA_SMILES = (
    "[H]C([H])=C(C(=O)OC([H])([H])C([H])([H])Oc1c([H])c([H])c([H])c([H])c1[H])C([H])([H])[H]"
)
POC_NHS_SMILES = "CC(=C)C(=O)ON1C(=O)CCC1=O"

DEFAULT_PROBABILITY_TOLERANCE = 1.0e-6


class PolymerMonomerRecipe(BaseModel):
    """Specification for one monomer in a stochastic polymer recipe."""

    label: str = Field(..., min_length=1, max_length=1)
    name: str = Field(..., min_length=1)
    residue_name: str = Field(..., min_length=3, max_length=3)
    smiles: str = Field(..., min_length=1)
    probability: float = Field(..., gt=0.0)

    @field_validator("label")
    @classmethod
    def normalize_label(cls, value: str) -> str:
        """Normalize a one-character sequence label."""
        normalized = value.strip().upper()
        if not normalized.isalnum():
            raise ValueError("Monomer labels must be alphanumeric")
        return normalized

    @field_validator("name")
    @classmethod
    def normalize_name(cls, value: str) -> str:
        """Normalize a monomer name for Polymerist fragment lookup."""
        normalized = value.strip()
        if not normalized:
            raise ValueError("Monomer names cannot be blank")
        return normalized

    @field_validator("residue_name")
    @classmethod
    def normalize_residue_name(cls, value: str) -> str:
        """Normalize and validate a PDB-safe residue name."""
        normalized = value.strip().upper()
        if len(normalized) != 3 or not normalized.isalnum():
            raise ValueError("Monomer residue names must be three alphanumeric PDB characters")
        return normalized


class PolymerRecipe(BaseModel):
    """Stochastic multi-residue polymer recipe for conjugation construction."""

    model_config = ConfigDict(populate_by_name=True)

    name: str = Field("polymer", min_length=1)
    monomers: tuple[PolymerMonomerRecipe, ...] = Field(..., min_length=1)
    length: int = Field(
        ...,
        ge=2,
        validation_alias=AliasChoices("length", "degree_of_polymerization"),
    )
    seed: int | None = None
    reactive_monomer_label: str | None = Field(
        None,
        validation_alias=AliasChoices(
            "reactive_monomer_label", "forced_reactive_monomer_label", "forced_label"
        ),
    )
    reactive_monomer_index: int | None = Field(
        None,
        ge=0,
        validation_alias=AliasChoices(
            "reactive_monomer_index", "forced_reactive_monomer_index", "forced_index"
        ),
    )
    probability_tolerance: float = Field(DEFAULT_PROBABILITY_TOLERANCE, gt=0.0)

    @field_validator("reactive_monomer_label")
    @classmethod
    def normalize_reactive_label(cls, value: str | None) -> str | None:
        """Normalize an optional forced reactive monomer label."""
        if value is None:
            return None
        normalized = value.strip().upper()
        if len(normalized) != 1 or not normalized.isalnum():
            raise ValueError("Reactive monomer label must be one alphanumeric character")
        return normalized

    @model_validator(mode="after")
    def validate_recipe(self) -> PolymerRecipe:
        """Validate probability, label, and forced reactive residue constraints."""
        labels = [monomer.label for monomer in self.monomers]
        if len(set(labels)) != len(labels):
            raise ValueError("Monomer labels must be unique")

        names = [monomer.name for monomer in self.monomers]
        if len(set(names)) != len(names):
            raise ValueError("Monomer names must be unique")

        total_probability = sum(monomer.probability for monomer in self.monomers)
        if abs(total_probability - 1.0) > self.probability_tolerance:
            raise ValueError(f"Monomer probabilities must sum to 1.0, got {total_probability}")

        if self.reactive_monomer_index is not None and self.reactive_monomer_index >= self.length:
            raise ValueError("Reactive monomer index must be within the polymer length")

        if self.reactive_monomer_label is not None and self.reactive_monomer_label not in labels:
            raise ValueError("Reactive monomer label must match a declared monomer label")

        return self

    @property
    def effective_reactive_index(self) -> int | None:
        """Return the configured reactive index or the centered forced-label index.

        Returns
        -------
        int or None
            Zero-based residue index used for forced reactive monomer placement.
        """
        if self.reactive_monomer_index is not None:
            return self.reactive_monomer_index
        if self.reactive_monomer_label is not None:
            return self.length // 2
        return None

    def monomer_by_label(self, label: str) -> PolymerMonomerRecipe:
        """Return a monomer specification by sequence label.

        Parameters
        ----------
        label : str
            One-character monomer sequence label.

        Returns
        -------
        PolymerMonomerRecipe
            Matching monomer recipe.

        Raises
        ------
        KeyError
            If no monomer has the requested label.
        """
        normalized = label.strip().upper()
        for monomer in self.monomers:
            if monomer.label == normalized:
                return monomer
        raise KeyError(f"Unknown monomer label: {label}")

    def generate_sequence(self, seed: int | None = None) -> str:
        """Generate a deterministic weighted monomer sequence.

        Parameters
        ----------
        seed : int or None, optional
            Override for the recipe seed, by default ``None``.

        Returns
        -------
        str
            Sequence string using monomer labels.
        """
        rng = random.Random(self.seed if seed is None else seed)
        labels = [monomer.label for monomer in self.monomers]
        weights = [monomer.probability for monomer in self.monomers]
        sequence = list(rng.choices(labels, weights=weights, k=self.length))

        reactive_index = self.effective_reactive_index
        if reactive_index is not None and self.reactive_monomer_label is not None:
            sequence[reactive_index] = self.reactive_monomer_label

        return "".join(sequence)

    def to_monomer_smiles(self) -> dict[str, str]:
        """Return the monomer-name to SMILES mapping expected by fragment generation."""
        return {monomer.name: monomer.smiles for monomer in self.monomers}

    def to_sequence_monomer_names(self) -> dict[str, str]:
        """Return the sequence-label to monomer-name mapping expected by Polymerist."""
        return {monomer.label: monomer.name for monomer in self.monomers}

    def to_polymerist_residue_names(self) -> dict[str, str]:
        """Return monomer-name to residue-name mapping for Polymerist PDB output."""
        return {monomer.name: monomer.residue_name for monomer in self.monomers}


class PolymeristGenerationSmokeResult(BaseModel):
    """Summary from the optional Polymerist recipe-generation boundary."""

    recipe_name: str
    sequence: str
    cache_directory: Path
    monomer_group_path: Path
    pdb_path: Path | None = None
    object_type: str
    atom_count: int | None = None


def sbma_egpma_nhs_recipe(
    *,
    length: int = 10,
    seed: int | None = 42,
    reactive_monomer_index: int | None = None,
) -> PolymerRecipe:
    """Build the SBMA/EGPMA/NHS recipe used by the conjugation POC.

    Parameters
    ----------
    length : int, optional
        Degree of polymerization, by default 10.
    seed : int or None, optional
        Random seed for deterministic sequence generation, by default 42.
    reactive_monomer_index : int or None, optional
        Zero-based NHS residue index. ``None`` centers NHS, by default ``None``.

    Returns
    -------
    PolymerRecipe
        Validated stochastic polymer recipe.
    """
    return PolymerRecipe(
        name="SBMA-EGPMA-NHS",
        monomers=(
            PolymerMonomerRecipe(
                label="A",
                name="SBMA",
                residue_name="SBM",
                smiles=POC_SBMA_SMILES,
                probability=0.945,
            ),
            PolymerMonomerRecipe(
                label="B",
                name="EGPMA",
                residue_name="EGP",
                smiles=POC_EGPMA_SMILES,
                probability=0.045,
            ),
            PolymerMonomerRecipe(
                label="C",
                name="NHS",
                residue_name="NHS",
                smiles=POC_NHS_SMILES,
                probability=0.01,
            ),
        ),
        length=length,
        seed=seed,
        reactive_monomer_label="C",
        reactive_monomer_index=reactive_monomer_index,
    )


def generate_polymerist_smoke_polymer(
    recipe: PolymerRecipe,
    cache_directory: Path | str,
    *,
    force_regenerate: bool = False,
    max_retries: int = 3,
) -> PolymeristGenerationSmokeResult:
    """Generate a small Polymerist-backed polymer artifact from a recipe.

    This boundary intentionally stops at Polymerist fragment and PDB generation.
    It does not perform placement, relaxation, charging, Interchange export, or
    simulation setup.

    Parameters
    ----------
    recipe : PolymerRecipe
        Validated polymer recipe containing monomer SMILES and probabilities.
    cache_directory : pathlib.Path or str
        Directory for Polymerist fragment and polymer PDB cache files.
    force_regenerate : bool, optional
        Rebuild cached Polymerist fragments when possible, by default ``False``.
    max_retries : int, optional
        Maximum Polymerist structure-building attempts, by default 3.

    Returns
    -------
    PolymeristGenerationSmokeResult
        Summary with generated sequence, cache paths, and object metadata.
    """
    from polyzymd.builders.conjugation.polymerist_compat import ensure_polymerist_py312_compat
    from polyzymd.data.reactions import get_atrp_reaction_paths

    ensure_polymerist_py312_compat()

    from polyzymd.builders.fragment_generator import FragmentGenerator
    from polyzymd.builders.polymer_generator import PolymerGenerator

    cache_path = Path(cache_directory)
    cache_path.mkdir(parents=True, exist_ok=True)
    reaction_paths = get_atrp_reaction_paths()

    fragment_generator = FragmentGenerator(
        initiation_rxn_path=reaction_paths["initiation"],
        polymerization_rxn_path=reaction_paths["polymerization"],
        termination_rxn_path=reaction_paths["termination"],
        cache_directory=cache_path,
    )
    monomer_group = fragment_generator.load_or_generate(
        monomer_smiles=recipe.to_monomer_smiles(),
        type_prefix=recipe.name,
        force_regenerate=force_regenerate,
    )

    sequence = recipe.generate_sequence()
    polymer_generator = PolymerGenerator(
        monomer_group=monomer_group,
        cache_directory=cache_path,
        max_retries=max_retries,
    )
    polymer_object, pdb_path = polymer_generator._build_polymer_structure(
        sequence=sequence,
        monomer_names=recipe.to_sequence_monomer_names(),
        residue_names=recipe.to_polymerist_residue_names(),
    )
    atom_count = _get_polymerist_atom_count(polymer_object)

    return PolymeristGenerationSmokeResult(
        recipe_name=recipe.name,
        sequence=sequence,
        cache_directory=cache_path,
        monomer_group_path=fragment_generator.get_cache_path(recipe.name),
        pdb_path=pdb_path,
        object_type=type(polymer_object).__name__,
        atom_count=atom_count,
    )


def _get_polymerist_atom_count(polymer_object: object) -> int | None:
    """Return a best-effort atom count from a Polymerist or mBuild object."""
    for attribute in ("n_particles", "n_atoms"):
        value = getattr(polymer_object, attribute, None)
        if isinstance(value, int):
            return value
    particles = getattr(polymer_object, "particles", None)
    if callable(particles):
        try:
            return sum(1 for _particle in particles())
        except TypeError:
            return None
    return None
