"""Polymer recipe models and Polymerist generation boundary."""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any

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
    fixed_sequence: str | None = Field(
        None,
        validation_alias=AliasChoices("fixed_sequence", "sequence", "deterministic_sequence"),
        description="Optional exact monomer-label sequence overriding stochastic generation",
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

    @field_validator("fixed_sequence")
    @classmethod
    def normalize_fixed_sequence(cls, value: str | None) -> str | None:
        """Normalize an optional fixed monomer-label sequence."""
        if value is None:
            return None
        normalized = "".join(value.split()).upper()
        if not normalized:
            raise ValueError("Fixed polymer sequence cannot be blank")
        if not normalized.isalnum():
            raise ValueError("Fixed polymer sequence labels must be alphanumeric")
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

        if self.fixed_sequence is not None:
            if len(self.fixed_sequence) != self.length:
                raise ValueError(
                    "Fixed polymer sequence length must match the configured polymer length"
                )
            unknown_labels = sorted(set(self.fixed_sequence) - set(labels))
            if unknown_labels:
                raise ValueError(
                    "Fixed polymer sequence labels must match declared monomer labels: "
                    f"{', '.join(unknown_labels)}"
                )
            if (
                self.reactive_monomer_index is not None
                and self.reactive_monomer_label is not None
                and self.fixed_sequence[self.reactive_monomer_index] != self.reactive_monomer_label
            ):
                raise ValueError(
                    "Fixed polymer sequence must contain the reactive monomer label at the "
                    "configured reactive index"
                )

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
        if self.fixed_sequence is not None:
            return self.fixed_sequence

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

    model_config = ConfigDict(arbitrary_types_allowed=True)

    recipe_name: str
    sequence: str
    cache_directory: Path
    monomer_group_path: Path
    pdb_path: Path | None = None
    sdf_path: Path | None = None
    charged_sdf_path: Path | None = None
    object_type: str
    atom_count: int | None = None
    rdkit_mol: Any | None = Field(default=None, exclude=True)


def sbma_egpma_nhs_recipe(
    *,
    length: int = 10,
    seed: int | None = 42,
    reactive_monomer_index: int | None = None,
    fixed_sequence: str | None = None,
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
    fixed_sequence : str or None, optional
        Exact monomer-label sequence overriding stochastic generation, by default
        ``None``.

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
        fixed_sequence=fixed_sequence,
    )


def sbma_nhs_egpma_acb_recipe() -> PolymerRecipe:
    """Build the deterministic v1 SBMA:NHS:EGPMA recipe.

    Returns
    -------
    PolymerRecipe
        Three-monomer recipe whose fixed sequence ``ACB`` maps to
        SBMA:NHS:EGPMA with the NHS monomer centered for Lys linkage.
    """
    return sbma_egpma_nhs_recipe(
        length=3,
        seed=None,
        reactive_monomer_index=1,
        fixed_sequence="ACB",
    )


def generate_polymerist_smoke_polymer(
    recipe: PolymerRecipe,
    cache_directory: Path | str,
    *,
    force_regenerate: bool = False,
    max_retries: int = 3,
    energy_minimize: bool = True,
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
    energy_minimize : bool, optional
        Whether Polymerist should run its built-in minimization during structure
        generation, by default ``True``.

    Returns
    -------
    PolymeristGenerationSmokeResult
        Summary with generated sequence, cache paths, and object metadata.
    """
    from polyzymd.builders.fragment_generator import FragmentGenerator
    from polyzymd.builders.polymer_generator import PolymerGenerator
    from polyzymd.data.reactions import get_atrp_reaction_paths

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
        energy_minimize=energy_minimize,
    )
    rdkit_mol = None
    sdf_path = None
    charged_sdf_path = None
    pdb_file = Path(pdb_path) if pdb_path is not None else None
    if pdb_file is not None and pdb_file.exists():
        sdf_path = pdb_file.with_suffix(".sdf")
        rdkit_mol = _polymerist_to_pdb_aligned_rdkit_mol(polymer_object, pdb_file)
        _write_rdkit_sdf_sidecar(rdkit_mol, sdf_path)
    generate_polymer = getattr(polymer_generator, "generate_polymer", None)
    make_filename = getattr(polymer_generator, "_make_polymer_filename", None)
    if callable(generate_polymer) and callable(make_filename):
        generate_polymer(
            sequence=sequence,
            monomer_names=recipe.to_sequence_monomer_names(),
            residue_names=recipe.to_polymerist_residue_names(),
            force_regenerate=force_regenerate,
        )
        charged_filename = make_filename(
            sequence,
            recipe.to_sequence_monomer_names(),
            charged=True,
        )
        expected_charged_sdf_path = cache_path / f"{charged_filename}.sdf"
        if expected_charged_sdf_path.exists():
            charged_sdf_path = expected_charged_sdf_path
    atom_count = _get_rdkit_atom_count(rdkit_mol)
    if atom_count is None:
        atom_count = _get_polymerist_atom_count(polymer_object)

    return PolymeristGenerationSmokeResult(
        recipe_name=recipe.name,
        sequence=sequence,
        cache_directory=cache_path,
        monomer_group_path=fragment_generator.get_cache_path(recipe.name),
        pdb_path=pdb_path,
        sdf_path=sdf_path,
        charged_sdf_path=charged_sdf_path,
        object_type=type(polymer_object).__name__,
        atom_count=atom_count,
        rdkit_mol=rdkit_mol,
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


def _get_rdkit_atom_count(rdkit_mol: object | None) -> int | None:
    """Return a best-effort atom count from an RDKit-like molecule."""
    if rdkit_mol is None:
        return None
    get_num_atoms = getattr(rdkit_mol, "GetNumAtoms", None)
    if not callable(get_num_atoms):
        return None
    try:
        atom_count = get_num_atoms()
    except (RuntimeError, TypeError, ValueError):
        return None
    if isinstance(atom_count, int):
        return atom_count
    return None


def _polymerist_to_pdb_aligned_rdkit_mol(
    polymer_object: object,
    pdb_path: Path,
    *,
    required: bool = True,
) -> Any | None:
    """Convert a Polymerist object to a PDB-aligned RDKit sidecar molecule.

    Parameters
    ----------
    polymer_object : object
        Polymerist or mBuild object that can provide an RDKit molecule.
    pdb_path : pathlib.Path
        PDB exported from the same final Polymerist/mBuild object. Its atom order
        is used as the authoritative sidecar atom order.
    required : bool, optional
        Raise on conversion failures when ``True``. When ``False``, failures are
        surfaced as warnings, by default ``True``.

    Returns
    -------
    Any or None
        RDKit molecule with atom ordering aligned to ``pdb_path``, or ``None``
        when optional conversion fails.

    Raises
    ------
    RuntimeError
        If required conversion fails.
    """
    to_rdkit = getattr(polymer_object, "to_rdkit", None)
    if not callable(to_rdkit):
        _handle_sdf_sidecar_failure(
            "Required RDKit SDF sidecar cannot be written because the polymer object does "
            "not expose to_rdkit()",
            required=required,
        )
        return None

    try:
        from rdkit import Chem
    except ImportError as exc:
        _handle_sdf_sidecar_failure(
            "RDKit is required to write the polymer SDF sidecar", required=required, cause=exc
        )
        return None

    try:
        mol = Chem.Mol(to_rdkit())
        _fix_tetravalent_neutral_nitrogens(mol)
        _set_unspecified_bonds_to_single(mol)
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol)
        aligned_mol = _align_rdkit_bond_orders_to_pdb_mol(mol, pdb_path)
        aligned_mol.UpdatePropertyCache(strict=False)
    except (AttributeError, RuntimeError, TypeError, ValueError) as exc:
        _handle_sdf_sidecar_failure(
            "Failed to convert Polymerist polymer to a PDB-aligned RDKit molecule",
            required=required,
            cause=exc,
        )
        return None

    return aligned_mol


def _set_unspecified_bonds_to_single(rdkit_mol: object) -> None:
    """Assign explicit single orders to mBuild/Polymerist topology-only bonds."""
    from rdkit import Chem

    for bond in rdkit_mol.GetBonds():
        if float(bond.GetBondTypeAsDouble()) <= 0.0:
            bond.SetBondType(Chem.BondType.SINGLE)


def _align_rdkit_bond_orders_to_pdb_mol(source_mol: object, pdb_path: Path) -> Any:
    """Return a PDB-ordered molecule carrying source-molecule bond orders."""
    from rdkit import Chem

    pdb_mol = Chem.MolFromPDBFile(
        str(pdb_path),
        sanitize=False,
        removeHs=False,
        proximityBonding=not _pdb_has_conect_records(pdb_path),
    )
    if pdb_mol is None or pdb_mol.GetNumAtoms() == 0:
        raise ValueError(f"RDKit could not read generated polymer PDB atoms from {pdb_path}")
    if source_mol.GetNumAtoms() != pdb_mol.GetNumAtoms():
        raise ValueError(
            "Polymerist RDKit/PDB atom count mismatch while writing SDF sidecar: "
            f"RDKit={source_mol.GetNumAtoms()}, PDB={pdb_mol.GetNumAtoms()}"
        )

    for atom_index, (source_atom, pdb_atom) in enumerate(
        zip(source_mol.GetAtoms(), pdb_mol.GetAtoms(), strict=True), start=1
    ):
        if source_atom.GetSymbol() != pdb_atom.GetSymbol():
            raise ValueError(
                "Polymerist RDKit/PDB atom order mismatch while writing SDF sidecar: "
                f"atom {atom_index} is {source_atom.GetSymbol()} in RDKit and "
                f"{pdb_atom.GetSymbol()} in PDB"
            )
        pdb_atom.SetFormalCharge(source_atom.GetFormalCharge())
        pdb_atom.SetNoImplicit(True)

    source_pairs = _bond_index_pairs(source_mol)
    pdb_pairs = _bond_index_pairs(pdb_mol)
    if source_pairs != pdb_pairs:
        missing = sorted(source_pairs - pdb_pairs)
        extra = sorted(pdb_pairs - source_pairs)
        raise ValueError(
            "Polymerist RDKit/PDB connectivity mismatch while writing SDF sidecar: "
            f"missing in PDB={missing[:5]}, extra in PDB={extra[:5]}"
        )

    editable = Chem.RWMol(pdb_mol)
    for source_bond in source_mol.GetBonds():
        pdb_bond = editable.GetBondBetweenAtoms(
            source_bond.GetBeginAtomIdx(), source_bond.GetEndAtomIdx()
        )
        if pdb_bond is None:
            raise ValueError("Polymerist PDB molecule lost a source bond during alignment")
        bond_order = float(source_bond.GetBondTypeAsDouble())
        if bond_order <= 0.0:
            raise ValueError("Polymerist source molecule contains a zero/unknown bond order")
        pdb_bond.SetBondType(source_bond.GetBondType())

    aligned = editable.GetMol()
    _fix_tetravalent_neutral_nitrogens(aligned)
    Chem.SanitizeMol(aligned)
    return aligned


def _bond_index_pairs(rdkit_mol: object) -> set[tuple[int, int]]:
    """Return zero-based RDKit bond-index pairs."""
    return {
        tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
        for bond in rdkit_mol.GetBonds()
    }


def _pdb_has_conect_records(path: Path) -> bool:
    """Return whether a generated PDB contains explicit CONECT records."""
    with path.open("r", encoding="utf-8") as handle:
        return any(line.startswith("CONECT") for line in handle)


def _fix_tetravalent_neutral_nitrogens(rdkit_mol: object) -> None:
    """Assign formal charge to tetravalent neutral nitrogen atoms in-place."""
    for atom in rdkit_mol.GetAtoms():
        if atom.GetSymbol() == "N" and atom.GetDegree() == 4 and atom.GetFormalCharge() == 0:
            atom.SetFormalCharge(1)


def _write_rdkit_sdf_sidecar(
    rdkit_mol: object,
    sdf_path: Path,
    *,
    required: bool = True,
) -> None:
    """Write an RDKit SDF sidecar carrying polymer bond orders.

    Parameters
    ----------
    rdkit_mol : object
        Explicit-hydrogen RDKit molecule to serialize.
    sdf_path : pathlib.Path
        Destination SDF path.
    required : bool, optional
        Raise on sidecar failures when ``True``. When ``False``, failures are
        surfaced as warnings, by default ``True``.

    Raises
    ------
    RuntimeError
        If the required sidecar cannot be generated or written.
    """
    try:
        from rdkit import Chem
    except ImportError as exc:
        _handle_sdf_sidecar_failure(
            "RDKit is required to write the polymer SDF sidecar", required=required, cause=exc
        )
        return

    try:
        Chem.MolToMolFile(rdkit_mol, str(sdf_path))
    except (AttributeError, OSError, RuntimeError, TypeError, ValueError) as exc:
        _handle_sdf_sidecar_failure(
            f"Failed to write required RDKit SDF sidecar to {sdf_path}",
            required=required,
            cause=exc,
        )
        return

    if not sdf_path.exists():
        _handle_sdf_sidecar_failure(
            f"Required RDKit SDF sidecar was not created at {sdf_path}",
            required=required,
        )


def _handle_sdf_sidecar_failure(
    message: str,
    *,
    required: bool,
    cause: BaseException | None = None,
) -> None:
    """Raise or warn for an SDF sidecar failure.

    Parameters
    ----------
    message : str
        Diagnostic message for the sidecar failure.
    required : bool
        Whether the caller requires the sidecar to continue.
    cause : BaseException or None, optional
        Original exception to chain when raising, by default ``None``.

    Raises
    ------
    RuntimeError
        If ``required`` is ``True``.
    """
    if required:
        raise RuntimeError(message) from cause

    import warnings

    warnings.warn(message, RuntimeWarning, stacklevel=2)
