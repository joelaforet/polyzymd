"""Internal moiety-source resolution for conjugation workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation._specs import _generated_fragment_from_moiety_plan
from polyzymd.builders.conjugation.polymer import (
    GeneratedMoietyFragment,
    GeneratedPolymerFragment,
    MultiResidueGenerationResult,
    PolymerRecipe,
    build_smiles_moiety_fragment,
    generate_multi_residue_molecule,
)
from polyzymd.builders.conjugation.polymer.polymerist import generated_fragment_from_polymerist_pdb


class ResolvedMoietySource(BaseModel):
    """Resolved conjugation moiety source ready for reaction planning."""

    model_config = {"arbitrary_types_allowed": True}

    fragment: GeneratedPolymerFragment = Field(exclude=True)
    source_fragment: Any | None = Field(default=None, exclude=True)
    source_kind: Literal["polymer", "smiles"]
    sidecars: dict[str, Path] = Field(default_factory=dict)
    generation: MultiResidueGenerationResult | None = None
    reactive_sequence_index: int | None = None
    reactive_selector: dict[str, int | str] | None = None
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)


def resolve_moiety_source(
    attachment: Any,
    *,
    attachment_index: int,
    output_dir: Path,
    protein_pdb_path: Path | str,
    force_regenerate: bool = False,
    max_retries: int = 3,
    energy_minimize: bool = True,
    random_seed: int | None = None,
) -> ResolvedMoietySource:
    """Resolve exactly one supported attachment moiety source.

    Parameters
    ----------
    attachment : Any
        Attachment config carrying a ``moiety`` object.
    attachment_index : int
        One-based attachment index used for deterministic sidecar paths.
    output_dir : pathlib.Path
        Directory for provider-generated artifacts.
    protein_pdb_path : pathlib.Path or str
        Prepared protein PDB used to map generated polymer reactive selectors.
    force_regenerate : bool, optional
        Regenerate polymer artifacts when using a polymer recipe, by default ``False``.
    max_retries : int, optional
        Multi-residue generation retry count, by default 3.
    energy_minimize : bool, optional
        Whether the generation backend should minimize generated polymers, by default ``True``.
    random_seed : int or None, optional
        Optional RDKit seed for SMILES moiety conformer generation, by default ``None``.

    Returns
    -------
    ResolvedMoietySource
        Provider-neutral fragment and provenance metadata.
    """
    moiety = getattr(attachment, "moiety", None)
    source_names = validate_moiety_source_config(moiety)
    if source_names[0] == "polymer_recipe":
        return _resolve_polymer_recipe_source(
            attachment,
            attachment_index=attachment_index,
            output_dir=output_dir,
            protein_pdb_path=protein_pdb_path,
            force_regenerate=force_regenerate,
            max_retries=max_retries,
            energy_minimize=energy_minimize,
        )
    return _resolve_smiles_source(
        attachment,
        attachment_index=attachment_index,
        output_dir=output_dir,
        random_seed=random_seed,
    )


def validate_moiety_source_config(moiety: Any) -> list[str]:
    """Validate a moiety source configuration without generating coordinates."""
    if moiety is None:
        raise ValueError("attachment.moiety is required")
    source_names = _configured_source_names(moiety)
    if len(source_names) != 1:
        joined = ", ".join(source_names) if source_names else "none"
        raise ValueError(
            "attachment.moiety must define exactly one supported source: "
            "polymer_recipe or smiles with residue_name. "
            f"Configured sources: {joined}"
        )
    if source_names[0] == "input_path":
        raise ValueError("attachment.moiety.input_path sources are not supported by this provider")
    return source_names


def _configured_source_names(moiety: Any) -> list[str]:
    """Return configured source kinds, including unsupported file inputs."""
    names = []
    if getattr(moiety, "polymer_recipe", None) is not None:
        names.append("polymer_recipe")
    if (
        getattr(moiety, "smiles", None) is not None
        or getattr(moiety, "residue_name", None) is not None
    ):
        if getattr(moiety, "smiles", None) is None or getattr(moiety, "residue_name", None) is None:
            raise ValueError(
                "SMILES moiety sources require both moiety.smiles and moiety.residue_name"
            )
        names.append("smiles")
    if getattr(moiety, "input_path", None) is not None:
        names.append("input_path")
    return names


def _resolve_polymer_recipe_source(
    attachment: Any,
    *,
    attachment_index: int,
    output_dir: Path,
    protein_pdb_path: Path | str,
    force_regenerate: bool,
    max_retries: int,
    energy_minimize: bool,
) -> ResolvedMoietySource:
    """Resolve a recipe-backed multi-residue polymer source."""
    recipe = getattr(attachment.moiety, "polymer_recipe", None)
    if not isinstance(recipe, PolymerRecipe):
        raise ValueError("attachment.moiety.polymer_recipe must define a PolymerRecipe")
    reactive_sequence_index = _reactive_sequence_index(recipe)
    generation = generate_multi_residue_molecule(
        recipe,
        output_dir / f"{attachment_index:02d}_{_safe_attachment_token(attachment.name)}",
        force_regenerate=force_regenerate,
        max_retries=max_retries,
        energy_minimize=energy_minimize,
    )
    if generation.pdb_path is None:
        raise RuntimeError("Generation backend did not produce a conjugate-polymer PDB")
    reactive_selector = _reactive_residue_selector(
        generation.pdb_path,
        sequence=generation.sequence,
        reactive_sequence_index=reactive_sequence_index,
    )
    fragment = generated_fragment_from_polymerist_pdb(
        generation.pdb_path,
        recipe=recipe,
        sequence=generation.sequence,
        name=attachment.moiety.name,
        reactive_residue_chain_id=_optional_str(reactive_selector.get("chain_id")),
        reactive_residue_name=str(reactive_selector["residue_name"]),
        reactive_residue_number=int(reactive_selector["residue_number"]),
    )
    sidecars = {}
    if generation.sdf_path is not None:
        sidecars["sdf"] = Path(generation.sdf_path)
        sidecars["bond_sdf"] = Path(generation.sdf_path)
    charged_sdf_path = getattr(generation, "charged_sdf_path", None)
    if charged_sdf_path is not None:
        sidecars["charged_sdf"] = Path(charged_sdf_path)
    return ResolvedMoietySource(
        fragment=fragment,
        source_fragment=fragment,
        source_kind="polymer",
        sidecars=sidecars,
        generation=generation,
        reactive_sequence_index=reactive_sequence_index,
        reactive_selector=reactive_selector,
        diagnostics=("Resolved polymer_recipe moiety source",),
    )


def _resolve_smiles_source(
    attachment: Any,
    *,
    attachment_index: int,
    output_dir: Path,
    random_seed: int | None,
) -> ResolvedMoietySource:
    """Resolve a SMILES-backed one-residue moiety source."""
    moiety = attachment.moiety
    source_fragment = build_smiles_moiety_fragment(
        moiety.smiles,
        moiety.residue_name,
        name=moiety.name,
        output_dir=output_dir / f"{attachment_index:02d}_{_safe_attachment_token(attachment.name)}",
        random_seed=random_seed,
    )
    fragment = _generated_fragment_from_moiety_plan_placeholder(source_fragment)
    sidecars = _moiety_sidecars(source_fragment)
    return ResolvedMoietySource(
        fragment=fragment,
        source_fragment=source_fragment,
        source_kind="smiles",
        sidecars=sidecars,
        diagnostics=("Resolved SMILES moiety source",),
    )


def _generated_fragment_from_moiety_plan_placeholder(
    fragment: GeneratedMoietyFragment,
) -> GeneratedPolymerFragment:
    """Adapt a moiety before the final resolved plan is available."""
    return GeneratedPolymerFragment(
        atoms=fragment.atoms,
        bonds=fragment.bonds,
        bond_orders=fragment.bond_orders,
        residues=fragment.residues,
        sequence=None,
        reactive_atom_serial=None,
        reactive_atom_index=None,
        reactive_atom_name=None,
        leaving_atom_serials=(),
        leaving_atom_indices=(),
        leaving_atom_names=(),
        name=fragment.name,
    )


def generated_fragment_for_resolved_source(
    source: ResolvedMoietySource,
    plan: Any,
) -> GeneratedPolymerFragment:
    """Return the construction fragment updated with resolved reactive atoms."""
    if isinstance(source.source_fragment, GeneratedMoietyFragment):
        return _generated_fragment_from_moiety_plan(source.source_fragment, plan)
    return source.fragment


def _moiety_sidecars(fragment: GeneratedMoietyFragment) -> dict[str, Path]:
    """Return generated SMILES sidecar paths."""
    sidecars = {}
    if fragment.pdb_path is not None:
        sidecars["pdb"] = fragment.pdb_path
    if fragment.sdf_path is not None:
        sidecars["sdf"] = fragment.sdf_path
    return sidecars


def _reactive_sequence_index(recipe: PolymerRecipe) -> int:
    """Return the configured reactive sequence index for a polymer recipe."""
    reactive_index = recipe.effective_reactive_index
    if reactive_index is None:
        raise ValueError("polymer recipe must define a reactive monomer index or label")
    return reactive_index


def _reactive_residue_selector(
    pdb_path: Path,
    *,
    sequence: str,
    reactive_sequence_index: int,
) -> dict[str, int | str]:
    """Return the PDB residue selector corresponding to a polymer sequence index."""
    from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord

    if reactive_sequence_index >= len(sequence):
        raise ValueError("reactive sequence index is outside the generated polymer sequence")
    residues: list[tuple[str, int, str, str]] = []
    seen: set[tuple[str, int, str, str]] = set()
    with Path(pdb_path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom = PdbAtomRecord.from_pdb_line(line, atom_index=0)
            key = (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)
            if key in seen:
                continue
            seen.add(key)
            residues.append(key)
    if not residues:
        raise ValueError(f"No ATOM/HETATM residues found in {pdb_path}")
    pdb_order_index = _polymerist_pdb_order_index(
        reactive_sequence_index,
        sequence_length=len(sequence),
    )
    if pdb_order_index >= len(residues):
        raise ValueError(
            "reactive sequence index maps outside the generated PDB residue list: "
            f"sequence_index={reactive_sequence_index}, pdb_order_index={pdb_order_index}, "
            f"residues={len(residues)}"
        )
    chain_id, residue_number, insertion_code, residue_name = residues[pdb_order_index]
    return {
        "sequence_index": reactive_sequence_index,
        "pdb_order_index": pdb_order_index,
        "label": sequence[reactive_sequence_index],
        "chain_id": chain_id,
        "residue_number": residue_number,
        "insertion_code": insertion_code,
        "residue_name": residue_name,
    }


def _polymerist_pdb_order_index(sequence_index: int, *, sequence_length: int) -> int:
    """Map user sequence order to Polymerist's PDB residue emission order."""
    if sequence_length <= 2:
        return sequence_index
    if sequence_index == 0:
        return sequence_length - 2
    if sequence_index == sequence_length - 1:
        return sequence_length - 1
    return sequence_index - 1


def _safe_attachment_token(name: str) -> str:
    """Return a conservative artifact-directory token."""
    token = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in name.strip())
    return token.strip("_-") or "attachment"


def _optional_str(value: object) -> str | None:
    """Return stripped text or None for empty values."""
    text = "" if value is None else str(value).strip()
    return text or None
