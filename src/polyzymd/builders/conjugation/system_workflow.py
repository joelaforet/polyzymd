"""Config-driven protein-polymer conjugate packing workflow."""

from __future__ import annotations

import copy
import inspect
import json
import logging
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
from pydantic import BaseModel, Field

from polyzymd.builders.conjugation._assembly import (
    ModifierConstructionResult,
    ModifierConstructionSettings,
    PackmolModifierPlacementSettings,
    place_modifiers_with_resolved_plans,
)
from polyzymd.builders.conjugation._linkage import (
    NhsLysModifierLinker,
    PabloCrosslinkRequirement,
    ResolvedAttachmentPlan,
    parse_pdb_atom_records,
    placed_fragment_from_resolved_plan,
    require_pablo_crosslink_requirement,
)
from polyzymd.builders.conjugation._relaxation import (
    VacuumSmokeSettings,
    run_post_crosslink_local_minimization,
    run_restrained_vacuum_smoke,
)
from polyzymd.builders.conjugation._specs import (
    AttachmentBuildSpec,
    attachment_spec_from_generated_polymer_plan,
    attachment_spec_from_moiety_plan,
)
from polyzymd.builders.conjugation.final_interchange import create_final_conjugated_interchange
from polyzymd.builders.conjugation.models import ConjugationResult
from polyzymd.builders.conjugation.pablo.ingestion import PabloIngestor
from polyzymd.builders.conjugation.pablo.parameterization import (
    InterchangeParameterizationSettings,
    build_formal_charge_smoke_template,
    create_interchange_from_pablo_topology,
)
from polyzymd.builders.conjugation.polymer import (
    GeneratedMoietyFragment,
    GeneratedPolymerFragment,
    PolymeristGenerationSmokeResult,
    PolymerRecipe,
    build_smiles_moiety_fragment,
    generate_polymerist_smoke_polymer,
    generated_fragment_from_polymerist_pdb,
)
from polyzymd.builders.conjugation.reactions._roles import (
    STRUCTURE_MATCHING_BLOCKER_MESSAGE,
    atom_mapped_reaction_from_mechanism_config,
    resolve_reaction_roles_from_identity_map,
)
from polyzymd.builders.conjugation.reactions.library import get_reaction
from polyzymd.builders.conjugation.structure.pdb import (
    CrosslinkedPdbAssemblyOptions,
    PdbAtomRecord,
    write_crosslinked_pdb,
)
from polyzymd.builders.conjugation.structure.preparation import (
    ProteinCanonicalizationResult,
    ProteinCanonicalizationSettings,
    canonicalize_protein_hydrogens,
)
from polyzymd.builders.system_builder import SystemBuilder
from polyzymd.config.schema import (
    ConjugationCcdCrosslinkConfig,
    ConjugationCcdPabloPolicyConfig,
    ConjugationChainPolicyConfig,
    ConjugationConfig,
    SimulationConfig,
)

_ATOM_RECORD_PREFIXES = ("ATOM", "HETATM")
_NHS_LYS_REACTION = get_reaction("nhs_lys")
_NHS_LYS_REACTION_NAME = _NHS_LYS_REACTION.name
_NHS_LYS_COORDINATE_BACKEND_MECHANISM = _NHS_LYS_REACTION.coordinate_backend_mechanism
LOGGER = logging.getLogger(__name__)


class ConjugatedPolymerSystemSettings(BaseModel):
    """Settings for config-driven conjugate construction, relaxation, and solvation."""

    model_config = {"arbitrary_types_allowed": True}

    force_regenerate_conjugate_polymer: bool = False
    conjugate_polymerist_max_retries: int = Field(3, ge=1)
    conjugate_polymerist_energy_minimize: bool = True
    conjugate_cache_dir_name: str = "conjugate-polymerist-cache"
    conjugate_artifact_dir_name: str = "conjugate-construction"
    solvated_pdb_name: str = "solvated_conjugate_free_polymers.pdb"
    workflow_json_name: str = "conjugated_polymer_system_workflow.json"
    create_final_interchange: bool = False
    preserve_reference_atom_names: bool = True
    canonicalize_source_protein_hydrogens: bool = True
    use_product_state_pablo_library: bool = True
    run_product_state_local_minimization: bool = True
    protein_canonicalization: ProteinCanonicalizationSettings = Field(
        default_factory=ProteinCanonicalizationSettings
    )
    local_minimization: Any = Field(default_factory=lambda: _default_local_minimization_settings())
    placement: PackmolModifierPlacementSettings = Field(
        default_factory=PackmolModifierPlacementSettings
    )
    conjugate_parameterization: InterchangeParameterizationSettings = Field(
        default_factory=InterchangeParameterizationSettings
    )
    vacuum_smoke: VacuumSmokeSettings = Field(default_factory=lambda: _protein_restrained_smoke())


class ConjugateConstructionResult(BaseModel):
    """Specs-first construction result with singular compatibility fields."""

    model_config = {"arbitrary_types_allowed": True}

    output_dir: Path
    resolved_plan: ResolvedAttachmentPlan
    resolved_plans: tuple[ResolvedAttachmentPlan, ...]
    attachment_specs: tuple[Any, ...] = Field(default_factory=tuple, exclude=True)
    crosslink_validation: Any
    crosslink_validations: tuple[Any, ...]
    placement: Any
    placements: tuple[Any, ...]
    assembly: Any
    pablo: Any
    parameterization: Any
    smoke: Any | None = None
    local_minimization: Any | None = None
    product_state_pablo_library: Any | None = Field(default=None, exclude=True)
    crosslinked_pdb_path: Path
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)


def build_conjugated_polymer_system_from_config_path(
    config_path: Path | str,
    *,
    output_dir: Path | str | None = None,
    settings: ConjugatedPolymerSystemSettings | None = None,
    free_polymer_seed: int | None = None,
) -> ConjugationResult:
    """Load a config YAML and build the relaxed, solvated conjugate system."""
    from polyzymd.config.loader import load_config

    path = Path(config_path)
    config = load_config(path)
    effective_output_dir = output_dir or path.parent / "artifacts" / path.stem
    LOGGER.info(
        "Starting config conjugation build from %s in %s",
        path,
        effective_output_dir,
    )
    return build_conjugated_polymer_system_from_config(
        config,
        output_dir=effective_output_dir,
        settings=settings,
        free_polymer_seed=free_polymer_seed,
    )


def build_conjugated_polymer_system_from_config(
    config: SimulationConfig,
    *,
    output_dir: Path | str,
    settings: ConjugatedPolymerSystemSettings | None = None,
    free_polymer_seed: int | None = None,
) -> ConjugationResult:
    """Build a relaxed protein-polymer conjugate, pack free polymers, and solvate.

    The workflow consumes the existing PolyzyMD config schema: ``conjugation``
    defines one covalent protein-polymer attachment, ``polymers`` defines free
    non-covalent polymer chains to pack around that conjugate, and ``solvent``
    defines the final solvent box.
    """
    workflow_settings = settings or ConjugatedPolymerSystemSettings()
    artifact_dir = Path(output_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    construction_dir = artifact_dir / workflow_settings.conjugate_artifact_dir_name

    attachments = _enabled_supported_nhs_lys_attachments(config.conjugation)
    LOGGER.info("Enabled conjugation attachment count: %d", len(attachments))
    _log_attachment_additions(attachments)
    LOGGER.info("Preparing and canonicalizing source protein")
    protein_pdb_path, protein_canonicalization = _prepared_protein_pdb_path(
        config.enzyme.pdb_path,
        output_dir=construction_dir,
        settings=workflow_settings,
    )
    spec_payloads = tuple(
        _build_nhs_lys_attachment_spec(
            attachment,
            attachment_index=index,
            protein_pdb_path=protein_pdb_path,
            artifact_dir=artifact_dir,
            workflow_settings=workflow_settings,
        )
        for index, attachment in enumerate(attachments, start=1)
    )
    specs = tuple(payload[0] for payload in spec_payloads)
    generations = tuple(payload[1] for payload in spec_payloads)
    reactive_sequence_indices = tuple(payload[2] for payload in spec_payloads)
    reactive_selectors = tuple(payload[3] for payload in spec_payloads)
    modifiers = tuple(spec.generated_fragment for spec in specs)
    resolved_plans = tuple(spec.resolved_plan for spec in specs)
    ccd_pablo_policy = _policy_with_resolved_crosslinks(
        config.conjugation.ccd_pablo,
        resolved_plans,
    )

    construction_settings = ModifierConstructionSettings(
        placement=workflow_settings.placement,
        parameterization=workflow_settings.conjugate_parameterization,
        smoke=workflow_settings.vacuum_smoke,
        run_smoke=True,
    )
    construction, construction_topology = _construct_conjugate_from_specs(
        protein_pdb_path=protein_pdb_path,
        specs=specs,
        ccd_pablo_policy=ccd_pablo_policy,
        output_dir=construction_dir,
        chain_policy=config.conjugation.chain_policy,
        settings=construction_settings,
        use_product_state_pablo_library=workflow_settings.use_product_state_pablo_library,
        run_product_state_local_minimization=(
            workflow_settings.run_product_state_local_minimization
        ),
        local_minimization_settings=workflow_settings.local_minimization,
    )
    if workflow_settings.preserve_reference_atom_names:
        _restore_smoke_pdb_atom_names(construction, construction.crosslinked_pdb_path)

    LOGGER.info("Preparing relaxed conjugate topology")
    relaxed_pdb = _relaxed_conjugate_pdb(construction)
    relaxed_topology = topology_with_pdb_positions(
        construction_topology,
        relaxed_pdb,
        atom_name_template_pdb=(
            construction.crosslinked_pdb_path
            if workflow_settings.preserve_reference_atom_names
            else None
        ),
    )

    builder = _build_solvated_system(
        config,
        relaxed_conjugate_topology=relaxed_topology,
        working_dir=artifact_dir,
        polymer_seed=free_polymer_seed,
        create_interchange=workflow_settings.create_final_interchange,
        product_state_pablo_library=getattr(construction, "product_state_pablo_library", None),
        parameterization_settings=workflow_settings.conjugate_parameterization,
    )
    solvated_pdb_path = artifact_dir / workflow_settings.solvated_pdb_name
    builder.save_topology(solvated_pdb_path)
    LOGGER.info("Wrote final solvated conjugate PDB to %s", solvated_pdb_path)
    if workflow_settings.preserve_reference_atom_names:
        _restore_pdb_atom_name_fields(solvated_pdb_path, construction.crosslinked_pdb_path)

    result = ConjugationResult(
        output_dir=artifact_dir,
        generated_sequence=generations[0].sequence,
        reactive_sequence_index=reactive_sequence_indices[0],
        reactive_residue_selector=reactive_selectors[0],
        conjugate_generation=generations[0],
        construction=construction,
        attachment_specs=specs,
        generated_sequences=tuple(generation.sequence for generation in generations),
        reactive_sequence_indices=reactive_sequence_indices,
        reactive_residue_selectors=reactive_selectors,
        conjugate_generations=generations,
        protein_canonicalization=protein_canonicalization,
        relaxed_conjugate_pdb_path=relaxed_pdb,
        solvated_pdb_path=solvated_pdb_path,
        final_interchange_created=builder.interchange is not None,
        modifier=modifiers[0],
        modifiers=modifiers,
        relaxed_conjugate_topology=relaxed_topology,
        solvated_topology=builder.solvated_topology,
        final_interchange=builder.interchange,
        system_builder=builder,
    )
    workflow_path = artifact_dir / workflow_settings.workflow_json_name
    result.workflow_json_path = workflow_path
    result.artifact_paths["workflow_json"] = workflow_path
    result.save(workflow_path)
    LOGGER.info("Saved conjugation workflow JSON to %s", workflow_path)
    LOGGER.info("Completed config conjugation build in %s", artifact_dir)
    return result


def build_direct_smiles_moiety_conjugate(
    *,
    protein_pdb_path: Path | str,
    attachments: tuple[Any, ...],
    output_dir: Path | str,
    ccd_pablo: Any | None = None,
    chain_policy: Any | None = None,
    settings: ConjugatedPolymerSystemSettings | None = None,
    random_seed: int | None = None,
) -> ConjugationResult:
    """Build a direct protein plus SMILES-moiety conjugate request.

    This public-engine path is intentionally narrower than the legacy config
    workflow: the input protein is already cleaned, each enabled attachment must
    define one SMILES moiety, and the combined product is parameterized once.
    """
    workflow_settings = settings or ConjugatedPolymerSystemSettings()
    artifact_dir = Path(output_dir)
    construction_dir = artifact_dir / workflow_settings.conjugate_artifact_dir_name
    moiety_dir = construction_dir / "moieties"
    construction_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Starting direct conjugation build in %s", artifact_dir)

    enabled_attachments = tuple(attachment for attachment in attachments if attachment.enabled)
    if not enabled_attachments:
        raise ValueError("Direct conjugation requests require at least one enabled attachment")
    LOGGER.info("Enabled conjugation attachment count: %d", len(enabled_attachments))
    _log_attachment_additions(enabled_attachments)

    LOGGER.info("Preparing and canonicalizing source protein")
    protein_path, protein_canonicalization = _prepared_protein_pdb_path(
        protein_pdb_path,
        output_dir=construction_dir,
        settings=workflow_settings,
    )
    specs: list[AttachmentBuildSpec] = []
    for index, attachment in enumerate(enabled_attachments, start=1):
        specs.append(
            _build_n_gly_direct_moiety_attachment_spec(
                attachment,
                attachment_index=index,
                protein_pdb_path=protein_path,
                moiety_dir=moiety_dir,
                random_seed=random_seed,
            )
        )
    resolved_plans = tuple(spec.resolved_plan for spec in specs)

    policy = _policy_with_resolved_crosslinks(
        ccd_pablo or ConjugationCcdPabloPolicyConfig(),
        resolved_plans,
    )
    chain_assignment = chain_policy or ConjugationChainPolicyConfig()
    construction_settings = ModifierConstructionSettings(
        placement=workflow_settings.placement,
        parameterization=workflow_settings.conjugate_parameterization,
        smoke=workflow_settings.vacuum_smoke,
        run_smoke=True,
    )

    construction, construction_topology = _construct_conjugate_from_specs(
        protein_pdb_path=protein_path,
        specs=tuple(specs),
        ccd_pablo_policy=policy,
        output_dir=construction_dir,
        chain_policy=chain_assignment,
        settings=construction_settings,
        use_product_state_pablo_library=workflow_settings.use_product_state_pablo_library,
        run_product_state_local_minimization=(
            workflow_settings.run_product_state_local_minimization
        ),
        local_minimization_settings=workflow_settings.local_minimization,
    )

    LOGGER.info("Preparing relaxed conjugate topology")
    relaxed_pdb = _relaxed_conjugate_pdb(construction)
    relaxed_topology = topology_with_pdb_positions(construction_topology, relaxed_pdb)
    builder = _build_direct_solvated_system(
        relaxed_conjugate_topology=relaxed_topology,
        working_dir=artifact_dir,
        create_interchange=workflow_settings.create_final_interchange,
        product_state_pablo_library=getattr(construction, "product_state_pablo_library", None),
        parameterization_settings=workflow_settings.conjugate_parameterization,
    )
    solvated_pdb_path = artifact_dir / workflow_settings.solvated_pdb_name
    builder.save_topology(solvated_pdb_path)
    LOGGER.info("Wrote final solvated conjugate PDB to %s", solvated_pdb_path)
    if workflow_settings.preserve_reference_atom_names:
        _restore_pdb_atom_name_fields(solvated_pdb_path, relaxed_pdb)

    result = ConjugationResult(
        output_dir=artifact_dir,
        construction=construction,
        protein_canonicalization=protein_canonicalization,
        relaxed_conjugate_pdb_path=relaxed_pdb,
        solvated_pdb_path=solvated_pdb_path,
        final_interchange_created=builder.interchange is not None,
        modifier=tuple(spec.generated_fragment for spec in specs),
        modifiers=tuple(spec.generated_fragment for spec in specs),
        attachment_specs=tuple(specs),
        relaxed_conjugate_topology=relaxed_topology,
        solvated_topology=builder.solvated_topology,
        final_interchange=builder.interchange,
        system_builder=builder,
    )
    workflow_path = artifact_dir / workflow_settings.workflow_json_name
    result.workflow_json_path = workflow_path
    result.artifact_paths["workflow_json"] = workflow_path
    result.save(workflow_path)
    _save_direct_workflow_summary(result, enabled_attachments, list(resolved_plans), workflow_path)
    LOGGER.info("Saved conjugation workflow JSON to %s", workflow_path)
    LOGGER.info("Completed direct conjugation build in %s", artifact_dir)
    return result


def topology_with_pdb_positions(
    topology: Any,
    pdb_path: Path | str,
    *,
    atom_name_template_pdb: Path | str | None = None,
) -> Any:
    """Return a topology copy whose atom positions come from a same-order PDB."""
    positioned_topology = copy.deepcopy(topology)
    coords_angstrom = _pdb_coordinates_angstrom(Path(pdb_path))
    n_atoms = getattr(positioned_topology, "n_atoms", None)
    if n_atoms is not None and len(coords_angstrom) != int(n_atoms):
        raise ValueError(
            "Relaxed PDB atom count does not match Pablo topology atom count: "
            f"PDB={len(coords_angstrom)}, topology={n_atoms}"
        )
    from openff.units import Quantity

    positioned_topology.set_positions(Quantity(coords_angstrom, "angstrom"))
    if atom_name_template_pdb is not None:
        _apply_pdb_atom_names_to_topology(positioned_topology, atom_name_template_pdb)
    return positioned_topology


def _protein_restrained_smoke() -> VacuumSmokeSettings:
    """Default short vacuum run that restrains protein heavy atoms only."""
    return VacuumSmokeSettings(
        minimization_max_iterations=100,
        nvt_steps=10,
        restrain_all_heavy_atoms=False,
    )


def _default_local_minimization_settings() -> Any:
    """Build default post-crosslink local minimization settings lazily."""
    from polyzymd.builders.conjugation._relaxation import LocalMinimizationSettings

    return LocalMinimizationSettings()


def _prepared_protein_pdb_path(
    protein_pdb_path: Path | str,
    *,
    output_dir: Path,
    settings: ConjugatedPolymerSystemSettings,
) -> tuple[Path, ProteinCanonicalizationResult | None]:
    """Return the protein PDB used for construction, canonicalizing hydrogens if enabled."""
    source_path = Path(protein_pdb_path)
    if not settings.canonicalize_source_protein_hydrogens:
        return source_path, None

    canonical_settings = settings.protein_canonicalization
    output_path = canonical_settings.output_path
    if output_path is None:
        output_path = output_dir / f"source_protein_canonical_pH{canonical_settings.ph:g}.pdb"
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    result = canonicalize_protein_hydrogens(
        source_path,
        output_path,
        settings=canonical_settings,
    )
    result.save(Path(result.output_path).with_suffix(".canonicalization.json"))
    return result.output_path, result


def _build_nhs_lys_attachment_spec(
    attachment: Any,
    *,
    attachment_index: int,
    protein_pdb_path: Path | str,
    artifact_dir: Path,
    workflow_settings: ConjugatedPolymerSystemSettings,
) -> tuple[AttachmentBuildSpec, PolymeristGenerationSmokeResult, int, dict[str, int | str], Any]:
    """Resolve one config-driven NHS-Lys polymer attachment into a build spec."""
    _require_supported_coordinate_backend(attachment)
    recipe = _polymer_recipe_from_attachment(attachment)
    reactive_sequence_index = _reactive_sequence_index(recipe)
    LOGGER.info(
        "Generating conjugate polymer/moiety for attachment %d (%s)",
        attachment_index,
        _attachment_moiety_name(attachment),
    )
    generation = generate_polymerist_smoke_polymer(
        recipe,
        artifact_dir
        / workflow_settings.conjugate_cache_dir_name
        / f"{attachment_index:02d}_{_safe_attachment_token(attachment.name)}",
        force_regenerate=workflow_settings.force_regenerate_conjugate_polymer,
        max_retries=workflow_settings.conjugate_polymerist_max_retries,
        energy_minimize=workflow_settings.conjugate_polymerist_energy_minimize,
    )
    if generation.pdb_path is None:
        raise RuntimeError("Polymerist did not produce a conjugate-polymer PDB")

    LOGGER.info("Resolving linkage for attachment %d", attachment_index)
    reactive_selector = _reactive_residue_selector(
        generation.pdb_path,
        sequence=generation.sequence,
        reactive_sequence_index=reactive_sequence_index,
    )
    modifier = generated_fragment_from_polymerist_pdb(
        generation.pdb_path,
        recipe=recipe,
        sequence=generation.sequence,
        name=attachment.moiety.name,
        reactive_residue_chain_id=_optional_str(reactive_selector.get("chain_id")),
        reactive_residue_name=str(reactive_selector["residue_name"]),
        reactive_residue_number=int(reactive_selector["residue_number"]),
    )
    linker = _nhs_lys_linker_from_attachment(attachment)
    resolved_plan = linker.resolve_plan(protein_pdb_path, modifier)
    return (
        attachment_spec_from_generated_polymer_plan(
            modifier,
            generation.sdf_path,
            resolved_plan,
            attachment_config=attachment,
            attachment_index=attachment_index,
            reaction_name=attachment.mechanism.name,
            charged_sdf_path=getattr(generation, "charged_sdf_path", None),
        ),
        generation,
        reactive_sequence_index,
        reactive_selector,
        linker,
    )


def _build_n_gly_direct_moiety_attachment_spec(
    attachment: Any,
    *,
    attachment_index: int,
    protein_pdb_path: Path | str,
    moiety_dir: Path,
    random_seed: int | None,
) -> AttachmentBuildSpec:
    """Resolve one direct SMILES N-gly moiety attachment into a build spec."""
    reaction_template = get_reaction(attachment.mechanism.name)
    if reaction_template.coordinate_backend_mechanism != "n_glycosylation":
        raise NotImplementedError(
            "Direct SMILES-moiety requests currently support mechanism "
            f"'n_glycosylation'; received '{attachment.mechanism.name}'"
        )
    moiety = attachment.moiety
    if moiety.smiles is None or moiety.residue_name is None:
        raise ValueError("Direct SMILES-moiety attachments require moiety.smiles and residue_name")
    LOGGER.info(
        "Generating conjugate polymer/moiety for attachment %d (%s)",
        attachment_index,
        _attachment_moiety_name(attachment),
    )
    fragment = build_smiles_moiety_fragment(
        moiety.smiles,
        moiety.residue_name,
        name=moiety.name,
        output_dir=moiety_dir / f"{attachment_index:02d}_{_safe_attachment_token(attachment.name)}",
        random_seed=random_seed,
    )
    settings_builder = getattr(reaction_template, "settings_from_attachment", None)
    reaction_settings = settings_builder(attachment) if callable(settings_builder) else None
    LOGGER.info("Resolving linkage for attachment %d", attachment_index)
    plan = reaction_template.resolve_plan(
        protein_pdb_path,
        attachment.site,
        fragment,
        settings=reaction_settings,
    )
    return attachment_spec_from_moiety_plan(
        fragment,
        plan,
        attachment_config=attachment,
        attachment_index=attachment_index,
        reaction_name=attachment.mechanism.name,
    )


def _enabled_supported_nhs_lys_attachments(
    conjugation: ConjugationConfig | None,
) -> tuple[Any, ...]:
    """Return all enabled attachments supported by the config NHS-Lys path.

    Parameters
    ----------
    conjugation : ConjugationConfig or None
        Conjugation configuration to inspect.

    Returns
    -------
    tuple of Any
        Enabled attachment configurations supported by this workflow.
    """
    if conjugation is None or not conjugation.enabled:
        raise ValueError("conjugation.enabled must be true for this workflow")
    attachments = tuple(attachment for attachment in conjugation.attachments if attachment.enabled)
    if not attachments:
        raise ValueError("conjugated polymer workflow requires at least one enabled attachment")
    for attachment in attachments:
        _require_supported_coordinate_backend(attachment)
        _polymer_recipe_from_attachment(attachment)
    return attachments


def _log_attachment_additions(attachments: tuple[Any, ...]) -> None:
    """Log the enabled attachment additions in user-facing site notation.

    Parameters
    ----------
    attachments : tuple of Any
        Enabled attachment-like objects with ``site`` and ``moiety`` attributes.
    """
    total = len(attachments)
    for index, attachment in enumerate(attachments, start=1):
        site = getattr(attachment, "site", None)
        LOGGER.info(
            "Addition %d/%d: adding %s to %s:%s%s:%s",
            index,
            total,
            _attachment_moiety_name(attachment),
            _site_value(site, "chain_id", "?"),
            _site_value(site, "residue_name", ""),
            _site_value(site, "residue_number", "?"),
            _site_value(site, "atom_name", "?"),
        )


def _attachment_moiety_name(attachment: Any) -> str:
    """Return a readable moiety name for logging."""
    moiety = getattr(attachment, "moiety", None)
    name = getattr(moiety, "name", None) or getattr(attachment, "name", None)
    return str(name or "moiety")


def _site_value(site: Any, field_name: str, default: str) -> str:
    """Return a non-empty site field value for logging."""
    value = getattr(site, field_name, None)
    if value is None or value == "":
        return default
    return str(value)


def _polymer_recipe_from_attachment(attachment: Any) -> PolymerRecipe:
    recipe = getattr(attachment.moiety, "polymer_recipe", None)
    if not isinstance(recipe, PolymerRecipe):
        raise ValueError("attachment.moiety.recipe must define a PolymerRecipe")
    return recipe


def _reactive_sequence_index(recipe: PolymerRecipe) -> int:
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
    if reactive_sequence_index >= len(sequence):
        raise ValueError("reactive sequence index is outside the generated polymer sequence")
    residues = _unique_pdb_residues(pdb_path)
    pdb_order_index = _polymerist_pdb_order_index(
        reactive_sequence_index,
        sequence_length=len(sequence),
    )
    if pdb_order_index >= len(residues):
        raise ValueError(
            "reactive sequence index maps outside the Polymerist PDB residue list: "
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
    """Map user sequence order to Polymerist's PDB residue emission order.

    ``PolymerGenerator`` gives Polymerist a middle sequence plus separately
    configured head/tail fragments. Polymerist emits middle residues first,
    followed by the head and tail terminal fragments. This mapper lets a config
    still use normal sequence indexing for the covalent reactive residue.
    """
    if sequence_length <= 2:
        return sequence_index
    if sequence_index == 0:
        return sequence_length - 2
    if sequence_index == sequence_length - 1:
        return sequence_length - 1
    return sequence_index - 1


def _require_supported_coordinate_backend(attachment: Any) -> None:
    """Gate config-driven coordinate construction to implemented backends."""
    mechanism = attachment.mechanism
    mechanism_name = mechanism.name.strip().lower()
    if (
        mechanism_name == _NHS_LYS_COORDINATE_BACKEND_MECHANISM
        and mechanism.reaction_smarts is None
    ):
        return

    if mechanism.reaction_smarts is not None:
        reaction = atom_mapped_reaction_from_mechanism_config(mechanism)
        preflight = resolve_reaction_roles_from_identity_map(
            reaction,
            {},
            require_required_identities=False,
        )
        added = len(preflight.bond_changes.added_bonds)
        removed = len(preflight.bond_changes.removed_bonds)
        order_changed = len(preflight.bond_changes.order_changes)
        raise NotImplementedError(
            f"{STRUCTURE_MATCHING_BLOCKER_MESSAGE} Mechanism '{mechanism.name}' passed "
            f"generic SMARTS preflight with {added} added, {removed} removed, and "
            f"{order_changed} order-changed mapped bonds, but the config-driven system "
            "workflow currently has coordinate surgery only for mechanism "
            f"'{_NHS_LYS_COORDINATE_BACKEND_MECHANISM}' without reaction_smarts."
        )

    raise NotImplementedError(
        "Config-driven conjugated polymer system construction currently implements coordinate "
        f"surgery only for mechanism '{_NHS_LYS_COORDINATE_BACKEND_MECHANISM}'. "
        f"Received mechanism '{mechanism.name}'. Provide reaction_smarts for a generic preflight "
        "or use the supported NHS-Lys mechanism."
    )


def _unique_pdb_residues(path: Path) -> tuple[tuple[str, int, str, str], ...]:
    residues: list[tuple[str, int, str, str]] = []
    seen: set[tuple[str, int, str, str]] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(_ATOM_RECORD_PREFIXES):
                continue
            atom = PdbAtomRecord.from_pdb_line(line, atom_index=0)
            key = (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)
            if key in seen:
                continue
            seen.add(key)
            residues.append(key)
    if not residues:
        raise ValueError(f"No ATOM/HETATM residues found in {path}")
    return tuple(residues)


def _nhs_lys_linker_from_attachment(attachment: Any) -> NhsLysModifierLinker:
    reaction_template: Any = get_reaction("nhs_lys")
    return reaction_template.create_linker_from_attachment(attachment)


def _policy_with_resolved_crosslink(
    policy: Any,
    resolved_plan: ResolvedAttachmentPlan,
) -> Any:
    requirement = _product_state_crosslink_requirement(resolved_plan)
    try:
        require_pablo_crosslink_requirement(policy, requirement)
        return policy
    except ValueError:
        pass

    crosslink = ConjugationCcdCrosslinkConfig(
        residues=requirement.residues,
        linking_atoms=requirement.linking_atoms,
        leaving_atoms=requirement.leaving_atoms,
        bond_order=requirement.bond_order,
    )
    existing = list(getattr(policy, "crosslinks", ()))
    if hasattr(policy, "model_copy"):
        return policy.model_copy(update={"crosslinks": [*existing, crosslink]})

    from types import SimpleNamespace

    data = getattr(policy, "__dict__", {}).copy() if policy is not None else {}
    data["crosslinks"] = [*existing, crosslink]
    return SimpleNamespace(**data)


def _product_state_crosslink_requirement(resolved_plan: ResolvedAttachmentPlan):
    """Return the Pablo crosslink requirement for an already-modified product PDB.

    The resolved attachment plan records reactant-state leaving atoms for graph
    surgery. By the time Pablo ingests the emitted crosslinked PDB, PolyzyMD has
    already removed those atoms and written the product-state ``CONECT`` bond, so
    the Pablo policy must not ask Pablo to remove them a second time.
    """
    requirement = resolved_plan.pablo_crosslink_requirement
    return requirement.model_copy(update={"leaving_atoms": ((), ())})


def _construct_nhs_lys_modifier_linked_protein(
    *,
    protein_pdb_path: Path | str,
    prepared_protein_pdb_path: Path | str,
    modifier: GeneratedPolymerFragment,
    polymer_sdf_path: Path | str | None,
    linker: NhsLysModifierLinker,
    resolved_plan: ResolvedAttachmentPlan,
    ccd_pablo_policy: Any,
    output_dir: Path | str,
    chain_policy: Any | None,
    settings: ModifierConstructionSettings,
    use_product_state_pablo_library: bool,
    run_product_state_local_minimization: bool,
    local_minimization_settings: Any,
) -> tuple[ModifierConstructionResult, Any]:
    """Compatibility shim for the single-site NHS-Lys private constructor."""
    _ = (
        protein_pdb_path,
        linker,
    )
    spec = SimpleNamespace(
        attachment_id="nhs_lys_attachment_01",
        attachment_index=1,
        reaction_name=_NHS_LYS_REACTION_NAME,
        generated_fragment=modifier,
        resolved_plan=resolved_plan,
        source_sidecars={"sdf": Path(polymer_sdf_path)} if polymer_sdf_path is not None else {},
        fragment=SimpleNamespace(source_kind="polymer"),
    )
    return _construct_conjugate_from_specs(
        protein_pdb_path=prepared_protein_pdb_path,
        specs=(spec,),
        ccd_pablo_policy=ccd_pablo_policy,
        output_dir=output_dir,
        chain_policy=chain_policy,
        settings=settings,
        use_product_state_pablo_library=use_product_state_pablo_library,
        run_product_state_local_minimization=run_product_state_local_minimization,
        local_minimization_settings=local_minimization_settings,
    )


def _construct_conjugate_from_specs(
    *,
    protein_pdb_path: Path | str,
    specs: tuple[AttachmentBuildSpec, ...],
    ccd_pablo_policy: Any,
    output_dir: Path | str,
    chain_policy: Any | None,
    settings: ModifierConstructionSettings,
    use_product_state_pablo_library: bool,
    run_product_state_local_minimization: bool = False,
    local_minimization_settings: Any | None = None,
) -> tuple[Any, Any]:
    """Construct, parameterize, and relax a conjugate from resolved attachment specs."""
    if not specs:
        raise ValueError("Conjugate construction requires at least one attachment spec")

    modifiers = tuple(spec.generated_fragment for spec in specs)
    resolved_plans = tuple(spec.resolved_plan for spec in specs)
    if len(modifiers) != len(resolved_plans):
        raise ValueError("Conjugate construction requires aligned attachment specs")

    artifact_dir = Path(output_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info("Resolving product-state Pablo crosslinks")
    product_state_requirements = tuple(
        _product_state_crosslink_requirement(plan) for plan in resolved_plans
    )
    crosslink_validations = tuple(
        require_pablo_crosslink_requirement(ccd_pablo_policy, requirement)
        for requirement in product_state_requirements
    )

    LOGGER.info("Placing conjugate polymer/moiety fragments")
    placements = place_modifiers_with_resolved_plans(
        protein_pdb_path,
        modifiers,
        resolved_plans,
        artifact_dir,
        settings=settings.placement,
    )
    placed_modifiers = tuple(
        placed_fragment_from_resolved_plan(
            placement.placed_modifier,
            plan,
        )
        for placement, plan in zip(placements, resolved_plans, strict=True)
    )

    crosslinked_pdb_path = artifact_dir / settings.crosslinked_pdb_name
    LOGGER.info("Writing crosslinked PDB to %s", crosslinked_pdb_path)
    assembly_result = write_crosslinked_pdb(
        protein_pdb_path,
        tuple(placed_modifiers),
        tuple(plan.to_pdb_linkage_attachment() for plan in resolved_plans),
        crosslinked_pdb_path,
        CrosslinkedPdbAssemblyOptions(),
    )

    product_state_pablo_library = None
    product_state_residue_library = None
    if use_product_state_pablo_library:
        LOGGER.info("Building product-state Pablo residue library")
        product_state_specs = _product_state_specs_with_assembly_mappings(
            specs,
            assembly_result=assembly_result,
        )
        product_state_pablo_library = _product_state_pablo_library_for_specs(
            product_pdb=crosslinked_pdb_path,
            source_protein_pdb=protein_pdb_path,
            specs=product_state_specs,
        )
        product_state_residue_library = product_state_pablo_library.residue_library

    LOGGER.info("Ingesting product-state PDB with Pablo")
    pablo_result = PabloIngestor(policy=ccd_pablo_policy).ingest_structure(
        crosslinked_pdb_path,
        chain_policy=chain_policy,
        output_dir=artifact_dir,
        residue_library=product_state_residue_library,
    )
    if not pablo_result.success or pablo_result.topology is None:
        raise RuntimeError(_pablo_failure_message(pablo_result))
    if product_state_pablo_library is not None:
        LOGGER.info("Building production product-state charge bridge")
        product_state_pablo_library = _product_state_library_with_charge_bridge(
            product_state_pablo_library,
            product_topology=pablo_result.topology,
            product_pdb=crosslinked_pdb_path,
            source_protein_pdb=protein_pdb_path,
            specs=product_state_specs,
            output_dir=artifact_dir,
            settings=settings.parameterization,
        )
        product_state_pablo_library = _product_state_library_with_charge_templates(
            product_state_pablo_library,
            pablo_result.topology,
        )

    LOGGER.info("Parameterizing conjugate with OpenFF Interchange")
    parameterization_result = create_interchange_from_pablo_topology(
        pablo_result.topology,
        settings=settings.parameterization,
        charge_from_molecules=_formal_charge_templates_from_topology(pablo_result.topology),
    )
    if not parameterization_result.success or parameterization_result.interchange is None:
        raise RuntimeError("OpenFF Interchange parameterization did not produce an interchange")

    site_label = "site" if len(resolved_plans) == 1 else "sites"
    diagnostics = [
        f"Conjugate construction completed ({len(resolved_plans)} attachment {site_label})",
    ]
    smoke_result = None
    local_minimization_result = None
    run_combined_smoke_relaxation = False
    if run_product_state_local_minimization and len(resolved_plans) == 1:
        if _supports_product_state_local_minimization(specs[0]):
            product_state_requirement = product_state_requirements[0]
            local_settings = _local_minimization_settings_for_product(
                crosslinked_pdb_path,
                base_settings=local_minimization_settings or _default_local_minimization_settings(),
                requirement=resolved_plans[0].pablo_crosslink_requirement,
                product_state_pablo_library=product_state_pablo_library,
            )
            LOGGER.info("Running product-state local minimization")
            local_minimization_result = run_post_crosslink_local_minimization(
                crosslinked_pdb_path,
                artifact_dir,
                settings=local_settings,
                pablo_crosslink_requirement=product_state_requirement,
                product_state_pablo_library=product_state_pablo_library,
                resolved_plan=resolved_plans[0],
            )
            if not local_minimization_result.success:
                blocker = getattr(local_minimization_result, "blocker", None)
                relaxed_path = getattr(local_minimization_result, "relaxed_pdb_path", None)
                if blocker is not None or relaxed_path is None:
                    detail = blocker or "local minimization did not produce a relaxed PDB"
                    raise RuntimeError(f"Product-state local minimization failed: {detail}")
                diagnostics.append(
                    "Product-state local minimization completed and wrote a relaxed PDB, "
                    "but post-minimization geometry validation did not pass. Continuing with "
                    "the local-minimized artifact without falling back to smoke relaxation."
                )
        else:
            mechanism_name = _attachment_mechanism_name(specs[0])
            diagnostics.append(
                "Product-state local minimization was requested for mechanism "
                f"'{mechanism_name}', but the NHS-Lys local selector machinery does not "
                "support this chemistry; using one combined restrained vacuum "
                "smoke/minimization for the assembled product."
            )
            run_combined_smoke_relaxation = True
    elif run_product_state_local_minimization and len(resolved_plans) > 1:
        diagnostics.append(
            "Product-state local minimization was requested for multiple attachments; "
            "using one combined restrained vacuum smoke/minimization for the assembled product."
        )
        run_combined_smoke_relaxation = True

    if (
        local_minimization_result is None
        and smoke_result is None
        and (settings.run_smoke or run_combined_smoke_relaxation)
    ):
        LOGGER.info("Running combined restrained vacuum smoke relaxation")
        smoke_result = _run_restrained_vacuum_smoke_with_diagnostics(
            parameterization_result.interchange,
            artifact_dir,
            settings=settings.smoke,
            crosslinked_pdb_path=crosslinked_pdb_path,
            attachment_specs=specs,
        )
        if not smoke_result.success:
            raise RuntimeError("OpenMM restrained vacuum smoke did not report success")

    return (
        ConjugateConstructionResult(
            output_dir=artifact_dir,
            resolved_plan=resolved_plans[0],
            resolved_plans=resolved_plans,
            attachment_specs=specs,
            crosslink_validation=crosslink_validations[0],
            crosslink_validations=crosslink_validations,
            placement=placements[0],
            placements=placements,
            assembly=assembly_result,
            pablo=pablo_result,
            parameterization=parameterization_result,
            smoke=smoke_result,
            local_minimization=local_minimization_result,
            product_state_pablo_library=product_state_pablo_library,
            crosslinked_pdb_path=crosslinked_pdb_path,
            diagnostics=tuple(diagnostics),
        ),
        pablo_result.topology,
    )


def _run_restrained_vacuum_smoke_with_diagnostics(
    interchange: Any,
    artifact_dir: Path,
    *,
    settings: VacuumSmokeSettings,
    crosslinked_pdb_path: Path,
    attachment_specs: tuple[Any, ...],
) -> Any:
    """Run smoke while preserving compatibility with legacy test doubles."""
    signature = inspect.signature(run_restrained_vacuum_smoke)
    if "crosslinked_pdb_path" not in signature.parameters:
        return run_restrained_vacuum_smoke(interchange, artifact_dir, settings=settings)
    return run_restrained_vacuum_smoke(
        interchange,
        artifact_dir,
        settings=settings,
        crosslinked_pdb_path=crosslinked_pdb_path,
        attachment_specs=attachment_specs,
    )


def _supports_product_state_local_minimization(spec: Any) -> bool:
    """Return whether the attachment can use the NHS-Lys local minimizer."""
    return _attachment_mechanism_name(spec) in {
        _NHS_LYS_REACTION_NAME,
        _NHS_LYS_COORDINATE_BACKEND_MECHANISM,
    }


def _attachment_mechanism_name(spec: Any) -> str:
    """Resolve the best available mechanism name for diagnostics and capability checks."""
    plan = getattr(spec, "resolved_plan", None)
    contract = getattr(plan, "contract", None)
    mechanism_name = getattr(contract, "mechanism_name", None) or getattr(
        spec,
        "reaction_name",
        None,
    )
    if mechanism_name is None:
        return "unknown"
    return str(mechanism_name).strip().lower() or "unknown"


def _construct_multi_modifier_linked_protein(
    *,
    protein_pdb_path: Path | str,
    modifiers: tuple[GeneratedPolymerFragment, ...],
    resolved_plans: tuple[ResolvedAttachmentPlan, ...],
    ccd_pablo_policy: Any,
    output_dir: Path | str,
    chain_policy: Any | None,
    settings: ModifierConstructionSettings,
    use_product_state_pablo_library: bool,
    attachment_specs: tuple[AttachmentBuildSpec, ...] | None = None,
    source_moieties: tuple[GeneratedMoietyFragment, ...] | None = None,
    run_product_state_local_minimization: bool = False,
    local_minimization_settings: Any | None = None,
) -> tuple[Any, Any]:
    """Compatibility shim for legacy multi-modifier private callers."""
    if not modifiers or len(modifiers) != len(resolved_plans):
        raise ValueError("Multi-attachment construction requires aligned modifiers and plans")
    if attachment_specs is not None and len(attachment_specs) != len(resolved_plans):
        raise ValueError("Multi-attachment construction requires aligned specs and plans")
    if source_moieties is not None and len(source_moieties) != len(resolved_plans):
        raise ValueError("Multi-attachment construction requires aligned source moieties and plans")

    if attachment_specs is not None:
        specs = attachment_specs
    elif source_moieties is not None:
        specs = _product_state_specs_for_inputs(
            attachment_specs=None,
            source_moieties=source_moieties,
            resolved_plans=resolved_plans,
        )
        specs = tuple(
            spec.model_copy(update={"generated_fragment": modifier})
            for spec, modifier in zip(specs, modifiers, strict=True)
        )
    else:
        specs = tuple(
            SimpleNamespace(
                attachment_id=f"attachment_{index:02d}",
                attachment_index=index,
                generated_fragment=modifier,
                resolved_plan=plan,
                source_sidecars={},
                fragment=SimpleNamespace(source_kind="polymer"),
            )
            for index, (modifier, plan) in enumerate(
                zip(modifiers, resolved_plans, strict=True), start=1
            )
        )

    return _construct_conjugate_from_specs(
        protein_pdb_path=protein_pdb_path,
        specs=specs,
        ccd_pablo_policy=ccd_pablo_policy,
        output_dir=output_dir,
        chain_policy=chain_policy,
        settings=settings,
        use_product_state_pablo_library=use_product_state_pablo_library,
        run_product_state_local_minimization=run_product_state_local_minimization,
        local_minimization_settings=local_minimization_settings,
    )


def _product_state_specs_with_assembly_mappings(
    specs: tuple[Any, ...],
    *,
    assembly_result: Any,
) -> tuple[Any, ...]:
    """Return specs whose plans point at concrete product-PDB modifier residues."""
    mappings = getattr(assembly_result, "residue_mappings", {}) or {}
    updated_specs = []
    for fragment_index, spec in enumerate(specs, start=1):
        plan = spec.resolved_plan
        fragment_prefix = f"fragment_{fragment_index}:"
        fragment_mappings = {
            key.removeprefix(fragment_prefix): value
            for key, value in mappings.items()
            if key.startswith(fragment_prefix)
        }
        source_key = _modifier_source_residue_key(plan.modifier_link_atom)
        mapping = fragment_mappings.get(source_key)
        if mapping is None:
            updated_specs.append(_copy_spec_with_product_mappings(spec, fragment_mappings))
            continue

        modifier_link_atom = plan.modifier_link_atom.model_copy(
            update={
                "chain_id": str(mapping.get("target_chain", "C")),
                "residue_number": int(mapping["target_residue_number"]),
                "residue_name": plan.modifier_product_residue_name,
            }
        )
        updated_plan = plan.model_copy(update={"modifier_link_atom": modifier_link_atom})
        updated_specs.append(
            _copy_spec_with_product_mappings(
                spec,
                fragment_mappings,
                resolved_plan=updated_plan,
            )
        )
    return tuple(updated_specs)


def _copy_spec_with_product_mappings(
    spec: Any,
    product_residue_mappings: dict[str, dict[str, int | str]],
    **updates: Any,
) -> Any:
    """Copy a spec while attaching assembly residue mappings for Pablo templating."""
    updates["product_residue_mappings"] = product_residue_mappings
    copier = getattr(spec, "model_copy", None)
    if callable(copier):
        return copier(update=updates)
    return _namespace_copy(spec, **updates)


def _modifier_source_residue_key(atom: PdbAtomRecord) -> str:
    """Return the residue-mapping key suffix used by the PDB assembly writer."""
    return f"{atom.residue_number}{atom.insertion_code or ''}"


def _namespace_copy(obj: Any, **updates: Any) -> Any:
    """Copy a simple object namespace while replacing selected attributes."""
    data = getattr(obj, "__dict__", {}).copy()
    data.update(updates)
    return SimpleNamespace(**data)


def _policy_with_resolved_crosslinks(
    policy: Any,
    resolved_plans: tuple[ResolvedAttachmentPlan, ...],
) -> Any:
    """Return a Pablo policy containing product-state crosslinks for all plans."""
    updated = policy
    for plan in resolved_plans:
        updated = _policy_with_resolved_crosslink(updated, plan)
    return updated


def _product_state_specs_for_inputs(
    *,
    attachment_specs: tuple[AttachmentBuildSpec, ...] | None,
    source_moieties: tuple[GeneratedMoietyFragment, ...] | None,
    resolved_plans: tuple[ResolvedAttachmentPlan, ...],
) -> tuple[AttachmentBuildSpec, ...]:
    """Return resolved attachment specs for product-state library generation."""
    if attachment_specs is not None:
        return attachment_specs
    if source_moieties is None:
        raise ValueError(
            "Product-state Pablo library generation requires attachment specs or legacy "
            "source moieties"
        )
    return tuple(
        attachment_spec_from_moiety_plan(
            moiety,
            plan,
            attachment_config=SimpleNamespace(name=f"attachment_{index:02d}"),
            attachment_index=index,
            reaction_name=getattr(plan.contract, "mechanism_name", None) or "unknown",
        )
        for index, (moiety, plan) in enumerate(
            zip(source_moieties, resolved_plans, strict=True), start=1
        )
    )


def _product_state_pablo_library_for_specs(
    *,
    product_pdb: Path,
    source_protein_pdb: Path | str,
    specs: tuple[AttachmentBuildSpec, ...],
) -> Any:
    """Build product-state residue definitions for one or more attachment specs."""
    from polyzymd.builders.conjugation.pablo.product import (
        build_product_state_pablo_library_for_specs,
    )

    return build_product_state_pablo_library_for_specs(
        product_pdb=product_pdb,
        source_protein_pdb=source_protein_pdb,
        specs=specs,
    )


def _product_state_pablo_library_for_plans(
    *,
    product_pdb: Path,
    source_protein_pdb: Path | str,
    source_moieties: tuple[GeneratedMoietyFragment, ...],
    resolved_plans: tuple[ResolvedAttachmentPlan, ...],
) -> Any:
    """Build product-state definitions for legacy moiety/plan callers."""
    return _product_state_pablo_library_for_specs(
        product_pdb=product_pdb,
        source_protein_pdb=source_protein_pdb,
        specs=_product_state_specs_for_inputs(
            attachment_specs=None,
            source_moieties=source_moieties,
            resolved_plans=resolved_plans,
        ),
    )


def _safe_attachment_token(name: str) -> str:
    """Return a conservative artifact-directory token."""
    token = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in name.strip())
    return token.strip("_-") or "attachment"


def _save_direct_workflow_summary(
    result: Any,
    attachments: tuple[Any, ...],
    resolved_plans: list[ResolvedAttachmentPlan],
    path: Path,
) -> None:
    """Write a JSON summary for direct public conjugation requests."""
    payload = result.model_dump(mode="json")
    payload.update(
        {
            "crosslinked_pdb_path": str(result.construction.crosslinked_pdb_path),
            "attachments": [attachment.name for attachment in attachments],
            "resolved_plans": [plan.model_dump(mode="json") for plan in resolved_plans],
            "diagnostics": tuple(getattr(result.construction, "diagnostics", ())),
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _pablo_failure_message(result: Any) -> str:
    diagnostics = [f"{diag.code}: {diag.message}" for diag in result.diagnostics]
    joined = "; ".join(diagnostics) if diagnostics else "no diagnostics were reported"
    return f"Pablo ingestion failed or returned no topology for {result.path}: {joined}"


def _relaxed_conjugate_pdb(construction: ModifierConstructionResult) -> Path:
    local_minimization = getattr(construction, "local_minimization", None)
    if local_minimization is not None:
        relaxed_path = getattr(local_minimization, "relaxed_pdb_path", None)
        if relaxed_path is not None:
            return relaxed_path
    if construction.smoke is None:
        raise RuntimeError("protein-restrained vacuum smoke did not run")
    if construction.smoke.equilibrated_pdb_path is not None:
        return construction.smoke.equilibrated_pdb_path
    if construction.smoke.minimized_pdb_path is not None:
        return construction.smoke.minimized_pdb_path
    raise RuntimeError("protein-restrained vacuum smoke did not write a relaxed PDB")


def _formal_charge_templates_from_topology(topology: Any) -> tuple[Any, ...]:
    """Build smoke-only formal-charge templates for product-state parameterization."""
    molecules = tuple(getattr(topology, "molecules", ()) or ())
    return tuple(build_formal_charge_smoke_template(molecule) for molecule in molecules)


def _product_state_library_with_charge_templates(
    product_state_pablo_library: Any,
    topology: Any,
) -> Any:
    """Attach genuinely charged product-state templates from a Pablo topology."""
    templates = _charged_product_state_molecules_from_topology(
        topology,
        product_residue_names=_product_state_library_residue_names(product_state_pablo_library),
    )
    if hasattr(product_state_pablo_library, "model_copy"):
        return product_state_pablo_library.model_copy(update={"charge_templates": templates})
    product_state_pablo_library.charge_templates = templates
    return product_state_pablo_library


def _product_state_library_with_charge_bridge(
    product_state_pablo_library: Any,
    *,
    product_topology: Any,
    product_pdb: Path,
    source_protein_pdb: Path | str,
    specs: tuple[Any, ...],
    output_dir: Path,
    settings: InterchangeParameterizationSettings,
) -> Any:
    """Attach production partial-charge bridge records to the product library."""
    if not _supports_charge_bridge(specs, product_state_pablo_library):
        return product_state_pablo_library
    from polyzymd.builders.conjugation.pablo.charge_bridge import (
        build_product_state_charge_bridge,
    )

    return build_product_state_charge_bridge(
        product_state_pablo_library=product_state_pablo_library,
        product_topology=product_topology,
        product_pdb=product_pdb,
        source_protein_pdb=source_protein_pdb,
        specs=specs,
        output_dir=output_dir,
        settings=settings,
    )


def _supports_charge_bridge(
    product_specs: tuple[Any, ...], product_state_pablo_library: Any
) -> bool:
    """Return whether the local NAGL bridge supports the product-state library."""
    names = tuple(getattr(product_state_pablo_library, "residue_names", ()) or ())
    definitions = tuple(getattr(product_state_pablo_library, "definitions", ()) or ())
    if not names and not definitions:
        return False
    return all(_supports_product_state_local_minimization(spec) for spec in product_specs)


def _charged_product_state_molecules_from_topology(
    topology: Any,
    *,
    product_residue_names: tuple[str, ...],
) -> tuple[Any, ...]:
    """Return charged topology molecules containing product-state residues."""
    if not product_residue_names:
        return ()
    product_names = {name.upper() for name in product_residue_names}
    templates: list[Any] = []
    for molecule in tuple(getattr(topology, "molecules", ()) or ()):
        if getattr(molecule, "partial_charges", None) is None:
            continue
        if _molecule_has_product_residue(molecule, product_names):
            templates.append(molecule)
    return tuple(templates)


def _product_state_library_residue_names(product_state_pablo_library: Any) -> tuple[str, ...]:
    """Resolve product-state residue names from a generated Pablo library."""
    names = tuple(
        str(name) for name in getattr(product_state_pablo_library, "residue_names", ()) or ()
    )
    if names:
        return names
    definitions = tuple(getattr(product_state_pablo_library, "definitions", ()) or ())
    return tuple(
        str(getattr(definition, "residue_name", "")).strip()
        for definition in definitions
        if str(getattr(definition, "residue_name", "")).strip()
    )


def _molecule_has_product_residue(molecule: Any, product_names: set[str]) -> bool:
    """Return whether a molecule contains any product-state residue metadata."""
    for atom in tuple(getattr(molecule, "atoms", ()) or ()):
        residue_name = _metadata_residue_name(getattr(atom, "metadata", None))
        if residue_name in product_names:
            return True
    properties = getattr(molecule, "properties", None)
    residue_name = _metadata_residue_name(properties)
    return residue_name in product_names


def _metadata_residue_name(metadata: Any) -> str:
    """Return an uppercase residue name from atom or molecule metadata."""
    if metadata is None:
        return ""
    if isinstance(metadata, dict):
        value = metadata.get("residue_name") or metadata.get("residue")
    else:
        value = getattr(metadata, "residue_name", None) or getattr(metadata, "residue", None)
    return str(value or "").strip().upper()


def _local_minimization_settings_for_product(
    product_pdb_path: Path,
    *,
    base_settings: Any,
    requirement: PabloCrosslinkRequirement,
    product_state_pablo_library: Any | None = None,
) -> Any:
    """Build local minimization selectors from emitted product atom identities."""
    from polyzymd.builders.conjugation._relaxation import CrosslinkAtomSelector

    explicit_fields = set(getattr(base_settings, "model_fields_set", set()))
    if {"nz_selector", "c047_selector", "o020_selector"}.issubset(explicit_fields):
        return base_settings

    atoms = parse_pdb_atom_records(product_pdb_path)
    updates = {}

    if "nz_selector" not in explicit_fields:
        protein_atom = _unique_product_atom(
            atoms,
            residue_name=requirement.residues[0],
            atom_name=requirement.linking_atoms[0],
        )
        updates["nz_selector"] = _local_selector(protein_atom, CrosslinkAtomSelector)

    modifier_atom = None
    if "c047_selector" not in explicit_fields or "o020_selector" not in explicit_fields:
        modifier_atom = _unique_product_atom(
            atoms,
            residue_name=requirement.residues[1],
            atom_name=requirement.linking_atoms[1],
        )
    if "c047_selector" not in explicit_fields and modifier_atom is not None:
        updates["c047_selector"] = _local_selector(modifier_atom, CrosslinkAtomSelector)
    if "o020_selector" not in explicit_fields and modifier_atom is not None:
        modifier_oxygen_atom = _product_modifier_carbonyl_oxygen_atom(
            product_pdb_path,
            atoms,
            modifier_atom=modifier_atom,
            product_state_pablo_library=product_state_pablo_library,
        )
        updates["o020_selector"] = _local_selector(modifier_oxygen_atom, CrosslinkAtomSelector)

    if not updates:
        return base_settings

    return base_settings.model_copy(update=updates)


def _local_selector(atom: PdbAtomRecord, selector_cls: Any) -> Any:
    return selector_cls(
        serial=None,
        chain_id=atom.chain_id,
        residue_name=atom.residue_name,
        residue_number=atom.residue_number,
        atom_name=atom.atom_name,
    )


def _unique_product_atom(
    atoms: tuple[PdbAtomRecord, ...],
    *,
    residue_name: str,
    atom_name: str,
    same_residue_as: PdbAtomRecord | None = None,
) -> PdbAtomRecord:
    matches = [
        atom
        for atom in atoms
        if atom.residue_name.upper() == residue_name.upper()
        and atom.atom_name.upper() == atom_name.upper()
    ]
    if same_residue_as is not None:
        matches = [
            atom
            for atom in matches
            if atom.chain_id == same_residue_as.chain_id
            and atom.residue_number == same_residue_as.residue_number
            and atom.insertion_code == same_residue_as.insertion_code
        ]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one product atom {residue_name}:{atom_name}, found {len(matches)}"
        )
    return matches[0]


def _product_modifier_carbonyl_oxygen_atom(
    product_pdb_path: Path,
    atoms: tuple[PdbAtomRecord, ...],
    *,
    modifier_atom: PdbAtomRecord,
    product_state_pablo_library: Any | None,
) -> PdbAtomRecord:
    """Return the product-state oxygen bonded to the modifier linking carbon."""
    conect = _parse_product_conect_records(product_pdb_path)
    bonded_oxygen_atoms = _oxygen_atoms_bonded_to_modifier_atom(atoms, modifier_atom, conect)
    if len(bonded_oxygen_atoms) == 1:
        return bonded_oxygen_atoms[0]

    definition_oxygen_atoms = _modifier_oxygen_atoms_from_product_definitions(
        atoms,
        modifier_atom=modifier_atom,
        product_state_pablo_library=product_state_pablo_library,
    )
    if len(definition_oxygen_atoms) == 1:
        return definition_oxygen_atoms[0]

    candidates = bonded_oxygen_atoms or definition_oxygen_atoms
    candidate_text = (
        ", ".join(
            f"{atom.chain_id}:{atom.residue_name}:{atom.residue_number}:{atom.atom_name}"
            for atom in candidates
        )
        or "none"
    )
    raise RuntimeError(
        "Could not identify a unique product-state modifier carbonyl oxygen bonded to "
        f"{modifier_atom.residue_name}:{modifier_atom.atom_name} in {product_pdb_path}; "
        f"candidates: {candidate_text}. Product-state local minimization requires either "
        "PDB CONECT records or product-state Pablo residue definitions that identify exactly "
        "one oxygen atom in the same modifier residue bonded to the modifier linking/carbonyl "
        "carbon. Reactant leaving atoms are removed from product PDBs and are not valid "
        "selectors. Provide explicit LocalMinimizationSettings.o020_selector to override."
    )


def _oxygen_atoms_bonded_to_modifier_atom(
    atoms: tuple[PdbAtomRecord, ...],
    modifier_atom: PdbAtomRecord,
    conect: dict[int, set[int]],
) -> tuple[PdbAtomRecord, ...]:
    if modifier_atom.serial is None:
        return ()
    bonded_serials = conect.get(modifier_atom.serial, set())
    if not bonded_serials:
        return ()
    return tuple(
        atom
        for atom in atoms
        if atom.serial in bonded_serials
        and _same_product_residue(atom, modifier_atom)
        and _is_oxygen_atom(atom)
    )


def _modifier_oxygen_atoms_from_product_definitions(
    atoms: tuple[PdbAtomRecord, ...],
    *,
    modifier_atom: PdbAtomRecord,
    product_state_pablo_library: Any | None,
) -> tuple[PdbAtomRecord, ...]:
    definitions = tuple(getattr(product_state_pablo_library, "definitions", ()) or ())
    if not definitions:
        return ()
    residue_atoms = tuple(atom for atom in atoms if _same_product_residue(atom, modifier_atom))
    residue_atom_by_name = {atom.atom_name: atom for atom in residue_atoms}
    oxygen_names: set[str] = set()
    for definition in definitions:
        if (
            str(getattr(definition, "residue_name", "")).upper()
            != modifier_atom.residue_name.upper()
        ):
            continue
        atom_defs = {str(atom.name): atom for atom in getattr(definition, "atoms", ())}
        if modifier_atom.atom_name not in atom_defs:
            continue
        for bond in getattr(definition, "bonds", ()):
            atom1 = str(getattr(bond, "atom1", ""))
            atom2 = str(getattr(bond, "atom2", ""))
            if modifier_atom.atom_name not in (atom1, atom2):
                continue
            other_name = atom2 if atom1 == modifier_atom.atom_name else atom1
            atom_def = atom_defs.get(other_name)
            if atom_def is not None and str(getattr(atom_def, "symbol", "")).upper() == "O":
                oxygen_names.add(other_name)
    return tuple(
        atom
        for name, atom in residue_atom_by_name.items()
        if name in oxygen_names and _is_oxygen_atom(atom)
    )


def _parse_product_conect_records(path: Path) -> dict[int, set[int]]:
    conect: dict[int, set[int]] = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("CONECT"):
                continue
            try:
                source = int(line[6:11])
            except ValueError:
                continue
            targets = conect.setdefault(source, set())
            for start in range(11, len(line), 5):
                field = line[start : start + 5].strip()
                if not field:
                    continue
                try:
                    target = int(field)
                except ValueError:
                    continue
                targets.add(target)
                conect.setdefault(target, set()).add(source)
    return conect


def _same_product_residue(atom: PdbAtomRecord, anchor: PdbAtomRecord) -> bool:
    return (
        atom.chain_id == anchor.chain_id
        and atom.residue_name.upper() == anchor.residue_name.upper()
        and atom.residue_number == anchor.residue_number
        and atom.insertion_code == anchor.insertion_code
    )


def _is_oxygen_atom(atom: PdbAtomRecord) -> bool:
    return atom.element.strip().upper() == "O"


def _build_solvated_system(
    config: SimulationConfig,
    *,
    relaxed_conjugate_topology: Any,
    working_dir: Path,
    polymer_seed: int | None,
    create_interchange: bool,
    product_state_pablo_library: Any | None = None,
    parameterization_settings: InterchangeParameterizationSettings | None = None,
) -> SystemBuilder:
    """Build the solvated system around a relaxed conjugate topology."""
    builder = SystemBuilder.from_config(config)
    builder._working_dir = working_dir
    builder._enzyme_topology = relaxed_conjugate_topology
    builder._n_enzyme_molecules = 1
    builder._preserve_enzyme_chain_ids = True

    if config.substrate is not None:
        LOGGER.info("Building substrate for conjugated system")
        builder.build_substrate(
            sdf_path=config.substrate.sdf_path,
            conformer_index=config.substrate.conformer_index,
            charge_method=config.substrate.charge_method.value,
            residue_name=config.substrate.residue_name,
        )

    LOGGER.info("Combining conjugate, substrate, and free polymer solutes")
    builder.combine_solutes()
    if config.polymers is not None and config.polymers.enabled:
        LOGGER.info("Building and packing free polymers")
        _build_and_pack_free_polymers(builder, config, polymer_seed=polymer_seed)

    LOGGER.info("Solvating conjugated system")
    builder._solvent_builder.solvate_from_config(builder._combined_topology, config.solvent)
    builder._solvated_topology = builder._solvent_builder.solvated_topology

    if create_interchange:
        LOGGER.info("Creating final solvated OpenFF Interchange")
        create_final_conjugated_interchange(
            builder,
            product_state_pablo_library=product_state_pablo_library,
            settings=parameterization_settings,
        )

    return builder


def _build_direct_solvated_system(
    *,
    relaxed_conjugate_topology: Any,
    working_dir: Path,
    create_interchange: bool,
    product_state_pablo_library: Any | None = None,
    parameterization_settings: InterchangeParameterizationSettings | None = None,
) -> SystemBuilder:
    """Solvate a direct SMILES-moiety conjugate without a full SimulationConfig."""
    builder = SystemBuilder()
    builder._working_dir = working_dir
    builder._enzyme_topology = relaxed_conjugate_topology
    builder._n_enzyme_molecules = 1
    builder._preserve_enzyme_chain_ids = True
    LOGGER.info("Combining direct conjugate solutes")
    builder.combine_solutes()
    LOGGER.info("Solvating direct conjugated system")
    builder.solvate(padding=0.8, box_shape="cube")
    if create_interchange:
        LOGGER.info("Creating final solvated OpenFF Interchange")
        create_final_conjugated_interchange(
            builder,
            product_state_pablo_library=product_state_pablo_library,
            settings=parameterization_settings,
        )
    return builder


def _build_and_pack_free_polymers(
    builder: SystemBuilder,
    config: SimulationConfig,
    *,
    polymer_seed: int | None,
) -> None:
    polymers = config.polymers
    if polymers is None:
        return

    generation_mode = polymers.generation_mode.value
    monomer_smiles = None
    monomer_names = None
    residue_names = None
    reactions = None
    if generation_mode == "dynamic":
        monomer_smiles = {monomer.name: monomer.smiles for monomer in polymers.monomers}
        monomer_names = {monomer.label: monomer.name for monomer in polymers.monomers}
        residue_names = {
            monomer.name: monomer.residue_name
            for monomer in polymers.monomers
            if monomer.residue_name is not None
        }
        reactions = polymers.reactions

    effective_seed = polymers.random_seed if polymers.random_seed is not None else polymer_seed
    builder.build_polymers(
        characters=[monomer.label for monomer in polymers.monomers],
        probabilities=[monomer.probability for monomer in polymers.monomers],
        length=polymers.length,
        count=polymers.count,
        type_prefix=polymers.type_prefix,
        sdf_directory=polymers.sdf_directory,
        seed=effective_seed,
        generation_mode=generation_mode,
        monomer_smiles=monomer_smiles,
        monomer_names=monomer_names,
        residue_names=residue_names,
        reactions=reactions,
        charger_type=polymers.charger.value,
        max_retries=polymers.max_retries,
        cache_directory=polymers.cache_directory,
    )
    packing = polymers.packing
    builder.pack_polymers(
        padding=packing.padding,
        tolerance=packing.tolerance,
        movebadrandom=packing.movebadrandom,
        working_directory=builder._working_dir,
        box_vectors_nm=packing.box_vectors,
    )


def _pdb_coordinates_angstrom(path: Path) -> np.ndarray:
    coords: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(_ATOM_RECORD_PREFIXES):
                continue
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not coords:
        raise ValueError(f"No ATOM/HETATM coordinates found in {path}")
    return np.asarray(coords, dtype=float)


def _restore_smoke_pdb_atom_names(
    construction: ModifierConstructionResult,
    atom_name_template_pdb: Path | str,
) -> None:
    """Restore reference atom-name fields on smoke PDB artifacts when present."""
    if construction.smoke is None:
        return
    for pdb_path in (
        construction.smoke.minimized_pdb_path,
        construction.smoke.equilibrated_pdb_path,
    ):
        if pdb_path is not None:
            _restore_pdb_atom_name_fields(pdb_path, atom_name_template_pdb)


def _restore_pdb_atom_name_fields(
    target_pdb_path: Path | str,
    template_pdb_path: Path | str,
) -> int:
    """Copy same-order PDB atom-name columns from a template into a target PDB.

    Only the first ``N`` atom records are rewritten, where ``N`` is the number
    of atom records in the template. This lets a linked conjugate template fix
    the conjugate atom names in a larger solvated PDB without touching solvent.
    """
    target_path = Path(target_pdb_path)
    template_fields = _pdb_atom_name_fields(Path(template_pdb_path))
    if not template_fields:
        raise ValueError(f"No ATOM/HETATM atom names found in {template_pdb_path}")

    lines = target_path.read_text(encoding="utf-8").splitlines(keepends=True)
    restored = 0
    updated_lines: list[str] = []
    for line in lines:
        if line.startswith(_ATOM_RECORD_PREFIXES) and restored < len(template_fields):
            updated_lines.append(_replace_pdb_atom_name_field(line, template_fields[restored]))
            restored += 1
        else:
            updated_lines.append(line)

    if restored != len(template_fields):
        raise ValueError(
            "Target PDB has fewer atom records than the atom-name template: "
            f"target={restored}, template={len(template_fields)}"
        )

    target_path.write_text("".join(updated_lines), encoding="utf-8")
    return restored


def _apply_pdb_atom_names_to_topology(topology: Any, template_pdb_path: Path | str) -> None:
    """Set OpenFF atom names from a same-order PDB template."""
    template_atoms = _pdb_atom_records(Path(template_pdb_path))
    topology_atoms = list(_iter_openff_topology_atoms(topology))
    if len(topology_atoms) != len(template_atoms):
        raise ValueError(
            "Topology atom count does not match atom-name template count: "
            f"topology={len(topology_atoms)}, template={len(template_atoms)}"
        )
    for topology_atom, pdb_atom in zip(topology_atoms, template_atoms, strict=True):
        topology_atom.name = pdb_atom.atom_name
        metadata = getattr(topology_atom, "metadata", None)
        if isinstance(metadata, dict):
            metadata["atom_name"] = pdb_atom.atom_name


def _pdb_atom_records(path: Path) -> tuple[PdbAtomRecord, ...]:
    records: list[PdbAtomRecord] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(_ATOM_RECORD_PREFIXES):
                records.append(PdbAtomRecord.from_pdb_line(line, atom_index=len(records)))
    return tuple(records)


def _pdb_atom_name_fields(path: Path) -> tuple[str, ...]:
    fields: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(_ATOM_RECORD_PREFIXES):
                fields.append(line[12:16].ljust(4)[:4])
    return tuple(fields)


def _replace_pdb_atom_name_field(line: str, atom_name_field: str) -> str:
    newline = "\n" if line.endswith("\n") else ""
    body = line[:-1] if newline else line
    if body.endswith("\r"):
        body = body[:-1]
        newline = "\r" + newline
    padded = body.ljust(16)
    return f"{padded[:12]}{atom_name_field[:4]}{padded[16:]}{newline}"


def _iter_openff_topology_atoms(topology: Any) -> tuple[Any, ...]:
    if hasattr(topology, "atoms"):
        atoms_attr = topology.atoms
        atoms = atoms_attr() if callable(atoms_attr) else atoms_attr
        return tuple(atoms)
    atoms: list[Any] = []
    for molecule in getattr(topology, "molecules", ()):
        atoms.extend(getattr(molecule, "atoms", ()))
    return tuple(atoms)


def _required_int(value: int | None, name: str) -> int:
    if value is None:
        raise ValueError(f"{name} is required")
    return int(value)


def _optional_str(value: object) -> str | None:
    text = "" if value is None else str(value).strip()
    return text or None
