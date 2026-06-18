"""Config-driven protein-polymer conjugate packing workflow."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import numpy as np
from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.construction import (
    ModifierConstructionResult,
    ModifierConstructionSettings,
)
from polyzymd.builders.conjugation.contracts import (
    ResolvedAttachmentPlan,
    placed_fragment_from_resolved_plan,
)
from polyzymd.builders.conjugation.crosslinks import (
    require_pablo_crosslink_requirement,
)
from polyzymd.builders.conjugation.direct_openff import build_direct_openff_linkage
from polyzymd.builders.conjugation.linkers import NhsLysModifierLinker
from polyzymd.builders.conjugation.pablo_adapter import PabloIngestor
from polyzymd.builders.conjugation.parameterization import (
    InterchangeParameterizationSettings,
    create_interchange_from_openff_topology,
    create_interchange_from_pablo_topology,
)
from polyzymd.builders.conjugation.pdb_assembly import (
    CrosslinkedPdbAssemblyOptions,
    PdbAtomRecord,
    write_crosslinked_pdb,
)
from polyzymd.builders.conjugation.placement import (
    PackmolModifierPlacementSettings,
    place_modifier_with_packmol,
)
from polyzymd.builders.conjugation.polymer_fragment import GeneratedPolymerFragment
from polyzymd.builders.conjugation.polymer_recipe import (
    PolymeristGenerationSmokeResult,
    PolymerRecipe,
    generate_polymerist_smoke_polymer,
)
from polyzymd.builders.conjugation.polymerist_pdb import generated_fragment_from_polymerist_pdb
from polyzymd.builders.conjugation.reaction_roles import (
    STRUCTURE_MATCHING_BLOCKER_MESSAGE,
    atom_mapped_reaction_from_mechanism_config,
    resolve_reaction_roles_from_identity_map,
)
from polyzymd.builders.conjugation.smoke import VacuumSmokeSettings, run_restrained_vacuum_smoke
from polyzymd.builders.system_builder import SystemBuilder
from polyzymd.config.schema import (
    ConjugationCcdCrosslinkConfig,
    ConjugationConfig,
    ConjugationMode,
    SimulationConfig,
)

_ATOM_RECORD_PREFIXES = ("ATOM", "HETATM")
_NHS_LYS_COORDINATE_BACKEND_MECHANISM = "nhs_lys_amide"


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
    allow_direct_openff_fallback: bool = False
    preserve_reference_atom_names: bool = True
    placement: PackmolModifierPlacementSettings = Field(
        default_factory=PackmolModifierPlacementSettings
    )
    conjugate_parameterization: InterchangeParameterizationSettings = Field(
        default_factory=InterchangeParameterizationSettings
    )
    vacuum_smoke: VacuumSmokeSettings = Field(default_factory=lambda: _protein_restrained_smoke())


class ConjugatedPolymerSystemResult(BaseModel):
    """Artifacts from the complete conjugate plus free-polymer solvation workflow."""

    model_config = {"arbitrary_types_allowed": True}

    output_dir: Path
    generated_sequence: str
    reactive_sequence_index: int
    reactive_residue_selector: dict[str, int | str]
    conjugate_generation: PolymeristGenerationSmokeResult
    construction: ModifierConstructionResult
    relaxed_conjugate_pdb_path: Path | None = None
    solvated_pdb_path: Path | None = None
    workflow_json_path: Path | None = None
    final_interchange_created: bool = False
    modifier: GeneratedPolymerFragment = Field(exclude=True)
    relaxed_conjugate_topology: Any = Field(default=None, exclude=True)
    solvated_topology: Any = Field(default=None, exclude=True)
    final_interchange: Any | None = Field(default=None, exclude=True)
    system_builder: SystemBuilder | None = Field(default=None, exclude=True)

    def save(self, path: Path | str) -> Path:
        """Write a JSON sidecar excluding heavy topology and Interchange objects."""
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(self.model_dump(mode="json"), indent=2) + "\n")
        return output_path


def build_conjugated_polymer_system_from_config_path(
    config_path: Path | str,
    *,
    output_dir: Path | str | None = None,
    settings: ConjugatedPolymerSystemSettings | None = None,
    free_polymer_seed: int | None = None,
) -> ConjugatedPolymerSystemResult:
    """Load a config YAML and build the relaxed, solvated conjugate system."""
    from polyzymd.config.loader import load_config

    path = Path(config_path)
    config = load_config(path)
    effective_output_dir = output_dir or path.parent / "artifacts" / path.stem
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
) -> ConjugatedPolymerSystemResult:
    """Build a relaxed protein-polymer conjugate, pack free polymers, and solvate.

    The workflow consumes the existing PolyzyMD config schema: ``conjugation``
    defines one covalent protein-polymer attachment, ``polymers`` defines free
    non-covalent polymer chains to pack around that conjugate, and ``solvent``
    defines the final solvent box.
    """
    workflow_settings = settings or ConjugatedPolymerSystemSettings()
    artifact_dir = Path(output_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)

    attachment = _single_enabled_attachment(config.conjugation)
    _require_supported_coordinate_backend(attachment)
    recipe = _polymer_recipe_from_attachment(attachment)
    reactive_sequence_index = _reactive_sequence_index(recipe)

    generation = generate_polymerist_smoke_polymer(
        recipe,
        artifact_dir / workflow_settings.conjugate_cache_dir_name,
        force_regenerate=workflow_settings.force_regenerate_conjugate_polymer,
        max_retries=workflow_settings.conjugate_polymerist_max_retries,
        energy_minimize=workflow_settings.conjugate_polymerist_energy_minimize,
    )
    if generation.pdb_path is None:
        raise RuntimeError("Polymerist did not produce a conjugate-polymer PDB")

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
    resolved_plan = linker.resolve_plan(config.enzyme.pdb_path, modifier)
    ccd_pablo_policy = _policy_with_resolved_crosslink(
        config.conjugation.ccd_pablo,
        resolved_plan,
    )

    construction_settings = ModifierConstructionSettings(
        placement=workflow_settings.placement,
        parameterization=workflow_settings.conjugate_parameterization,
        smoke=workflow_settings.vacuum_smoke,
        run_smoke=True,
    )
    construction, construction_topology = _construct_nhs_lys_modifier_linked_protein(
        protein_pdb_path=config.enzyme.pdb_path,
        modifier=modifier,
        linker=linker,
        resolved_plan=resolved_plan,
        ccd_pablo_policy=ccd_pablo_policy,
        chain_policy=config.conjugation.chain_policy,
        output_dir=artifact_dir / workflow_settings.conjugate_artifact_dir_name,
        settings=construction_settings,
        allow_direct_openff_fallback=workflow_settings.allow_direct_openff_fallback,
    )
    if workflow_settings.preserve_reference_atom_names:
        _restore_smoke_pdb_atom_names(construction, construction.crosslinked_pdb_path)

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
    )
    solvated_pdb_path = artifact_dir / workflow_settings.solvated_pdb_name
    builder.save_topology(solvated_pdb_path)
    if workflow_settings.preserve_reference_atom_names:
        _restore_pdb_atom_name_fields(solvated_pdb_path, construction.crosslinked_pdb_path)

    result = ConjugatedPolymerSystemResult(
        output_dir=artifact_dir,
        generated_sequence=generation.sequence,
        reactive_sequence_index=reactive_sequence_index,
        reactive_residue_selector=reactive_selector,
        conjugate_generation=generation,
        construction=construction,
        relaxed_conjugate_pdb_path=relaxed_pdb,
        solvated_pdb_path=solvated_pdb_path,
        final_interchange_created=builder.interchange is not None,
        modifier=modifier,
        relaxed_conjugate_topology=relaxed_topology,
        solvated_topology=builder.solvated_topology,
        final_interchange=builder.interchange,
        system_builder=builder,
    )
    result.workflow_json_path = result.save(artifact_dir / workflow_settings.workflow_json_name)
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


def _single_enabled_attachment(conjugation: ConjugationConfig | None) -> Any:
    if conjugation is None or not conjugation.enabled:
        raise ValueError("conjugation.enabled must be true for this workflow")
    if conjugation.mode != ConjugationMode.CONSTRUCT:
        raise ValueError("this workflow requires conjugation.mode: construct")
    attachments = [attachment for attachment in conjugation.attachments if attachment.enabled]
    if len(attachments) != 1:
        raise ValueError("v1 conjugated polymer workflow requires exactly one enabled attachment")
    return attachments[0]


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
    site = attachment.site
    product_residues = attachment.mechanism.product_residues
    return NhsLysModifierLinker(
        target_chain=site.chain_id,
        target_residue_name=site.residue_name or "LYS",
        target_residue_number=_required_int(site.residue_number, "site.residue_number"),
        target_insertion_code=site.insertion_code,
        target_atom_name=site.atom_name or "NZ",
        lysine_target_resname=product_residues.site or "LYX",
        modifier_target_resname=product_residues.moiety or "NHX",
    )


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
    modifier: GeneratedPolymerFragment,
    linker: NhsLysModifierLinker,
    resolved_plan: ResolvedAttachmentPlan,
    ccd_pablo_policy: Any,
    output_dir: Path | str,
    chain_policy: Any | None,
    settings: ModifierConstructionSettings,
    allow_direct_openff_fallback: bool,
) -> tuple[ModifierConstructionResult, Any]:
    """Run the NHS-Lys construction path using Pablo product-residue crosslinks."""
    artifact_dir = Path(output_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)

    product_state_requirement = _product_state_crosslink_requirement(resolved_plan)
    crosslink_validation = require_pablo_crosslink_requirement(
        ccd_pablo_policy,
        product_state_requirement,
    )

    placement_result = place_modifier_with_packmol(
        protein_pdb_path,
        modifier,
        linker,
        artifact_dir,
        settings=settings.placement,
    )
    placed_modifier = placed_fragment_from_resolved_plan(
        placement_result.placed_modifier,
        resolved_plan,
    )
    crosslinked_pdb_path = artifact_dir / settings.crosslinked_pdb_name
    assembly_result = write_crosslinked_pdb(
        protein_pdb_path,
        placed_modifier,
        resolved_plan.to_nhs_lys_pdb_attachment(),
        crosslinked_pdb_path,
        CrosslinkedPdbAssemblyOptions(),
    )

    pablo_result = PabloIngestor(policy=ccd_pablo_policy).ingest_structure(
        crosslinked_pdb_path,
        chain_policy=chain_policy,
        output_dir=artifact_dir,
    )
    if not pablo_result.success or pablo_result.topology is None:
        if not allow_direct_openff_fallback:
            raise RuntimeError(_pablo_failure_message(pablo_result))
        return _construct_with_direct_openff_fallback(
            protein_pdb_path=protein_pdb_path,
            placed_modifier=placed_modifier,
            resolved_plan=resolved_plan,
            output_dir=artifact_dir,
            pablo_result=pablo_result,
            crosslink_validation=crosslink_validation,
            placement_result=placement_result,
            assembly_result=assembly_result,
            settings=settings,
        )

    parameterization_result = create_interchange_from_pablo_topology(
        pablo_result.topology,
        settings=settings.parameterization,
    )
    if not parameterization_result.success or parameterization_result.interchange is None:
        raise RuntimeError("OpenFF Interchange parameterization did not produce an interchange")

    smoke_result = None
    if settings.run_smoke:
        smoke_result = run_restrained_vacuum_smoke(
            parameterization_result.interchange,
            artifact_dir,
            settings=settings.smoke,
        )
        if not smoke_result.success:
            raise RuntimeError("OpenMM restrained vacuum smoke did not report success")

    return (
        ModifierConstructionResult(
            output_dir=artifact_dir,
            resolved_plan=resolved_plan,
            crosslink_validation=crosslink_validation,
            placement=placement_result,
            assembly=assembly_result,
            pablo=pablo_result,
            parameterization=parameterization_result,
            smoke=smoke_result,
            crosslinked_pdb_path=crosslinked_pdb_path,
            diagnostics=("NHS-Lys modifier-linked protein construction completed",),
        ),
        pablo_result.topology,
    )


def _construct_with_direct_openff_fallback(
    *,
    protein_pdb_path: Path | str,
    placed_modifier: Any,
    resolved_plan: Any,
    output_dir: Path,
    pablo_result: Any,
    crosslink_validation: Any,
    placement_result: Any,
    assembly_result: Any,
    settings: ModifierConstructionSettings,
) -> tuple[ModifierConstructionResult, Any]:
    """Parameterize and relax the linked conjugate without Pablo topology ingestion."""
    direct_result = build_direct_openff_linkage(
        protein_pdb_path=protein_pdb_path,
        modifier=placed_modifier,
        resolved_plan=resolved_plan,
        output_dir=output_dir / "direct-openff-fallback",
        build_openff_topology=True,
    )
    if direct_result.topology is None or direct_result.charge_template is None:
        raise RuntimeError(
            _pablo_failure_message(pablo_result)
            + "; direct OpenFF fallback could not build a connected topology"
        )

    parameterization_result = create_interchange_from_openff_topology(
        direct_result.topology,
        settings=settings.parameterization,
        charge_from_molecules=[direct_result.charge_template],
        require_charge_templates=True,
        success_diagnostic="OpenFF Interchange was created from the direct conjugate topology",
        failure_subject="direct OpenFF conjugate topology",
    )
    if not parameterization_result.success or parameterization_result.interchange is None:
        raise RuntimeError("Direct OpenFF fallback did not produce an interchange")

    smoke_result = None
    if settings.run_smoke:
        smoke_result = run_restrained_vacuum_smoke(
            parameterization_result.interchange,
            output_dir,
            settings=settings.smoke,
        )
        if not smoke_result.success:
            raise RuntimeError("OpenMM restrained vacuum smoke did not report success")

    return (
        ModifierConstructionResult(
            output_dir=output_dir,
            resolved_plan=resolved_plan,
            crosslink_validation=crosslink_validation,
            placement=placement_result,
            assembly=assembly_result,
            pablo=pablo_result,
            parameterization=parameterization_result,
            smoke=smoke_result,
            crosslinked_pdb_path=direct_result.linked_pdb_path,
            diagnostics=(
                "Pablo ingestion failed; direct OpenFF fallback produced the relaxed conjugate",
                *direct_result.limitations,
            ),
        ),
        direct_result.topology,
    )


def _pablo_failure_message(result: Any) -> str:
    diagnostics = [f"{diag.code}: {diag.message}" for diag in result.diagnostics]
    joined = "; ".join(diagnostics) if diagnostics else "no diagnostics were reported"
    return f"Pablo ingestion failed or returned no topology for {result.path}: {joined}"


def _relaxed_conjugate_pdb(construction: ModifierConstructionResult) -> Path:
    if construction.smoke is None:
        raise RuntimeError("protein-restrained vacuum smoke did not run")
    if construction.smoke.equilibrated_pdb_path is not None:
        return construction.smoke.equilibrated_pdb_path
    if construction.smoke.minimized_pdb_path is not None:
        return construction.smoke.minimized_pdb_path
    raise RuntimeError("protein-restrained vacuum smoke did not write a relaxed PDB")


def _build_solvated_system(
    config: SimulationConfig,
    *,
    relaxed_conjugate_topology: Any,
    working_dir: Path,
    polymer_seed: int | None,
    create_interchange: bool,
) -> SystemBuilder:
    builder = SystemBuilder.from_config(config)
    builder._working_dir = working_dir
    builder._enzyme_topology = relaxed_conjugate_topology
    builder._n_enzyme_molecules = 1
    builder._preserve_enzyme_chain_ids = True

    if config.substrate is not None:
        builder.build_substrate(
            sdf_path=config.substrate.sdf_path,
            conformer_index=config.substrate.conformer_index,
            charge_method=config.substrate.charge_method.value,
            residue_name=config.substrate.residue_name,
        )

    builder.combine_solutes()
    if config.polymers is not None and config.polymers.enabled:
        _build_and_pack_free_polymers(builder, config, polymer_seed=polymer_seed)

    builder._solvent_builder.solvate_from_config(builder._combined_topology, config.solvent)
    builder._solvated_topology = builder._solvent_builder.solvated_topology

    if create_interchange:
        use_optimized = config.polymers is not None and config.polymers.enabled
        builder.create_interchange(use_optimized_combining=use_optimized)

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
