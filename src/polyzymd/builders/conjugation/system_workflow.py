"""Config-driven protein-polymer conjugate packing workflow."""

from __future__ import annotations

import copy
import json
import logging
from collections.abc import Iterable, Mapping, MutableMapping
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Literal

import numpy as np
from pydantic import BaseModel, Field

from polyzymd.builders.conjugation._linkage import (
    ReactionProduct,
    require_pablo_crosslink_requirement,
)
from polyzymd.builders.conjugation._moiety_provider import (
    attachment_uses_pdb_fragment,
    prepare_resolved_moiety_source,
    resolve_moiety_source,
    validate_moiety_source_config,
)
from polyzymd.builders.conjugation.construction import (
    ModifierConstructionResult,
    ModifierConstructionSettings,
)
from polyzymd.builders.conjugation.final_interchange import create_final_conjugated_interchange
from polyzymd.builders.conjugation.force_fields import resolve_conjugate_force_fields
from polyzymd.builders.conjugation.models import ConjugationResult
from polyzymd.builders.conjugation.native_openmm_glycam import (
    create_native_openmm_glycam_handoff,
    native_glycam_enabled,
)
from polyzymd.builders.conjugation.pablo.charge_templates import (
    build_conjugate_charge_templates,
)
from polyzymd.builders.conjugation.pablo.ingestion import PabloIngestor
from polyzymd.builders.conjugation.pablo.parameterization import (
    InterchangeParameterizationSettings,
    create_interchange_from_pablo_topology,
)
from polyzymd.builders.conjugation.pablo.product_state import (
    metadata_residue_name,
    molecule_contains_product_residue,
    product_residue_names,
    product_state_library_has_provenance,
)
from polyzymd.builders.conjugation.placement import (
    PackmolModifierPlacementSettings,
    PackmolOutputValidationError,
    place_modifier_with_resolved_plan,
    place_modifiers_with_resolved_plans,
)
from polyzymd.builders.conjugation.polymer import MultiResidueGenerationResult
from polyzymd.builders.conjugation.reactions._roles import (
    STRUCTURE_MATCHING_BLOCKER_MESSAGE,
    atom_mapped_reaction_from_mechanism_config,
    resolve_reaction_roles_from_identity_map,
)
from polyzymd.builders.conjugation.reactions.library import get_reaction
from polyzymd.builders.conjugation.relaxation import (
    ConjugateRelaxationSettings,
    relax_conjugate,
)
from polyzymd.builders.conjugation.structure.parsing import (
    ATOM_RECORD_PREFIXES,
    parse_pdb_conect_pairs,
    pdb_coordinates,
)
from polyzymd.builders.conjugation.structure.parsing import (
    parse_pdb_atom_records as parse_structure_pdb_atom_records,
)
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
from polyzymd.builders.conjugation.validation import (
    ValidationStatus,
    build_conjugate_validation_report,
    summarize_nonbonded_heavy_clashes,
    validate_conect_graph,
)
from polyzymd.builders.system_builder import SystemBuilder
from polyzymd.config.schema import (
    BuildScope,
    ConjugationCcdCrosslinkConfig,
    ConjugationCcdPabloPolicyConfig,
    ConjugationChainPolicyConfig,
    ConjugationConfig,
    SimulationConfig,
)
from polyzymd.exporters.exact_openmm import create_exact_export_bundle

_ATOM_RECORD_PREFIXES = ATOM_RECORD_PREFIXES
_NHS_LYS_REACTION = get_reaction("nhs_lys")
_NHS_LYS_COORDINATE_BACKEND_MECHANISM = _NHS_LYS_REACTION.coordinate_backend_mechanism
_PDB_FRAGMENT_COORDINATE_ONLY_MIN_PACKMOL_TOLERANCE_ANGSTROM = 2.5
_PDB_FRAGMENT_COORDINATE_ONLY_MAX_PACKMOL_ATTEMPTS = 20
_PDB_FRAGMENT_COORDINATE_ONLY_MIN_HEAVY_NONBONDED_ANGSTROM = 2.0
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
    run_relaxation: bool = True
    direct_solvation_padding: float = Field(0.8, gt=0.0)
    direct_solvation_box_shape: str = "cube"
    protein_canonicalization: ProteinCanonicalizationSettings = Field(
        default_factory=ProteinCanonicalizationSettings
    )
    placement: PackmolModifierPlacementSettings = Field(
        default_factory=PackmolModifierPlacementSettings
    )
    conjugate_parameterization: InterchangeParameterizationSettings = Field(
        default_factory=InterchangeParameterizationSettings
    )
    relaxation: ConjugateRelaxationSettings = Field(default_factory=ConjugateRelaxationSettings)
    pdb_fragment_output_mode: Literal["coordinate_only", "experimental_pablo"] = "coordinate_only"
    build_scope: BuildScope = BuildScope.SYSTEM


def _settings_with_config_defaults(
    settings: ConjugatedPolymerSystemSettings | None,
    config: Any,
) -> ConjugatedPolymerSystemSettings:
    """Return workflow settings with missing runtime defaults copied from config."""
    workflow_settings = settings or ConjugatedPolymerSystemSettings()
    if settings is None:
        conjugation_placement = getattr(getattr(config, "conjugation", None), "placement", None)
        if conjugation_placement is not None:
            placement = workflow_settings.placement.model_copy(
                update={"timeout_seconds": conjugation_placement.timeout_seconds}
            )
            workflow_settings = workflow_settings.model_copy(update={"placement": placement})
    platform_name = getattr(getattr(config, "openmm", None), "platform", None)
    if platform_name is None or workflow_settings.relaxation.platform_name is not None:
        return workflow_settings
    relaxation = workflow_settings.relaxation.model_copy(
        update={"platform_name": str(platform_name)}
    )
    return workflow_settings.model_copy(update={"relaxation": relaxation})


class ConjugateConstructionResult(BaseModel):
    """Specs-first construction result for product-state conjugation."""

    model_config = {"arbitrary_types_allowed": True}

    output_dir: Path
    reaction_product: ReactionProduct
    reaction_products: tuple[ReactionProduct, ...]
    crosslink_validation: Any
    crosslink_validations: tuple[Any, ...]
    placement: Any
    placements: tuple[Any, ...]
    assembly: Any
    pablo: Any
    parameterization: Any
    relaxation: Any | None = None
    product_state_pablo_library: Any | None = Field(default=None, exclude=True)
    crosslinked_pdb_path: Path
    validation_report_path: Path | None = None
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
    workflow_settings = _settings_with_config_defaults(settings, config)
    artifact_dir = Path(output_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    construction_dir = artifact_dir / workflow_settings.conjugate_artifact_dir_name

    attachments = _enabled_supported_attachments(config.conjugation)
    LOGGER.info("Enabled conjugation attachment count: %d", len(attachments))
    _log_attachment_additions(attachments)
    LOGGER.info("Preparing and canonicalizing source protein")
    protein_pdb_path, protein_canonicalization = _prepared_protein_pdb_path(
        config.enzyme.pdb_path,
        output_dir=construction_dir,
        settings=workflow_settings,
    )
    spec_payloads = tuple(
        _build_attachment_spec(
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
    modifiers = tuple(spec.fragment for spec in specs)
    resolved_plans = tuple(specs)
    if workflow_settings.build_scope is BuildScope.STRUCTURE:
        result = _build_pdb_fragment_coordinate_only_result(
            protein_pdb_path=protein_pdb_path,
            specs=specs,
            output_dir=artifact_dir,
            construction_dir=construction_dir,
            protein_canonicalization=protein_canonicalization,
            placement_settings=workflow_settings.placement,
            output_name="assembled_crosslinked.pdb",
            status="structure",
        )
        workflow_path = artifact_dir / workflow_settings.workflow_json_name
        result.workflow_json_path = workflow_path
        result.artifact_paths["workflow_json"] = workflow_path
        result.save(workflow_path)
        LOGGER.info("Saved structure-scope conjugation workflow JSON to %s", workflow_path)
        return result
    resolved_force_fields = resolve_conjugate_force_fields(config)
    use_native_glycam = native_glycam_enabled(config)
    use_mixed_overlay = resolved_force_fields.route == "mixed_overlay"
    if workflow_settings.build_scope is BuildScope.SOLUTE and use_mixed_overlay:
        raise NotImplementedError(
            "build scope 'solute' fails closed for mixed-overlay parameterization because "
            "isolated scope cannot yet be propagated without full-system recursion"
        )
    if (
        _uses_pdb_fragment_sources(attachments)
        and workflow_settings.pdb_fragment_output_mode == "coordinate_only"
        and not use_native_glycam
        and resolved_force_fields.route != "mixed_overlay"
    ):
        result = _build_pdb_fragment_coordinate_only_result(
            protein_pdb_path=protein_pdb_path,
            specs=specs,
            output_dir=artifact_dir,
            construction_dir=construction_dir,
            protein_canonicalization=protein_canonicalization,
            placement_settings=workflow_settings.placement,
        )
        workflow_path = artifact_dir / workflow_settings.workflow_json_name
        result.workflow_json_path = workflow_path
        result.artifact_paths["workflow_json"] = workflow_path
        result.save(workflow_path)
        LOGGER.info(
            "Saved coordinate-only PDB-fragment conjugation workflow JSON to %s", workflow_path
        )
        return result
    if _uses_pdb_fragment_sources(attachments):
        pdb_fragment_coordinate_artifact_path, _, _ = _write_pdb_fragment_coordinate_artifacts(
            protein_pdb_path=protein_pdb_path,
            specs=specs,
            construction_dir=construction_dir,
            placement_settings=workflow_settings.placement,
        )
        LOGGER.warning(
            "Continuing PDB-fragment input into experimental Pablo/OpenFF mode; "
            "coordinate-only artifact path is %s",
            pdb_fragment_coordinate_artifact_path,
        )
    else:
        pdb_fragment_coordinate_artifact_path = None
    ccd_pablo_policy = _policy_with_resolved_crosslinks(
        config.conjugation.ccd_pablo,
        resolved_plans,
    )

    construction_settings = ModifierConstructionSettings(
        placement=workflow_settings.placement,
        parameterization=workflow_settings.conjugate_parameterization,
        relaxation=workflow_settings.relaxation,
        run_relaxation=workflow_settings.run_relaxation,
    )
    try:
        construction, construction_topology = _construct_conjugate_from_specs(
            protein_pdb_path=protein_pdb_path,
            specs=specs,
            ccd_pablo_policy=ccd_pablo_policy,
            output_dir=construction_dir,
            chain_policy=config.conjugation.chain_policy,
            settings=construction_settings,
            use_product_state_pablo_library=workflow_settings.use_product_state_pablo_library,
            use_conjugate_relaxation=workflow_settings.run_relaxation,
        )
    except RuntimeError as exc:
        if pdb_fragment_coordinate_artifact_path is not None:
            raise RuntimeError(
                "Experimental PDB-fragment Pablo/OpenFF continuation failed after coordinate-only "
                f"artifact was written to {pdb_fragment_coordinate_artifact_path}"
            ) from exc
        raise
    if workflow_settings.preserve_reference_atom_names:
        _restore_relaxed_pdb_atom_names(construction, construction.crosslinked_pdb_path)

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

    if workflow_settings.build_scope is BuildScope.SOLUTE:
        builder = _build_isolated_conjugate_system(
            config,
            relaxed_conjugate_topology=relaxed_topology,
            working_dir=artifact_dir,
            product_state_pablo_library=getattr(construction, "product_state_pablo_library", None),
            parameterization_settings=workflow_settings.conjugate_parameterization,
            create_interchange=not use_native_glycam,
        )
        solute_pdb_path = artifact_dir / "solute.pdb"
        exact_export_bundle = None
        if use_native_glycam:
            exact_export_bundle = create_native_openmm_glycam_handoff(
                builder,
                config=config,
                construction=construction,
                output_dir=artifact_dir,
                solute_scope=True,
            )
            exact_export_bundle.write_pdb(solute_pdb_path)
        else:
            builder.save_topology(solute_pdb_path)
        result = ConjugationResult(
            status="solute",
            output_dir=artifact_dir,
            construction=construction,
            protein_canonicalization=protein_canonicalization,
            relaxed_conjugate_pdb_path=relaxed_pdb,
            solvated_pdb_path=solute_pdb_path,
            final_interchange_created=True,
            modifier=modifiers[0],
            modifiers=modifiers,
            attachment_specs=specs,
            relaxed_conjugate_topology=relaxed_topology,
            solvated_topology=builder.solvated_topology,
            final_interchange=None if exact_export_bundle is not None else builder.interchange,
            exact_export_bundle=exact_export_bundle,
            system_builder=builder,
        )
        if exact_export_bundle is not None:
            result.artifact_paths["native_openmm_glycam_audit"] = exact_export_bundle.audit_path
            result.artifact_paths["exact_openmm_exceptions"] = exact_export_bundle.sidecar_path
        return result

    builder = _build_solvated_system(
        config,
        relaxed_conjugate_topology=relaxed_topology,
        working_dir=artifact_dir,
        polymer_seed=free_polymer_seed,
        create_interchange=workflow_settings.create_final_interchange and not use_native_glycam,
        product_state_pablo_library=getattr(construction, "product_state_pablo_library", None),
        parameterization_settings=workflow_settings.conjugate_parameterization,
    )
    solvated_pdb_path = artifact_dir / workflow_settings.solvated_pdb_name
    exact_export_bundle = None
    overlay_artifact_paths: dict[str, Path] = {}
    if use_native_glycam:
        LOGGER.info("Creating native OpenMM GLYCAM handoff")
        exact_export_bundle = create_native_openmm_glycam_handoff(
            builder,
            config=config,
            construction=construction,
            output_dir=artifact_dir,
        )
        exact_export_bundle.write_pdb(solvated_pdb_path)
    elif use_mixed_overlay and workflow_settings.create_final_interchange:
        LOGGER.info("Creating mixed GLYCAM/OpenFF overlay exact handoff")
        exact_export_bundle, overlay_artifact_paths = _create_mixed_overlay_exact_handoff(
            builder,
            config=config,
            construction=construction,
            output_dir=artifact_dir,
            workflow_settings=workflow_settings,
        )
        exact_export_bundle.write_pdb(solvated_pdb_path)
    else:
        builder.save_topology(solvated_pdb_path)
    LOGGER.info("Wrote final solvated conjugate PDB to %s", solvated_pdb_path)
    if workflow_settings.preserve_reference_atom_names and exact_export_bundle is None:
        _restore_pdb_atom_name_fields(solvated_pdb_path, construction.crosslinked_pdb_path)

    result = ConjugationResult(
        output_dir=artifact_dir,
        generated_sequence=getattr(generations[0], "sequence", None) if generations else None,
        reactive_sequence_index=reactive_sequence_indices[0],
        reactive_residue_selector=reactive_selectors[0],
        conjugate_generation=generations[0] if generations else None,
        construction=construction,
        attachment_specs=specs,
        generated_sequences=tuple(
            generation.sequence for generation in generations if generation is not None
        ),
        reactive_sequence_indices=reactive_sequence_indices,
        reactive_residue_selectors=reactive_selectors,
        conjugate_generations=tuple(
            generation for generation in generations if generation is not None
        ),
        protein_canonicalization=protein_canonicalization,
        relaxed_conjugate_pdb_path=relaxed_pdb,
        solvated_pdb_path=solvated_pdb_path,
        final_interchange_created=(
            exact_export_bundle is not None
            or (False if use_native_glycam else builder.interchange is not None)
        ),
        modifier=modifiers[0],
        modifiers=modifiers,
        relaxed_conjugate_topology=relaxed_topology,
        solvated_topology=builder.solvated_topology,
        final_interchange=None if exact_export_bundle is not None else builder.interchange,
        exact_export_bundle=exact_export_bundle,
        system_builder=builder,
    )
    if exact_export_bundle is not None:
        result.artifact_paths["native_openmm_glycam_audit"] = exact_export_bundle.audit_path
        sidecar_path = getattr(exact_export_bundle, "sidecar_path", None)
        if sidecar_path is not None:
            result.artifact_paths["exact_openmm_exceptions"] = sidecar_path
        result.artifact_paths.update(overlay_artifact_paths)
        _refresh_final_parameter_validation(
            construction=construction,
            resolved_plans=resolved_plans,
            exact_export_bundle=exact_export_bundle,
            artifact_paths=result.artifact_paths,
        )
    workflow_path = artifact_dir / workflow_settings.workflow_json_name
    result.workflow_json_path = workflow_path
    result.artifact_paths["workflow_json"] = workflow_path
    result.save(workflow_path)
    LOGGER.info("Saved conjugation workflow JSON to %s", workflow_path)
    LOGGER.info("Completed config conjugation build in %s", artifact_dir)
    return result


def _refresh_final_parameter_validation(
    *,
    construction: Any,
    resolved_plans: tuple[Any, ...],
    exact_export_bundle: Any,
    artifact_paths: Mapping[str, Path],
) -> None:
    """Rewrite only validation evidence after the final OpenMM System exists."""

    final_system = exact_export_bundle.to_openmm()
    evidence_paths = {
        key: value
        for key, value in artifact_paths.items()
        if key
        in {
            "ownership_manifest",
            "overlay_diagnostics",
            "mixed_overlay_charge_audit",
            "exact_openmm_exceptions",
            "native_openmm_glycam_audit",
        }
    }
    validation_report = build_conjugate_validation_report(
        product_pdb_path=construction.crosslinked_pdb_path,
        resolved_plans=resolved_plans,
        assembly=construction.assembly,
        output_dir=construction.output_dir,
        openmm_system=final_system,
        expected_particle_count=final_system.getNumParticles(),
        parameter_evidence_paths=evidence_paths,
        write=True,
    )
    construction.validation_report_path = validation_report.report_path


def build_direct_moiety_conjugate(
    *,
    protein_pdb_path: Path | str,
    attachments: tuple[Any, ...],
    output_dir: Path | str,
    ccd_pablo: Any | None = None,
    chain_policy: Any | None = None,
    settings: ConjugatedPolymerSystemSettings | None = None,
    random_seed: int | None = None,
) -> ConjugationResult:
    """Build a direct protein plus moiety conjugate request.

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
    specs = [
        _build_attachment_spec(
            attachment,
            attachment_index=index,
            protein_pdb_path=protein_path,
            artifact_dir=moiety_dir,
            workflow_settings=workflow_settings,
            random_seed=random_seed,
            use_cache_dir=False,
        )[0]
        for index, attachment in enumerate(enabled_attachments, start=1)
    ]
    resolved_plans = tuple(specs)

    policy = _policy_with_resolved_crosslinks(
        ccd_pablo or ConjugationCcdPabloPolicyConfig(),
        resolved_plans,
    )
    chain_assignment = chain_policy or ConjugationChainPolicyConfig()
    construction_settings = ModifierConstructionSettings(
        placement=workflow_settings.placement,
        parameterization=workflow_settings.conjugate_parameterization,
        relaxation=workflow_settings.relaxation,
        run_relaxation=workflow_settings.run_relaxation,
    )

    construction, construction_topology = _construct_conjugate_from_specs(
        protein_pdb_path=protein_path,
        specs=tuple(specs),
        ccd_pablo_policy=policy,
        output_dir=construction_dir,
        chain_policy=chain_assignment,
        settings=construction_settings,
        use_product_state_pablo_library=workflow_settings.use_product_state_pablo_library,
        use_conjugate_relaxation=workflow_settings.run_relaxation,
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
        padding=workflow_settings.direct_solvation_padding,
        box_shape=workflow_settings.direct_solvation_box_shape,
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
        modifier=specs[0].fragment,
        modifiers=tuple(spec.fragment for spec in specs),
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


def _build_attachment_spec(
    attachment: Any,
    *,
    attachment_index: int,
    protein_pdb_path: Path | str,
    artifact_dir: Path,
    workflow_settings: ConjugatedPolymerSystemSettings,
    random_seed: int | None = None,
    use_cache_dir: bool = True,
) -> tuple[
    ReactionProduct,
    MultiResidueGenerationResult | None,
    int,
    dict[str, int | str],
    Any,
]:
    """Resolve one config-driven polymer attachment into a generic build spec."""
    _require_supported_coordinate_backend(attachment)
    LOGGER.info(
        "Generating conjugate polymer/moiety for attachment %d (%s)",
        attachment_index,
        _attachment_moiety_name(attachment),
    )
    provider_dir = (
        artifact_dir / workflow_settings.conjugate_cache_dir_name if use_cache_dir else artifact_dir
    )
    source = resolve_moiety_source(
        attachment,
        attachment_index=attachment_index,
        output_dir=provider_dir,
        protein_pdb_path=protein_pdb_path,
        force_regenerate=workflow_settings.force_regenerate_conjugate_polymer,
        max_retries=workflow_settings.conjugate_polymerist_max_retries,
        energy_minimize=workflow_settings.conjugate_polymerist_energy_minimize,
        random_seed=random_seed,
    )
    LOGGER.info("Resolving linkage for attachment %d", attachment_index)
    reaction_template = get_reaction(attachment.mechanism.name)
    settings_builder = getattr(reaction_template, "settings_from_attachment", None)
    reaction_settings = settings_builder(attachment) if callable(settings_builder) else None
    prepared_fragment = prepare_resolved_moiety_source(source)
    reaction_product = reaction_template.resolve_plan(
        protein_pdb_path,
        attachment.site,
        source.reaction_fragment,
        prepared_fragment=prepared_fragment,
        attachment_id=str(attachment.name or f"attachment_{attachment_index:02d}"),
        attachment_index=attachment_index,
        attachment_config=attachment,
        attachment_force_field_domain=str(attachment.moiety.force_field or ""),
        diagnostics=(*source.diagnostics, f"Resolved {source.source_kind} reaction product"),
        settings=reaction_settings,
    )
    return (
        reaction_product,
        source.generation,
        source.reactive_sequence_index if source.reactive_sequence_index is not None else 0,
        source.reactive_selector or {},
        reaction_template,
    )


def _enabled_supported_attachments(
    conjugation: ConjugationConfig | None,
) -> tuple[Any, ...]:
    """Return all enabled attachments supported by the generic config path.

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
        validate_moiety_source_config(
            getattr(attachment, "moiety", None),
            mechanism_name=getattr(getattr(attachment, "mechanism", None), "name", None),
        )
    return attachments


def _uses_pdb_fragment_sources(attachments: tuple[Any, ...]) -> bool:
    """Return whether any enabled attachment uses PDB-fragment input."""
    return any(attachment_uses_pdb_fragment(attachment) for attachment in attachments)


def _build_pdb_fragment_coordinate_only_result(
    *,
    protein_pdb_path: Path | str,
    specs: tuple[ReactionProduct, ...],
    output_dir: Path,
    construction_dir: Path,
    protein_canonicalization: ProteinCanonicalizationResult | None,
    placement_settings: PackmolModifierPlacementSettings | None = None,
    run_packmol_func: Any | None = None,
    output_name: str = "pdb_fragment_coordinate_only_conjugate.pdb",
    status: str = "coordinate_only",
) -> ConjugationResult:
    """Write a coordinate-only residue-resolved PDB-fragment conjugate artifact."""
    output_path, assembly, placement = _write_pdb_fragment_coordinate_artifacts(
        protein_pdb_path=protein_pdb_path,
        specs=specs,
        construction_dir=construction_dir,
        placement_settings=placement_settings,
        run_packmol_func=run_packmol_func,
        output_name=output_name,
    )
    construction = SimpleNamespace(
        crosslinked_pdb_path=output_path,
        validation_report_path=None,
        assembly=assembly,
        placement=placement,
        placements=placement if isinstance(placement, tuple) else (placement,),
        pablo=None,
        parameterization=None,
        relaxation=None,
        diagnostics=("Coordinate-only PDB-fragment artifact; Pablo/OpenFF not run",),
    )
    result = ConjugationResult(
        status=status,
        output_dir=output_dir,
        crosslinked_conjugate_pdb_path=output_path,
        construction=construction,
        attachment_specs=specs,
        generated_sequence=specs[0].fragment.sequence,
        reactive_sequence_index=0,
        reactive_residue_selector={
            "chain_id": "C",
            "atom_serial": specs[0].modifier_link_atom.serial or 0,
            "atom_index": specs[0].modifier_link_atom.atom_index or 0,
            "atom_name": specs[0].modifier_link_atom.atom_name,
        },
        generated_sequences=tuple(
            fragment.sequence
            for fragment in (spec.fragment for spec in specs)
            if fragment.sequence is not None
        ),
        reactive_sequence_indices=tuple(range(len(specs))),
        reactive_residue_selectors=(*(_reactive_selector_for_product(spec) for spec in specs),),
        conjugate_generations=(),
        protein_canonicalization=protein_canonicalization,
        modifier=specs[0].fragment,
        modifiers=tuple(spec.fragment for spec in specs),
        final_interchange_created=False,
    )
    result.artifact_paths.update(
        {
            "crosslinked_conjugate_pdb": output_path,
            **(
                {"pdb_fragment_coordinate_only_pdb": output_path}
                if status == "coordinate_only"
                else {"structure_scope_pdb": output_path}
            ),
            **_pdb_fragment_sidecar_artifacts(specs),
        }
    )
    return result


def _reactive_selector_for_product(product: ReactionProduct) -> dict[str, int | str]:
    """Return a compact reactive atom selector for coordinate-only artifacts."""
    atom = product.modifier_link_atom
    return {
        "chain_id": atom.chain_id,
        "atom_serial": atom.serial or 0,
        "atom_index": atom.atom_index or 0,
        "atom_name": atom.atom_name,
    }


def _pdb_fragment_sidecar_artifacts(specs: tuple[ReactionProduct, ...]) -> dict[str, Path]:
    """Return namespaced PDB-fragment sidecar artifact paths."""
    artifacts: dict[str, Path] = {}
    for index, spec in enumerate(specs, start=1):
        for name, path in spec.fragment.sidecars.items():
            artifacts[f"pdb_fragment_{index}_{name}"] = path
            if index == 1:
                artifacts.setdefault(f"pdb_fragment_{name}", path)
    return artifacts


def _write_pdb_fragment_coordinate_artifact(
    *,
    protein_pdb_path: Path | str,
    spec: ReactionProduct,
    construction_dir: Path,
    placement_settings: PackmolModifierPlacementSettings | None = None,
    run_packmol_func: Any | None = None,
    output_name: str = "pdb_fragment_coordinate_only_conjugate.pdb",
) -> tuple[Path, Any, Any]:
    """Write the residue-preserved PDB-fragment coordinate artifact through Packmol placement."""
    output_path = construction_dir / output_name
    placement, assembly = _place_pdb_fragment_coordinate_only_with_packmol(
        protein_pdb_path,
        spec.fragment,
        spec,
        construction_dir,
        output_path=output_path,
        settings=_pdb_fragment_coordinate_only_placement_settings(placement_settings),
        run_packmol_func=run_packmol_func,
    )
    sidecar_path = spec.fragment.sidecars.get("pdb_fragment_ingestion")
    if sidecar_path is not None:
        _annotate_pdb_fragment_sidecar(sidecar_path, coordinate_artifact_path=output_path)
    return output_path, assembly, placement


def _write_pdb_fragment_coordinate_artifacts(
    *,
    protein_pdb_path: Path | str,
    specs: tuple[ReactionProduct, ...],
    construction_dir: Path,
    placement_settings: PackmolModifierPlacementSettings | None = None,
    run_packmol_func: Any | None = None,
    output_name: str = "pdb_fragment_coordinate_only_conjugate.pdb",
) -> tuple[Path, Any, Any]:
    """Write one coordinate artifact for one or more independent PDB fragments."""
    if len(specs) == 1:
        return _write_pdb_fragment_coordinate_artifact(
            protein_pdb_path=protein_pdb_path,
            spec=specs[0],
            construction_dir=construction_dir,
            placement_settings=placement_settings,
            run_packmol_func=run_packmol_func,
            output_name=output_name,
        )
    output_path = construction_dir / output_name
    placements = place_modifiers_with_resolved_plans(
        protein_pdb_path,
        tuple(spec.fragment for spec in specs),
        specs,
        construction_dir,
        settings=_pdb_fragment_coordinate_only_placement_settings(placement_settings),
        run_packmol_func=run_packmol_func,
    )
    placed_fragments = tuple(
        placement.placed_modifier.model_copy(update={"name": spec.fragment.name})
        for index, (placement, spec) in enumerate(zip(placements, specs, strict=True), start=1)
    )
    assembly = write_crosslinked_pdb(
        protein_pdb_path,
        placed_fragments,
        tuple(spec.to_pdb_linkage_attachment() for spec in specs),
        output_path,
        CrosslinkedPdbAssemblyOptions(
            protein_chain="A",
            polymer_chain="C",
            include_link_records=True,
        ),
    )
    try:
        summary = _summarize_pdb_fragment_true_nonbonded_contacts(
            output_path,
            assembly=assembly,
        )
    except RuntimeError as error:
        for placement in placements:
            _record_final_packmol_validation(
                placement,
                accepted=False,
                diagnostics={
                    "final_conect_graph_valid": False,
                    "final_validation_error": str(error),
                },
            )
        raise
    pair = _format_nonbonded_contact_pair(summary.min_contact)
    accepted = summary.contact_count == 0
    diagnostics = {
        "final_conect_graph_valid": True,
        "true_nonbonded_heavy_contact_count_below_2_angstrom": summary.contact_count,
        "min_true_nonbonded_heavy_distance_angstrom": summary.min_distance_angstrom,
        "min_true_nonbonded_heavy_pair": pair,
    }
    for placement in placements:
        _record_final_packmol_validation(
            placement,
            accepted=accepted,
            diagnostics=diagnostics,
        )
    if not accepted:
        raise RuntimeError(
            "Joint Packmol placement produced severe final protein-fragment or "
            "fragment-fragment heavy-atom clashes below "
            f"{_PDB_FRAGMENT_COORDINATE_ONLY_MIN_HEAVY_NONBONDED_ANGSTROM:.1f} A: "
            f"minimum {summary.min_distance_angstrom:.3f} A for {pair} with "
            f"{summary.contact_count} true contacts"
        )
    placements = tuple(
        placement.model_copy(
            update={
                "min_true_nonbonded_heavy_distance_angstrom": summary.min_distance_angstrom,
                "min_true_nonbonded_heavy_pair": pair,
                "true_nonbonded_heavy_contact_count_below_2_angstrom": summary.contact_count,
                "final_conect_graph_valid": True,
            }
        )
        for placement in placements
    )
    for spec in specs:
        sidecar_path = spec.fragment.sidecars.get("pdb_fragment_ingestion")
        if sidecar_path is not None:
            _annotate_pdb_fragment_sidecar(sidecar_path, coordinate_artifact_path=output_path)
    return output_path, assembly, placements


def _pdb_fragment_coordinate_only_placement_settings(
    settings: PackmolModifierPlacementSettings | None,
) -> PackmolModifierPlacementSettings:
    """Return Packmol settings with a conservative fragment steric tolerance."""
    placement_settings = settings or PackmolModifierPlacementSettings()
    if (
        placement_settings.tolerance_angstrom
        >= _PDB_FRAGMENT_COORDINATE_ONLY_MIN_PACKMOL_TOLERANCE_ANGSTROM
    ):
        return placement_settings
    return placement_settings.model_copy(
        update={
            "tolerance_angstrom": _PDB_FRAGMENT_COORDINATE_ONLY_MIN_PACKMOL_TOLERANCE_ANGSTROM,
        }
    )


def _place_pdb_fragment_coordinate_only_with_packmol(
    protein_pdb_path: Path | str,
    modifier: Any,
    plan: ReactionProduct,
    construction_dir: Path,
    *,
    output_path: Path,
    settings: PackmolModifierPlacementSettings,
    run_packmol_func: Any | None,
) -> tuple[Any, Any]:
    """Place a PDB fragment through Packmol until final-graph clashes are acceptable."""
    failures: list[str] = []
    for attempt_index in range(1, _PDB_FRAGMENT_COORDINATE_ONLY_MAX_PACKMOL_ATTEMPTS + 1):
        attempt_settings = settings
        if attempt_index > 1:
            attempt_settings = settings.model_copy(
                update={
                    "work_dir_name": f"{settings.work_dir_name}_attempt_{attempt_index:02d}",
                    "tolerance_angstrom": settings.tolerance_angstrom
                    + 0.25 * float(attempt_index - 1),
                    "nloop": max(settings.nloop, 1000),
                }
            )
        try:
            placement = place_modifier_with_resolved_plan(
                protein_pdb_path,
                modifier,
                plan,
                construction_dir,
                settings=attempt_settings,
                run_packmol_func=run_packmol_func,
            )
        except PackmolOutputValidationError as error:
            failures.append(error.diagnostic(attempt_index))
            continue
        attempt_path = output_path
        if attempt_index > 1:
            attempt_path = output_path.with_name(
                f"{output_path.stem}_attempt_{attempt_index:02d}{output_path.suffix}"
            )
        assembly = _write_pdb_fragment_placed_artifact(
            protein_pdb_path=protein_pdb_path,
            placed_modifier=placement.placed_modifier,
            plan=plan,
            output_path=attempt_path,
            attachment_id=getattr(modifier, "name", "pdb_fragment"),
        )
        try:
            summary = _summarize_pdb_fragment_true_nonbonded_contacts(attempt_path)
        except RuntimeError as error:
            _record_final_packmol_validation(
                placement,
                accepted=False,
                diagnostics={
                    "final_conect_graph_valid": False,
                    "final_validation_error": str(error),
                },
            )
            raise
        min_distance = summary.min_distance_angstrom or float("inf")
        pair = _format_nonbonded_contact_pair(summary.min_contact)
        if summary.contact_count == 0:
            _record_final_packmol_validation(
                placement,
                accepted=True,
                diagnostics={
                    "final_conect_graph_valid": True,
                    "true_nonbonded_heavy_contact_count_below_2_angstrom": 0,
                    "min_true_nonbonded_heavy_distance_angstrom": summary.min_distance_angstrom,
                },
            )
            placement = placement.model_copy(
                update={
                    "min_true_nonbonded_heavy_distance_angstrom": summary.min_distance_angstrom,
                    "min_true_nonbonded_heavy_pair": pair,
                    "true_nonbonded_heavy_contact_count_below_2_angstrom": summary.contact_count,
                    "final_conect_graph_valid": True,
                }
            )
            if attempt_path != output_path:
                assembly = _write_pdb_fragment_placed_artifact(
                    protein_pdb_path=protein_pdb_path,
                    placed_modifier=placement.placed_modifier,
                    plan=plan,
                    output_path=output_path,
                    attachment_id=getattr(modifier, "name", "pdb_fragment"),
                )
            return placement, assembly
        _record_final_packmol_validation(
            placement,
            accepted=False,
            diagnostics={
                "final_conect_graph_valid": True,
                "true_nonbonded_heavy_contact_count_below_2_angstrom": summary.contact_count,
                "min_true_nonbonded_heavy_distance_angstrom": summary.min_distance_angstrom,
                "min_true_nonbonded_heavy_pair": pair,
            },
        )
        failures.append(
            _format_pdb_fragment_packmol_attempt_diagnostic(
                attempt_index=attempt_index,
                placement=placement,
                min_distance=min_distance,
                pair=pair,
                contact_count=summary.contact_count,
            )
        )
    raise RuntimeError(
        "Packmol could not place the coordinate-only PDB fragment without severe "
        "protein-fragment heavy-atom clashes below "
        f"{_PDB_FRAGMENT_COORDINATE_ONLY_MIN_HEAVY_NONBONDED_ANGSTROM:.1f} A. "
        + "; ".join(failures)
    )


def _record_final_packmol_validation(
    placement: Any,
    *,
    accepted: bool,
    diagnostics: dict[str, object],
) -> Path:
    """Persist final-graph and graph-aware clash validation for a candidate."""
    from polyzymd.utils.packmol import record_packmol_validation

    work_dir = Path(placement.packmol_input_path).parent
    return record_packmol_validation(
        work_dir,
        accepted=accepted,
        diagnostics=diagnostics,
    )


def _write_pdb_fragment_placed_artifact(
    *,
    protein_pdb_path: Path | str,
    placed_modifier: Any,
    plan: ReactionProduct,
    output_path: Path,
    attachment_id: str,
) -> Any:
    """Write one candidate final PDB-fragment coordinate-only artifact."""
    placed_fragment = placed_modifier
    return write_crosslinked_pdb(
        protein_pdb_path,
        placed_fragment.model_copy(update={"name": attachment_id}),
        plan.to_pdb_linkage_attachment(),
        output_path,
        CrosslinkedPdbAssemblyOptions(
            protein_chain="A",
            polymer_chain="C",
            include_link_records=True,
        ),
    )


def _summarize_pdb_fragment_true_nonbonded_contacts(
    output_path: Path,
    *,
    assembly: Any | None = None,
) -> Any:
    """Return graph-distance-aware clash metrics for a final PDB-fragment PDB."""
    atoms = tuple(parse_structure_pdb_atom_records(output_path))
    bonds = _pdb_fragment_final_clash_graph_bonds(output_path, atoms)
    if not validate_conect_graph(atoms, bonds):
        raise RuntimeError(
            f"PDB-fragment coordinate-only CONECT graph has unknown endpoints: {output_path}"
        )
    include_pair = _is_protein_fragment_pair
    if assembly is not None:
        include_pair = _multi_fragment_clash_pair_filter(atoms, assembly)
    return summarize_nonbonded_heavy_clashes(
        atoms,
        bonds,
        cutoff_angstrom=_PDB_FRAGMENT_COORDINATE_ONLY_MIN_HEAVY_NONBONDED_ANGSTROM,
        excluded_bond_depth=3,
        include_pair=include_pair,
    )


def _multi_fragment_clash_pair_filter(atoms: tuple[PdbAtomRecord, ...], assembly: Any) -> Any:
    """Return a filter covering protein-fragment and cross-fragment atom pairs."""
    fragment_by_serial: dict[int, str] = {}
    for source_key, mapping in assembly.atom_mappings.items():
        if not source_key.startswith("fragment_"):
            continue
        target_serial = mapping.get("target_serial")
        if isinstance(target_serial, int):
            fragment_by_serial[target_serial] = source_key.split(":", maxsplit=1)[0]
    fragment_serials = {
        atom.serial for atom in atoms if atom.chain_id.strip().upper() == "C" and atom.serial
    }
    if fragment_serials != set(fragment_by_serial):
        raise RuntimeError(
            "Final multi-fragment clash validation could not map every chain C atom to its "
            "source fragment"
        )

    def include_pair(left: PdbAtomRecord, right: PdbAtomRecord) -> bool:
        if _is_protein_fragment_pair(left, right):
            return True
        if left.chain_id.strip().upper() != "C" or right.chain_id.strip().upper() != "C":
            return False
        return fragment_by_serial[left.serial] != fragment_by_serial[right.serial]

    return include_pair


def _format_nonbonded_contact_pair(contact: Any | None) -> str:
    """Format a nonbonded contact pair for Packmol diagnostics."""
    if contact is None:
        return "none"
    return f"{contact.left_identity}-{contact.right_identity}"


def _format_pdb_fragment_packmol_attempt_diagnostic(
    *,
    attempt_index: int,
    placement: Any,
    min_distance: float,
    pair: str,
    contact_count: int,
) -> str:
    """Return a retry diagnostic for a Packmol placement with final clashes."""
    output_path = Path(placement.packmol_output_path)
    output_exists = output_path.exists()
    atom_count = _safe_pdb_atom_count(output_path) if output_exists else None
    error_log_path = Path(placement.packmol_input_path).with_name("packmol_error.log")
    exit_code = _packmol_exit_code_from_log(error_log_path)
    atom_count_text = "unreadable" if atom_count is None else str(atom_count)
    exit_code_text = "unknown" if exit_code is None else str(exit_code)
    return (
        f"attempt {attempt_index}: exit code {exit_code_text}, "
        f"output exists={output_exists}, atoms={atom_count_text}, "
        f"input={placement.packmol_input_path}, error log={error_log_path}, "
        f"min true heavy nonbonded distance {min_distance:.3f} A for {pair} "
        f"with {contact_count} true contacts"
    )


def _safe_pdb_atom_count(path: Path) -> int | None:
    """Return a PDB atom count, or ``None`` when the output is unreadable."""
    try:
        return len(parse_structure_pdb_atom_records(path))
    except OSError:
        return None


def _packmol_exit_code_from_log(error_log_path: Path) -> int | None:
    """Return the Packmol exit code recorded in a retained error log."""
    from polyzymd.utils.packmol import load_packmol_status

    status = load_packmol_status(error_log_path.parent)
    return_code = status.get("return_code")
    if isinstance(return_code, int):
        return return_code
    if not error_log_path.exists():
        return 0
    text = error_log_path.read_text(encoding="utf-8", errors="replace")
    if "173" in text:
        return 173
    return None


def _pdb_fragment_final_clash_graph_bonds(
    output_path: Path,
    atoms: tuple[PdbAtomRecord, ...],
) -> tuple[tuple[int, int], ...]:
    """Return final graph bonds for PDB-fragment protein-fragment clash classification."""
    conect_bonds = set(parse_pdb_conect_pairs(output_path))
    local_protein_bonds = {
        tuple(sorted((left.serial, right.serial)))
        for index, left in enumerate(atoms)
        for right in atoms[index + 1 :]
        if _is_same_residue_protein_pair(left, right)
        and left.serial is not None
        and right.serial is not None
        and _protein_local_bond_distance(left, right) <= 1.9
    }
    return tuple(sorted(conect_bonds | local_protein_bonds))


def _is_protein_fragment_pair(left: PdbAtomRecord, right: PdbAtomRecord) -> bool:
    """Return whether a pair spans protein chain A and fragment chain C."""
    chains = {left.chain_id.strip().upper(), right.chain_id.strip().upper()}
    return chains == {"A", "C"}


def _is_same_residue_protein_pair(left: PdbAtomRecord, right: PdbAtomRecord) -> bool:
    """Return whether two atoms are in the same protein residue."""
    return (
        left.chain_id.strip().upper() == "A"
        and right.chain_id.strip().upper() == "A"
        and left.residue_number == right.residue_number
        and left.insertion_code == right.insertion_code
    )


def _protein_local_bond_distance(left: PdbAtomRecord, right: PdbAtomRecord) -> float:
    """Return a local protein atom-pair distance for inferred graph edges."""
    return float(
        np.linalg.norm(
            np.asarray([left.x - right.x, left.y - right.y, left.z - right.z], dtype=float)
        )
    )


def _annotate_pdb_fragment_sidecar(sidecar_path: Path, *, coordinate_artifact_path: Path) -> None:
    """Add the coordinate-only artifact path to an existing PDB-fragment sidecar."""
    payload = json.loads(sidecar_path.read_text(encoding="utf-8"))
    payload["coordinate_artifact_path"] = str(coordinate_artifact_path)
    sidecar_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


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


def _require_supported_coordinate_backend(attachment: Any) -> None:
    """Gate config-driven coordinate construction to implemented backends."""
    mechanism = attachment.mechanism
    mechanism_name = mechanism.name.strip().lower()
    try:
        reaction_template = get_reaction(mechanism_name)
    except KeyError:
        reaction_template = None
    if reaction_template is None:
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
                f"'{_NHS_LYS_COORDINATE_BACKEND_MECHANISM}' or registered generic "
                "attachment mechanisms."
            )
        raise NotImplementedError(
            "Config-driven conjugated polymer system construction currently implements coordinate "
            "surgery only for registered generic attachment mechanisms. "
            f"Received mechanism '{mechanism.name}'. Provide a built-in mechanism with "
            "resolve_plan() support."
        )
    if getattr(reaction_template, "supports_coordinate_assembly", False) or callable(
        getattr(reaction_template, "resolve_plan", None)
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
            "workflow currently requires a registered reaction template with coordinate "
            "assembly capability."
        )

    raise NotImplementedError(
        "Config-driven conjugated polymer system construction currently implements coordinate "
        "surgery only for registered generic attachment mechanisms. "
        f"Received mechanism '{mechanism.name}'. Provide a built-in mechanism with "
        "resolve_plan() support."
    )


def _policy_with_resolved_crosslink(
    policy: Any,
    resolved_plan: ReactionProduct,
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


def _product_state_crosslink_requirement(resolved_plan: ReactionProduct):
    """Return the Pablo crosslink requirement for an already-modified product PDB.

    The resolved attachment plan records reactant-state leaving atoms for graph
    surgery. By the time Pablo ingests the emitted crosslinked PDB, PolyzyMD has
    already removed those atoms and written the product-state ``CONECT`` bond, so
    the Pablo policy must not ask Pablo to remove them a second time.
    """
    requirement = resolved_plan.pablo_crosslink_requirement
    return requirement.model_copy(update={"leaving_atoms": ((), ())})


def _construct_conjugate_from_specs(
    *,
    protein_pdb_path: Path | str,
    specs: tuple[ReactionProduct, ...],
    ccd_pablo_policy: Any,
    output_dir: Path | str,
    chain_policy: Any | None,
    settings: ModifierConstructionSettings,
    use_product_state_pablo_library: bool,
    use_conjugate_relaxation: bool = True,
) -> tuple[Any, Any]:
    """Construct, parameterize, and relax a conjugate from resolved attachment specs."""
    if not specs:
        raise ValueError("Conjugate construction requires at least one attachment spec")

    modifiers = tuple(spec.fragment for spec in specs)
    resolved_plans = specs
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
        placement.placed_modifier
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
    product_state_specs = _product_state_specs_with_assembly_mappings(
        specs,
        assembly_result=assembly_result,
        product_pdb_path=crosslinked_pdb_path,
    )
    pablo_product_pdb_path = _write_scoped_pablo_ingestion_pdb(
        crosslinked_pdb_path,
        product_state_specs,
        artifact_dir / "assembled_crosslinked.pablo_scoped.pdb",
    )

    product_state_pablo_library = None
    product_state_residue_library = None
    if use_product_state_pablo_library:
        LOGGER.info("Building product-state Pablo residue library")
        product_state_pablo_library = _product_state_pablo_library_for_specs(
            product_pdb=pablo_product_pdb_path,
            source_protein_pdb=protein_pdb_path,
            specs=product_state_specs,
        )
        product_state_residue_library = product_state_pablo_library.residue_library

    LOGGER.info("Ingesting product-state PDB with Pablo")
    pablo_result = PabloIngestor(policy=ccd_pablo_policy).ingest_structure(
        pablo_product_pdb_path,
        chain_policy=chain_policy,
        output_dir=artifact_dir,
        residue_library=product_state_residue_library,
    )
    if not pablo_result.success or pablo_result.topology is None:
        raise RuntimeError(_pablo_failure_message(pablo_result, specs=product_state_specs))
    if product_state_pablo_library is not None:
        _apply_pdb_atom_identity_to_topology(pablo_result.topology, pablo_product_pdb_path)
        LOGGER.info("Building preproduction product-state charge bridge")
        product_state_pablo_library = _product_state_library_with_charge_bridge(
            product_state_pablo_library,
            product_topology=pablo_result.topology,
            product_pdb=pablo_product_pdb_path,
            source_protein_pdb=protein_pdb_path,
            specs=product_state_specs,
            output_dir=artifact_dir,
            settings=settings.parameterization,
        )
        product_state_pablo_library = _product_state_library_with_charge_templates(
            product_state_pablo_library,
            pablo_result.topology,
        )
    else:
        _apply_pdb_atom_identity_to_topology(pablo_result.topology, pablo_product_pdb_path)

    charge_templates, require_charge_templates = _intermediate_conjugate_charge_templates(
        pablo_result.topology,
        product_state_pablo_library,
    )
    _apply_pdb_atom_identity_to_topology(pablo_result.topology, crosslinked_pdb_path)
    _restore_scoped_alias_metadata(pablo_result.topology, product_state_specs)
    LOGGER.info("Parameterizing conjugate with OpenFF Interchange")
    parameterization_result = create_interchange_from_pablo_topology(
        pablo_result.topology,
        settings=settings.parameterization,
        charge_from_molecules=charge_templates,
        require_charge_templates=require_charge_templates,
    )
    if not parameterization_result.success or parameterization_result.interchange is None:
        raise RuntimeError("OpenFF Interchange parameterization did not produce an interchange")

    site_label = "site" if len(resolved_plans) == 1 else "sites"
    diagnostics = [
        f"Conjugate construction completed ({len(resolved_plans)} attachment {site_label})",
    ]
    relaxation_result = None
    if use_conjugate_relaxation and settings.run_relaxation:
        LOGGER.info("Running conjugate relaxation")
        relaxation_result = relax_conjugate(
            parameterization_result.interchange,
            artifact_dir,
            product_pdb_path=crosslinked_pdb_path,
            attachment_specs=product_state_specs,
            assembly=assembly_result,
            settings=settings.relaxation,
        )
        if not relaxation_result.success:
            raise RuntimeError("Conjugate relaxation did not report success")

    validation_report = build_conjugate_validation_report(
        product_pdb_path=crosslinked_pdb_path,
        resolved_plans=resolved_plans,
        assembly=assembly_result,
        output_dir=artifact_dir,
        expected_particle_count=getattr(pablo_result.topology, "n_atoms", None),
        write=True,
    )
    relaxation_was_requested = use_conjugate_relaxation and settings.run_relaxation
    relaxation_evidence_status = validation_report.relaxation_evidence.status
    if relaxation_was_requested and relaxation_evidence_status != ValidationStatus.PASS:
        raise RuntimeError(
            "Conjugate validation report relaxation evidence did not pass "
            f"after requested relaxation: {relaxation_evidence_status.value}"
        )
    if not relaxation_was_requested and relaxation_evidence_status == ValidationStatus.FAIL:
        raise RuntimeError("Conjugate validation report relaxation evidence failed")

    return (
        ConjugateConstructionResult(
            output_dir=artifact_dir,
            reaction_product=product_state_specs[0],
            reaction_products=product_state_specs,
            crosslink_validation=crosslink_validations[0],
            crosslink_validations=crosslink_validations,
            placement=placements[0],
            placements=placements,
            assembly=assembly_result,
            pablo=pablo_result,
            parameterization=parameterization_result,
            relaxation=relaxation_result,
            product_state_pablo_library=product_state_pablo_library,
            crosslinked_pdb_path=crosslinked_pdb_path,
            validation_report_path=validation_report.report_path,
            diagnostics=tuple(diagnostics),
        ),
        pablo_result.topology,
    )


def _product_state_specs_with_assembly_mappings(
    specs: tuple[Any, ...],
    *,
    assembly_result: Any,
    product_pdb_path: Path | str,
) -> tuple[Any, ...]:
    """Return specs whose plans point at exact product-PDB linkage atoms."""
    mappings = assembly_result.residue_mappings
    added_pairs = tuple(assembly_result.added_conect_pairs)
    if len(added_pairs) != len(specs):
        raise ValueError(
            "Product-state spec mapping requires one ordered assembly CONECT pair per "
            f"attachment spec (pairs={len(added_pairs)}, specs={len(specs)})"
        )
    product_atom_by_serial = _product_pdb_atoms_by_serial(product_pdb_path)
    _validate_added_conect_pairs(product_pdb_path, added_pairs, assembly_result=assembly_result)
    updated_specs = []
    alias_cursor = 1
    for fragment_index, (spec, pair) in enumerate(zip(specs, added_pairs, strict=True), start=1):
        plan = spec
        fragment_prefix = f"fragment_{fragment_index}:"
        fragment_mappings = {
            key.removeprefix(fragment_prefix): value
            for key, value in mappings.items()
            if key.startswith(fragment_prefix)
        }
        protein_atom = _product_atom_for_serial(product_atom_by_serial, pair[0], "protein")
        modifier_atom = _product_atom_for_serial(product_atom_by_serial, pair[1], "modifier")
        expected_protein_chain = _assembly_protein_chain(assembly_result, fragment_index)
        _validate_product_pair_for_spec(
            spec,
            protein_atom=protein_atom,
            modifier_atom=modifier_atom,
            fragment_index=fragment_index,
            fragment_mappings=fragment_mappings,
            expected_protein_chain=expected_protein_chain,
        )
        fragment_mappings, alias_map, alias_cursor = _fragment_mappings_with_scoped_aliases(
            fragment_mappings,
            product_atom_by_serial,
            alias_cursor,
        )
        endpoint_provenance = _endpoint_provenance_for_spec(
            assembly_result,
            fragment_index=fragment_index,
            protein_atom=protein_atom,
            modifier_atom=modifier_atom,
            conect_pair=pair,
        )
        updated_specs.append(
            _copy_spec_with_product_mappings(
                spec,
                fragment_mappings,
                endpoint_provenance=endpoint_provenance,
                scoped_residue_aliases=alias_map,
            )
        )
    return tuple(updated_specs)


def _validate_added_conect_pairs(
    product_pdb_path: Path | str,
    added_pairs: tuple[tuple[int, int], ...],
    *,
    assembly_result: Any,
) -> None:
    """Verify every assembly endpoint pair is present in the emitted CONECT graph."""
    if not assembly_result.attachment_endpoint_records:
        return
    conect_pairs = {frozenset(pair) for pair in parse_pdb_conect_pairs(Path(product_pdb_path))}
    for pair in added_pairs:
        if frozenset(pair) not in conect_pairs:
            raise ValueError(
                "Crosslinked PDB assembly endpoint provenance is inconsistent: emitted "
                f"CONECT records do not contain pair {pair} in {product_pdb_path}"
            )


def _fragment_mappings_with_scoped_aliases(
    fragment_mappings: dict[str, dict[str, int | str]],
    product_atom_by_serial: dict[int, PdbAtomRecord],
    alias_cursor: int,
) -> tuple[dict[str, dict[str, int | str]], dict[str, str], int]:
    """Attach attachment-local Pablo residue aliases to product residue mappings."""
    atoms_by_residue = _product_atoms_by_residue(product_atom_by_serial.values())
    updated: dict[str, dict[str, int | str]] = {}
    alias_map: dict[str, str] = {}
    for source_key, mapping in sorted(fragment_mappings.items()):
        chain = str(mapping.get("target_chain", "C") or "C")
        residue_number = int(mapping["target_residue_number"])
        insertion_code = str(mapping.get("target_insertion_code", "") or "")
        residue_atoms = atoms_by_residue.get((chain, residue_number, insertion_code), ())
        if not residue_atoms:
            raise ValueError(
                "Product-state scoped Pablo aliasing could not find emitted residue "
                f"{chain}:{residue_number}{insertion_code} for source {source_key}"
            )
        canonical_name = residue_atoms[0].residue_name.strip().upper()
        scoped_name = _scoped_pablo_residue_name(alias_cursor)
        alias_cursor += 1
        residue_key = _product_residue_alias_key(chain, residue_number, insertion_code)
        alias_map[residue_key] = scoped_name
        enriched = dict(mapping)
        enriched["canonical_residue_name"] = canonical_name
        enriched["scoped_residue_name"] = scoped_name
        updated[source_key] = enriched
    return updated, alias_map, alias_cursor


def _product_atoms_by_residue(
    atoms: Iterable[PdbAtomRecord],
) -> dict[tuple[str, int, str], tuple[PdbAtomRecord, ...]]:
    """Group product atoms by chain, residue number, and insertion code."""
    grouped: dict[tuple[str, int, str], list[PdbAtomRecord]] = {}
    for atom in atoms:
        key = (atom.chain_id, atom.residue_number, atom.insertion_code or "")
        grouped.setdefault(key, []).append(atom)
    return {key: tuple(value) for key, value in grouped.items()}


def _scoped_pablo_residue_name(index: int) -> str:
    """Return a deterministic three-character residue alias for Pablo matching."""
    if index < 1 or index > 899:
        raise ValueError("Attachment-scoped Pablo aliases exceeded the supported residue count")
    if index <= 99:
        return f"Z{index:02d}"
    block, suffix = divmod(index - 100, 100)
    prefix = chr(ord("Y") - block)
    return f"{prefix}{suffix:02d}"


def _product_residue_alias_key(chain: str, residue_number: int, insertion_code: str) -> str:
    """Return the scoped-alias map key for one emitted product residue."""
    return f"{chain}:{residue_number}{insertion_code or ''}"


def _endpoint_provenance_for_spec(
    assembly_result: Any,
    *,
    fragment_index: int,
    protein_atom: PdbAtomRecord,
    modifier_atom: PdbAtomRecord,
    conect_pair: tuple[int, int],
) -> dict[str, Any]:
    """Return serial-first endpoint provenance for one attachment spec."""
    endpoint_records = tuple(assembly_result.attachment_endpoint_records)
    record = next(
        (
            item
            for item in endpoint_records
            if int(item.get("attachment_index", 0)) == fragment_index
        ),
        None,
    )
    return {
        "attachment_index": fragment_index,
        "conect_pair": {"protein_serial": conect_pair[0], "modifier_serial": conect_pair[1]},
        "protein_endpoint": (
            record.get("protein_endpoint") if record else _pdb_atom_payload(protein_atom)
        ),
        "modifier_endpoint": (
            record.get("modifier_endpoint") if record else _pdb_atom_payload(modifier_atom)
        ),
        "atom_mappings": dict(assembly_result.atom_mappings),
    }


def _pdb_atom_payload(atom: PdbAtomRecord) -> dict[str, int | str]:
    """Return JSON-safe PDB atom identity payload."""
    return {
        "serial": int(atom.serial) if atom.serial is not None else "",
        "atom_name": atom.atom_name,
        "residue_name": atom.residue_name,
        "chain_id": atom.chain_id,
        "residue_number": atom.residue_number,
        "insertion_code": atom.insertion_code or "",
    }


def _product_pdb_atoms_by_serial(product_pdb_path: Path | str) -> dict[int, PdbAtomRecord]:
    """Return product PDB atom records keyed by unique serial."""
    atoms_by_serial: dict[int, PdbAtomRecord] = {}
    for atom in parse_structure_pdb_atom_records(Path(product_pdb_path), require_atoms=True):
        serial = atom.serial
        if serial is None:
            raise ValueError(f"Product PDB atom {atom.atom_name} is missing a serial")
        if serial in atoms_by_serial:
            raise ValueError(f"Product PDB contains duplicate atom serial {serial}")
        atoms_by_serial[serial] = atom
    return atoms_by_serial


def _product_atom_for_serial(
    product_atom_by_serial: dict[int, PdbAtomRecord], serial: int, role: str
) -> PdbAtomRecord:
    """Return the product atom record for an assembly CONECT serial."""
    atom = product_atom_by_serial.get(int(serial))
    if atom is None:
        raise ValueError(f"Assembly {role} CONECT serial {serial} is missing from product PDB")
    return atom


def _validate_product_pair_for_spec(
    spec: Any,
    *,
    protein_atom: PdbAtomRecord,
    modifier_atom: PdbAtomRecord,
    fragment_index: int,
    fragment_mappings: dict[str, dict[str, int | str]],
    expected_protein_chain: str | None = None,
) -> None:
    """Validate an ordered product linkage pair against one attachment spec."""
    plan = spec
    _validate_product_protein_atom(
        plan,
        protein_atom,
        fragment_index,
        expected_chain_id=expected_protein_chain,
    )
    _validate_product_modifier_atom(plan, modifier_atom, fragment_index, fragment_mappings)


def _validate_product_protein_atom(
    plan: Any,
    atom: PdbAtomRecord,
    fragment_index: int,
    *,
    expected_chain_id: str | None = None,
) -> None:
    """Validate a product PDB protein endpoint for one resolved plan."""
    selector = getattr(getattr(plan, "contract", None), "protein_endpoint", None)
    selector = getattr(selector, "selector", None)
    source_atom = plan.protein_link_atom
    expected_chain = expected_chain_id or getattr(selector, "chain_id", source_atom.chain_id)
    expected_residue_number = getattr(selector, "residue_number", source_atom.residue_number)
    expected_insertion_code = getattr(selector, "insertion_code", source_atom.insertion_code)
    expected_atom_name = source_atom.atom_name
    expected_residue_names = {
        str(source_atom.residue_name).upper(),
        str(plan.protein_product_residue_name).upper(),
    }
    if (
        atom.chain_id.upper() != str(expected_chain).upper()
        or atom.residue_number != int(expected_residue_number)
        or (atom.insertion_code or "").upper() != str(expected_insertion_code or "").upper()
        or atom.atom_name.upper() != str(expected_atom_name).upper()
        or atom.residue_name.upper() not in expected_residue_names
    ):
        raise ValueError(
            "Assembly CONECT pair protein endpoint does not match attachment "
            f"{fragment_index}: serial={atom.serial} {atom.chain_id}:{atom.residue_name}:"
            f"{atom.residue_number}{atom.insertion_code}:{atom.atom_name}"
        )


def _assembly_protein_chain(assembly_result: Any, fragment_index: int) -> str | None:
    """Return the canonical protein chain emitted for one source attachment."""

    for record in tuple(assembly_result.attachment_endpoint_records):
        if int(record.get("attachment_index", 0)) != fragment_index:
            continue
        endpoint = record.get("protein_endpoint", {})
        chain_id = endpoint.get("chain_id") if isinstance(endpoint, dict) else None
        return str(chain_id) if chain_id else None
    return None


def _validate_product_modifier_atom(
    plan: Any,
    atom: PdbAtomRecord,
    fragment_index: int,
    fragment_mappings: dict[str, dict[str, int | str]],
) -> None:
    """Validate a product PDB modifier endpoint for one resolved plan."""
    source_atom = plan.modifier_link_atom
    expected_atom_name = source_atom.atom_name.upper()
    expected_residue_names = {
        str(source_atom.residue_name).upper(),
        str(plan.modifier_product_residue_name).upper(),
    }
    mapping = fragment_mappings.get(_modifier_source_residue_key(source_atom), {})
    has_mapping = bool(mapping)
    mismatches = []
    if atom.atom_name.upper() != expected_atom_name:
        mismatches.append(f"atom name expected={expected_atom_name!r} observed={atom.atom_name!r}")
    if atom.residue_name.upper() not in expected_residue_names:
        mismatches.append(
            "residue name expected one of "
            f"{sorted(expected_residue_names)!r} observed={atom.residue_name!r}"
        )
    expected_chain = str(mapping.get("target_chain", "") or "").upper()
    expected_residue_number = mapping.get("target_residue_number")
    expected_insertion_code = str(mapping.get("target_insertion_code", "") or "").upper()
    if expected_chain and atom.chain_id.upper() != expected_chain:
        mismatches.append(f"chain expected={expected_chain!r} observed={atom.chain_id!r}")
    if expected_residue_number not in {None, ""} and atom.residue_number != int(
        expected_residue_number
    ):
        mismatches.append(
            f"residue number expected={expected_residue_number} observed={atom.residue_number}"
        )
    if has_mapping and (atom.insertion_code or "").upper() != expected_insertion_code:
        mismatches.append(
            f"insertion code expected={expected_insertion_code!r} observed={atom.insertion_code!r}"
        )
    if mismatches:
        raise ValueError(
            "Assembly CONECT pair modifier endpoint does not match attachment "
            f"{fragment_index}: serial={atom.serial} {atom.chain_id}:{atom.residue_name}:"
            f"{atom.residue_number}{atom.insertion_code}:{atom.atom_name}; " + "; ".join(mismatches)
        )


def _copy_spec_with_product_mappings(
    spec: Any,
    product_residue_mappings: dict[str, dict[str, int | str]],
    **updates: Any,
) -> Any:
    """Copy a spec while attaching assembly residue mappings for Pablo templating."""
    updates["product_residue_mappings"] = product_residue_mappings
    return spec.model_copy(update=updates)


def _modifier_source_residue_key(atom: PdbAtomRecord) -> str:
    """Return the residue-mapping key suffix used by the PDB assembly writer."""
    return f"{atom.residue_number}{atom.insertion_code or ''}"


def _policy_with_resolved_crosslinks(
    policy: Any,
    resolved_plans: tuple[ReactionProduct, ...],
) -> Any:
    """Return a Pablo policy containing product-state crosslinks for all plans."""
    updated = policy
    for plan in resolved_plans:
        updated = _policy_with_resolved_crosslink(updated, plan)
    return updated


def _product_state_pablo_library_for_specs(
    *,
    product_pdb: Path,
    source_protein_pdb: Path | str,
    specs: tuple[ReactionProduct, ...],
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


def _write_scoped_pablo_ingestion_pdb(
    canonical_pdb_path: Path | str,
    specs: tuple[Any, ...],
    scoped_pdb_path: Path,
) -> Path:
    """Write a Pablo-only PDB with attachment-local residue-name aliases."""
    alias_by_residue: dict[str, str] = {}
    for spec in specs:
        alias_by_residue.update(getattr(spec, "scoped_residue_aliases", {}) or {})
    if not alias_by_residue:
        return Path(canonical_pdb_path)

    updated_lines: list[str] = []
    with Path(canonical_pdb_path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(_ATOM_RECORD_PREFIXES):
                key = _pdb_line_residue_alias_key(line)
                alias = alias_by_residue.get(key)
                if alias is not None:
                    line = _replace_pdb_residue_name_field(line, alias)
            updated_lines.append(line)
    scoped_pdb_path.write_text("".join(updated_lines), encoding="utf-8")
    return scoped_pdb_path


def _pdb_line_residue_alias_key(line: str) -> str:
    """Return the scoped alias key from one fixed-width PDB atom line."""
    chain = line[21:22].strip() or " "
    residue_number = int(line[22:26])
    insertion_code = line[26:27].strip()
    return _product_residue_alias_key(chain, residue_number, insertion_code)


def _replace_pdb_residue_name_field(line: str, residue_name: str) -> str:
    """Replace the fixed-width PDB residue name field."""
    newline = "\n" if line.endswith("\n") else ""
    body = line[:-1] if newline else line
    padded = body.ljust(20)
    return f"{padded[:17]}{residue_name.strip().upper():>3}{padded[20:]}{newline}"


def _restore_scoped_alias_metadata(topology: Any, specs: tuple[Any, ...]) -> None:
    """Restore canonical GLYCAM residue names after Pablo-only alias matching."""
    canonical_by_alias: dict[str, str] = {}
    for spec in specs:
        for mapping in (getattr(spec, "product_residue_mappings", {}) or {}).values():
            alias = str(mapping.get("scoped_residue_name", "") or "").strip().upper()
            canonical = str(mapping.get("canonical_residue_name", "") or "").strip().upper()
            if alias and canonical:
                canonical_by_alias[alias] = canonical
    if not canonical_by_alias:
        return
    for atom in _iter_openff_topology_atoms(topology):
        metadata = getattr(atom, "metadata", {}) or {}
        residue_name = str(metadata.get("residue_name", "") or "").strip().upper()
        canonical = canonical_by_alias.get(residue_name)
        if canonical is None:
            continue
        atom.name = str(metadata.get("atom_name", getattr(atom, "name", ""))).strip()
        _update_atom_metadata(
            atom,
            {
                "residue_name": canonical,
                "canonical_residue_name": canonical,
                "pablo_scoped_residue_name": residue_name,
                "product_identity_source": "canonical_product_pdb_after_pablo_alias",
            },
        )


def _safe_attachment_token(name: str) -> str:
    """Return a conservative artifact-directory token."""
    token = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in name.strip())
    return token.strip("_-") or "attachment"


def _save_direct_workflow_summary(
    result: Any,
    attachments: tuple[Any, ...],
    resolved_plans: list[ReactionProduct],
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


def _pablo_failure_message(result: Any, *, specs: tuple[Any, ...] = ()) -> str:
    """Return an actionable Pablo failure message with attachment context."""
    diagnostics = [f"{diag.code}: {diag.message}" for diag in result.diagnostics]
    joined = "; ".join(diagnostics) if diagnostics else "no diagnostics were reported"
    details = []
    for spec in specs:
        endpoint = getattr(spec, "endpoint_provenance", {}) or {}
        sidecars = spec.fragment.sidecars
        details.append(
            {
                "attachment_id": getattr(spec, "attachment_id", ""),
                "attachment_index": getattr(spec, "attachment_index", ""),
                "source_glycan": str(
                    sidecars.get("pdb") or sidecars.get("pdb_fragment_ingestion") or ""
                ),
                "conect_pair": endpoint.get("conect_pair", {}),
                "protein_endpoint": endpoint.get("protein_endpoint", {}),
                "modifier_endpoint": endpoint.get("modifier_endpoint", {}),
            }
        )
    context = f" attachment_context={json.dumps(details, sort_keys=True)}" if details else ""
    return f"Pablo ingestion failed or returned no topology for {result.path}: {joined}{context}"


def _relaxed_conjugate_pdb(construction: ModifierConstructionResult) -> Path:
    relaxation = getattr(construction, "relaxation", None)
    if relaxation is not None:
        if relaxation.relaxed_pdb_path is not None:
            return relaxation.relaxed_pdb_path
    crosslinked = getattr(construction, "crosslinked_pdb_path", None)
    if crosslinked is not None:
        return crosslinked
    raise RuntimeError("Conjugate construction did not produce a downstream-ready PDB")


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


def _intermediate_conjugate_charge_templates(
    topology: Any,
    product_state_pablo_library: Any | None,
) -> tuple[tuple[Any, ...], bool]:
    """Build validated charge templates for intermediate conjugate parameterization.

    Parameters
    ----------
    topology : Any
        Current product-state Pablo topology to parameterize.
    product_state_pablo_library : Any or None
        Product-state Pablo library carrying validated/explicit partial-charge provenance.

    Returns
    -------
    tuple of tuple of Any and bool
        Charged molecule templates and whether OpenFF parameterization must
        reject empty templates before implicit charge assignment can occur.

    Raises
    ------
    RuntimeError
        If product-state provenance is active but validated charge templates
        cannot be built.
    """
    if product_state_pablo_library is None:
        return (), False

    try:
        templates = build_conjugate_charge_templates(topology, product_state_pablo_library)
    except ValueError as exc:
        raise RuntimeError(
            "Product-state conjugate parameterization requires validated charge templates from "
            "the product-state Pablo library. Refusing to let OpenFF fall back to whole-conjugate "
            f"implicit charge assignment. Original error: {exc}"
        ) from exc
    return templates, True


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
    """Attach validated/explicit partial-charge bridge records to the product library."""
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
    _ = product_specs
    return product_state_library_has_provenance(product_state_pablo_library)


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
    return product_residue_names(product_state_pablo_library)


def _molecule_has_product_residue(molecule: Any, product_names: set[str]) -> bool:
    """Return whether a molecule contains any product-state residue metadata."""
    return molecule_contains_product_residue(molecule, product_names)


def _metadata_residue_name(metadata: Any) -> str:
    """Return an uppercase residue name from atom or molecule metadata."""
    return metadata_residue_name(metadata)


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
    builder._n_enzyme_molecules = _topology_molecule_count(relaxed_conjugate_topology)
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

    # Exact OpenMM handoffs snapshot the parameterized topology, so canonical
    # PDB identities must exist before parameterization rather than only at export.
    builder._assign_pdb_identifiers()

    if create_interchange:
        LOGGER.info("Creating final solvated OpenFF Interchange")
        create_final_conjugated_interchange(
            builder,
            product_state_pablo_library=product_state_pablo_library,
            settings=parameterization_settings,
        )

    return builder


def _build_isolated_conjugate_system(
    config: SimulationConfig,
    *,
    relaxed_conjugate_topology: Any,
    working_dir: Path,
    product_state_pablo_library: Any | None = None,
    parameterization_settings: InterchangeParameterizationSettings | None = None,
    create_interchange: bool = True,
) -> SystemBuilder:
    """Parameterize the assembled conjugate without secondary components."""
    builder = SystemBuilder.from_config(config)
    builder._working_dir = working_dir
    builder._enzyme_topology = relaxed_conjugate_topology
    builder._n_enzyme_molecules = _topology_molecule_count(relaxed_conjugate_topology)
    builder._preserve_enzyme_chain_ids = True
    builder.combine_solutes()
    builder._solvated_topology = builder._combined_topology
    builder._assign_pdb_identifiers()
    if create_interchange:
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
    padding: float = 0.8,
    box_shape: str = "cube",
) -> SystemBuilder:
    """Solvate a direct SMILES-moiety conjugate without a full SimulationConfig."""
    builder = SystemBuilder()
    builder._working_dir = working_dir
    builder._enzyme_topology = relaxed_conjugate_topology
    builder._n_enzyme_molecules = _topology_molecule_count(relaxed_conjugate_topology)
    builder._preserve_enzyme_chain_ids = True
    LOGGER.info("Combining direct conjugate solutes")
    builder.combine_solutes()
    LOGGER.info("Solvating direct conjugated system")
    builder.solvate(padding=padding, box_shape=box_shape)
    builder._assign_pdb_identifiers()
    if create_interchange:
        LOGGER.info("Creating final solvated OpenFF Interchange")
        create_final_conjugated_interchange(
            builder,
            product_state_pablo_library=product_state_pablo_library,
            settings=parameterization_settings,
        )
    return builder


def _topology_molecule_count(topology: Any) -> int:
    """Return a validated molecule count for a relaxed conjugate topology."""

    count = int(topology.n_molecules)
    if count < 1:
        raise ValueError("Relaxed conjugate topology must contain at least one molecule")
    return count


def _create_mixed_overlay_exact_handoff(
    builder: SystemBuilder,
    *,
    config: Any,
    construction: Any,
    output_dir: Path | str,
    workflow_settings: ConjugatedPolymerSystemSettings,
) -> tuple[Any, dict[str, Path]]:
    """Create an exact bundle from a merged mixed GLYCAM/OpenFF OpenMM System.

    The baseline Interchange is fully parameterized first and supplies all generic
    protein/polymer/linkage terms. A native GLYCAM reference is then generated for
    the same solvated topology, scoped GLYCAM terms are overlaid, and the final
    exact bundle is created from the merged System rather than the stale baseline.
    """

    from openmm import unit

    from polyzymd.builders.conjugation.system_overlay import merge_openmm_system_overlay

    if builder.interchange is None:
        raise RuntimeError("Mixed overlay requires a complete baseline OpenFF Interchange")

    baseline_topology = builder.interchange.to_openmm_topology()
    baseline_positions = builder.interchange.positions.to_openmm()
    baseline_system = builder.interchange.to_openmm_system()
    native_audit: dict[str, Any]
    native_audit_path: Path | None
    native_bundle = _build_mixed_overlay_native_reference(
        config=config,
        output_dir=Path(output_dir) / "mixed_overlay_native_reference",
        free_polymer_seed=None,
        workflow_settings=workflow_settings,
    )
    native_system = native_bundle.system
    native_topology = native_bundle.topology
    native_audit = native_bundle.audit
    native_audit_path = native_bundle.audit_path
    glycam_particles = _glycam_particles_from_native_audit(native_topology, native_audit)
    if not glycam_particles:
        raise RuntimeError("Mixed overlay could not infer any GLYCAM-owned particles")
    atom_mapping = _build_native_to_baseline_atom_mapping(
        baseline_topology,
        native_topology,
        required_native_indices=glycam_particles,
        allow_scoped_fallback=True,
    )
    attachments = tuple(
        resolved.to_dict() for resolved in resolve_conjugate_force_fields(config).attachments
    )
    overlay = merge_openmm_system_overlay(
        baseline_system=baseline_system,
        native_system=native_system,
        glycam_particles=glycam_particles,
        atom_mapping=atom_mapping,
        attachments=attachments,
    )
    charge_audit_path = _write_mixed_overlay_charge_audit(
        overlay.system,
        Path(output_dir) / "mixed_overlay_charge_audit.json",
        unit_module=unit,
    )
    overlay_paths = overlay.save_artifacts(output_dir)
    overlay_paths["mixed_overlay_charge_audit"] = charge_audit_path
    audit = {
        "route": "mixed_overlay",
        "native_reference_audit": str(native_audit_path) if native_audit_path is not None else None,
        "ownership_manifest": str(overlay_paths["ownership_manifest"]),
        "overlay_diagnostics": str(overlay_paths["overlay_diagnostics"]),
        "charge_audit": str(charge_audit_path),
        "glycam_particle_count": len(glycam_particles),
        "native_audit": native_audit,
        "overlay": overlay.diagnostics,
    }
    _stringify_openmm_topology_ids(baseline_topology)
    exact_bundle = create_exact_export_bundle(
        topology=baseline_topology,
        system=overlay.system,
        positions=baseline_positions,
        output_dir=output_dir,
        audit_path=overlay_paths["overlay_diagnostics"],
        audit=audit,
    )
    builder._exact_export_bundle = exact_bundle
    return exact_bundle, overlay_paths


def _stringify_openmm_topology_ids(topology: Any) -> None:
    """Normalize OpenMM topology IDs to strings for PDB writing."""

    for chain in topology.chains():
        chain.id = str(chain.id)
    for residue in topology.residues():
        residue.id = str(residue.id)


def _write_mixed_overlay_charge_audit(
    system: Any,
    path: Path,
    *,
    unit_module: Any,
) -> Path:
    """Write final charge audit for a mixed overlay System."""

    total_charge = 0.0
    particle_count = system.getNumParticles()
    for force in system.getForces():
        if force.__class__.__name__ != "NonbondedForce":
            continue
        for index in range(force.getNumParticles()):
            charge, _sigma, _epsilon = force.getParticleParameters(index)
            total_charge += float(charge.value_in_unit(unit_module.elementary_charge))
        payload = {
            "particle_count": particle_count,
            "nonbonded_particle_count": force.getNumParticles(),
            "total_charge_e": total_charge,
            "rounded_total_charge_e": round(total_charge),
            "charge_mismatch_e": total_charge - round(total_charge),
        }
        if abs(payload["charge_mismatch_e"]) > 1e-4:
            raise RuntimeError(
                f"Mixed overlay final charge is not integral: {payload['total_charge_e']:.8f} e"
            )
        path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        return path
    raise RuntimeError("Mixed overlay charge audit requires a NonbondedForce")


def _build_mixed_overlay_native_reference(
    *,
    config: Any,
    output_dir: Path,
    free_polymer_seed: int | None,
    workflow_settings: ConjugatedPolymerSystemSettings,
) -> Any:
    """Build an isolated native GLYCAM reference from GLYCAM attachments only."""

    reference_config = _mixed_overlay_native_reference_config(config)
    reference_settings = ConjugatedPolymerSystemSettings(
        create_final_interchange=False,
        preserve_reference_atom_names=True,
        run_relaxation=False,
        pdb_fragment_output_mode="experimental_pablo",
        relaxation=workflow_settings.relaxation,
        conjugate_parameterization=workflow_settings.conjugate_parameterization,
    )
    reference_result = build_conjugated_polymer_system_from_config(
        reference_config,
        output_dir=output_dir,
        settings=reference_settings,
        free_polymer_seed=free_polymer_seed,
    )
    if reference_result.exact_export_bundle is None:
        raise RuntimeError("Mixed overlay native reference did not produce an exact OpenMM bundle")
    return reference_result.exact_export_bundle


def _mixed_overlay_native_reference_config(config: Any) -> Any:
    """Return config with only GLYCAM attachments enabled for native reference build."""

    resolved = resolve_conjugate_force_fields(config)
    glycam_names = {
        attachment.attachment_name for attachment in resolved.attachments if attachment.is_glycam
    }
    if not glycam_names:
        raise RuntimeError("Mixed overlay native reference requires at least one GLYCAM attachment")
    reference_config = (
        config.model_copy(deep=True) if hasattr(config, "model_copy") else copy.deepcopy(config)
    )
    conjugation = getattr(reference_config, "conjugation", None)
    if conjugation is None:
        raise RuntimeError("Mixed overlay native reference requires conjugation config")
    filtered_attachments = []
    for attachment in getattr(conjugation, "attachments", ()) or ():
        enabled = str(getattr(attachment, "name", "attachment")) in glycam_names
        if hasattr(attachment, "model_copy"):
            filtered_attachments.append(
                attachment.model_copy(update={"enabled": enabled}, deep=True)
            )
        else:
            cloned = copy.deepcopy(attachment)
            cloned.enabled = enabled
            filtered_attachments.append(cloned)
    reference_config.conjugation = conjugation.model_copy(
        update={"attachments": filtered_attachments}, deep=True
    )
    polymers = getattr(reference_config, "polymers", None)
    if polymers is not None and hasattr(polymers, "model_copy"):
        reference_config.polymers = polymers.model_copy(update={"enabled": False}, deep=True)
    return reference_config


def _build_native_to_baseline_atom_mapping(
    baseline_topology: Any,
    native_topology: Any,
    *,
    required_native_indices: set[int] | frozenset[int] | None = None,
    allow_scoped_fallback: bool = False,
) -> dict[int, int]:
    """Build native-reference to baseline atom mapping from stable PDB identities."""

    baseline_by_key: dict[tuple[str, str, str, str], list[int]] = {}
    baseline_scoped = _scoped_residue_atom_index(baseline_topology)["key_to_atom"]
    native_atom_to_scoped_key = _scoped_residue_atom_index(native_topology)["atom_to_key"]
    for atom in baseline_topology.atoms():
        for key in _overlay_atom_identity_keys(atom):
            baseline_by_key.setdefault(key, []).append(atom.index)
    mapping: dict[int, int] = {}
    used_baseline_indices: set[int] = set()
    unmapped: list[str] = []
    required = set(required_native_indices) if required_native_indices is not None else None
    for atom in native_topology.atoms():
        is_required = required is None or atom.index in required
        candidates: set[int] = set()
        for key in _overlay_atom_identity_keys(atom):
            candidates.update(baseline_by_key.get(key, ()))
        if len(candidates) > 1 and is_required:
            raise ValueError(
                "Ambiguous native-to-baseline overlay atom mapping for "
                f"{_openmm_atom_label(atom)}: candidates={sorted(candidates)}"
            )
        if len(candidates) > 1:
            continue
        mapped = next(iter(candidates), None)
        if mapped is None and allow_scoped_fallback:
            scoped_key = native_atom_to_scoped_key[atom.index]
            scoped_candidates = baseline_scoped.get(scoped_key, ())
            if len(scoped_candidates) > 1:
                raise ValueError(
                    "Ambiguous scoped native-to-baseline overlay atom mapping for "
                    f"{_openmm_atom_label(atom)}: candidates={sorted(scoped_candidates)}"
                )
            mapped = scoped_candidates[0] if scoped_candidates else None
        if mapped is None:
            if is_required:
                unmapped.append(_openmm_atom_label(atom))
            continue
        if mapped in used_baseline_indices:
            if not is_required:
                continue
            raise ValueError(
                "Duplicate native-to-baseline overlay atom mapping targets baseline atom "
                f"{mapped} while mapping {_openmm_atom_label(atom)}"
            )
        mapping[atom.index] = mapped
        used_baseline_indices.add(mapped)
    if unmapped:
        raise ValueError(
            f"Native-to-baseline overlay atom mapping could not map native atoms: {unmapped[:20]}"
        )
    return mapping


def _scoped_residue_atom_index(topology: Any) -> dict[str, dict[Any, Any]]:
    """Return residue-occurrence scoped atom identity indexes."""

    residue_ordinals: dict[tuple[str, str], int] = {}
    residue_keys: dict[Any, tuple[str, str, int]] = {}
    for residue in topology.residues():
        residue_name = str(residue.name).strip().upper()
        base_key = (str(residue.chain.id).strip(), residue_name)
        residue_ordinals[base_key] = residue_ordinals.get(base_key, 0) + 1
        residue_keys[residue] = (*base_key, residue_ordinals[base_key])
    atom_to_key = {}
    key_to_atom = {}
    for atom in topology.atoms():
        scoped_key = (*residue_keys[atom.residue], str(atom.name).strip().upper())
        atom_to_key[atom.index] = scoped_key
        key_to_atom.setdefault(scoped_key, []).append(atom.index)
    return {"atom_to_key": atom_to_key, "key_to_atom": key_to_atom}


def _overlay_atom_identity_keys(atom: Any) -> tuple[tuple[str, str, str, str], ...]:
    """Return equivalent stable keys for cross-route atom mapping."""

    residue = atom.residue
    chain_id = str(residue.chain.id).strip()
    residue_id = str(residue.id).strip()
    residue_name = str(residue.name).strip().upper()
    atom_name = str(atom.name).strip().upper()
    keys = [(chain_id, residue_id, residue_name, atom_name)]
    if residue_name in {"ASX", "NLN", "ASN"}:
        residue_aliases = ("ASX", "NLN", "ASN")
        atom_aliases = ("HD21", "HD22") if atom_name in {"HD21", "HD22"} else (atom_name,)
        keys.extend(
            (chain_id, residue_id, alias, atom_alias)
            for alias in residue_aliases
            for atom_alias in atom_aliases
        )
    if residue_name in {"LYS", "LYX"}:
        keys.extend((chain_id, residue_id, alias, atom_name) for alias in ("LYS", "LYX"))
    return tuple(dict.fromkeys(keys))


def _glycam_particles_from_native_audit(
    native_topology: Any, native_audit: dict[str, Any]
) -> frozenset[int]:
    """Return native GLYCAM-owned particles from native audit domain assignments."""

    assignments = (
        native_audit.get("domain_assignments", {}).get("residues", ())
        if isinstance(native_audit.get("domain_assignments"), dict)
        else ()
    )
    labels = {
        str(entry.get("residue"))
        for entry in assignments
        if str(entry.get("domain", "")).strip().lower() in {"glycan", "protein_modified_glycam"}
    }
    if not labels:
        raise ValueError("Native GLYCAM audit contains no explicit GLYCAM-owned residue labels")
    topology_labels = {_openmm_residue_label(residue) for residue in native_topology.residues()}
    missing_labels = sorted(labels - topology_labels)
    if missing_labels:
        raise ValueError(
            "Native GLYCAM audit residue labels are absent from the native topology: "
            f"{missing_labels}"
        )
    particles = {
        atom.index
        for atom in native_topology.atoms()
        if _openmm_residue_label(atom.residue) in labels
    }
    if not particles:
        raise ValueError("Native GLYCAM audit labels resolved to zero topology particles")
    return frozenset(particles)


def _openmm_atom_label(atom: Any) -> str:
    """Return a concise OpenMM atom label for diagnostics."""

    return f"{_openmm_residue_label(atom.residue)}:{atom.name}#{atom.index}"


def _openmm_residue_label(residue: Any) -> str:
    """Return the residue label format used by native GLYCAM audit records."""

    insertion = getattr(residue, "insertionCode", "") or ""
    return f"{residue.chain.id}:{residue.name}{residue.id}{insertion}"


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
    return np.asarray(pdb_coordinates(path), dtype=float)


def _restore_relaxed_pdb_atom_names(
    construction: ModifierConstructionResult,
    atom_name_template_pdb: Path | str,
) -> None:
    """Restore reference atom-name fields on relaxed PDB artifacts when present."""
    relaxation = getattr(construction, "relaxation", None)
    if relaxation is None:
        return
    for pdb_path in (getattr(relaxation, "relaxed_pdb_path", None),):
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
    _apply_pdb_atom_identity_to_topology(topology, template_pdb_path)


def _apply_pdb_atom_identity_to_topology(topology: Any, template_pdb_path: Path | str) -> None:
    """Attach authoritative product PDB identity metadata to topology atoms.

    Parameters
    ----------
    topology : Any
        OpenFF-like topology whose atom order matches the product PDB.
    template_pdb_path : pathlib.Path or str
        Product PDB containing serial, chain, residue, and atom-name identity.

    Raises
    ------
    ValueError
        If the topology and PDB atom counts differ or any PDB atom lacks the
        identity fields required for deterministic charge transfer.
    """
    template_atoms = _pdb_atom_records(Path(template_pdb_path))
    topology_atoms = list(_iter_openff_topology_atoms(topology))
    if not topology_atoms:
        return
    if len(topology_atoms) != len(template_atoms):
        raise ValueError(
            "Topology atom count does not match product identity template count: "
            f"topology={len(topology_atoms)}, template={len(template_atoms)}"
        )
    for topology_atom, pdb_atom in zip(topology_atoms, template_atoms, strict=True):
        _validate_pdb_atom_identity(pdb_atom, template_pdb_path)
        topology_atom.name = pdb_atom.atom_name
        _update_atom_metadata(topology_atom, _pdb_atom_identity_metadata(pdb_atom))


def _update_atom_metadata(atom: Any, values: dict[str, int | str]) -> None:
    """Update OpenFF atom metadata across mutable and property-backed atoms."""
    metadata = getattr(atom, "metadata", None)
    if isinstance(metadata, MutableMapping):
        metadata.update(values)
        return
    private_metadata = getattr(atom, "_metadata", None)
    if isinstance(private_metadata, dict):
        private_metadata.update(values)
        return
    try:
        atom.metadata = dict(values)
    except AttributeError as exc:
        raise ValueError(
            "Could not attach product PDB identity metadata to topology atom; "
            f"atom type {type(atom).__name__} exposes no mutable metadata mapping"
        ) from exc


def _pdb_atom_identity_metadata(atom: PdbAtomRecord) -> dict[str, int | str]:
    """Return product PDB atom identity metadata for charge transfer."""
    metadata: dict[str, int | str] = {
        "chain_id": atom.chain_id.strip(),
        "residue_name": atom.residue_name.strip().upper(),
        "residue_number": int(atom.residue_number),
        "insertion_code": atom.insertion_code.strip(),
        "atom_name": atom.atom_name.strip(),
        "product_identity_source": "product_pdb",
    }
    if atom.atom_index is not None:
        metadata["product_atom_index"] = int(atom.atom_index)
    if atom.serial is not None:
        metadata["serial"] = int(atom.serial)
        metadata["product_atom_serial"] = int(atom.serial)
    return metadata


def _validate_pdb_atom_identity(atom: PdbAtomRecord, source_path: Path | str) -> None:
    """Fail clearly when a product PDB atom cannot identify a charge target."""
    missing = []
    if not atom.atom_name.strip():
        missing.append("atom_name")
    if not atom.residue_name.strip():
        missing.append("residue_name")
    if atom.residue_number is None:
        missing.append("residue_number")
    if atom.serial is None:
        missing.append("serial")
    if missing:
        raise ValueError(
            "Product PDB atom identity metadata is incomplete for charge template transfer: "
            f"serial={atom.serial} missing {', '.join(missing)} in {source_path}"
        )


def _pdb_atom_records(path: Path) -> tuple[PdbAtomRecord, ...]:
    return parse_structure_pdb_atom_records(path)


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
    atoms: list[Any] = []
    for molecule in getattr(topology, "molecules", ()):
        atoms.extend(getattr(molecule, "atoms", ()))
    if atoms:
        return tuple(atoms)
    if hasattr(topology, "atoms"):
        atoms_attr = topology.atoms
        atoms = list(atoms_attr() if callable(atoms_attr) else atoms_attr)
    return tuple(atoms)


def _required_int(value: int | None, name: str) -> int:
    if value is None:
        raise ValueError(f"{name} is required")
    return int(value)
