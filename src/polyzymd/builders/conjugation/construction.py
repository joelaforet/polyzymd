"""Integrated modifier-linking construction workflow for conjugation POCs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.contracts import (
    ExplicitLinkageContract,
    ResolvedAttachmentPlan,
    parse_pdb_atom_records,
    placed_fragment_from_resolved_plan,
    resolve_explicit_linkage_contract,
)
from polyzymd.builders.conjugation.crosslinks import (
    CrosslinkValidationResult,
    require_pablo_crosslink_requirement,
)
from polyzymd.builders.conjugation.linkers import ModifierLinker, NhsLysModifierLinker
from polyzymd.builders.conjugation.pablo_adapter import PabloIngestionResult, PabloIngestor
from polyzymd.builders.conjugation.parameterization import (
    InterchangeParameterizationResult,
    InterchangeParameterizationSettings,
    create_interchange_from_pablo_topology,
)
from polyzymd.builders.conjugation.pdb_assembly import (
    CrosslinkedPdbAssemblyOptions,
    CrosslinkedPdbAssemblyResult,
    NhsLysPdbAttachment,
    write_crosslinked_pdb,
)
from polyzymd.builders.conjugation.placement import (
    PackmolModifierPlacementResult,
    PackmolModifierPlacementSettings,
    place_modifier_with_packmol,
    place_modifier_with_resolved_plan,
)
from polyzymd.builders.conjugation.polymer_fragment import GeneratedPolymerFragment
from polyzymd.builders.conjugation.smoke import (
    VacuumSmokeResult,
    VacuumSmokeSettings,
    run_restrained_vacuum_smoke,
)


class ModifierConstructionSettings(BaseModel):
    """Settings for integrated modifier-linking construction."""

    crosslinked_pdb_name: str = "assembled_crosslinked.pdb"
    placement: PackmolModifierPlacementSettings = Field(
        default_factory=PackmolModifierPlacementSettings
    )
    parameterization: InterchangeParameterizationSettings = Field(
        default_factory=InterchangeParameterizationSettings
    )
    smoke: VacuumSmokeSettings = Field(default_factory=VacuumSmokeSettings)
    run_smoke: bool = True


class ModifierConstructionResult(BaseModel):
    """Integrated construction result and artifact summary."""

    model_config = {"arbitrary_types_allowed": True}

    output_dir: Path
    resolved_plan: ResolvedAttachmentPlan
    crosslink_validation: CrosslinkValidationResult
    placement: PackmolModifierPlacementResult
    assembly: CrosslinkedPdbAssemblyResult
    pablo: PabloIngestionResult
    parameterization: InterchangeParameterizationResult
    smoke: VacuumSmokeResult | None = None
    crosslinked_pdb_path: Path
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)


def construct_modifier_linked_protein(
    *,
    protein_pdb_path: Path | str,
    modifier: GeneratedPolymerFragment,
    linker: ModifierLinker,
    ccd_pablo_policy: Any,
    output_dir: Path | str,
    chain_policy: Any | None = None,
    settings: ModifierConstructionSettings | None = None,
    run_packmol_func: Any | None = None,
    pablo_ingestor: PabloIngestor | None = None,
    parameterizer: Any | None = None,
    smoke_runner: Any | None = None,
) -> ModifierConstructionResult:
    """Construct a modifier-linked protein through the POC workflow.

    Parameters
    ----------
    protein_pdb_path : pathlib.Path or str
        Protein PDB containing the target linker atom.
    modifier : GeneratedPolymerFragment
        Generated modifier or polymer fragment to attach.
    linker : ModifierLinker
        Linker strategy. ``NhsLysModifierLinker`` is the first concrete
        realization supported by the PDB assembly step.
    ccd_pablo_policy : Any
        Explicit Pablo policy with a matching ``crosslinks`` entry.
    output_dir : pathlib.Path or str
        Directory for placement, assembly, Pablo, parameterization, and smoke artifacts.
    chain_policy : Any or None, optional
        Chain policy forwarded to Pablo ingestion, by default ``None``.
    settings : ModifierConstructionSettings or None, optional
        Workflow settings, by default ``None``.
    run_packmol_func : callable or None, optional
        Optional Packmol executor injection for tests, by default ``None``.
    pablo_ingestor : PabloIngestor or None, optional
        Optional Pablo ingestor injection for tests, by default ``None``.
    parameterizer : callable or None, optional
        Optional parameterization function injection for tests, by default ``None``.
    smoke_runner : callable or None, optional
        Optional smoke runner injection for tests, by default ``None``.

    Returns
    -------
    ModifierConstructionResult
        Integrated workflow result and artifacts.
    """
    construction_settings = settings or ModifierConstructionSettings()
    artifact_dir = Path(output_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)

    resolved_plan = _resolved_plan_from_linker(linker, protein_pdb_path, modifier)
    crosslink_validation = require_pablo_crosslink_requirement(
        ccd_pablo_policy,
        resolved_plan.pablo_crosslink_requirement,
    )

    placement_result = place_modifier_with_packmol(
        protein_pdb_path,
        modifier,
        linker,
        artifact_dir,
        settings=construction_settings.placement,
        run_packmol_func=run_packmol_func,
    )

    crosslinked_pdb_path = artifact_dir / construction_settings.crosslinked_pdb_name
    attachment = _attachment_from_plan_or_linker(resolved_plan, linker, protein_pdb_path)
    placed_modifier = placed_fragment_from_resolved_plan(
        placement_result.placed_modifier,
        resolved_plan,
    )
    assembly_result = write_crosslinked_pdb(
        protein_pdb_path,
        placed_modifier,
        attachment,
        crosslinked_pdb_path,
        CrosslinkedPdbAssemblyOptions(),
    )

    ingestor = pablo_ingestor or PabloIngestor(policy=ccd_pablo_policy)
    pablo_result = ingestor.ingest_structure(
        crosslinked_pdb_path,
        chain_policy=chain_policy,
        output_dir=artifact_dir,
    )
    if not pablo_result.success or pablo_result.topology is None:
        raise RuntimeError(_pablo_failure_message(pablo_result))

    parameterizer_func = parameterizer or create_interchange_from_pablo_topology
    parameterization_result = parameterizer_func(
        pablo_result.topology,
        settings=construction_settings.parameterization,
    )
    if not parameterization_result.success or parameterization_result.interchange is None:
        raise RuntimeError("OpenFF Interchange parameterization did not produce an interchange")

    smoke_result = None
    if construction_settings.run_smoke:
        smoke_func = smoke_runner or run_restrained_vacuum_smoke
        smoke_result = smoke_func(
            parameterization_result.interchange,
            artifact_dir,
            settings=construction_settings.smoke,
        )
        if not smoke_result.success:
            raise RuntimeError("OpenMM restrained vacuum smoke did not report success")

    return ModifierConstructionResult(
        output_dir=artifact_dir,
        resolved_plan=resolved_plan,
        crosslink_validation=crosslink_validation,
        placement=placement_result,
        assembly=assembly_result,
        pablo=pablo_result,
        parameterization=parameterization_result,
        smoke=smoke_result,
        crosslinked_pdb_path=crosslinked_pdb_path,
        diagnostics=("Modifier-linked protein construction POC completed",),
    )


def construct_explicit_pdb_linkage(
    *,
    protein_pdb_path: Path | str,
    modifier_pdb_path: Path | str,
    contract: ExplicitLinkageContract,
    ccd_pablo_policy: Any,
    output_dir: Path | str,
    chain_policy: Any | None = None,
    settings: ModifierConstructionSettings | None = None,
    run_packmol_func: Any | None = None,
    pablo_ingestor: PabloIngestor | None = None,
    parameterizer: Any | None = None,
    smoke_runner: Any | None = None,
) -> ModifierConstructionResult:
    """Construct a protein-modifier conjugate from an explicit PDB contract.

    Parameters
    ----------
    protein_pdb_path : pathlib.Path or str
        Protein PDB containing the selected protein endpoint atom.
    modifier_pdb_path : pathlib.Path or str
        Modifier PDB containing the selected modifier endpoint atom.
    contract : ExplicitLinkageContract
        PDB-only explicit linkage contract.
    ccd_pablo_policy : Any
        Explicit Pablo policy with a matching ``crosslinks`` entry.
    output_dir : pathlib.Path or str
        Directory for placement, assembly, Pablo, parameterization, and smoke artifacts.
    chain_policy : Any or None, optional
        Chain policy forwarded to Pablo ingestion, by default ``None``.
    settings : ModifierConstructionSettings or None, optional
        Workflow settings, by default ``None``.
    run_packmol_func : callable or None, optional
        Optional Packmol executor injection for tests, by default ``None``.
    pablo_ingestor : PabloIngestor or None, optional
        Optional Pablo ingestor injection for tests, by default ``None``.
    parameterizer : callable or None, optional
        Optional parameterization function injection for tests, by default ``None``.
    smoke_runner : callable or None, optional
        Optional smoke runner injection for tests, by default ``None``.

    Returns
    -------
    ModifierConstructionResult
        Integrated workflow result and artifacts.
    """
    construction_settings = settings or ModifierConstructionSettings()
    artifact_dir = Path(output_dir)
    artifact_dir.mkdir(parents=True, exist_ok=True)

    modifier_path = Path(modifier_pdb_path)
    resolved_plan = resolve_explicit_linkage_contract(
        protein_pdb_path,
        modifier_path,
        contract,
    )
    crosslink_validation = require_pablo_crosslink_requirement(
        ccd_pablo_policy,
        resolved_plan.pablo_crosslink_requirement,
    )

    modifier = _generated_fragment_from_explicit_pdb(modifier_path, resolved_plan)
    placement_result = place_modifier_with_resolved_plan(
        protein_pdb_path,
        modifier,
        resolved_plan,
        artifact_dir,
        settings=construction_settings.placement,
        run_packmol_func=run_packmol_func,
    )

    crosslinked_pdb_path = artifact_dir / construction_settings.crosslinked_pdb_name
    placed_modifier = placed_fragment_from_resolved_plan(
        placement_result.placed_modifier,
        resolved_plan,
    )
    assembly_result = write_crosslinked_pdb(
        protein_pdb_path,
        placed_modifier,
        resolved_plan.to_pdb_linkage_attachment(),
        crosslinked_pdb_path,
        CrosslinkedPdbAssemblyOptions(),
    )

    ingestor = pablo_ingestor or PabloIngestor(policy=ccd_pablo_policy)
    pablo_result = ingestor.ingest_structure(
        crosslinked_pdb_path,
        chain_policy=chain_policy,
        output_dir=artifact_dir,
    )
    if not pablo_result.success or pablo_result.topology is None:
        raise RuntimeError(_pablo_failure_message(pablo_result))

    parameterizer_func = parameterizer or create_interchange_from_pablo_topology
    parameterization_result = parameterizer_func(
        pablo_result.topology,
        settings=construction_settings.parameterization,
    )
    if not parameterization_result.success or parameterization_result.interchange is None:
        raise RuntimeError("OpenFF Interchange parameterization did not produce an interchange")

    smoke_result = None
    if construction_settings.run_smoke:
        smoke_func = smoke_runner or run_restrained_vacuum_smoke
        smoke_result = smoke_func(
            parameterization_result.interchange,
            artifact_dir,
            settings=construction_settings.smoke,
        )
        if not smoke_result.success:
            raise RuntimeError("OpenMM restrained vacuum smoke did not report success")

    return ModifierConstructionResult(
        output_dir=artifact_dir,
        resolved_plan=resolved_plan,
        crosslink_validation=crosslink_validation,
        placement=placement_result,
        assembly=assembly_result,
        pablo=pablo_result,
        parameterization=parameterization_result,
        smoke=smoke_result,
        crosslinked_pdb_path=crosslinked_pdb_path,
        diagnostics=("Explicit PDB linkage construction completed",),
    )


def _resolved_plan_from_linker(
    linker: ModifierLinker,
    protein_pdb_path: Path | str,
    modifier: GeneratedPolymerFragment,
) -> ResolvedAttachmentPlan:
    """Resolve a generic attachment plan from a linker."""
    plan_builder = getattr(linker, "resolve_plan", None)
    if callable(plan_builder):
        return plan_builder(protein_pdb_path, modifier)
    raise TypeError("The integrated POC requires a linker with a generic resolved plan")


def _generated_fragment_from_explicit_pdb(
    modifier_pdb_path: Path,
    plan: ResolvedAttachmentPlan,
) -> GeneratedPolymerFragment:
    """Parse a modifier PDB into a generated fragment using resolved selectors."""
    atoms = parse_pdb_atom_records(modifier_pdb_path)
    bonds = _conect_bonds_from_pdb(modifier_pdb_path)
    return GeneratedPolymerFragment.from_atom_records(
        atoms,
        bonds=bonds,
        reactive_atom_serial=plan.modifier_link_atom.serial,
        reactive_atom_index=plan.modifier_link_atom.atom_index,
        reactive_atom_name=None,
        leaving_atom_serials=tuple(
            atom.serial for atom in plan.modifier_leaving_atoms if atom.serial is not None
        ),
        leaving_atom_indices=tuple(
            atom.atom_index for atom in plan.modifier_leaving_atoms if atom.atom_index is not None
        ),
        leaving_atom_names=(),
        name=modifier_pdb_path.stem,
    )


def _conect_bonds_from_pdb(path: Path) -> tuple[tuple[int, int], ...]:
    """Parse unique serial-based CONECT bonds from a PDB file."""
    bonds: set[tuple[int, int]] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("CONECT"):
                continue
            source = _parse_conect_serial(line[6:11])
            if source is None:
                continue
            for start in range(11, len(line), 5):
                target = _parse_conect_serial(line[start : start + 5])
                if target is None or target == source:
                    continue
                bonds.add(tuple(sorted((source, target))))
    return tuple(sorted(bonds))


def _parse_conect_serial(value: str) -> int | None:
    """Parse a CONECT serial field."""
    stripped = value.strip()
    if not stripped:
        return None
    try:
        return int(stripped)
    except ValueError:
        return None


def _attachment_from_plan_or_linker(
    plan: ResolvedAttachmentPlan, linker: ModifierLinker, protein_pdb_path: Path | str
) -> NhsLysPdbAttachment:
    """Build the PDB assembly attachment for supported linkers."""
    try:
        return plan.to_nhs_lys_pdb_attachment()
    except ValueError:
        pass
    if isinstance(linker, NhsLysModifierLinker):
        return linker.attachment(protein_pdb_path)
    attachment_builder = getattr(linker, "attachment", None)
    if callable(attachment_builder):
        attachment = attachment_builder(protein_pdb_path)
        if isinstance(attachment, NhsLysPdbAttachment):
            return attachment
    raise TypeError(
        "The integrated POC currently requires a linker that can provide "
        "NhsLysPdbAttachment for write_crosslinked_pdb()"
    )


def _pablo_failure_message(result: PabloIngestionResult) -> str:
    """Return a compact Pablo ingestion failure message."""
    diagnostics = [f"{diag.code}: {diag.message}" for diag in result.diagnostics]
    joined = "; ".join(diagnostics) if diagnostics else "no diagnostics were reported"
    return f"Pablo ingestion failed or returned no topology for {result.path}: {joined}"
