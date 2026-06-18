"""Covalent modification builder skeleton."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Mapping

from polyzymd.builders.conjugation.diagnostics import (
    ConjugationDiagnosticsReport,
    DiagnosticCode,
    DiagnosticSeverity,
    write_diagnostics_report,
)
from polyzymd.builders.conjugation.exceptions import (
    ConjugationError,
    ConjugationNotImplementedError,
    PabloIngestionError,
)
from polyzymd.builders.conjugation.execution import (
    ExplicitNhsReactiveGroup,
    RdkitGraphEditExecutionRequest,
    RdkitGraphEditExecutionResult,
    RdkitGraphEditExecutionSummary,
    extract_explicit_rdkit_execution_request,
)
from polyzymd.builders.conjugation.mechanism_library import get_builtin_mechanism
from polyzymd.builders.conjugation.metadata import (
    ConjugationMetadata,
    chain_policy_from_config,
    save_metadata,
)
from polyzymd.builders.conjugation.models import ConjugationBuildResult
from polyzymd.builders.conjugation.moieties import normalize_moiety_descriptor
from polyzymd.builders.conjugation.nhs_lys import (
    LysineReactiveSite,
    NhsReactiveGroup,
    detect_nhs_reactive_group,
    execute_nhs_lys_amide_rdkit_graph_edit,
    extract_lysine_reactive_site,
    plan_nhs_lys_amide,
)
from polyzymd.builders.conjugation.pablo_adapter import PabloIngestionResult, PabloIngestor
from polyzymd.builders.conjugation.sites import (
    AttachmentSite,
    match_site_rule,
    normalize_attachment_site,
)
from polyzymd.config.schema import ConjugationConfig, ConjugationMode

LOGGER = logging.getLogger(__name__)


class CovalentModificationBuilder:
    """No-op/skeleton builder for covalent modifications.

    Phase 0-1 establishes configuration, metadata, diagnostics, and the Pablo
    adapter boundary. Full chemistry graph surgery and parameterization are
    intentionally deferred to later phases.
    """

    def __init__(self, config: ConjugationConfig | None, output_dir: Path | None = None) -> None:
        """Initialize the covalent modification builder.

        Parameters
        ----------
        config : ConjugationConfig or None
            Top-level covalent modification config. ``None`` is treated as
            disabled for backwards compatibility.
        output_dir : Path or None, optional
            Working directory for metadata and diagnostics sidecars, by default
            ``None``.
        """
        self.config = config or ConjugationConfig(enabled=False)
        self.output_dir = output_dir

    @classmethod
    def from_config(
        cls,
        config: Any,
        output_dir: Path | str | None = None,
    ) -> "CovalentModificationBuilder":
        """Create a builder from a simulation or conjugation config.

        Parameters
        ----------
        config : Any
            Either a :class:`SimulationConfig` with a ``conjugation`` attribute
            or a :class:`ConjugationConfig` directly.
        output_dir : Path, str, or None, optional
            Working directory for sidecar files, by default ``None``.

        Returns
        -------
        CovalentModificationBuilder
            Configured builder instance.
        """
        conjugation_config = getattr(config, "conjugation", config)
        output_path = Path(output_dir) if output_dir is not None else None
        return cls(conjugation_config, output_dir=output_path)

    def build(
        self,
        topology: Any,
        *,
        protein_topology: Any | None = None,
        context: Mapping[str, Any] | None = None,
    ) -> ConjugationBuildResult:
        """Apply covalent modifications to the current pre-solvation topology.

        Parameters
        ----------
        topology : Any
            Current OpenFF topology or topology-like object.
        protein_topology : Any or None, optional
            Protein topology context for future construction workflows, by
            default ``None``.
        context : Mapping[str, Any] or None, optional
            Additional build context supplied by :class:`SystemBuilder`, by
            default ``None``.

        Returns
        -------
        ConjugationBuildResult
            The unchanged topology plus empty metadata and diagnostics when the
            workflow is disabled.

        Raises
        ------
        ConjugationNotImplementedError
            If an enabled Phase 0-1 workflow requests chemistry that is not yet
            implemented.
        """
        if not self.config.enabled:
            LOGGER.debug("Covalent modification disabled; returning topology unchanged")
            return ConjugationBuildResult(
                topology=topology,
                metadata=ConjugationMetadata(),
                diagnostics=ConjugationDiagnosticsReport(),
            )

        report = self._build_enabled_report()
        metadata = self._build_enabled_metadata()
        self._write_sidecars(metadata, report)
        try:
            execution_request = extract_explicit_rdkit_execution_request(context)
        except ValueError as exc:
            report.add(
                DiagnosticCode.GRAPH_EDIT_EXECUTION,
                "Explicit RDKit graph edit context validation failed",
                severity=DiagnosticSeverity.ERROR,
                details={"error": str(exc)},
            )
            self._write_sidecars(metadata, report)
            raise ConjugationError(str(exc)) from exc

        if execution_request is not None and self.config.mode != ConjugationMode.CONSTRUCT:
            report.add(
                DiagnosticCode.GRAPH_EDIT_EXECUTION,
                "Explicit RDKit graph edit execution requires construct mode",
                severity=DiagnosticSeverity.ERROR,
                details={"mode": self.config.mode.value},
            )
            self._write_sidecars(metadata, report)
            raise ConjugationError(
                "Explicit RDKit graph edit execution requires conjugation mode 'construct'"
            )

        if _is_source_backed_config(self.config):
            ingestion = _cached_pablo_ingestion(context)
            if ingestion is None:
                ingestor = PabloIngestor(self.config.ccd_pablo)
                try:
                    ingestion = ingestor.ingest_structure(
                        self.config.source_pdb_path,
                        chain_policy=self.config.chain_policy,
                        output_dir=self.output_dir,
                    )
                except PabloIngestionError as exc:
                    report.add(
                        DiagnosticCode.PABLO_INGESTION,
                        "Pablo ingestion input validation failed",
                        severity=DiagnosticSeverity.ERROR,
                        details={"error": str(exc)},
                    )
                    self._write_sidecars(metadata, report)
                    raise ConjugationNotImplementedError(
                        "Pablo ingestion could not produce a usable topology because the source "
                        f"structure path is invalid: {exc}"
                    ) from exc

            planned_attachments: list[dict[str, Any]] = []
            if self.config.attachments and self.config.mode in {
                ConjugationMode.CONSTRUCT,
                ConjugationMode.MIXED,
            }:
                try:
                    planned_attachments = self._validate_construct_plan(report)
                except ConjugationNotImplementedError:
                    self._write_sidecars(metadata, report)
                    raise
            elif self.config.attachments:
                planned_attachments = [
                    attachment.model_dump(mode="json") for attachment in self.config.attachments
                ]

            metadata = ingestion.metadata.model_copy(
                update={
                    "attachments": planned_attachments,
                    "notes": [
                        *ingestion.metadata.notes,
                        "Phase 2 Pablo ingestion did not run Interchange parameterization or solvation",
                        "Source-backed attachment construction is recorded for a later chemistry phase",
                    ],
                }
            )
            report.diagnostics.extend(ingestion.diagnostics)
            self._write_sidecars(metadata, report)
            if ingestion.success and ingestion.topology is not None:
                return ConjugationBuildResult(
                    topology=ingestion.topology,
                    metadata=metadata,
                    diagnostics=report,
                )

            report.add(
                DiagnosticCode.UNSUPPORTED_OPERATION,
                "Pablo ingestion did not return a usable topology for downstream parameterization",
                severity=DiagnosticSeverity.ERROR,
                details={"pablo_success": ingestion.success},
            )
            self._write_sidecars(metadata, report)
            raise ConjugationNotImplementedError(
                "Pablo ingestion was attempted but did not return a usable topology. Review "
                "conjugation_diagnostics.json and pablo_ingestion_result.json for parser details."
            )

        if self.config.mode == ConjugationMode.CONSTRUCT:
            try:
                planned_attachments = self._validate_construct_plan(report)
            except ConjugationNotImplementedError:
                self._write_sidecars(metadata, report)
                raise
            metadata = metadata.model_copy(
                update={
                    "attachments": planned_attachments,
                    "notes": [
                        *metadata.notes,
                        "Phase 3 construct planning validated mechanisms, sites, and moieties",
                        "Graph surgery and parameterization are intentionally deferred",
                    ],
                }
            )
            if execution_request is not None:
                execution = self._execute_explicit_rdkit_graph_edit(
                    execution_request,
                    report,
                    metadata,
                )
                metadata = metadata.model_copy(
                    update={
                        "notes": [
                            *metadata.notes,
                            "Explicit NHS-Lys RDKit graph edit executed in memory",
                            "ConjugationBuildResult.topology remains unchanged",
                        ],
                    }
                )
                self._write_sidecars(metadata, report)
                return ConjugationBuildResult(
                    topology=topology,
                    metadata=metadata,
                    diagnostics=report,
                    graph_edit_results=[execution],
                    graph_edit_summaries=[execution.summary],
                )
            report.add(
                DiagnosticCode.UNSUPPORTED_OPERATION,
                "Construct mode validation completed, but covalent graph surgery is not implemented",
                severity=DiagnosticSeverity.ERROR,
                details={"planned_attachments": len(planned_attachments)},
            )
            self._write_sidecars(metadata, report)
            raise ConjugationNotImplementedError(
                "Construct mode validated requested attachments, but covalent graph surgery is not "
                "implemented in this phase. Review conjugation diagnostics for the validated plan."
            )

        report.add(
            DiagnosticCode.UNSUPPORTED_OPERATION,
            f"Conjugation mode '{self.config.mode.value}' is not implemented in Phase 0-1",
            severity=DiagnosticSeverity.ERROR,
        )
        self._write_sidecars(metadata, report)
        raise ConjugationNotImplementedError(
            f"Conjugation mode '{self.config.mode.value}' is not implemented in the Phase 0-1 "
            "skeleton. This pass only adds configuration, metadata, diagnostics, and the "
            "Pablo adapter boundary."
        )

    def _execute_explicit_rdkit_graph_edit(
        self,
        request: RdkitGraphEditExecutionRequest,
        report: ConjugationDiagnosticsReport,
        metadata: ConjugationMetadata,
    ) -> RdkitGraphEditExecutionResult:
        """Execute the single supported explicit RDKit graph edit path.

        Parameters
        ----------
        request : RdkitGraphEditExecutionRequest
            Explicit in-memory molecules and optional atom-index overrides.
        report : ConjugationDiagnosticsReport
            Diagnostics report updated with graph edit execution details.
        metadata : ConjugationMetadata
            Current metadata, used only for serializable summary context.

        Returns
        -------
        RdkitGraphEditExecutionResult
            In-memory RDKit result plus a JSON-safe summary.

        Raises
        ------
        ConjugationError
            If the explicit execution request is incompatible with the controlled
            single NHS-Lys execution boundary.
        """
        enabled_attachments = [
            attachment for attachment in self.config.attachments if attachment.enabled
        ]
        if len(enabled_attachments) != 1:
            message = "Explicit RDKit graph edit execution requires exactly one enabled attachment"
            report.add(
                DiagnosticCode.GRAPH_EDIT_EXECUTION,
                message,
                severity=DiagnosticSeverity.ERROR,
                details={"enabled_attachments": len(enabled_attachments)},
            )
            self._write_sidecars(metadata, report)
            raise ConjugationError(message)

        attachment = enabled_attachments[0]
        mechanism_name = attachment.mechanism.name.strip().lower()
        if mechanism_name != "nhs_lys_amide":
            message = "Explicit RDKit graph edit execution requires mechanism 'nhs_lys_amide'"
            report.add(
                DiagnosticCode.GRAPH_EDIT_EXECUTION,
                message,
                severity=DiagnosticSeverity.ERROR,
                details={"attachment": attachment.name, "mechanism": mechanism_name},
            )
            self._write_sidecars(metadata, report)
            raise ConjugationError(message)

        try:
            mechanism = get_builtin_mechanism(mechanism_name)
            site = normalize_attachment_site(attachment.site)
            match_site_rule(site, mechanism)
            moiety = normalize_moiety_descriptor(attachment.moiety)
            lys_site = self._resolve_nhs_lys_site(site, request)
            reactive_group = self._resolve_nhs_reactive_group(request)
            plan = plan_nhs_lys_amide(
                lys_site,
                reactive_group,
                site_hydrogen_indices_to_remove=request.explicit_site_hydrogen_indices,
            )
            result = execute_nhs_lys_amide_rdkit_graph_edit(
                protein_mol=request.protein_mol,
                moiety_mol=request.moiety_mol,
                site_atom_index=plan.site_atom_index,
                reactive_carbon_index=plan.reactive_carbon_index,
                leaving_atom_indices=plan.leaving_group_atom_indices,
                site_hydrogen_indices=plan.site_hydrogen_indices_to_remove,
                sanitize=request.sanitize,
            )
        except Exception as exc:
            report.add(
                DiagnosticCode.GRAPH_EDIT_EXECUTION,
                "Explicit NHS-Lys RDKit graph edit execution failed",
                severity=DiagnosticSeverity.ERROR,
                details={
                    "attachment": attachment.name,
                    "mechanism": mechanism_name,
                    "error": str(exc),
                },
            )
            self._write_sidecars(metadata, report)
            raise ConjugationError(
                f"Explicit NHS-Lys RDKit graph edit execution failed: {exc}"
            ) from exc

        warnings = (
            *lys_site.warnings,
            *reactive_group.diagnostics,
            *plan.warnings,
            *result.warnings,
        )
        summary = RdkitGraphEditExecutionSummary(
            attachment=attachment.name,
            mechanism=mechanism.identifier,
            site=site.model_dump(mode="json"),
            moiety=moiety.model_dump(mode="json"),
            added_bond=result.added_bond,
            removed_protein_atom_indices=result.removed_protein_atom_indices,
            removed_moiety_atom_indices=result.removed_moiety_atom_indices,
            removed_atoms_count=(
                len(result.removed_protein_atom_indices) + len(result.removed_moiety_atom_indices)
            ),
            product_atom_count=result.product_mol.GetNumAtoms(),
            warnings=warnings,
            topology_unchanged=True,
        )
        report.add(
            DiagnosticCode.GRAPH_EDIT_EXECUTION,
            "Executed explicit NHS-Lys RDKit graph edit; topology remains unchanged",
            details={
                "attachment": attachment.name,
                "mechanism": mechanism.identifier,
                "site": site.model_dump(mode="json"),
                "moiety": moiety.model_dump(mode="json"),
                "added_bond": result.added_bond.model_dump(mode="json"),
                "removed_atoms_count": summary.removed_atoms_count,
                "removed_protein_atom_indices": list(result.removed_protein_atom_indices),
                "removed_moiety_atom_indices": list(result.removed_moiety_atom_indices),
                "warnings": list(warnings),
                "topology_unchanged": True,
            },
        )
        return RdkitGraphEditExecutionResult(
            plan=plan,
            graph_edit_result=result,
            summary=summary,
        )

    def _resolve_nhs_lys_site(
        self,
        site: AttachmentSite,
        request: RdkitGraphEditExecutionRequest,
    ) -> LysineReactiveSite:
        """Resolve a lysine NZ site from explicit indices or topology metadata.

        Parameters
        ----------
        site : AttachmentSite
            Normalized declared attachment site.
        request : RdkitGraphEditExecutionRequest
            Explicit execution request with optional site index overrides.

        Returns
        -------
        LysineReactiveSite
            Resolved NHS-Lys reactive site.
        """
        if request.explicit_site_atom_index is not None:
            return LysineReactiveSite(
                chain_id=site.chain_id,
                residue_name=site.residue_name,
                residue_number=site.residue_number,
                atom_name=site.atom_name,
                nz_atom_index=request.explicit_site_atom_index,
                nz_hydrogen_indices=request.explicit_site_hydrogen_indices or (),
                evidence={"source": "explicit_execution_context"},
            )

        if request.protein_topology_atoms is None:
            raise ValueError(
                "Explicit RDKit graph edit execution requires either explicit_site_atom_index "
                "or protein_topology_atoms for lysine NZ extraction"
            )

        return extract_lysine_reactive_site(
            site,
            request.protein_topology_atoms,
            bonds=request.protein_topology_bonds,
            positions=request.protein_topology_positions,
        )

    def _resolve_nhs_reactive_group(
        self,
        request: RdkitGraphEditExecutionRequest,
    ) -> NhsReactiveGroup:
        """Resolve an NHS ester reactive group from explicit indices or autodetection.

        Parameters
        ----------
        request : RdkitGraphEditExecutionRequest
            Explicit execution request with optional NHS group indices.

        Returns
        -------
        NhsReactiveGroup
            Resolved NHS reactive group.
        """
        if isinstance(request.explicit_nhs_group, NhsReactiveGroup):
            return request.explicit_nhs_group
        if isinstance(request.explicit_nhs_group, ExplicitNhsReactiveGroup):
            return request.explicit_nhs_group.to_reactive_group(request.moiety_mol)

        return detect_nhs_reactive_group(
            request.moiety_mol,
            candidate_atom_indices=request.nhs_candidate_atom_indices,
        )

    def _validate_construct_plan(
        self, report: ConjugationDiagnosticsReport
    ) -> list[dict[str, Any]]:
        """Validate declarative construct-mode attachments before graph surgery.

        Parameters
        ----------
        report : ConjugationDiagnosticsReport
            Diagnostics report to append validation events to.

        Returns
        -------
        list[dict[str, Any]]
            Serializable planned attachment records.

        Raises
        ------
        ConjugationNotImplementedError
            If a requested mechanism, site, or moiety declaration is invalid for
            this planning phase.
        """
        planned_attachments: list[dict[str, Any]] = []
        enabled_attachments = [
            attachment for attachment in self.config.attachments if attachment.enabled
        ]
        if not enabled_attachments:
            report.add(
                DiagnosticCode.MECHANISM_VALIDATION,
                "Construct mode has no enabled attachments to validate",
                severity=DiagnosticSeverity.WARNING,
            )
            return planned_attachments

        for attachment in enabled_attachments:
            mechanism_name = attachment.mechanism.name
            try:
                mechanism = get_builtin_mechanism(mechanism_name)
            except KeyError as exc:
                report.add(
                    DiagnosticCode.MECHANISM_VALIDATION,
                    f"Unknown conjugation mechanism '{mechanism_name}'",
                    severity=DiagnosticSeverity.ERROR,
                    details={"attachment": attachment.name, "mechanism": mechanism_name},
                )
                raise ConjugationNotImplementedError(
                    f"Unknown conjugation mechanism '{mechanism_name}'. Add a mechanism definition "
                    "before construct-mode graph surgery can proceed."
                ) from exc

            try:
                site = normalize_attachment_site(attachment.site)
                site_rule = match_site_rule(site, mechanism)
            except (ConjugationNotImplementedError, ValueError) as exc:
                report.add(
                    DiagnosticCode.SITE_SELECTION,
                    f"Attachment site validation failed for '{attachment.name}'",
                    severity=DiagnosticSeverity.ERROR,
                    details={"attachment": attachment.name, "error": str(exc)},
                )
                raise ConjugationNotImplementedError(
                    f"Attachment site validation failed for '{attachment.name}': {exc}"
                ) from exc

            try:
                moiety = normalize_moiety_descriptor(attachment.moiety)
            except ValueError as exc:
                report.add(
                    DiagnosticCode.MOIETY_NORMALIZATION,
                    f"Moiety normalization failed for '{attachment.name}'",
                    severity=DiagnosticSeverity.ERROR,
                    details={"attachment": attachment.name, "error": str(exc)},
                )
                raise ConjugationNotImplementedError(
                    f"Moiety normalization failed for '{attachment.name}': {exc}"
                ) from exc

            report.add(
                DiagnosticCode.MECHANISM_VALIDATION,
                f"Mechanism '{mechanism.identifier}' validated for attachment '{attachment.name}'",
                details={
                    "attachment": attachment.name,
                    "mechanism": mechanism.model_dump(mode="json"),
                },
            )
            report.add(
                DiagnosticCode.SITE_SELECTION,
                f"Explicit site validated for attachment '{attachment.name}'",
                details={
                    "attachment": attachment.name,
                    "site": site.model_dump(mode="json"),
                    "matched_rule": site_rule.model_dump(mode="json"),
                },
            )
            report.add(
                DiagnosticCode.MOIETY_NORMALIZATION,
                f"Moiety descriptor normalized for attachment '{attachment.name}'",
                details={
                    "attachment": attachment.name,
                    "moiety": moiety.model_dump(mode="json"),
                },
            )
            primitive_support = mechanism.identifier == "nhs_lys_amide"
            if primitive_support:
                report.add(
                    DiagnosticCode.MECHANISM_VALIDATION,
                    "NHS-Lys RDKit graph edit primitive is available behind an explicit-object "
                    "executor boundary",
                    details={
                        "attachment": attachment.name,
                        "mechanism": mechanism.identifier,
                        "executor": "execute_nhs_lys_amide_rdkit_graph_edit",
                        "requires_explicit_rdkit_molecules": True,
                    },
                )
            planned_attachments.append(
                {
                    "attachment": attachment.name,
                    "mechanism": mechanism.model_dump(mode="json"),
                    "site": site.model_dump(mode="json"),
                    "site_rule": site_rule.model_dump(mode="json"),
                    "moiety": moiety.model_dump(mode="json"),
                    "placement": attachment.placement.model_dump(mode="json"),
                    "config_overrides": attachment.mechanism.model_dump(mode="json"),
                    "executable_primitives": {
                        "nhs_lys_rdkit_graph_edit": primitive_support,
                        "invoked_from_config": False,
                    },
                }
            )

        return planned_attachments

    def _build_enabled_report(self) -> ConjugationDiagnosticsReport:
        """Create a diagnostics report for an enabled skeleton run."""
        report = ConjugationDiagnosticsReport(enabled=True, mode=self.config.mode.value)
        report.add(
            DiagnosticCode.ENABLED_MODE,
            f"Covalent modification requested with mode '{self.config.mode.value}'",
            details={"attachments": len(self.config.attachments)},
        )
        report.add(
            DiagnosticCode.CCD_POLICY,
            "CCD/Pablo policy recorded for future ingestion",
            details={
                "enabled": self.config.ccd_pablo.enabled,
                "lookup_policy": self.config.ccd_pablo.lookup_policy.value,
                "use_canonical_atom_names": self.config.ccd_pablo.use_canonical_atom_names,
            },
        )
        report.add(
            DiagnosticCode.CHARGE_PARAMETERIZATION,
            "Local junction charge patching is reserved for future chemistry phases",
            severity=DiagnosticSeverity.WARNING,
            details={
                "local_junction_patching": self.config.charge.local_junction_patching,
                "patch_radius_bonds": self.config.charge.patch_radius_bonds,
                "preserve_total_charge": self.config.charge.preserve_total_charge,
            },
        )
        return report

    def _build_enabled_metadata(self) -> ConjugationMetadata:
        """Create metadata for an enabled skeleton run."""
        return ConjugationMetadata(
            enabled=True,
            mode=self.config.mode.value,
            chain_policy=chain_policy_from_config(self.config.chain_policy),
            attachments=[
                attachment.model_dump(mode="json") for attachment in self.config.attachments
            ],
            notes=["Phase 0-1 skeleton did not modify topology or chemistry"],
        )

    def _write_sidecars(
        self,
        metadata: ConjugationMetadata,
        report: ConjugationDiagnosticsReport,
    ) -> None:
        """Write metadata and diagnostics sidecars when configured."""
        if self.output_dir is None or not self.config.diagnostics.enabled:
            return

        save_metadata(metadata, self.output_dir / self.config.diagnostics.metadata_filename)
        write_diagnostics_report(report, self.output_dir / self.config.diagnostics.output_filename)


def _cached_pablo_ingestion(context: Mapping[str, Any] | None) -> PabloIngestionResult | None:
    """Return a cached Pablo ingestion result from build context.

    Parameters
    ----------
    context : Mapping[str, Any] or None
        Build context supplied by :class:`SystemBuilder`.

    Returns
    -------
    PabloIngestionResult or None
        Cached result when present and type-compatible, otherwise ``None``.
    """
    if context is None:
        return None
    ingestion = context.get("pablo_ingestion_result")
    if ingestion is None:
        return None
    if not isinstance(ingestion, PabloIngestionResult):
        raise ConjugationError("Cached Pablo ingestion context has an invalid result type")
    return ingestion


def _is_source_backed_config(config: ConjugationConfig) -> bool:
    """Return whether a conjugation config should load a prepared source PDB.

    Parameters
    ----------
    config : ConjugationConfig
        Conjugation workflow settings.

    Returns
    -------
    bool
        ``True`` when the workflow is ingest-existing or explicitly supplies a
        prepared source PDB path.
    """
    return config.mode == ConjugationMode.INGEST_EXISTING or config.source_pdb_path is not None
