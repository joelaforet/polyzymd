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
    ConjugationNotImplementedError,
    PabloIngestionError,
)
from polyzymd.builders.conjugation.mechanism_library import get_builtin_mechanism
from polyzymd.builders.conjugation.metadata import (
    ConjugationMetadata,
    chain_policy_from_config,
    save_metadata,
)
from polyzymd.builders.conjugation.models import ConjugationBuildResult
from polyzymd.builders.conjugation.moieties import normalize_moiety_descriptor
from polyzymd.builders.conjugation.pablo_adapter import PabloIngestor
from polyzymd.builders.conjugation.polymerist_compat import polymerist_py312_compat_status
from polyzymd.builders.conjugation.sites import match_site_rule, normalize_attachment_site
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

        if self.config.mode == ConjugationMode.INGEST_EXISTING:
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

            metadata = ingestion.metadata.model_copy(
                update={
                    "attachments": [
                        attachment.model_dump(mode="json") for attachment in self.config.attachments
                    ],
                    "notes": [
                        *ingestion.metadata.notes,
                        "Phase 2 Pablo ingestion did not run Interchange parameterization or solvation",
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

    def _validate_construct_plan(self, report: ConjugationDiagnosticsReport) -> list[dict[str, Any]]:
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
        enabled_attachments = [attachment for attachment in self.config.attachments if attachment.enabled]
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
            planned_attachments.append(
                {
                    "attachment": attachment.name,
                    "mechanism": mechanism.model_dump(mode="json"),
                    "site": site.model_dump(mode="json"),
                    "site_rule": site_rule.model_dump(mode="json"),
                    "moiety": moiety.model_dump(mode="json"),
                    "placement": attachment.placement.model_dump(mode="json"),
                    "config_overrides": attachment.mechanism.model_dump(mode="json"),
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
        shim_status = polymerist_py312_compat_status()
        if shim_status["relevant"]:
            report.add(
                DiagnosticCode.POLYMERIST_COMPAT,
                "Polymerist Python 3.12 compatibility shim is available for lazy imports",
                severity=(
                    DiagnosticSeverity.WARNING
                    if shim_status["shim_required"]
                    else DiagnosticSeverity.INFO
                ),
                details=shim_status,
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
