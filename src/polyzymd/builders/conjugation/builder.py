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
from polyzymd.builders.conjugation.exceptions import ConjugationNotImplementedError
from polyzymd.builders.conjugation.metadata import (
    ConjugationMetadata,
    chain_policy_from_config,
    save_metadata,
)
from polyzymd.builders.conjugation.models import ConjugationBuildResult
from polyzymd.builders.conjugation.pablo_adapter import PabloIngestor
from polyzymd.builders.conjugation.polymerist_compat import polymerist_py312_compat_status
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
            availability = ingestor.probe_available()
            report.add(
                DiagnosticCode.PABLO_ADAPTER,
                "OpenFF Pablo availability checked for future ingestion",
                severity=(
                    DiagnosticSeverity.INFO
                    if availability.available
                    else DiagnosticSeverity.WARNING
                ),
                details=availability.model_dump(mode="json"),
            )
            report.add(
                DiagnosticCode.UNSUPPORTED_OPERATION,
                "OpenFF Pablo ingestion is not implemented in the Phase 0-1 skeleton",
                severity=DiagnosticSeverity.ERROR,
            )
            self._write_sidecars(metadata, report)
            return ingestor.ingest_existing(self.config.source_pdb_path)

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
