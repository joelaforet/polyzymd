"""Internal diagnostics and metadata artifacts for conjugation workflows."""

from __future__ import annotations

import json
from enum import Enum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field


class DiagnosticSeverity(str, Enum):
    """Diagnostic severity levels."""

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


class DiagnosticCode(str, Enum):
    """Stable diagnostic codes for conjugation reports."""

    CCD_POLICY = "ccd_policy"
    ENABLED_MODE = "enabled_mode"
    UNSUPPORTED_OPERATION = "unsupported_operation"
    CHARGE_PARAMETERIZATION = "charge_parameterization"
    PABLO_ADAPTER = "pablo_adapter"
    PABLO_INGESTION = "pablo_ingestion"
    COMPONENT_METADATA = "component_metadata"
    LINK_DISCOVERY = "link_discovery"
    MECHANISM_VALIDATION = "mechanism_validation"
    SITE_SELECTION = "site_selection"
    MOIETY_NORMALIZATION = "moiety_normalization"
    GRAPH_EDIT_EXECUTION = "graph_edit_execution"
    PDB_STRUCTURE_INSPECTION = "pdb_structure_inspection"
    PDB_CHAIN_NORMALIZATION = "pdb_chain_normalization"


class ConjugationDiagnostic(BaseModel):
    """One covalent modification diagnostic event."""

    code: DiagnosticCode
    severity: DiagnosticSeverity = DiagnosticSeverity.INFO
    message: str
    details: dict[str, Any] = Field(default_factory=dict)


class ConjugationDiagnosticsReport(BaseModel):
    """Collection of covalent modification diagnostic events."""

    enabled: bool = False
    mode: str | None = None
    diagnostics: list[ConjugationDiagnostic] = Field(default_factory=list)

    def add(
        self,
        code: DiagnosticCode,
        message: str,
        *,
        severity: DiagnosticSeverity = DiagnosticSeverity.INFO,
        details: dict[str, Any] | None = None,
    ) -> None:
        """Append a diagnostic event to this report.

        Parameters
        ----------
        code : DiagnosticCode
            Stable diagnostic code.
        message : str
            Human-readable diagnostic message.
        severity : DiagnosticSeverity, optional
            Event severity, by default ``DiagnosticSeverity.INFO``.
        details : dict[str, Any] or None, optional
            Optional structured diagnostic details, by default ``None``.
        """
        self.diagnostics.append(
            ConjugationDiagnostic(
                code=code,
                severity=severity,
                message=message,
                details=details or {},
            )
        )


def write_diagnostics_report(report: ConjugationDiagnosticsReport, path: Path | str) -> Path:
    """Write a covalent modification diagnostics report as JSON.

    Parameters
    ----------
    report : ConjugationDiagnosticsReport
        Report to serialize.
    path : Path or str
        Destination JSON path.

    Returns
    -------
    Path
        Path that was written.
    """
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report.model_dump(mode="json"), indent=2) + "\n")
    return output_path


class ComponentRole(str, Enum):
    """Canonical component roles used in conjugation metadata."""

    PROTEIN = "protein"
    SUBSTRATE = "substrate"
    MOIETY = "moiety"
    POLYMER = "polymer"
    SOLVENT = "solvent"


class ChainPolicyMetadata(BaseModel):
    """Serializable chain policy used for component metadata."""

    protein_chain: str = "A"
    substrate_chain: str = "B"
    moiety_chain: str = "C"
    solvent_start_chain: str = "D"


class ComponentMetadata(BaseModel):
    """Metadata for one molecular component in a covalently modified system."""

    component_id: str
    role: ComponentRole
    chain_id: str
    atom_indices: list[int] = Field(default_factory=list)
    residue_numbers: list[int] = Field(default_factory=list)
    source: str | None = None
    details: dict[str, Any] = Field(default_factory=dict)


class ConjugationMetadata(BaseModel):
    """Top-level metadata for a covalent modification build pass."""

    enabled: bool = False
    mode: str | None = None
    chain_policy: ChainPolicyMetadata = Field(default_factory=ChainPolicyMetadata)
    components: list[ComponentMetadata] = Field(default_factory=list)
    attachments: list[dict[str, Any]] = Field(default_factory=list)
    notes: list[str] = Field(default_factory=list)


def chain_policy_from_config(config: Any | None) -> ChainPolicyMetadata:
    """Build serializable chain metadata from a chain policy config.

    Parameters
    ----------
    config : Any or None
        Conjugation chain policy config. ``None`` yields the PolyzyMD defaults.

    Returns
    -------
    ChainPolicyMetadata
        Serializable chain policy metadata.
    """
    if config is None:
        return ChainPolicyMetadata()
    return ChainPolicyMetadata(
        protein_chain=config.protein_chain,
        substrate_chain=config.substrate_chain,
        moiety_chain=config.moiety_chain,
        solvent_start_chain=config.solvent_start_chain,
    )


def save_metadata(metadata: ConjugationMetadata, path: Path | str) -> Path:
    """Save conjugation metadata as JSON.

    Parameters
    ----------
    metadata : ConjugationMetadata
        Metadata model to serialize.
    path : Path or str
        Destination JSON path.

    Returns
    -------
    Path
        Path that was written.
    """
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(metadata.model_dump(mode="json"), indent=2) + "\n")
    return output_path
