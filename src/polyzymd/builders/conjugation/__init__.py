"""Covalent modification builders for PolyzyMD."""

from polyzymd.builders.conjugation.builder import CovalentModificationBuilder
from polyzymd.builders.conjugation.diagnostics import (
    ConjugationDiagnostic,
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
from polyzymd.builders.conjugation.graph import AddedBond, RdkitGraphEditResult
from polyzymd.builders.conjugation.mechanism_library import (
    get_builtin_mechanism,
    list_builtin_mechanisms,
    load_builtin_mechanisms,
)
from polyzymd.builders.conjugation.mechanisms import (
    BondSpec,
    ChargePatchHint,
    GraphEditPlan,
    MoietyReactiveGroup,
    ReactionMechanism,
    SiteAtomRule,
)
from polyzymd.builders.conjugation.metadata import (
    ChainPolicyMetadata,
    ComponentMetadata,
    ComponentRole,
    ConjugationMetadata,
    save_metadata,
)
from polyzymd.builders.conjugation.models import ConjugationBuildResult
from polyzymd.builders.conjugation.moieties import MoietyDescriptor, normalize_moiety_descriptor
from polyzymd.builders.conjugation.nhs_lys import (
    LysineReactiveSite,
    NhsLysGraphEditPlan,
    NhsReactiveGroup,
    detect_nhs_reactive_group,
    execute_nhs_lys_amide_rdkit_graph_edit,
    extract_lysine_reactive_site,
    plan_nhs_lys_amide,
)
from polyzymd.builders.conjugation.pablo_adapter import (
    PabloAvailability,
    PabloIngestionResult,
    PabloIngestor,
    PabloLinkCandidate,
    PabloResidueSummary,
    PabloStructureCounts,
    PabloStructurePreflight,
)
from polyzymd.builders.conjugation.polymerist_compat import (
    ensure_polymerist_py312_compat,
    import_polymerist_building,
)
from polyzymd.builders.conjugation.sites import (
    AttachmentSite,
    match_site_rule,
    normalize_attachment_site,
)

__all__ = [
    "AddedBond",
    "AttachmentSite",
    "BondSpec",
    "ChainPolicyMetadata",
    "ChargePatchHint",
    "ComponentMetadata",
    "ComponentRole",
    "ConjugationBuildResult",
    "ConjugationDiagnostic",
    "ConjugationDiagnosticsReport",
    "ConjugationError",
    "ConjugationMetadata",
    "ConjugationNotImplementedError",
    "CovalentModificationBuilder",
    "DiagnosticCode",
    "DiagnosticSeverity",
    "ExplicitNhsReactiveGroup",
    "GraphEditPlan",
    "LysineReactiveSite",
    "MoietyDescriptor",
    "MoietyReactiveGroup",
    "NhsLysGraphEditPlan",
    "NhsReactiveGroup",
    "PabloAvailability",
    "PabloIngestionError",
    "PabloIngestor",
    "PabloIngestionResult",
    "PabloLinkCandidate",
    "PabloResidueSummary",
    "PabloStructureCounts",
    "PabloStructurePreflight",
    "ReactionMechanism",
    "RdkitGraphEditResult",
    "RdkitGraphEditExecutionRequest",
    "RdkitGraphEditExecutionResult",
    "RdkitGraphEditExecutionSummary",
    "SiteAtomRule",
    "save_metadata",
    "ensure_polymerist_py312_compat",
    "detect_nhs_reactive_group",
    "execute_nhs_lys_amide_rdkit_graph_edit",
    "extract_explicit_rdkit_execution_request",
    "extract_lysine_reactive_site",
    "get_builtin_mechanism",
    "import_polymerist_building",
    "list_builtin_mechanisms",
    "load_builtin_mechanisms",
    "match_site_rule",
    "normalize_attachment_site",
    "normalize_moiety_descriptor",
    "plan_nhs_lys_amide",
    "write_diagnostics_report",
]
