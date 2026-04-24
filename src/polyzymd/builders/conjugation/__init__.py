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
    "GraphEditPlan",
    "MoietyDescriptor",
    "MoietyReactiveGroup",
    "PabloAvailability",
    "PabloIngestionError",
    "PabloIngestor",
    "PabloIngestionResult",
    "PabloLinkCandidate",
    "PabloResidueSummary",
    "PabloStructureCounts",
    "PabloStructurePreflight",
    "ReactionMechanism",
    "SiteAtomRule",
    "save_metadata",
    "ensure_polymerist_py312_compat",
    "get_builtin_mechanism",
    "import_polymerist_building",
    "list_builtin_mechanisms",
    "load_builtin_mechanisms",
    "match_site_rule",
    "normalize_attachment_site",
    "normalize_moiety_descriptor",
    "write_diagnostics_report",
]
