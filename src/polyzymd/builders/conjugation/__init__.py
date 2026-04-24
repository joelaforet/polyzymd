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
from polyzymd.builders.conjugation.pdb_assembly import (
    CrosslinkedPdbAssemblyOptions,
    CrosslinkedPdbAssemblyResult,
    NhsLysPdbAttachment,
    PdbAtomRecord,
    PlacedPolymerFragment,
    canonicalize_poc_residue_name,
    write_crosslinked_pdb,
)
from polyzymd.builders.conjugation.polymer_fragment import (
    GeneratedPolymerFragment,
    PolymerFragmentAtom,
    PolymerFragmentResidue,
)
from polyzymd.builders.conjugation.polymer_recipe import (
    POC_EGPMA_SMILES,
    POC_NHS_SMILES,
    POC_SBMA_SMILES,
    PolymeristGenerationSmokeResult,
    PolymerMonomerRecipe,
    PolymerRecipe,
    generate_polymerist_smoke_polymer,
    sbma_egpma_nhs_recipe,
)
from polyzymd.builders.conjugation.polymerist_compat import (
    ensure_polymerist_py312_compat,
    import_polymerist_building,
)
from polyzymd.builders.conjugation.polymerist_pdb import generated_fragment_from_polymerist_pdb
from polyzymd.builders.conjugation.sites import (
    AttachmentSite,
    match_site_rule,
    normalize_attachment_site,
)
from polyzymd.builders.conjugation.structure_normalization import (
    PDBChainNormalizationAction,
    PDBCleanlinessIssue,
    PDBNormalizationPlan,
    plan_pdb_chain_normalization,
    write_normalized_pdb,
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
    "CrosslinkedPdbAssemblyOptions",
    "CrosslinkedPdbAssemblyResult",
    "CovalentModificationBuilder",
    "DiagnosticCode",
    "DiagnosticSeverity",
    "ExplicitNhsReactiveGroup",
    "GeneratedPolymerFragment",
    "GraphEditPlan",
    "LysineReactiveSite",
    "MoietyDescriptor",
    "MoietyReactiveGroup",
    "NhsLysGraphEditPlan",
    "NhsLysPdbAttachment",
    "NhsReactiveGroup",
    "PabloAvailability",
    "PDBChainNormalizationAction",
    "PDBCleanlinessIssue",
    "PDBNormalizationPlan",
    "PabloIngestionError",
    "PabloIngestor",
    "PabloIngestionResult",
    "PabloLinkCandidate",
    "PabloResidueSummary",
    "PabloStructureCounts",
    "PabloStructurePreflight",
    "PdbAtomRecord",
    "PlacedPolymerFragment",
    "POC_EGPMA_SMILES",
    "POC_NHS_SMILES",
    "POC_SBMA_SMILES",
    "PolymerFragmentAtom",
    "PolymerFragmentResidue",
    "PolymeristGenerationSmokeResult",
    "PolymerMonomerRecipe",
    "PolymerRecipe",
    "ReactionMechanism",
    "RdkitGraphEditResult",
    "RdkitGraphEditExecutionRequest",
    "RdkitGraphEditExecutionResult",
    "RdkitGraphEditExecutionSummary",
    "SiteAtomRule",
    "save_metadata",
    "canonicalize_poc_residue_name",
    "ensure_polymerist_py312_compat",
    "detect_nhs_reactive_group",
    "execute_nhs_lys_amide_rdkit_graph_edit",
    "extract_explicit_rdkit_execution_request",
    "extract_lysine_reactive_site",
    "generate_polymerist_smoke_polymer",
    "generated_fragment_from_polymerist_pdb",
    "get_builtin_mechanism",
    "import_polymerist_building",
    "list_builtin_mechanisms",
    "load_builtin_mechanisms",
    "match_site_rule",
    "normalize_attachment_site",
    "normalize_moiety_descriptor",
    "plan_pdb_chain_normalization",
    "plan_nhs_lys_amide",
    "sbma_egpma_nhs_recipe",
    "write_normalized_pdb",
    "write_crosslinked_pdb",
    "write_diagnostics_report",
]
