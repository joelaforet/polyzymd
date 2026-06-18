"""Compatibility inventory for pre-refactor conjugation imports."""

from __future__ import annotations

import importlib
import importlib.util
from types import ModuleType

import pytest

OPTIONAL_HEAVY_IMPORT_ROOTS = {
    "MDAnalysis",
    "openff",
    "openmm",
    "pablo",
    "pdbfixer",
    "polymerist",
    "rdkit",
}


CURRENT_HELPER_MODULE_EXPORTS = {
    "polyzymd.builders.conjugation.product_pablo": (
        "ProductStatePabloDefinitionSummary",
        "ProductStatePabloLibrary",
        "build_product_state_pablo_library",
    ),
    "polyzymd.builders.conjugation.protein_preparation": (
        "ProteinCanonicalizationResult",
        "ProteinCanonicalizationSettings",
        "canonicalize_protein_hydrogens",
    ),
    "polyzymd.builders.conjugation.local_minimization": (
        "LocalMinimizationResult",
        "LocalMinimizationSettings",
        "product_state_pablo_crosslink_requirement",
        "run_post_crosslink_local_minimization",
    ),
    "polyzymd.builders.conjugation.workflow.preparation": (
        "ProteinCanonicalizationResult",
        "ProteinCanonicalizationSettings",
        "canonicalize_protein_hydrogens",
    ),
    "polyzymd.builders.conjugation.workflow.minimization": (
        "LocalMinimizationResult",
        "LocalMinimizationSettings",
        "product_state_pablo_crosslink_requirement",
        "run_post_crosslink_local_minimization",
    ),
    "polyzymd.builders.conjugation.workflow.pablo": (
        "ProductStatePabloDefinitionSummary",
        "ProductStatePabloLibrary",
        "PabloIngestionResult",
        "PabloIngestor",
        "build_product_state_pablo_library",
    ),
    "polyzymd.builders.conjugation.workflow.parameterization": (
        "InterchangeParameterizationResult",
        "InterchangeParameterizationSettings",
        "create_interchange_from_pablo_topology",
    ),
    "polyzymd.builders.conjugation.reaction_roles": (
        "AtomMappedReaction",
        "GenericResolvedReactionPlan",
        "PdbAtomIdentity",
        "resolve_reaction_roles_from_identity_map",
    ),
    "polyzymd.builders.conjugation.pablo_reaction": (
        "PabloProductDefinitionDiagnostic",
        "PabloReactionDiagnostic",
        "PabloReactionRequest",
        "explore_pablo_residue_reaction",
    ),
}

LEGACY_WORKFLOW_EXPORTS = (
    "ConjugatedPolymerSystemSettings",
    "build_conjugated_polymer_system_from_config",
    "build_conjugated_polymer_system_from_config_path",
)


def _import_current_module(module_name: str) -> ModuleType:
    """Import a current compatibility target, skipping only optional heavy deps."""
    try:
        return importlib.import_module(module_name)
    except ImportError as exc:
        missing_root = (exc.name or "").split(".", maxsplit=1)[0]
        if missing_root in OPTIONAL_HEAVY_IMPORT_ROOTS:
            pytest.skip(f"{module_name} currently requires optional dependency {exc.name!r}: {exc}")
        raise


@pytest.mark.parametrize(
    ("module_name", "expected_exports"),
    sorted(CURRENT_HELPER_MODULE_EXPORTS.items()),
)
def test_current_conjugation_helper_modules_import(module_name, expected_exports):
    """Working POC helper module paths should remain import-compatible."""
    module = _import_current_module(module_name)

    assert module.__name__ == module_name
    for name in expected_exports:
        assert hasattr(module, name), f"{module_name} missing {name}"


def test_workflow_preparation_imports_match_legacy_path():
    """The new workflow preparation path should expose the legacy objects."""
    legacy = _import_current_module("polyzymd.builders.conjugation.protein_preparation")
    workflow = _import_current_module("polyzymd.builders.conjugation.workflow.preparation")

    for name in CURRENT_HELPER_MODULE_EXPORTS["polyzymd.builders.conjugation.protein_preparation"]:
        assert getattr(workflow, name) is getattr(legacy, name)


def test_workflow_minimization_imports_match_legacy_path():
    """The new workflow minimization path should expose the legacy objects."""
    legacy = _import_current_module("polyzymd.builders.conjugation.local_minimization")
    workflow = _import_current_module("polyzymd.builders.conjugation.workflow.minimization")

    for name in CURRENT_HELPER_MODULE_EXPORTS["polyzymd.builders.conjugation.local_minimization"]:
        assert getattr(workflow, name) is getattr(legacy, name)


def test_workflow_pablo_imports_match_legacy_paths():
    """The workflow Pablo path should expose product library and ingestion helpers."""
    product_pablo = _import_current_module("polyzymd.builders.conjugation.product_pablo")
    pablo_adapter = _import_current_module("polyzymd.builders.conjugation.pablo_adapter")
    workflow = _import_current_module("polyzymd.builders.conjugation.workflow.pablo")

    for name in CURRENT_HELPER_MODULE_EXPORTS["polyzymd.builders.conjugation.product_pablo"]:
        assert getattr(workflow, name) is getattr(product_pablo, name)
    assert workflow.PabloIngestor is pablo_adapter.PabloIngestor
    assert workflow.PabloIngestionResult is pablo_adapter.PabloIngestionResult


def test_workflow_parameterization_imports_match_legacy_path():
    """The workflow parameterization path should expose Interchange helpers."""
    legacy = _import_current_module("polyzymd.builders.conjugation.parameterization")
    workflow = _import_current_module("polyzymd.builders.conjugation.workflow.parameterization")

    for name in CURRENT_HELPER_MODULE_EXPORTS[
        "polyzymd.builders.conjugation.workflow.parameterization"
    ]:
        assert getattr(workflow, name) is getattr(legacy, name)


def test_legacy_conjugated_system_workflow_exports_import():
    """Config-driven conjugated-polymer workflow entry points should remain available."""
    module = _import_current_module("polyzymd.builders.conjugation.system_workflow")

    for name in LEGACY_WORKFLOW_EXPORTS:
        assert hasattr(module, name), f"system_workflow missing {name}"


def test_current_conjugation_package_facade_reexports_legacy_workflow():
    """The current builders.conjugation facade should keep legacy workflow exports."""
    facade = _import_current_module("polyzymd.builders.conjugation")

    for name in LEGACY_WORKFLOW_EXPORTS:
        assert hasattr(facade, name), f"builders.conjugation facade missing {name}"

    exported_names = set(getattr(facade, "__all__", ()))
    for name in LEGACY_WORKFLOW_EXPORTS:
        assert name in exported_names, f"builders.conjugation __all__ missing {name}"


def test_new_conjugation_package_facade_exports_legacy_workflow_if_present():
    """A future polyzymd.conjugation facade should expose legacy workflow wrappers."""
    module_name = "polyzymd.conjugation"
    if importlib.util.find_spec(module_name) is None:
        pytest.skip("polyzymd.conjugation facade is not present in the current tree")

    facade = _import_current_module(module_name)
    for name in LEGACY_WORKFLOW_EXPORTS:
        assert hasattr(facade, name), f"polyzymd.conjugation facade missing {name}"
