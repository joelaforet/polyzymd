"""Tests for Phase 3 declarative conjugation mechanism planning."""

from __future__ import annotations

import json

import pytest
from pydantic import ValidationError

from polyzymd.builders.conjugation.builder import CovalentModificationBuilder
from polyzymd.builders.conjugation.diagnostics import DiagnosticCode
from polyzymd.builders.conjugation.exceptions import ConjugationNotImplementedError
from polyzymd.builders.conjugation.mechanism_library import (
    get_builtin_mechanism,
    list_builtin_mechanisms,
)
from polyzymd.builders.conjugation.mechanisms import GraphEditPlan, ReactionMechanism
from polyzymd.builders.conjugation.moieties import normalize_moiety_descriptor
from polyzymd.builders.conjugation.sites import normalize_attachment_site
from polyzymd.config.schema import (
    ConjugationConfig,
    ConjugationMoietyConfig,
    ConjugationSiteConfig,
)


def _valid_mechanism_definition() -> dict:
    """Build a minimal valid mechanism definition for validation tests.

    Returns
    -------
    dict
        Mechanism definition suitable for :class:`ReactionMechanism` validation.
    """
    return {
        "identifier": "test_mechanism",
        "display_name": "Test mechanism",
        "allowed_sites": [{"residue_name": "LYS", "atom_name": "NZ"}],
        "moiety_reactive_group": {"group": "test_group", "anchor_atom": "C1"},
        "graph_edits": {"add_bonds": [{"site_atom": "NZ", "moiety_atom": "C1"}]},
        "charge_patch_hint": {"strategy": "defer"},
        "rationale": "Used by tests to validate declarative mechanism models",
    }


def _plausible_nhs_lys_config(tmp_path) -> ConjugationConfig:
    """Build a plausible construct-mode NHS-Lys attachment config."""
    return ConjugationConfig(
        enabled=True,
        mode="construct",
        attachments=[
            {
                "name": "lys23-peg",
                "site": {
                    "chain_id": "A",
                    "residue_name": "LYS",
                    "residue_number": 23,
                    "atom_name": "NZ",
                },
                "moiety": {
                    "name": "PEG-NHS",
                    "role": "polymer",
                    "input_path": tmp_path / "peg_nhs.sdf",
                    "smiles": "CC(=O)ON1C(=O)CCC1=O",
                },
                "mechanism": {"name": "nhs_lys_amide"},
            }
        ],
    )


def test_builtin_mechanism_library_loads_all_required_mechanisms():
    """Built-in mechanism lookup should expose all Phase 3 declarations."""
    names = set(list_builtin_mechanisms())

    assert {"nhs_lys_amide", "n_glycosidic_asn", "residue_replacement"} <= names
    assert get_builtin_mechanism("NHS_LYS_AMIDE").identifier == "nhs_lys_amide"


def test_duplicate_site_rules_are_rejected():
    """Duplicate residue and atom rules should fail mechanism validation."""
    definition = _valid_mechanism_definition()
    definition["allowed_sites"] = [
        {"residue_name": "LYS", "atom_name": "NZ"},
        {"residue_name": "lys", "atom_name": "nz"},
    ]

    with pytest.raises(ValidationError, match="Duplicate allowed site rules"):
        ReactionMechanism.model_validate(definition)


def test_duplicate_graph_edit_bonds_are_rejected():
    """Duplicate bond placeholders should fail graph edit validation."""
    with pytest.raises(ValidationError, match="Duplicate bond specifications"):
        GraphEditPlan.model_validate(
            {
                "add_bonds": [
                    {"site_atom": "NZ", "moiety_atom": "C1"},
                    {"site_atom": "nz", "moiety_atom": "c1"},
                ]
            }
        )


def test_missing_required_mechanism_fields_are_rejected():
    """Mechanism declarations should require rationale and core chemistry fields."""
    definition = _valid_mechanism_definition()
    definition.pop("moiety_reactive_group")

    with pytest.raises(ValidationError, match="moiety_reactive_group"):
        ReactionMechanism.model_validate(definition)


def test_explicit_atom_site_normalization_works():
    """Explicit site configs should normalize chain, residue, and atom labels."""
    site = normalize_attachment_site(
        ConjugationSiteConfig(
            chain_id="a",
            residue_name="lys",
            residue_number=23,
            atom_name="nz",
        )
    )

    assert site.chain_id == "A"
    assert site.residue_name == "LYS"
    assert site.residue_number == 23
    assert site.atom_name == "NZ"
    assert site.selector_mode == "explicit_atom"


def test_missing_site_fields_produce_clear_errors():
    """Missing explicit site atom or residue fields should fail clearly."""
    with pytest.raises(ValueError, match="residue_name, residue_number, and atom_name"):
        normalize_attachment_site(ConjugationSiteConfig(chain_id="A", residue_name="LYS"))


def test_moiety_descriptor_normalizes_file_metadata(tmp_path):
    """Moiety descriptors should preserve names and file metadata."""
    input_path = tmp_path / "peg_nhs.sdf"
    moiety = normalize_moiety_descriptor(
        ConjugationMoietyConfig(
            name="PEG-NHS",
            role="Polymer",
            input_path=input_path,
            smiles="CCO",
            residue_name="peg",
        )
    )

    assert moiety.name == "PEG-NHS"
    assert moiety.role == "polymer"
    assert moiety.kind == "file"
    assert moiety.input_path == input_path
    assert moiety.smiles == "CCO"
    assert moiety.residue_name == "PEG"
    assert moiety.metadata["input_path"] == str(input_path)
    assert moiety.metadata["smiles"] == "CCO"


def test_moiety_descriptor_classifies_smiles_metadata():
    """SMILES-only moieties should remain declarative but identifiable."""
    moiety = normalize_moiety_descriptor(
        ConjugationMoietyConfig(name="NHS-linker", role="moiety", smiles="CC(=O)O")
    )

    assert moiety.kind == "smiles"
    assert moiety.metadata["has_smiles"] is True
    assert moiety.metadata["has_input_path"] is False


def test_construct_mode_validates_nhs_lys_and_writes_diagnostics(tmp_path):
    """Construct mode should plan plausible NHS-Lys attachments before failing."""
    builder = CovalentModificationBuilder(_plausible_nhs_lys_config(tmp_path), output_dir=tmp_path)

    with pytest.raises(ConjugationNotImplementedError, match="graph surgery"):
        builder.build(object())

    diagnostics = json.loads((tmp_path / "conjugation_diagnostics.json").read_text())
    codes = [event["code"] for event in diagnostics["diagnostics"]]
    assert DiagnosticCode.MECHANISM_VALIDATION.value in codes
    assert DiagnosticCode.SITE_SELECTION.value in codes
    assert DiagnosticCode.MOIETY_NORMALIZATION.value in codes
    assert DiagnosticCode.UNSUPPORTED_OPERATION.value in codes

    metadata = json.loads((tmp_path / "conjugation_metadata.json").read_text())
    planned = metadata["attachments"][0]
    assert planned["mechanism"]["identifier"] == "nhs_lys_amide"
    assert planned["site"]["residue_name"] == "LYS"
    assert planned["moiety"]["role"] == "polymer"


def test_unknown_construct_mechanism_fails_clearly(tmp_path):
    """Unknown mechanism identifiers should produce diagnostics before failure."""
    config = _plausible_nhs_lys_config(tmp_path)
    config.attachments[0].mechanism.name = "unknown_mechanism"
    builder = CovalentModificationBuilder(config, output_dir=tmp_path)

    with pytest.raises(ConjugationNotImplementedError, match="Unknown conjugation mechanism"):
        builder.build(object())

    diagnostics = json.loads((tmp_path / "conjugation_diagnostics.json").read_text())
    mechanism_errors = [
        event
        for event in diagnostics["diagnostics"]
        if event["code"] == DiagnosticCode.MECHANISM_VALIDATION.value
        and event["severity"] == "error"
    ]
    assert mechanism_errors
    assert mechanism_errors[0]["details"]["mechanism"] == "unknown_mechanism"
