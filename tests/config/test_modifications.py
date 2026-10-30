"""The mutations:, modifications: and modification_parameters: keys of a protein solute."""

from __future__ import annotations

import copy

import pytest
from pydantic import ValidationError

from polyzymd.config.modifications import declared_modifications
from polyzymd.config.schema import ProteinSolute

GLYCAM = {"force_field": ["GLYCAM_06j-1.xml"], "charges": {"method": "library"}}
SAGE = {
    "force_field": ["openff-2.3.0.offxml"],
    "charges": {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"},
}
BASE = {
    "type": "protein",
    "pdb_path": "trpcage.pdb",
    "force_field": ["ff14sb_off_impropers_0.0.4.offxml"],
    "charges": {"method": "library"},
}
N_GLYCAN = {
    "type": "glycan",
    "reaction": "N-glycosylation",
    "atom_name": "ND2",
    "fragment_atom_name": "C1",
    "fragment_resnum": 2,
}


def _protein(**keys) -> ProteinSolute:
    return ProteinSolute.model_validate({**copy.deepcopy(BASE), **keys})


def _error(**keys) -> str:
    with pytest.raises(ValidationError) as caught:
        _protein(**keys)
    return str(caught.value)


def test_an_unmodified_protein_dumps_as_before() -> None:
    dumped = _protein().model_dump(mode="json")
    for key in ("mutations", "modifications", "modification_parameters", "max_charge_shift"):
        assert key not in dumped


def test_an_unmodified_protein_keeps_its_build_hash() -> None:
    from polyzymd.simulation.artifact_integrity import _digest

    explicit = _protein(mutations=[], modifications=[], max_charge_shift=0.005, relax="joint")
    assert _digest(explicit.model_dump(mode="json")) == _digest(_protein().model_dump(mode="json"))


def test_a_glycosylated_protein_validates() -> None:
    protein = _protein(
        mutations=[{"site": "A:Q5", "to": "ASN"}],
        modifications=[{"attach": {**N_GLYCAN, "site": "A:N5", "fragment": {"pdb_path": "g.pdb"}}}],
        modification_parameters={"glycan": GLYCAM},
    )
    assert protein.is_modified
    assert protein.modifications[0].verb == "attach"
    assert protein.mutations[0].is_canonical


def test_a_mutation_inside_modifications_points_to_mutations() -> None:
    message = _error(modifications=[{"mutate": {"site": "A:Q5", "to": "ASN"}}])
    assert "mutations go under `mutations:`" in message


def test_an_unknown_verb_is_named() -> None:
    assert "unknown modification glycosylate" in _error(
        modifications=[{"glycosylate": {"site": "A:N5"}}]
    )


@pytest.mark.parametrize(
    ("keys", "missing"),
    [
        ({"mutations": [{"site": "A:Q5", "to": "4II"}]}, "noncanonical sites (the mutation to 4II)"),
        (
            {
                "modifications": [
                    {
                        "attach": {
                            "type": "generic",
                            "site": "A:K8",
                            "atom_name": "NZ",
                            "fragment": {"smiles": "*C(=O)C", "resname": "ACE"},
                        }
                    }
                ]
            },
            "generic sites (a generic attachment)",
        ),
    ],
)
def test_a_site_type_without_a_parameter_set_is_named(keys, missing) -> None:
    message = _error(**keys)
    assert f"no parameter set for {missing}" in message
    assert "modification_parameters" in message


def test_a_canonical_or_d_mutation_needs_no_parameter_set() -> None:
    protein = _protein(mutations=[{"site": "A:L2", "to": "DLE"}, {"site": "A:Q5", "to": "asn"}])
    assert all(mutation.is_canonical for mutation in protein.mutations)


def test_glycam_is_refused_for_other_site_types() -> None:
    assert "GLYCAM06j parameterizes glycans only; generic" in _error(
        modification_parameters={"generic": GLYCAM}
    )


def test_an_openff_parameter_set_needs_nagl_charges() -> None:
    wrong = {"force_field": ["openff-2.3.0.offxml"], "charges": {"method": "library"}}
    assert "NAGL on a capped model" in _error(modification_parameters={"generic": wrong})


def test_per_entry_parameters_are_refused() -> None:
    entry = {**N_GLYCAN, "site": "A:N5", "fragment": {"pdb_path": "g.pdb"}, "parameters": "x"}
    message = _error(modifications=[{"attach": entry}], modification_parameters={"glycan": GLYCAM})
    assert "set once per site type" in message


def test_bad_site_notation_is_explained() -> None:
    assert "Write chain:residue" in _error(mutations=[{"site": "A5", "to": "ASN"}])


def test_a_group_entry_cannot_take_a_name() -> None:
    entry = {**N_GLYCAN, "chains": ["A"], "glycans": {"G1": ["N5"]}, "name": "G001"}
    assert "cannot set `name`" in _error(
        modifications=[{"attach": entry}], modification_parameters={"glycan": GLYCAM}
    )


def test_glycans_are_for_glycan_attachments_only() -> None:
    entry = {**N_GLYCAN, "type": "generic", "glycans": {"G1": ["A:N5"]}}
    assert "`glycans` is for type: glycan" in _error(
        modifications=[{"attach": entry}], modification_parameters={"generic": SAGE}
    )


def test_an_attachment_needs_a_fragment() -> None:
    entry = {"type": "generic", "site": "A:K8", "atom_name": "NZ"}
    assert "an attachment needs `fragment`" in _error(
        modifications=[{"attach": entry}], modification_parameters={"generic": SAGE}
    )


def test_one_site_form_per_entry() -> None:
    assert "exactly one of" in _error(
        mutations=[{"site": "A:Q5", "sites": ["A:Q5"], "to": "ASN"}]
    )


def test_declared_modifications_expand_every_group_form(tmp_path) -> None:
    library = tmp_path / "glycans"
    library.mkdir()
    (library / "glycam_G42666HT.pdb").write_text("REMARK test\n")
    sites = tmp_path / "sites.csv"
    sites.write_text("﻿site, glycan\nA:S13, G42666HT\n")
    protein = _protein(
        mutations=[{"chains": ["A", "B"], "residues": ["Q5"], "to": "ASN"}],
        modifications=[
            {"attach": {**N_GLYCAN, "chains": ["A"], "glycans": {"G42666HT": "N5"},
                        "glycan_library": str(library)}},
            {"attach": {**N_GLYCAN, "atom_name": "OG", "sites_file": str(sites),
                        "glycan_library": str(library)}},
        ],
        modification_parameters={"glycan": GLYCAM},
    )
    declared = declared_modifications(
        protein.mutations, protein.modifications, protein.modification_parameters
    )
    assert [(d["kind"], d["site"]) for d in declared] == [
        ("mutate", "A:Q5"),
        ("mutate", "B:Q5"),
        ("attach", "A:N5"),
        ("attach", "A:S13"),
    ]
    glycan = declared[2]
    assert glycan["fragment"]["sha256"] and glycan["parameters"]["force_field"] == [
        "GLYCAM_06j-1.xml"
    ]
