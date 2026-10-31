"""The mutations:, modifications: and modification_parameters: keys of a protein solute."""

from __future__ import annotations

import copy
from pathlib import Path

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
        (
            {"mutations": [{"site": "A:Q5", "to": "4II"}]},
            "noncanonical sites (the mutation to 4II)",
        ),
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
    assert "one parameter set per attachment type in 1.4" in message


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
    assert "exactly one of" in _error(mutations=[{"site": "A:Q5", "sites": ["A:Q5"], "to": "ASN"}])


def test_declared_modifications_expand_every_group_form(tmp_path) -> None:
    library = tmp_path / "glycans"
    library.mkdir()
    (library / "glycam_G42666HT.pdb").write_text("REMARK test\n")
    sites = tmp_path / "sites.csv"
    sites.write_text("﻿site, glycan\nA:S13, G42666HT\n")
    protein = _protein(
        mutations=[{"chains": ["A", "B"], "residues": ["Q5"], "to": "ASN"}],
        modifications=[
            {
                "attach": {
                    **N_GLYCAN,
                    "chains": ["A"],
                    "glycans": {"G42666HT": "N5"},
                    "glycan_library": str(library),
                }
            },
            {
                "attach": {
                    **N_GLYCAN,
                    "atom_name": "OG",
                    "sites_file": str(sites),
                    "glycan_library": str(library),
                }
            },
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


def test_attach_at_must_be_a_position_of_the_chain() -> None:
    chain = {
        "monomers": [{"label": "A", "name": "SBMA", "smiles": "[*]CC([*])C", "probability": 1.0}],
        "end_groups": {
            "head": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
            "tail": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
        },
        "length": 9,
        "random_seed": 1,
    }
    graft = {
        "type": "polymer",
        "site": "A:K8",
        "atom_name": "NZ",
        "fragment": {"polymer": chain},
    }
    message = _error(
        modifications=[{"attach": {**graft, "attach_at": 10}}],
        modification_parameters={"polymer": SAGE},
    )
    assert "attach_at 10 is outside the chain's positions 1..9" in message
    _protein(
        modifications=[{"attach": {**graft, "attach_at": 5}}],
        modification_parameters={"polymer": SAGE},
    )


# --- What `validate` checks before a build ------------------------------------------------

QUICKSTART = Path(__file__).resolve().parents[2] / "examples" / "quickstart"
DATA = Path(__file__).resolve().parents[1] / "data" / "modifying_proteins"


def _buildable(tmp_path: Path, protein: dict, sites: str | None = None, glycan: bool = True):
    """Load the quickstart with ``protein`` added and run the checks before a build."""
    import shutil

    import yaml

    from polyzymd.config.loader import load_config

    folder = tmp_path / "study"
    (folder / "structures" / "glycans").mkdir(parents=True)
    shutil.copy(QUICKSTART / "trpcage.pdb", folder)
    if glycan:
        shutil.copy(DATA / "glycam_G42666HT.pdb", folder / "structures" / "glycans")
    if sites is not None:
        (folder / "sites.csv").write_text(sites)
    config = yaml.safe_load((QUICKSTART / "config.yaml").read_text())
    config["solutes"][0].update(protein)
    (folder / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    load_config(folder / "config.yaml").require_buildable()


def _refused(tmp_path: Path, protein: dict, **kwargs) -> str:
    with pytest.raises(ValueError) as caught:
        _buildable(tmp_path, protein, **kwargs)
    message = str(caught.value)
    assert "\nfix: " in message, message
    return message


GLYCAN_GROUP = {
    "modifications": [{"attach": {**N_GLYCAN, "chains": ["A"], "glycans": {"G42666HT": ["N1"]}}}],
    "modification_parameters": {"glycan": GLYCAM},
}


def test_a_complete_glycan_group_is_buildable(tmp_path) -> None:
    _buildable(tmp_path, GLYCAN_GROUP)


def test_a_missing_glycan_file_is_named(tmp_path) -> None:
    message = _refused(tmp_path, GLYCAN_GROUP, glycan=False)
    assert "no file for glycan G42666HT" in message and "glycan_library" in message


def test_a_missing_sites_file_is_named(tmp_path) -> None:
    protein = {
        "modifications": [{"attach": {**N_GLYCAN, "sites_file": "absent.csv"}}],
        "modification_parameters": {"glycan": GLYCAM},
    }
    assert "sites file" in _refused(tmp_path, protein)


def test_a_glycan_sites_file_needs_a_glycan_column(tmp_path) -> None:
    protein = {
        "modifications": [{"attach": {**N_GLYCAN, "sites_file": "sites.csv"}}],
        "modification_parameters": {"glycan": GLYCAM},
    }
    message = _refused(tmp_path, protein, sites="site\nA:N1\n")
    assert "no 'glycan' column" in message


def test_a_missing_fragment_or_ccd_file_is_named(tmp_path) -> None:
    fragment = {
        "modifications": [
            {"attach": {**N_GLYCAN, "site": "A:N1", "fragment": {"pdb_path": "absent.pdb"}}}
        ],
        "modification_parameters": {"glycan": GLYCAM},
    }
    assert "fragment PDB" in _refused(tmp_path / "a", fragment)
    ccd = {
        "mutations": [{"site": "A:Y3", "to": {"ccd_file": "absent.cif"}}],
        "modification_parameters": {"noncanonical": SAGE},
    }
    assert "CCD file" in _refused(tmp_path / "b", ccd)


def test_an_unknown_site_force_field_or_nagl_model_is_named(tmp_path) -> None:
    acyl = {
        "modifications": [
            {
                "attach": {
                    "type": "generic",
                    "site": "A:K8",
                    "atom_name": "NZ",
                    "fragment": {"smiles": "*C(C)=O", "resname": "ACE"},
                }
            }
        ],
    }
    wrong_ff = {**SAGE, "force_field": ["openff-9.9.9.offxml"]}
    message = _refused(tmp_path / "a", {**acyl, "modification_parameters": {"generic": wrong_ff}})
    assert "openff-9.9.9.offxml" in message
    wrong_model = {**SAGE, "charges": {"method": "nagl", "model": "no-such-model.pt"}}
    message = _refused(
        tmp_path / "b", {**acyl, "modification_parameters": {"generic": wrong_model}}
    )
    assert "no-such-model.pt" in message


def test_two_attachments_at_one_atom_name_both_entries(tmp_path) -> None:
    protein = {
        "modifications": [
            {"attach": {**N_GLYCAN, "chains": ["A"], "glycans": {"G42666HT": ["N1"]}}},
            {
                "attach": {
                    **N_GLYCAN,
                    "site": "A:N1",
                    "fragment": {"pdb_path": "structures/glycans/glycam_G42666HT.pdb"},
                }
            },
        ],
        "modification_parameters": {"glycan": GLYCAM},
    }
    message = _refused(tmp_path, protein)
    assert "modifications.1.attach" in message and "modifications.0.attach" in message


def test_residue_number_selections_without_a_segment_are_flagged() -> None:
    from polyzymd.config.modifications import unscoped_residue_selection

    assert unscoped_residue_selection("resid 141-148")
    assert unscoped_residue_selection("chainid A and resid 5")
    assert not unscoped_residue_selection("segid PROT and resid 141-148")
    assert not unscoped_residue_selection("protein and resid 5")
    assert not unscoped_residue_selection("name CA")


def test_validate_warns_about_a_restraint_by_residue_number_on_a_glycosylated_protein(
    tmp_path,
) -> None:
    import shutil

    import yaml

    from polyzymd.config.loader import load_config
    from polyzymd.config.validation import collect_protocol_warnings

    folder = tmp_path / "study"
    (folder / "structures" / "glycans").mkdir(parents=True)
    shutil.copy(QUICKSTART / "trpcage.pdb", folder)
    shutil.copy(DATA / "glycam_G42666HT.pdb", folder / "structures" / "glycans")
    config = yaml.safe_load((QUICKSTART / "config.yaml").read_text())
    config["solutes"][0].update(GLYCAN_GROUP)

    def restraint(selection: str) -> dict:
        return {
            "type": "flat_bottom",
            "name": "pin",
            "atom1": {"component": "trpcage", "selection": selection},
            "atom2": {"component": "trpcage", "selection": "resid 20 and name CA"},
            "distance": 20.0,
        }

    config["restraints"] = [restraint("resid 2 and name CA")]
    (folder / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    warnings = collect_protocol_warnings(load_config(folder / "config.yaml"))
    assert any("restraint pin atom1" in w and "may also match" in w for w in warnings)
    config["restraints"] = [restraint("protein and resid 2 and name CA")]
    config["restraints"][0]["atom2"]["selection"] = "protein and resid 20 and name CA"
    (folder / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    assert not any(
        "restraint" in w for w in collect_protocol_warnings(load_config(folder / "config.yaml"))
    )


def test_graft_fields_are_part_of_the_modification_identity() -> None:
    """Grafts that differ only in attach_at or in their sequence are different systems."""
    from polyzymd.config.modifications import modification_identity

    chain = {
        "monomers": [
            {"label": "A", "name": "SBMA", "smiles": "[*]CC([*])C", "probability": 1.0},
            {"label": "N", "name": "NHS", "smiles": "[*]CC([*])CO", "probability": 0.0},
        ],
        "end_groups": {
            "head": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
            "tail": [{"label": "H", "smiles": "[*][H]", "probability": 1.0}],
        },
        "length": 9,
        "sequence": {"fixed": {5: "NHS"}},
    }

    def identity(attach_at: int = 5, sequence: dict | None = None) -> str:
        polymer = {**chain, **({"sequence": sequence} if sequence else {})}
        return modification_identity(
            _protein(
                modifications=[
                    {
                        "attach": {
                            "type": "polymer",
                            "site": "A:K8",
                            "atom_name": "NZ",
                            "attach_at": attach_at,
                            "fragment": {"polymer": polymer},
                        }
                    }
                ],
                modification_parameters={"polymer": SAGE},
            )
        )

    base = identity()
    assert base and base == identity()
    assert identity(attach_at=4) != base
    assert identity(sequence={"fixed": {4: "NHS"}}) != base


def test_a_ccd_file_fragment_needs_a_link_atom_and_an_existing_file(tmp_path) -> None:
    entry = {"type": "generic", "site": "A:K8", "atom_name": "NZ"}
    assert "needs `link_atom`" in _error(
        modifications=[{"attach": {**entry, "fragment": {"ccd_file": "x.cif"}}}],
        modification_parameters={"generic": SAGE},
    )
    protein = {
        "modifications": [
            {"attach": {**entry, "fragment": {"ccd_file": "absent.cif", "link_atom": "CB"}}}
        ],
        "modification_parameters": {"generic": SAGE},
    }
    assert "fragment CCD file" in _refused(tmp_path, protein)


def test_a_glycan_attachment_may_leave_its_bonding_atoms_to_the_defaults() -> None:
    entry = {
        "type": "glycan",
        "reaction": "N-glycosylation",
        "chains": ["A"],
        "glycans": {"G1": ["N5"]},
    }
    protein = _protein(
        modifications=[{"attach": entry}], modification_parameters={"glycan": GLYCAM}
    )
    attachment = protein.modifications[0].arguments
    assert attachment.atom_name is None and attachment.fragment_atom_name is None


def test_other_attachments_still_name_the_protein_atom() -> None:
    entry = {"type": "generic", "site": "A:K8", "fragment": {"smiles": "*C(=O)C", "resname": "ACE"}}
    message = _error(modifications=[{"attach": entry}], modification_parameters={"generic": SAGE})
    assert "a generic attachment needs atom_name" in message
