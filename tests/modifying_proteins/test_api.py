"""Modifying real proteins through the PolyzyMD API and recording the changes."""

from pathlib import Path

import pytest

from polyzymd import modifying_proteins as mp
from polyzymd.modifying_proteins import AtomKey

ROOT = Path(__file__).resolve().parents[2]
TRPCAGE = ROOT / "examples" / "quickstart" / "trpcage.pdb"
DATA = ROOT / "tests" / "data" / "modifying_proteins"
CHITOBIOSE = DATA / "glycam_G42666HT.pdb"
AZIDOPHENYLALANINE = DATA / "4II.cif"


def test_amide_coupling_records_the_acyl_group_and_the_lost_protons() -> None:
    protein = mp.load_protein(TRPCAGE)
    n_before = protein.mbuild.n_particles
    deprotonated = mp.deprotonate(protein, "A:K8", "NZ")
    acyl = mp.fragment_from_smiles("*C(=O)CCCCCCC", "OC8")
    record = mp.attach(
        protein,
        acyl,
        type="generic",
        site="A:K8",
        atom_name="NZ",
        relax=False,
    )

    assert deprotonated.formal_charge_change == -1
    assert len(deprotonated.removed) == 1
    # Octanoyl, C8H15O, bonds NZ in place of one more proton.
    assert len(record.added) == 24
    assert [key.atom_name for key in record.removed] == ["HZ1"]
    assert record.new_residues == (("A", 21, "", "OC8"),)
    assert record.formal_charge_change == 0
    assert record.name == "M001"
    assert record.arguments["fragment"] == {"smiles": "*C(=O)CCCCCCC", "resname": "OC8"}
    assert protein.mbuild.n_particles == n_before - 2 + 24


def test_glycan_from_a_glycam_file_adds_its_residues() -> None:
    protein = mp.load_protein(TRPCAGE)
    record = mp.attach(
        protein,
        {"pdb_path": str(CHITOBIOSE)},
        type="glycan",
        site="A:N1",
        atom_name="ND2",
        reaction="N-glycosylation",
        fragment_atom_name="C1",
        fragment_resnum=2,
        relax=False,
    )

    assert record.name == "G001"
    assert [residue[3] for residue in record.new_residues] == ["4YB", "0YB"]
    # The ROH residue and the amide proton HD22 leave.
    assert [key.atom_name for key in record.removed] == ["HD22"]
    assert record.arguments["fragment"]["sha256"]
    assert record.residue_after == "ASN"


def test_canonical_and_noncanonical_mutations_are_recorded() -> None:
    protein = mp.load_protein(TRPCAGE)
    alanine = mp.mutate(protein, "A:L2", "ALA", relax=False)
    azido = mp.mutate(protein, "A:Y3", {"ccd_file": str(AZIDOPHENYLALANINE)}, relax=False)
    side_chain = mp.mutate(
        protein, "A:I4", {"smiles": "*CCCCN=[N+]=[N-]", "resname": "AZK"}, relax=False
    )

    assert (alanine.residue_before, alanine.residue_after) == ("LEU", "ALA")
    assert (azido.residue_before, azido.residue_after) == ("TYR", "4II")
    assert azido.arguments["to"]["sha256"]
    assert side_chain.residue_after == "AZK"
    # Backbone atoms keep their identity through a mutation.
    backbone = {AtomKey("A", 3, "", name) for name in ("N", "CA", "C", "O")}
    assert not backbone & set(azido.removed)


def test_click_merged_into_a_mutated_residue_stays_in_the_residue() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.mutate(protein, "A:Y3", {"ccd_file": str(AZIDOPHENYLALANINE)}, relax=False)
    alkyne = mp.fragment_from_smiles("C#Cc1ccccc1", "ALK")
    record = mp.attach(
        protein,
        alkyne,
        type="generic",
        site="A:4II3",
        atom_name="N3",
        reaction="azide-alkyne triazole",
        merge=True,
        relax=False,
    )

    assert record.new_residues == ()
    assert len(record.added) == 14
    assert {key.resnum for key in record.added} == {3}
    assert record.formal_charge_change == 0


def test_mutation_after_an_attachment_is_refused() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.deprotonate(protein, "A:K8", "NZ")
    with pytest.raises(ValueError, match="Mutations come before"):
        mp.mutate(protein, "A:L2", "ALA", relax=False)


def test_attachments_are_named_by_type_in_order() -> None:
    protein = mp.load_protein(TRPCAGE)
    names = [
        mp.attach(
            protein,
            {"smiles": "*C(=O)C", "resname": "ACY"},
            type="generic",
            site=site,
            atom_name="OG",
            relax=False,
        ).name
        for site in ("A:S13", "A:S14")
    ]
    glycan = mp.attach(
        protein,
        {"pdb_path": str(CHITOBIOSE)},
        type="glycan",
        site="A:S20",
        atom_name="OG",
        reaction="O-glycosylation",
        fragment_atom_name="C1",
        fragment_resnum=2,
        relax=False,
    )
    assert names == ["M001", "M002"]
    assert glycan.name == "G001"
    with pytest.raises(ValueError, match="already exists"):
        mp.attach(
            protein,
            {"smiles": "*C", "resname": "MET"},
            type="generic",
            site="A:K8",
            atom_name="NZ",
            name="M001",
            relax=False,
        )


def test_plan_shows_the_change_before_it_is_applied() -> None:
    protein = mp.load_protein(TRPCAGE)
    plan = mp.plan_deprotonate(protein, "A:K8", "NZ")
    assert plan.kind == "deprotonate"
    assert len(plan.leaving) == 1 and plan.leaving[0].startswith("LYS8:HZ")
    assert "NZ" in str(plan)
    assert protein.records == []
    record = mp.apply(protein, plan)
    assert record.kind == "deprotonate"


def test_unknown_attachment_type_is_refused() -> None:
    protein = mp.load_protein(TRPCAGE)
    with pytest.raises(ValueError, match="glycan, polymer, generic"):
        mp.plan_attach(
            protein, {"smiles": "*C", "resname": "MET"}, type="dye", site="A:K8", atom_name="NZ"
        )


@pytest.mark.parametrize(
    ("verb", "site", "atom_name"),
    [("protonate", "A:S13", "OG"), ("protonate", "A:D9", "OD1"), ("protonate", "A:K8", "NZ")],
)
def test_protonation_that_changes_nothing_is_refused(verb, site, atom_name) -> None:
    protein = mp.load_protein(TRPCAGE)
    with pytest.raises(ValueError, match="would change nothing"):
        getattr(mp, verb)(protein, site, atom_name)
    assert protein.records == []


def test_repeated_deprotonation_is_refused() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.deprotonate(protein, "A:K8", "NZ")
    with pytest.raises(ValueError, match="would change nothing"):
        mp.deprotonate(protein, "A:K8", "NZ")
    assert len(protein.records) == 1


def test_a_planned_mutation_is_refused_after_another_change() -> None:
    protein = mp.load_protein(TRPCAGE)
    plan = mp.plan_mutate(protein, "A:L2", "ALA", relax=False)
    mp.deprotonate(protein, "A:K8", "NZ")
    with pytest.raises(ValueError, match="Mutations come before"):
        mp.apply(protein, plan)


def test_a_plan_cannot_be_applied_twice() -> None:
    protein = mp.load_protein(TRPCAGE)
    plan = mp.plan_mutate(protein, "A:L2", "ALA", relax=False)
    mp.apply(protein, plan)
    with pytest.raises(ValueError, match="already applied"):
        mp.apply(protein, plan)


def test_automatic_names_skip_names_already_taken() -> None:
    protein = mp.load_protein(TRPCAGE)
    acetyl = {"smiles": "*C(=O)C", "resname": "ACY"}
    names = [
        mp.attach(
            protein, acetyl, type="generic", site=site, atom_name="OG", relax=False, **name
        ).name
        for site, name in (("A:S13", {"name": "M002"}), ("A:S14", {}), ("A:S20", {}))
    ]
    assert names == ["M002", "M001", "M003"]


def test_attachment_names_must_fit_a_segment_id() -> None:
    protein = mp.load_protein(TRPCAGE)
    with pytest.raises(ValueError, match="one to four letters or digits"):
        mp.attach(
            protein,
            {"smiles": "*C(=O)C", "resname": "ACY"},
            type="generic",
            site="A:S13",
            atom_name="OG",
            name="ACETYL",
            relax=False,
        )
    assert protein.records == []


def test_site_code_ending_in_a_digit_is_read_against_the_structure() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.mutate(protein, "A:Y3", {"smiles": "*CCO", "resname": "HY3"}, relax=False)
    # A:HY33 reads as HY 33 or HY3 3; only residue 3 (HY3) exists.
    plan = mp.plan_attach(
        protein,
        {"smiles": "*C(=O)C", "resname": "ACY"},
        type="generic",
        site="A:HY33",
        atom_name="O1",
    )
    assert (plan.site.code, plan.site.resnum) == ("HY3", 3)


def test_formal_charge_change_comes_from_the_chemistry() -> None:
    protein = mp.load_protein(TRPCAGE)
    # A C-terminal N-methylamide removes the carboxylate's charged OXT.
    amide = mp.attach(
        protein,
        {"smiles": "*NC", "resname": "NME"},
        type="generic",
        site="A:S20",
        atom_name="C",
        leaving_atom_names=["OXT"],
        relax=False,
    )
    assert [key.atom_name for key in amide.removed] == ["OXT"]
    assert amide.formal_charge_change == +1


def test_formal_charge_change_of_a_lysine_acylation() -> None:
    protein = mp.load_protein(TRPCAGE)
    deprotonated = mp.deprotonate(protein, "A:K8", "NZ")
    acylated = mp.attach(
        protein,
        {"smiles": "*C(=O)C", "resname": "ACY"},
        type="generic",
        site="A:K8",
        atom_name="NZ",
        relax=False,
    )
    assert (deprotonated.formal_charge_change, acylated.formal_charge_change) == (-1, 0)
