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
