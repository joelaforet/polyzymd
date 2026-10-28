"""The atoms each modification owns, and the segment of each atom."""

from pathlib import Path

from polyzymd import modifying_proteins as mp

ROOT = Path(__file__).resolve().parents[2]
TRPCAGE = ROOT / "examples" / "quickstart" / "trpcage.pdb"
DATA = ROOT / "tests" / "data" / "modifying_proteins"
CHITOBIOSE = DATA / "glycam_G42666HT.pdb"
AZIDOPHENYLALANINE = DATA / "4II.cif"

LYS_SIDE_CHAIN = {"CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3", "CE", "HE2", "HE3"}


def _glycosylate(protein, site, atom_name, reaction):
    return mp.attach(
        protein,
        {"pdb_path": str(CHITOBIOSE)},
        type="glycan",
        site=site,
        atom_name=atom_name,
        reaction=reaction,
        fragment_atom_name="C1",
        fragment_resnum=2,
        relax=False,
    )


def test_acylated_lysine_site_is_its_side_chain_and_the_acyl_group() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.deprotonate(protein, "A:K8", "NZ")
    record = mp.attach(
        protein,
        {"smiles": "*C(=O)CCCCCCC", "resname": "OC8"},
        type="generic",
        site="A:K8",
        atom_name="NZ",
        relax=False,
    )
    (site,) = mp.modification_sites(protein)

    lysine = {atom.atom_name for atom in site.atoms if atom.resnum == 8}
    assert lysine == LYS_SIDE_CHAIN | {"NZ", "HZ2"}
    assert len(site.atoms) == 14 + 24
    assert set(record.added) <= set(site.atoms)
    assert (site.name, site.type, site.input_residues) == ("M001", "generic", ("LYS",))
    keys = protein.atom_keys()
    assert [keys[index] for index in site.atom_indices] == list(site.atoms)

    segments = mp.segments(protein)
    assert segments.count("M001") == 24
    assert {segments[i] for i, key in enumerate(keys) if key.resnum == 8} == {"PROT"}


def test_glycan_site_spans_the_asparagine_side_chain_and_the_glycan() -> None:
    protein = mp.load_protein(TRPCAGE)
    record = _glycosylate(protein, "A:N1", "ND2", "N-glycosylation")
    (site,) = mp.modification_sites(protein)

    asparagine = {atom.atom_name for atom in site.atoms if atom.resnum == 1}
    assert asparagine == {"CB", "HB2", "HB3", "CG", "OD1", "ND2", "HD21"}
    assert len(site.atoms) == 7 + len(record.added)
    assert site.type == "glycan"
    assert mp.segments(protein).count("G001") == len(record.added)


def test_noncanonical_residue_is_a_site_and_a_click_on_it_grows_it() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.mutate(protein, "A:L2", "ALA", relax=False)
    mp.mutate(protein, "A:Y3", {"ccd_file": str(AZIDOPHENYLALANINE)}, relax=False)
    (site,) = mp.modification_sites(protein)
    assert (site.name, site.type, site.input_residues) == ("R001", "noncanonical", ("TYR",))
    n_side_chain = len(site.atoms)

    click = mp.attach(
        protein,
        {"smiles": "C#Cc1ccccc1", "resname": "ALK"},
        type="generic",
        site="A:4II3",
        atom_name="N3",
        reaction="azide-alkyne triazole",
        merge=True,
        relax=False,
    )
    (grown,) = mp.modification_sites(protein)
    assert grown.name == "R001"
    assert len(grown.atoms) == n_side_chain + len(click.added)
    assert set(mp.segments(protein)) == {"PROT"}


def test_protonation_changes_make_sites_only_outside_the_ff14sb_variants() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.deprotonate(protein, "A:K8", "NZ")  # LYN: an ff14SB residue
    mp.deprotonate(protein, "A:S13", "OG")  # no ff14SB variant
    mp.protonate(protein, "A:S14", "OG")  # claimed by the attachment below
    mp.attach(
        protein,
        {"smiles": "*C(=O)C", "resname": "ACY"},
        type="generic",
        site="A:S14",
        atom_name="OG",
        relax=False,
    )
    sites = {site.name: site for site in mp.modification_sites(protein)}
    assert set(sites) == {"M001", "R001"}
    assert sites["R001"].residues == (("A", 13, ""),)
    assert sites["M001"].residues == (("A", 14, ""),)


def test_backbone_anchor_brings_its_atom_and_hydrogens_into_the_site() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.attach(
        protein,
        {"smiles": "*C(=O)C", "resname": "ACY"},
        type="generic",
        site="A:N1",
        atom_name="N",
        relax=False,
    )
    (site,) = mp.modification_sites(protein)
    names = {atom.atom_name for atom in site.atoms if atom.resnum == 1}
    assert "N" in names and "CA" not in names and "C" not in names
    assert names & {"H", "H2", "H3"}


def test_two_glycans_are_two_sites_and_two_segments() -> None:
    protein = mp.load_protein(TRPCAGE)
    _glycosylate(protein, "A:S13", "OG", "O-glycosylation")
    _glycosylate(protein, "A:S20", "OG", "O-glycosylation")
    assert [site.name for site in mp.modification_sites(protein)] == ["G001", "G002"]
    assert {"G001", "G002", "PROT"} == set(mp.segments(protein))
