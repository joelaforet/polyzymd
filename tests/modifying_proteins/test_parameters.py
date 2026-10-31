"""Parameters and charges of modified proteins: protein backbone, site side chain."""

import shutil
from pathlib import Path

import numpy as np
import pytest

from polyzymd import modifying_proteins as mp
from polyzymd.modifying_proteins.parameters import (
    ParameterSet,
    parameterize,
    reference_protein,
    to_openff,
)

ROOT = Path(__file__).resolve().parents[2]
TRPCAGE = ROOT / "examples" / "quickstart" / "trpcage.pdb"
DATA = ROOT / "tests" / "data" / "modifying_proteins"
CHITOBIOSE = DATA / "glycam_G42666HT.pdb"
AZIDOPHENYLALANINE = DATA / "4II.cif"

SAGE = ParameterSet(
    ("openff-2.3.0.offxml",), {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"}
)
#: ff14SB library charges of the lysine backbone (amber14 protein.ff14SB, LYS).
LYS_BACKBONE = {"N": -0.3479, "H": 0.2747, "CA": -0.24, "HA": 0.1426, "C": 0.7341, "O": -0.5894}


def _acylated_lysine():
    protein = mp.load_protein(TRPCAGE)
    mp.deprotonate(protein, "A:K8", "NZ")
    mp.attach(
        protein,
        {"smiles": "*C(=O)CCCCCCC", "resname": "OC8"},
        type="generic",
        site="A:K8",
        atom_name="NZ",
        relax=False,
    )
    return protein


@pytest.fixture(scope="module")
def acylated():
    protein = _acylated_lysine()
    interchange, record = parameterize(protein, {"generic": SAGE}, max_charge_shift=0.005)
    return protein, interchange, record


def _charges(interchange) -> np.ndarray:
    from openff.units import unit

    charges = interchange["Electrostatics"].charges
    values = np.zeros(len(charges))
    for key, charge in charges.items():
        values[key.atom_indices[0]] = charge.m_as(unit.elementary_charge)
    return values


def _index(protein, resnum: int, atom_name: str) -> int:
    keys = protein.atom_keys()
    return next(
        i for i, key in enumerate(keys) if (key.resnum, key.atom_name) == (resnum, atom_name)
    )


def _source(interchange, collection: str, atoms: tuple[int, ...]) -> set[str]:
    """Return 'protein' and/or 'site' for the terms on exactly these atoms."""
    wanted = {atoms, atoms[::-1]}
    return {
        "protein" if potential.id.startswith("protein::") else "site"
        for key, potential in interchange[collection].key_map.items()
        if key.atom_indices in wanted
    }


def test_site_charges_reach_the_target_and_the_protein_stays_integral(acylated) -> None:
    protein, interchange, record = acylated
    (site,) = record.sites
    charges = _charges(interchange)

    assert site.name == "M001" and len(site.atom_indices) == 38
    assert charges[site.atom_indices].sum() == pytest.approx(site.target, abs=1e-6)
    # The target leaves the site and its residue's ff14SB backbone at formal charge 0.
    assert site.target == pytest.approx(-sum(LYS_BACKBONE.values()), abs=1e-4)
    assert abs(site.shift_per_atom) < 0.005
    assert site.shift_total == pytest.approx(site.target - site.nagl_sum)
    assert charges.sum() == pytest.approx(0.0, abs=1e-6)
    assert record.total_charge == pytest.approx(0.0, abs=1e-6)


def test_modified_residue_backbone_keeps_the_input_residue_charges(acylated) -> None:
    protein, interchange, _ = acylated
    charges = _charges(interchange)
    for name, value in LYS_BACKBONE.items():
        assert charges[_index(protein, 8, name)] == pytest.approx(value, abs=1e-4)


def test_backbone_terms_are_the_proteins_and_side_chain_terms_the_sites(acylated) -> None:
    protein, interchange, _ = acylated
    n, ca, c = (_index(protein, 8, name) for name in ("N", "CA", "C"))
    cb, cg = _index(protein, 8, "CB"), _index(protein, 8, "CG")
    previous_c = _index(protein, 7, "C")
    next_n = _index(protein, 9, "N")

    assert _source(interchange, "Bonds", (previous_c, n)) == {"protein"}
    assert _source(interchange, "Bonds", (c, next_n)) == {"protein"}
    assert _source(interchange, "Bonds", (n, ca)) == {"protein"}
    assert _source(interchange, "ProperTorsions", (previous_c, n, ca, c)) == {"protein"}  # phi
    assert _source(interchange, "ProperTorsions", (n, ca, c, next_n)) == {"protein"}  # psi
    assert _source(interchange, "Bonds", (ca, cb)) == {"site"}
    assert _source(interchange, "Angles", (n, ca, cb)) == {"site"}
    assert _source(interchange, "ProperTorsions", (n, ca, cb, cg)) == {"site"}  # chi1


def test_record_names_segments_and_force_fields(acylated) -> None:
    protein, _, record = acylated
    data = record.to_dict()
    assert data["protein_force_field"] == ["ff14sb_off_impropers_0.0.4.offxml"]
    assert data["sites"][0]["force_field"] == ["openff-2.3.0.offxml"]
    assert data["sites"][0]["backbone_residues"] == ["LYS"]
    assert data["segments"].count("M001") == 24
    assert len(data["segments"]) == len(protein.atom_keys())


@pytest.mark.skipif(shutil.which("gmx") is None, reason="GROMACS is not installed")
def test_openmm_and_gromacs_agree_on_the_energy(acylated, tmp_path, monkeypatch) -> None:
    from openff.interchange.drivers import get_gromacs_energies, get_openmm_energies
    from openff.units import unit

    protein, interchange, _ = acylated
    interchange.positions = to_openff(protein).conformers[0]
    interchange.box = np.eye(3) * 6.0 * unit.nanometer
    monkeypatch.chdir(tmp_path)
    gromacs = get_gromacs_energies(interchange)
    openmm = get_openmm_energies(interchange, combine_nonbonded_forces=True)
    for term in ("Bond", "Angle", "Torsion"):
        a = openmm[term].m_as(unit.kilojoule_per_mole)
        b = gromacs[term].m_as(unit.kilojoule_per_mole)
        assert a == pytest.approx(b, rel=1e-4, abs=1e-2), term
    # PME and cutoff details differ slightly between the engines' nonbonded sums.
    nonbonded = (gromacs["vdW"] + gromacs["Electrostatics"]).m_as(unit.kilojoule_per_mole)
    assert openmm["Nonbonded"].m_as(unit.kilojoule_per_mole) == pytest.approx(nonbonded, rel=1e-3)


def test_noncanonical_residue_is_a_site_and_a_canonical_mutation_is_not() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.mutate(protein, "A:L2", "ALA", relax=False)
    mp.mutate(protein, "A:Y3", {"ccd_file": str(AZIDOPHENYLALANINE)}, relax=False)
    interchange, record = parameterize(protein, {"noncanonical": SAGE}, max_charge_shift=0.005)

    (site,) = record.sites
    assert (site.type, site.backbone_residues) == ("noncanonical", ["TYR"])
    assert all(protein.atom_keys()[i].resnum == 3 for i in site.atom_indices)
    charges = _charges(interchange)
    # Ala2 is a canonical ff14SB residue: library charge on its methyl carbon.
    assert charges[_index(protein, 2, "CB")] == pytest.approx(-0.1825, abs=1e-4)
    assert charges.sum() == pytest.approx(round(charges.sum()), abs=1e-6)


def test_reference_replays_only_changes_that_keep_ff14sb_residues() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.mutate(protein, "A:L2", "ALA", relax=False)
    mp.deprotonate(protein, "A:K8", "NZ")
    mp.attach(
        protein,
        {"smiles": "*C(=O)CCCCCCC", "resname": "OC8"},
        type="generic",
        site="A:K8",
        atom_name="NZ",
        relax=False,
    )
    reference = reference_protein(protein)
    names = {(key.resnum, key.atom_name) for key in reference.atom_keys()}
    assert (2, "CG") not in names  # Leu2 became Ala
    assert (8, "HZ3") in names  # Lys8 keeps its charged side chain, not LYN
    assert not any(key.resnum > 20 for key in reference.atom_keys())  # no acyl residue


def test_multi_residue_glycan_on_an_openff_force_field() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.attach(
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
    interchange, record = parameterize(protein, {"glycan": SAGE}, max_charge_shift=0.005)
    (site,) = record.sites
    residues = {protein.atom_keys()[i].resnum for i in site.atom_indices}
    assert len(residues) == 3  # Asn1 side chain and both GlcNAc residues
    charges = _charges(interchange)
    assert charges[site.atom_indices].sum() == pytest.approx(site.target, abs=1e-6)
    assert site.capped_model_atoms > len(site.atom_indices)


def test_a_shift_beyond_the_bound_is_refused() -> None:
    with pytest.raises(ValueError, match="charge shift"):
        parameterize(_acylated_lysine(), {"generic": SAGE}, max_charge_shift=1e-7)


def test_a_site_type_without_a_parameter_set_is_refused() -> None:
    with pytest.raises(ValueError, match="No parameter set for the site type"):
        parameterize(_acylated_lysine(), {"glycan": SAGE}, max_charge_shift=0.005)


def test_openff_parameter_sets_need_nagl_charges() -> None:
    with pytest.raises(ValueError, match="NAGL"):
        ParameterSet(("openff-2.3.0.offxml",), {"method": "am1bcc"})


def test_a_mutated_then_modified_residue_keeps_the_mutants_backbone_charges() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.mutate(protein, "A:S13", "CYS", relax=False)
    mp.attach(
        protein,
        {"smiles": "CCN1C(=O)C=CC1=O", "resname": "NEM"},
        type="generic",
        site="A:C13",
        atom_name="SG",
        reaction="thiol-maleimide",
        relax=False,
    )
    interchange, record = parameterize(protein, {"generic": SAGE}, max_charge_shift=0.005)
    (site,) = record.sites
    assert site.backbone_residues == ["CYS"]
    charges = _charges(interchange)
    # amber14 protein.ff14SB CYS backbone.
    cys = {"N": -0.4157, "H": 0.2719, "CA": 0.0213, "HA": 0.1124, "C": 0.5973, "O": -0.5679}
    for name, value in cys.items():
        assert charges[_index(protein, 13, name)] == pytest.approx(value, abs=1e-4)
    assert charges.sum() == pytest.approx(1.0, abs=1e-6)


def test_a_protein_without_sites_is_plain_ff14sb() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.mutate(protein, "A:L2", "ALA", relax=False)
    interchange, record = parameterize(protein, {}, max_charge_shift=0.005)
    assert record.sites == []
    charges = _charges(interchange)
    assert charges[_index(protein, 2, "CB")] == pytest.approx(-0.1825, abs=1e-4)
    assert all(not p.id.startswith("protein::") for p in interchange["Bonds"].key_map.values())


#: Plumbing tests on terminal sites widen the shift bound; the bound itself is
#: tested per case below.
WIDE = 0.05


def _neutral_n_terminus(protein):
    mp.deprotonate(protein, "A:N1", "N")


def _protonated_c_terminus(protein):
    # Ser20 follows Pro19: the capped model crosses the proline ring.
    mp.protonate(protein, "A:S20", "OXT")


def _n_terminal_acetyl(protein):
    mp.deprotonate(protein, "A:N1", "N")
    mp.attach(
        protein,
        {"smiles": "*C(C)=O", "resname": "ACE"},
        type="generic",
        site="A:N1",
        atom_name="N",
        relax=False,
    )


def _c_terminal_methylamide(protein):
    mp.attach(
        protein,
        {"smiles": "*NC", "resname": "NMA"},
        type="generic",
        site="A:S20",
        atom_name="C",
        leaving_atom_names=["OXT"],
        relax=False,
    )


TERMINAL_CHANGES = {
    "N-terminal acetyl": (_n_terminal_acetyl, "generic", 1, 0.0),
    "C-terminal methylamide": (_c_terminal_methylamide, "generic", 20, 2.0),
}

ASN_SIDE_CHAIN = {"CB", "HB2", "HB3", "CG", "OD1", "ND2", "HD21", "HD22"}
SER_SIDE_CHAIN = {"CB", "HB2", "HB3", "OG", "HG"}
#: amber14 protein.ff14SB main-chain (not terminal) CA charges.
MAIN_CHAIN_CA = {1: 0.0143, 20: -0.0249}  # ASN, SER


def _terminal(change, max_charge_shift):
    apply_change, site_type, resnum, total = TERMINAL_CHANGES[change]
    protein = mp.load_protein(TRPCAGE)
    apply_change(protein)
    interchange, record = parameterize(
        protein, {site_type: SAGE}, max_charge_shift=max_charge_shift
    )
    return protein, interchange, record, resnum, total


@pytest.mark.parametrize("change", sorted(TERMINAL_CHANGES))
def test_a_terminal_site_leaves_the_rest_of_the_residue_on_the_main_chain_template(change) -> None:
    from polyzymd.modifying_proteins.parameters import library_charges, reference_protein

    protein, interchange, record, resnum, total = _terminal(change, WIDE)
    (site,) = record.sites
    charges = _charges(interchange)
    reference = library_charges(to_openff(reference_protein(protein)))
    keys = protein.atom_keys()
    side_chain = ASN_SIDE_CHAIN if resnum == 1 else SER_SIDE_CHAIN
    for name in side_chain:
        index = _index(protein, resnum, name)
        assert index not in site.atom_indices
        assert charges[index] == pytest.approx(reference[keys[index]], abs=1e-10)
    ca, cb = _index(protein, resnum, "CA"), _index(protein, resnum, "CB")
    if ca not in site.atom_indices:
        # The terminus is modified, so the untouched CA takes the main-chain template.
        assert charges[ca] == pytest.approx(MAIN_CHAIN_CA[resnum], abs=1e-4)
    assert charges[site.atom_indices].sum() == pytest.approx(site.target, abs=1e-6)
    assert charges.sum() == pytest.approx(total, abs=1e-6)
    if resnum == 1:
        n = _index(protein, 1, "N")
        assert n in site.atom_indices and ca not in site.atom_indices
        assert _source(interchange, "Bonds", (n, ca)) == {"site"}
        assert _source(interchange, "Bonds", (ca, cb)) == {"protein"}
    else:
        c = _index(protein, resnum, "C")
        assert c in site.atom_indices
        assert _source(interchange, "Bonds", (ca, c)) == {"site"}


@pytest.mark.parametrize("change", ["C-terminal methylamide", "N-terminal acetyl"])
def test_a_terminal_label_stays_within_the_shift_bound(change) -> None:
    _, _, record, _, _ = _terminal(change, 0.005)
    assert abs(record.sites[0].shift_per_atom) < 0.002


@pytest.mark.parametrize(
    "change", [_neutral_n_terminus, _protonated_c_terminus], ids=["neutral N", "protonated C"]
)
def test_a_terminus_changed_by_protonation_alone_is_refused(change) -> None:
    protein = mp.load_protein(TRPCAGE)
    change(protein)
    with pytest.raises(ValueError, match="does not model a neutral N-terminus"):
        parameterize(protein, {"noncanonical": SAGE})


def test_only_a_touched_terminus_is_capped_in_the_reference() -> None:
    from polyzymd.modifying_proteins.parameters import modified_termini, reference_protein

    protein = mp.load_protein(TRPCAGE)
    _n_terminal_acetyl(protein)
    sites = mp.modification_sites(protein)
    assert modified_termini(mp.load_protein(TRPCAGE), sites) == {("A", 1, ""): "N"}

    reference = reference_protein(protein, sites)
    caps = [
        (r.site.resnum, r.arguments["atom_name"]) for r in reference.records if r.kind == "attach"
    ]
    assert caps == [(1, "N")]
    # Ser20's C-terminus is untouched: it keeps OXT and its terminal template.
    assert any(key.resnum == 20 and key.atom_name == "OXT" for key in reference.atom_keys())


def _constrained(interchange) -> dict[tuple[int, int], float]:
    from openff.units import unit

    constraints = interchange["Constraints"]
    return {
        tuple(sorted(key.atom_indices)): constraints.potentials[potential]
        .parameters["distance"]
        .m_as(unit.angstrom)
        for key, potential in constraints.key_map.items()
    }


@pytest.mark.parametrize("modified", [False, True])
def test_every_bond_to_hydrogen_is_constrained_at_its_own_bond_length(modified) -> None:
    from openff.units import unit

    protein = _acylated_lysine() if modified else mp.load_protein(TRPCAGE)
    interchange, _ = parameterize(protein, {"generic": SAGE})
    molecule = to_openff(protein)
    to_hydrogen = {
        tuple(sorted((bond.atom1_index, bond.atom2_index)))
        for bond in molecule.bonds
        if 1 in (bond.atom1.atomic_number, bond.atom2.atomic_number)
    }
    constrained = _constrained(interchange)
    assert set(constrained) == to_hydrogen
    bonds = interchange["Bonds"]
    for key, potential in bonds.key_map.items():
        atoms = tuple(sorted(key.atom_indices))
        if atoms in constrained:
            length = bonds.potentials[potential].parameters["length"].m_as(unit.angstrom)
            assert constrained[atoms] == pytest.approx(length, abs=1e-6)


def test_two_site_force_fields_each_keep_their_own_terms() -> None:
    from openff.toolkit import ForceField, Topology
    from openff.units import unit

    protein = _acylated_lysine()
    mp.attach(
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
    older = ParameterSet(("openff-2.2.1.offxml",), SAGE.charges)
    interchange, record = parameterize(protein, {"generic": SAGE, "glycan": older})
    sites = {site.type: set(site.atom_indices) for site in record.sites}
    molecule = to_openff(protein)
    molecule.partial_charges = _charges(interchange) * unit.elementary_charge

    for site_type, parameter_set in (("generic", SAGE), ("glycan", older)):
        label = "+".join(parameter_set.force_field)
        force_field = ForceField(*parameter_set.force_field)
        own = force_field.create_interchange(
            Topology.from_molecules([molecule]), charge_from_molecules=[molecule]
        )
        for name in ("Bonds", "Angles", "ProperTorsions"):
            expected = {
                (key.atom_indices, getattr(key, "mult", None)): own[name]
                .potentials[potential]
                .parameters
                for key, potential in own[name].key_map.items()
                if set(key.atom_indices) & sites[site_type]
            }
            spliced = {
                (key.atom_indices, getattr(key, "mult", None)): interchange[name]
                .potentials[potential]
                .parameters
                for key, potential in interchange[name].key_map.items()
                if potential.id.startswith(f"site:{label}::")
            }
            assert spliced.keys() == expected.keys(), (site_type, name)
            for term, parameters in expected.items():
                assert spliced[term] == parameters, (site_type, name, term)


def _clicked_noncanonical():
    protein = mp.load_protein(TRPCAGE)
    mp.mutate(protein, "A:Y3", {"ccd_file": str(AZIDOPHENYLALANINE)}, relax=False)
    mp.attach(
        protein,
        {"smiles": "C#Cc1ccccc1", "resname": "ALK"},
        type="generic",
        site="A:4II3",
        atom_name="N3",
        reaction="azide-alkyne triazole",
        merge=True,
        relax=False,
    )
    return protein


def test_a_site_merging_two_types_with_one_parameter_set_uses_it() -> None:
    protein = _clicked_noncanonical()
    interchange, record = parameterize(protein, {"noncanonical": SAGE, "generic": SAGE})
    (site,) = record.sites
    assert site.force_field == ["openff-2.3.0.offxml"]
    assert _charges(interchange)[site.atom_indices].sum() == pytest.approx(site.target, abs=1e-6)


def test_a_site_merging_two_types_with_different_parameter_sets_is_refused() -> None:
    older = ParameterSet(("openff-2.2.1.offxml",), SAGE.charges)
    with pytest.raises(ValueError, match="merges a generic and a noncanonical modification"):
        parameterize(_clicked_noncanonical(), {"noncanonical": SAGE, "generic": older})
