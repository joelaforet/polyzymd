"""The audit of covalently modified proteins: correct builds pass, each known mistake fails.

Every system is trp-cage with real modifications (an acylated Lys8; N- and
O-linked chitobiose on GLYCAM06j; both together). Each is realized for OpenMM
or GROMACS through one helper and checked against the protein the audit
rebuilds from the input PDB and the modification blocks.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("openff.interchange")

import openmm  # noqa: E402

from polyzymd import modifying_proteins as mp  # noqa: E402
from polyzymd.audit.findings import run_checks  # noqa: E402
from polyzymd.audit.realized import RealizedSystem  # noqa: E402
from polyzymd.modifying_proteins import glycam  # noqa: E402
from polyzymd.modifying_proteins.components import protein_component  # noqa: E402
from polyzymd.modifying_proteins.parameters import (  # noqa: E402
    ParameterSet,
    parameterize,
    to_openff,
)

ROOT = Path(__file__).resolve().parents[2]
TRPCAGE = ROOT / "examples" / "quickstart" / "trpcage.pdb"
CHITOBIOSE = ROOT / "tests" / "data" / "modifying_proteins" / "glycam_G42666HT.pdb"
AZIDOPHENYLALANINE = ROOT / "tests" / "data" / "modifying_proteins" / "4II.cif"

MODEL = "openff-gnn-am1bcc-1.0.0.pt"
SAGE = ParameterSet(("openff-2.3.0.offxml",), {"method": "nagl", "model": MODEL})
GLYCAM = ParameterSet(["GLYCAM_06j-1.xml"], {"method": "library"})
ENGINES = ("openmm", "gromacs")

#: The parameter checks of every domain and the checks of modified proteins.
CHECKS = (
    "parameters.bonded",
    "parameters.nonbonded",
    "parameters.one_four",
    "parameters.constraints",
    "charges.source",
    "masses",
    "parameters.coverage",
    "charges.site_integer",
    "charges.site_shift",
    "modifications.applied",
    "parameters.modified_backbone",
    "parameters.glycan_force_field",
    "modifications.protonation_state",
    "modifications.sites",
    "modifications.declared",
    "modifications.chemistry",
)


def _acylate(protein, site="A:K8"):
    mp.deprotonate(protein, site, "NZ")
    mp.attach(
        protein,
        {"smiles": "*C(=O)CCCCCCC", "resname": "OC8"},
        type="generic",
        site=site,
        atom_name="NZ",
        relax=False,
    )


def _glycosylate(protein):
    for site, atom, reaction in (
        ("A:N5", "ND2", "N-glycosylation"),
        ("A:S13", "OG", "O-glycosylation"),
    ):
        mp.attach(
            protein,
            {"pdb_path": str(CHITOBIOSE)},
            type="glycan",
            site=site,
            atom_name=atom,
            reaction=reaction,
            fragment_atom_name="C1",
            fragment_resnum=2,
            relax=False,
        )


def _build(kind: str, sets: dict, build_sets: dict | None = None, max_charge_shift=0.005):
    """A modified trp-cage, its Interchange and record, and its component record."""
    protein = mp.load_protein(TRPCAGE)
    if kind in ("glycans", "mixed"):
        mp.mutate(protein, "A:Q5", "ASN", relax=False, minimize=False)
    if kind in ("acyl", "mixed"):
        _acylate(protein)
    if kind in ("glycans", "mixed"):
        _glycosylate(protein)
    if kind == "asp":
        mp.mutate(protein, "A:S13", "ASP", relax=False, minimize=False)
    if kind == "click":
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
    result = parameterize(protein, build_sets or sets, max_charge_shift=max_charge_shift)
    return protein, result, result.record, protein_component(protein, result.record, sets)


def _realize(built, engine: str, folder: Path, edit=None, pairs: bool = True) -> RealizedSystem:
    """Realize a build for ``engine`` through its public export.

    ``pairs=False`` exports the bare Interchange instead, which scales every 1-4
    pair uniformly; ``edit`` breaks the OpenMM system before it is written.
    """
    from openff.units import unit

    protein, result, _, _ = built
    bare = result.interchange_without_glycam_pairs
    folder.mkdir(parents=True, exist_ok=True)
    if engine == "openmm":
        if pairs:
            system = result.to_openmm()
        else:
            system = bare.to_openmm(combine_nonbonded_forces=True)
        if edit is not None:
            edit(system, protein)
        (folder / "system.xml").write_text(openmm.XmlSerializer.serialize(system))
        protein.save_pdb(folder / "protein.pdb")
        return RealizedSystem.from_openmm(folder / "system.xml", folder / "protein.pdb")
    bare.box = np.eye(3) * 6 * unit.nanometer
    if pairs:
        result.to_gromacs(prefix=str(folder / "sys"), monolithic=True)
    else:
        bare.to_gromacs(prefix=str(folder / "sys"), monolithic=True)
    (folder / "run.mdp").write_text(
        "coulombtype = PME\nrcoulomb = 1.0\nrvdw = 1.0\nvdwtype = cut-off\n"
        "vdw-modifier = None\nDispCorr = no\nconstraints = h-bonds\ndt = 0.002\n"
    )
    return RealizedSystem.from_gromacs(folder / "sys.top", folder / "sys.gro", folder / "run.mdp")


def _status(realized, component, declared=None) -> dict[str, tuple[str, str]]:
    return {
        f.check: (f.status, f.detail)
        for f in run_checks(realized, [component], CHECKS, declared=declared)
        if f.component in (component.name, "system")
    }


def _not_passing(status: dict, allowed=()) -> dict:
    """Findings other than PASS (and the INFO of a comparison not asked for)."""
    return {k: v for k, v in status.items() if v[0] not in ("PASS", "INFO") and k not in allowed}


def _declared(component, **changes) -> dict:
    """The component's own blocks and parameter sets, as a config would declare them."""
    import copy

    recipe = component.domains[0].reference["modified_protein"]
    declared = {
        "spec": copy.deepcopy(recipe["spec"]),
        "parameter_sets": copy.deepcopy(recipe["parameter_sets"]),
        "base_dir": recipe["base_dir"],
    }
    for key, edit in changes.items():
        edit(declared)
    return {component.name: declared}


@pytest.fixture(scope="module")
def acylated():
    return _build("acyl", {"generic": SAGE})


@pytest.fixture(scope="module")
def glycosylated():
    return _build("glycans", {"glycan": GLYCAM})


@pytest.fixture(scope="module")
def mixed():
    return _build("mixed", {"glycan": GLYCAM, "generic": SAGE})


@pytest.mark.parametrize("engine", ENGINES)
def test_an_acylated_lysine_passes_every_check(acylated, engine, tmp_path) -> None:
    status = _status(_realize(acylated, engine, tmp_path), acylated[3])
    assert _not_passing(status) == {}


@pytest.mark.parametrize("engine", ENGINES)
def test_glycam_glycans_and_a_sage_site_pass_every_check(mixed, engine, tmp_path) -> None:
    status = _status(_realize(mixed, engine, tmp_path), mixed[3])
    assert _not_passing(status) == {}
    assert status["parameters.glycan_force_field"][1] == "2 glycans on GLYCAM06j"


def test_glycam_glycans_carry_glycam_masses(glycosylated, tmp_path) -> None:
    status = _status(_realize(glycosylated, "openmm", tmp_path), glycosylated[3])
    assert status["masses"][0] == "PASS", status["masses"][1]


def test_a_click_on_a_noncanonical_residue_passes_every_check(tmp_path) -> None:
    """Azidophenylalanine (CCD 4II) clicked to an alkyne inside its residue (merge)."""
    built = _build("click", {"noncanonical": SAGE, "generic": SAGE})
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert _not_passing(status) == {}


def test_a_chemistry_the_build_perceives_wrongly_fails(monkeypatch, tmp_path) -> None:
    """The build reads the acyl amide as its zwitterionic resonance form; Pablo does not."""
    from openff.toolkit import Molecule

    from polyzymd.modifying_proteins import parameters

    read = parameters.to_openff

    def zwitterion(protein):
        molecule = read(protein)
        keys = protein.atom_keys()
        if not any(k.resnum != 8 and k.atom_name == "C1" for k in keys):
            return molecule  # the ff14SB reference, which has no acyl group
        nz = _atom(protein, 8, "NZ")
        c1 = next(i for i, k in enumerate(keys) if k.resnum != 8 and k.atom_name == "C1")
        o1 = next(i for i, k in enumerate(keys) if k.resnum != 8 and k.atom_name == "O1")
        rdmol = molecule.to_rdkit()
        rdmol.GetBondBetweenAtoms(c1, o1).SetBondType(
            rdmol.GetBondBetweenAtoms(c1, nz).GetBondType()
        )
        from rdkit import Chem

        rdmol.GetBondBetweenAtoms(c1, nz).SetBondType(Chem.BondType.DOUBLE)
        rdmol.GetAtomWithIdx(o1).SetFormalCharge(-1)
        rdmol.GetAtomWithIdx(nz).SetFormalCharge(1)
        changed = Molecule.from_rdkit(
            rdmol, allow_undefined_stereo=True, hydrogens_are_explicit=True
        )
        for atom, original in zip(changed.atoms, molecule.atoms):
            atom.name = original.name
            atom.metadata.update(original.metadata)
        return changed

    monkeypatch.setattr(parameters, "to_openff", zwitterion)
    built = _build("acyl", {"generic": SAGE}, max_charge_shift=0.05)
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert status["modifications.chemistry"][0] == "FAIL"
    assert "formal charges differ" in status["modifications.chemistry"][1]
    assert "bonds differ" in status["modifications.chemistry"][1]


def test_another_kekule_form_of_an_aromatic_ring_is_the_same_chemistry(
    monkeypatch, tmp_path
) -> None:
    """The build and Pablo may place a ring's double bonds differently; that is not a
    difference (the FRET donor's DBCO and 4II rings did so)."""
    from polyzymd.modifying_proteins import parameters

    read = parameters.to_openff

    def other_kekule_form(protein):
        molecule = read(protein)
        rdmol = molecule.to_rdkit()
        ring = next(
            bonds
            for bonds in rdmol.GetRingInfo().BondRings()
            if len(bonds) == 6
            and all(rdmol.GetBondWithIdx(b).GetIsAromatic() for b in bonds)
            and all(
                molecule.atom(rdmol.GetBondWithIdx(b).GetBeginAtomIdx()).metadata["residue_number"]
                == "3"
                for b in bonds
            )
        )
        ends = {
            tuple(
                sorted(
                    (
                        rdmol.GetBondWithIdx(b).GetBeginAtomIdx(),
                        rdmol.GetBondWithIdx(b).GetEndAtomIdx(),
                    )
                )
            )
            for b in ring
        }
        for bond in molecule.bonds:
            if tuple(sorted((bond.atom1_index, bond.atom2_index))) in ends:
                bond._bond_order = 3 - bond.bond_order
        return molecule

    monkeypatch.setattr(parameters, "to_openff", other_kekule_form)
    built = _build("click", {"noncanonical": SAGE, "generic": SAGE})
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert status["modifications.chemistry"][0] == "PASS", status["modifications.chemistry"]


def test_an_attached_amino_acid_is_read_as_an_attachment(tmp_path) -> None:
    """A CCD amino acid (4II) attached by its CB to Lys8: mBuild numbers it after the
    chain, and it has N, CA and C, but it is an attachment, not a chain residue."""
    protein = mp.load_protein(TRPCAGE, ccd_paths=[AZIDOPHENYLALANINE.parent])
    _acylate_with = {"ccd": "4II", "link_atom": "CB"}
    mp.deprotonate(protein, "A:K8", "NZ")
    mp.attach(protein, _acylate_with, type="generic", site="A:K8", atom_name="NZ", relax=False)
    sets = {"generic": SAGE}
    result = parameterize(protein, sets)
    built = (protein, result, result.record, protein_component(protein, result.record, sets))
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert _not_passing(status) == {}


def test_pablo_reads_a_protein_of_two_chains(tmp_path) -> None:
    """Two chains no bond joins are two molecules to Pablo; the audit reads them as one
    protein, atom by atom (acylated Lys8 on chain A, a glycan on chain B)."""
    from polyzymd.audit import _pablo
    from polyzymd.audit.modified import chemistry_differences

    lines = [
        line for line in TRPCAGE.read_text().splitlines() if line.startswith(("ATOM", "HETATM"))
    ]
    chain_a = [line[:21] + "A" + line[22:] for line in lines]
    chain_b = [
        line[:21] + "B" + line[22:30] + f"{float(line[30:38]) + 30.0:8.3f}" + line[38:]
        for line in lines
    ]
    pdb = tmp_path / "two_chains.pdb"
    pdb.write_text("\n".join(chain_a + ["TER"] + chain_b + ["TER", "END"]) + "\n")
    protein = mp.load_protein(pdb)
    mp.mutate(protein, "B:Q5", "ASN", relax=False, minimize=False)
    _acylate(protein)
    mp.attach(
        protein,
        {"pdb_path": str(CHITOBIOSE)},
        type="glycan",
        site="B:N5",
        atom_name="ND2",
        reaction="N-glycosylation",
        fragment_atom_name="C1",
        fragment_resnum=2,
        relax=False,
    )
    read = _pablo.perceive(protein)
    assert read.n_atoms == len(protein.atom_keys())
    assert chemistry_differences(read, to_openff(protein), protein.atom_keys()) == []


def test_the_component_record_names_segments_domains_and_the_rebuild(mixed) -> None:
    component = mixed[3]
    data = component.to_dict()
    assert [d["name"] for d in data["domains"]] == ["protein", "M001", "G001", "G002"]
    assert set(data["segments"]) == {"M001", "G001", "G002"}
    recipe = component.domains[0].reference["modified_protein"]
    assert recipe["spec"]["mutations"][0]["site"] == "A:Q5"
    assert recipe["parameter_sets"]["glycan"]["force_field"] == ["GLYCAM_06j-1.xml"]


@pytest.mark.parametrize("engine", ENGINES)
def test_glycam_pairs_scaled_uniformly_fail_one_four(glycosylated, engine, tmp_path) -> None:
    """The mistake an export through Interchange makes: one global 1-4 scale."""
    status = _status(_realize(glycosylated, engine, tmp_path, pairs=False), glycosylated[3])
    assert status["parameters.one_four"][0] == "FAIL"
    assert "G001" in status["parameters.one_four"][1]
    # The bare Interchange also lacks GLYCAM's atom-type masses, which only the export adds.
    assert status["masses"][0] == "FAIL"
    assert _not_passing(status, allowed=("parameters.one_four", "masses")) == {}


def test_nagl_charges_on_the_backbone_fail(acylated, tmp_path) -> None:
    """The earlier rule: the whole residue (backbone included) charged by NAGL."""
    from polyzymd.audit.modified import capped_model_charges

    protein, _, record, component = acylated
    keys = protein.atom_keys()
    residue = [i for i, k in enumerate(keys) if k.resnum == 8]
    site = sorted(set(residue) | set(record.sites[0].atom_indices))
    nagl = capped_model_charges(to_openff(protein), site, [], {"charges": {"model": MODEL}})

    def edit(system, _protein):
        force = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
        total = sum(nagl.values())
        shift = (round(total) - total) / len(nagl)
        for i, q in nagl.items():
            _, sigma, epsilon = force.getParticleParameters(i)
            force.setParticleParameters(i, q + shift, sigma, epsilon)

    status = _status(_realize(acylated, "openmm", tmp_path, edit), component)
    assert status["charges.source"][0] == "FAIL"
    assert status["parameters.modified_backbone"][0] == "FAIL"
    assert "backbone charges" in status["parameters.modified_backbone"][1]


def test_sage_phi_and_psi_on_a_modified_residue_fail(acylated, tmp_path) -> None:
    from openff.toolkit import ForceField

    protein, interchange, _, component = acylated
    keys = protein.atom_keys()
    at = {(k.resnum, k.atom_name): i for i, k in enumerate(keys)}
    dihedrals = {
        (at[(7, "C")], at[(8, "N")], at[(8, "CA")], at[(8, "C")]),  # phi
        (at[(8, "N")], at[(8, "CA")], at[(8, "C")], at[(9, "N")]),  # psi
    }
    molecule = to_openff(protein)
    sage = ForceField("openff-2.3.0.offxml").create_interchange(molecule.to_topology())
    sage_system = sage.to_openmm(combine_nonbonded_forces=True)
    sage_torsions = next(
        f for f in sage_system.getForces() if isinstance(f, openmm.PeriodicTorsionForce)
    )

    def edit(system, _protein):
        force = next(f for f in system.getForces() if isinstance(f, openmm.PeriodicTorsionForce))
        for n in range(force.getNumTorsions()):
            *atoms, periodicity, phase, k = force.getTorsionParameters(n)
            if tuple(atoms) in dihedrals or tuple(atoms[::-1]) in dihedrals:
                force.setTorsionParameters(n, *atoms, periodicity, phase, 0.0)
        for n in range(sage_torsions.getNumTorsions()):
            *atoms, periodicity, phase, k = sage_torsions.getTorsionParameters(n)
            if tuple(atoms) in dihedrals or tuple(atoms[::-1]) in dihedrals:
                force.addTorsion(*atoms, periodicity, phase, k)

    status = _status(_realize(acylated, "openmm", tmp_path, edit), component)
    assert status["parameters.bonded"][0] == "FAIL"
    assert status["parameters.modified_backbone"][0] == "FAIL"


def _replace_bonds(system, keep=lambda i, j: True, add=()):
    for n, force in enumerate(system.getForces()):
        if isinstance(force, openmm.HarmonicBondForce):
            bonds = openmm.HarmonicBondForce()
            for b in range(force.getNumBonds()):
                i, j, r, k = force.getBondParameters(b)
                if keep(i, j):
                    bonds.addBond(i, j, r, k)
            for i, j in add:
                bonds.addBond(i, j, 0.15, 1000.0)
            system.removeForce(n)
            system.addForce(bonds)
            return


def _atom(protein, resnum, name):
    return next(
        i for i, k in enumerate(protein.atom_keys()) if (k.resnum, k.atom_name) == (resnum, name)
    )


def test_a_dropped_crosslink_bond_fails_modifications_applied(acylated, tmp_path) -> None:
    protein, _, _, component = acylated
    nz = _atom(protein, 8, "NZ")
    c1 = next(i for i, k in enumerate(protein.atom_keys()) if k.resnum != 8 and k.atom_name == "C1")

    def edit(system, _protein):
        _replace_bonds(system, keep=lambda i, j: {i, j} != {nz, c1})

    status = _status(_realize(acylated, "openmm", tmp_path, edit), component)
    assert status["modifications.applied"][0] == "FAIL"
    assert "planned bonds are missing" in status["modifications.applied"][1]


def test_a_bond_no_change_planned_fails_modifications_applied(acylated, tmp_path) -> None:
    """A bond inside the molecule that no change made: only the new logic can see it."""
    protein, _, _, component = acylated
    nz = _atom(protein, 8, "NZ")
    cb = _atom(protein, 8, "CB")

    def edit(system, _protein):
        _replace_bonds(system, add=[(nz, cb)])

    status = _status(_realize(acylated, "openmm", tmp_path, edit), component)
    assert status["modifications.applied"][0] == "FAIL"
    assert "bonds are not planned" in status["modifications.applied"][1]
    assert status["parameters.coverage"][0] == "PASS"


def test_a_term_touching_two_sites_fails_coverage(mixed, tmp_path) -> None:
    protein, _, record, component = mixed
    first = {site.name: site.atom_indices[0] for site in record.sites}

    def edit(system, _protein):
        _replace_bonds(system, add=[(first["M001"], first["G001"])])

    status = _status(_realize(mixed, "openmm", tmp_path, edit), component)
    assert status["parameters.coverage"][0] == "FAIL"
    assert "touch two modification sites" in status["parameters.coverage"][1]


def test_a_site_that_swallows_another_residue_fails(monkeypatch, tmp_path) -> None:
    """The build and the audit share the site derivation; the audit's own rule catches it."""
    from polyzymd.modifying_proteins import parameters, topology

    derive = topology.modification_sites

    def swallowing(protein):
        keys = protein.atom_keys()
        extra = {i for i, k in enumerate(keys) if k.resnum == 2}
        out = []
        for site in derive(protein):
            if site.name == "M001":
                atoms = tuple(sorted(set(site.atom_indices) | extra))
                site = dataclasses.replace(
                    site, atom_indices=atoms, atoms=tuple(keys[i] for i in atoms)
                )
            out.append(site)
        return out

    monkeypatch.setattr(topology, "modification_sites", swallowing)
    monkeypatch.setattr(parameters, "modification_sites", swallowing)
    built = _build("acyl", {"generic": SAGE})
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert status["modifications.sites"][0] == "FAIL"
    assert "not atoms any change defines" in status["modifications.sites"][1]
    assert status["parameters.modified_backbone"][0] == "FAIL"


def test_a_formal_charge_the_records_do_not_explain_fails_site_integer(
    monkeypatch, tmp_path
) -> None:
    """The chemistry's formal charge must equal the input residues' plus the recorded changes."""
    from polyzymd.audit import modified

    rebuild = modified.rebuild

    def misrecorded(recipe):
        protein = rebuild(recipe)
        attach = next(r for r in protein.records if r.kind == "attach")
        index = protein.records.index(attach)
        protein.records[index] = dataclasses.replace(attach, formal_charge_change=1)
        return protein

    monkeypatch.setattr(modified, "rebuild", misrecorded)
    built = acylated_build()
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert status["charges.site_integer"][0] == "FAIL"
    assert "from the input residues and the recorded changes" in status["charges.site_integer"][1]


def acylated_build():
    return _build("acyl", {"generic": SAGE})


def _terminus(change: str, max_charge_shift: float = 0.005):
    """Trp-cage with its N-terminus neutralized or acetylated, or its C-terminus amidated."""
    protein = mp.load_protein(TRPCAGE)
    if change in ("neutral", "acetyl"):
        mp.deprotonate(protein, "A:N1", "N")
    if change == "acetyl":
        mp.attach(
            protein,
            {"smiles": "*C(C)=O", "resname": "ACE"},
            type="generic",
            site="A:N1",
            atom_name="N",
            relax=False,
        )
    if change == "amide":
        mp.attach(
            protein,
            {"smiles": "*NC", "resname": "NME"},
            type="generic",
            site="A:S20",
            atom_name="C",
            leaving_atom_names=["OXT"],
            relax=False,
        )
    sets = {"generic": SAGE, "noncanonical": SAGE}
    result = parameterize(protein, sets, max_charge_shift=max_charge_shift)
    return protein, result, result.record, protein_component(protein, result.record, sets)


@pytest.mark.parametrize("change", ["acetyl", "amide"])
def test_a_modified_terminus_passes_every_check(change, tmp_path) -> None:
    """A label or amide on a terminus: the residue's other atoms take the main-chain template."""
    built = _terminus(change)
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert _not_passing(status) == {}


@pytest.mark.parametrize("change", ["acetyl", "amide"])
def test_a_labelled_terminus_keeps_the_rest_of_its_backbone(change, tmp_path) -> None:
    """The touched backbone atoms join the site by design; the rest stay ff14SB's."""
    built = _terminus(change)
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert status["parameters.modified_backbone"][0] == "PASS", status


def test_a_neutral_terminus_fails_modifications_sites(monkeypatch, tmp_path) -> None:
    """A terminus changed by protonation alone, built by bypassing the builder's refusal."""
    from polyzymd.modifying_proteins import parameters

    monkeypatch.setattr(parameters, "_refuse_neutral_termini", lambda protein, termini: None)
    built = _terminus("neutral", max_charge_shift=0.05)
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert status["modifications.sites"][0] == "FAIL"
    assert "neutral terminus is not supported" in status["modifications.sites"][1]
    assert "residue 1 (N-terminus)" in status["modifications.sites"][1]


def test_a_modified_terminus_left_uncapped_in_the_reference_fails(monkeypatch, tmp_path) -> None:
    """A reference that leaves an acetylated N-terminus on NASN's charges (not ASN's).

    The build and the audit share ``reference_protein``; the audit finds the
    modified termini itself and fails when one is not capped.
    """
    from polyzymd.modifying_proteins import parameters

    monkeypatch.setattr(parameters, "modified_termini", lambda reference, sites: {})
    # On the terminal template the site needs a larger shift than the default bound.
    built = _terminus("acetyl", max_charge_shift=0.05)
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert status["modifications.sites"][0] == "FAIL"
    assert "does not cap these modified termini" in status["modifications.sites"][1]
    assert "residue 1 (N-terminus)" in status["modifications.sites"][1]


def test_uncapped_termini_reads_the_caps_of_the_reference() -> None:
    from polyzymd.audit.modified import uncapped_termini
    from polyzymd.modifying_proteins.parameters import reference_protein, to_openff

    protein, *_ = _terminus("acetyl")
    sites = mp.modification_sites(protein)
    capped = to_openff(reference_protein(protein, sites))
    assert uncapped_termini(capped, {("A", 1, ""): "N"}) == []
    # Ser20's C-terminus is not modified, so the reference does not cap it.
    assert uncapped_termini(capped, {("A", 20, ""): "C"}) == ["chain A residue 20 (C-terminus)"]


def test_the_declared_blocks_are_compared_with_the_build(acylated, tmp_path) -> None:
    realized = _realize(acylated, "openmm", tmp_path)
    component = acylated[3]
    same = _status(realized, component, _declared(component))
    assert same["modifications.declared"][0] == "PASS"

    def other_fragment(declared):
        entry = declared["spec"]["modifications"][1]["attach"]
        entry["fragment"] = {"smiles": "*C(=O)CCCCCC", "resname": "OC7"}

    def other_charges(declared):
        declared["parameter_sets"]["generic"]["charges"]["model"] = "other.pt"

    for edit, words in ((other_fragment, "change 2"), (other_charges, "site M001")):
        status = _status(realized, component, _declared(component, edit=edit))
        assert status["modifications.declared"][0] == "FAIL"
        assert words in status["modifications.declared"][1]


def test_the_audit_site_rule_tables_match_the_modification_package() -> None:
    from polyzymd.audit import modified
    from polyzymd.modifying_proteins import topology

    assert modified.BACKBONE == topology.BACKBONE
    for atom, hydrogens in topology._BACKBONE_HYDROGENS.items():
        assert modified.TOUCHED_BACKBONE[atom] == (atom, *hydrogens)
    assert set(topology._TEMPLATE_PROTONATION) <= modified.FF14SB_VARIANTS


def test_a_shift_beyond_the_recorded_bound_fails(acylated, tmp_path) -> None:
    protein, interchange, record, _ = acylated
    tight = protein_component(
        protein, dataclasses.replace(record, max_charge_shift=1e-6), {"generic": SAGE}
    )
    status = _status(_realize(acylated, "openmm", tmp_path), tight)
    assert status["charges.site_shift"][0] == "FAIL"
    assert "M001" in status["charges.site_shift"][1]


def test_charges_from_another_nagl_model_fail(tmp_path) -> None:
    # 0.1.0-rc.2, rc.3 and 1.0.0 give the same charges; rc.1 differs.
    other = ParameterSet(
        ("openff-2.3.0.offxml",), {"method": "nagl", "model": "openff-gnn-am1bcc-0.1.0-rc.1.pt"}
    )
    built = _build("acyl", {"generic": SAGE}, build_sets={"generic": other})
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert status["charges.source"][0] == "FAIL"
    assert "M001" in status["charges.source"][1]


def test_a_glycan_on_sage_warns(tmp_path) -> None:
    built = _build("glycans", {"glycan": SAGE})
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert status["parameters.glycan_force_field"][0] == "WARN"
    assert "G001" in status["parameters.glycan_force_field"][1]
    assert _not_passing(status, allowed=("parameters.glycan_force_field",)) == {}


def test_a_neutral_aspartate_from_a_ccd_mutation_warns(tmp_path) -> None:
    built = _build("asp", {})
    status = _status(_realize(built, "openmm", tmp_path), built[3])
    assert status["modifications.protonation_state"][0] == "WARN"
    assert (
        "deprotonate: {site: A:ASP13, atom_name: OD2}"
        in (status["modifications.protonation_state"][1])
    )
    assert _not_passing(status, allowed=("modifications.protonation_state",)) == {}


def test_planned_names_are_compared_as_openmms_pdb_reader_names_them() -> None:
    from polyzymd.audit.checks.modifications import _as_read

    # An input NME cap's CCD amide hydrogen HN2 is read back as H2.
    assert _as_read("NME", "HN2") == "H2"
    assert _as_read("NME", "H1") == "H1"
    assert _as_read("ALA", "CA") == "CA"


def test_a_declared_grafted_chain_is_compared_without_its_draw() -> None:
    """The declared chain has no drawn sequence or replicate; the rest must match."""
    from pathlib import Path

    from polyzymd.audit.modified import _fragment

    declared = {"polymer": {"monomers": [{"label": "A"}], "length": 9, "force_field": []}}
    built = {
        "polymer": {
            "monomers": [{"label": "A"}],
            "length": 9,
            "force_field": ["openff-2.3.0.offxml"],
            "drawn": "HAAAAAAAH",
            "replicate": 2,
        }
    }
    assert _fragment(declared, Path(".")) == _fragment(built, Path("."))
    longer = {"polymer": {**built["polymer"], "length": 11}}
    assert _fragment(declared, Path(".")) != _fragment(longer, Path("."))
    assert _fragment(declared, Path(".")) != _fragment({"smiles": None}, Path("."))


def test_a_fragment_placed_through_the_protein_fails_the_clash_check(acylated, tmp_path) -> None:
    from polyzymd.audit.findings import run_checks

    protein, interchange, record, component = acylated
    system = interchange.to_openmm(combine_nonbonded_forces=True)
    # A box wider than the protein: OpenMM's default 2 nm box would bring the
    # protein's ends onto each other's images.
    system.setDefaultPeriodicBoxVectors(*(openmm.Vec3(*row) for row in 5.0 * np.eye(3)))
    (tmp_path / "system.xml").write_text(openmm.XmlSerializer.serialize(system))
    protein.save_pdb(tmp_path / "protein.pdb")
    keys = protein.atom_keys()
    acyl = next(
        i for i, s in enumerate(record.segments) if s != "PROT" and keys[i].atom_name[0] == "C"
    )
    target = next(i for i, k in enumerate(keys) if k.resnum == 20 and k.atom_name == "CA")
    lines = (tmp_path / "protein.pdb").read_text().splitlines()
    atoms = [k for k, line in enumerate(lines) if line.startswith(("ATOM", "HETATM"))]
    x, y, z = (float(lines[atoms[target]][c : c + 8]) for c in (30, 38, 46))
    line = lines[atoms[acyl]]
    lines[atoms[acyl]] = f"{line[:30]}{x + 0.5:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
    clashing = tmp_path / "clashing.pdb"
    clashing.write_text("\n".join(lines) + "\n")

    def status(pdb):
        realized = RealizedSystem.from_openmm(tmp_path / "system.xml", pdb)
        (finding,) = [
            f
            for f in run_checks(realized, [component], ["geometry.modified_protein_clash"])
            if f.check == "geometry.modified_protein_clash"
        ]
        return finding

    passed = status(tmp_path / "protein.pdb")
    assert passed.status == "PASS", passed.detail
    failed = status(clashing)
    assert failed.status == "FAIL" and "1 heavy-heavy pairs closer than 1.0 Å" in failed.detail
    assert ":CA" in failed.detail

    # A fragment hydrogen 0.6 Å from a protein hydrogen, and 1.7 Å from its carbon,
    # leaves every heavy atom clear: the check fails on the hydrogens alone.
    hydrogen = next(
        i for i, s in enumerate(record.segments) if s != "PROT" and keys[i].atom_name[0] == "H"
    )
    ha = next(i for i, k in enumerate(keys) if k.resnum == 20 and k.atom_name in ("HA", "HA2"))
    lines = (tmp_path / "protein.pdb").read_text().splitlines()
    ca_xyz, ha_xyz = (
        np.array([float(lines[atoms[i]][c : c + 8]) for c in (30, 38, 46)]) for i in (target, ha)
    )
    x, y, z = ha_xyz + 0.6 * (ha_xyz - ca_xyz) / np.linalg.norm(ha_xyz - ca_xyz)
    line = lines[atoms[hydrogen]]
    lines[atoms[hydrogen]] = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
    clashing = tmp_path / "hydrogen_clash.pdb"
    clashing.write_text("\n".join(lines) + "\n")
    failed = status(clashing)
    assert failed.status == "FAIL", failed.detail
    assert "1 hydrogen-hydrogen pairs closer than 0.8 Å" in failed.detail
    assert "heavy-heavy" not in failed.detail and ":HA" in failed.detail
