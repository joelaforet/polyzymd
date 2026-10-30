"""GLYCAM06j glycan sites: native terms and per-pair 1-4 scales."""

import re
from pathlib import Path

import numpy as np
import pytest

from polyzymd import modifying_proteins as mp
from polyzymd.modifying_proteins import glycam
from polyzymd.modifying_proteins.parameters import (
    ParameterizedProtein,
    ParameterSet,
    parameterize,
)

ROOT = Path(__file__).resolve().parents[2]
TRPCAGE = ROOT / "examples" / "quickstart" / "trpcage.pdb"
CHITOBIOSE = ROOT / "tests" / "data" / "modifying_proteins" / "glycam_G42666HT.pdb"

GLYCAM = ParameterSet(["GLYCAM_06j-1.xml"], {"method": "library"})


def _glycosylate(protein, site, atom_name, reaction, pdb_path=CHITOBIOSE):
    mp.attach(
        protein,
        {"pdb_path": str(pdb_path)},
        type="glycan",
        site=site,
        atom_name=atom_name,
        reaction=reaction,
        fragment_atom_name="C1",
        fragment_resnum=2,
    )


def _glycoprotein():
    """Trp-cage with an N-glycan on Asn5 (from Q5N) and an O-glycan on Ser13."""
    protein = mp.load_protein(TRPCAGE)
    mp.mutate(protein, "A:Q5", "ASN")
    _glycosylate(protein, "A:N5", "ND2", "N-glycosylation")
    _glycosylate(protein, "A:S13", "OG", "O-glycosylation")
    return protein


def _parameterized(protein):
    result = parameterize(protein, {"glycan": GLYCAM}, max_charge_shift=0.005)
    return result.interchange_without_glycam_pairs, result.record, result.to_openmm()


@pytest.fixture(scope="module")
def glycosylated():
    return _parameterized(_glycoprotein())


def _native_system(interchange, record):
    """ff14SB + GLYCAM06j straight from OpenMM's force field, on the whole glycoprotein."""
    from openmm import app

    topology = interchange.to_openmm_topology()
    linkages = {}
    for site in record.glycam:
        residue, template = site.linkage.split(" -> ")
        chain, name = residue.split(":")
        linkages[(chain, re.sub(r"\D", "", name))] = template
    names = glycam.glycam_residue_names()
    templates = {}
    for residue in topology.residues():
        key = (residue.chain.id, str(residue.id))
        if residue.name in ("ASN", "SER", "THR") and key in linkages:
            residue.name = linkages[key]
            templates[residue] = residue.name
        elif residue.name in names:
            templates[residue] = residue.name
    return app.ForceField(*glycam.OPENMM_FILES).createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        removeCMMotion=False,
        residueTemplates=templates,
    )


def _energies(system, positions, atoms=None) -> dict[str, float]:
    """Energy of each force type; with ``atoms``, only bonded terms touching them."""
    import openmm
    from openmm import unit

    copy = openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(system))
    for index in reversed(range(copy.getNumConstraints())):
        copy.removeConstraint(index)
    forces = {}
    for group, force in enumerate(copy.getForces()):
        force.setForceGroup(group)
        if atoms is not None:
            _keep_terms_touching(force, atoms)
        forces[group] = type(force).__name__
    context = openmm.Context(
        copy, openmm.VerletIntegrator(1e-3), openmm.Platform.getPlatformByName("Reference")
    )
    context.setPositions(positions)
    energies: dict[str, float] = {}
    for group, name in forces.items():
        state = context.getState(getEnergy=True, groups={group})
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        energies[name] = energies.get(name, 0.0) + energy
    return energies


def _keep_terms_touching(force, atoms) -> None:
    import openmm

    if isinstance(force, openmm.HarmonicAngleForce):
        for n in range(force.getNumAngles()):
            *indices, angle, k = force.getAngleParameters(n)
            if not set(indices) & atoms:
                force.setAngleParameters(n, *indices, angle, 0.0)
    elif isinstance(force, openmm.PeriodicTorsionForce):
        for n in range(force.getNumTorsions()):
            *indices, periodicity, phase, k = force.getTorsionParameters(n)
            if not set(indices) & atoms:
                force.setTorsionParameters(n, *indices, periodicity, phase, 0.0)


def test_glycam_sites_record_linkage_templates_and_unscaled_pairs(glycosylated) -> None:
    _, record, _ = glycosylated
    assert [site.linkage for site in record.glycam] == ["A:ASN5 -> NLN", "A:SER13 -> OLS"]
    assert all(site.residue_templates[1:] == ["4YB", "0YB"] for site in record.glycam)
    for site in record.glycam:
        # GLYCAM leaves about 70 1-4 pairs of a disaccharide unscaled.
        assert site.n_unscaled_pairs > 50
        assert site.n_pairs == len(site.pairs) > site.n_unscaled_pairs
    assert [site.shift_per_atom for site in record.sites] == [None, None]
    assert abs(record.total_charge - round(record.total_charge)) < 1e-6
    assert "pairs" not in record.to_dict()["glycam"][0]


def test_glycam_sites_reproduce_openmms_own_glycam06j_system(glycosylated) -> None:
    _assert_matches_native(*glycosylated)


def _assert_matches_native(interchange, record, ours) -> None:
    native = _native_system(interchange, record)
    # GLYCAM-owned atoms carry GLYCAM's atom-type masses (N is 14.01, not 14.0067).
    for site in record.glycam:
        for index in site.atom_indices:
            assert ours.getParticleMass(index) == native.getParticleMass(index)
    positions = interchange.positions.to_openmm()
    mine, theirs = _energies(ours, positions), _energies(native, positions)
    assert mine["HarmonicBondForce"] == pytest.approx(theirs["HarmonicBondForce"], rel=1e-9)
    # Charges, Lennard-Jones and every 1-4 pair, GLYCAM's unscaled ones included.
    assert mine["NonbondedForce"] == pytest.approx(theirs["NonbondedForce"], rel=1e-7)
    # GLYCAM's angles and torsions are copied exactly. Away from the glycans the
    # SMIRNOFF port of ff14SB differs from OpenMM's XML: rounded angle constants,
    # and impropers averaged over three orderings (an Asn side-chain amide, for one).
    site = {index for site in record.glycam for index in site.atom_indices}
    mine_site = _energies(ours, positions, site)
    theirs_site = _energies(native, positions, site)
    for name in ("HarmonicAngleForce", "PeriodicTorsionForce"):
        assert mine_site[name] == pytest.approx(theirs_site[name], rel=1e-10)
    assert mine["HarmonicAngleForce"] == pytest.approx(theirs["HarmonicAngleForce"], rel=1e-4)


def test_a_glycam_result_exports_only_with_its_pairs(glycosylated) -> None:
    interchange, record, ours = glycosylated
    result = ParameterizedProtein(interchange, record)
    with pytest.raises(ValueError, match="to_openmm"):
        _ = result.interchange
    with pytest.raises(ValueError, match="to_openmm"):
        _, _ = result
    with pytest.raises(ValueError, match="repartitioning") as refused:
        result.to_openmm(hydrogen_mass=3.0)
    assert "not supported yet" in str(refused.value)
    assert "https://github.com/joelaforet/polyzymd/issues" in str(refused.value)
    positions = interchange.positions.to_openmm()
    exported = _energies(result.to_openmm(), positions)["NonbondedForce"]
    assert exported == pytest.approx(_energies(ours, positions)["NonbondedForce"], rel=1e-12)


def test_thr_sites_proline_neighbours_and_adjacent_glycans_match_glycam06j() -> None:
    protein = mp.load_protein(TRPCAGE)
    mp.mutate(protein, "A:L7", "THR")
    mp.mutate(protein, "A:G11", "ASN")  # Pro12 follows
    _glycosylate(protein, "A:T7", "OG1", "O-glycosylation")
    _glycosylate(protein, "A:N11", "ND2", "N-glycosylation")
    _glycosylate(protein, "A:S13", "OG", "O-glycosylation")
    _glycosylate(protein, "A:S14", "OG", "O-glycosylation")  # next to Ser13's glycan
    interchange, record, system = _parameterized(protein)
    assert [site.linkage.split(" -> ")[1] for site in record.glycam] == [
        "OLT",
        "NLN",
        "OLS",
        "OLS",
    ]
    _assert_matches_native(interchange, record, system)


def test_glycans_on_two_chains_with_clashing_residue_numbers_match_glycam06j(tmp_path) -> None:
    lines = [line for line in TRPCAGE.read_text().splitlines() if line.startswith("ATOM")]
    second = []
    for line in lines:
        x = float(line[30:38]) + 30.0
        second.append(f"{line[:21]}B{line[22:30]}{x:8.3f}{line[38:]}")
    pdb = tmp_path / "two_chains.pdb"
    pdb.write_text("\n".join(lines + ["TER"] + second + ["TER", "END"]) + "\n")
    protein = mp.load_protein(pdb)
    mp.mutate(protein, "B:Q5", "ASN")
    _glycosylate(protein, "B:N5", "ND2", "N-glycosylation")
    _glycosylate(protein, "A:S13", "OG", "O-glycosylation")
    interchange, record, system = _parameterized(protein)
    # Each chain numbers its glycan residues after its own last residue.
    assert [site.linkage for site in record.glycam] == ["B:ASN5 -> NLN", "A:SER13 -> OLS"]
    _assert_matches_native(interchange, record, system)


def test_without_the_pair_overrides_the_glycan_1_4_energy_is_wrong(glycosylated) -> None:
    interchange, record, ours = glycosylated
    uniform = interchange.to_openmm(combine_nonbonded_forces=True)
    positions = interchange.positions.to_openmm()
    with_pairs = _energies(ours, positions)["NonbondedForce"]
    without = _energies(uniform, positions)["NonbondedForce"]
    assert abs(with_pairs - without) > 10.0  # kJ/mol; GLYCAM's unscaled pairs


# ---------------------------------------------------------------------------
# Inputs GLYCAM06j cannot use
# ---------------------------------------------------------------------------
def _edited_glycan(tmp_path: Path, edit) -> Path:
    lines = CHITOBIOSE.read_text().splitlines(keepends=True)
    path = tmp_path / "glycan.pdb"
    path.write_text("".join(edit(lines)))
    return path


def test_a_glycan_not_named_for_glycam_is_refused_by_glycam(tmp_path) -> None:
    path = _edited_glycan(tmp_path, lambda lines: [line.replace("0YB", "XYZ") for line in lines])
    protein = mp.load_protein(TRPCAGE)
    _glycosylate(protein, "A:S13", "OG", "O-glycosylation", path)
    with pytest.raises(ValueError, match="GLYCAM06j does not define: XYZ"):
        parameterize(protein, {"glycan": GLYCAM}, max_charge_shift=0.005)


def test_a_glycan_on_a_terminal_residue_is_refused_by_glycam() -> None:
    protein = mp.load_protein(TRPCAGE)
    _glycosylate(protein, "A:S20", "OG", "O-glycosylation")
    with pytest.raises(ValueError) as refused:
        parameterize(protein, {"glycan": GLYCAM}, max_charge_shift=0.005)
    message = str(refused.value)
    assert message.startswith("Refused: glycan G001 on A:SER20")
    assert "linkage residue (OLS)" in message
    assert "a neighbouring residue on each side" in message
    assert "https://github.com/joelaforet/polyzymd/issues" in message


def test_glycam_parameterizes_glycans_only() -> None:
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
    with pytest.raises(ValueError, match="glycan sites only"):
        parameterize(protein, {"generic": GLYCAM}, max_charge_shift=0.005)


def test_a_glycam_parameter_set_takes_library_charges() -> None:
    with pytest.raises(ValueError, match="charges=\\{'method': 'library'\\}"):
        ParameterSet(["GLYCAM_06j-1.xml"], {"method": "nagl", "model": "x.pt"})


# ---------------------------------------------------------------------------
# What the workarounds above rely on. When one of these tests fails, the
# upstream gap is closed and the workaround can go (issue #281).
# ---------------------------------------------------------------------------
def test_interchange_holds_one_1_4_scale_per_collection() -> None:
    from openff.toolkit import ForceField, Molecule

    interchange = ForceField("openff-2.3.0.offxml").create_interchange(
        Molecule.from_smiles("CCO").to_topology()
    )
    assert isinstance(interchange["Electrostatics"].scale_14, float)
    assert not hasattr(interchange["Electrostatics"], "pair_scales")


def test_interchange_gromacs_export_drops_an_improper_out_of_its_own_order(tmp_path) -> None:
    from openff.interchange.components.potentials import Potential
    from openff.interchange.models import ImproperTorsionKey, PotentialKey
    from openff.toolkit import ForceField, Molecule
    from openff.units import unit

    molecule = Molecule.from_smiles("CC(=O)NC")
    molecule.generate_conformers(n_conformers=1)
    interchange = ForceField("openff-2.3.0.offxml").create_interchange(molecule.to_topology())
    interchange.box = [3, 3, 3] * unit.nanometer
    impropers = interchange["ImproperTorsions"]
    central_first = next(iter(impropers.key_map))
    center, i, k, l_ = central_first.atom_indices  # SMIRNOFF keys put the center first
    amber_order = ImproperTorsionKey(atom_indices=(i, k, center, l_), mult=0)
    potential_key = PotentialKey(id="amber-order", mult=0, associated_handler="ImproperTorsions")
    impropers.key_map[amber_order] = potential_key
    impropers.potentials[potential_key] = Potential(
        parameters={
            "k": 1.234 * unit.kilocalorie_per_mole,  # 5.163056 kJ/mol, unique here
            "periodicity": 2 * unit.dimensionless,
            "phase": 180.0 * unit.degree,
            "idivf": 1.0 * unit.dimensionless,
        }
    )
    interchange.to_gromacs(prefix=str(tmp_path / "m"), monolithic=True)
    assert "5.163056" not in (tmp_path / "m.top").read_text()
