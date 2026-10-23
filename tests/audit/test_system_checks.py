"""System, geometry and restraint checks: pass a correct build, fail each way one goes wrong.

The small all-component system is built once, with a flat-bottom restraint
between a protein atom and a substrate atom; each test breaks a copy of its
system, its coordinates or its config and checks that the right check fails.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

pytest.importorskip("openff.interchange")
pytest.importorskip("MDAnalysis")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import openmm  # noqa: E402
from _support.small_system import write_inputs  # noqa: E402

from polyzymd.audit.findings import run_checks  # noqa: E402
from polyzymd.audit.realized import RealizedSystem  # noqa: E402

SYSTEM_CHECKS = [
    "components.declared",
    "composition.counts",
    "composition.polymer_sequences",
    "composition.polymer_mass_fraction",
    "ions.neutral",
    "ions.ionic_strength",
    "ions.salt",
    "constraints.hbonds",
    "nonbonded.settings",
    "geometry.water_overlap",
    "geometry.polymer_overlap",
    "geometry.interchain_overlap",
    "geometry.buried_water",
    "geometry.periodic_image",
    "geometry.image_bridging",
    "restraint.parameters",
    "restraint.atom_identity",
]
RESTRAINT = {
    "type": "flat_bottom",
    "name": "tether",
    "atom1": {"component": "trpcage", "selection": "resid 6 and name CA"},
    "atom2": {"component": "substrate", "selection": "name C1x"},
    "distance": 8.0,
    "force_constant": 1000.0,
}


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.config.loader import load_config
    from polyzymd.core.restraints import RestraintFactory, apply_restraints

    folder = tmp_path_factory.mktemp("system_checks")
    path = write_inputs(folder)
    data = yaml.safe_load(path.read_text())
    data["restraints"] = [RESTRAINT]
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    config = load_config(path)
    builder = SystemBuilder.from_config(config)
    builder.build_from_config(config, working_dir=folder / "build", polymer_seed=2)
    topology, system, _ = builder.get_openmm_components()
    scopes = {c.name: c.atom_indices() for c in builder.components()}
    apply_restraints([RestraintFactory.from_config(RESTRAINT)], topology, system, scopes)
    return builder, config, folder, system


def _realize(folder: Path, system, name: str, positions=None, box=None) -> RealizedSystem:
    xml = folder / f"{name}.xml"
    xml.write_text(openmm.XmlSerializer.serialize(system))
    realized = RealizedSystem.from_openmm(xml, folder / "build" / "solvated_system.pdb", 2.0)
    if positions is not None:
        realized.positions_nm = positions
    if box is not None:
        realized.box_nm = box
    return realized


def _copy(system):
    return openmm.XmlSerializer.deserialize(openmm.XmlSerializer.serialize(system))


def _run(built, system=None, config=None, positions=None, box=None, name="case"):
    builder, base_config, folder, base_system = built
    realized = _realize(folder, system or base_system, name, positions, box)
    return {
        (f.check, f.component): f
        for f in run_checks(realized, builder.components(), SYSTEM_CHECKS, config or base_config)
    }


def _status(findings, check):
    return {f.status for (c, _), f in findings.items() if c == check}


def _atoms(builder, type_):
    return next(c for c in builder.components() if c.type == type_).atom_indices()


def test_a_correct_build_passes_every_system_check(built) -> None:
    findings = _run(built, name="correct")
    bad = [f for f in findings.values() if f.status not in ("PASS", "INFO")]
    assert not bad, bad
    assert {c for c, _ in findings} >= set(SYSTEM_CHECKS) - {"ions.ionic_strength"}
    assert "trpcage" in findings[("restraint.atom_identity", "system")].detail


def test_a_config_that_disagrees_with_the_build_fails_declared(built) -> None:
    """The records say what was built; a config naming another model is caught."""
    config = built[1].model_copy(deep=True)
    substrate = next(s for s in config.solutes if s.type == "small_molecule")
    object.__setattr__(substrate.charges, "model", "openff-gnn-am1bcc-0.1.0-rc.3.pt")
    assert _status(_run(built, config=config), "components.declared") == {"FAIL"}


def test_a_chain_count_other_than_the_config_fails(built) -> None:
    config = built[1].model_copy(deep=True)
    object.__setattr__(config.polymers, "count", config.polymers.count + 1)
    assert _status(_run(built, config=config), "composition.counts") == {"FAIL"}


def test_salt_and_ionic_strength_are_checked_against_the_config(built) -> None:
    config = built[1].model_copy(deep=True)
    object.__setattr__(config.solvent.ions, "nacl_concentration", 0.5)
    assert _status(_run(built, config=config), "ions.salt") == {"FAIL"}
    config = built[1].model_copy(deep=True)
    object.__setattr__(config.solvent.ions, "nacl_concentration", 0.0)
    object.__setattr__(config.solvent.ions, "ionic_strength", 0.5)
    findings = _run(built, config=config)
    assert _status(findings, "ions.ionic_strength") == {"FAIL"}
    assert ("ions.salt", "system") not in findings


def test_a_net_charge_fails_neutral(built) -> None:
    system = _copy(built[3])
    force = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    q, sigma, epsilon = force.getParticleParameters(0)
    force.setParticleParameters(0, q + 0.5 * q.unit, sigma, epsilon)
    assert _status(_run(built, system=system), "ions.neutral") == {"FAIL"}


def test_a_hydrogen_without_its_constraint_fails_hbonds(built) -> None:
    system = _copy(built[3])
    system.removeConstraint(0)
    assert _status(_run(built, system=system), "constraints.hbonds") == {"FAIL"}


def test_a_cutoff_other_than_the_force_fields_fails(built) -> None:
    system = _copy(built[3])
    force = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    force.setCutoffDistance(1.2)
    assert _status(_run(built, system=system), "nonbonded.settings") == {"FAIL"}


def test_a_water_placed_inside_the_protein_fails_and_is_reported_buried(built) -> None:
    builder = built[0]
    realized = _realize(built[2], built[3], "positions")
    positions = np.array(realized.positions_nm)
    protein = _atoms(builder, "protein")
    solvent = next(c for c in builder.components() if c.type == "solvent")
    water = next(d for d in solvent.domains if d.name == "water").atoms[0]
    centre = positions[protein].mean(axis=0)
    positions[list(water)] += centre - positions[water.start]
    findings = _run(built, positions=positions)
    assert _status(findings, "geometry.water_overlap") == {"FAIL"}


def test_a_polymer_atom_on_the_protein_fails_polymer_overlap(built) -> None:
    builder = built[0]
    positions = np.array(_realize(built[2], built[3], "p2").positions_nm)
    polymer = next(c for c in builder.components() if c.type == "polymer")
    chain = list(polymer.molecules[0])
    protein = _atoms(builder, "protein")
    positions[chain] += positions[protein[10]] - positions[chain[0]] + 0.05
    assert _status(_run(built, positions=positions), "geometry.polymer_overlap") == {"FAIL"}


def test_two_chains_on_top_of_each_other_warn(built) -> None:
    builder = built[0]
    positions = np.array(_realize(built[2], built[3], "p3").positions_nm)
    polymer = next(c for c in builder.components() if c.type == "polymer")
    first, second = list(polymer.molecules[0]), list(polymer.molecules[1])
    positions[second] += positions[first[0]] - positions[second[0]] + 0.05
    assert _status(_run(built, positions=positions), "geometry.interchain_overlap") == {"WARN"}


def test_a_box_too_small_for_its_padding_fails_periodic_image(built) -> None:
    realized = _realize(built[2], built[3], "p4")
    box = np.array(realized.box_nm) * 0.55
    assert _status(_run(built, box=box), "geometry.periodic_image") == {"FAIL"}


def test_a_wrong_restraint_constant_or_atom_fails(built) -> None:
    builder = built[0]
    system = _copy(built[3])
    force = next(f for f in system.getForces() if f.getName() == "polyzymd_restraint:tether")
    a1, a2, (r0, k) = force.getBondParameters(0)
    force.setBondParameters(0, a1, a2, [r0, k * 2])
    assert _status(_run(built, system=system), "restraint.parameters") == {"FAIL"}

    system = _copy(built[3])
    force = next(f for f in system.getForces() if f.getName() == "polyzymd_restraint:tether")
    force.setBondParameters(0, a1 + 1, a2, [r0, k])  # the atom after the selected one
    findings = _run(built, system=system)
    assert _status(findings, "restraint.atom_identity") == {"FAIL"}
    assert "selects 1 atom" in findings[("restraint.atom_identity", "system")].detail


def test_two_restraints_with_different_distances_can_run_together() -> None:
    """Restraint parameters are per bond: global ones collided in one Context."""
    from polyzymd.core.restraints import RestraintFactory

    system = openmm.System()
    for _ in range(4):
        system.addParticle(12.0)
    for k, (name, distance) in enumerate((("a", 3.0), ("b", 10.0))):
        restraint = RestraintFactory.from_config(
            {"type": "flat_bottom", "name": name, "distance": distance, "atom1": {}, "atom2": {}}
        )
        force = restraint._create_flat_bottom_force(2 * k, 2 * k + 1)
        system.addForce(force)
    openmm.Context(
        system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference")
    )


def test_a_chain_touching_the_protein_and_its_image_fails_bridging() -> None:
    """A hand-made box: a straight chain reaching from the protein to the protein's image."""
    from openmm import app

    from polyzymd.audit.findings import AuditContext
    from polyzymd.audit.checks.geometry import image_bridging
    from polyzymd.core.components import AtomRange, ChargeSource, ComponentRecord, Domain

    n_chain = 7
    system = openmm.System()
    topology = app.Topology()
    chain = topology.addChain()
    bonds = openmm.HarmonicBondForce()
    nonbonded = openmm.NonbondedForce()
    for i in range(1 + n_chain):
        system.addParticle(12.0)
        nonbonded.addParticle(0.0, 0.3, 0.1)
        topology.addAtom("C", app.element.carbon, topology.addResidue("X", chain))
    for i in range(1, n_chain):
        bonds.addBond(i, i + 1, 0.4, 1000.0)
    system.addForce(bonds)
    system.addForce(nonbonded)
    positions = np.zeros((1 + n_chain, 3)) + 1.5
    positions[0, 0] = 0.3  # the protein
    positions[1:, 0] = 0.6 + 0.4 * np.arange(n_chain)  # 0.6 ... 3.0: touches 0.3 and 3.3
    realized = RealizedSystem("openmm", system, topology, positions, np.eye(3) * 3.0)

    def record(name, type_, segid, ranges):
        return ComponentRecord(
            name,
            type_,
            segid,
            "fixed",
            ranges,
            [Domain(name, tuple(ranges), ("ff",), ChargeSource("library"))],
        )

    components = [
        record("protein", "protein", "PROT", [AtomRange(0, 1)]),
        record("polymer", "polymer", "POLY", [AtomRange(1, 1 + n_chain)]),
    ]
    findings = image_bridging(AuditContext(realized, components))
    assert [f.status for f in findings] == ["FAIL"]
