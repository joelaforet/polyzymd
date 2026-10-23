"""The parameter checks pass a correct build and fail each way a build has gone wrong.

One small system (trp-cage, a substrate, polymers, DMSO, water, NaCl) is built
once; each test breaks a copy of its OpenMM system the way a real build went
wrong and checks that the right check, on the right component, fails.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

pytest.importorskip("openff.interchange")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import openmm  # noqa: E402
from _support.small_system import build  # noqa: E402

from polyzymd.audit.findings import run_checks  # noqa: E402
from polyzymd.audit.realized import RealizedSystem  # noqa: E402

PARAMETER_CHECKS = (
    "parameters.bonded",
    "parameters.nonbonded",
    "parameters.constraints",
    "charges.source",
    "masses",
)


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    folder = tmp_path_factory.mktemp("parameter_checks")
    builder, config = build(folder, replicate=2)
    builder.create_interchange()
    return builder, config, folder


def _realize(folder: Path, system, name: str) -> RealizedSystem:
    xml = folder / f"{name}.xml"
    xml.write_text(openmm.XmlSerializer.serialize(system))
    return RealizedSystem.from_openmm(xml, folder / "build" / "solvated_system.pdb")


def _openmm_system(builder):
    return builder._interchange.to_openmm(combine_nonbonded_forces=True)


def _status(findings) -> dict[tuple[str, str], str]:
    return {(f.component, f.check): f.status for f in findings}


def _failed(findings) -> set[tuple[str, str]]:
    return {(f.component, f.check) for f in findings if f.status == "FAIL"}


@pytest.fixture(scope="module")
def correct(built):
    builder, _, folder = built
    realized = _realize(folder, _openmm_system(builder), "correct")
    return realized, run_checks(realized, builder.components(), PARAMETER_CHECKS)


def test_a_correct_build_passes_every_parameter_check(built, correct) -> None:
    builder, _, _ = built
    _, findings = correct
    names = [c.name for c in builder.components()]
    assert {(f.component, f.check) for f in findings} == {
        (name, check) for name in names for check in PARAMETER_CHECKS
    }
    assert all(f.status == "PASS" for f in findings), [f for f in findings if f.status != "PASS"]


def test_the_gromacs_files_pass_the_same_checks(built) -> None:
    """A GROMACS build is checked on its .top, .gro and production .mdp."""
    if shutil.which("gmx") is None and not Path.home().joinpath(".pixi/bin/gmx").exists():
        pytest.skip("GROMACS is not installed")
    from polyzymd.exporters.gromacs import GromacsExporter

    builder, config, folder = built
    files = GromacsExporter(builder._interchange, config, builder.get_component_info()).export(
        folder / "gromacs"
    )
    realized = RealizedSystem.from_gromacs(files["top"], files["gro"], files["prod_mdp"])
    findings = run_checks(realized, builder.components(), PARAMETER_CHECKS)
    assert all(f.status == "PASS" for f in findings), [f for f in findings if f.status != "PASS"]


def test_the_force_field_order_bug_fails_the_protein(built, monkeypatch) -> None:
    """The 1.3 bug: Sage loaded after ff14SB overrode the protein's bonded and LJ terms."""
    from polyzymd.builders.system_builder import SystemBuilder

    builder, _, folder = built
    monkeypatch.setattr(SystemBuilder, "_apply_protein_force_field", lambda self, ic: None)
    buggy = builder._create_interchange_single_call(
        __import__("openff.toolkit").toolkit.ForceField(builder._protein_ff, builder._sm_ff),
        __import__("polyzymd.data.solvent_molecules").data.solvent_molecules.get_solvent_molecule(
            "tip3p"
        ),
    )
    monkeypatch.undo()
    realized = _realize(folder, buggy.to_openmm(combine_nonbonded_forces=True), "order_bug")
    failed = _failed(run_checks(realized, builder.components(), PARAMETER_CHECKS))
    protein = builder.components()[0].name
    assert {(protein, "parameters.bonded"), (protein, "parameters.nonbonded")} <= failed
    assert {component for component, _ in failed} == {protein}


def test_ff14sb_loaded_last_fails_a_small_molecule(built) -> None:
    """Loading ff14SB last lets its cap and amide SMIRKS match the substrate's amide."""
    from openff.toolkit import ForceField

    from polyzymd.data.solvent_molecules import get_solvent_molecule

    builder, _, folder = built
    reversed_ff = ForceField(builder._sm_ff, builder._protein_ff)
    wrong = builder._create_interchange_single_call(reversed_ff, get_solvent_molecule("tip3p"))
    realized = _realize(folder, wrong.to_openmm(combine_nonbonded_forces=True), "ff14sb_last")
    failed = _failed(run_checks(realized, builder.components(), PARAMETER_CHECKS))
    names = {c.type: c.name for c in builder.components()}
    assert any(
        component in (names["small_molecule"], names["polymer"], "solvent")
        for component, check in failed
        if check.startswith("parameters.")
    ), failed


def _broken(built, name, edit) -> list:
    builder, _, folder = built
    system = _openmm_system(builder)
    edit(system, builder)
    return run_checks(_realize(folder, system, name), builder.components(), PARAMETER_CHECKS)


def _nonbonded(system):
    return next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))


def _atoms_of(builder, type_):
    return next(c for c in builder.components() if c.type == type_).atom_indices()


def test_a_wrong_charge_fails_charges_source(built) -> None:
    """One substrate atom charged 0.05 e off its NAGL charge fails that component."""

    def edit(system, builder):
        atom = _atoms_of(builder, "small_molecule")[3]
        force = _nonbonded(system)
        q, sigma, epsilon = force.getParticleParameters(atom)
        force.setParticleParameters(atom, q + 0.05 * q.unit, sigma, epsilon)

    findings = _broken(built, "wrong_charge", edit)
    builder = built[0]
    substrate = next(c.name for c in builder.components() if c.type == "small_molecule")
    status = _status(findings)
    assert status[(substrate, "charges.source")] == "FAIL"
    assert status[(substrate, "parameters.bonded")] == "PASS"


def test_a_missing_constraint_fails_parameters_constraints(built) -> None:
    """A protein bond to hydrogen left unconstrained at a 2 fs step fails."""

    def edit(system, builder):
        protein = set(_atoms_of(builder, "protein"))
        for k in range(system.getNumConstraints()):
            a, b, _ = system.getConstraintParameters(k)
            if a in protein:
                system.removeConstraint(k)
                return

    status = _status(_broken(built, "missing_constraint", edit))
    protein = built[0].components()[0].name
    assert status[(protein, "parameters.constraints")] == "FAIL"
    assert status[(protein, "masses")] == "PASS"


def test_a_wrong_mass_fails_masses(built) -> None:
    """A hydrogen given 3 Da (as mass repartitioning does) fails the masses check."""

    def edit(system, builder):
        protein = _atoms_of(builder, "protein")
        elements = [a.element.symbol for a in builder.get_openmm_components()[0].atoms()]
        hydrogen = next(i for i in protein if elements[i] == "H")
        system.setParticleMass(hydrogen, 3.024)

    status = _status(_broken(built, "wrong_mass", edit))
    protein = built[0].components()[0].name
    assert status[(protein, "masses")] == "FAIL"
    assert status[(protein, "charges.source")] == "PASS"


def test_parameters_a_change_in_openff_would_move_are_pinned(correct, built) -> None:
    """ff14SB's backbone C-N bond and Sage's C-H constraint length, as the build realizes them.

    If OpenFF Interchange ever applies these force fields differently, this fails in CI
    before a build does.
    """
    realized, _ = correct
    builder = built[0]
    atoms = list(realized.topology.atoms())
    protein = _atoms_of(builder, "protein")
    by_residue = {}
    for i in protein:
        by_residue[(atoms[i].residue.index, atoms[i].name)] = i
    first = atoms[protein[0]].residue.index
    carbon, nitrogen = by_residue[(first, "C")], by_residue[(first + 1, "N")]
    bonds = realized.index.extract(protein)["bonds"]
    r0, _ = bonds[tuple(sorted((carbon, nitrogen)))]
    assert r0 == pytest.approx(0.1335, abs=1e-6)  # ff14SB C-N peptide bond, 1.335 Å
    substrate = _atoms_of(builder, "small_molecule")
    elements = realized.elements
    constraints = realized.index.extract(substrate)["constraints"]
    lengths = {
        round(d, 6) for (a, b), d in constraints.items() if {elements[a], elements[b]} == {"C", "H"}
    }
    assert lengths == SAGE_CH_CONSTRAINTS_NM


#: The C-H constraint lengths (aromatic and aliphatic) Sage 2.3.0 gives the substrate,
#: as OpenFF Interchange realized them when this test was written.
SAGE_CH_CONSTRAINTS_NM = {0.108664, 0.109245}
