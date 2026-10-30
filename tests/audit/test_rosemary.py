"""The audit of a protein on OpenFF Rosemary alpha: NAGL charges on the whole molecule."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from polyzymd.audit.findings import run_checks
from polyzymd.audit.realized import RealizedSystem
from polyzymd.core.components import AtomRange, ChargeSource, ComponentRecord, Domain

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _support.rosemary import newest_rosemary  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
TRPCAGE = ROOT / "examples" / "quickstart" / "trpcage.pdb"
ROSEMARY = newest_rosemary()
MODEL = "openff-gnn-am1bcc-1.0.0.pt"
CHECKS = [
    "charges.source",
    "parameters.bonded",
    "parameters.nonbonded",
    "parameters.one_four",
    "masses",
    "parameters.rosemary_alpha",
]

pytestmark = pytest.mark.filterwarnings("ignore")


def _realized(folder: Path, library_charges: bool, model: str = MODEL) -> RealizedSystem:
    """Trp-cage on Rosemary alpha in vacuum; optionally with ff14SB's library charges."""
    import openmm
    from openff.toolkit import ForceField, Molecule, Topology
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

    topology = Topology.from_pdb(str(TRPCAGE))
    molecule = Molecule(topology.molecule(0))
    NAGLToolkitWrapper().assign_partial_charges(molecule, partial_charge_method=model)
    system = (
        ForceField(ROSEMARY)
        .create_interchange(topology, charge_from_molecules=[molecule])
        .to_openmm(combine_nonbonded_forces=True)
    )
    if library_charges:
        amber = ForceField("ff14sb_off_impropers_0.0.4.offxml").create_interchange(topology)
        library = amber["Electrostatics"].get_charge_array()
        force = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
        for i, charge in enumerate(library):
            _, sigma, epsilon = force.getParticleParameters(i)
            force.setParticleParameters(i, charge.m, sigma, epsilon)
    (folder / "system.xml").write_text(openmm.XmlSerializer.serialize(system))
    return RealizedSystem.from_openmm(folder / "system.xml", TRPCAGE)


def _component(
    n_atoms: int, model: str = MODEL, charges: ChargeSource | None = None, force_field=ROSEMARY
) -> ComponentRecord:
    atoms = (AtomRange(0, n_atoms),)
    return ComponentRecord(
        name="trpcage",
        type="protein",
        segid="PROT",
        placement="fixed",
        molecules=list(atoms),
        domains=[
            Domain(
                "protein",
                atoms,
                (force_field,),
                charges or ChargeSource("nagl", model=model),
                reference={"pdb": str(TRPCAGE), "custom_substructures": None},
            )
        ],
    )


def _status(findings) -> dict[str, str]:
    return {finding.check: finding.status for finding in findings}


def test_a_rosemary_protein_with_nagl_charges_passes_and_is_flagged_as_alpha(tmp_path) -> None:
    realized = _realized(tmp_path, library_charges=False)
    status = _status(run_checks(realized, [_component(realized.n_atoms)], CHECKS))
    assert status.pop("parameters.rosemary_alpha") == "WARN"
    assert set(status.values()) == {"PASS"}, status


def test_library_charges_declared_as_rosemary_nagl_fail_charges_source(tmp_path) -> None:
    realized = _realized(tmp_path, library_charges=True)
    status = _status(run_checks(realized, [_component(realized.n_atoms)], CHECKS))
    assert status["charges.source"] == "FAIL"
    assert status["parameters.bonded"] == "PASS"


def test_the_declared_nagl_model_is_the_one_the_audit_recomputes(tmp_path) -> None:
    other = "openff-gnn-am1bcc-0.1.0-rc.1.pt"
    realized = _realized(tmp_path, library_charges=False, model=other)
    declared_other = _status(run_checks(realized, [_component(realized.n_atoms, other)], CHECKS))
    declared_default = _status(run_checks(realized, [_component(realized.n_atoms)], CHECKS))
    assert declared_other["charges.source"] == "PASS"
    assert declared_default["charges.source"] == "FAIL"


def test_library_charges_declared_on_a_nagl_force_field_fail(tmp_path) -> None:
    """Rosemary would charge the reference itself; declaring library must not pass."""
    realized = _realized(tmp_path, library_charges=False)
    component = _component(realized.n_atoms, charges=ChargeSource("library"))
    findings = run_checks(realized, [component], ["charges.source"])
    (finding,) = findings
    assert finding.status == "FAIL" and "no protein library charges" in finding.detail


def test_an_ff14sb_protein_gets_no_rosemary_finding(tmp_path) -> None:
    realized = _realized(tmp_path, library_charges=False)
    component = _component(
        realized.n_atoms,
        charges=ChargeSource("library"),
        force_field="ff14sb_off_impropers_0.0.4.offxml",
    )
    assert run_checks(realized, [component], ["parameters.rosemary_alpha"]) == []


@pytest.mark.parametrize(
    ("verb", "site", "atom"), [("deprotonate", "A:N1", "N"), ("protonate", "A:S20", "OXT")]
)
def test_the_audit_passes_a_neutral_terminus_on_rosemary(tmp_path, verb, site, atom) -> None:
    """Rosemary needs no residue templates, so a terminus changed by protonation
    alone passes the audit of a protein parameterized whole."""
    import openmm

    from polyzymd import modifying_proteins as mp
    from polyzymd.modifying_proteins.components import protein_component

    protein = mp.load_protein(TRPCAGE)
    getattr(mp, verb)(protein, site, atom)
    result = mp.parameterize_whole(protein, (ROSEMARY,), {"method": "nagl", "model": MODEL})
    (tmp_path / "system.xml").write_text(openmm.XmlSerializer.serialize(result.to_openmm()))
    protein.save_pdb(tmp_path / "protein.pdb")
    realized = RealizedSystem.from_openmm(tmp_path / "system.xml", tmp_path / "protein.pdb")
    component = protein_component(protein, result.record, {}, name="trpcage")
    (finding,) = run_checks(realized, [component], ["modifications.sites"])
    assert finding.status == "PASS", finding.detail
