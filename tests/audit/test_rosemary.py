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


SAGE = "openff-2.3.0.offxml"


def _solvated(water: str, co_solvent: str | None = None) -> list[ComponentRecord]:
    """A Rosemary protein of 304 atoms, one water of ``water`` and optionally a co-solvent."""
    water_atoms = (AtomRange(304, 307),)
    domains = [
        Domain(
            "water",
            water_atoms,
            (SAGE, "opc3.offxml") if water == "opc3" else (SAGE,),
            ChargeSource("water_model", water),
            reference={"smiles": "O"},
        )
    ]
    molecules = list(water_atoms)
    if co_solvent is not None:
        atoms = (AtomRange(307, 317),)
        molecules += atoms
        domains.append(
            Domain(
                "co_solvent:dmso",
                atoms,
                (co_solvent,),
                ChargeSource("preset", source_file="dmso.sdf"),
                reference={"smiles": "CS(=O)C"},
            )
        )
    solvent = ComponentRecord(
        name="solvent",
        type="solvent",
        segid="SOL",
        placement="packed",
        molecules=molecules,
        domains=domains,
    )
    return [_component(304), solvent]


def _audit(components, exceptions=()):
    from types import SimpleNamespace

    from polyzymd.audit.report import audit
    from polyzymd.config.schema import ModelingException

    config = SimpleNamespace(
        modeling_exceptions=[ModelingException(**entry) for entry in exceptions]
    )
    realized = SimpleNamespace(n_atoms=317)
    ids = ["solvent.rosemary_water", "parameters.rosemary_components"]
    return {f.check: f for f in audit(realized, components, config, ids).findings}


def test_a_rosemary_protein_in_tip3p_fails_rosemary_water() -> None:
    finding = _audit(_solvated("tip3p"))["solvent.rosemary_water"]
    assert finding.status == "FAIL"
    assert "OPC3" in finding.detail and "solvent.water.model: opc3" in finding.detail


def test_a_modeling_exception_accepts_rosemary_in_tip3p() -> None:
    reason = "Comparing with earlier TIP3P runs"
    findings = _audit(_solvated("tip3p"), [{"check": "solvent.rosemary_water", "reason": reason}])
    finding = findings["solvent.rosemary_water"]
    assert (finding.status, finding.accepted_from, finding.reason) == ("ACCEPTED", "FAIL", reason)
    assert "modeling_exceptions.unused" not in findings


def test_a_rosemary_protein_in_opc3_passes_rosemary_water() -> None:
    assert _audit(_solvated("opc3"))["solvent.rosemary_water"].status == "PASS"


def test_an_ff14sb_protein_in_tip3p_gets_no_rosemary_water_finding() -> None:
    components = _solvated("tip3p")
    components[0] = _component(
        304, charges=ChargeSource("library"), force_field="ff14sb_off_impropers_0.0.4.offxml"
    )
    assert _audit(components) == {}


def test_soluble_components_off_rosemary_get_an_info_recommendation() -> None:
    finding = _audit(_solvated("opc3", co_solvent=SAGE))["parameters.rosemary_components"]
    assert finding.status == "INFO"
    assert "co-solvent dmso (openff-2.3.0.offxml) is not" in finding.detail
    assert "recommended for them too" in finding.detail and ROSEMARY in finding.detail
    on_rosemary = _audit(_solvated("opc3", co_solvent=ROSEMARY))
    assert on_rosemary["parameters.rosemary_components"].status == "PASS"


def test_glycans_on_a_rosemary_protein_are_an_info_line(monkeypatch) -> None:
    from types import SimpleNamespace

    from polyzymd.audit.checks import modifications

    reference = SimpleNamespace(
        recipe={"protein_charges": {"method": "nagl", "model": MODEL}},
        site_types={},
        records=[
            SimpleNamespace(type="glycan", site="A:N5"),
            SimpleNamespace(type=None, site="A:K8"),
        ],
    )
    monkeypatch.setattr(modifications, "is_modified", lambda component: True)
    monkeypatch.setattr(modifications, "modified_reference", lambda context, c: reference)
    (finding,) = run_checks(
        SimpleNamespace(n_atoms=304), [_component(304)], ["parameters.glycan_force_field"]
    )
    assert finding.status == "INFO"
    assert finding.detail == "glycans parameterized by Rosemary (whole-protein mode): 1 at A:N5"
