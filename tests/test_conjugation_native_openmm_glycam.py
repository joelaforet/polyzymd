"""Tests for the native OpenMM GLYCAM handoff seam."""

from __future__ import annotations

from types import SimpleNamespace

import pytest
from pydantic import ValidationError

from polyzymd.builders.conjugation.force_fields import resolve_conjugate_force_fields
from polyzymd.builders.conjugation.native_openmm_glycam import (
    _annotate_force_field_domains_from_config,
    _build_native_glycam_audit,
    _domain_assignment_audit,
    _openff_topology_to_openmm_for_glycam,
    _register_disconnected_sage_template_generator,
    _require_essential_linkage_terms,
    create_native_openmm_glycam_handoff,
    native_glycam_enabled,
)
from polyzymd.config.schema import (
    ConjugationConfig,
    EquilibrationStageConfig,
    ForceFieldConfig,
    SimulationConfig,
    SimulationPhaseConfig,
    SimulationPhasesConfig,
    SolventConfig,
    ThermodynamicsConfig,
)
from polyzymd.exporters.exact_openmm import ExactExportBundle
from tests.test_exact_openmm_export import _sidecar


def test_force_field_native_glycam_is_attachment_scoped() -> None:
    """Keep the OpenFF route as default unless an attachment requests GLYCAM."""
    config = _native_config(
        force_field=ForceFieldConfig(),
        conjugation=_conjugation_with_domain("openff-2.0.0.offxml"),
    )
    resolved = resolve_conjugate_force_fields(config)
    assert resolved.route == "standard_interchange"
    config = _native_config(conjugation=_conjugation_with_domain("glycam06"))
    resolved = resolve_conjugate_force_fields(config)
    assert resolved.route == "native_exact"
    assert resolved.attachments[0].source == "glycam06"


def test_force_field_rejects_unknown_moiety_label() -> None:
    """Unknown moiety labels must never silently fall back to Sage."""
    with pytest.raises(ValidationError, match="Unknown moiety.force_field"):
        _native_config(conjugation=_conjugation_with_domain("glycam-special"))

    with pytest.raises(ValidationError, match="explicitly declare moiety.force_field"):
        _conjugation_with_missing_owner()


def test_simulation_config_rejects_removed_routing_fields() -> None:
    """Removed mixed-force-field UX fields are forbidden by the schema."""
    with pytest.raises(ValueError, match="force_field_domain"):
        _conjugation_with_removed_domain_field()
    with pytest.raises(ValueError, match="glycan_policy"):
        ForceFieldConfig.model_validate({"glycan_policy": "strict_glycam"})


def test_native_glycam_enabled_reads_force_field_mode() -> None:
    """Detect only the attachment-scoped native GLYCAM-only route."""
    assert not native_glycam_enabled(
        _native_config(conjugation=_conjugation_with_domain("openff-2.0.0.offxml"))
    )
    assert native_glycam_enabled(_native_config(conjugation=_conjugation_with_domain("glycam06")))
    assert not native_glycam_enabled(_native_config(conjugation=_mixed_conjugation()))


def test_exact_bundle_exposes_authoritative_openmm_methods() -> None:
    """Return stored OpenMM objects without pretending to be Interchange."""
    bundle = ExactExportBundle(
        topology=object(),
        system=object(),
        positions=object(),
        private_baseline_interchange=object(),
        sidecar=_sidecar(),
    )

    assert bundle.to_openmm() is bundle.system
    assert bundle.to_openmm_topology() is bundle.topology
    assert bundle.to_openmm_positions() is bundle.positions
    assert bundle.is_exact_export_bundle is True


def test_create_native_handoff_maps_asx_to_nln_and_writes_audit(monkeypatch, tmp_path) -> None:
    """Create native handoff through the narrow seam with NLN residue mapping."""
    asx_residue = object()
    converted = SimpleNamespace(
        topology=SimpleNamespace(),
        positions=[object()],
        asx_residue=asx_residue,
        crosslink_atoms=(object(), object()),
        renamed_atoms=({"atom_name": "HD22", "native_atom_name": "HD21"},),
    )
    calls = {}

    class FakeForceField:
        def createSystem(self, topology, **kwargs):  # noqa: N802 - OpenMM API name
            calls["topology"] = topology
            calls["kwargs"] = kwargs
            return SimpleNamespace()

    monkeypatch.setattr(
        "polyzymd.builders.conjugation.native_openmm_glycam._openff_topology_to_openmm_for_glycam",
        lambda topology, construction, **kwargs: converted,
    )
    monkeypatch.setattr(
        "polyzymd.builders.conjugation.native_openmm_glycam._load_native_glycam_force_field",
        lambda: FakeForceField(),
    )
    monkeypatch.setattr(
        "polyzymd.builders.conjugation.native_openmm_glycam._build_native_glycam_audit",
        lambda *args, **kwargs: {
            "adjacent_terms": {
                "bonds": [{"category": "crosslink_bond"}],
                "angles": [
                    {"category": "asx_side_angle"},
                    {"category": "glycan_side_angle"},
                ],
                "torsions": [{"category": "proper_crosslink_torsion"}],
            },
            "local_exclusions_and_14_exceptions": [
                {"category": "zero_12_exclusion"},
                {"category": "zero_13_exclusion"},
                {"category": "scaled_14_exception"},
            ],
        },
    )
    monkeypatch.setattr(
        "polyzymd.builders.conjugation.native_openmm_glycam.create_exact_export_bundle",
        lambda **kwargs: ExactExportBundle(
            topology=kwargs["topology"],
            system=kwargs["system"],
            positions=kwargs["positions"],
            private_baseline_interchange=object(),
            sidecar=_sidecar(),
            audit_path=kwargs["audit_path"],
            audit=kwargs["audit"],
        ),
    )
    builder = SimpleNamespace(_solvated_topology=SimpleNamespace())

    handoff = create_native_openmm_glycam_handoff(
        builder,
        config=_native_config(),
        construction=SimpleNamespace(),
        output_dir=tmp_path,
    )

    assert calls["kwargs"]["residueTemplates"] == {asx_residue: "NLN"}
    assert handoff.audit_path.exists()
    assert builder._exact_export_bundle is handoff


def test_native_conversion_renames_source_hd22_to_native_hd21_without_mutating_source() -> None:
    """Keep Pablo/OpenFF source metadata unchanged while native topology uses NLN HD21."""
    topology = _fake_native_source_topology(asx_hydrogens=("HD22",))

    converted = _openff_topology_to_openmm_for_glycam(topology, construction=SimpleNamespace())

    source_names = [atom.metadata["atom_name"] for atom in topology.molecules[0].atoms]
    native_names = [atom.name for atom in converted.asx_residue.atoms()]
    assert "HD22" in source_names
    assert "HD21" not in source_names
    assert "HD21" in native_names
    assert "HD22" not in native_names


@pytest.mark.parametrize("asx_hydrogens", [(), ("HD21",), ("HD21", "HD22")])
def test_native_conversion_rejects_missing_or_ambiguous_asx_hydrogens(asx_hydrogens) -> None:
    """Require exactly one source HD22 for native-only HD21 remapping."""
    topology = _fake_native_source_topology(asx_hydrogens=asx_hydrogens)

    with pytest.raises(ValueError, match="exactly one retained amide hydrogen"):
        _openff_topology_to_openmm_for_glycam(topology, construction=SimpleNamespace())


def test_native_audit_validation_identifies_missing_categories() -> None:
    """Report the exact missing linkage category from audit validation."""
    audit = {
        "adjacent_terms": {
            "bonds": [{"category": "crosslink_bond"}],
            "angles": [{"category": "asx_side_angle"}],
            "torsions": [{"category": "proper_crosslink_torsion"}],
        },
        "local_exclusions_and_14_exceptions": [
            {"category": "zero_12_exclusion"},
            {"category": "zero_13_exclusion"},
            {"category": "scaled_14_exception"},
        ],
    }

    with pytest.raises(RuntimeError, match="glycan_side_angle"):
        _require_essential_linkage_terms(audit)


def test_disconnected_precharged_sage_molecule_registers_after_glycan_exclusion(
    monkeypatch,
) -> None:
    """Only complete disconnected Sage molecules enter SMIRNOFFTemplateGenerator."""
    calls = {}

    class FakeGenerator:
        def __init__(self, *, molecules, forcefield):
            calls["molecules"] = molecules
            calls["forcefield"] = forcefield
            self.generator = object()

    class FakeForceField:
        def registerTemplateGenerator(self, generator):  # noqa: N802 - OpenMM API name
            calls["registered"] = generator

    monkeypatch.setitem(
        __import__("sys").modules,
        "openmmforcefields.generators",
        SimpleNamespace(SMIRNOFFTemplateGenerator=FakeGenerator),
    )
    glycan = _FakeSageMolecule("glycan", domain="glycan", charges=[0.0])
    dmso = _FakeSageMolecule("dmso", domain="sage", charges=[-0.2, 0.2])

    audit = _register_disconnected_sage_template_generator(
        FakeForceField(),
        SimpleNamespace(molecules=[glycan, dmso]),
        _native_config(),
        SimpleNamespace(),
    )

    assert calls["molecules"] == [dmso]
    assert calls["forcefield"] == "openff-2.0.0.offxml"
    assert any(entry["domain"] == "glycan" and not entry["routed_to_sage"] for entry in audit)


def test_native_route_classifies_dmso_cosolvent_as_disconnected_sage() -> None:
    """DMSO co-solvent residue names should be Sage before GLYCAM template matching."""
    dmso = SimpleNamespace(
        atoms=[
            _FakeDomainAtom("DMS", "C1"),
            _FakeDomainAtom("DMS", "S1"),
            _FakeDomainAtom("DMS", "O1"),
        ]
    )
    glycan = SimpleNamespace(atoms=[_FakeDomainAtom("4YB", "C1")])
    topology = SimpleNamespace(molecules=[glycan, dmso])
    config = _native_config(
        solvent=SolventConfig.model_validate(
            {
                "primary": {"type": "water", "model": "tip3p"},
                "co_solvents": [{"name": "dmso", "mole_fraction": 0.01, "residue_name": "DMS"}],
            }
        )
    )

    _annotate_force_field_domains_from_config(topology, SimpleNamespace(), config)

    assert not hasattr(glycan, "force_field_domain")
    assert dmso.force_field_domain == "sage"
    assert {atom.metadata["force_field_domain"] for atom in dmso.atoms} == {"sage"}


def test_sage_molecule_missing_charges_or_covalent_boundary_rejected() -> None:
    """Missing Sage charges and covalent Amber/GLYCAM boundaries fail clearly."""
    with pytest.raises(ValueError, match="missing trusted assigned partial charges"):
        _register_disconnected_sage_template_generator(
            SimpleNamespace(),
            SimpleNamespace(molecules=[_FakeSageMolecule("dmso", domain="sage", charges=None)]),
            _native_config(),
            SimpleNamespace(),
        )
    with pytest.raises(ValueError, match="crosses an Amber/GLYCAM boundary"):
        _register_disconnected_sage_template_generator(
            SimpleNamespace(),
            SimpleNamespace(
                molecules=[
                    _FakeSageMolecule(
                        "polymer",
                        domain="sage",
                        charges=[0.0],
                        metadata={"cross_force_field_boundary": True},
                    )
                ]
            ),
            _native_config(),
            SimpleNamespace(),
        )


def test_multi_linkage_audit_validates_each_linkage() -> None:
    """Each NLN/glycan linkage must independently carry essential terms."""
    audit = {
        "linkages": (
            _complete_linkage_audit(0),
            {
                **_complete_linkage_audit(1),
                "adjacent_terms": {"bonds": [], "angles": [], "torsions": []},
            },
        )
    }

    with pytest.raises(RuntimeError, match="linkage 1 crosslink_bond"):
        _require_essential_linkage_terms(audit)


def test_audit_records_domains_linkages_sage_and_no_glycan_proof(monkeypatch) -> None:
    """Native audit exposes domain routing and unsupported boundary provenance."""
    monkeypatch.setattr(
        "polyzymd.builders.conjugation.native_openmm_glycam._nonbonded_force",
        lambda system: _FakeNonbonded(),
    )
    topology = _FakeAuditTopology()
    system = _FakeSystem()
    nd2, c1 = topology.crosslink

    audit = _build_native_glycam_audit(
        topology,
        system,
        (nd2.residue,),
        ((nd2, c1),),
        renamed_atoms=({"atom_name": "HD22", "native_atom_name": "HD21"},),
        template_matches=({"residue": "A:ASX60", "template": "NLN"},),
        sage_components=({"domain": "glycan", "routed_to_sage": False},),
        construction=SimpleNamespace(attachments=()),
    )

    assert audit["residue_templates"] == {"A:ASX60": "NLN"}
    assert audit["glycam_template_matches"]
    assert audit["sage_template_generator"]["proof_no_glycan_entered_sage"] is True
    assert (
        audit["unsupported_boundary_diagnostics"][
            "covalently_attached_sage_polymer_across_amber_glycam_boundary"
        ]
        == "unsupported"
    )


def test_disconnected_multi_residue_sage_polymer_converts_as_one_template_unit() -> None:
    """A complete disconnected Sage molecule should be one OpenMM residue for matching."""
    topology = _fake_native_source_topology_with_sage(_multi_residue_sage_molecule())

    converted = _openff_topology_to_openmm_for_glycam(topology, construction=SimpleNamespace())

    sage_residues = [residue for residue in converted.topology.residues() if residue.name == "SBM"]
    assert len(sage_residues) == 1
    assert len(list(sage_residues[0].atoms())) == 4
    assert sage_residues[0].force_field_domain == "sage"
    assert [segment["residue_name"] for segment in sage_residues[0].source_residue_segments] == [
        "SBM",
        "EGP",
    ]


def test_repeated_disconnected_sage_copies_get_deterministic_unique_residue_ids() -> None:
    """Repeated complete Sage copies should not reuse OpenMM residue identity."""
    topology = _fake_native_source_topology_with_sage(
        _multi_residue_sage_molecule(),
        _multi_residue_sage_molecule(start_residue_number=3),
    )

    converted = _openff_topology_to_openmm_for_glycam(topology, construction=SimpleNamespace())

    sage_residues = [residue for residue in converted.topology.residues() if residue.name == "SBM"]
    assert [residue.id for residue in sage_residues] == ["S1", "S2"]


def test_glycan_residue_names_are_not_collapsed_into_sage_template_units() -> None:
    """Strict glycan residues must remain GLYCAM-owned even when Sage is present."""
    topology = _fake_native_source_topology_with_sage(_multi_residue_sage_molecule())

    converted = _openff_topology_to_openmm_for_glycam(topology, construction=SimpleNamespace())

    glycan_residues = [
        residue for residue in converted.topology.residues() if residue.name == "4YB"
    ]
    sage_residues = [residue for residue in converted.topology.residues() if residue.name == "SBM"]
    assert len(glycan_residues) == 1
    assert len(sage_residues) == 1


def test_sage_source_residue_provenance_is_available_in_domain_audit() -> None:
    """Collapsed Sage template units should preserve original monomer segmentation."""
    topology = _fake_native_source_topology_with_sage(_multi_residue_sage_molecule())
    converted = _openff_topology_to_openmm_for_glycam(topology, construction=SimpleNamespace())

    audit = _domain_assignment_audit(list(converted.topology.residues()), ())

    sage_entries = [entry for entry in audit["residues"] if entry["domain"] == "sage"]
    assert len(sage_entries) == 1
    assert sage_entries[0]["source_component_id"] == "polymer-copy"
    assert [segment["residue_name"] for segment in sage_entries[0]["source_residue_segments"]] == [
        "SBM",
        "EGP",
    ]


def test_sage_molecule_with_explicit_glycam_boundary_bond_is_rejected() -> None:
    """A Sage molecule bonded to glycan/Amber atoms must fail before createSystem."""
    sage_atom = _FakeAtom("C1", "SBM", "1", 6)
    sage_atom.metadata["force_field_domain"] = "sage"
    glycan_atom = _FakeAtom("O1", "4YB", "1", 8)
    glycan_atom.metadata["force_field_domain"] = "glycan"
    molecule = _FakeMolecule([sage_atom, glycan_atom], [_FakeBond(sage_atom, glycan_atom)])
    molecule.force_field_domain = "sage"
    molecule.partial_charges = [0.0, 0.0]

    with pytest.raises(ValueError, match="crosses an Amber/GLYCAM boundary"):
        _openff_topology_to_openmm_for_glycam(
            _fake_native_source_topology_with_sage(molecule), construction=SimpleNamespace()
        )


def test_native_conversion_reuses_one_openmm_chain_per_chain_id_for_waters(tmp_path) -> None:
    """Native conversion should not create one OpenMM chain and TER per water."""

    from openmm.app import PDBFile

    base = _fake_native_source_topology(asx_hydrogens=("HD22",)).molecules[0]
    waters = [_water_molecule(index) for index in range(1, 4)]
    topology = _FakeMultiMoleculeTopology([base, *waters])

    converted = _openff_topology_to_openmm_for_glycam(topology, construction=SimpleNamespace())
    pdb_path = tmp_path / "native_waters.pdb"
    with pdb_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(converted.topology, converted.positions, handle, keepIds=True)
    reloaded = PDBFile(str(pdb_path)).topology

    chain_ids = [chain.id for chain in converted.topology.chains()]
    water_residues = [residue for residue in reloaded.residues() if residue.name == "HOH"]

    assert chain_ids.count("D") == 1
    assert len(water_residues) == 3
    assert all(len(list(residue.atoms())) == 3 for residue in water_residues)
    assert pdb_path.read_text(encoding="utf-8").count("TER") == len(chain_ids)


def _native_config(**updates):
    """Return a minimal native GLYCAM-compatible config double."""
    solvent = SolventConfig()
    solvent.ions.neutralize = False
    data = {
        "name": "native-glycam",
        "enzyme": {"name": "E", "pdb_path": "enzyme.pdb"},
        "conjugation": ConjugationConfig(enabled=False),
        "solvent": solvent,
        "thermodynamics": ThermodynamicsConfig(temperature=300.0),
        "simulation_phases": SimulationPhasesConfig(
            equilibration_stages=[
                EquilibrationStageConfig(
                    name="equilibration",
                    duration=0.1,
                    samples=1,
                    ensemble="NVT",
                    temperature=300.0,
                )
            ],
            production=SimulationPhaseConfig(
                ensemble="NVT",
                duration=0.1,
                samples=1,
                report_interval=1,
                checkpoint_interval=60.0,
            ),
        ),
        "force_field": ForceFieldConfig(),
        "engine": "openmm",
    }
    data.update(updates)
    return SimulationConfig.model_validate(data)


def _conjugation_with_domain(force_field: str | None) -> ConjugationConfig:
    """Return a minimal enabled conjugation config with one force-field-labeled moiety."""
    moiety = {
        "name": force_field or "generic",
        "input_path": "glycan.pdb",
        "link_site": {
            "chain_id": "C",
            "residue_name": "4YB",
            "residue_number": 1,
            "atom_name": "C1",
        },
    }
    if force_field is not None:
        moiety["force_field"] = force_field
    return ConjugationConfig.model_validate(
        {
            "enabled": True,
            "attachments": [
                {
                    "name": f"{force_field or 'generic'}-attachment",
                    "site": {
                        "chain_id": "A",
                        "residue_name": "ASN",
                        "residue_number": 60,
                        "atom_name": "ND2",
                    },
                    "moiety": moiety,
                    "mechanism": {
                        "name": "explicit_linkage",
                        "product_residues": {"site": "ASX", "moiety": "4YB"},
                    },
                }
            ],
        }
    )


def _conjugation_with_removed_domain_field() -> ConjugationConfig:
    """Return a config payload using a removed moiety routing field."""

    data = _conjugation_with_domain("glycam06").model_dump(mode="json")
    data["attachments"][0]["moiety"].pop("force_field")
    data["attachments"][0]["moiety"]["force_field_domain"] = "glycan"
    return ConjugationConfig.model_validate(data)


def _conjugation_with_missing_owner() -> ConjugationConfig:
    """Return an enabled attachment payload without force-field ownership."""

    data = _conjugation_with_domain("openff-2.0.0.offxml").model_dump(mode="json")
    data["attachments"][0]["moiety"].pop("force_field")
    return ConjugationConfig.model_validate(data)


def _mixed_conjugation() -> ConjugationConfig:
    """Return one GLYCAM and one generic attachment config."""

    data = _conjugation_with_domain("glycam06").model_dump(mode="json")
    generic = data["attachments"][0].copy()
    generic["name"] = "generic-attachment"
    generic["moiety"] = generic["moiety"].copy()
    generic["moiety"]["name"] = "generic"
    generic["moiety"]["force_field"] = "openff-2.0.0.offxml"
    data["attachments"].append(generic)
    return ConjugationConfig.model_validate(data)


class _FakeAtom:
    """Minimal OpenFF atom double for native GLYCAM conversion tests."""

    def __init__(self, atom_name: str, residue_name: str, residue_number: str, atomic_number: int):
        self.atomic_number = atomic_number
        self.name = atom_name
        self.metadata = {
            "chain_id": "A" if residue_name == "ASX" else "C",
            "residue_name": residue_name,
            "residue_number": residue_number,
            "atom_name": atom_name,
        }


class _FakeBond:
    """Minimal OpenFF bond double for native GLYCAM conversion tests."""

    def __init__(self, atom1: _FakeAtom, atom2: _FakeAtom):
        self.atom1 = atom1
        self.atom2 = atom2


class _FakeMolecule:
    """Minimal OpenFF molecule double for native GLYCAM conversion tests."""

    def __init__(self, atoms: list[_FakeAtom], bonds: list[_FakeBond]):
        self.atoms = atoms
        self.bonds = bonds


class _FakeTopology:
    """Minimal OpenFF topology double for native GLYCAM conversion tests."""

    box_vectors = None

    def __init__(self, molecule: _FakeMolecule):
        self.molecules = [molecule]

    def get_positions(self):
        """Return source positions in nanometers."""
        return [[float(index), 0.0, 0.0] for index, _atom in enumerate(self.molecules[0].atoms)]


class _FakeMultiMoleculeTopology:
    """Minimal OpenFF topology double carrying multiple molecules."""

    box_vectors = None

    def __init__(self, molecules: list[_FakeMolecule]):
        self.molecules = molecules

    def get_positions(self):
        """Return flattened source positions in nanometers."""
        atoms = [atom for molecule in self.molecules for atom in molecule.atoms]
        return [[float(index), 0.0, 0.0] for index, _atom in enumerate(atoms)]


class _FakeSageMolecule:
    """Minimal molecule carrying domain and partial-charge provenance."""

    def __init__(
        self,
        name: str,
        *,
        domain: str,
        charges: list[float] | None,
        metadata: dict[str, object] | None = None,
    ):
        self.metadata = {"component_id": name, **(metadata or {})}
        self.force_field_domain = domain
        self.partial_charges = charges


class _FakeDomainAtom:
    """Minimal OpenFF atom double carrying residue metadata."""

    def __init__(self, residue_name: str, atom_name: str):
        self.metadata = {"residue_name": residue_name, "atom_name": atom_name}


def _complete_linkage_audit(index: int) -> dict[str, object]:
    """Return a complete linkage audit dictionary for validation tests."""
    return {
        "index": index,
        "adjacent_terms": {
            "bonds": [{"category": "crosslink_bond"}],
            "angles": [{"category": "asx_side_angle"}, {"category": "glycan_side_angle"}],
            "torsions": [{"category": "proper_crosslink_torsion"}],
        },
        "local_exclusions_and_14_exceptions": [
            {"category": "zero_12_exclusion"},
            {"category": "zero_13_exclusion"},
            {"category": "scaled_14_exception"},
        ],
    }


class _FakeAuditResidue:
    """Minimal OpenMM residue double for audit tests."""

    def __init__(self, name: str, residue_id: str, chain_id: str):
        self.name = name
        self.id = residue_id
        self.chain = SimpleNamespace(id=chain_id)
        self.insertionCode = ""


class _FakeAuditAtom:
    """Minimal OpenMM atom double for audit tests."""

    def __init__(self, index: int, name: str, residue: _FakeAuditResidue):
        self.index = index
        self.name = name
        self.residue = residue


class _FakeAuditTopology:
    """Minimal OpenMM topology double for native audit tests."""

    def __init__(self):
        asx = _FakeAuditResidue("ASX", "60", "A")
        glycan = _FakeAuditResidue("4YB", "1", "C")
        self.atoms_ = [_FakeAuditAtom(0, "ND2", asx), _FakeAuditAtom(1, "C1", glycan)]
        self.crosslink = (self.atoms_[0], self.atoms_[1])

    def atoms(self):
        """Return topology atoms."""
        return iter(self.atoms_)

    def residues(self):
        """Return topology residues."""
        return iter([self.atoms_[0].residue, self.atoms_[1].residue])

    def bonds(self):
        """Return topology bonds."""
        return iter([self.crosslink])

    def getNumAtoms(self):  # noqa: N802 - OpenMM API name
        """Return atom count."""
        return len(self.atoms_)

    def getNumBonds(self):  # noqa: N802 - OpenMM API name
        """Return bond count."""
        return 1


class _FakeNonbonded:
    """Minimal NonbondedForce double for native audit tests."""

    def getNumParticles(self):  # noqa: N802 - OpenMM API name
        """Return particle count."""
        return 2

    def getParticleParameters(self, index):  # noqa: N802 - OpenMM API name
        """Return zero charge parameters with unit-compatible values."""
        from openmm import unit

        _ = index
        return (0.0 * unit.elementary_charge, 1.0 * unit.nanometer, 0.0 * unit.kilojoule_per_mole)

    def getNonbondedMethod(self):  # noqa: N802 - OpenMM API name
        """Return PME method code."""
        return 4

    def getCutoffDistance(self):  # noqa: N802 - OpenMM API name
        """Return cutoff distance."""
        from openmm import unit

        return 1.0 * unit.nanometer

    def getEwaldErrorTolerance(self):  # noqa: N802 - OpenMM API name
        """Return Ewald tolerance."""
        return 0.0005

    def getUseDispersionCorrection(self):  # noqa: N802 - OpenMM API name
        """Return dispersion correction flag."""
        return True

    def getNumExceptions(self):  # noqa: N802 - OpenMM API name
        """Return exception count."""
        return 0


class _FakeSystem:
    """Minimal OpenMM system double for native audit tests."""

    def getNumParticles(self):  # noqa: N802 - OpenMM API name
        """Return particle count."""
        return 2

    def getNumConstraints(self):  # noqa: N802 - OpenMM API name
        """Return constraint count."""
        return 0

    def getForces(self):  # noqa: N802 - OpenMM API name
        """Return no bonded forces for this audit-only test."""
        return []


def _fake_native_source_topology(*, asx_hydrogens: tuple[str, ...]) -> _FakeTopology:
    """Build a minimal ASX--4YB source topology with configurable ASX H names."""
    nd2 = _FakeAtom("ND2", "ASX", "60", 7)
    asx_c = _FakeAtom("C", "ASX", "60", 6)
    hydrogens = [_FakeAtom(name, "ASX", "60", 1) for name in asx_hydrogens]
    c1 = _FakeAtom("C1", "4YB", "1", 6)
    o4 = _FakeAtom("O4", "4YB", "1", 8)
    atoms = [nd2, asx_c, *hydrogens, c1, o4]
    bonds = [_FakeBond(nd2, asx_c), _FakeBond(nd2, c1), _FakeBond(c1, o4)]
    for hydrogen in hydrogens:
        bonds.append(_FakeBond(nd2, hydrogen))
    return _FakeTopology(_FakeMolecule(atoms, bonds))


def _fake_native_source_topology_with_sage(
    *sage_molecules: _FakeMolecule,
) -> _FakeMultiMoleculeTopology:
    """Build an ASX--glycan source topology plus disconnected Sage molecules."""
    base = _fake_native_source_topology(asx_hydrogens=("HD22",)).molecules[0]
    return _FakeMultiMoleculeTopology([*sage_molecules, base])


def _water_molecule(residue_number: int) -> _FakeMolecule:
    """Return a canonical chain-D water molecule double."""

    oxygen = _FakeAtom("O", "HOH", str(residue_number), 8)
    h1 = _FakeAtom("H1", "HOH", str(residue_number), 1)
    h2 = _FakeAtom("H2", "HOH", str(residue_number), 1)
    for atom in (oxygen, h1, h2):
        atom.metadata["chain_id"] = "D"
    return _FakeMolecule([oxygen, h1, h2], [_FakeBond(oxygen, h1), _FakeBond(oxygen, h2)])


def _multi_residue_sage_molecule(start_residue_number: int = 1) -> _FakeMolecule:
    """Return a two-residue disconnected Sage molecule double."""
    c1 = _FakeAtom("C1", "SBM", str(start_residue_number), 6)
    h1 = _FakeAtom("H1", "SBM", str(start_residue_number), 1)
    c2 = _FakeAtom("C2", "EGP", str(start_residue_number + 1), 6)
    h2 = _FakeAtom("H2", "EGP", str(start_residue_number + 1), 1)
    for atom in (c1, h1, c2, h2):
        atom.metadata["force_field_domain"] = "sage"
    molecule = _FakeMolecule(
        [c1, h1, c2, h2], [_FakeBond(c1, h1), _FakeBond(c1, c2), _FakeBond(c2, h2)]
    )
    molecule.force_field_domain = "sage"
    molecule.partial_charges = [0.0, 0.0, 0.0, 0.0]
    molecule.metadata = {"component_id": "polymer-copy"}
    return molecule
