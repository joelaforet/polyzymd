"""Tests for restrained vacuum smoke validation helpers."""

from __future__ import annotations

import json
import math
import os
import shutil
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from pydantic import ValidationError

import polyzymd.builders.conjugation._relaxation as relaxation
from polyzymd.builders.conjugation._relaxation import (
    FrozenProteinRelaxationDiagnostics,
    FrozenProteinRelaxationSettings,
    ProductLinkagePair,
    VacuumSmokeSettings,
    _add_linkage_anchor_restraints,
    _positions_to_numpy,
    _remove_barostats,
    _resolve_restrained_indices,
    _run_frozen_product_md,
    _write_vacuum_smoke_failure,
    analyze_pre_smoke_geometry,
    freeze_chain_a_masses,
    freeze_protein_chain_masses,
    resolve_product_linkage_pairs,
    restore_particle_masses,
    run_frozen_protein_product_relaxation,
    validate_finite_energy,
    validate_finite_positions,
)


class _TopologyDouble:
    """Minimal OpenMM topology-like object for restraint selection tests."""

    def __init__(self, atoms: tuple[object, ...]) -> None:
        """Store atom-like objects returned by ``atoms()``."""
        self._atoms = atoms

    def atoms(self) -> tuple[object, ...]:
        """Return atom-like objects in topology order."""
        return self._atoms

    def bonds(self) -> tuple[tuple[object, object], ...]:
        """Return atom-like topology bonds."""
        return getattr(self, "_bonds", ())


class _AtomDouble:
    """Minimal OpenMM atom-like object for restraint selection tests."""

    def __init__(self, index: int, name: str, element: str, chain_id: str) -> None:
        """Store the atom attributes used by restraint inference."""
        self.index = index
        self.name = name
        self.element = type("ElementDouble", (), {"symbol": element})()
        chain = type("ChainDouble", (), {"id": chain_id})()
        self.residue = type(
            "ResidueDouble",
            (),
            {"chain": chain, "name": "RES", "id": "1"},
        )()


def _pdb_line(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: float,
) -> str:
    """Build a compact PDB atom line for product-linkage tests."""
    element = atom_name[0]
    return (
        f"ATOM  {serial:5d} {atom_name:<4s} {residue_name:>3s} {chain_id:1s}"
        f"{residue_number:4d}    {x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00"
        f"          {element:>2s}\n"
    )


def _topology_atom(
    index: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_id: str,
) -> object:
    """Create a topology atom double with product identity fields."""
    chain = type("ChainDouble", (), {"id": chain_id})()
    residue = type("ResidueDouble", (), {"chain": chain, "name": residue_name, "id": residue_id})()
    return type("TopologyAtomDouble", (), {"index": index, "name": atom_name, "residue": residue})()


def test_validate_finite_energy_accepts_finite_values():
    """Finite energies should pass smoke validation."""
    validate_finite_energy(-123.4, label="test_energy")


def test_resolve_product_linkage_pairs_uses_generic_assembly_metadata(tmp_path):
    """Arbitrary residue and atom names should resolve from CONECT metadata."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(1, "Q1", "ABC", "A", 7, 0.0),
                _pdb_line(2, "Z9", "MNO", "C", 42, 1.4),
                "CONECT    1    2\nEND\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "Q1", "ABC", "A", "7"),
            _topology_atom(1, "Z9", "MNO", "C", "42"),
        )
    )
    plan = type("PlanDouble", (), {"target_bond_length_angstrom": 1.4})()
    spec = type(
        "SpecDouble",
        (),
        {"resolved_plan": plan, "attachment_id": "x", "attachment_index": 1},
    )()
    assembly = type("AssemblyDouble", (), {"added_conect_pairs": ((1, 2),)})()

    pairs = resolve_product_linkage_pairs(
        topology,
        product_pdb_path=product,
        attachment_specs=(spec,),
        assembly=assembly,
    )

    assert pairs[0].protein_atom_index == 0
    assert pairs[0].modifier_atom_index == 1
    assert pairs[0].target_bond_length_angstrom == pytest.approx(1.4)


def test_resolve_product_linkage_pairs_disambiguates_duplicate_moieties(tmp_path):
    """Identical duplicated moieties should resolve distinct product linkage pairs."""
    product = tmp_path / "product.pdb"
    product.write_text(
        "".join(
            (
                _pdb_line(1, "Q1", "ABC", "A", 7, 0.0),
                _pdb_line(2, "Z9", "MNO", "C", 42, 1.4),
                _pdb_line(3, "Q1", "ABC", "A", 8, 3.0),
                _pdb_line(4, "Z9", "MNO", "C", 43, 4.4),
                "END\n",
            )
        ),
        encoding="utf-8",
    )
    topology = _TopologyDouble(
        (
            _topology_atom(0, "Q1", "ABC", "A", "7"),
            _topology_atom(1, "Z9", "MNO", "C", "42"),
            _topology_atom(2, "Q1", "ABC", "A", "8"),
            _topology_atom(3, "Z9", "MNO", "C", "43"),
        )
    )
    plan = type("PlanDouble", (), {"target_bond_length_angstrom": 1.4})()
    specs = tuple(
        type(
            "SpecDouble",
            (),
            {"resolved_plan": plan, "attachment_id": f"x{index}", "attachment_index": index},
        )()
        for index in (1, 2)
    )
    assembly = type("AssemblyDouble", (), {"added_conect_pairs": ((1, 2), (3, 4))})()

    pairs = resolve_product_linkage_pairs(
        topology,
        product_pdb_path=product,
        attachment_specs=specs,
        assembly=assembly,
    )

    assert [(pair.protein_atom_index, pair.modifier_atom_index) for pair in pairs] == [
        (0, 1),
        (2, 3),
    ]


def test_freeze_chain_a_masses_restores_only_protein_chain():
    """Temporary zero masses should affect only chain A particles."""

    class SystemDouble:
        """Minimal system double with mutable particle masses."""

        def __init__(self) -> None:
            self.masses = [12.0, 14.0, 16.0]

        def getParticleMass(self, index: int) -> float:
            """Return one particle mass."""
            return self.masses[index]

        def setParticleMass(self, index: int, mass: float) -> None:
            """Set one particle mass."""
            self.masses[index] = mass

    topology = _TopologyDouble(
        (
            _AtomDouble(0, "CA", "C", "A"),
            _AtomDouble(1, "CB", "C", "A"),
            _AtomDouble(2, "C1", "C", "C"),
        )
    )
    system = SystemDouble()

    frozen, original = freeze_chain_a_masses(system, topology, type("Unit", (), {"dalton": 1.0})())
    assert frozen == (0, 1)
    assert system.masses == [0.0, 0.0, 16.0]

    restore_particle_masses(system, original)
    assert system.masses == [12.0, 14.0, 16.0]


def test_freeze_protein_chain_masses_supports_configurable_chains():
    """Temporary zero masses should support generic configured protein chains."""

    class SystemDouble:
        """Minimal system double with mutable particle masses."""

        def __init__(self) -> None:
            """Store particle masses."""
            self.masses = [12.0, 14.0, 16.0]

        def getParticleMass(self, index: int) -> float:
            """Return one particle mass."""
            return self.masses[index]

        def setParticleMass(self, index: int, mass: float) -> None:
            """Set one particle mass."""
            self.masses[index] = mass

    topology = _TopologyDouble(
        (
            _AtomDouble(0, "CA", "C", "A"),
            _AtomDouble(1, "CB", "C", "B"),
            _AtomDouble(2, "C1", "C", "C"),
        )
    )
    system = SystemDouble()

    frozen, original = freeze_protein_chain_masses(
        system,
        topology,
        type("Unit", (), {"dalton": 1.0})(),
        chain_ids=("B",),
    )

    assert frozen == (1,)
    assert system.masses == [12.0, 0.0, 16.0]
    restore_particle_masses(system, original)
    assert system.masses == [12.0, 14.0, 16.0]


def test_temporary_anchor_restraint_count_matches_linkages():
    """Temporary anchor restraints should add one bond per attachment."""

    class ForceDouble:
        """CustomBondForce test double."""

        def __init__(self, expression: str) -> None:
            self.expression = expression
            self.bonds = []

        def addPerBondParameter(self, name: str) -> None:
            """Accept per-bond parameter declarations."""

        def addGlobalParameter(self, name: str, value: float) -> None:
            """Accept global parameter declarations."""

        def addBond(self, atom_i: int, atom_j: int, parameters: list[float]) -> None:
            """Store one anchor bond."""
            self.bonds.append((atom_i, atom_j, parameters))

    class SystemDouble:
        """System test double collecting added forces."""

        def __init__(self) -> None:
            self.forces = []

        def addForce(self, force: ForceDouble) -> None:
            """Store one added force."""
            self.forces.append(force)

    openmm = type("OpenMMDouble", (), {"CustomBondForce": ForceDouble})()
    unit = type("UnitDouble", (), {"kilojoule_per_mole": 1.0, "nanometer": 1.0})()
    pairs = (
        ProductLinkagePair(
            protein_atom_index=0, modifier_atom_index=1, target_bond_length_angstrom=1.3
        ),
        ProductLinkagePair(
            protein_atom_index=2, modifier_atom_index=3, target_bond_length_angstrom=1.4
        ),
    )
    system = SystemDouble()

    count = _add_linkage_anchor_restraints(system, pairs, 100.0, openmm, unit)

    assert count == 2
    assert len(system.forces[0].bonds) == 2


def test_frozen_relaxation_stage_b_uses_fresh_frozen_system(monkeypatch, tmp_path):
    """Stage A should use normal masses and Stage B should use a fresh frozen system."""
    calls: list[tuple[str, str, tuple[float, ...]]] = []
    openmm, openmm_app = _install_fake_openmm(monkeypatch, calls)
    monkeypatch.setattr(
        relaxation, "_select_platform", lambda *_args: SimpleNamespace(getName=lambda: "CPU")
    )
    monkeypatch.setattr(relaxation, "_assign_force_groups", lambda _system: {})
    monkeypatch.setattr(relaxation, "_force_group_labels", lambda _system, *, existing_labels: {})
    monkeypatch.setattr(relaxation, "_force_group_energies", lambda *_args: {})
    monkeypatch.setattr(relaxation, "_add_linkage_anchor_restraints", lambda *_args: 0)
    monkeypatch.setattr(relaxation, "_write_openmm_pdb", lambda *_args: None)

    topology = _relaxation_topology()
    interchange = _RelaxationInterchange(topology)

    result = run_frozen_protein_product_relaxation(
        interchange,
        tmp_path,
        product_pdb_path=tmp_path / "product.pdb",
        attachment_specs=(),
        settings=FrozenProteinRelaxationSettings(md_steps=5),
    )

    assert result.success is True
    assert len(interchange.systems) == 2
    assert interchange.systems[0] is not interchange.systems[1]
    assert calls == [
        ("stage_a", "init", (12.0, 16.0)),
        ("stage_b", "init", (0.0, 16.0)),
        ("stage_b", "step", (0.0, 16.0)),
    ]
    assert interchange.systems[0].masses == [12.0, 16.0]
    assert interchange.systems[1].masses == [12.0, 16.0]
    assert sys.modules["openmm"] is openmm
    assert sys.modules["openmm.app"] is openmm_app


def test_frozen_relaxation_restores_masses_after_stage_b_error(monkeypatch, tmp_path):
    """Temporary Stage B zero masses should be restored when MD fails."""
    calls: list[tuple[str, str, tuple[float, ...]]] = []
    _install_fake_openmm(monkeypatch, calls)
    monkeypatch.setattr(
        relaxation, "_select_platform", lambda *_args: SimpleNamespace(getName=lambda: "CPU")
    )
    monkeypatch.setattr(relaxation, "_assign_force_groups", lambda _system: {})
    monkeypatch.setattr(relaxation, "_force_group_labels", lambda _system, *, existing_labels: {})
    monkeypatch.setattr(relaxation, "_force_group_energies", lambda *_args: {})
    monkeypatch.setattr(relaxation, "_add_linkage_anchor_restraints", lambda *_args: 0)
    monkeypatch.setattr(relaxation, "_write_openmm_pdb", lambda *_args: None)

    def fail_md(*_args, **_kwargs):
        raise RuntimeError("backend instability")

    monkeypatch.setattr(relaxation, "_run_frozen_product_md", fail_md)
    interchange = _RelaxationInterchange(_relaxation_topology())

    with pytest.raises(RuntimeError, match="backend instability"):
        run_frozen_protein_product_relaxation(
            interchange,
            tmp_path,
            product_pdb_path=tmp_path / "product.pdb",
            attachment_specs=(),
            settings=FrozenProteinRelaxationSettings(md_steps=5),
        )

    assert interchange.systems[1].masses == [12.0, 16.0]


def test_frozen_relaxation_settings_reject_zero_md_steps():
    """Frozen-protein relaxation must run a real Stage B MD segment."""
    with pytest.raises(ValidationError, match="greater than 0"):
        FrozenProteinRelaxationSettings(md_steps=0)


def test_frozen_product_md_propagates_programming_errors(monkeypatch):
    """Retry handling should not mask programming errors as instability."""

    class BuggySimulation:
        def __init__(self, *_args):
            raise TypeError("bad fake API")

    warnings: list[str] = []
    openmm = SimpleNamespace(
        LangevinMiddleIntegrator=lambda *_args: object(),
        OpenMMException=RuntimeError,
    )
    app = SimpleNamespace(Simulation=BuggySimulation)
    unit = SimpleNamespace(kelvin=1.0, picosecond=1.0, femtosecond=1.0, nanometer=1.0)

    with pytest.raises(TypeError, match="bad fake API"):
        _run_frozen_product_md(
            _relaxation_topology(),
            _RelaxationSystem("stage_b"),
            np.zeros((2, 3)),
            FrozenProteinRelaxationSettings(md_steps=5),
            openmm,
            app,
            unit,
            object(),
            warnings,
        )

    assert warnings == []


def test_remove_barostats_removes_only_barostat_forces():
    """Vacuum relaxation should discard barostat forces from transient systems."""

    class MonteCarloBarostat:
        """Barostat marker double."""

    class HarmonicBondForce:
        """Non-barostat marker double."""

    class SystemDouble:
        """System double with removable forces."""

        def __init__(self) -> None:
            self.forces = [HarmonicBondForce(), MonteCarloBarostat(), HarmonicBondForce()]

        def getNumForces(self) -> int:
            """Return force count."""
            return len(self.forces)

        def getForce(self, index: int) -> object:
            """Return force by index."""
            return self.forces[index]

        def removeForce(self, index: int) -> None:
            """Remove force by index."""
            self.forces.pop(index)

    system = SystemDouble()

    assert _remove_barostats(system) == 1
    assert [type(force).__name__ for force in system.forces] == [
        "HarmonicBondForce",
        "HarmonicBondForce",
    ]


def test_frozen_protein_diagnostics_serializes(tmp_path):
    """Frozen-protein diagnostics should write JSON evidence."""
    path = tmp_path / "diagnostics.json"
    diagnostics = FrozenProteinRelaxationDiagnostics(
        success=True,
        stage_a_success=True,
        stage_b_success=True,
        frozen_atom_count=2,
        temporary_anchor_count=1,
        stage_a_energy_after_min_kj_mol=-10.0,
        stage_b_energy_after_md_kj_mol=-11.0,
        stage_a_protein_rmsd_from_initial_angstrom=0.2,
        stage_b_protein_rmsd_from_stage_a_angstrom=0.0,
        linkage_distances_angstrom=(1.4,),
        stage_b_linkage_distances_angstrom=(1.4,),
        final_relaxed_pdb_path=tmp_path / "relaxed.pdb",
    )

    diagnostics.write_json(path)

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["success"] is True
    assert payload["stage_a_success"] is True
    assert payload["stage_b_success"] is True
    assert payload["temporary_anchor_count"] == 1
    assert payload["barostat_used"] is False
    assert payload["final_relaxed_pdb_path"] == str(tmp_path / "relaxed.pdb")
    assert payload["json_path"] == str(path)


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_validate_finite_energy_rejects_nonfinite_values(value: float):
    """NaN and infinite energies should fail hard."""
    with pytest.raises(RuntimeError, match="non-finite test_energy"):
        validate_finite_energy(value, label="test_energy")


def test_validate_finite_positions_accepts_numpy_arrays():
    """Finite coordinate arrays should pass smoke validation."""
    span_nm = validate_finite_positions(np.zeros((3, 3)), label="positions")

    assert span_nm == 0.0


def test_positions_to_numpy_logs_expected_conversion_fallback(caplog):
    """Expected unit conversion failures should fall back with a warning."""

    class RawArrayPositions:
        """Position container with a broken unit API and raw array fallback."""

        def value_in_unit(self, unit):
            """Raise an expected conversion error for the unit-aware path."""
            raise TypeError(f"cannot convert to {unit}")

        def __array__(self, dtype=None):
            """Return raw nanometer coordinates for NumPy coercion."""
            return np.asarray([[0.0, 0.0, 0.0]], dtype=dtype)

    unit_module = type("UnitModule", (), {"nanometer": "nanometer"})()

    with caplog.at_level("WARNING"):
        array = _positions_to_numpy(RawArrayPositions(), unit_module)

    np.testing.assert_allclose(array, [[0.0, 0.0, 0.0]])
    assert "Falling back to raw np.asarray()" in caplog.text
    assert "RawArrayPositions" in caplog.text


def test_validate_finite_positions_propagates_unexpected_conversion_errors():
    """Unexpected position conversion failures should not be hidden."""

    class BrokenPositions:
        """Position container with an unexpected conversion failure."""

        def value_in_unit(self, unit):
            """Raise an unexpected error from the unit-aware path."""
            raise RuntimeError(f"unexpected conversion failure for {unit}")

    unit_module = type("UnitModule", (), {"nanometer": "nanometer"})()

    with pytest.raises(RuntimeError, match="unexpected conversion failure"):
        validate_finite_positions(BrokenPositions(), unit_module, label="positions")


def test_validate_finite_positions_rejects_nonfinite_arrays():
    """NaN coordinates should fail hard before downstream use."""
    positions = np.array([[0.0, 0.0, 0.0], [np.nan, 1.0, 2.0]])

    with pytest.raises(RuntimeError, match="non-finite positions"):
        validate_finite_positions(positions, label="positions")


def test_validate_finite_positions_rejects_unrealistic_span():
    """Blown-up post-MD coordinates should fail before solvation."""
    positions = np.array([[0.0, 0.0, 0.0], [51.0, 1.0, 1.0]])

    with pytest.raises(RuntimeError, match="unrealistic coordinate span"):
        validate_finite_positions(positions, label="equilibrated_positions", max_span_nm=50.0)


def test_vacuum_smoke_settings_require_positive_nvt_steps():
    """Vacuum smoke settings should allow minimization-only smoke."""
    assert VacuumSmokeSettings(nvt_steps=0).nvt_steps == 0


def test_vacuum_smoke_settings_reject_negative_nvt_steps():
    """Vacuum smoke settings should reject negative MD step counts."""
    with pytest.raises(ValidationError, match="greater than or equal to 0"):
        VacuumSmokeSettings(nvt_steps=-1)


def test_vacuum_smoke_settings_use_conservative_md_defaults():
    """Default smoke MD settings should be conservative for vacuum stability."""
    settings = VacuumSmokeSettings()

    assert settings.nvt_steps > 0
    assert settings.timestep_femtoseconds <= 0.25
    assert settings.temperature_kelvin <= 50.0
    assert settings.friction_per_picosecond >= 10.0
    assert settings.restrain_all_heavy_atoms is True
    assert settings.max_position_span_nm == 50.0


def test_analyze_pre_smoke_geometry_reports_contacts_and_bond_outliers():
    """Pre-smoke diagnostics should identify close contacts and stretched bonds."""
    topology = _TopologyDouble(
        (
            _AtomDouble(0, "C1", "C", "C"),
            _AtomDouble(1, "C2", "C", "C"),
            _AtomDouble(2, "H1", "H", "C"),
            _AtomDouble(3, "O1", "O", "C"),
        )
    )
    topology._bonds = ((topology._atoms[0], topology._atoms[3]),)
    positions = np.array(
        [
            [0.00, 0.00, 0.00],
            [0.10, 0.00, 0.00],
            [0.02, 0.00, 0.00],
            [0.40, 0.00, 0.00],
        ]
    )

    diagnostics = analyze_pre_smoke_geometry(topology, positions)

    assert diagnostics.coordinate_span_nm == pytest.approx(0.40)
    assert diagnostics.min_heavy_heavy_distance_nm == pytest.approx(0.10)
    assert diagnostics.min_h_heavy_distance_nm == pytest.approx(0.02)
    assert [pair.category for pair in diagnostics.close_contacts] == [
        "h-heavy-close-contact",
        "heavy-heavy-close-contact",
    ]
    assert diagnostics.bonded_distance_outliers[0].distance_nm == pytest.approx(0.40)


def test_pre_smoke_geometry_writes_json_sidecar(tmp_path):
    """Pre-smoke geometry diagnostics should persist JSON artifacts."""
    topology = _restraint_selection_topology()
    diagnostics = analyze_pre_smoke_geometry(topology, np.zeros((4, 3)))

    path = diagnostics.write_json(tmp_path / "pre_smoke_geometry.json")

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["atom_count"] == 4
    assert payload["json_path"].endswith("pre_smoke_geometry.json")


def test_vacuum_smoke_failure_artifact_records_span_and_energies(tmp_path):
    """Smoke failures should write structured diagnostics before raising upstream."""
    topology = _restraint_selection_topology()
    pre_smoke = analyze_pre_smoke_geometry(
        topology,
        np.array([[0.0, 0.0, 0.0], [60.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]]),
    )

    path = _write_vacuum_smoke_failure(
        tmp_path / "vacuum_smoke_failure.json",
        exc=RuntimeError("span validation failed"),
        settings=VacuumSmokeSettings(nvt_steps=0),
        pre_smoke=pre_smoke,
        energy_before_min=10.0,
        energy_after_min=5.0,
        energy_before_nvt=5.0,
        energy_after_nvt=math.nan,
        failed_pdb_path=tmp_path / "failed.pdb",
    )

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["success"] is False
    assert payload["energies_kj_mol"]["after_min"] == pytest.approx(5.0)
    assert payload["energies_kj_mol"]["after_nvt"] is None
    assert payload["pre_smoke_geometry"]["coordinate_span_nm"] == pytest.approx(60.0)


def test_resolve_restrained_indices_defaults_to_all_heavy_atoms():
    """All-heavy mode should restrain non-protein conjugate atoms by default."""
    topology = _restraint_selection_topology()

    indices = _resolve_restrained_indices(
        topology,
        protein_heavy_atom_indices=None,
        restrain_all_heavy_atoms=True,
    )

    assert indices == (0, 2)


def test_resolve_restrained_indices_ignores_protein_only_selection_in_all_heavy_mode():
    """Supplying chain-A atoms must not restrict all-heavy vacuum restraints."""
    topology = _restraint_selection_topology()

    indices = _resolve_restrained_indices(
        topology,
        protein_heavy_atom_indices=(0,),
        restrain_all_heavy_atoms=True,
    )

    assert indices == (0, 2)


def test_resolve_restrained_indices_allows_legacy_protein_only_mode():
    """Protein-only restraints remain available when all-heavy mode is disabled."""
    topology = _restraint_selection_topology()

    indices = _resolve_restrained_indices(
        topology,
        protein_heavy_atom_indices=None,
        restrain_all_heavy_atoms=False,
    )

    assert indices == (0,)


def _restraint_selection_topology() -> _TopologyDouble:
    """Build a small topology with protein, conjugate, and hydrogen atoms."""
    return _TopologyDouble(
        (
            _AtomDouble(0, "CA", "C", "A"),
            _AtomDouble(1, "HA", "H", "A"),
            _AtomDouble(2, "C1", "C", "C"),
            _AtomDouble(3, "H1", "H", "C"),
        )
    )


class _RelaxationSystem:
    """Mutable system double for staged frozen-protein relaxation tests."""

    def __init__(self, label: str) -> None:
        """Store a label and particle masses."""
        self.label = label
        self.masses = [12.0, 16.0]

    def getNumForces(self) -> int:  # noqa: N802 - OpenMM API compatibility
        """Return no forces for the minimal relaxation path."""
        return 0

    def getParticleMass(self, index: int) -> float:  # noqa: N802 - OpenMM API compatibility
        """Return one particle mass."""
        return self.masses[index]

    def setParticleMass(self, index: int, mass: float) -> None:  # noqa: N802
        """Set one particle mass."""
        self.masses[index] = float(mass)


class _RelaxationInterchange:
    """Interchange double returning fresh OpenMM systems."""

    def __init__(self, topology: _TopologyDouble) -> None:
        """Store topology and emitted systems."""
        self.topology = topology
        self.positions = np.array([[0.0, 0.0, 0.0], [0.15, 0.0, 0.0]])
        self.systems: list[_RelaxationSystem] = []

    def to_openmm_topology(self) -> _TopologyDouble:
        """Return the fake topology."""
        return self.topology

    def to_openmm_system(self) -> _RelaxationSystem:
        """Return a fresh system for each staged OpenMM conversion."""
        label = "stage_a" if not self.systems else "stage_b"
        system = _RelaxationSystem(label)
        self.systems.append(system)
        return system


class _RelaxationState:
    """OpenMM state double with energy and positions."""

    def __init__(self, positions: np.ndarray) -> None:
        """Store state positions."""
        self._positions = np.asarray(positions, dtype=float)

    def getPotentialEnergy(self) -> float:  # noqa: N802 - OpenMM API compatibility
        """Return a finite energy."""
        return -1.0

    def getPositions(self, asNumpy: bool = False) -> np.ndarray:  # noqa: N802, FBT001, FBT002
        """Return the stored positions."""
        return self._positions.copy()


class _RelaxationContext:
    """OpenMM context double storing positions."""

    def __init__(self) -> None:
        """Initialize empty positions."""
        self.positions = np.zeros((2, 3))

    def setPositions(self, positions: np.ndarray) -> None:  # noqa: N802
        """Store context positions."""
        self.positions = np.asarray(positions, dtype=float)

    def setVelocities(self, _velocities: np.ndarray) -> None:  # noqa: N802
        """Accept velocity initialization."""

    def getState(self, **_kwargs) -> _RelaxationState:  # noqa: N802
        """Return a finite state for requested data."""
        return _RelaxationState(self.positions)


def _relaxation_topology() -> _TopologyDouble:
    """Build a minimal two-atom protein-conjugate topology."""
    return _TopologyDouble(
        (
            _AtomDouble(0, "CA", "C", "A"),
            _AtomDouble(1, "C1", "C", "C"),
        )
    )


def _install_fake_openmm(monkeypatch, calls: list[tuple[str, str, tuple[float, ...]]]):
    """Install fake OpenMM modules for staged relaxation behavior tests."""

    class FakeSimulation:
        """OpenMM Simulation double recording stage masses."""

        def __init__(self, _topology, system: _RelaxationSystem, _integrator, _platform) -> None:
            self.system = system
            self.context = _RelaxationContext()
            calls.append((system.label, "init", tuple(system.masses)))

        def step(self, _steps: int) -> None:
            """Record masses used during MD integration."""
            calls.append((self.system.label, "step", tuple(self.system.masses)))

    openmm = SimpleNamespace(
        VerletIntegrator=lambda *_args: object(),
        LangevinMiddleIntegrator=lambda *_args: object(),
        LocalEnergyMinimizer=SimpleNamespace(minimize=lambda *_args, **_kwargs: None),
        OpenMMException=RuntimeError,
    )
    openmm_app = SimpleNamespace(
        Simulation=FakeSimulation,
        PDBFile=SimpleNamespace(writeFile=lambda *_args, **_kwargs: None),
    )
    openmm_unit = SimpleNamespace(
        picoseconds=1.0,
        picosecond=1.0,
        femtosecond=1.0,
        kelvin=1.0,
        kilojoule_per_mole=1.0,
        nanometer=1.0,
        dalton=1.0,
    )
    openmm.app = openmm_app
    openmm.unit = openmm_unit
    monkeypatch.setitem(sys.modules, "openmm", openmm)
    monkeypatch.setitem(sys.modules, "openmm.app", openmm_app)
    monkeypatch.setitem(sys.modules, "openmm.unit", openmm_unit)
    return openmm, openmm_app


@pytest.mark.slow
@pytest.mark.conjugation_stack
def test_opt_in_conjugation_stack_smoke_requirements_are_available():
    """Opt-in guard for stack availability before the integrated smoke test.

    The physics acceptance test lives in
    ``tests/test_conjugation_integrated_smoke.py`` and should be run under::

        module load slurm/blanca
        salloc ...
        PYTHONNOUSERSITE=1 POLYZYMD_RUN_CONJUGATION_PABLO_SMOKE=1 \
          pixi run -e sim-cuda-12-4 pytest \
          tests/test_conjugation_integrated_smoke.py -v
    """
    if os.environ.get("POLYZYMD_RUN_CONJUGATION_PABLO_SMOKE") != "1":
        pytest.skip("Set POLYZYMD_RUN_CONJUGATION_PABLO_SMOKE=1 to check the Pablo smoke stack")
    if shutil.which("packmol") is None:
        pytest.skip("Packmol binary is not available on PATH")
    pytest.importorskip("polymerist")
    pytest.importorskip("openff.toolkit")
    pytest.importorskip("openff.interchange")
    pytest.importorskip("openff.pablo")
    pytest.importorskip("openmm")
