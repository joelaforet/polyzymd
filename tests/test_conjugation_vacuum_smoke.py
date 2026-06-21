"""Tests for restrained vacuum smoke validation helpers."""

from __future__ import annotations

import json
import math
import os
import shutil

import numpy as np
import pytest
from pydantic import ValidationError

from polyzymd.builders.conjugation._relaxation import (
    VacuumSmokeSettings,
    _positions_to_numpy,
    _resolve_restrained_indices,
    _write_vacuum_smoke_failure,
    analyze_pre_smoke_geometry,
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
        self.residue = type("ResidueDouble", (), {"chain": chain})()


def test_validate_finite_energy_accepts_finite_values():
    """Finite energies should pass smoke validation."""
    validate_finite_energy(-123.4, label="test_energy")


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
