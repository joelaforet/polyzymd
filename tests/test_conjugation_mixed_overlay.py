"""Tests for mixed GLYCAM/OpenFF force-field overlay helpers."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.force_fields import resolve_conjugate_force_fields
from polyzymd.builders.conjugation.ownership import build_particle_ownership_manifest
from polyzymd.builders.conjugation.system_overlay import merge_openmm_system_overlay
from polyzymd.config.schema import ConjugationMoietyConfig


def test_resolver_mixed_overlay_and_no_fallback() -> None:
    """Resolve mixed attachment sources and reject unknown labels."""

    config = _config(
        _attachment("glycan", "glycam_06"),
        _attachment("polymer", "openff-2.0.0.offxml"),
    )
    resolved = resolve_conjugate_force_fields(config)
    assert resolved.route == "mixed_overlay"
    assert resolved.attachments[0].source == "glycam06"
    assert resolved.attachments[1].source == "openff-2.0.0.offxml"
    assert "inherited" not in resolved.attachments[1].to_dict()
    assert "default_small_molecule_force_field" not in resolved.to_dict()

    with pytest.raises(ValueError, match="must explicitly declare moiety.force_field"):
        resolve_conjugate_force_fields(_config(_attachment("missing", None)))

    with pytest.raises(ValueError, match="Unknown conjugation moiety force_field"):
        resolve_conjugate_force_fields(_config(_attachment("bad", "glycan-special")))


@pytest.mark.parametrize(
    ("sources", "expected_route"),
    [
        (("glycam06", "glycam-06"), "native_exact"),
        (("openff-2.0.0.offxml", "/tmp/custom.offxml"), "standard_interchange"),
        (("glycam06", "openff-2.0.0.offxml"), "mixed_overlay"),
    ],
)
def test_route_is_derived_only_from_explicit_attachment_owners(
    sources: tuple[str, str], expected_route: str
) -> None:
    """Explicit all-GLYCAM, all-OpenFF, and mixed owners select the route."""

    config = _config(
        *(_attachment(f"attachment-{index}", source) for index, source in enumerate(sources))
    )

    assert resolve_conjugate_force_fields(config).route == expected_route


def test_disabled_conjugation_is_inert_without_attachment_owner() -> None:
    """A disabled conjugation block does not resolve attachment ownership."""

    config = _config(_attachment("disabled", None))
    config.conjugation.enabled = False

    resolved = resolve_conjugate_force_fields(config)

    assert resolved.route == "standard_interchange"
    assert resolved.attachments == ()


def test_schema_canonicalizes_moiety_force_field_and_forbids_legacy_fields() -> None:
    """Schema validation should normalize approved values and reject unknown fields."""

    inert = ConjugationMoietyConfig.model_validate({"name": "polymer"})
    assert inert.force_field is None
    glycam = ConjugationMoietyConfig.model_validate({"name": "glycan", "force_field": "GLYCAM_06"})
    assert glycam.force_field == "glycam06"
    offxml = ConjugationMoietyConfig.model_validate(
        {"name": "polymer", "force_field": "openff-2.0.0.offxml"}
    )
    assert offxml.force_field == "openff-2.0.0.offxml"
    with pytest.raises(ValueError, match="Unknown moiety.force_field"):
        ConjugationMoietyConfig.model_validate({"name": "bad", "force_field": "glycan-special"})
    with pytest.raises(ValueError, match="must not be blank"):
        ConjugationMoietyConfig.model_validate({"name": "bad", "force_field": "  "})
    with pytest.raises(ValueError, match="Extra inputs are not permitted"):
        ConjugationMoietyConfig.model_validate({"name": "bad", "glycan_force_field": "glycam06"})


def test_ownership_manifest_requires_one_owner_per_particle() -> None:
    """Particle ownership must be complete and conflict-free."""

    manifest = build_particle_ownership_manifest(particle_count=3, glycam_particles={1})
    assert manifest.domains == {"glycam": [1], "generic": [0, 2]}
    with pytest.raises(ValueError, match="out of range"):
        build_particle_ownership_manifest(particle_count=3, glycam_particles={3})


def test_overlay_replaces_particles_terms_constraints_and_exceptions() -> None:
    """Mixed overlay should keep generic terms and replace GLYCAM-touching terms."""

    baseline = _tiny_system(charge_shift=0.0, bond_length=0.10)
    native = _tiny_system(charge_shift=1.0, bond_length=0.20)

    result = merge_openmm_system_overlay(
        baseline_system=baseline,
        native_system=native,
        glycam_particles=frozenset({1}),
    )
    system = result.system
    nonbonded = _force(system, "NonbondedForce")
    assert _charge(nonbonded, 0) == pytest.approx(0.1)
    assert _charge(nonbonded, 1) == pytest.approx(1.2)
    assert _force(system, "HarmonicBondForce").getNumBonds() == 2
    assert _force(system, "HarmonicAngleForce").getNumAngles() == 1
    assert _force(system, "PeriodicTorsionForce").getNumTorsions() == 1
    assert system.getNumConstraints() == 2
    assert result.ownership_manifest.domains["glycam"] == [1]
    assert result.diagnostics["forces"]["NonbondedForce"]["added_native_glycam_exceptions"] == 2
    assert result.diagnostics["parity"]["glycam_vs_native"]["particles"]["max_abs_delta"] == 0.0
    assert result.diagnostics["parity"]["generic_vs_baseline"]["bond"]["max_abs_delta"] == 0.0


def test_overlay_preserves_baseline_nonbonded_globals_and_force_metadata() -> None:
    """Merged handled forces should retain baseline PME and force-level metadata."""

    import openmm
    from openmm import unit

    baseline = _tiny_system(charge_shift=0.0, bond_length=0.10)
    native = _tiny_system(charge_shift=1.0, bond_length=0.20)
    baseline_nb = _force(baseline, "NonbondedForce")
    native_nb = _force(native, "NonbondedForce")
    baseline_nb.setName("baseline nonbonded")
    baseline_nb.setForceGroup(3)
    baseline_nb.setNonbondedMethod(openmm.NonbondedForce.PME)
    baseline_nb.setCutoffDistance(1.1 * unit.nanometer)
    baseline_nb.setUseSwitchingFunction(True)
    baseline_nb.setSwitchingDistance(0.9 * unit.nanometer)
    baseline_nb.setUseDispersionCorrection(False)
    baseline_nb.setEwaldErrorTolerance(2.0e-5)
    baseline_nb.setPMEParameters(0.31 / unit.nanometer, 16, 18, 20)
    baseline_nb.setReciprocalSpaceForceGroup(6)
    native_nb.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
    bond_force = _force(baseline, "HarmonicBondForce")
    bond_force.setName("baseline bonds")
    bond_force.setForceGroup(4)
    bond_force.setUsesPeriodicBoundaryConditions(True)

    result = merge_openmm_system_overlay(
        baseline_system=baseline,
        native_system=native,
        glycam_particles=frozenset({1}),
    )

    merged_nb = _force(result.system, "NonbondedForce")
    assert merged_nb.getName() == "baseline nonbonded"
    assert merged_nb.getForceGroup() == 3
    assert merged_nb.getNonbondedMethod() == openmm.NonbondedForce.PME
    assert merged_nb.getUseSwitchingFunction() is True
    assert merged_nb.getUseDispersionCorrection() is False
    assert merged_nb.getReciprocalSpaceForceGroup() == 6
    assert tuple(merged_nb.getPMEParameters()[1:]) == (16, 18, 20)
    assert result.diagnostics["forces"]["NonbondedForce"]["global_settings"]["method"] == int(
        openmm.NonbondedForce.PME
    )
    assert result.diagnostics["forces"]["NonbondedForce"]["force_metadata"] == {
        "name": "baseline nonbonded",
        "force_group": 3,
        "uses_pbc": True,
    }
    merged_bonds = _force(result.system, "HarmonicBondForce")
    assert merged_bonds.getName() == "baseline bonds"
    assert merged_bonds.getForceGroup() == 4
    assert merged_bonds.usesPeriodicBoundaryConditions() is True
    assert result.diagnostics["forces"]["HarmonicBondForce"]["force_metadata"] == {
        "name": "baseline bonds",
        "force_group": 4,
        "uses_pbc": True,
    }


def test_overlay_rejects_custom_baseline_force_touching_glycam() -> None:
    """Unsupported baseline custom forces may not cross into GLYCAM-owned atoms."""

    import openmm

    baseline = _tiny_system(charge_shift=0.0, bond_length=0.10)
    native = _tiny_system(charge_shift=1.0, bond_length=0.20)
    custom = openmm.CustomBondForce("k*r")
    custom.addPerBondParameter("k")
    custom.addBond(1, 2, [1.0])
    baseline.addForce(custom)

    with pytest.raises(ValueError, match="touches GLYCAM-owned atoms"):
        merge_openmm_system_overlay(
            baseline_system=baseline,
            native_system=native,
            glycam_particles=frozenset({1}),
        )


def test_overlay_allows_generic_only_and_global_custom_baseline_forces() -> None:
    """Provably generic-only and global baseline custom forces should be retained."""

    import openmm

    baseline = _tiny_system(charge_shift=0.0, bond_length=0.10)
    native = _tiny_system(charge_shift=1.0, bond_length=0.20)
    generic_custom = openmm.CustomBondForce("k*r")
    generic_custom.addPerBondParameter("k")
    generic_custom.addBond(0, 2, [1.0])
    global_custom = openmm.CustomExternalForce("g")
    global_custom.addGlobalParameter("g", 0.0)
    baseline.addForce(generic_custom)
    baseline.addForce(global_custom)

    result = merge_openmm_system_overlay(
        baseline_system=baseline,
        native_system=native,
        glycam_particles=frozenset({1}),
    )

    names = [force.__class__.__name__ for force in result.system.getForces()]
    assert names.count("CustomBondForce") == 1
    assert names.count("CustomExternalForce") == 1
    assert result.diagnostics["forces"]["CustomBondForce"]["action"] == "copy"


def test_overlay_rejects_count_only_particle_fallback() -> None:
    """Overlay mapping must be explicit if particle counts differ."""

    baseline = _tiny_system(charge_shift=0.0, bond_length=0.10)
    native = _tiny_system(charge_shift=1.0, bond_length=0.20, particle_count=4)
    with pytest.raises(ValueError, match="atom mapping is required"):
        merge_openmm_system_overlay(
            baseline_system=baseline,
            native_system=native,
            glycam_particles=frozenset({1}),
        )


def test_overlay_accepts_native_subset_mapping_and_preserves_baseline_only() -> None:
    """Native subset mapping should preserve baseline-only Sage or solvent particles."""

    baseline = _tiny_system(charge_shift=0.0, bond_length=0.10, particle_count=4)
    native = _tiny_system(charge_shift=1.0, bond_length=0.20, particle_count=3)

    result = merge_openmm_system_overlay(
        baseline_system=baseline,
        native_system=native,
        glycam_particles=frozenset({1}),
        atom_mapping={0: 0, 1: 1, 2: 2},
    )

    assert result.system.getNumParticles() == 4
    nonbonded = _force(result.system, "NonbondedForce")
    assert nonbonded.getNumParticles() == 4
    assert _charge(nonbonded, 1) == pytest.approx(1.2)
    assert _charge(nonbonded, 3) == pytest.approx(0.4)
    assert result.ownership_manifest.atom_mapping == {0: 0, 1: 1, 2: 2}
    assert result.ownership_manifest.domains["glycam"] == [1]
    assert result.diagnostics["parity"]["unmapped_glycam_count"] == 0


def test_overlay_parity_fails_when_glycam_term_is_tampered(monkeypatch: pytest.MonkeyPatch) -> None:
    """Parity audit should reject a copied GLYCAM term that differs from native."""

    from openmm import unit

    import polyzymd.builders.conjugation.system_overlay as overlay_module

    baseline = _tiny_system(charge_shift=0.0, bond_length=0.10)
    native = _tiny_system(charge_shift=1.0, bond_length=0.20)
    original = overlay_module._add_terms_from

    def tampering_add_terms(*args, **kwargs):
        count = original(*args, **kwargs)
        owner = args[5]
        kind = args[6]
        if owner == "glycam" and kind == "bond":
            force = args[0]
            i, j, _, k = force.getBondParameters(force.getNumBonds() - 1)
            force.setBondParameters(force.getNumBonds() - 1, i, j, 0.21 * unit.nanometer, k)
        return count

    monkeypatch.setattr(overlay_module, "_add_terms_from", tampering_add_terms)
    with pytest.raises(ValueError, match="GLYCAM bond terms parity failed"):
        merge_openmm_system_overlay(
            baseline_system=baseline,
            native_system=native,
            glycam_particles=frozenset({1}),
        )


def test_overlay_parity_fails_when_generic_term_is_tampered(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Parity audit should reject a retained Sage term that differs from baseline."""

    from openmm import unit

    import polyzymd.builders.conjugation.system_overlay as overlay_module

    baseline = _tiny_system(charge_shift=0.0, bond_length=0.10)
    native = _tiny_system(charge_shift=1.0, bond_length=0.20)
    original = overlay_module._add_terms_from

    def tampering_add_terms(*args, **kwargs):
        count = original(*args, **kwargs)
        owner = args[5]
        kind = args[6]
        force = args[0]
        if owner == "generic" and kind == "bond" and force.getNumBonds() > 0:
            i, j, _, k = force.getBondParameters(0)
            force.setBondParameters(0, i, j, 0.11 * unit.nanometer, k)
        return count

    monkeypatch.setattr(overlay_module, "_add_terms_from", tampering_add_terms)
    with pytest.raises(ValueError, match="generic bond terms parity failed"):
        merge_openmm_system_overlay(
            baseline_system=baseline,
            native_system=native,
            glycam_particles=frozenset({2}),
        )


def _attachment(name: str, force_field: str | None) -> SimpleNamespace:
    """Return a minimal attachment object."""

    return SimpleNamespace(
        name=name,
        enabled=True,
        moiety=SimpleNamespace(name=name, force_field=force_field),
    )


def _config(*attachments: SimpleNamespace) -> SimpleNamespace:
    """Return a minimal force-field resolver config."""

    return SimpleNamespace(
        force_field=SimpleNamespace(
            protein="ff14sb_off_impropers_0.0.4.offxml",
            small_molecule="openff-2.0.0.offxml",
        ),
        conjugation=SimpleNamespace(attachments=attachments),
    )


def _tiny_system(
    *,
    charge_shift: float,
    bond_length: float,
    particle_count: int = 3,
) -> object:
    """Build a tiny OpenMM System with all overlay-supported force types."""

    import openmm
    from openmm import unit

    system = openmm.System()
    for _ in range(particle_count):
        system.addParticle(12.0 * unit.dalton)
    if particle_count >= 3:
        system.addConstraint(0, 1, bond_length * unit.nanometer)
        system.addConstraint(1, 2, bond_length * unit.nanometer)

    nonbonded = openmm.NonbondedForce()
    nonbonded.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
    for index in range(particle_count):
        nonbonded.addParticle(
            (0.1 * (index + 1) + charge_shift) * unit.elementary_charge,
            0.3 * unit.nanometer,
            0.5 * unit.kilojoule_per_mole,
        )
    if particle_count >= 3:
        nonbonded.addException(0, 1, 0.0, 0.3 * unit.nanometer, 0.0)
        nonbonded.addException(1, 2, 0.0, 0.3 * unit.nanometer, 0.0)
    system.addForce(nonbonded)

    bonds = openmm.HarmonicBondForce()
    if particle_count >= 3:
        bond_k = 100.0 * unit.kilojoule_per_mole / unit.nanometer**2
        bonds.addBond(0, 1, bond_length * unit.nanometer, bond_k)
        bonds.addBond(1, 2, bond_length * unit.nanometer, bond_k)
    system.addForce(bonds)

    angles = openmm.HarmonicAngleForce()
    if particle_count >= 3:
        angle_k = 10.0 * unit.kilojoule_per_mole / unit.radian**2
        angles.addAngle(0, 1, 2, 1.0 * unit.radian, angle_k)
    system.addForce(angles)

    torsions = openmm.PeriodicTorsionForce()
    if particle_count >= 3:
        torsions.addTorsion(0, 1, 2, 0, 1, 0.0 * unit.radian, 1.0 * unit.kilojoule_per_mole)
    system.addForce(torsions)
    return system


def _force(system: object, name: str) -> object:
    """Return one force by class name."""

    for force in system.getForces():
        if force.__class__.__name__ == name:
            return force
    raise AssertionError(f"Missing force {name}")


def _charge(nonbonded: object, index: int) -> float:
    """Return a particle charge in elementary charge."""

    from openmm import unit

    charge, _, _ = nonbonded.getParticleParameters(index)
    return float(charge.value_in_unit(unit.elementary_charge))
